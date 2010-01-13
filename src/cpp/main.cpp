// Copyright (C) 2008-2009 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-07-16
// Last changed: 2009-08-25
//

#include <dolfin.h>
#include <boost/scoped_array.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/if.hpp>

#include "DarcyFlow.h"
#include "SaturationEquation.h"
#include "VelocityP0Projection.h"
#include "VelocityP1Projection.h"
#include "ScalarP1Projection.h"

using namespace dolfin;
using namespace boost::lambda;

// Velocity
/*
class Velocity : public Expression
{
  public:

    Velocity() : Expression(2) {}

    void eval(double* values, const std::vector<double>& x) const
    { 
      values[0] =  1.0; 
      values[1] =  0.0; 
    }
};
*/

// Source term
class Source : public Expression
{
  public:

    void eval(double* values, const std::vector<double>& x) const
    { 
      values[0] =  0.0; 
    }
};

// Snake Permeability 
class SnakePermeabilty : public Expression
{
  public:

    void eval(double* values, const std::vector<double>& x) const
    { 
      // Snake-like crack
      const double l = (x[1] - 0.5 - 0.1*sin(10*x[0]))/0.1;
      values[0] = std::max(exp(-(l*l)), 0.01); 
    }
};

// Snake Permeability 
class RandomPermeabilty : public Expression
{
  public:

    RandomPermeabilty() : x_r(80, std::vector<double>(2))

    {
      for(dolfin::uint i = 0; i < x_r.size(); ++i)
        for(dolfin::uint j = 0; j < 2; ++j)
          x_r[i][j] = dolfin::rand(); 
    }

    void eval(double* values, const std::vector<double>& x) const
    { 
      const dolfin::uint N = 80;
      double l;
      double k1 = 0.01;
      const double k2 = 4.0;
      double sigma= 0.0;
      for(dolfin::uint i = 0; i < N; ++i)
      {
        l = (x[0]-x_r[i][0])*(x[0]-x_r[i][0]) + (x[1]-x_r[i][1])*(x[1]-x_r[i][1]); 
        sigma = std::max(sigma, exp(-l/0.0025));         
      }
      k1 = std::max(sigma, k1);      
      values[0] = std::min(k1, k2); 
    }

  private:
    std::vector< std::vector<double> > x_r;
};

// Pressure boundary condition
class Gp : public Expression
{
  public:

    void eval(double* values, const std::vector<double>& x) const
    { 
      values[0] =  1.0 - x[0]; 
    }
};

// Saturation boundary condition
class Gs : public Expression
{
  public:

    Gs(const double& t, const Constant& dt) : t(t), dt(dt) {}

    void eval(double* values, const std::vector<double>& x) const
    {       
      values[0] = 0.0;
      if(std::abs(x[0]) < DOLFIN_EPS)
      {
        /*
        const double tmp0 = 2.0*sin(DOLFIN_PI*std::pow(0.5, 1.0/3.0));
        double tmp1 = sin(DOLFIN_PI*std::pow(x[1], 1.0/3.0))/tmp0;
        tmp1       += sin(DOLFIN_PI*std::pow(1.0- x[1], 1.0/3.0))/tmp0;
        if(t < 10*dt)
          values[0] = (t/(10.0*dt))*tmp1;
        else 
          values[0] = tmp1;
        */
        if(t < 100*dt)
          values[0] = (t/(100.0*dt));
        else 
          values[0] = 1.0;
      }
    }
  private:
    const double& t;
    const Constant& dt;
};

// Sub domain for inflow boundary condition
class InflowBoundary : public SubDomain
{
  bool inside(const double* x, bool on_boundary) const
  {
    return std::abs(x[0]) < DOLFIN_EPS && on_boundary;
  }
};

// Sub domain for top and bottom boundaries
class SideBoundary : public SubDomain
{
  bool inside(const double* x, bool on_boundary) const
  {
    if (std::abs(x[1]) < DOLFIN_EPS && on_boundary)
      return true;
    else if (std::abs(x[1] - 1.0)  < DOLFIN_EPS && on_boundary)
      return true;
    else
      return false;
  }
};


int main()
{
  dolfin::seed(3);

  // Create mesh and source term
  double hmesh = 32;
  UnitSquare mesh(hmesh, hmesh);
  //UnitSquare mesh(hmesh, hmesh, "crossed");
  CellSize h(mesh);
  
  // Time
  double t  = 0.0;
  Constant dt(0.0);

  // Create mesh function over the cell facets
  MeshFunction<unsigned int> facet_domains(mesh, mesh.topology().dim() - 1);
  facet_domains = 0;
  InflowBoundary inflow_boundary;
  SideBoundary side_boundary;
  inflow_boundary.mark(facet_domains, 1);
  side_boundary.mark(facet_domains, 2);

  // Create Darcy function space and extract velocity subspace
  DarcyFlow::FunctionSpace Vdarcy(mesh);
  std::vector<dolfin::uint> sub_space = boost::assign::list_of(0);
  boost::shared_ptr<const FunctionSpace> Vu(Vdarcy.extract_sub_space(sub_space));

  SaturationEquation::FunctionSpace Vs(mesh);

  // Parameters and functions
  //RandomPermeabilty permeability(mesh);   // permeability
  //SnakePermeabilty permeability;   
  Constant permeability(0.01);   
  Function U(Vdarcy);                    // velocity, pressure
  Function u(Vu);                    // velocity
  Source f;                           // Darcy flow source term
  Function S(Vs), S0(Vs);                           // saturation
  Function F(Vs);                               // Fractional flow function
  Function dFdS(Vs);                            // dF/dS
  Function mobility(Vs);                        // mobility

  // Presssure and saturation boundary conditions
  Gp gp;
  Gs gs(t, dt);

  std::vector<const BoundaryCondition*> empty_bcs;

  // Define Darcy flow PDE
  DarcyFlow::BilinearForm a_darcy(Vdarcy, Vdarcy);
  DarcyFlow::LinearForm L_darcy(Vdarcy);
  a_darcy.k = permeability; a_darcy.mobility = mobility; 
  L_darcy.g = gp; L_darcy.f = f;
  VariationalProblem pde_darcy(a_darcy, L_darcy, empty_bcs, 0, &facet_domains, 0);

  // Saturation PDE
  SaturationEquation::BilinearForm a_s(Vs, Vs);
  SaturationEquation::LinearForm L_s(Vs);
  L_s.w0 = U; L_s.s0 = S0, L_s.F = F; L_s.g = gs; L_s.dt = dt; 
  VariationalProblem pde_saturation(a_s, L_s, empty_bcs, 0, &facet_domains, 0);


  // Projection forms (u onto P1 functions for post-processing)
  VelocityP1Projection::FunctionSpace Vp1(mesh);
  VelocityP1Projection::BilinearForm a_u_p1(Vp1, Vp1);
  VelocityP1Projection::LinearForm L_u_p1(Vp1, U);
  VariationalProblem pde_u_p1(a_u_p1, L_u_p1);

  // Project u onto P0 functions to compute maximum time step
  VelocityP0Projection::FunctionSpace Vp0(mesh);
  VelocityP0Projection::BilinearForm a_dt(Vp0, Vp0);
  VelocityP0Projection::LinearForm L_dt(Vp0, U, dFdS, h);

  // Perform some sanity checks
  assert(mobility.vector().size() == S0.vector().size());
  assert(F.vector().size() == S0.vector().size());
  assert(dFdS.vector().size() == S0.vector().size());

  // Project u onto P0 functions to compute maximum time step
  VariationalProblem pde_dt(a_dt, L_dt);
  pde_dt.parameters["linear_solver"] = "iterative";

  // Project saturation for post-processing
  ScalarP1Projection::FunctionSpace Vsp1(mesh);
  Function saturation_projected(Vsp1);
  ScalarP1Projection::BilinearForm a_sat(Vsp1, Vsp1);
  ScalarP1Projection::LinearForm L_sat(Vsp1, S0);
  VariationalProblem pde_s(a_sat, L_sat);
  pde_s.solve(saturation_projected);

  // Project mobility for post-processing
  Function mobility_projected(Vsp1);
  ScalarP1Projection::BilinearForm a_mob(Vsp1, Vsp1);
  ScalarP1Projection::LinearForm L_mob(Vsp1, mobility);
  VariationalProblem pde_m(a_mob, L_mob);
  pde_m.solve(mobility_projected);

  // Project permeabilty for post-processing
  //ScalarP1Projection::FunctionSpace Vsp1(mesh);
  //Function permeability_p1;
  //ScalarP1Projection::BilinearForm a_perm(Vsp1, Vsp1);
  //ScalarP1Projection::LinearForm L_perm(Vsp1, permeability);
  //VariationalProblem pde_perm(a_perm, L_perm);
  //pde_perm.solve(permeability_p1);

  // Output files
  File f0("results/velocity.pvd",     "compressed");
  File f1("results/presure.pvd",      "compressed");
  File f2("results/saturation.pvd",   "compressed");
  File f3("results/permeability.pvd", "compressed");
  File f4("results/mobility.pvd", "compressed");
  f2 << S0;
  //f3 << permeability_p1;

  // Allocate scratch arrays for updating terms
  std::vector<double> _S0(S0.vector().size());
  std::vector<double> _F(S0.vector().size());
  std::vector<double> _dFdS(S0.vector().size());

  // Start solution procedure
  tic();
  for(dolfin::uint i = 0; i < 1000; ++i)
  {
    cout << "Step, time " << i << "  " << t << endl;

    // Compute mobility 
    //mobility.vector().lambda( _1 = 1.0);
    mobility.vector().lambda(S0.vector(), (1.0/0.2)*_1*_1 + (1.0-_1)*(1.0-_1));

    // Solve Darcy flow problem for velocity and pressure
    pde_darcy.solve(U);

    // Compute fractional flow function and its derivative (derivative is used to compute dt)
    //F.vector().lambda(S0.vector(), _1);
    //dFdS.vector().lambda(_1 = 1.0);

    S0.vector().get_local(&_S0[0]);
    F.vector().get_local(&_F[0]);
    dFdS.vector().get_local(&_dFdS[0]);
    for(dolfin::uint i = 0; i < F.vector().size(); ++i)
    {
      // _F[i] = _S0[i]*_S0[i]/( _S0[i]*_S0[i] + 0.2*(1.0-_S0[i])*(1.0-_S0[i]) );

      double ftmp  = ( _S0[i]*_S0[i] + 0.2*(1.0-_S0[i])*(1.0-_S0[i]) );
      _F[i]    = _S0[i]*_S0[i]/ftmp;      
      if( _S0[i] > DOLFIN_EPS)
        _dFdS[i] = 2.0*_S0[i]/ftmp - _S0[i]*_S0[i]*(2.0*_S0[i] - 0.4*(1.0-_S0[i]))/(ftmp*ftmp); 
      else
        _dFdS[i]  = 1.0;        
    }
    F.vector().set_local(&_F[0]);  
    dFdS.vector().set_local(&_dFdS[0]);  

    // Compute dt and update timne
    u = U[0];
    Function unorm(Vp0);
    pde_dt.solve(unorm);
    if(std::abs(unorm.vector().max()) < DOLFIN_EPS)
      error("Cannot compute time step (dividing by zero)");
    dt = 0.25*0.9/(4.0*sqrt(unorm.vector().max()));
    t += dt;
    cout << "Computed dt " << dt << endl;

    // Solve saturation equation
    pde_saturation.solve(S);

    // Update saturation
    S0 = S;

    // Check saturation bounds and bring back to [0, 1]
    S0.vector().lambda(if_then_else(_1 < 0.0, _1 = 0.0, if_then(_1 > 1.0, _1 = 1.0))); 

    // Project velocity for post-processing
    Function u_p1(Vp1);
    pde_u_p1.solve(u_p1);

    //pde_s.solve(saturation_projected);
    //pde_m.solve(mobility_projected);

    // Output to file
    f0 << u_p1;
    f1 << U[1];
    f2 << S0;
    //f2 << saturation_projected;
    //f4 << mobility_projected;
  }
  cout << "Finished. Solve time " << toc() << endl;
}
