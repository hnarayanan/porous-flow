"""
This program solves pressure-driven, time-dependent flow of two phases
through porous media.

Strong form:

    (lambda(s)*K)^(-1)*u + grad(p) = 0
                            div(u) = 0
              ds/dt + u.grad(F(s)) = 0,

where,

    lambda(s) = 1.0/mu_rel*s^2 + (1 - s)^2
         F(s) = k_rw(s)/mu_w/(k_rw(s)/mu_w + k_ro(s)/mu_o)
              = s^2/(s^2 + mu_rel*(1 - s)^2).

One can then can post-calculate the velocity of each phase using the
relation: u_j = - (k_rj(s)/mu_j)*K*grad(p).

Weak form:

Find u, p, s in V such that,

   (v, (lambda*K)^(-1)*u) - (div(v), p) = - (v, pbar*n)_N       (1)
                            (q, div(u)) = 0                     (2)
            (r, ds/dt) - (grad(r), F*u) = - (r, F*u.n)_N        (3)

for all v, q, r in V'.

Adjoint problem:

Find z_u, z_p, z_s in V' such that,

(z_u, (lambda*K)^(-1)*v) - (div(z_u), q) + (z_p, div(v)) + (z_s, r)
- (z_s, s0) - dt*(grad(z_s), F*v) + dt*(z_s, F*v.n)_N = M(v)

for all v, q, r in V - V_h.

Error estimate (per time step):

Using the dual weighted residual method (ignoring jump terms and and
constant multipliers),

  M(u) - M(u_h) ~=  \sum_T |<R_T, z - z_h>_T + <R_dT, z - z_h>_dT|

Model problem:

 -----4-----
 |         |
 1         2
 |         |
 -----3-----

Initial conditions:
u(x, 0) = 0
p(x, 0) = 0
s(x, 0) = 0 in \Omega

Boundary conditions:
p(x, t) = 1 - x on \Gamma_{1, 2, 3, 4}
s(x, t) = 1 on \Gamma_1 if u.n < 0
s(x, t) = 0 on \Gamma_{2, 3, 4} if u.n > 0

Goal functionals:
M(v) = inner(grad(u_h), grad(v))*dx
M(v) = inner(u_h, v)*dx
M(v) = inner(v, n)*ds(2)

Parameters:
mu_rel, Kinv, lmbdainv, F, dt, T

This implementation includes functional forms from the deal.II demo
available at: http://www.dealii.org/6.2.1/doxygen/deal.II/step_21.html
"""

__author__    = "Harish Narayanan (harish@simula.no) and Garth N. Wells (gnw20@cam.ac.uk)"
__date__      = "2010-01-30"
__copyright__ = "Copyright (C) 2010 Harish Narayanan and Garth N. Wells"
__license__   = "GNU GPL Version 3.0"

from dolfin import *
from numpy import array, sort, zeros, max, abs
import sys
import os

# This program does not run in parallel
not_working_in_parallel("This program")

# This program does not work without CGAL
if not has_cgal():
    print "DOLFIN must be compiled with CGAL to run this program."
    exit(0)

# Optimise compilation of forms
#parameters["form_compiler"]["log_level"] = INFO
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True

# Parameters related to the adaptivity
TOL = 0.005          # Desired error tolerance
REFINE_RATIO = 0.05  # Fraction of cells to refine in each iteration
MAX_ITER = 10        # Maximum number of iterations
MIN_SIZE = 0.015625  # Minimum cell diameter

# Physical parameters, functional forms and boundary conditions
# Relative viscosity of water w.r.t. crude oil
mu_rel = 0.2

# Spatially-varying permeability matrix (inverse)
kinv = Expression("1.0/std::max(exp(-pow((x[1] - 0.5 - 0.1*sin(10*x[0]))/0.1, 2.0)), 0.01)")
zero = Expression("0.0")
Kinv = as_matrix(((kinv, zero), (zero, kinv)))

# Total mobility
def lmbdainv(s):
    return 1.0/((1.0/mu_rel)*s**2 + (1.0 - s)**2)

# Fractional flow function
def F(s):
    return s**2/(s**2 + mu_rel*(1.0 - s)**2)

# Pressure boundary condition
class PressureBC(Expression):
    def eval(self, values, x):
        values[0] = 1.0 - x[0]

# Saturation boundary condition
class SaturationBC(Expression):
    def eval(self, values, x):
        if x[0] < DOLFIN_EPS:
            values[0] =  1.0

class TwoPhaseFlow(NonlinearProblem):
    def __init__(self, a, L, ffc_parameters):
        NonlinearProblem.__init__(self)
        self.L = L
        self.a = a
        self.reset_sparsity = True
        self.ffc_parameters = ffc_parameters
    def F(self, b, x):
        assemble(self.L, tensor=b, form_compiler_parameters=self.ffc_parameters)
    def J(self, A, x):
        assemble(self.a, tensor=A, reset_sparsity=self.reset_sparsity,
                 form_compiler_parameters=self.ffc_parameters)
        self.reset_sparsity = False

# Print in colour
def print_good(output):
    print '\033[1;32m%s\033[1;m' % output

def print_better(output):
    print '\033[1;42m%s\033[1;m' % output

# Order of the elements
order = 1

# Parameters to the form compiler
#ffc_parameters = {"quadrature_degree": order + 1, "representation": "quadrature"}
ffc_parameters = {"representation": "quadrature"}

if not os.path.exists("./results"):
    os.makedirs("./results")

u_file = File("./results/velocity.pvd")
p_file = File("./results/pressure.pvd")
s_file = File("./results/saturation.pvd")

# Computational domain and geometry information
mesh0 = UnitSquare(8, 8, "crossed")
mesh_new = Mesh(mesh0)

# Function spaces and functions on the intial mesh
BDM_init = FunctionSpace(mesh0, "Brezzi-Douglas-Marini", order)
DG_init  = FunctionSpace(mesh0, "Discontinuous Lagrange", order - 1)
ME       = MixedFunctionSpace([BDM_init, DG_init, DG_init])
U0       = Function(ME)

# Time step size and number of steps
dt = Constant(0.005)
N = 1500
t = 0.0
T = N*float(dt)

while t < T:

    t += float(dt)
    print_good("Solving at time t = %f" % t)

    # Start the adaptive algorithm
    for level in xrange(MAX_ITER):

        # Update to new mesh
        mesh = mesh_new
        n = FacetNormal(mesh)

        # Function spaces
        BDM = FunctionSpace(mesh, "Brezzi-Douglas-Marini", order)
        DG  = FunctionSpace(mesh, "Discontinuous Lagrange", order - 1)
        ME  = MixedFunctionSpace([BDM, DG, DG])

        # Spaces to project to
        P0s = FunctionSpace(mesh, "Discontinuous Lagrange", 0)
        P0v = VectorFunctionSpace(mesh, "Discontinuous Lagrange", 0)

        # Functions
        V   = TestFunction(ME)
        dU  = TrialFunction(ME)
        U   = Function(ME)

        v, q, r = split(V)
        u, p, s = split(U)
        u0, p0, s0 = split(U0)

        s_mid = 0.5*(s0 + s)

        pbar = PressureBC(degree=1)
        sbar = SaturationBC(degree=1)

        # Variational forms for the primal problem
        L1 = inner(v, lmbdainv(s_mid)*Kinv*u)*dx - div(v)*p*dx \
            + inner(v, pbar*n)*ds
        L2 = q*div(u)*dx

        # Upwind normal velocity: (inner(v, n) + |inner(v, n)|)/2.0
        # (using velocity from previous step on facets)
        un   = 0.5*(dot(u0, n) + sqrt(dot(u0, n)*dot(u0, n)))
        un_h = 0.5*(dot(u0, n) - sqrt(dot(u0, n)*dot(u0, n)))
        L_facet = dt('+')*dot(jump(r), un('+')*F(s_mid)('+') - un('-')*F(s_mid)('-'))*dS \
                + dt*r*un_h*sbar*ds
        L3 = r*(s - s0)*dx - dt*dot(grad(r), F(s_mid)*u)*dx \
           + dt*r*F(s_mid)*un*ds + L_facet

        # Total L
        L = L1 + L2 + L3

        # Jacobian
        a = derivative(L, U, dU)

        # Setup and solve the primal problem
        problem = TwoPhaseFlow(a, L, ffc_parameters)
        solver  = NewtonSolver()
        solver.parameters["absolute_tolerance"] = 1e-14
        solver.parameters["relative_tolerance"] = 1e-9
        solver.parameters["maximum_iterations"] = 10

        print "Solving primal problem"
        solver.solve(problem, U.vector())
        u, p, s = U.split()

        # Variational forms for the adjoint problem
        a_adjoint = adjoint(a)
        L_adjoint = inner(grad(u), grad(v))*dx

        print "Create and solve adjoint problem"
        problem_adjoint = VariationalProblem(a_adjoint, L_adjoint)
        (z_u, z_p, z_s) = problem_adjoint.solve().split()

        R1 = project(lmbdainv(s)*Kinv*u + grad(p), P0v)
        R2 = project(div(u), P0s)
        R3 = project((s - s0)/float(dt) + inner(u, grad(F(s))), P0s)

        # Compute the derivatives of the solutions of the adjoint problem
        Dz_u = project(div(z_u), P0s)

        # FIXME: These are fake derivatives
        Dz_p = project(Constant((1.0, 1.0)), P0v)
        Dz_s = project(Constant((1.0, 1.0)), P0v)

        # Estimate the error
        E1 = zeros(mesh.num_cells()) # From ||Dz_u|| ||h R1||
        E2 = zeros(mesh.num_cells()) # From ||Dz_s|| ||h R2||
        E3 = zeros(mesh.num_cells()) # From ||Dz_p|| ||h R3||
        E = zeros(mesh.num_cells())  # Total

        # FIXME: The following can be improved by evaluation at cells,
        # rather than points. This will be cleaned up for efficiency and
        # pythonic style after the error estimators begin to make sense.

        i = 0
        for c in cells(mesh):
            h = c.diameter()
            K = c.volume()
            x = array((c.midpoint().x(), c.midpoint().y()))
            E1[i] = abs(Dz_u(x))*h*sqrt(R1(x)[0]**2 + R1(x)[1]**2)*sqrt(K)
            E2[i] = sqrt(Dz_p(x)[0]**2 + Dz_p(x)[1]**2)*h*abs(R2(x))*sqrt(K)
            E3[i] = sqrt(Dz_s(x)[0]**2 + Dz_s(x)[1]**2)*h*abs(R3(x))*sqrt(K)
            E[i] = E1[i] + E2[i] + E3[i]
            i = i + 1

        E1_norm = sqrt(sum([e1*e1 for e1 in E1]))
        E2_norm = sqrt(sum([e2*e2 for e2 in E2]))
        E3_norm = sqrt(sum([e3*e3 for e3 in E3]))
        E_norm  = sqrt(sum([e*e for e in E]))

        print_good("Level %d: E = %g (TOL = %g)" % (level, E_norm, TOL))

        # Check convergence
        if E_norm < TOL:
            print_better("Success, solution converged after %d iterations" % level)
            break

        # Mark cells for refinement
        cell_markers = MeshFunction("bool", mesh, mesh.topology().dim())
        E_0 = sorted(E, reverse=True)[int(len(E)*REFINE_RATIO)]
        for c in cells(mesh):
            cell_markers[c] = E[c.index()] > E_0 and c.diameter() > MIN_SIZE

        # Refine
        print "Refining mesh"
        mesh_new = refine(mesh, cell_markers)
        print "Finished refining mesh"

        print "----------------------------------"
        print memory_usage()
        print "----------------------------------"


    # Plot and store interesting quantities
    # plot(mesh1, title="Mesh at time t = %f" % t)
    # plot(u, title="Velocity")
    # plot(p, title="Pressure")
    # plot(s, title="Saturation at time t = %f" % t)
    u_file << u
    p_file << p
    s_file << s

    # Note: does this affect conservation properties???
    # Interpolate solution to new mesh
    BDM0 = FunctionSpace(mesh, "Brezzi-Douglas-Marini", order)
    DG0 = FunctionSpace(mesh, "Discontinuous Lagrange", order - 1)
    mixed_space0 = MixedFunctionSpace([BDM0, DG0, DG0])
    U0 = Function(mixed_space0)
    print "Interpolate U0"
    U0.interpolate(U)
    print "Finished interpolate U0"

