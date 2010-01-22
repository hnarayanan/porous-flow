"""
This program solves pressure-driven, steady-state Darcy's flow on a
square plate with spatially varying permeability.

        Kinv*u + grad(p) = 0
                  div(u) = 0

Which, in weak form reads:

 (v, Kinv*u) - (div(v), p) = - (v, p*n)_N  for all v
               (q, div(u)) = 0             for all q

It then proceeds to construct and solve the adjoint problem, with a
suitable goal functional of interest as its right hand-side.

 (w, Kinv*v) - (div(w), q) + (r, div(v)) = M(v) for all v, q

We solve this adjoint problem for the variable z = (w, r). After this,
the error is estimated (rather crudely for now: ignoring jump terms
and constant multipliers) using the dual weighted residual method:

  M(u) - M(u_h) ~=  \sum_T |<R_T, z - z_h>_T + <R_dT, z - z_h>_dT|

The error estimate is used to suitably refine the mesh, and the above
process is repeated until a certain tolerance is reached.
"""

__author__    = "Harish Narayanan (harish@simula.no)"
__date__      = "2010-01-20"
__copyright__ = "Copyright (C) 2009 Harish Narayanan"
__license__   = "GNU GPL Version 3.0"

from dolfin import *
from numpy import array, sort, zeros, max, abs

# This program does not run in parallel
not_working_in_parallel("This program")

# This program does not work without CGAL
if not has_cgal():
    print "DOLFIN must be compiled with CGAL to run this program."
    exit(0)

# Parameters related to the adaptivity
TOL = 1.e-15         # Desired error tolerance
REFINE_RATIO = 0.25  # Fraction of cells to refine in each iteration
MAX_ITER = 10        # Maximum number of iterations

# Parameters and boundary conditions related to the physics
# Spatially-varying permeability matrix (inverse)
k = "std::max(exp(-(((x[1] - 0.5 - 0.1*sin(10*x[0]))/0.1)*((x[1] - 0.5 - 0.1*sin(10*x[0]))/0.1))), 0.01) + 1.0"
# k = "cos(4*pi*x[1]*x[0])/5.0 + 1.0"
kinv11 = Expression(k)
kinv12 = Constant(0.0)
kinv21 = Constant(0.0)
kinv22 = Expression(k)
Kinv = as_matrix(((kinv11, kinv12), (kinv21, kinv22)))

# Pressure boundary condition
class PressureBC(Expression):
    def eval(self, values, x):
        values[0] = 1.0 - x[0]

# Create initial mesh
mesh = UnitSquare(8, 8)
n = FacetNormal(mesh)

# Plot the permeability inverse to make some sense of subsequent
# solutions
# plot(kinv11, title="Inverse permeability magnitude", mesh=mesh)

# Start the adaptive algorithm
for level in xrange(MAX_ITER):

    # FIXME: Continuous function spaces used in order to drive
    # residual to 0. The other BDM + DG set works as well, just that
    # the indicators won't be driven to very small values.
    BDM = FunctionSpace(mesh, "BDM", 2)
    DG = FunctionSpace(mesh, "DG", 1)
    V = BDM + DG

    # P2 = VectorFunctionSpace(mesh, "CG", 2)
    # P1 = FunctionSpace(mesh, "CG", 1)
    # V  = P2 + P1

    # Define the required functions
    (u, p) = TrialFunctions(V)
    (v, q) = TestFunctions(V)
    pbar = PressureBC()

    # Pose primal problem
    a_primal = dot(v, Kinv*u)*dx - div(v)*p*dx + q*div(u)*dx
    L_primal = - inner(v, pbar*n)*ds

    # Compute primal solution
    print "Solve primal problem"
    problem_primal = VariationalProblem(a_primal, L_primal)
    (u_h, p_h) = problem_primal.solve().split()

    # Construct continuous spaces to project to
    P1s = FunctionSpace(mesh, "CG", 1)
    P1v = VectorFunctionSpace(mesh, "CG", 1)

    # Project u_h for post-processing
    print "Projecting solution to P1 for post-processing"
    u_h_proj = project(u_h, P1v)

    # Plot the primal solution fields
    # plot(u_h_proj, title="Primal velocity")
    # plot(p_h, title="Primal pressure")

    # Pose adjoint problem with some arbitrary goal functional. Note
    # that playing with different forms of the goal results in
    # different interesting results
    a_adjoint = adjoint(a_primal)
    L_adjoint = inner(grad(u_h), grad(v))*dx

    print "Solve adjoint problem"
    problem_adjoint = VariationalProblem(a_adjoint, L_adjoint)
    (w_h, r_h) = problem_adjoint.solve().split()

    # Project w_h for post-processing
    print "Projecting solution to P1 for post-processing"
    w_h_proj = project(w_h, P1v)
    r_h_proj = project(r_h, P1s)

    # Plot the adjoint solution fields
    # plot(w_h_proj, title="Adjoint velocity")
    # plot(r_h, title="Adjoint pressure")

    # Construct discontinuous spaces to project to
    P0s = FunctionSpace(mesh, "DG", 0)
    P0v = VectorFunctionSpace(mesh, "DG", 0)

    # Calculate the residuals and project them obtain constant values
    # on cells. FIXME: Perform integration by parts by hand instead of
    # using project.
    print "Calculating and projecting residuals"
    P1 = VectorFunctionSpace(mesh, "CG", 1)
    R1 = project(Kinv*u_h + grad(p_h), P0v)
    R2 = project(div(u_h), P0s)

    # plot(R1, title="Residual of the Darcy flow equation")
    # plot(R2, title="Residual of the mass conservation equation")

    # Compute the derivatives of the solutions of the adjoint problem
    # FIXME: The following projects might do something strange
    Dw_h = project(div(w_h), P0s)
    Dr_h = project(grad(r_h), P0v)

    # plot(Dw_h, title="Divergence of the adjoint velocity")
    # plot(Dr_h, title="Gradient of the adjoint pressure")

    # Estimate the error
    E1 = zeros(mesh.num_cells()) # From ||Dw|| ||h R1||
    E2 = zeros(mesh.num_cells()) # From ||Dr|| ||h R2||
    E = zeros(mesh.num_cells())  # Total
    
    # FIXME: The following can be improved by evaluation at cells,
    # rather than points. This will be cleaned up for efficiency and
    # pythonic style after the error estimators begin to make sense.

    i = 0
    for c in cells(mesh):
        h = c.diameter()
        K = c.volume()
        x = array((c.midpoint().x(), c.midpoint().y()))
        E1[i] = abs(Dw_h(x))*h*sqrt(R1(x)[0]**2 + R1(x)[1]**2)*sqrt(K)
        E2[i] = sqrt(Dr_h(x)[0]**2 + Dr_h(x)[1]**2)*h*abs(R2(x))*sqrt(K)
        E[i] = E1[i] + E2[i]
        i = i + 1

    E1_norm = sqrt(sum([e1*e1 for e1 in E1]))
    E2_norm = sqrt(sum([e2*e2 for e2 in E2]))
    E_norm  = sqrt(sum([e*e for e in E]))
    
    print "Level %d: E = %g (TOL = %g)" % (level, E_norm, TOL)

    # Check convergence
    if E_norm < TOL:
        print "Success, solution converged after %d iterations" % level
        break

    # Mark cells for refinement
    cell_markers = MeshFunction("bool", mesh, mesh.topology().dim())
    E_0 = sorted(E, reverse=True)[int(len(E)*REFINE_RATIO)]
    for c in cells(mesh):
        cell_markers[c] = E[c.index()] > E_0

    # Refine mesh
    mesh.refine(cell_markers)

    # Plot mesh
    plot(mesh)

# Hold plot
interactive()
