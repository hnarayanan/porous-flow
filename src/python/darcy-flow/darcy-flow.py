"""
This program solves pressure-driven, steady-state Darcy's flow on a
square plate with spatially varying permeability.

        Kinv*u + grad(p) = 0
                  div(u) = 0

Which, in weak form reads:

 (v, Kinv*u) - (div(v), p) = - (v, p*n)_N  for all v
               (q, div(u)) = 0             for all q
"""

__author__    = "Harish Narayanan (harish@simula.no)"
__date__      = "2010-01-15"
__copyright__ = "Copyright (C) 2009 Harish Narayanan"
__license__   = "GNU GPL Version 3.0"

from dolfin import *

# Construct a spatially-varying permeability matrix (inverse)
kinv11 = Expression("(cos(4*pi*x[1])/5.0 + 1.0)")
kinv12 = Constant(0.0)
kinv21 = Constant(0.0)
kinv22 = Expression("(cos(4*pi*x[1])/5.0 + 1.0)")
Kinv = as_matrix(((kinv11, kinv12), (kinv21, kinv22)))

# Pressure boundary condition
class PressureBC(Expression):
    def eval(self, values, x):
        values[0] = 1.0 - x[0]

# Define the bilinear form
#def a(v, q, u, p):
#    return dot(Kinv*v, u)*dx - div(v)*p*dx + q*div(u)*dx

# Create mesh and define function spaces
mesh = UnitSquare(32, 32)
n = FacetNormal(mesh)
BDM = FunctionSpace(mesh, "BDM", 1)
DG = FunctionSpace(mesh, "DG", 0)
V = BDM + DG

# Define variational problem
(u, p) = TrialFunctions(V)
(v, q) = TestFunctions(V)

pbar = PressureBC()
f = Constant(0.0)

# Pose primal problem
#a_primal = a(v, q, u, p)
a_primal = dot(Kinv*v, u)*dx - div(v)*p*dx + q*div(u)*dx
L_primal = q*f*dx - inner(v, pbar*n)*ds

# Compute primal solution
print "Solve primal problem"
problem_primal = VariationalProblem(a_primal, L_primal)
(U, P) = problem_primal.solve().split()

# Pose adjoint problem
#a_dual = a(u, p, v, q)
a_dual = adjoint(a_primal)
L_dual = v[0]*1.0*ds # Some goal: Average of one component of the velocity
                 #            over the boundary

print "Solve dual problem"
problem_dual = VariationalProblem(a_dual, L_dual)
(W, R) = problem_dual.solve().split()

# Project U, W for post-processing
print "Post-processing projections"
P1 = VectorFunctionSpace(mesh, "CG", 1)
U_proj = project(U, P1)
W_proj = project(W, P1)

# Plot interesting fields
#plot(kinv11, title="Inverse permeability magnitude", mesh=mesh)
#plot(U_proj, title="Velocity")
#plot(P, title="Pressure")

plot(W_proj, title="Dual velocity")
#plot(R, title="Dual pressure")
interactive()


# Save solutions in pvd format
u_file = File("primal_velocity.pvd")
p_file = File("primal_pressure.pvd")
w_file = File("dual_velocity.pvd")
r_file = File("dual_pressure.pvd")

u_file << U_proj
p_file << P

w_file << W_proj
r_file << R
