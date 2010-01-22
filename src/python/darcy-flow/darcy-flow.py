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
k = "1.0/(exp(-(((x[1] - 0.5 - 0.1*sin(10*x[0]))/0.1)*((x[1] - 0.5 - 0.1*sin(10*x[0]))/0.1))) + 1.0)"
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

# Define the bilinear form
def a(v, q, u, p):
    return dot(Kinv*v, u)*dx - div(v)*p*dx + q*div(u)*dx

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
a_primal = a(v, q, u, p)
L_primal = q*f*dx - inner(v, pbar*n)*ds

# Compute primal solution
print "Solve primal problem"
problem_primal = VariationalProblem(a_primal, L_primal)
(u_h, p_h) = problem_primal.solve().split()

# Post-process and plot interesting fields
print "Projecting velocity to a continuous space"
P1 = VectorFunctionSpace(mesh, "CG", 1)
u_h_proj = project(u_h, P1)

plot(kinv11, title="Inverse permeability magnitude", mesh=mesh)
plot(u_h_proj, title="Velocity")
plot(p_h, title="Pressure")
interactive()

# Save solutions in pvd format
u_file = File("primal_velocity.pvd")
p_file = File("primal_pressure.pvd")

u_file << u_h_proj
p_file << p_h
