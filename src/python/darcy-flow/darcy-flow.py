"""
This program solves pressure-driven, steady-state Darcy's flow on a
square plate with spatially varying permeability.

        u + K*grad(p) = 0
               div(u) = 0

Which, in weak form reads:

 (v, u) - (div(K*v), p) = - (K*v, p*n)_N  for all v
            (q, div(u)) = 0               for all q
"""

__author__    = "Harish Narayanan (harish@simula.no)"
__date__      = "2010-01-15"
__copyright__ = "Copyright (C) 2009 Harish Narayanan"
__license__   = "GNU GPL Version 3.0"

from dolfin import *

# Pressure boundary condition
class PressureBC(Expression):
    def eval(self, values, x):
        values[0] = 1.0 - x[0]

# Create mesh and define function spaces
mesh = UnitSquare(32, 32)
n = FacetNormal(mesh)
BDM = FunctionSpace(mesh, "BDM", 1)
DG = FunctionSpace(mesh, "DG", 0)
V = BDM + DG

# Define variational problem
(u, p) = TrialFunctions(V)
(v, q) = TestFunctions(V)

# Construct a spatially-varying permeability matrix
# FIXME: The following doesn't work for spatially-varying expressions
k11 = Constant(1.0)
k12 = Constant(0.0)
k21 = Constant(0.0)
k22 = Constant(1.0)
f = Constant(0.0)

K = as_matrix(((k11, k12), (k21, k22)))
pbar = PressureBC()

a = dot(v, u)*dx - div(K*v)*p*dx + q*div(u)*dx
L = q*f*dx - inner(K*v, pbar*n)*ds

# Compute solution
problem = VariationalProblem(a, L)
(u, p) = problem.solve().split()

# Project u for post-processing
P1 = VectorFunctionSpace(mesh, "CG", 1)
u_proj = project(u, P1)

# Plot solution fields
plot(u_proj, title="Velocity", interactive=True)
plot(p, title="Pressure", interactive=True)

# Save solution to pvd format
u_file = File("u_darcy.pvd")
p_file = File("p_darcy.pvd")

u_file << u_proj
p_file << p
