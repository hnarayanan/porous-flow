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
from numpy import array, sort, zeros, max, abs

# Construct a spatially-varying permeability matrix (inverse)
kinv11 = Expression("(cos(4*pi*x[1])*cos(4*pi*x[0])/5.0 + 1.0)")
kinv12 = Constant(0.0)
kinv21 = Constant(0.0)
kinv22 = Expression("(cos(4*pi*x[1])*cos(4*pi*x[0])/5.0 + 1.0)")
Kinv = as_matrix(((kinv11, kinv12), (kinv21, kinv22)))

# The number of cells refined in each iteration
REFINE_RATIO = 0.10

# Pressure boundary condition
class PressureBC(Expression):
    def eval(self, values, x):
        values[0] = 1.0 - x[0]

# Define the bilinear form
def a(v, q, u, p):
    return dot(Kinv*v, u)*dx - div(v)*p*dx + q*div(u)*dx

# Create mesh and define function spaces
mesh = UnitSquare(64, 64)
n = FacetNormal(mesh)
BDM = FunctionSpace(mesh, "BDM", 1) # VectorFunctionSpace(mesh, "CG", 2)
DG = FunctionSpace(mesh, "DG", 0)   # FunctionSpace(mesh, "CG", 1)
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

# Pose adjoint problem
a_adjoint = adjoint(a_primal)
L_adjoint = v[0]*1.0*ds # Some goal: Average of one component of the velocity
                        #            over the boundary

print "Solve adjoint problem"
problem_adjoint = VariationalProblem(a_adjoint, L_adjoint)
(w_h, r_h) = problem_adjoint.solve().split()

# Estimate the error using the dual weighted residual method
#  M(u) - M(u_h) ~=  \sum_T |<R_T, z - z_h>_T + <R_dT, z - z_h>_dT|

# Calculate residuals
scalarDG = FunctionSpace(mesh, "DG", 1)
vectorDG = VectorFunctionSpace(mesh, "DG", 0)

R1 = project(Kinv*u_h + grad(p_h), vectorDG)
R2 = project(div(u_h), scalarDG)

# Calculate derivatives of dual fields
# FIXME: Perform integration by parts by hand instead of using project
Dw_h = project(div(w_h), scalarDG)
#Dr_h = project(grad(r_h), vectorDG)

# Determine error indicators
h = array([c.diameter() for c in cells(mesh)])
K = array([c.volume() for c in cells(mesh)])

E = zeros(mesh.num_cells())
E1 = zeros(mesh.num_cells())
E2 = zeros(mesh.num_cells())

# Variables to hold scalars and vectors
vectorval = array((0.0, 0.0))
scalarval = array((0.0))

# Calculate components of the error in each cell
i = 0
for c in cells(mesh):  
    x = array((c.midpoint().x(), c.midpoint().y()))
#    R1.eval(vectorval, x) FIXME: Needs CGAL
#    Dw_h.eval(scalarval, x)
    E1[i] = sqrt(vectorval[0]**2 + vectorval[1]**2)*abs(scalarval)
    E[i] = h[i]*E1[i]
    i = i + 1

# Calculate the error norm    
Enorm = 0
for i2 in range(mesh.num_cells()):
    Enorm  = Enorm + abs(E[i2])*h[i2]*h[i2]

cell_markers = MeshFunction("bool", mesh, mesh.topology().dim())
E_0 = sorted(E, reverse=True)[int(len(E)*REFINE_RATIO)]
for c in cells(mesh):
    cell_markers[c] = E[c.index()] > E_0

# plot(cell_markers, interactive=True)
	        
# mesh.refine(cell_markers) #FIXME: Needs CGAL

# Project u_h, w_h for post-processing
print "Post-processing projections"
P1 = VectorFunctionSpace(mesh, "CG", 1)
u_h_proj = project(u_h, P1)
w_h_proj = project(w_h, P1)

# Plots for debugging
# plot(div(u_h), title="div(u_h)")
# plot(div(u_h_proj), title="div(u_h_proj)")
# div_u_h_proj_dg = project(div(u_h), scalarDG)
# plot(div_u_h_proj_dg, title="div(u_h_proj_dg)")
# plot(Kinv*u_h + grad(project(p_h, FunctionSpace(mesh, "CG", 1))), title="R1")

# Plot interesting fields
plot(kinv11, title="Inverse permeability magnitude", mesh=mesh)
plot(u_h_proj, title="Velocity")
plot(p_h, title="Pressure")

plot(w_h_proj, title="Dual velocity")
plot(r_h, title="Dual pressure")
interactive()


# Save solutions in pvd format
u_file = File("primal_velocity.pvd")
p_file = File("primal_pressure.pvd")
w_file = File("dual_velocity.pvd")
r_file = File("dual_pressure.pvd")

u_file << u_h_proj
p_file << p_h

w_file << w_h_proj
r_file << r_h
