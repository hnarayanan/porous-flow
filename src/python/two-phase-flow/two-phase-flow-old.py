""" . 
"""

__author__    = "Garth N. Wells (gnw20@cam.ac.uk)"
__date__      = "2009-10-28"
__copyright__ = "Copyright (C) 2009 Garth N. Wells"
__license__   = "GNU GPL Version 3.0"

from dolfin import *

# Optimize compilation of the form
parameters.optimize = False

# Sub domain for inflow (right)
class Inflow(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS and on_boundary

# Pressure boundary condition
class Gp(Expression):
    def eval(self, values, x):
        values[0] = 1.0 - x[0] 

# Saturation boundary condition
class Gsw(Expression):
    def eval(self, values, x):
        if x[0] < DOLFIN_EPS:
            values[0] =  0.3

# Velocity boundary condition
class Velocity(Expression):
    def eval(self, values, x):
        values[0] =  1.0
        values[1] =  0.0

# Define boundary of domain
class Boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

mesh = UnitSquare(16, 16)
n = FacetNormal(mesh)

# Create mesh function over the cell facets
sub_domains = MeshFunction("uint", mesh, mesh.topology().dim() - 1)
sub_domains.set_all(3)
inflow = Inflow()
inflow.mark(sub_domains, 1)

# Create function spaces
Pv1 = VectorFunctionSpace(mesh, "Lagrange", 1)
BDM = FunctionSpace(mesh, "BDM", 1)
P0  = FunctionSpace(mesh, "DG", 0)
P1  = FunctionSpace(mesh, "CG", 1)

mixed_space = MixedFunctionSpace([BDM, P0, P0])

# Boundary function
gsw = Gsw(degree=1)
gp = Gp(degree=1)

dt = 0.01
f = 0.0
k = 1.0
mobility = 1.0

V   = TestFunction(mixed_space)
dU  = TrialFunction(mixed_space)
U   = Function(mixed_space)
U0  = Function(mixed_space)

w, q, v    = split(V)
u, p, s    = split(U)
u0, p0, s0 = split(U0)

## Darcy equations

# Mass conservation eqn
L_darcy_mass = q*div(u)*dx - q*f*dx
a_darcy_mass = derivative(L_darcy_mass, U, dU)

# Pressure eqn
L_darcy_p = (1.0/k)*(1.0/mobility)*dot(w, u)*dx - div(w)*p*dx \
           + dot(w, n)*gp*ds
a_darcy_p = derivative(L_darcy_p, U, dU)

L_darcy = L_darcy_mass + L_darcy_p
a_darcy = a_darcy_mass + a_darcy_p

# Upwind normal velocity: ( dot(v, n) + |dot(v, n)| )/2.0 
# (using velocity from previous step on facets)
un   = (dot(u0, n) + sqrt(dot(u0, n)*dot(u0, n)))/2.0
un_h = (dot(u0, n) - sqrt(dot(u0, n)*dot(u0, n)))/2.0

### Saturation equation

# Mid-point value of saturation
s_mid = 0.5*(s + s0)

# Fractional flow function
F = s_mid*s_mid/( s_mid*s_mid + 0.2*(1.0-s_mid)*(1.0-s_mid) );

# Forms

# DG form
L_s =  dot(v, s)*dx - dot(v, s0)*dx - dt*dot(grad(v), u*F)*dx \
     + dt*dot(jump(v), un('+')*F('+') - un('-')*F('-'))*dS \
     + dt*dot(v, un*F)*ds + dt*dot(v, un_h*gsw)*ds

# CG form (crude with artificial diffusion)
# L_s =  dot(v, s)*dx - dot(v, s0)*dx - dt*dot(grad(v), u*F)*dx \
#       + dt*dot(grad(v), grad(s_mid))*dx \
#       + dt*dot(v, un*F)*ds + dt*dot(v, un_h*gsw)*ds

a_s = derivative(L_s, U, dU)

# Sum forms
L = L_darcy + L_s
a = a_darcy + a_s

pde= VariationalProblem(a, L, nonlinear=True)

fs = File("s_phase.pvd")
fu = File("u.pvd")

t = 0.0
T = 100*dt
while t < T:
    t += dt
    U0.assign(U)
    pde.solve(U)
    u, p, s = U.split() 
    fs << s
    uh = project(u, Pv1)
    fu << uh

uh = project(u, Pv1)
plot(uh)
interactive() 

