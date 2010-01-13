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

# Saturation boundary condition
class Gss(Expression):
    def eval(self, values, x):
        if x[0] < DOLFIN_EPS:
            values[0] =  0.1

# Saturation boundary condition
class Velocity(Expression):
    def eval(self, values, x):
        values[0] =  1.0
        values[1] =  0.0

# Define boundary of domain
class Boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

# Create mesh
#mesh = UnitSquare(36, 36, "crossed")

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

mixed_space = MixedFunctionSpace([BDM, P0, P0, P0])
#mixed_space = MixedFunctionSpace([BDM, P0])

# Boundary function
gsw = Gsw(degree=1)
gss = Gss(degree=1)
gp = Gp(degree=1)
#u = Velocity(V = Pv1)

dt = 0.01
f = 0.0
k = 1.0
mobility = 1.0

V   = TestFunction(mixed_space)
dU  = TrialFunction(mixed_space)
U   = Function(mixed_space)
U0  = Function(mixed_space)

w, q, v0, v1       = split(V)
u, p, s0, s1       = split(U)
u0, p0, s0_0, s1_0 = split(U0)

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
un   = (dot(u0, n) + sqrt(dot(u0, n)*dot(u0, n)))/2.0
un_h = (dot(u0, n) - sqrt(dot(u0, n)*dot(u0, n)))/2.0

### Saturation equations

## Phase 0
# Fractional flow function
s0_mid = 0.5*(s0 + s0_0)
F0 = s0_mid*s0_mid/( s0_mid*s0_mid + 0.2*(1.0-s0_mid)*(1.0-s0_mid) );
L_s0 =  dot(v0, s0)*dx - dot(v0, s0_0)*dx - dt*dot(grad(v0), u*F0)*dx \
      + dt*dot(jump(v0), un('+')*F0('+') - un('-')*F0('-'))*dS \
      + dt*dot(v0, un*F0)*ds + dt*dot(v0, un_h*gsw)*ds
a_s0 = derivative(L_s0, U, dU)

## Phase 1
# Fractional flow function
s1_mid = 0.5*(s1 + s1_0)
F1 = s1_mid*s1_mid/( s1_mid*s1_mid + 0.2*(1.0-s1_mid)*(1.0-s1_mid) );
L_s1 =  dot(v1, s1)*dx - dot(v1, s1_0)*dx - dt*dot(grad(v1), u*F1)*dx \
      + dt*dot(jump(v1), un('+')*F0('+') - un('-')*F1('-'))*dS \
      + dt*dot(v1, un*F1)*ds + dt*dot(v1, un_h*gss)*ds
a_s1 = derivative(L_s1, U, dU)

L_s = L_s0 + L_s1
a_s = a_s0 + a_s1

# Sum forms
L = L_darcy + L_s
a = a_darcy + a_s

pde= VariationalProblem(a, L, nonlinear=True)

fs0 = File("results/s0_phase.pvd")
fs1 = File("results/s1_phase.pvd")
fu = File("results/u.pvd")

t = 0.0
T = 50*dt
while t < T:
    t += dt
    U0.assign(U)
    pde.solve(U)
    u, p, s0, s1 = U.split() 
    fs0 << s0
    fs1 << s1
    uh = project(u, Pv1)
    fu << uh

uh = project(u, Pv1)
plot(uh)
interactive() 

