from dolfin import *

mesh = UnitSquare(3, 3)
n = FacetNormal(mesh)

PV = VectorFunctionSpace(mesh, "CG", 1)
P0 = FunctionSpace(mesh, "Discontinuous Lagrange", 0)
mixed_space = MixedFunctionSpace([PV, P0])

# Function spaces and functions
r  = TestFunction(P0)
U0 = Function(mixed_space)
u0, s0 = split(U0)

U0.vector()[:] = 1.0

F    = 1.0/(s0**2 + 1.0)

# This breaks FFC
un   = (dot(u0, n) + sqrt(dot(u0, n)*dot(u0, n)))/2.0
un_h = (dot(u0, n) - sqrt(dot(u0, n)*dot(u0, n)))/2.0

# This works
#un   = 0.5*(dot(u0, n) + sqrt(dot(u0, n)*dot(u0, n)))
#un_h = 0.5*(dot(u0, n) - sqrt(dot(u0, n)*dot(u0, n)))

L = jump(r)*(un('+')*F('+') - un('-')*F('-'))*dS 

param0 = {"representation": "quadrature", "optimize": False}
param1 = {"representation": "quadrature", "optimize": True}
print "Vec norm (no opt):   ", assemble(L,
                             form_compiler_parameters=param0).norm("l2")
print "Vec norm (with opt): ", assemble(L,
                             form_compiler_parameters=param1).norm("l2")

