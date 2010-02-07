from dolfin import *
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True

def lmbdainv(s): return 1.0/((1.0/0.2)*s**2 + (1.0 - s)**2)
def F(s): return s**2/(s**2 + 0.2*(1.0 - s)**2)

mesh = UnitSquare(8, 8, "crossed")
n = FacetNormal(mesh)

t  = 0.0
dt = Constant(0.005)
T  = 100*float(dt) 
while t < T:
    t += float(dt)
    for level in xrange(5):
        BDM = FunctionSpace(mesh, "Brezzi-Douglas-Marini", 1)
        DG  = FunctionSpace(mesh, "Discontinuous Lagrange", 0)
        ME  = MixedFunctionSpace([BDM, DG, DG])

        V   = TestFunction(ME)
        dU  = TrialFunction(ME)
        U   = Function(ME)
        U0  = Function(ME)

        v, q, r = split(V)
        u, p, s = split(U)
        u0, p0, s0 = split(U0)
        s_mid = 0.5*(s0 + s)

        un   = 0.5*(dot(u0, n) + sqrt(dot(u0, n)*dot(u0, n)))
        un_h = 0.5*(dot(u0, n) - sqrt(dot(u0, n)*dot(u0, n)))
        L1 = inner(v, lmbdainv(s_mid)*u)*dx - div(v)*p*dx
        L2 = q*div(u)*dx
        L3 = r*(s - s0)*dx - dt*dot(grad(r), F(s_mid)*u)*dx \
            + dt*r*F(s_mid)*un*ds \
            + dt('+')*dot(jump(r), un('+')*F(s_mid)('+') \
            - un('-')*F(s_mid)('-'))*dS

        L = L1 + L2 + L3
        a = derivative(L, U, dU)

        pde = VariationalProblem(a, L)
        pde.solve()
    
    raw_input("Check memory use and press ENTER to continue")

