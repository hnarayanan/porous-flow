from dolfin import *
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True

# Computational domain and geometry information
mesh_init = UnitSquare(2, 2, "crossed")
mesh_new = Mesh(mesh_init)

V0 = FunctionSpace(mesh_init, "CG", 1)
ME0 = V0 + V0
U0 = Function(ME0)
#u0 = Function(V0)
#p0 = Function(V0)
#s0 = Function(V0)

#BDM0 = FunctionSpace(mesh_init, "Brezzi-Douglas-Marini", 1)
#DG0 = FunctionSpace(mesh_init, "Discontinuous Lagrange", 0)
#ME0 = MixedFunctionSpace([BDM0, DG0, DG0])
#U0 = Function(ME0)

t  = 0.0
dt = Constant(0.005)
T  = 100*float(dt) 
while t < T:

    t += float(dt)

    for level in xrange(2):

        mesh = mesh_new

        #BDM = FunctionSpace(mesh, "Brezzi-Douglas-Marini", 1)
        #DG = FunctionSpace(mesh, "Discontinuous Lagrange", 0)
        #ME = MixedFunctionSpace([BDM, DG, DG])

        V  = FunctionSpace(mesh, "CG", 1)
        ME = V + V
        ME = MixedFunctionSpace([V, V])

        (v, q) = TestFunctions(ME)
        dU     = TrialFunction(ME)
        U      = Function(ME)

        #(v, q, r) = TestFunctions(ME)
        #dU        = TrialFunction(ME)
        #U         = Function(ME)

        #u0,  p0, s0  = split(U0)
        #u,   p, s    = split(U)
        #du, dp, ds   = split(dU)

        u0,  p0 = split(U0)
        u,   p  = split(U)
        du, dp  = split(dU)

        print "Verify function dim (a)"
        print U.function_space().mesh().num_cells()
        print "Verify function dim (b)"
        print U0.function_space().mesh().num_cells()
        print "End verify functions dims"
  
        L = dot(v, u-u0)*dx + q*(p-p0)*dx #+ r*(s-s0)*dx
        a = derivative(L, U, dU)

        pde = VariationalProblem(a, L)
        U = pde.solve()

        mesh_new = Mesh(mesh)
        mesh_new.refine()

    #U0 = U
    V1 = FunctionSpace(mesh_new, "CG", 1)
    ME1 = V1 + V1
    U0 = Function(ME1)
    #U0 = Function(ME)
    U0.interpolate(U)

    raw_input("Check memory use and press ENTER to continue")

