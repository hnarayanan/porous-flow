from dolfin import *
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True

# Computational domain and geometry information
mesh_init = UnitSquare(2, 2, "crossed")
mesh_new = Mesh(mesh_init)

V0 = FunctionSpace(mesh_init, "CG", 1)
ME0 = V0
U0 = Function(ME0)

# Mixed
#ME0 = V0 + V0
#U0 = Function(ME0)
#u0, p0 = split(U0)

t  = 0.0
dt = Constant(0.005)
T  = 100*float(dt) 
while t < T:

    t += float(dt)

    for level in xrange(2):

        mesh = mesh_new

        V  = FunctionSpace(mesh, "CG", 1)

        ME = V
        v  = TestFunction(ME)
        dU = TrialFunction(ME)
        U  = Function(ME)

        #ME = MixedFunctionSpace([V, V])
        #(v, q) = TestFunctions(ME)
        #dU     = TrialFunction(ME)
        #U      = Function(ME)

        #u0,  p0 = split(U0)
        #u,   p  = split(U)
        #du, dp  = split(dU)

        print "Verify function dim (a)"
        print U.function_space().mesh().num_cells()
        print "Verify function dim (b)"
        print U0.function_space().mesh().num_cells()
        print "End verify functions dims"
  
        #L = dot(v, u-u0)*dx + q*(p-p0)*dx
        L = v*(U-U0)*dx
        a = derivative(L, U, dU)

        pde = VariationalProblem(a, L)
        U = pde.solve()

        mesh_new = Mesh(mesh)
        mesh_new.refine()

    # This breaks for mixed element
    U0 = U

    # This works for mixed element
    #V1 = FunctionSpace(mesh_new, "CG", 1)
    #ME1 = V1 + V1
    #U0 = Function(ME1)
    #U0.interpolate(U)

    #U0 = Function(ME)
    #U0.interpolate(U)

    raw_input("Check memory use and press ENTER to continue")

