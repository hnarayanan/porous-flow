from dolfin import *
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True

# Computational domain and geometry information
mesh_init = UnitSquare(2, 2, "crossed")
mesh_new = Mesh(mesh_init)

V0 = FunctionSpace(mesh_init, "CG", 1)
u0 = Function(V0)

t  = 0.0
dt = Constant(0.005)
T  = 100*float(dt) 
while t < T:

    t += float(dt)

    for level in xrange(5):

        mesh = mesh_new
        V   = FunctionSpace(mesh, "CG", 1)
        v   = TestFunction(V)
        du  = TrialFunction(V)
        u   = Function(V)

        print "Verify function dim (a)"
        print u.function_space().mesh().num_cells()
        print "Verify function dim (b)"
        print u0.function_space().mesh().num_cells()
        print "End verify functions dims"
  
        L = v*(u-u0)*dx
        a = derivative(L, u, du)

        #pde = VariationalProblem(a, L)
        #pde.solve()

    mesh_new = Mesh(mesh)
    mesh_new.refine()

    u0 = Function(V)
    u0.interpolate(u)

    raw_input("Check memory use and press ENTER to continue")

