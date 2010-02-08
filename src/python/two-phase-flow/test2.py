from dolfin import *

# Computational domain and geometry information
mesh_init = UnitSquare(2, 2, "crossed")
mesh_new = Mesh(mesh_init)
V0 = FunctionSpace(mesh_init, "CG", 1)
U0 = Function(V0)

# Test refinement for single-field element
for t in xrange(3):
    for level in xrange(2):
        mesh = mesh_new
        V  = FunctionSpace(mesh, "CG", 1)
        v  = TestFunction(V)
        dU = TrialFunction(V)
        U  = Function(V)

        L = v*(U-U0)*dx
        a = v*dU*dx
        pde = VariationalProblem(a, L)
        U = pde.solve()

        # Refine mesh
        mesh_new = Mesh(mesh)
        mesh_new.refine()

    # Update U0
    U0 = U

print "\n \n ****Testing mixed element \n"

# Test refinement for mixed element
mesh_init = UnitSquare(2, 2, "crossed")
mesh_new = Mesh(mesh_init)
V0 = FunctionSpace(mesh_init, "CG", 1)
ME0 = V0 + V0
U0 = Function(ME0)
for t in xrange(3):
    for level in xrange(2):

        mesh = mesh_new
        V  = FunctionSpace(mesh, "CG", 1)
        ME = MixedFunctionSpace([V, V])
        (v, q) = TestFunctions(ME)
        dU     = TrialFunction(ME)
        U      = Function(ME)

        u0,  p0 = split(U0)
        u,   p  = split(U)
        du, dp  = split(dU)

        print "Verify function dim (a)"
        print U.function_space().mesh().num_cells()
        print "Verify function dim (b)"
        print U0.function_space().mesh().num_cells()
        print "End verify functions dims"
  
        L = dot(v, u-u0)*dx + q*(p-p0)*dx
        a = v*du*dx + q*dp*dx

        pde = VariationalProblem(a, L)
        U = pde.solve()

        # Refine mesh
        mesh_new = Mesh(mesh)
        mesh_new.refine()

    # This breaks for mixed element
    #U0 = U

    # This works for mixed element
    V1 = FunctionSpace(mesh_new, "CG", 1)
    ME1 = V1 + V1
    U0 = Function(ME1)
    U0.interpolate(U)

print "Finsihed successfully"
