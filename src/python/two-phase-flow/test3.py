from dolfin import *
from numpy import array, zeros, abs

kinv = Expression("1.0/std::max(exp(-pow((x[1] - 0.5 - 0.1*sin(10*x[0]))/0.1, 2.0)), 0.01)")
REFINE_RATIO = 0.24

mesh_0 = UnitSquare(8, 8, "crossed")
CG_0 = FunctionSpace(mesh_0, "Lagrange", 1)
U_0 = project(Expression("1 - x[0]*x[1]"), CG_0)

plot(U_0, title="0")

for time_step in range(3):

    mesh = Mesh(mesh_0)

    for refinement_iterator in range(3 + 1 - time_step):

        CG = FunctionSpace(mesh, "Lagrange", 1)
        CG6 = FunctionSpace(mesh, "Lagrange", 6)
        DG = VectorFunctionSpace(mesh, "Discontinuous Lagrange", 0)

        U_0 = interpolate(U_0, CG)

        u = TrialFunction(CG)
        v = TestFunction(CG)
        f = Constant(1.0)
        a = u*v*dx
        L = (U_0 + f)*v*dx
        problem = VariationalProblem(a, L)
        U = problem.solve()

        kinv = project(kinv, CG6)

        R = project(grad(kinv), DG)
        E = zeros(mesh.num_cells())

        i = 0
        for c in cells(mesh):
            h = c.diameter()
            K = c.volume()
            x = array((c.midpoint().x(), c.midpoint().y()))
            E[i] = abs(h*sqrt(R(x)[0]**2 + R(x)[1]**2)*sqrt(K))
            i += 1

        E_norm  = sqrt(sum([e*e for e in E]))

        cell_markers = MeshFunction("bool", mesh, mesh.topology().dim())
        E_0 = sorted(E, reverse=True)[int(len(E)*REFINE_RATIO)]
        for c in cells(mesh):
            cell_markers[c] = E[c.index()] > E_0
    
        mesh.refine(cell_markers)

    plot(U, title="%d" % (time_step + 1))
    U_0 = U


interactive()


