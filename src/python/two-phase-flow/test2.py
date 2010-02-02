from dolfin import *
from numpy import array, zeros, abs

kinv = Expression("1.0/std::max(exp(-pow((x[1] - 0.5 - 0.1*sin(10*x[0]))/0.1, 2.0)), 0.01)")
REFINE_RATIO = 0.24

mesh = UnitSquare(8, 8)

for iter in range(7):

    CG = FunctionSpace(mesh, "Lagrange", 6)
    DG = VectorFunctionSpace(mesh, "Discontinuous Lagrange", 0)

    kinv = project(kinv, CG)

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
    plot(mesh)

file = File("mesh.xml")
file << mesh

interactive()


