from dolfin import *

mesh = UnitSquare(32, 32, "crossed")

for i in xrange(10):
    for j in xrange(100):
        E = FunctionSpace(mesh, "CG", 1)


    raw_input("Check memory use and press ENTER to continue")
    print ". . . .", i

