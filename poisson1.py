
from dolfin import *
import math

error_prev = 0 
for N in [10, 20, 40, 80, 160]: 
    mesh = UnitIntervalMesh(N)
    V = FunctionSpace(mesh, "Lagrange", 3) 
    u = TrialFunction(V)
    v = TestFunction(V)

    def boundary(x, on_boundary): return on_boundary
    bc = DirichletBC(V, 0, boundary)

    f = Expression("M_PI*M_PI*sin(M_PI*x[0])", element = V.ufl_element())
    a = -div(grad(u))*v*dx 
    a = inner(grad(u),grad(v))*dx 
    L = f*v*dx 

    U = Function(V)

    solve(a == L, U, bc) 

    U_exact = Expression("sin(M_PI*x[0])", degree=10, PI = math.pi)
    error = errornorm(U_exact, U, "H1")
    print (" N ", N, " error ", error)
    if error_prev > 0:  
        print ("rate ", math.log(error_prev / error, 2))
    error_prev = error 


