
from dolfin import *

import math

prev_L2_error = 0 
prev_H1_error = 0 
for N in [8,16,32,64]:

  mesh = UnitSquareMesh(N, N)
  V = FunctionSpace(mesh, "CG", 1)
  V3 = FunctionSpace(mesh, "CG", 3)

  u = TrialFunction(V)
  v = TestFunction(V)

  mu_value = 1.0e-2 
  mu_value = 1.0 
  mu = Constant(mu_value)
  f = Constant(0)
  h = mesh.hmin()

  a = (-u.dx(0)*v + mu*inner(grad(u),grad(v)))*dx  
  L = f*v*dx  

  u_analytical = Expression("(exp(-x[0]/%e) - 1)/ (exp(-1/%e) - 1)" % (mu_value, mu_value), degree=3 )  
  def boundary(x, on_boundary):
      return on_boundary and (x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS)

  bc = DirichletBC(V, u_analytical, boundary) 

  U = Function(V)
  solve(a == L, U, bc) 

  U_analytical = project(u_analytical, V3)

  L2_error = errornorm(U, u_analytical, mesh = V.mesh()) 
  print (" N = %d, L2 error = %e " % (N, L2_error))
  if prev_L2_error > 0: 
    print ("rate of L2 error ", math.log(prev_L2_error / L2_error, 2)) 
  prev_L2_error = L2_error 

  H1_error = errornorm(U, u_analytical, "H1", mesh = V.mesh()) 
  print (" N = %d, H1 error = %e " %  (N, H1_error))
  if prev_H1_error > 0: 
    print ("rate of H1 error ", math.log(prev_H1_error / H1_error, 2)) 
  prev_H1_error = H1_error 



