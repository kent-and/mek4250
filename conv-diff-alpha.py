
from dolfin import *

import math

prev_L2_error = 0 
prev_H1_error = 0 
N = 32
for alpha_value in [0.3**p for p in [4, 3, 2, 1, 0]]:

  print ("alpha ", alpha_value) 

  mesh = UnitSquareMesh(N, N)
  V = FunctionSpace(mesh, "CG", 1)
  V3 = FunctionSpace(mesh, "CG", 3)

  u = TrialFunction(V)
  v = TestFunction(V)

  alpha = Constant(alpha_value)
  f = Constant(0)
  h = mesh.hmin()

  a = (-u.dx(0)*v + alpha*inner(grad(u),grad(v)))*dx  
  L = f*v*dx  

  u_analytical = Expression("(exp(-x[0]/%e) - 1)/ (exp(-1/%e) - 1)" % (alpha_value, alpha_value), degree=3 )  
  def boundary(x, on_boundary):
      return on_boundary and (x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS)

  bc = DirichletBC(V, u_analytical, boundary) 

  U = Function(V)
  solve(a == L, U, bc) 

  U_analytical = project(u_analytical, V3)

  # notice that when computing the rate vrt alpha we should not use log! 

  L2_error = errornorm(U, u_analytical) 
  print (" N = %d, L2 error = %e " % (N, L2_error))
  if prev_L2_error > 0: 
    print "rate of L2 error ", prev_L2_error / L2_error 
  prev_L2_error = L2_error 

  H1_error = errornorm(U, u_analytical, "H1") 
  print (" N = %d, H1 error = %e " %  (N, H1_error))
  if prev_H1_error > 0: 
    print "rate of H1 error ", prev_H1_error / H1_error 
  prev_H1_error = H1_error 



