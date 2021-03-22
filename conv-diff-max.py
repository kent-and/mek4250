
from dolfin import *

import math

prev_Linf_error = 0 
prev_dLinf_error = 0 
for N in [8,16,32,64]:

  mesh = UnitSquareMesh(N, N)
  V = FunctionSpace(mesh, "CG", 2)
  V3 = FunctionSpace(mesh, "CG", 4)

  u = TrialFunction(V)
  v = TestFunction(V)

  alpha_value = 1.0e-2 
  alpha_value = 1.0 
  alpha = Constant(alpha_value)
  f = Constant(0)
  h = mesh.hmin()

  a = (-u.dx(0)*v + alpha*inner(grad(u),grad(v)))*dx  
  L = f*v*dx  

  u_analytical = Expression("(exp(-x[0]/%e) - 1)/ (exp(-1/%e) - 1)" % (alpha_value, alpha_value), degree=3 )  
  du_analytical = Expression("(-(1.0/%e)*exp(-x[0]/%e))/ (exp(-1/%e) - 1)" % (alpha_value, alpha_value, alpha_value), degree=3 )  
  def boundary(x, on_boundary):
      return on_boundary and (x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS)

  bc = DirichletBC(V, u_analytical, boundary) 

  U = Function(V)
  solve(a == L, U, bc) 

  U_analytical = project(u_analytical, V3)

  difference = project(u_analytical - U, V3) 
  Linf_error = max(difference.vector().max(), -difference.vector().min())
  print (" N = %d, Linf error = %e " % (N, Linf_error))
  if prev_Linf_error > 0: 
    print "rate of Linf error ", math.log(prev_Linf_error / Linf_error, 2) 
  prev_Linf_error = Linf_error 

  
  difference = project(du_analytical - inner(Constant((1,0)), grad(U)), V3) 
  print (difference)
  dLinf_error = max(difference.vector().max(), -difference.vector().min())
  print (" N = %d, dLinf error = %e " % (N, dLinf_error))
  if prev_dLinf_error > 0: 
    print "rate of dLinf error ", math.log(prev_dLinf_error / dLinf_error, 2) 
  prev_dLinf_error = dLinf_error 


