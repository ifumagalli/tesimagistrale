from optimization import *
from geometry import *

# mesh and spaces definitions
Nq = 10
I = UnitInterval(Nq) # mesh for controls
Nup = 20
Omega0 = UnitSquare(Nup,Nup) # mesh for velocity and pressure on reference domain
degq = 1	# polynomial degree for controls FE space Q
degu = 2	# polynomial degree for velocity FE space V
degp = 1	# polynomial degree for pressure FE space P
Q = FunctionSpace(I,'CG',degq)
V = VectorFunctionSpace(Omega0,'CG',degu)
P = FunctionSpace(Omega0,'CG',degp)

# parameters for the optimization method
tol = 1.0e-5
maxit = 100
alpha0 = 0.5 # for descent methods

# data for the state problem
nu = Constant(1.0e-6)	# viscosity
eta = Constant(0)		# 0-order parameter for generalized Stokes problem
f = Expression(('0','0'))		# external force
# for boundary conditions
u_in = Expression(('y','0'))		# inflow velocity
stress_out = Expression(('0','0'))	# outflow stress

# functional


# initializations
solver = StokesSolver(V,P)
optmethod = GradientMethod(f,x0,tol,maxit,alpha0)