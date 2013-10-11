from optimization import *
from geometry import *

my_debug = True
# my_debug = False

# mesh and spaces definitions
Nq = 10
I = UnitInterval(Nq) # mesh for controls
Nup = 20
Omega = UnitSquare(Nup,Nup) # mesh for velocity and pressure on reference domain
degq = 1	# polynomial degree for controls FE space Q
degu = 2	# polynomial degree for velocity FE space V
degp = 1	# polynomial degree for pressure FE space P
Q = FunctionSpace(I,'CG',degq)
Q_Omega = FunctionSpace(Omega,'CG',degq)
V = VectorFunctionSpace(Omega,'CG',degu)
P = FunctionSpace(Omega,'CG',degp)
Vexpr = V	# used by Expressions: maybe a high order is better

# parameters for the optimization method
tol = 1.0e-5
maxit = 100
alpha0 = 0.5 # for descent methods
# penalty coefficients
lambda1 = 1		# for volume conservation
lambda2 = 1		# for q'' penalization

# data for the state problem
nu = Constant(1.0e-6)	# viscosity
eta = Constant(0)		# 0-order parameter for generalized Stokes problem
f = Expression(('0','0'))		# external force
# for boundary conditions
u_in = Expression(('x[1]','0'))		# inflow velocity
stress_out = Expression(('0','0'))	# outflow stress

# functions and initial settings
qsym = TrialFunction(Q)
usym = TrialFunction(V)
q = Function(Q)
u = Function(V)
#q0 = Expression('0.2-0.2*4*(x[0]-0.5)*(x[0]-0.5)')
q0 = Expression('sin(2*pi*x[0])')
u0 = Expression('0')
u.vector()[:] = 0.0
#q = q0
Vbar_form = abs(q)*dx	
Vbar = assemble(Vbar_form)	# volume under q(x)
# TODO GESTIRE q_Omega(x,y) = q(x)

# boundary managing
boundary_parts = MeshFunction("uint", Omega, Omega.topology().dim()-1)
Gamma_moving = MovingBoundary()
Gamma_moving.mark(boundary_parts,0)

# coefficients due to transformation on reference domain
# DEVONO ESSERE ESPRESSIONI CHE DIPENDONO DAL CONTROLLO
Tq = MapTq(element=Vexpr.ufl_element())
	# ??? passandogli questo element ho il grado di interpolazione dato da Q: non e' che facendo il gradiente poi perdo un grado?
Tq.q = q0
# DTq = grad(Tq)
# gammaq = det(DTq)
# DTq_invT = transpose(inv(DTq))
# Aq = gammaq*DTq*DTq_invT

if my_debug == True:
	#u.vector()[:] = as_vector([Tq(Omega.coordinates()[i]) for i in range(0,len(Omega.coordinates()))])
	plot(Tq,mesh=Omega,title="Tq(x,y)")
	plot(Tq-Expression(('x[0]','x[1]')),mesh=Omega,title="T\_q(x,y) - I")
	interactive()

# # mesh utilities
# n = FacetNormal(Omega)

# cost functional J and its gradient
# FOR DERIVATIVES WE COULD USE UFL
cost_functional_form = 0.5*inner(grad(u)*Aq,grad(u))*dx
penalty_terms = 0.5*lambda1*((q*dx) - Vbar)**2 \
				+ 0.5*lambda2*Dx(Dx(q,0),0)*Dx(Dx(q,0),0)*dx
	# ??? non so se il primo pezzo funziona, scritto cosi'...
J_form = cost_functional_form + penalty_terms
cost_functional_gradient = - 0.5*inner(grad(u)*DTq,grad(u)*mapTq.Jac)*ds(0)\
							 - dot(grad(u)*DTq*DTq_invT*n, (grad(u)-nu*grad(z))*DTq*DTq_invT*n + s*n)*ds(0)
penalty_gradient = lambda1*((q*dx)-Vbar) + lambda2*Dx(Dx(Dx(Dx(q,0),0),0),0)*dx
	# ??? va bene con la derivata 4'
gradJ = cost_functional_gradient + penalty_gradient

# initializations
optmethod = GradientMethod(tol,maxit)
solver = StokesSolver(V,P)
solver.assemble(nu,f,Expression('0','0'))
	# TODO per ora lasciamo questa Expression: poi ci andranno le bc

# functional object definition
J = CostFunctional(J_form,solver)

# apply optimization method
optmethod.set_gradf(gradJ)	# TODO per ora lo mettiamo qui, ma dev'esserci un metodo gradiente nella classe di J
q_opt = optmethod.apply(J,q0)
u_opt,p_opt = J.last_state()

# post-processing


raw_input("Press ENTER to continue")
print ""
print "Bien: On a fini Normalement!"