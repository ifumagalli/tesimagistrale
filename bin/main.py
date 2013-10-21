from dolfin import *
from optimization import *
from geometry import *
from numpy import array
from myschur import StokesSolver
from myutilities import *

# info(parameters,True)
parameters["allow_extrapolation"] = True	# useful for evaluating on boundaries (see footnote on p.116 of the manual)

# mesh and spaces definitions
Nq = 20
I = UnitInterval(Nq) # mesh for controls
# I = UnitIntervalMesh(Nq) # mesh for controls !!! better: previous is deprecated
Nup = Nq
    # ??? ci vuole un legame tra Nq e Nup? Non direi...
    # ??? idem per i gradi polinomiali
Omega = UnitSquare(Nup,Nup) # mesh for velocity and pressure on reference domain
# Omega = UnitSquareMesh(Nup,Nup) # mesh for velocity and pressure on reference domain !!! better: previous is deprecated
degq = 2	# polynomial degree for controls FE space Q
	# every q in Q is supposed to be in H2, then we guarantee not to have necessarily q''=0 setting degq>=2
degu = 2	# polynomial degree for velocity FE space V
degp = 1	# polynomial degree for pressure FE space P
Q = FunctionSpace(I,'CG',degq)
Q_Omega = FunctionSpace(Omega,'CG',degq)
V = VectorFunctionSpace(Omega,'CG',degu)
P = FunctionSpace(Omega,'CG',degp)
Vexpr = V	# used by Expressions: maybe a high order is better

# parameters for the optimization method
# tol = 1.0e-3
tol = DOLFIN_EPS
maxit = 2
alpha0 = 0.5 # for descent methods

# penalty coefficients of the functional
penalty_coeffs=[1.0,1.0]
	# first:	for volume conservation
	# second:	for q'' penalization

# data for the state problem
eta = Constant(0.0)		# 0-order parameter for generalized Stokes problem
nu = Constant(1.e-6)	# viscosity
f = Expression(('0.0','0.0'))		# external force
g = Expression(('0.0','0.0'))		# Neumann datum
# for boundary conditions
u_in = Expression(('sqrt(x[1])','0.0'))		# inflow velocity
stress_out = Expression(('0.0','0.0'))	        # outflow stress

# functions and initial settings
q = Function(Q)
u = Function(V)
z = Function(V)
s = Function(P)
# q0 = Expression('0.0')
# q0 = Expression('sin(2*pi*x[0])')
q0 = Expression('0.2-0.2*4*(x[0]-0.5)*(x[0]-0.5)')
# q.vector()[:] = array([q0(x) for x in I.coordinates()])
	# !!! the line above works only for polynomial degree = 1
q = project(q0,Q)

# boundary managing
Gamma_moving = MovingBoundary()
# Dirichlet bc
noslip = DirichletBC(V,Constant(('0.0','0.0')),lambda x, on_boundary : Gamma_moving.inside(x,on_boundary))
inflow = DirichletBC(V,u_in,lambda x, on_boundary : (not Gamma_moving.inside(x,on_boundary)) and near(x[0],0.0))
symmetry = DirichletBC(V.sub(1),Expression('0.0'),lambda x, on_boundary : (not Gamma_moving.inside(x,on_boundary)) and near(x[1],1.0))
bc = [noslip,inflow,symmetry]
# Neumann bc
class NeuBoundary(SubDomain):
	def inside(self,x,on_boundary):
		return on_boundary and near(x[0],1.0)
			# ??? c'e' un modo per far chiamare DirBoundary.inside senza creare un'istanza?
			# perche' altrimenti dovrei fare l'override di qualche metodo affinche' NeuBoundary abbia un oggetto DirBoundary...
boundaries = FacetFunction("uint",Omega)
neu_boundary = NeuBoundary()
neu_boundary.mark(boundaries,10)

# DEBUG
if MYDEBUG == True:
	# coefficients due to transformation on reference domain
	# DEVONO ESSERE ESPRESSIONI CHE DIPENDONO DAL CONTROLLO
	Tq = MapTq(element=Vexpr.ufl_element(),q=q0)
		# ??? passandogli questo element ho il grado di interpolazione dato da Vexpr: non e' che facendo il gradiente poi perdo un grado?
	DTq = grad(Tq)
	gammaq = det(DTq)
	DTq_invT = transpose(inv(DTq))
	Aq = gammaq*DTq*DTq_invT

	plot(q,mesh=I,title="q(x)")
	plot(Tq,mesh=Omega,title="Tq(x,y)")
	plot(Tq-Expression(('x[0]','x[1]')),mesh=Omega,title="Tq(x,y) - I")
	plot(gammaq,mesh=Omega,title="gamma_q(x,y)")
	interactive()
# END DEBUG

# initializations
optmethod = GradientMethod(tol,maxit)
solver = StokesSolver(V,P)
adjoint_solver = StokesSolver(V,P)
solver.set_data(eta,nu,f,bc,g,boundaries)
	# adjoint_solver data cannot be set up a priori: they depend on the state variables

# functional object definition
J = CostFunctional(nu,solver,q0=q0)

# apply optimization method
q_opt = optmethod.apply(J,q)
	# ! Be sure to have set q to its initial value
u_opt,p_opt = J.last_state()

# post-processing
print "\nPost-processing"
# DEBUG
if MYDEBUG == True:
	plot(q_opt,mesh=I,title="q(x)")
	plot(Tq,mesh=Omega,title="Tq(x,y)")
	plot(Tq-Expression(('x[0]','x[1]')),mesh=Omega,title="Tq(x,y) - I")
	plot(gammaq,mesh=Omega,title="gamma_q(x,y)")
	interactive()
# END DEBUG
raw_input("Press ENTER to continue")

print ""
print "Bien: On a fini Normalement!"
