from dolfin import *
from optimization import *
from geometry import *
from numpy import array
from myschur import StokesSolver
	# ??? dove viene usata la PETSc

# my_debug = True
my_debug = False

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
# penalty coefficients
lambda1 = 1		# for volume conservation
lambda2 = 1		# for q'' penalization

# data for the state problem
nu = Constant(1)	# viscosity
eta = Constant(0)		# 0-order parameter for generalized Stokes problem
f = Expression(('0.0','0.0'))		# external force
# for boundary conditions
u_in = Expression(('sqrt(x[1])','0'))		# inflow velocity
stress_out = Expression(('0','0'))	        # outflow stress

# functions and initial settings
qt = TrialFunction(Q)
ut = TrialFunction(V)
pt = TrialFunction(P)
zt = TrialFunction(V)
st = TrialFunction(P)
q = Function(Q)
u = Function(V)
z = Function(V)
s = Function(P)
#q0 = Expression('0.2-0.2*4*(x[0]-0.5)*(x[0]-0.5)')
q0 = Expression('sin(2*pi*x[0])')
u0 = Expression(('0.0','0.0'))
# q.vector()[:] = array([q0(x) for x in I.coordinates()])
	# !!! the line above works only for polynomial degree = 1
q = project(q0,Q)

# boundary managing
# boundary_parts = MeshFunction("uint", Omega, Omega.topology().dim()-1)
# boundary_parts = MeshFunction("size_t", Omega, Omega.topology().dim()-1) # uint as MeshFunction type has been deprecated since DOLFIN 1.1.0
	# ??? sarebbe meglio usare questa MeshFunction, invece dei SubDomain sotto definiti?
Gamma_moving = MovingBoundary()
# Gamma_still = StillBoundary()
# Dirichlet bc
noslip = DirichletBC(V,Expression(('0','0')),lambda x, on_boundary : Gamma_moving.inside(x,on_boundary))
inflow = DirichletBC(V,u_in,lambda x, on_boundary : (not Gamma_moving.inside(x,on_boundary)) and near(x[0],0))#x[0] < DOLFIN_EPS)
# noslip2 = DirichletBC(V,Expression(('0','0')),lambda x, on_boundary : Gamma_still.inside(x,on_boundary) and x[1] > 1-DOLFIN_EPS)
symmetry = DirichletBC(V.sub(1),Expression('0'),lambda x, on_boundary : (not Gamma_moving.inside(x,on_boundary)) and near(x[1],1))#x[1] > 1-DOLFIN_EPS)
	# ??? DA VERIFICARE
bc = [noslip,inflow,symmetry]


# coefficients due to transformation on reference domain
# DEVONO ESSERE ESPRESSIONI CHE DIPENDONO DAL CONTROLLO
Tq = MapTq(element=Vexpr.ufl_element())
	# ??? passandogli questo element ho il grado di interpolazione dato da Vexpr: non e' che facendo il gradiente poi perdo un grado?
Tq.q = q0
DTq = grad(Tq)
gammaq = det(DTq)
DTq_invT = transpose(inv(DTq))
Aq = gammaq*DTq*DTq_invT

if my_debug == True:
	#u.vector()[:] = as_vector([Tq(Omega.coordinates()[i]) for i in range(0,len(Omega.coordinates()))])
	plot(q,mesh=I,title="q(x)")
	plot(Tq,mesh=Omega,title="Tq(x,y)")
	plot(Tq-Expression(('x[0]','x[1]')),mesh=Omega,title="Tq(x,y) - I")
	plot(gammaq,mesh=Omega,title="gamma_q(x,y)")
	interactive()

# mesh utilities
n = FacetNormal(Omega)

# cost functional J and its gradient
# FOR DERIVATIVES WE COULD USE UFL
cost_functional_form = 0.5*inner(grad(u)*Aq,grad(u))*dx
# penalty_terms = 0.5*lambda1*((q*dx) - Vbar)**2 \
#				+ 0.5*lambda2*Dx(Dx(q,0),0)*Dx(Dx(q,0),0)*dx
	# ??? non so se il primo pezzo funziona, scritto cosi'... per ora proviamo il seguente
    # quando decommenti sopra, decommenta anche penalty_gradient!!
penalty_terms = 0.5*lambda2*Dx(Dx(q,0),0)*Dx(Dx(q,0),0)*dx
J_form = cost_functional_form #+ penalty_terms
# J_form = (u - u0)*(u - u0)*dx
### J_form = dot(u-u0,u-u0)*dx
    # TODO per ora teniamo questo funzionale elementare, per vedere se le cose vanno...
cost_functional_gradient = - 0.5*inner(grad(u)*DTq,grad(u)*DTq)*ds(0)\
							 - dot(grad(u)*DTq*DTq_invT*n, (grad(u)-nu*grad(z))*DTq*DTq_invT*n + s*n)*ds(0)
# penalty_gradient = lambda1*((q*dx)-Vbar) + lambda2*Dx(Dx(Dx(Dx(q,0),0),0),0)*dx
penalty_gradient = lambda2*Dx(Dx(Dx(Dx(qt,0),0),0),0)*dx
	# ??? va bene con la derivata 4'?
gradJ_form = cost_functional_gradient #+ penalty_gradient
### gradJ_form = -inner(grad(z), w0*grad(u))*dx

# initializations
optmethod = GradientMethod(tol,maxit)
solver = StokesSolver(V,P)
#solver.assemble(nu,f,bc,q)
	# TODO per ora lasciamo questa Expression: poi ci andranno le bc

# functional object definition
J = CostFunctional(nu,f,bc,J_form,solver,grad_form=gradJ_form,q0=q0)
# apply optimization method
# optmethod.set_gradf(gradJ)	# TODO per ora lo mettiamo qui, ma dev'esserci un metodo gradiente nella classe di J
				# DONE messo in CostFunctional
q_opt = optmethod.apply(J,q)
	# ! Be sure to have set q to its initial value
u_opt,p_opt = J.last_state()

# post-processing
print "\nPost-processing"
# plot(q_opt,mesh=I,title="q(x)")
# plot(Tq,mesh=Omega,title="Tq(x,y)")
# plot(Tq-Expression(('x[0]','x[1]')),mesh=Omega,title="Tq(x,y) - I")
# plot(gammaq,mesh=Omega,title="gamma_q(x,y)")
# interactive()
raw_input("Press ENTER to continue")

print ""
print "Bien: On a fini Normalement!"
