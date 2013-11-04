from dolfin import *
from numpy import *
from optimization import *
from geometry import *
from myschur import StokesSolver
from myutilities import *
import sys

# info(parameters,True)
parameters["allow_extrapolation"] = True	# useful for evaluating on boundaries (see footnote on p.116 of the manual)

# ======================================================================
# mesh and spaces definitions
# ======================================================================
Nq = 20
# I = UnitInterval(Nq) # mesh for controls
# I = UnitIntervalMesh(Nq) # mesh for controls !!! better: previous is deprecated
	# !!! usiamo RestrictFunction
Nup = Nq
    # ??? ci vuole un legame tra Nq e Nup? Non direi...
    # ??? idem per i gradi polinomiali
Omega = UnitSquareMesh(Nup,Nup) # mesh for velocity and pressure on reference domain
# Omega = UnitSquareMesh(Nup,Nup) # mesh for velocity and pressure on reference domain !!! better: previous is deprecated
restrictor = RestrictFunction(Omega)	# restriction operator "from state domain to control domain"
I = restrictor.get_restr_mesh()
degq = 2	# polynomial degree for controls FE space Q
	# every q in Q is supposed to be in H2, then we guarantee not to have necessarily q''=0 setting degq>=2
degu = 2	# polynomial degree for velocity FE space V
degp = 1	# polynomial degree for pressure FE space P
# Q = FunctionSpace(I,'CG',degq)
	# !!! usiamo RestrictFunction
Q = restrictor.get_restr_space()
# Q_Omega = FunctionSpace(Omega,'CG',degq) # mai usato...
V = VectorFunctionSpace(Omega,'CG',degu)
P = FunctionSpace(Omega,'CG',degp)
Vexpr = V	# used by Expressions: maybe a high order is better
W = V*P		# used by StokesSolver

# ======================================================================
# parameters for the optimization method
# ======================================================================
# tol = 1.0e-3
tol = DOLFIN_EPS
maxit = 7
alpha0 = 0.1 # for descent methods

# penalty coefficients of the functional
#penalty_coeffs=[1.e-3,1.e-5]
penalty_coeffs=[0.0,0.0]
	# first:	for volume conservation
	# second:	for q'' penalization

# ======================================================================
# data for the state problem
# ======================================================================
eta = Constant('0.0')		# 0-order parameter for generalized Stokes problem
nu = Constant('1.0')#.e-6)	# viscosity
f = Expression(('0.0','0.0'))		# external force
# for boundary conditions
u_in = Expression(('2*x[1]*(1-x[1]/2.0)','0.0'))
g = Expression(('0.0','0.0'))	        	# outflow stress (Neumann datum)

# ======================================================================
# boundary managing
# ======================================================================
Gamma_moving = MovingBoundary()
# Dirichlet bc
noslip = DirichletBC(W.sub(0),Constant(('0.0','0.0')),lambda x, on_boundary : Gamma_moving.inside(x,on_boundary))
inflow = DirichletBC(W.sub(0),u_in,lambda x, on_boundary : near(x[0],0))
symmetry = DirichletBC(W.sub(0).sub(1),Expression('0.0'),lambda x, on_boundary : near(x[1],1))
# noslip = DirichletBC(V,Constant(('0.0','0.0')),lambda x, on_boundary : Gamma_moving.inside(x,on_boundary))
# symmetry = DirichletBC(V.sub(1),Expression('0.0'),lambda x, on_boundary : (not Gamma_moving.inside(x,on_boundary)) and near(x[1],1.0))
# inflow = DirichletBC(V,u_in,lambda x, on_boundary : (not Gamma_moving.inside(x,on_boundary)) and near(x[0],0.0))
bc = [noslip,symmetry,inflow]
# Neumann bc
class NeuBoundary(SubDomain):
	def inside(self,x,on_boundary):
		return on_boundary and near(x[0],1.0)
			# ??? c'e' un modo per far chiamare DirBoundary.inside senza creare un'istanza?
			# perche' altrimenti dovrei fare l'override di qualche metodo affinche' NeuBoundary abbia un oggetto DirBoundary...
neu_boundary = NeuBoundary()

boundaries = FacetFunction("size_t",Omega)
boundaries.set_all(1)
neu_boundary.mark(boundaries,10)
Gamma_moving.mark(boundaries,0)
ds = Measure('ds')[boundaries]

# ======================================================================
# initial settings
# ======================================================================
print "# ======================================================================\n\
# initial settings\n\
# ======================================================================"
q = Function(Q)
# q0 = Expression('0.0')
q0 = Expression('0.2*sin(2*pi*x[0])*sin(2*pi*x[0])')
# q0 = Expression('0.2-0.2*4*(x[0]-0.5)*(x[0]-0.5)')
# q0 = Expression('(x[0]<0.5)*0.5*x[0] + (x[0]>=0.5)*0.5*(1-x[0])')
# q0 = Expression('0.5*x[0]')
# q.vector()[:] = array([q0(x) for x in I.coordinates()])
	# !!! the line above works only for polynomial degree = 1
q = project(q0,Q)
file = File("q0.pvd")
file << q
pfig = plot(q,mesh=I,title="q0(x)")
interactive()
# pfig.write_png("q0.png")
# pfig.write_ps("q0", format="eps")

# DEBUG
#Tq = MapTq(element=Vexpr.ufl_element(),q=q0)
#plot(Tq,mesh=Omega,title="Tq(x,y)")
#plot(Tq-Expression(('x[0]','x[1]')),mesh=Omega,title="Tq(x,y) - I")
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

# methods initialization
optmethod = GradientMethod(tol,maxit)
solver = StokesSolver(W)#(V,P)
adjoint_solver = StokesSolver(W)#(V,P)
solver.set_data(eta,nu,f,bc,g,boundaries)
	# adjoint_solver data cannot be set up a priori: they depend on the state variables

# functional object definition
J = CostFunctional(nu,solver,adjoint_solver,q0=q0,N=Nq)

# ======================================================================
# optimization method applying
# ======================================================================
print "# ======================================================================\n\
# optimization method applying\n\
# ======================================================================"
q_opt = optmethod.apply(J,q,alpha0)
	# ! Be sure to have set q to its initial value
u_opt,p_opt = J.last_state()
z_opt,s_opt = J.last_adjoint_state()

# ======================================================================
# post-processing
# ======================================================================
print "# ======================================================================\n\
# post-processing\n\
# ======================================================================"
pfig = plot(q_opt,mesh=I,title="q_opt(x)",basename="q_opt")
# pfig.write_png("q_opt.png")
# pfig.write_ps("q_opt", format="eps")
pfig = plot(q_opt-q0,mesh=I,title="q_opt(x)-q0(x)")
# pfig.write_png("q_opt-q0.png")
# pfig.write_ps("q_opt-q0", format="eps")
file = File("q_opt.pvd")
file << (q_opt,0.0) 
file << (project(q_opt-q0,Q),10.0)

u_x,u_y = u_opt.split()
pfig = plot(u_opt,title="Optimal velocity vector u")
# pfig.write_png("u.png")
# pfig.write_ps("u", format="eps")
# plot(u_x,title=("Velocity x component"))
# plot(u_y,title=("Velocity y component"))
# pfig = plot(p_opt,title=("Pressure"))
# # pfig.write_png("p.png")
# # pfig.write_ps("p", format="eps")
file = File("state_opt.pvd")
file << (u_opt,0.0)
file << (u_x,10.0)
file << (u_y,20.0)
file << (p_opt,100.0)
z_x,z_y = z_opt.split()
pfig = plot(z_opt,title="Optimal adjoint velocity vector z")
# pfig.write_png("z.png")
# pfig.write_ps("z", format="eps")
# pfig = plot(s_opt,title=("Adjoint pressure"))
# # pfig.write_png("s.png")
# # pfig.write_ps("s", format="eps")
file = File("adjoint_opt.pvd")
file << (z_opt,0.0)
file << (z_x,10.0)
file << (z_y,20.0)
file << (s_opt,100.0)


interactive()
raw_input("Press ENTER to continue")

print ""
print "Bien: On a fini Normalement!"
