#!/usr/bin/env python

from dolfin import *
from geometry import *
from myutilities import *

# parameters['linear_algebra_backend'] = 'PETSc'
parameters['linear_algebra_backend'] = 'uBLAS'

class StokesSolver(object):
  # TODO vedere se qualcosa va messo private
  # constructor
  def __init__(self, W):
	# self.V = V
	# self.P = P
	self.W = W

  def set_data(self,eta,nu,f,bc,g=None,boundaries=None,f_gradv=False):
	# TODO forse meglio un generico passaggio per kwargs?
	# si potrebbe tenere un self.kwargs da cui estrarre al momento dell'assemble
	self.eta = eta
	self.nu = nu
	self.f = f
	self.bc = bc
	self.g = g
	self.boundaries = boundaries
	self.f_gradv = f_gradv

  # assemble matrix
  def assemble(self, q):
	if self.eta==None or self.nu==None or self.f==None or self.bc==None:
		exit("@@@ EXITING: data setting is required before assembling")
	if self.g == None and self.boundaries != None:
		print("@@@ WARNING: 'boundaries' FacetFunction passed won't be used")
	if self.g != None and self.boundaries == None:
		exit("@@@ EXITING: Neumann data need a FacetFunction ('boundaries')")
	
	u,p = TrialFunctions(self.W)
	v,r = TestFunctions(self.W)

	# coefficients due to control
	Tq = MapTq(element=self.W.sub(0).ufl_element(),q=q)
		# ??? meglio fare un membro Tq anche in questa classe? in CostFunctional c'e' gia'
	DTq = grad(Tq)
		# ??? meglio grad o nabla_grad?
	gammaq = det(DTq)
	DTq_invT = transpose(inv(DTq))
	Aq = gammaq*DTq*DTq_invT

# DEBUG
#	plot(Tq,mesh=self.W.mesh(),title="Tq(x,y) in myschur",key="Tq_schur")
#	plot(Tq-Expression(('x[0]','x[1]')),mesh=self.W.mesh(),title="Tq(x,y) - I in myschur",key="Tq_schur")
# END DEBUG

	# assemble A matrix
	ds = Measure('ds')[self.boundaries]
	a = self.nu*inner(grad(u)*Aq,grad(v)) * dx \
	  + inner(grad(u),DTq_invT)*r*gammaq * dx \
	  - inner(grad(v),DTq_invT)*p*gammaq * dx \
	  + Constant('0.0')*dot(u,v) * ds(10)
		# ??? altrimenti da' un warning e non usa ds in L
	#a = self.nu*inner(grad(u),grad(v)) * dx \
	#  - p * div(v) * dx \
	#  + r * div(u) * dx \
	#  + Constant('0.0')*dot(u,v) * ds(10)
	#L = None
	if self.f_gradv == False:
		L = dot(self.f, v)*gammaq * dx + dot(self.g,v) * ds(10)
		#L = dot(self.f, v) * dx + dot(self.g,v) * ds(10)
	else:
		if isinstance(self.f,Expression):
			self.f = project(self.f,self.V)
			# grad() can't be applied on a generical Expression: it needs gdim
		L = inner(self.f*Aq, grad(v)) * dx + dot(self.g,v) * ds(10)
		#L = inner(self.f, grad(v)) * dx + dot(self.g,v) * ds(10)
	# self.A, self.b = assemble_system(a, L, self.bc, exterior_facet_domains=self.boundaries)
	self.A, self.b = assemble_system(a, L, self.bc)
	# self.A = assemble(a)#, exterior_facet_domains=self.boundaries)
	# self.b = Vector(self.A.size(0))
	# self.b.zero()
	# self.b = assemble(L)#, exterior_facet_domains=self.boundaries)
	# for item in self.bc:
	#	item.apply(self.A,self.b)

	# assemble A preconditioner
	# a = self.nu * inner(grad(u),grad(v)) * dx \
	#	+ p * r / self.nu * dx
	# self.Ap = assemble(a)

  def solve(self):
	x = Function(self.W)
	
	#solver = KrylovSolver('gmres', 'ilu')
	#solver.set_operators(self.A, self.Ap)
	solver = LUSolver()
	solver.set_operator(self.A)
	solver.solve(x.vector(), self.b)

	#return x.split()
	return x
	# Se vuoi che u,p siano modificabili all'esterno senza dover fare deep copies, usa invece le righe seguenti
	# u = Function(self.V)
	# p = Function(self.P)
	# u = assign(x.split()[0])
	# p = assign(x.split()[1])
	# return u,p
	
  def mesh(self):
	return self.W.mesh()

# =====================


# ===== Test code =====

if __name__=="__main__":

  N = 10
  mesh = UnitSquare(N,N,'crossed')
  plot(mesh)
  V = VectorFunctionSpace(mesh, 'CG', 2)
  # W = FunctionSpace(mesh, 'CG', 1)
  # B = FunctionSpace(mesh, 'B', 3)
  # V = W + B
  P = FunctionSpace(mesh, 'CG', 1)
  # W = (V*V)*P
  W = V*P

  # problem data
  eta = Constant(0.0)
  nu = Constant(1)#.e-6)
  f = Expression(('0.0','0.0'))

  # Dirichlet bc
  #u_in = Expression(('sqrt(x[1])','0.0'))       # inflow velocity
  u_in = Expression(('100*2*x[1]*(1-x[1]/2.0)','0.0'))
  Gamma_moving = MovingBoundary()
  noslip = DirichletBC(W.sub(0),Constant(('0.0','0.0')),lambda x, on_boundary : Gamma_moving.inside(x,on_boundary))
  #inflow = DirichletBC(W.sub(0),u_in,lambda x, on_boundary : (not Gamma_moving.inside(x,on_boundary)) and near(x[0],0))
  #symmetry = DirichletBC(W.sub(0).sub(1),Expression('0.0'),lambda x, on_boundary : (not Gamma_moving.inside(x,on_boundary)) and near(x[1],1))
  # PRESSIONE RISOLTA!
  inflow = DirichletBC(W.sub(0),u_in,lambda x, on_boundary : near(x[0],0))
  symmetry = DirichletBC(W.sub(0).sub(1),Expression('0.0'),lambda x, on_boundary : near(x[1],1))
  bc = [noslip,symmetry,inflow]
  # Neumann bc
  g = Expression(('0.0','0.0'))# boundary conditions
  class NeuBoundary(SubDomain):
	def inside(self,x,on_boundary):
		return on_boundary and near(x[0],1.0)
			# ??? c'e' un modo per far chiamare DirBoundary.inside senza creare un'istanza?
			# perche' altrimenti dovrei fare l'override di qualche metodo affinche' NeuBoundary abbia un oggetto DirBoundary...
  neu_boundary = NeuBoundary()
  # FacetFunction and ds settings
  boundaries = FacetFunction("uint",mesh)
  boundaries.set_all(1)
  neu_boundary.mark(boundaries,10)
  Gamma_moving.mark(boundaries,0)
  ds = Measure('ds')[boundaries]

  # control
  restrictor = RestrictFunction(mesh)	# restriction operator "from state domain to control domain"
  Q = restrictor.get_restr_space()
  #Q = FunctionSpace(UnitInterval(N),'CG',2)
  q = Function(Q)
  # q0 = Expression('0.0')
  # q0 = Expression('sin(2*pi*x[0])')
  # q0 = Expression('0.2-0.2*4*(x[0]-0.5)*(x[0]-0.5)')
  q0 = Expression('0.2*x[0]')
  q = project(q0,Q)
  
  # solver initializing
  stokes = StokesSolver(W)
  stokes.set_data(eta, nu, f, bc, g, boundaries)
  
  # solution
  stokes.assemble(q)
  x = Function(W)
  x = stokes.solve()
  u,p = x.split()

  # post-processing
  u_x,u_y = u.split()
  plot(u,title=("Velocity vector u"))
  # plot(u_x,title=("Velocity x component"))
  # plot(u_y,title=("Velocity y component"))
  plot(p,title=("Pressure"))

	#  interactive()
#  raw_input("Press ENTER")
#  restrictor0 = RestrictFunction(mesh)#geometry.moving_boundary):
#  class Gamma1(SubDomain):
#	def inside(self,x,on_boundary):
#		return near(x[0],0.0)
#  gamma1 = Gamma1()
#  restrictor1 = RestrictFunction(mesh,gamma1)#geometry.moving_boundary):
#
#  file = File("u.pvd")
#  file << u
#
#  u_x,u_y = u.split()
#  u_x_copy = Function(FunctionSpace(mesh,'CG',2))
#  u_x_copy.vector()[:] = u_x.vector()
#  u_x_rest0 = restrictor0(u_x_copy)
# 
#  u_rest0, p_rest0 = restrictor0(x).split()
#  u_rest1, p_rest1 = restrictor1(x).split()
#  u_x_rest0,u_y_rest0 = u_rest0.split()
#  u_x_rest1,u_y_rest1 = u_rest1.split()
#
#  file = File("ux0.pvd")
#  file << u_x_rest0
#  file = File("uy0.pvd")
#  file << u_y_rest0
#  file = File("ux1.pvd")
#  file << u_x_rest1
#  file = File("uy1.pvd")
#  file << u_y_rest1

  
  # Tq = MapTq(element=V.ufl_element(),q=q)
  # u = TildeFunction(fun=u,map=Tq)
  # p = TildeFunction(fun=p,map=Tq)
  # plot(u,mesh=mesh,title=("Velocity vector u"))
  # # plot(u_x,title=("Velocity x component"))
  # # plot(u_y,title=("Velocity y component"))
  # plot(p,mesh=mesh,title=("Pressure"))

  # ========================
  # aggiunto
  # ========================
  n = FacetNormal(mesh)
  f_adj = grad(u)
  noslip = DirichletBC(W.sub(0),Constant(('0.0','0.0')),lambda x, on_boundary : Gamma_moving.inside(x,on_boundary))
  # symmetry = DirichletBC(W.sub(0).sub(1),Expression('0.0'),lambda x, on_boundary : (not Gamma_moving.inside(x,on_boundary)) and near(x[1],1.0))
  # inflow = DirichletBC(W.sub(0),Constant(('0.0','0.0')),lambda x, on_boundary : (not Gamma_moving.inside(x,on_boundary)) and near(x[0],0.0))
  # PRESSIONE RISOLTA!
  inflow = DirichletBC(W.sub(0),Constant(('0.0','0.0')),lambda x, on_boundary : near(x[0],0))
  symmetry = DirichletBC(W.sub(0).sub(1),Expression('0.0'),lambda x, on_boundary : near(x[1],1))
  bc_adj = [noslip,symmetry,inflow]
  #g_adj = p*n-g
  g_adj = Constant(('0.0','0.0'))
	# !!! perche' e' annullata dal True

  stokes.set_data(eta,nu,f_adj,bc_adj,g_adj,boundaries,True)
  stokes.assemble(q)
  z,s = stokes.solve().split()

  plot(z,title="Adjoint velocity")
  plot(s,title="Adjoint pressure")

  # ========================
  # traccia
  # ========================
  dq_uzs = TrialFunction(P)
  phi = TestFunction(P)
  Tq = MapTq(element=V.ufl_element(),q=q)
  DTq = grad(Tq)
  DTq_invT = transpose(inv(DTq))
  D = lambda u : Dx(u,0)

  #boundaries = FacetFunction("uint",mesh)
  #boundary = MovingBoundary()
  #boundaries.set_all(1)
  #boundary.mark(boundaries,0)
  #ds = Measure('ds')[boundaries]

  # trace_form = dq_uzs * phi * ds(0)
  trace_form = dot(grad(dq_uzs),grad(phi))*dx + Constant('100.0')*dq_uzs * phi * ds(0)
	# ??? ha senso questa sorta di rilevamento?
  trace_rhs = Constant('100.0')*(-0.5)*inner(grad(u)*DTq,grad(u)*DTq) * phi * ds(0) \
		+ Constant('100.0')*dot(grad(u)*DTq*DTq_invT*n, (grad(u)+nu*grad(z))*DTq*DTq_invT*n - s*n) * phi * ds(0)
  value_from_state_and_adjoint = Function(P)
  trace_A = assemble(trace_form)
  trace_A.ident_zeros()
  trace_b = assemble(trace_rhs)

  print "\nsolving trace problem\n"
  solver = LUSolver()
  solver.set_operator(trace_A)
  solver.solve(value_from_state_and_adjoint.vector(), trace_b)
  #solve(trace_A, value_from_state_and_adjoint.vector(), trace_b)
  plot(value_from_state_and_adjoint,title="value_from_state_and_adjoint",key="adjoint")

#  for_value_from_state_and_adjoint = -0.5*inner(grad(u)*DTq,grad(u)*DTq) \
#			+ dot(grad(u)*DTq*DTq_invT*n, (grad(u)+nu*grad(z))*DTq*DTq_invT*n - s*n)
#  for_restr_form = dot(D(dq_uzs),D(phi)) * dx
#  for_restr_rhs = Constant('0.0')*phi * dx
#  for_restr_bc = DirichletBC(P,for_value_from_state_and_adjoint,lambda x,on_boundary : on_boundary)
#  solve(for_restr_form==for_restr_rhs, value_from_state_and_adjoint, for_restr_bc)

  func_value_from_state_and_adjoint = restrictor(value_from_state_and_adjoint)
	# ! Using myutilities.MyRestrictedFunction to return a Function of the same type as the control q
  pfig = plot(func_value_from_state_and_adjoint,title="func_value_from_state_and_adjoint",key="func_adjoint")


  interactive()
  raw_input("Press ENTER to finish")
