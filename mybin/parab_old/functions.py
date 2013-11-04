""" Functions to be optimized by an OptimizationMethod.
	TODO nel costruttore c'e' troppa roba: snellire
	TODO rivedere se c'e' un modo per snellire __call__ e gradient
"""

from dolfin import *
import numpy
import copy
from geometry import *
from myutilities import *

class CostFunctional(object):
	def __init__(self,nu,solver,adjoint_solver=None,penalty_coeffs=[0.0,0.0],**kwargs):
		self._nu = nu
		self._solver = solver
		if (adjoint_solver == None):
			print "\n@@@ adjoint_solver is None"
			self._adjoint_solver = solver
		else:
			self._adjoint_solver = adjoint_solver
		up = Function(solver.W)
		zs = Function(self._adjoint_solver.W)
		self._u, self._p = up.split()
		self._z, self._s = zs.split()
		# self._u = Function(solver.W.sub(0))	# state velocity
		# self._p = Function(solver.W.sub(1))	# state pressure
		# self._z = Function(self._adjoint_solver.W.sub(0))	# adjoint state velocity
		# self._s = Function(self._adjoint_solver.W.sub(1))	# adjoint state pressure
		self._penalty_coeffs = penalty_coeffs
		# Vbar assembling:	a value can be passed, or an Expression q whose volume will be Vbar
		#					if both are passed, the numerical value will be used
		coeffN = kwargs.get('N')
		if coeffN == None:
			coeffN = 100
		# else:
			# coeffN *= 2
		self._coeffspace = FunctionSpace(UnitIntervalMesh(coeffN),'CG',5)
		self._coeffVspace = VectorFunctionSpace(UnitSquareMesh(coeffN,coeffN),'CG',5)
			# Space to which Expressions will be projected if a mesh is needed
			# TODO per ora usiamo FunctionSpace(UnitInterval(100),'CG',10)
			
		# @@@ setting Vbar
		if MYVERBOSE == True:
			print "\n@@@ setting Vbar"
		self._Vbar = kwargs.get('Vbar',None)
		q0 = kwargs.get('q0',None)
		if self._Vbar==None:
			if q==None:
				exit("@@@ EXITING: Required fixed volume value or initial control to initialize Vbar ")
			else:
				qf = Function(self._coeffspace)
				qf = project(q0,self._coeffspace)
				self._Vbar = assemble(qf * dx)
		print("Vbar = {0}".format(self._Vbar))
		self._Tq = MapTq(element=self._coeffVspace.ufl_element(),q=q0)
		
		
	def __call__(self,q):
		""" Given a control Function q, it returns the functional evaluated in q.
		"""
		if MYVERBOSE == True:
			print "\n@@@@@@@ functions.CostFunctional.__call__"
		
		# @@@ state problem assembling
		if MYVERBOSE == True:
			print "\n@@@ state problem assembling"
		self._solver.assemble(q)
		
		# @@@ state solution computing
		if MYVERBOSE == True:
			print "\n@@@ state solution computing"
		self._u,self._p = self._solver.solve().split()
		
	# DEBUG
		if MYDEBUG == True:
			u_x,u_y = self._u.split()
			plot(self._u,title=("Velocity vector u"))
			# plot(u_x,title=("Velocity x component"))
			# plot(u_y,title=("Velocity y component"))
			# plot(self._p,title=("Pressure"))
			interactive()
	# END DEBUG
		
		# self._z = None			# to evaluate the functional we don't need the adjoint state
		# self._s = None
		
		# @@@ Functional assembling
		if MYVERBOSE == True:
			print "\n@@@ Functional assembling"
		self._Tq.q = q
		DTq = grad(self._Tq)
		gammaq = det(DTq)
		DTq_inv = inv(DTq)
		Aq = gammaq * DTq_inv * transpose(DTq_inv)
		D = lambda u : Dx(u,0)	# convenient for code readability
		#ux,uy = self._u.split()
		#gradux = grad(ux)
		#graduy = grad(uy)
		#plot(gradux,mesh=self._solver.mesh(),title="gradux")
		#plot(graduy,mesh=self._solver.mesh(),title="graduy")
		value = 0.5*( \
				assemble(inner(grad(self._u)*Aq,grad(self._u)) * dx) \
				+ self._penalty_coeffs[0]*(assemble(q * dx) - self._Vbar)**2 \
				+ self._penalty_coeffs[1]*norm(project(D(q),q.function_space()),'H10')**2 
				)
		return value

	def gradient(self,q):
		""" Returns the gradient of the functional, intended as the
			L2-Riesz representative of the functional Fre'chet derivative computed in q.
		"""
		if MYVERBOSE == True:
			print "\n@@@@@@@ functions.CostFunctional.gradient"
		
		# @@@ solver assembling
		if MYVERBOSE == True:
			print "\n@@@ solver assembling"
		self._solver.assemble(q)
		
		# @@@ state solution computing
		if MYVERBOSE == True:
			print "\n@@@ state solution computing"
		self._u,self._p = self._solver.solve().split()
		
		# @@@ adjoint state problem setting
		if MYVERBOSE == True:
			print "\n@@@ adjoint state problem setting"
		n = FacetNormal(self._solver.mesh())
		f_adj = grad(self._u)
		Gamma_moving = MovingBoundary()
		noslip = DirichletBC(self._adjoint_solver.W.sub(0),Constant(('0.0','0.0')),lambda x, on_boundary : Gamma_moving.inside(x,on_boundary))
		# symmetry = DirichletBC(self._adjoint_solver.W.sub(0).sub(1),Expression('0.0'),lambda x, on_boundary : (not Gamma_moving.inside(x,on_boundary)) and near(x[1],1.0))
		# inflow = DirichletBC(self._adjoint_solver.W.sub(0),Constant(('0.0','0.0')),lambda x, on_boundary : (not Gamma_moving.inside(x,on_boundary)) and near(x[0],0.0))
		# PRESSIONE RISOLTA!
		inflow = DirichletBC(self._adjoint_solver.W.sub(0),Constant(('0.0','0.0')),lambda x, on_boundary : near(x[0],0))
		symmetry = DirichletBC(self._adjoint_solver.W.sub(0).sub(1),Expression('0.0'),lambda x, on_boundary : near(x[1],1))
		bc_adj = [noslip,symmetry,inflow]
		# ??? vorrei fare una cosa come nella riga seguente
		# bc_adj = [copy.deepcopy(item).homogenize() for item in self._solver.bc]
		# bc_adj = []
		# bc_adj = self._solver.bc[:]
		# for item in bc_adj:
			# item.homogenize()
		#g_adj = self._p*n-self._solver.g
		g_adj = Constant(('0.0','0.0'))
		self._adjoint_solver.set_data(self._solver.eta,self._nu,f_adj,bc_adj,g_adj,self._solver.boundaries,True)
			# !!! True serve per scaricare la derivata del laplaciano di u
		
		# @@@ adjoint state problem assembling 
		if MYVERBOSE == True:
			print "\n@@@ adjoint state problem assembling"
		self._adjoint_solver.assemble(q)
		
		# @@@ adjoint state solution computing
		if MYVERBOSE == True:
			print "\n@@@ adjoint state solution computing"
		self._z,self._s = self._adjoint_solver.solve().split()

	# DEBUG	
		plot(self._u,title="u in grad",key="u_grad")
		plot(self._p,title="p in grad",key="p_grad")
		plot(self._z,title="z in grad",key="z_grad")
		plot(self._s,title="s in grad",key="s_grad")
		#interactive()
	# END DEBUG

		# @@@ Gradient assembling
		if MYVERBOSE == True:
			print "\n@@@ Gradient assembling"
		DTq = grad(self._Tq)
		DTq_inv = inv(DTq)
		DTq_invT = transpose(DTq_inv)
		D4 = lambda u :  Dx(Dx(Dx(Dx(u,0),0),0),0) # convenient for code readability

		# Trace operator on terms depending on u,z,s
		Q = q.function_space()
		P = FunctionSpace(self._solver.mesh(),'CG',1)
		#P = self._solver.W.sub(0).sub(0)
			# !!! non si possono creare Function da sottospazi
		dq_uzs = TrialFunction(P)
		phi = TestFunction(P)

		#b# boundaries = FacetFunction("uint",mesh)
		#b# boundary = MovingBoundary()
		#b# boundaries.set_all(1)
		#b# boundary.mark(boundaries,0)
		ds = Measure('ds')[self._solver.boundaries]

#		# trace_form = dq_uzs * phi * ds(0)
#		trace_form = dot(grad(dq_uzs),grad(phi))*dx + Constant('1000000')*dq_uzs * phi * ds(0)
#			# ??? ha senso questa sorta di rilevamento?
#		trace_rhs = - Constant('1000000')*0.5*inner(grad(self._u)*DTq,grad(self._u)*DTq) * phi * ds(0) \
#				- Constant('1000000')*dot(grad(self._u)*DTq*DTq_invT*n, \
#				 (grad(self._u)-self._nu*grad(self._z))*DTq*DTq_invT*n) * phi * ds(0) # + self._s*DTq_invT*n) * phi * ds(0)
#		trace_bc = DirichletBC(P,Constant('0.0'),lambda x,on_boundary : on_boundary and ( near(x[0],0.0) or near(x[0],1.0) ))
		# Dal Manzoni, p.75
		trace_form = dot(grad(dq_uzs),grad(phi))*dx + dq_uzs*phi*dx
		trace_rhs = 0.5*inner(grad(self._u)*DTq_inv,grad(self._u)*DTq_inv) * phi * ds(0) \
				+ dot(grad(self._u)*DTq_inv*DTq_invT*n, \
				 (grad(self._u)-self._nu*grad(self._z))*DTq_inv*DTq_invT*n) * phi * ds(0) # + self._s*DTq_invT*n) * phi * ds(0)
		trace_bc = DirichletBC(P,Constant('0.0'),lambda x,on_boundary : on_boundary and ( near(x[0],0.0) or near(x[0],1.0) or near(x[1],1.0) ))
		value_from_state_and_adjoint = Function(P)
		trace_A = assemble(trace_form)
		trace_A.ident_zeros()
		trace_b = assemble(trace_rhs)
		trace_bc.apply(trace_A,trace_b)

#		for_value_from_state_and_adjoint = -0.5*inner(grad(self._u)*DTq,grad(self._u)*DTq) \
#					+ dot(grad(self._u)*DTq*DTq_invT*n, (grad(self._u)+self._nu*grad(self._z))*DTq*DTq_invT*n - self._s*n)
#		D = lambda u : Dx(u,0)
#		for_restr_form = dot(D(dq_uzs),D(phi)) * dx
#		for_restr_rhs = Constant('0.0')*phi * dx
#		for_restr_bc = DirichletBC(P,for_value_from_state_and_adjoint,lambda x,on_boundary : on_boundary)
#		solve(for_restr_form==for_restr_rhs, value_from_state_and_adjoint, for_restr_bc)
	# DEBUG
	#	plot(self._Tq,mesh=self._solver.mesh(),title="Tq(x,y) in grad")
	#	plot(self._Tq-Expression(('x[0]','x[1]')),mesh=self._solver.mesh(),title="Tq(x,y) - I in grad")
	# END DEBUG

		# @ gradient trace problem solving
		print "\n@ gradient trace problem solving"
		# solve(trace_form==trace_rhs, value_from_state_and_adjoint)#, trace_bc)
		solve(trace_A, value_from_state_and_adjoint.vector(), trace_b)#, trace_bc)
	# DEBUG
		# if MYDEBUG == True:
		pfig = plot(value_from_state_and_adjoint,mesh=self._solver.mesh(),title="value_from_state_and_adjoint",key="adjoint")
		# pfig.write_png("adjoint_only_border.png")
	# END DEBUG
		# func_value_from_state_and_adjoint = MyRestrictedFunction(to_be_restricted=value_from_state_and_adjoint,mesh_in=q.function_space().mesh())#,V_in=self._solver.P)
		restrictor = RestrictFunction(self._solver.mesh())
		func_value_from_state_and_adjoint = restrictor(value_from_state_and_adjoint)
			# ! Using myutilities.MyRestrictedFunction to return a Function of the same type as the control q
		pfig = plot(func_value_from_state_and_adjoint,title="func_value_from_state_and_adjoint",key="func_adjoint")
		#interactive()
		#func_value_from_state_and_adjoint = Function(q.function_space())
		#for i in range(1,len(func_value_from_state_and_adjoint.vector())):
		#	func_value_from_state_and_adjoint.vector()[i] = value_from_state_and_adjoint.vector()[i]
		if MYVERBOSE == True:
			if q.function_space().mesh().geometry().dim() == 1:
				print "gradiente (senza penalizzazione) in 0.7", func_value_from_state_and_adjoint(0.7)
			else:
				print "gradiente (senza penalizzazione) in 0.7", func_value_from_state_and_adjoint(0.7,0.0)
		# value_from_state_and_adjoint = - 0.5*inner(grad(self._u)*DTq,grad(self._u)*DTq) - dot(grad(self._u)*DTq*DTq_invT*n, (grad(self._u)-self._nu*grad(self._z))*DTq*DTq_invT*n + self._s*n)
		# func_value_from_state_and_adjoint = MyRestrictedFunction(Q,to_be_restricted=project(value_from_state_and_adjoint,self._solver.V))
		
#@#@#@		Omega = P.mesh()
#@#@#@		bOmega = BoundaryMesh(Omega,"exterior")
#@#@#@		mfunc = MeshFunction("uint",bOmega,1)
#@#@#@		#ffunc = FacetFunction("uint",bOmega)
#@#@#@		mfunc.set_all(1)
#@#@#@		class MyGamma0(SubDomain):
#@#@#@			def inside(self,x,on_boundary):
#@#@#@				return abs(x[1]-0)<DOLFIN_EPS#near(x[1],0)
#@#@#@		Gamma_moving = MyGamma0()
#@#@#@		Gamma_moving.mark(mfunc,0)
#@#@#@		I = SubMesh(bOmega,mfunc,0)
#@#@#@		Q = FunctionSpace(I,'CG',1)
#@#@#@		parent_vertex = I.data().mesh_function("parent_vertex_indices")
#@#@#@		func_value_from_state_and_adjoint = Function(Q)
#@#@#@		for i in range(0,I.size(0)):
#@#@#@			func_value_from_state_and_adjoint.vector()[i] = \
#@#@#@				value_from_state_and_adjoint.vector()[P.dofmap().dof_to_vertex_map(Omega)[\
#@#@#@					bOmega.entity_map(0)[parent_vertex[Q.dofmap().vertex_to_dof_map(I)[i]]]]]

		func_value = func_value_from_state_and_adjoint \
				+ self._penalty_coeffs[0]*Constant(assemble(q * dx)-self._Vbar) \
				+ self._penalty_coeffs[1]*project(D4(q),Q)
				# !!! qui c'e' una derivata quarta di una funzione EF...
		# return func_value
		
		# @ gradient projecting
		if MYVERBOSE == True:
			print "\n@ gradient projecting"
		value = project(func_value,q.function_space())
		return value
	
	def last_state(self):
		return self._u, self._p
	
	def last_adjoint_state(self):
		return self._z, self._s
	
	def last_variables(self):
		return q,u,p,z,s
		
#=====================


class CostFunction(Expression):

	def eval(self,value,x):
		# value[0] = 0.5-x[0]**2-x[1]**2
		# value[0] = 0.25-x[0]**2
		value[0] = sin(x[0])*cos(x[1])

	def gradient(self,x):
		# return numpy.array([2*x[0],2*x[1]])
		# return 2*x[0]
		return numpy.array([cos(x[0])*cos(x[1]), -sin(x[0])*sin(x[1])])
		# ??? sarebbe meglio se fosse un oggetto derivato da gradient?
