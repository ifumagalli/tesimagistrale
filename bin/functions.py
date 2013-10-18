""" Functions to be optimized by an OptimizationMethod
	NB: classes beginning with MY need to be reviewed in Fenics environment

	TODO: far dipendere da q anche i problemi di stato e aggiunto
"""

from dolfin import *
import numpy	# finora serve solo per CostFunction
from geometry import *
from myutilities import *

class CostFunctional(object):
	def __init__(self,nu,f,bc,form,solver,adjoint_solver=None,grad_form=None,penalty_coeffs=[1,1],**kwargs):
		# ??? per ora tengo dentro anche nu,f,bc che servono per l'assemble del solver
		# q0 is needed to establish Vbar
		# TODO: completare
		self.nu = nu
		self.f = f
		self.bc = bc
		self._form = form
		self._solver = solver
		if (adjoint_solver == None):
			print "@@@ adjoint_solver is None"
			self._adjoint_solver = solver
		else:
			self._adjoint_solver = adjoint_solver
		if (grad_form == None):
			print "@@@ grad_form is None"
		self._u = Function(solver.Vh)	# state velocity
		self._p = Function(solver.Qh)	# state pressure
		self._z = Function(solver.Vh)	# adjoint state velocity
		self._s = Function(solver.Qh)	# adjoint state pressure
			# !!! per ora ho messo Vh,Qh per non dover cambiare il solutore: poi lo cambieremo per poter mettere V,P
		self._grad_form = grad_form
		self._penalty_coeffs = penalty_coeffs
		# Vbar assembling:	a value can be passed, or an Expression q whose volume will be Vbar
		#					if both are passed, the numerical value will be used
		self._coeffspace = FunctionSpace(UnitInterval(100),'CG',10)
			# Space to which Expressions will be projected if a mesh is needed
			# TODO per ora usiamo FunctionSpace(UnitInterval(100),'CG',10)
		self._Vbar = kwargs.get('Vbar',None)
		q0 = kwargs.get('q0',None)
		if self._Vbar==None:
			if q==None:
				exit("@@@ EXITING: Required fixed volume value or initial control to initialize Vbar ")
			else:
				qf = Function(self._coeffspace)
				qf = project(q0,self._coeffspace)
				self._Vbar = assemble(qf * dx)
		
	def __call__(self,q):
		""" Given a control function q, it returns J(q): the functional evaluated in q
		"""
		self._solver.assemble(self.nu,self.f,self.bc,q)
		self._u,self._p = self._solver.solve()
		self._z = None			# to evaluate the functional we don't need the adjoint state
		self._s = None
		# Functional assembling
		# Tq = MapTq(element=self._coeffspace.ufl_element())
		Tq = MapTq(element=self._solver.Vh.ufl_element())
		Tq.q = q
		DTq = grad(Tq)
		gammaq = det(DTq)
		DTq_invT = transpose(inv(DTq))
		Aq = gammaq*DTq*DTq_invT
		D = lambda u : Dx(u,0)	# convenient for code readability
		value = 0.5*( \
				assemble(inner(grad(self._u)*Aq,grad(self._u)) * dx) \
				+ self._penalty_coeffs[0]*(assemble(q * dx) - self._Vbar)**2 \
#				+ self._penalty_coeffs[1]*norm(grad(q),'H10')**2 
				)
		return value

	def gradient(self,q):
		""" Returns the gradient of the functional, intended as the Riesz representative
			of the functional J'(q)(.) expressed as an integral on the interval I
		"""
		self._solver.assemble(self.nu,self.f,self.bc,q)
		self._u,self._p = self._solver.solve()
		self._solver.assemble(self.nu,self.f,self.bc,q)
		self._z,self._s = self._adjoint_solver.solve()
		# Gradient assembling
		# Tq = MapTq(element=self._coeffspace.ufl_element())
		Tq = MapTq(element=self._solver.Vh.ufl_element())
		Tq.q = q
		DTq = grad(Tq)
		DTq_invT = transpose(inv(DTq))
		D4 = lambda u :  Dx(Dx(Dx(Dx(u,0),0),0),0) # convenient for code readability
		n = FacetNormal(self._solver.mesh())
		# value = - 0.5*assemble(inner(grad(self._u)*DTq,grad(self._u)*DTq) * ds(0)) \
				# - assemble(dot(grad(self._u)*DTq*DTq_invT*n, (grad(self._u)-self.nu*grad(self._z))*DTq*DTq_invT*n + self._s*n) * ds(0)) \
				# + self._penalty_coeffs[0]*(assemble(q * dx)-self._Vbar) \
				# + self._penalty_coeffs[1]*assemble(D4(q) * dx)
		# Trace operator on terms depending from u,z,s
		Q = q.function_space()
		P = self._solver.Qh
		dq_uzs = TrialFunction(P)
		phi = TestFunction(P)
		trace_form = dq_uzs * phi * ds(0)
		trace_rhs = - 0.5*inner(grad(self._u)*DTq,grad(self._u)*DTq) * phi * ds(0) \
				- dot(grad(self._u)*DTq*DTq_invT*n, (grad(self._u)-self.nu*grad(self._z))*DTq*DTq_invT*n + self._s*n) * phi * ds(0)
		trace_bc = DirichletBC(P,Constant('0.0'),lambda x,on_boundary : on_boundary)
		value_from_state_and_adjoint = Function(P)
		solve(trace_form==trace_rhs, value_from_state_and_adjoint, trace_bc)
		func_value_from_state_and_adjoint = MyRestrictedFunction(Q,to_be_restricted=value_from_state_and_adjoint)
		func_value = func_value_from_state_and_adjoint \
				+ self._penalty_coeffs[0]*Constant(assemble(q * dx)-self._Vbar) #\
				# + self._penalty_coeffs[1]*project(D4(q),Q)
				# !!! qui c'e' una derivata quarta di una funzione EF...
		# func_value = Function(Q)
		# print("Q.mesh().coordinates():")
		# print Q.mesh().coordinates()
		# boundary_coordinates = numpy.array(zip([item for sublist in Q.mesh().coordinates().tolist() for item in sublist], \
										# [0.0 for i in range(1,len(Q.mesh().coordinates()))]))
		# print("boundary_coordinates:")
		# print boundary_coordinates
		# # print("P.mesh().coordinates():")
		# # print P.mesh().coordinates()
		# func_value.vector()[:] = value(boundary_coordinates)
			# # to return a Function of the same kind as q
		return func_value
	
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

		
#############################################################################
# class MYGenericFunction(object):
	# def __init__(self,f,df):
		# self._compute = f
		# self._derivative = df
	
	# def __call__(self,x):
		# return f(x)

# class CostFunctionalDerivative(NonlinearProblem):
	# """ Reduced cost functional derivative j'(q), computed on control variation deltaq
	# """
	# """	From FEniCS manual, p.98
	# def __init__(self, L, a, bc):
		# NonlinearProblem.__init__(self)
		# self.L = L
		# self.a = a
		# self.bc = bc
	# def F(self, b, x):
		# assemble(self.L, tensor=b)
		# self.bc.apply(b, x)
	# def J(self, A, x):
		# assemble(self.a, tensor=A)
		# self.bc.apply(A)
	# """
	# def __init__(self,form,penalty_terms,bc):
		# NonlinearProblem.__init__(self)
		# self._form = form + penalty_terms	# form representing j'(q) (as a functional on variations)
		# self._bc = bc						# boundary conditions on control variations
		
	# def F(self,assembled_form,deltaq):
		# assemble(self._form,tensor=assembled_form)
		# self._bc.apply(assembled_form,deltaq)
	
	# # def J(self,assembled_jacobian,deltaq):
		# # TODO se serve
		# # per il momento non la facciamo, che' per il metodo del gradiente non credo serva...
		

# class CostFunctional(Form):
	# """ Reduced cost functional defined on control q in Q
	# """
	# # inheriting from Form class, applied derivative can be computed Fenics-way
	# def __init__(Q,state_solver,lambda1=1,volume,lambda2=1):
		# self._Q = Q # controls set
		# self._state_solver = state_solver
		# self._form = Expression(form,p1=???)
		# self._lambda1 = lambda1 # penalty coefficient for volume
		# self._volume = volume	# fixed volume to be penalized
		# self._lambda2 = lambda2	# penalty coefficient for q''
		# TODO
	
	# def gradient(q):
		# """ Given a control q, it returns the Riesz representative of the Frechet derivative computed in q
		# """
		# # TODO: per il momento lo scriviamo esplicitamente rispetto al nostro specifico funzionale
		# int_q = q*dx
		# from_penalty = self._lambda1*(assemble(int_q)-self._volume) + self._lambda2*grad(grad(grad(grad(q))))
		# TODO: DTq = ...(q)
		# from_moving_boundary = 0.5*inner(grad(u)*DTq,grad(u)*DTq)
		
		# return from_penalty + from_moving_boundary + from_state
	

# class SolutionOperator(MYGenericFunction):
	# # TODO: prima fai il funzionale, che' gia' c'e' un solutore in schur.py
	# def __init__(self,V,P,solver): #altri argomenti possono essere S,dS...
		# TODO
		# self._solver = solver
		# self._V = V		# velocity space
		# self._P = P		# pressure space
	
	# def __call__(q)
		# # Takes a control function q and returns the solution sol of the state problem
		# return solver(
