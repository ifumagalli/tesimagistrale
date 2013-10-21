""" Functions to be optimized by an OptimizationMethod.
	TODO nel costruttore c'e' troppa roba: snellire
	TODO rivedere se c'e' un modo per snellire __call__ e gradient
"""

from dolfin import *
import numpy
from geometry import *
from myutilities import *

class CostFunctional(object):
	def __init__(self,nu,solver,adjoint_solver=None,penalty_coeffs=[1.0,1.0],**kwargs):
		self.nu = nu
		self._solver = solver
		if (adjoint_solver == None):
			print "@@@ adjoint_solver is None"
			self._adjoint_solver = solver
		else:
			self._adjoint_solver = adjoint_solver
		self._u = Function(solver.V)	# state velocity
		self._p = Function(solver.P)	# state pressure
		self._z = Function(solver.V)	# adjoint state velocity
		self._s = Function(solver.P)	# adjoint state pressure
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
		self._Tq = MapTq(element=self._solver.V.ufl_element(),q=q0)
		
		
	def __call__(self,q):
		""" Given a control Function q, it returns the functional evaluated in q.
		"""
		# self._solver.assemble(self.nu,self.f,self.bc,q)
		self._solver.assemble(q)
		self._u,self._p = self._solver.solve()
	# DEBUG
		if MYDEBUG == True:
			u_x,u_y = self._u.split()
			plot(self._u,title=("Velocity vector u"))
			plot(u_x,title=("Velocity x component"))
			plot(u_y,title=("Velocity y component"))
			plot(self._p,title=("Pressure"))
			interactive()
	# END DEBUG

		self._z = None			# to evaluate the functional we don't need the adjoint state
		self._s = None
		# Functional assembling
		self._Tq.q = q
		DTq = grad(self._Tq)
		gammaq = det(DTq)
		Aq = gammaq * DTq * transpose(inv(DTq))
		D = lambda u : Dx(u,0)	# convenient for code readability
		value = 0.5*( \
				assemble(inner(grad(self._u)*Aq,grad(self._u)) * dx) \
				+ self._penalty_coeffs[0]*(assemble(q * dx) - self._Vbar)**2 \
#				+ self._penalty_coeffs[1]*norm(grad(q),'H10')**2 
				)
		return value

	def gradient(self,q):
		""" Returns the gradient of the functional, intended as the
			L2-Riesz representative of the functional Fre'chet derivative computed in q.
		"""
		# self._solver.assemble(self.nu,self.f,self.bc,q)
		self._solver.assemble(q)
		self._u,self._p = self._solver.solve()
		# self._adjoint_solver.assemble(self.nu,self.f,self.bc,q)
		self._adjoint_solver.assemble(q)
		self._z,self._s = self._adjoint_solver.solve()
		# Gradient assembling
		DTq = grad(self._Tq)
		DTq_invT = transpose(inv(DTq))
		D4 = lambda u :  Dx(Dx(Dx(Dx(u,0),0),0),0) # convenient for code readability
		n = FacetNormal(self._solver.mesh())
		# Trace operator on terms depending from u,z,s
		Q = q.function_space()
		P = self._solver.P
		dq_uzs = TrialFunction(P)
		phi = TestFunction(P)
		trace_form = dq_uzs * phi * ds(0)
		trace_rhs = - 0.5*inner(grad(self._u)*DTq,grad(self._u)*DTq) * phi * ds(0) \
				- dot(grad(self._u)*DTq*DTq_invT*n, (grad(self._u)-self.nu*grad(self._z))*DTq*DTq_invT*n + self._s*n) * phi * ds(0)
		trace_bc = DirichletBC(P,Constant('0.0'),lambda x,on_boundary : on_boundary)
		value_from_state_and_adjoint = Function(P)
		solve(trace_form==trace_rhs, value_from_state_and_adjoint, trace_bc)
		func_value_from_state_and_adjoint = MyRestrictedFunction(Q,to_be_restricted=value_from_state_and_adjoint)
			# ! Using myutilities.MyRestrictedFunction to return a Function of the same type as the control q
		func_value = func_value_from_state_and_adjoint \
				+ self._penalty_coeffs[0]*Constant(assemble(q * dx)-self._Vbar) #\
				# + self._penalty_coeffs[1]*project(D4(q),Q)
				# !!! qui c'e' una derivata quarta di una funzione EF...
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
