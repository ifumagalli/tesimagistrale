""" Functions to be optimized by an OptimizationMethod
	NB: classes beginning with MY need to be reviewed in Fenics environment

	TODO: far dipendere da q anche i problemi di stato e aggiunto
"""

from dolfin import *
import numpy	# finora serve solo per CostFunction

class CostFunctional(object):
	def __init__(self,nu,f,bc,form,solver,adjoint_solver=None,grad_form=None):
		# ??? per ora tengo dentro anche nu,f,bc che servono per l'assemble del solver
		# TODO: completare
		self.nu = nu
		self.f = f
		self.bc = bc
		self._form = form
		self._solver = solver
		if (adjoint_solver == None):
			print "***adjoint_solver is None***"
			self._adjoint_solver = solver
		else:
			self._adjoint_solver = adjoint_solver
		if (grad_form == None):
			print "***grad_form is None***"
		self._u = Function(solver.Vh)	# state velocity
		self._p = Function(solver.Qh)	# state pressure
		self._z = Function(solver.Vh)	# adjoint state velocity
		self._s = Function(solver.Qh)	# adjoint state pressure
			# !!! per ora ho messo Vh,Qh per non dover cambiare il solutore: poi lo cambieremo per poter mettere V,P
		self._grad_form = grad_form
	
	def __call__(self,q):
		""" Given a control function q, it returns J(q): the functional evaluated in q
		"""
		self._solver.assemble(self.nu,self.f,self.bc,q)
		self._u,self._p = self._solver.solve()
		self._z = None			# to evaluate the functional we don't need the adjoint state
		self._s = None
		value = assemble(self._form)	# TODO controllare che calcoli effettivamente J(q)(u,p,z,s)
		return value

	def gradient(self,q):
		""" Returns the gradient of the functional, intended as the Riesz representative
			of the functional J'(q)(.) expressed as an integral on the interval I
		"""
		# per ora glielo passiamo al costruttore: magari poi glielo si puo' far calcolare
		self._solver.assemble(self.nu,self.f,self.bc,q)
		self._u,self._p = self._solver.solve()
		self._solver.assemble(self.nu,self.f,self.bc,q)
		self._z,self._s = self._adjoint_solver.solve()
		grad_assembled = assemble(self._grad_form)
		return grad_assembled
	
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
