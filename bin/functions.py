""" Functions to be optimized by an OptimizationMethod
	NB: classes beginning with MY need to be reviewed in Fenics environment
"""

class MYGenericFunction(object):
	def __init__(self,f,df):
		self._compute = f
		self._derivative = df
	
	def __call__(self,x):
		return f(x)


class CostFunctional(Form):
	""" Reduced cost functional defined on control q in Q
	"""
	# inheriting from Form class, applied derivative can be computed Fenics-way
	def __init__(Q,state_solver,lambda1=1,volume,lambda2=1):
		self._Q = Q # controls set
		self._state_solver = state_solver
		self._form = Expression(form,p1=???)
		self._lambda1 = lambda1 # penalty coefficient for volume
		self._volume = volume	# fixed volume to be penalized
		self._lambda2 = lambda2	# penalty coefficient for q''
		TODO
	
	def gradient(q):
		""" Given a control q, it returns the Riesz representative of the Frechet derivative computed in q
		"""
		# TODO: per il momento lo scriviamo esplicitamente rispetto al nostro specifico funzionale
		int_q = q*dx
		from_penalty = self._lambda1*(assemble(int_q)-self._volume) + self._lambda2*grad(grad(grad(grad(q))))
		from_moving_boundary = 0.5*inner(grad(u)*DTq,grad(u)*DTq)
		return from_penalty + from_moving_boundary + from_state
	

class SolutionOperator(MYGenericFunction):
	# TODO: prima fai il funzionale, ché già c'è un solutore in schur.py
	def __init__(self,V,P,solver): #altri argomenti possono essere S,dS...
		TODO
		self._solver = solver
		self._V = V		# velocity space
		self._P = P		# pressure space
	
	def __call__(q)
		# Takes a control function q and returns the solution sol of the state problem
		return solver(