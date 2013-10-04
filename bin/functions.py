""" Functions to be optimized by an OptimizationMethod
	NB: classes beginning with MY need to be reviewed in Fenics environment
"""

class MYGenericFunction(object):
	def __init__(self,f,df):
		self._compute = f
		self._derivative = df
	
	def __call__(self,x):
		return f(x)


class CostFunctional(MYGenericFunction):
	""" Reduced cost functional defined on control q in Q
	"""
	# TODO sarebbe utile derivare da un funzionale fenics
	def __init__(Q,state_solver):
		self._Q = Q # controls set
		self._state_solver = state_solver
		TODO
	
	def derivative(q):
		""" Given a control q, it returns the Riesz representative of the Frechet derivative computed in q
		"""
		
	

class SolutionOperator(MYGenericFunction):
	# TODO: prima fai il funzionale, ché già c'è un solutore in schur.py
	def __init__(self,solver): #altri argomenti possono essere S,dS...
		TODO
		self._solver = solver
	
	def __call__(q)
		# Takes a control function q and returns the solution sol of the state problem
		return solver(