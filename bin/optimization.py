### Inizio codice Ivan Fumagalli
""" Commenti vari
"""

from math import *
from functions import *

"""class OptimizationAlgorithm:
	def __init__(self,optTol,stateSolver,optMethod,linearSystemTol=1.e-3*optTol):
		self._optTol = optTol
			# serviranno anche altre tolleranze come linearSystemTol
		self._stateSolver = stateSolver
		self._optMethod = optMethod
			# Has to be an OptimizationMethod
		
	def setMethod(self,optMethod):
		self._optMethod = optMethod
	def getMethod(self):
		return self._optMethod
	def infoMethod(self):
		return self._optMethod.info
	# ...
	# analoghi metodi per gli altri attributi della classe

#=====================
"""
class OptimizationMethod(object):
	""" Abstract class defining an iterative method for optimize a "function" f(x)
		starting from x0, with a tolerance and a maximum number of iterations
	"""
	def __init__(self,f,x0,tol,maxit):
		self._f = f
		self._x0 = x0
		self._tol = tol
		self._maxit = maxit
		# serviranno altre tolleranze in altri oggetti
		self._iter = 0
		self._xk = x0
		self._fxk = self._f(x0)
		self._xold = None
		self._fxold = None

	def apply(self):
		pass
	
	def stop_normx(self):
		return (self._iter > self._maxit) or (abs(self.xk-self.xold) < self.tol)
		# returns True if stop is needed
	
	_stop_criterion = stop_normx

	# TODO
	# ...set,get,info...

#=======
	
class GradientMethod(OptimizationMethod):
	""" The "function" f(x) to be optimized need to be a class with a method derivative(x) that returns the derivative of f applied in x.
		With 'derivative' we mean the gradient or, in general, the Riesz element of the (Fréchet-)differential
	"""
	def __init__(self,f,x0,tol,maxit,alpha0):
		OptimizationMethod.__init__(self,f,x0,tol,maxit)
		#?# non si può evitare di chiamarlo?
		self._alphak = _alpha0

	def apply(self):
		while !self._stop_criterion():
			self._xold = self._xk
			self._alphak = self._update_alphak()
			self._xk = self._xk + self._alphak*self._f.gradient(xk)
			self._fxold = self._fxk
			self._fxk = self._f(xk)
			self._iter += 1
		return xk
	
	# TODO updating techniques
	# def linear_search(self,...):
	#	...
	#	return alphak
	# per il momento non aggiorniamo, con il seguente
	def _update_alphak():
		return self._alphak
	
	# ...set,get,info...

#=====================