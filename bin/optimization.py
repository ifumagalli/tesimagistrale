""" Classes implementing optimization methods.
    ??? Per ora GradientMethod eredita da OptimizationMethod, ma essendo una policy forse sara' meglio cambiare il legame tra le due classi
"""

from dolfin import *
from math import *
from functions import *
import copy

#=====================

class OptimizationMethod(object):
	""" Abstract class defining an iterative method to optimize a "function" f(x)
		starting from x0, with a tolerance and a maximum number of iterations.
		No hypotheses are set on f or x0.
	"""
	def __init__(self,tol,maxit):
		self._tol = tol
		self._maxit = maxit
		self._xk = None
		self._xold = None
		self._fxk = None
		self._fxold = None

	def apply(self,f,x0):
		pass
	
	# stop criteria
		# le righe commentate sono da usare se x e' un vettore, non una Function
	def _stop_normx(self):
		# return (sum(abs(self._xk-self._xold)) < self._tol)
		return (errornorm(self._xk,self._xold,'l2') < self._tol)
	def _stop_normf(self):
		# return (sum(abs(self._fxk-self._fxold)) < self._tol)
		return (errornorm(self._fxk,self._fxold,'l2') < self._tol)

	_criterion = _stop_normx
	
	def _stop_criterion(self):
		""" Returns True if stop is needed. False otherwise.
			It use one of the above stop criteria.
		"""
		iter_check = (self._iter >= self._maxit)
		if iter_check:
			print "@@@ WARNING: Maximum number of optimization iterations exceeded"
			return iter_check
		return self._criterion()
		# ??? appesantisco troppo a fare una chiamata di funzione in piu', oppure in Python non cambia nulla?
	
	def set_stop_criterion(self,criterion):
		""" To set possible user defined stop criterion.
		"""
		self._criterion = criterion

	# TODO
	# ...set,get,info...

#=====================
	
class GradientMethod(OptimizationMethod):
	""" The "function" f(x) to be optimized need to be a class with a method gradient(x) that returns the gradient of f applied in x.
		By 'gradient', in general we mean the Riesz element of the (Frechet-)differential.
	"""
	def __init__(self,tol,maxit):
		OptimizationMethod.__init__(self,tol,maxit)
		
	def apply(self,f,x0,alpha0=1.0):
		self._xk = x0
		self._fxk = f(x0)
		self._alphak = alpha0
		# self._xold = x0
		# self._fxold = self._fxk
		self._iter = 0
		while self._iter==0 or (not self._stop_criterion()):
			# !!! At first iteration, self._stop_criterion is not called
			self._iter += 1
			print("iter = {0}".format(self._iter))
			self._xold = self._xk				# shallow copy, but it's fine
			self._fxold = self._fxk				# shallow copy, but it's fine
			print("xold = {0} \t f(xold) = {1}".format(self._xold,self._fxold))
			print("f'(xold) = {0}".format(f.gradient(self._xold)))
			self._alphak = self._update_alphak()
			gradfxk = f.gradient(self._xk)
			self._xk = project(self._xk - (self._alphak*gradfxk), self._xold.function_space())
			self._fxk = f(self._xk)
			print("xk = {0} \t f(xk) = {1}".format(self._xk,self._fxk))
			print("f'(xk) = {0}".format(gradfxk))
	# DEBUG	
			if MYDEBUG == True:
				# plot(self._xk,mesh=UnitIntervalMesh(20),title="qk(x)")
				print [self._xk.vector()[i] for i in range(0,len(self._xk.vector()))]
				print [self._xold.vector()[i] for i in range(0,len(self._xold.vector()))]
				print numpy.abs(self._xk.vector()-self._xold.vector())
	# END DEBUG
			raw_input("Press ENTER to continue")
			print "-------"
		return self._xk
	
	# TODO updating techniques
	# def linear_search(self,...):
	#	...
	#	return alphak
	# per il momento non aggiorniamo: usiamo il seguente
	def _update_alphak(self):
		return self._alphak
	
	# TODO
	# ...set,get,info...
