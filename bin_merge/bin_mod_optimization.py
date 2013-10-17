""" Classes implementing optimization methods.
    NB ??? Per ora GradientMethod eredita da OptimizationMethod, ma essendo una policy forse sara' meglio cambiare il legame tra le due classi
"""

from math import *
from functions import *
import copy

#=====================

class OptimizationMethod(object):
	""" Abstract class defining an iterative method for optimize a "function" f(x)
		starting from x0, with a tolerance and a maximum number of iterations
	"""
	def __init__(self,tol,maxit):
		# self._f = f
		# self._x0 = x0
		self._tol = tol
		self._maxit = maxit
		# serviranno altre tolleranze in altri oggetti
		# self._iter = 0
		# self._xk = x0
		# self._fxk = self._f(x0)
		# self._xold = None
		# self._fxold = None

	def apply(self,f,x0):
		pass
	
	# stop criteria
	def stop_normx(self):
		iter_check = (self._iter > self._maxit)
		if iter_check:
			print "  ***MAXIMUM NUMBER OF ITERATIONS EXCEEDED***"
			print "  ***in OptimizationMethod.stop_normx()***"
			return iter_check
		# return iter_check or (sum(abs(self._xk-self._xold)) < self._tol)
		return iter_check or (abs(self._xk-self._xold) < self._tol)
			# ??? come fare per tenere generalizzazione al caso di x vettoriale? Perche' gli scalari non hanno un metodo length
		# returns True if stop is needed
	
	_stop_criterion = stop_normx

	# TODO
	# ...set,get,info...

#=====================
	
class GradientMethod(OptimizationMethod):
	""" The "function" f(x) to be optimized need to be a class with a method derivative(x) that returns the derivative of f applied in x.
		With 'derivative' we mean the gradient or, in general, the Riesz element of the (Frechet-)differential
	"""
	def __init__(self,tol,maxit):
		OptimizationMethod.__init__(self,tol,maxit)
		# self._alphak = _alpha0
		# self._gradf = gradf		# TODO il gradiente sara' un metodo dell'oggetto f

	def apply(self,f,x0,alpha0=1.0):
		self._xk = x0
		self._fxk = f(x0)
		#f.eval(self._fxk,x0)
		self._alphak = alpha0
		self._xold = 0*x0			# ??? metto 0*stesso_oggetto affinche' funzioni il criterio di arresto, ma preferirei un bel None
		self._fxold = 0*self._fxk	# idem
		self._iter = 1
		# while not self._stop_criterion():
		while conditional(Not(self._stop_criterion()),True,False):
			# ??? non c'e' un modo per far si' che _stop_criterion restituisca un bool di Python?
			# perche' usando la riga sopra commentata da' l'errore seguente: 
			# " UFL conditions cannot be evaluated as bool in a Python context "
			print "iter =",self._iter
			print "xk =",self._xk,"xold = ",self._xold
			print "f(xk) =",self._fxk,"f(xold) =",self._fxold
			print "f'(xk) =",f.gradient(self._xk),"f'(xold) =",f.gradient(self._xold)
			# print abs(self._xk-self._xold),self._stop_criterion()
			# print "-------"
			# self._xold = copy.deepcopy(self._xk)	# ??? DEEPCOPY
			self._xold = self._xk					# shallow copy, but it's fine
			self._alphak = self._update_alphak()
			self._xk = self._xk - (self._alphak*f.gradient(self._xk))	# DO NOT use -= because so far _xold and _xk share the same object
			self._fxold = self._fxk					# shallow copy, but it's fine
			self._fxk = f(self._xk)
			self._iter += 1
			# print self._xk,self._xold,abs(self._xk-self._xold),self._stop_criterion()
		# DEBUG	
			# plot(self._xk,mesh=UnitIntervalMesh(20),title="qk(x)")
			print [self._xk.vector()[i] for i in range(0,len(self._xk.vector()))]
			raw_input("Press ENTER to continue")
			print "-------"
		return self._xk
	
	# TODO updating techniques
	# def linear_search(self,...):
	#	...
	#	return alphak
	# per il momento non aggiorniamo, con il seguente
	def _update_alphak(self):
		return self._alphak
	
	# ...set,get,info...


#############################################################################
# class OptimizationAlgorithm:
	# def __init__(self,optTol,stateSolver,optMethod,linearSystemTol=1.e-3*optTol):
		# self._optTol = optTol
			# # serviranno anche altre tolleranze come linearSystemTol
		# self._stateSolver = stateSolver
		# self._optMethod = optMethod
			# # Has to be an OptimizationMethod
		
	# def setMethod(self,optMethod):
		# self._optMethod = optMethod
	# def getMethod(self):
		# return self._optMethod
	# def infoMethod(self):
		# return self._optMethod.info
	# # ...
	# # analoghi metodi per gli altri attributi della classe
