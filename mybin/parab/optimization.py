""" Classes implementing optimization methods.
    ??? Per ora GradientMethod eredita da OptimizationMethod, ma essendo una policy forse sara' meglio cambiare il legame tra le due classi
"""

from dolfin import *
from math import *
from functions import *
from numpy import *
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
		return (errornorm(self._xk,self._xold,'L2') < self._tol)
	def _stop_normf(self):
		# return (sum(abs(self._fxk-self._fxold)) < self._tol)
		return (errornorm(self._fxk,self._fxold,'L2') < self._tol)

	_criterion = _stop_normx
	
	def _stop_criterion(self):
		""" Returns True if stop is needed. False otherwise.
			It use one of the above stop criteria.
		"""
		iter_check = (self._iter >= self._maxit)
		if iter_check:
			print "@@@ WARNING: Maximum number of optimization iterations exceeded"
			return True
		if self._alphak <= self._alphamin:
			print "@@@ WARNING: alpha_k reached alpha_min"
			return True
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
		self._alphamin = 1.e-5
		
	def apply(self,f,x0,alpha0=1.0):
		self._xk = x0
		self._fxk = f(x0)
		print("Initial functional value f(x0) = {0}".format(self._fxk))
		self._fxks = [f(x0)]
		self._alphak = alpha0
		self._gradfxk = f.gradient(x0)
		self._errx = [0.0]
		self._errf = [0.0]
		# self._xold = x0
		# self._fxold = self._fxk
		self._iter = 0
		file_xk = File("xk.pvd")
		file_xk << (x0,0.0)
		file_gradk = File("gradk.pvd")
		file_gradk << (self._gradfxk,0.0)
		print "-------------------------------------------------------"
		while self._iter==0 or (not self._stop_criterion()):
			# !!! At first iteration, self._stop_criterion is not called
			self._iter += 1
			print("iter = {0}".format(self._iter))
			self._xold = self._xk				# shallow copy, but it's fine
			self._fxold = self._fxk				# shallow copy, but it's fine
			self._gradfxold = self._gradfxk
			print("xold = {0} \t f(xold) = {1}".format(self._xold,self._fxold))
			print("f'(xold) = {0}".format((self._gradfxold)))
			self._alphak = alpha0
			self._gradfxk = f.gradient(self._xk)
			if MYVERBOSE == True:
				if self._gradfxk.function_space().mesh().geometry().dim() == 1:
					print "======== uscita gradiente in 0.7", self._gradfxk(0.7)
				else:
					print "======== uscita gradiente in 0.7", self._gradfxk(0.7,0.0)
			plot(self._gradfxk,mesh=self._xold.function_space().mesh(),title="gradfxk",key="gradfxk")
			file = File("gradfxk.pvd")
			file_gradk << (self._gradfxk,float(self._iter))
			# self._radfxk = self._xold
			self._xk = project(self._xold - (self._alphak*self._gradfxk), self._xold.function_space())
				# !!! altrimenti non va eval() nella prossima riga
			self._fxk = f(self._xk)
			print "@@@ Applying back tracking on descent pace"
			while (self._fxk > self._fxold) and (self._alphak > self._alphamin) :
				self._alphak = 0.5*self._alphak
				self._xk = project(self._xold - (self._alphak*self._gradfxk), self._xold.function_space())
				self._fxk = f(self._xk)
				plot(self._xk,title="xk on backtracking",key="backtrack")
				print(" alpha = {0} ,   f(xk) = {1}".format(self._alphak,self._fxk))
			print "  @ End of back tracking"
			self._fxks.append(self._fxk)
			print("xk = {0} \t f(xk) = {1}".format(self._xk,self._fxk))
			print("f'(xk) = {0}".format(self._gradfxk))
			plot(self._xk,mesh=self._xold.function_space().mesh(),title="xk",key="xk")
			file_xk << (self._xk,float(self._iter))
			# print("f'(xk) = {0} \t ||f'(xk)|| = {1}".format(self._gradfxk,norm(self._gradfxk,'L2')))
			Vincr = assemble(self._xk*dx)-f._Vbar
			if f._Vbar > DOLFIN_EPS :
				print("Volume increase = {0} (relative increase = {1} %)".format(Vincr,100.0*Vincr/f._Vbar))
			else: 
				print("Volume increase = {0}".format(Vincr))
	# DEBUG	
			if MYDEBUG == True:
				# plot(self._xk,mesh=UnitIntervalMesh(20),title="qk(x)")
				print [self._xk.vector()[i] for i in range(0,len(self._xk.vector()))]
				print [self._xold.vector()[i] for i in range(0,len(self._xold.vector()))]
				print numpy.abs(self._xk.vector()-self._xold.vector())
				raw_input("Press ENTER to continue")
	# END DEBUG
			self._errx.append(errornorm(self._xk,self._xold,'L2'))
			self._errf.append(abs(self._fxk-self._fxold))
			print "-------------------------------------------------------"
			#raw_input("Press ENTER to continue\n")
			print "Press ENTER to continue\n"

		fun_tmp = Function(FunctionSpace(UnitIntervalMesh(self._iter+1),'DG',0))
		fun_tmp.vector()[:] = array(self._fxks)
		file = File("fxks.pvd")
		file << fun_tmp
		fun_tmp.vector()[:] = array(self._errx)
		file = File("errx.pvd")
		file << fun_tmp
		fun_tmp.vector()[:] = array(self._errf)
		file = File("errf.pvd")
		file << fun_tmp

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
