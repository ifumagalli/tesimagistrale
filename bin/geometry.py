""" Classes and functions for geometry managing
	Everything is done for a 2D domain.
"""

from dolfin import Expression

class MapTq(Expression):	# inheriting from class dolfin.functions.expression.Expression
	def __init__(self,q)
		# NB Arguments have to be passed using **kwargs (e.g.   qh=Function(Qh)
		#														Tq = MapTq(q=qh)) 
		#	 and Expression.__init__(...) has not to be called directly: see the link below
		# http://fenicsproject.org/documentation/dolfin/1.2.0/python/programmers-reference/functions/
		#	expression/Expression.html#dolfin.functions.expression.Expression
		self.q = q	# control: a scalar Function on an Interval
		TODO

	def eval(self,values,x):
		""" Given a point x=[x[0] x[1]] in the REFERENCE domain Omega_0
			it returns in 'values' the point Tq(x) in the ACTUAL domain Omega_q
		"""
		values = (x[0], x[1]+(1-x[1])*self.q(x[0]))
	
	def inverse(self,values,x):
		""" Given a point x=[x[0] x[1]] in the ACTUAL domain Omega_q
			it returns in 'values' the point Tq^{-1}(x) in the REFERENCE domain Omega_0
		"""
		values = (x[0], (x[1]-self.q(x[0]))/(1-self.q(x[0])))
	
	def move_back(self,u):
		""" Given a function u defined on the ACTUAL domain Omega_q
			it returns u composed with Tq^{-1}
		"""
		TODO
		MA FORSE NON E' DA FARE: E' SEMPLICEMENTE UNA COMPOSIZIONE (v. Function in Fenics)

#=======
# other objects deriving from MapTq
# ??? vorrei fare l'overloading dell'operatore nabla_grad, ma non so come è fatto, e poi in Python non esiste overloading

# class DTq(LinearOperator):
	# def __init__(self,q)
		# vector = Vector(2)
		# LinearOperator.__init__(self, vector, vector)

#	# matrix-vector product
	# def mult(self, src, dst):
		# self.B.transpmult(src, self.tmp1)
		# self.Ainv.mult(self.tmp1, self.tmp2)
		# self.B.mult(self.tmp2, dst)
	# MEGLIO FARLO DERIVARE DA Expression
#class DTq(Expression):
class DTq(Expr): # from ufl
	def __init__(self,q):
		Expression.__init__(self,2,2)
		self.q = q
	
	def eval(self,values,x):
		values = np.array([[1 , 0] , [(1-x[1])*grad(q)(x[0]) , 1-q(x[0])]])
	
	def inverse(self,values,x):
		values = np.array([[1 , 0] , [-(1-x[1])/(1-q(x[0]))*grad(q)(x[0]) , 1.0/(1-q(x[0]))]])
	
	def transpose