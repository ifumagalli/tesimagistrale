""" Classes and functions for geometry managing.
	Everything is done for a 2D domain.
    NB !!! ho messo q(x[0],0), invece di q(x[0])
"""

from dolfin import *
from numpy import array

# ===== Boundaries =====

class MovingBoundary(SubDomain):
	def __init__(self):
		SubDomain.__init__(self)

	def inside(self,x,on_boundary):
		#return on_boundary and near(x[1],0)
		return near(x[1],0)
			# !!! altrimenti non funziona RestrictFunction

moving_boundary = MovingBoundary()

class StillBoundary(SubDomain):
    def __init__(self,moving_boundary):
        SubDomain.__init__(self)
        self._moving_boundary=moving_boundary
	def inside(self,x,on_boundary):
		return on_boundary and not self._moving_boundary.inside(x,on_boundary)


# ===== MapTq: Omega_0->Omega_q =====

class MapTq(Expression):
	""" Class implementing the map Tq: Omega_0->Omega_q
		which maps the reference domain into the actual one.
		It depends on the control function q, which identifies Omega_q.
		! Cartesian 2D domains are considered
	"""
	def __init__(self,*args,**kwargs):
		# NB Arguments have to be passed using **kwargs (e.g.   qh=Function(Qh)
		#														Tq = MapTq(q=qh)) 
		#	 and Expression.__init__(...) has not to be called directly: see the link below
		# http://fenicsproject.org/documentation/dolfin/1.2.0/python/programmers-reference/functions/
		#	expression/Expression.html#dolfin.functions.expression.Expression
		self.q = kwargs.get('q')
			# control: a scalar Function on an Interval

	def eval(self,values,x):
		""" Given a point x=[x[0] x[1]] in the REFERENCE domain Omega_0
			it returns in 'values' the point Tq(x) in the ACTUAL domain Omega_q.
			As far, Omega_q is similar to Omega_0, with the bottom border substituted by q's graph.
		"""
		values[0] = x[0]
		tmp = self.q(x[0],0)
		values[1] = x[1]+(1-x[1])*tmp
		#values[1] = x[1]
		
	def value_shape(self):
		return (2,)
			# we consider 2D domains
		
#	def inverse(self,value,x):
	def inverse(self,x):
		""" Given a point x=[x[0] x[1]] in the ACTUAL domain Omega_q
			it returns in 'values' the point Tq^{-1}(x) in the REFERENCE domain Omega_0
		"""
		value = (x[0], (x[1]-self.q(x[0],0))/(1-self.q(x[0],0)))
		return value
