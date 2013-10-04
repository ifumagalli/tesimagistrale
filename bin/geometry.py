""" Classes and functions for geometry managing
"""

class map_Tq(object):
	def __init__(q):
		self.q = q	# control
		TODO

	def __call__(x,y):
		""" Given a point (x,y) in the reference domain Omega_0
			it returns the point Tq(x,y) in the actual domain Omega_q
		"""
		return x, y+(1-y)*self.q(x)
	
	def inverse(x,y):
		""" Given a point (x,y) in the ACTUAL domain Omega_q
			it returns the point Tq^{-1}(x,y) in the REFERENCE domain Omega_0
		"""
		return x, (y-self.q(x))/(1-self.q(x))
	
	def move_back(u):
		""" Given a function u defined on the ACTUAL domain Omega_q
			it returns u composed with Tq^{-1}
		"""
		TODO
		MA FORSE NON E' DA FARE: E' SEMPLICEMENTE UNA COMPOSIZIONE (v. Function in Fenics)