""" In this file we collect some hand-made utilities for the program
"""

from dolfin import *

# ===== "macro" for debugging =====
MYDEBUG = True

# ===== Class to restrict a Function =====
class MyRestrictedFunction(Expression):
	""" Class to allow restriction of a Function to the bottom boundary
	"""
	# DEVO USARE Expression PERCHE' CON Function NON POSSO SOVRASCRIVERE eval
	# TODO cosi' com'e' mi sa di espediente sporco
	def __init__(self,*args,**kwargs):
		self._to_be_restricted = kwargs.get('to_be_restricted',None)
		if self._to_be_restricted == None:
			exit("@@@ Missing the function to restrict")
		
	def eval(self,value,x):
		value = self._to_be_restricted((x,0.0))
