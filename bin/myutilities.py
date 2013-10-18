""" In this file we collect some hand-made utilities for the program
"""

from dolfin import *

class MyRestrictedFunction(Expression):
	""" Class to allow restriction of a Function to the bottom boundary
	"""
	# DEVO USARE Expression PERCHE' SE NO NON POSSO SOVRASCRIVERE eval
	# def __init__(self,*args,**kwargs):
		# Function.__init__(self,*args)#,**kwargs)
	# def __init__(self,func_space,**kwargs):
	def __init__(self,*args,**kwargs):
		# # first argument MUST be passable to a Function
		# self._function = Function.__init__(self,args[0])
		self._to_be_restricted = kwargs.get('to_be_restricted',None)
		if self._to_be_restricted == None:
			exit("@@@ Missing the function to restrict")
		
	def eval(self,value,x):
		value = self._to_be_restricted((x,0.0))
