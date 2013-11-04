""" Per testare l'algoritmo nel caso di funzioni reali:
	modificare anche CostFunction in functions.py
"""

from optimization import *
import numpy

tol = 1.e-3
maxit = 50
optmethod = GradientMethod(tol,maxit)

# x0 = numpy.array([0.5,-0.5])
# x0 = numpy.array([0.25])
x0 = numpy.array([0.5,-1.2])
f = CostFunction()
x_opt = optmethod.apply(f,x0,0.7)

print x_opt