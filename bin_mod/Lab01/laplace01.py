#!/usr/bin/env python

from dolfin import *

# mesh generation
N = 10
mesh = UnitSquareMesh(N,N,'crossed')

# finite element space
V = FunctionSpace(mesh, 'CG', 1)

# variational formulation
u = TrialFunction(V)
v = TestFunction(V)

n = 2
f = Expression('(1+2*n*n*pi*pi)*cos(n*pi*x[0])*cos(n*pi*x[1])', n=n)

a = (dot(grad(u),grad(v)) + u*v) * dx
L = f * v * dx

# assemble
A,b = assemble_system(a,L)

# solve
uh = Function(V)

solver = LUSolver()
solver.set_operator(A)
solver.solve(uh.vector(), b)

# error
u0 = Expression('cos(n*pi*x[0])*cos(n*pi*x[1])', n=n)

print '-- L2 error norm:',  errornorm(u0,uh,'l2')
print '-- H1 error norm:',  errornorm(u0,uh,'h1')

# plot
plot(mesh)
plot(uh)
plot(grad(uh))
interactive()

