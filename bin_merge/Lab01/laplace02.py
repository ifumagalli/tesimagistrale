#!/usr/bin/env python

from dolfin import *

# mesh generation
N = 10
mesh = UnitSquareMesh(N,N,'crossed')

# finite element space
V = FunctionSpace(mesh, 'CG', 1)

# boundaries
boundaries = FacetFunction('size_t', mesh)
boundaries.set_all(0)

class Bottom(SubDomain):
  def inside(self, x, on_boundary):
    return near(x[1],0)
class Right(SubDomain):
  def inside(self, x, on_boundary):
    return near(x[0],1)
class Top(SubDomain):
  def inside(self, x, on_boundary):
    return near(x[1],1)
class Left(SubDomain):
  def inside(self, x, on_boundary):
    return near(x[0],0)

Bottom().mark(boundaries, 1)
Right().mark(boundaries, 2)
Top().mark(boundaries, 3)
Left().mark(boundaries, 4)

ds = Measure('ds')[boundaries]

# Dirichlet bcs
bcs = [ DirichletBC(V, Constant(0),                 boundaries, 2),
        DirichletBC(V, Expression('x[1]*(1-x[1])'), boundaries, 4) ]

# variational formulation
u = TrialFunction(V)
v = TestFunction(V)

f = Constant(0)
g = Constant(1)

# https://bitbucket.org/fenics-project/dolfin/issue/78/
a = dot(grad(u),grad(v))*dx + f*u*v*ds(1)
L = g*v*ds(1)

# assemble
A,b = assemble_system(a,L,bcs)

# solve
uh = Function(V)

solver = LUSolver()
solver.set_operator(A)
solver.solve(uh.vector(), b)

# save
File('solution.pvd') << uh

# plot
plot(uh)
interactive()
