#!/usr/bin/env python

from dolfin import *

# mesh
N = 10
#mesh = RectangleMesh(0,0,3,1,3*N,N,'crossed')
mesh = UnitSquare(N,N,'crossed')

# finite element spaces
V = VectorFunctionSpace(mesh, 'CG', 2)
#B = VectorFunctionSpace(mesh, 'B', 3)
Q = FunctionSpace(mesh, 'CG', 1)

W = MixedFunctionSpace([V,Q])

# solution
u0 = Expression(('x[1]','0'))
u_in = Expression(('100*2*x[1]*(1-x[1]/2.0)','0.0'))
p0 = Expression('0')

# boundaries
boundaries = FacetFunction('size_t', mesh)
boundaries.set_all(0)

class Bottom(SubDomain):
  def inside(self, x, on_boundary):
    return near(x[1],0)
class Right(SubDomain):
  def inside(self, x, on_boundary):
    return near(x[0],3)
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

# dirichlet bcs
bcs = [ DirichletBC(W.sub(0), Constant(('0.0','0.0')), boundaries, 1),
        DirichletBC(W.sub(0).sub(1), Constant('0.0'), boundaries, 3),
	DirichletBC(W.sub(0), u_in, boundaries, 4) ]

# variational formulation
u,p = TrialFunctions(W)
v,q = TestFunctions(W)

f = Constant(0)

a = (inner(grad(u),grad(v)) - p * div(v) + q * div(u)) * dx

# assemble
A = assemble(a)

b = Vector(A.size(0))
b.zero()
for bc in bcs:
  bc.apply(A,b)

# solve
xh = Function(W)

solver = LUSolver()#'gmres', 'none')
solver.set_operator(A)
solver.solve(xh.vector(), b)

uh,ph = xh.split()

# plot
plot(uh)
plot(ph)
interactive()
