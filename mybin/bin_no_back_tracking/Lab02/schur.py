#!/usr/bin/env python

from dolfin import *
from geometry import *

parameters['linear_algebra_backend'] = 'PETSc'

class IterativeInverse(LinearOperator):
  # constructor
  def __init__(self, vector, matrix, solver, preconditioner=None):
    LinearOperator.__init__(self, vector, vector)
    self.matrix = matrix
    self.solver = solver
    self.preconditioner = preconditioner

    if preconditioner is None:
      self.solver.set_operator(self.matrix)
    else:
      self.solver.set_operators(self.matrix, self.preconditioner)

  def size(self, dim):
    return self.matrix.size(dim)

  # matrix-vector product
  def mult(self, src, dst):
    self.solver.solve(dst, src)

class SchurComplement(LinearOperator):
  # constructor
  def __init__(self, vector, Ainv, B):
    LinearOperator.__init__(self, vector, vector)
    self.Ainv = Ainv
    self.B = B

    self.tmp1 = Vector(Ainv.size(1));
    self.tmp2 = Vector(Ainv.size(1));

  # matrix-vector product
  def mult(self, src, dst):
    self.B.transpmult(src, self.tmp1)
    self.Ainv.mult(self.tmp1, self.tmp2)
    self.B.mult(self.tmp2, dst)

class StokesSolver(object):
  # constructor
  def __init__(self, Vh, Qh):
    self.Vh = Vh
    self.Qh = Qh

  # assemble matrix
  def assemble(self, mu, f, u0, bc, boundaries):
    # Dirichlet bc
    u0_boundary = lambda x, on_boundary : on_boundary
    # bc = DirichletBC(self.Vh, u0, u0_boundary)

    u = TrialFunction(self.Vh)
    p = TrialFunction(self.Qh)

    v = TestFunction(self.Vh)
    q = TestFunction(self.Qh)

    ds = Measure('ds')[boundaries]

    # assemble A matrix
    eps = lambda u: sym(grad(u))
    a = 2 * mu * inner(eps(u), eps(v)) * dx
    L = dot(f, v) * dx
    self.A, self.b = assemble_system(a, L, bc)

    # assemble A preconditioner
    a = mu * inner(grad(u),grad(v)) * dx
    self.Ap = assemble_system(a, L, bc)[0]

    # assemble B matrix
    a = - div(u) * q * dx
    self.B = assemble(a)

    # assemble Schur preconditioner
    a = p * q / mu * dx
    self.M = assemble(a)

  def solve(self):
    u = Function(self.Vh)
    p = Function(self.Qh)

    #solver1 = KrylovSolver('gmres', 'none')
    #Ainv = IterativeInverse(u.vector(), self.A, solver1, self.Ap)
    solver1 = LUSolver()
    Ainv = IterativeInverse(u.vector(), self.A, solver1)

    S = SchurComplement(p.vector(), Ainv, self.B)
    #solver2 = KrylovSolver('gmres', 'ilu')
    solver2 = LUSolver()
    solver2.set_operator(S)#, self.M)

    tmp1 = Vector(u.vector().size())
    tmp2 = Vector(p.vector().size())

    Ainv.mult(self.b, tmp1)
    self.B.mult(tmp1, tmp2)

    solver2.solve(p.vector(), tmp2)

    self.B.transpmult(p.vector(), tmp1)
    tmp1 -= self.b
    tmp1 *= -1

    Ainv.mult(tmp1, u.vector())

    return u,p


N = 10
mesh = UnitSquareMesh(N,N,'crossed')

Vh = VectorFunctionSpace(mesh, 'CG', 2)
Qh = FunctionSpace(mesh, 'CG', 1)

stokes = StokesSolver(Vh, Qh)

# viscosity and force
mu = Constant(1)
f = Expression(('0.5 * sin(2*x[0]) - 2 * cos(x[0]) * sin(x[1])',
                '0.5 * sin(2*x[1]) + 2 * sin(x[0]) * cos(x[1])'))

# solution
u0 = Expression(('- cos(x[0]) * sin(x[1])',
                 'sin(x[0]) * cos(x[1])'))
p0 = Expression('0.25 * (cos(2*x[0]) + cos(2*x[1]))')

# Dirichlet bc
#u_in = Expression(('sqrt(x[1])','0.0'))       # inflow velocity
u_in = Expression(('100*2*x[1]*(1-x[1]/2.0)','0.0'))
Gamma_moving = MovingBoundary()
noslip = DirichletBC(Vh,Constant(('0.0','0.0')),lambda x, on_boundary : Gamma_moving.inside(x,on_boundary))
inflow = DirichletBC(Vh,u_in,lambda x, on_boundary : (not Gamma_moving.inside(x,on_boundary)) and near(x[0],0))
symmetry = DirichletBC(Vh.sub(1),Expression('0.0'),lambda x, on_boundary : (not Gamma_moving.inside(x,on_boundary)) and near(x[1],1))
bc = [noslip,inflow,symmetry]
# Neumann bc
g = Expression(('0.0','0.0'))# boundary conditions
class NeuBoundary(SubDomain):
	def inside(self,x,on_boundary):
		return on_boundary and near(x[0],1.0)
			# ??? c'e' un modo per far chiamare DirBoundary.inside senza creare un'istanza?
			# perche' altrimenti dovrei fare l'override di qualche metodo affinche' NeuBoundary abbia un oggetto DirBoundary...
neu_boundary = NeuBoundary()
# FacetFunction and ds settings
boundaries = FacetFunction("uint",mesh)
boundaries.set_all(1)
neu_boundary.mark(boundaries,10)
Gamma_moving.mark(boundaries,0)
ds = Measure('ds')[boundaries]


stokes.assemble(mu, f, u0, bc, boundaries)
u,p = stokes.solve()

plot(u)
plot(p)
interactive()

