#!/usr/bin/env python

from dolfin import *
from geometry import *

# parameters['linear_algebra_backend'] = 'PETSc'
parameters['linear_algebra_backend'] = 'uBLAS'

class StokesSolver(object):
  # constructor
  def __init__(self, V, P):
    self.V = V
    self.P = P
    self.W = V * P

  def set_data(self,eta,nu,f,bc,g=None,boundaries=None):
    # TODO forse meglio un generico passaggio per kwargs?
    # si potrebbe tenere un self._kwargs da cui estrarre al momento dell'assemble
    self._eta = eta
    self._nu = nu
    self._f = f
    self._bc = bc
    self._g = g
    self._boundaries = boundaries

  # assemble matrix
  def assemble(self, q):
    if self._eta==None or self._nu==None or self._f==None or self._bc==None:
        exit("@@@ EXITING: data setting is required before assembling")
    if self._g == None and self._boundaries != None:
        print("@@@ WARNING: 'boundaries' FacetFunction passed won't be used")
    if self._g != None and self._boundaries == None:
        exit("@@@ EXITING: Neumann data need a FacetFunction ('boundaries')")
    
    u,p = TrialFunctions(self.W)
    v,r = TestFunctions(self.W)

    # coefficients due to control
    Tq = MapTq(element=self.V.ufl_element(),q=q)
        # ??? meglio fare un membro Tq anche in questa classe? in CostFunctional c'e' gia'
    DTq = grad(Tq)
    gammaq = det(DTq)
    DTq_invT = transpose(inv(DTq))
    Aq = gammaq*DTq*DTq_invT

    # assemble A matrix
    a = self._nu*inner(grad(u)*Aq,grad(v)) * dx \
      - inner(grad(u),DTq_invT)*r*gammaq * dx \
      - inner(grad(v),DTq_invT)*p*gammaq * dx
    L = dot(self._f, v)*gammaq * dx #+ dot(self._g,v) * ds(10)
    self.A, self.b = assemble_system(a, L, self._bc, exterior_facet_domains=self._boundaries)

    # assemble A preconditioner
    a = self._nu * inner(grad(u),grad(v)) * dx \
        + p * r / self._nu * dx
    self.Ap = assemble_system(a, L, self._bc, exterior_facet_domains=self._boundaries)[0]

  def solve(self):
    x = Function(self.W)
    
    #solver = KrylovSolver('gmres', 'ilu')
    #solver.set_operators(self.A, self.Ap)
    solver = LUSolver()
    solver.set_operator(self.A)
    solver.solve(x.vector(), self.b)

    return x.split()
	# Se vuoi che u,p siano modificabili all'esterno senza dover fare deep copies, usa invece le righe seguenti
    # u = Function(self.V)
    # p = Function(self.P)
    # u = assign(x.split()[0])
    # p = assign(x.split()[1])
    # return u,p
    
  def mesh(self):
    return self.V.mesh()

# =====================


# ===== Test code =====

if __name__=="__main__":

  N = 10
  mesh = UnitSquare(N,N,'crossed')
  V = VectorFunctionSpace(mesh, 'CG', 2)
  P = FunctionSpace(mesh, 'CG', 1)

  # problem data
  eta = Constant(0.0)
  nu = Constant(1.e-6)
  f = Expression(('0.0','0.0'))

  # Dirichlet bc
  u_in = Expression(('sqrt(x[1])','0.0'))       # inflow velocity
  Gamma_moving = MovingBoundary()
  noslip = DirichletBC(V,Constant(('0.0','0.0')),lambda x, on_boundary : Gamma_moving.inside(x,on_boundary))
  inflow = DirichletBC(V,u_in,lambda x, on_boundary : (not Gamma_moving.inside(x,on_boundary)) and near(x[0],0))
  symmetry = DirichletBC(V.sub(1),Expression('0.0'),lambda x, on_boundary : (not Gamma_moving.inside(x,on_boundary)) and near(x[1],1))
  bc = [noslip,inflow,symmetry]

  # Neumann bc
  g = Expression(('0.0','0.0'))# boundary conditions
  class NeuBoundary(SubDomain):
    def inside(self,x,on_boundary):
		return on_boundary and near(x[0],1.0)
            # ??? c'e' un modo per far chiamare DirBoundary.inside senza creare un'istanza?
            # perche' altrimenti dovrei fare l'override di qualche metodo affinche' NeuBoundary abbia un oggetto DirBoundary...
  boundaries = FacetFunction("uint",mesh)
  neu_boundary = NeuBoundary()
  neu_boundary.mark(boundaries,10)

  # control
  Q = FunctionSpace(UnitInterval(N),'CG',2)
  q = Function(Q)
  # q0 = Expression('0.0')
  # q0 = Expression('sin(2*pi*x[0])')
  q0 = Expression('0.2-0.2*4*(x[0]-0.5)*(x[0]-0.5)')
  q = project(q0,Q)
  
  # solver initializing
  stokes = StokesSolver(V,P)
  stokes.set_data(eta, nu, f, bc, g, boundaries)
  
  # solution
  stokes.assemble(q)
  u,p = stokes.solve()

  # post-processing
  u_x,u_y = u.split()
  plot(u,title=("Velocity vector u"))
  plot(u_x,title=("Velocity x component"))
  plot(u_y,title=("Velocity y component"))
  plot(p,title=("Pressure"))
  interactive()
  raw_input("Press ENTER to finish")
