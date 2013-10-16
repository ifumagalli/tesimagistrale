from dolfin import *

N = 5
# due spazi diversi in tutto (ma entrambi CG)
V = VectorFunctionSpace(UnitSquare(N,N),'CG',2)
Q = FunctionSpace(UnitInterval(N),'CG',1)

v = Function(V)
q = Function(Q)

v_expr = Expression(('sin(2*pi*x[0])','5*exp(x[1])+x[0]'))
q_expr = Expression('sin(2*pi*x[0])')

v = project(v_expr,V)
q = project(q_expr,Q)

pv = Point(0.4,0.2)
pq = Point(0.4)
print("v_expr(0.5,0.2) = {0},   q_expr(0.5) = {1}".format(v_expr(pv), q_expr(pq)))

val_q = q(pq)
val_v = v(pv)

print("v(0.5,0.2) = {0},   q(0.5) = {1}".format(val_v, val_q))

print("Test superato!")
