from dolfin import *
from geometry import *
from myutilities import *

N = 25
Omega = UnitSquare(N,N)#,"crossed")
V = FunctionSpace(Omega,'CG',1)

r = Function(V)
r_expr = Expression('sin(4*pi*x[0])*cos(6*pi*x[1])')
#r_expr = Expression('1-x[1]')#+25.0*x[1]')
r = project(r_expr,V)
#for i in range(0,len(r.vector())):
#	r.vector()[i] = i
pfig = plot(r,title="Funzione di partenza",axes=True)
pfig.elevate(-15)
# interactive()

# ===================================
# funzione definita su tutto il bordo
# ===================================
bOmega = BoundaryMesh(Omega,"exterior")
bV = FunctionSpace(bOmega,'CG',1)
br = Function(bV)

emap = bOmega.entity_map(0)
#@v# emap = bOmega.vertex_map()
# print emap
# print emap[0],emap[7],emap[19],emap[39],emap[41],emap[56]
print "=== br.vector()"
for i in range(0,bOmega.size(0)):
	br.vector()[i] = r.vector()[emap[i]]
	print br.vector()[i]
print "==="
# #plot(br,title="Funzione sul bordo completo - NO")
# #interactive()
	# # !!! il plot non rispecchia la realta'
# print br.vector().str(True)#[1:len(br.vector())]

## ===================================
## Proviamo con Dirichlet
## ===================================
#Gamma_moving = MovingBoundary()
#bc = DirichletBC(V,r,lambda x,on_boundary : Gamma_moving.inside(x,on_boundary))
#values = bc.get_boundary_values()
#I = UnitInterval(N)
#Q = FunctionSpace(I,'CG',1)
#q = Function(Q)
#print values[1]
#for i in range(0,I.size(0)):
#	q.vector()[i] = values[i]
#plot(q,title="Con DiricletBC")
#print q.vector().str(True)


# ===================================
# funzione definita sul bordo inferiore
# ===================================
Gamma_moving = MovingBoundary()
Gamma_still = StillBoundary(Gamma_moving)
mfunc = MeshFunction("uint",bOmega,1)
ffunc = FacetFunction("uint",bOmega)
#Gamma_still.mark(mfunc,0)
#Gamma_moving.mark(mfunc,1)
#Gamma_still.mark(ffunc,0)
mfunc.set_all(1)
class MyGamma0(SubDomain):
	def inside(self,x,on_boundary):
		return near(x[1],0)
		#return near(x[0],1)# or near(x[1],0)
Gamma_moving = MyGamma0()
Gamma_moving.mark(mfunc,0)
print Omega,bOmega
print "============="
I = SubMesh(bOmega,mfunc,0)
#I = SubMesh(Omega,Gamma_moving)
print I
plot(Omega,title="Omega")
plot(bOmega,title="bOmega")
plot(I,title="I")
interactive()
print "============="
print I.geometry().dim(),I.topology().dim()
#I.clean()
#print I.geometry().dim(),I.topology().dim()
#for i in range(0,I.size(0)):
#	print br.vector()[I.parent_vertex_indices[i]]
#I = UnitInterval(N)
parent_vertex = I.data().mesh_function("parent_vertex_indices")
print "#@#@"
for i in range(0,I.size(0)):
	#print br.vector()[map[parent_vertex[i]]]#[I.data().array("parent_vertex_indices")[i]]
	print br.vector()[parent_vertex[i]]
	#print r.vector()[i]
	#print r.vector()[emap[parent_vertex[i]]]
	#print br(bOmega.coordinates()[i])
#print("r = {0}".format(r.vector().str(True)))
# print I.__dict__
#I = Mesh(I)
#I = UnitInterval(N)
Q = FunctionSpace(I,'CG',1)
q = Function(Q)
# dof_map = bV.dofmap()
# map = bV.dofmap().vertex_to_dof_map(bOmega)
for i in range(0,I.size(0)):
	#q.vector()[i] = br.vector()[map[parent_vertex[i]]]#[I.data().array("parent_vertex_indices")[i]]
	#q.vector()[i] = r.vector()[emap[parent_vertex[i]]]
	#q.vector()[i] = r.vector()[emap[parent_vertex[i]]]
	#q.vector()[i] = br(bOmega.coordinates()[i])
	q.vector()[i] = r.vector()[V.dofmap().dof_to_vertex_map(Omega)[bOmega.entity_map(0)[parent_vertex[Q.dofmap().vertex_to_dof_map(I)[i]]]]]
#q = MyRestrictedFunction(to_be_restricted=r,mesh_in=Omega)
#q = project(q,Q)
plot(q,title="Funzione sul bordo inferiore")
interactive()
print q.vector().str(True)

# ==================================
# prova con class
# =================================
restrictor = RestrictFunction(Omega,MyGamma0())
other_q = restrictor(r,V,)
plot(other_q,title="Prova con classe")
print "q diff", (q.vector()-other_q.vector()).str(True)


##print "Omega v2d ",	dof_map.vertex_to_dof_map(Omega)
##print "Omega d2v ",	dof_map.dof_to_vertex_map(Omega)
#print "bOmega v2d ",	dof_map.vertex_to_dof_map(bOmega)
#print "bOmega d2v ",	dof_map.dof_to_vertex_map(bOmega)
#print "bOmega e ",	bOmega.entity_map(0)
#dof_map  = Q.dofmap()
#print "I v2d ",		dof_map.vertex_to_dof_map(I)
#print "I d2v ",		dof_map.dof_to_vertex_map(I)
##print "I e ",		I.entity_map(0)
#print "bOmega e(I d2v): "
#for i in range(0,len(dof_map.vertex_to_dof_map(I))):
#	print bOmega.entity_map(0)[i]
	



#print dof_map.dofs()
#print "=============="
#print dof_map.str(True)
# print dof_map.__dict__
#plot(bI)
#interactive()
#dof_map = bQ.dofmap()
#print dof_map.vertex_to_dof_map(bI)
#print dof_map.dof_to_vertex_map(bI)
#print emap
#print emap[39]
#print emap(dof_map.vertex_to_dof_map(bI))
#print emap(dof_map.dof_to_vertex_map(bI))
#for x in I.coordinates():
#	dof = dof_map.vertex_to_dof_map(I,x)
#	q.vector()[dof] = r.vector()[dof]
#plot(q,title="Funzione rimappata")
#interactive()
