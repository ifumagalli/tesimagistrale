""" In this file we collect some hand-made utilities for the program
"""

from dolfin import *
from copy import deepcopy, copy
import numpy
import geometry


# ===== "macro" for debugging =====
MYDEBUG = False


# ===== Functor to restrict a Function =====
class RestrictFunction(object):
	""" Functor which returns a Function restricted on a part of the boundary
	"""
	def __init__(self,big_mesh,restr_line = None):#geometry.moving_boundary):
		self._big_mesh = big_mesh
		self._big_boundary = BoundaryMesh(big_mesh,"exterior")
		mfunc = MeshFunction("uint",self._big_boundary,1)
		# ffunc = FacetFunction("uint",bOmega)
			# con questa non va
		if restr_line == None:
			restr_line = geometry.MovingBoundary()
		mfunc.set_all(1)
		restr_line.mark(mfunc,0)
			# ??? non si modifica restriction_line, giusto?
#		print "====A"
		self._restr_mesh = SubMesh(self._big_boundary,mfunc,0)
#		print self._restr_mesh.geometry().str(True)
#		print self._restr_mesh.topology().str(True)
#		print "====B"
		self._restr_space = FunctionSpace(self._restr_mesh,'CG',1)
			# !!! per ora usiamo Lagrange P1: si puo' generalizzare?
	#d#	self._restr_map = self._big_boundary.entity_map(0)[\
	#d#			self._restr_mesh.data().mesh_function("parent_vertex_indices")[\
	#d#			self._restr_space.dofmap().vertex_to_dof_map(self._restr_mesh)]]

	def __call__(self,fun):
		# !!! mi sa che e' necessario che fun_space sia Lagrange P1, per ora...
		restr = Function(self._restr_space)
		for i in range(0,self._restr_mesh.size(0)):
			# ??? c'e' un modo piu' veloce del for per fare cio' che segue?
			restr.vector()[i] = fun.vector()[\
					fun.function_space().dofmap().dof_to_vertex_map(self._big_mesh)[\
				#d#	self._restr_map[\
						# ??? vorrei poter usare una sintassi del tipo #d# riportata nel costruttore
						# per creare una mappa al momento della costruzione, cosi' da accorciare i tempi di valutazione
					self._big_boundary.entity_map(0)[\
					self._restr_mesh.data().mesh_function("parent_vertex_indices")[\
					self._restr_space.dofmap().vertex_to_dof_map(self._restr_mesh)[\
					i]]]]]
		return restr

	def get_restr_mesh(self):
		return self._restr_mesh

	def get_restr_space(self):
		return self._restr_space


# ===== Class to restrict a Function =====
class MyRestrictedFunction(Expression):
	""" Class to allow restriction of a Function to the bottom boundary
	"""
	# DEVO USARE Expression PERCHE' CON Function NON POSSO SOVRASCRIVERE eval
	# TODO cosi' com'e' mi sa di espediente sporco
	def __init__(self,*args,**kwargs):
		self._to_be_restricted = kwargs.get('to_be_restricted',None)
		self._mesh_in = kwargs.get('mesh_in',None)
		if self._to_be_restricted == None or self._mesh_in == None:
			exit("@@@ Missing argument(s)")
		# self._V_in = kwargs.get('V_in')
		# self._to_be_restricted = project(self._to_be_restricted,self._V_in)
		
	def eval(self,value,x):
		#print self._to_be_restricted(0.5,0.0)
		point = Point(x[0],0.0)
		point = self._mesh_in.closest_point(point)
		# print self._to_be_restricted(0.5,0.0), self._to_be_restricted(point[0],point[1]), point[0], point[1], value
		value = self._to_be_restricted(point[0],point[1])
		# value = []
		# self._to_be_restricted.eval(value,numpy.array([point[0],point[1]],'d'))


# ===== Class to compound with Tq^-1 =====
class TildeFunction(Expression):
	# DEVO USARE Expression PERCHE' CON Function NON POSSO SOVRASCRIVERE eval
	# TODO cosi' com'e' mi sa di espediente sporco
	def __init__(self,*args,**kwargs):
		self._fun = kwargs.get('fun')
		self._map = kwargs.get('map')
		if self._fun == None:
			exit("@@@ Missing the function to remap")
		if self._map == None:
			exit("@@@ Missing the map")
		
	def eval(self,value,x):
		# TODO
		xtilde = self._map.inverse(x)
		value = []
		if xtilde : # se e' dentro il dominio
			value = self._fun(xtilde)
		else:
			value = [0.0 for i in range(1,self._fun.geometric_dimension())]


# ===================================
# programmini di prova
# ===================================
if __name__=="__main__":
	class MyClass(object):
		def __init__(self,other_my_class = None):
			self.id = 0
			if other_my_class != None:
				self.id = other_my_class.id + 1
		# def __copy__(self,memo):
			# cls = self.__class__
			# result = cls.__new__(cls)
			# memo[id(self)] = result
			# for k,v in self.__dict__.items():
				# setattr(result,k,deepcopy(v,memo))
			# return result

	e = MyClass()
	ee = MyClass(e)
	a = [e,ee,MyClass(ee)]
	b = deepcopy(a)
	# tmp = tuple(a)
	c = []
	for item in b:
		c.append(item)
	print "a.id =", [item.id for item in a]
	print "b.id =", [item.id for item in b]
	print "c.id =", [item.id for item in c]
	a[0].id = 6
	b[1].id = 7
	c[2].id = 8
	print "a.id =", [item.id for item in a]
	print "b.id =", [item.id for item in b]
	print "c.id =", [item.id for item in c]
	
	#======================================
	print "\n=== Test for RestrictFunction ==="
	
	parameters["allow_extrapolation"] = True	# useful for evaluating on boundaries (see footnote on p.116 of the manual)

	N = 25
	Omega = UnitSquare(N,N,"crossed")
	V = FunctionSpace(Omega,'CG',1)
	
	class MyGamma0(SubDomain):
		def inside(self,x,on_boundary):
			#return near(x[1],0)
			return near(x[0],1) #or near(x[1],0)

	r = Function(V)
	r_expr = Expression('sin(4*pi*x[0])*cos(6*pi*x[1])')
	#r_expr = Expression('1-x[1]')#+25.0*x[1]')
	r = project(r_expr,V)
	#for i in range(0,len(r.vector())):
	#	r.vector()[i] = i
	pfig = plot(r,title="Funzione di partenza",axes=True)

	restrictor = RestrictFunction(Omega)#,MyGamma0())
	q = restrictor(r)
	
	plot(q,title="Prova con classe")

	Q = FunctionSpace(UnitInterval(N),'CG',1)
	pq = project(q,Q)
	plot(pq,title="q proiettata - NO")
		# !!! proietta in maniera sbagliata
	# pr = project(r,Q)
	# plot(pr,title="r proiettata")
		# !!! le righe precedenti ovviamente non vanno: da' "matrices not aligned" 

	interactive()

#	# ==============
#	# definizione di sottomesh e sottospazi
#	# ==============
#	bOmega = BoundaryMesh(Omega,"exterior")
#	bV = FunctionSpace(bOmega,'CG',1)
#	br = Function(bV)
#	
#	emap = bOmega.entity_map(0)
#	#@v# emap = bOmega.vertex_map()
#	# print emap
#	# print emap[0],emap[7],emap[19],emap[39],emap[41],emap[56]
#	print "=== br.vector()"
#	for i in range(0,bOmega.size(0)):
#		br.vector()[i] = r.vector()[emap[i]]
#		print br.vector()[i]
#	
#	# ==============
#	# vera restrizione
#	# ==============
#	class MyGamma0(SubDomain):
#		def inside(self,x,on_boundary):
#			return near(x[1],0)
#			#return near(x[0],1)# or near(x[1],0)
#	Gamma_moving = MyGamma0()
#	Gamma_moving.mark(mfunc,0)
#	mfunc = MeshFunction("uint",bOmega,1)
#	#ffunc = FacetFunction("uint",bOmega)
#	mfunc.set_all(1)
#
#	I = SubMesh(bOmega,mfunc,0)
#	plot(Omega,title="Omega")
#	plot(bOmega,title="bOmega")
#	plot(I,title="I")
#	interactive()
#	parent_vertex = I.data().mesh_function("parent_vertex_indices")
#
#	Q = FunctionSpace(I,'CG',1)
#	q = Function(Q)
#	for i in range(0,I.size(0)):
#		#q.vector()[i] = br.vector()[map[parent_vertex[i]]]#[I.data().array("parent_vertex_indices")[i]]
#		#q.vector()[i] = r.vector()[emap[parent_vertex[i]]]
#		#q.vector()[i] = r.vector()[emap[parent_vertex[i]]]
#		#q.vector()[i] = br(bOmega.coordinates()[i])
#		q.vector()[i] = r.vector()[V.dofmap().dof_to_vertex_map(Omega)[bOmega.entity_map(0)[parent_vertex[Q.dofmap().vertex_to_dof_map(I)[i]]]]]
#	#q = MyRestrictedFunction(to_be_restricted=r,mesh_in=Omega)
#	#q = project(q,Q)
#	plot(q,title="Funzione sul bordo inferiore")
#	interactive()
#	print q.vector().str(True)


	ans = raw_input("Vuoi continuare con il test vecchio?  ")
	if ans=='n' or ans=='no':
		exit(0)

	#======================================
	print " === Old test for MyRestrictedFunction ==="
	mesh = UnitSquare(10,10)
	V = FunctionSpace(mesh,'CG',2)
	u = TrialFunction(V)
	v = TestFunction(V)
	form = dot(grad(u),grad(v)) * dx
	rhs = Constant('0.0')*v * dx
	bc1 = DirichletBC(V,Constant('0.0'),lambda x, on_boundary : on_boundary and not near(x[1],0.0))
	bc2 = DirichletBC(V,Constant('1.0'),lambda x, on_boundary : on_boundary and near(x[1],0.0))
	bc = [bc1,bc2]
	u = Function(V)
	solve(form==rhs, u, bc)
	pfig = plot(u,title="Solution")
	pfig.write_png("prova.png")
	interactive()
	
	point = Point(0.5,0.5)
	print mesh.closest_point(point)[0], mesh.closest_point(point)[1]
	
	u_restr = MyRestrictedFunction(to_be_restricted=u,mesh_in=mesh)
	print "uscita ", u(0.5,0.5),u_restr(0.5)
	u_restr = Constant('1.0')
	plot(u_restr,mesh=UnitInterval(10))
	interactive()
