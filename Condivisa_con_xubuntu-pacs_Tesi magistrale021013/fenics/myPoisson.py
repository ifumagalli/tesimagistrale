from dolfin import *
from numpy import *
import sys

### domain and mesh
n = 5
mesh = UnitSquare(n,n)
deg = 1
# or from command line
# IN WINDOWS NON VA: PROVARE CON MACCHINA VIRTUALE
# #print dir(sys)
# info(sys)
# if len(sys.argv)>=2:
	# domain_type = [UnitInterval,UnitSquare,UnitCube]
	# divisions = [int(arg) for arg in sys.argv[:len(argv)-2]]
	# mesh = domain_type[len(divisions)-1](*divisions)
	# deg = int(sys.argv[len(argv)-1])

### data
f = Expression("-6.0")
u0 = Expression("1+x[0]*x[0]+2*x[1]*x[1]")
uex = u0

### weak form
V = FunctionSpace(mesh,"CG",deg)
u = TrialFunction(V)
v = TestFunction(V)
a = inner(nabla_grad(u),nabla_grad(v))*dx
	# NB grad and nabla_grad are the same for scalars,
	#    but for vectors nabla_grad performs in a more intuitive way
L = f*v*dx

### boundary conditions
def u0_boundary(x,on_boundary):
	return on_boundary
bc = DirichletBC(V,u0,u0_boundary)

### settings for solution
parameters["linear_algebra_backend"] = "uBLAS"
prm = parameters["krylov_solver"] # short form
prm["absolute_tolerance"] = 1E-10
prm["relative_tolerance"] = 1E-6
prm["maximum_iterations"] = 1000
set_log_level(PROGRESS)
	# DEBUG for more info, comment line to set the lowest verbosity

### solution
u = Function(V)
solve(a==L,u,bc,
	solver_parameters={ "linear_solver":"cg",
						"preconditioner":"ilu"} )
"""Instead of the above line, the following statements can be used
	(NB: instead of using the variable 'parameters', 'solver.parameters' has to be used):
	problem = LinearVariationalProblem(a, L, u, bc)
	solver = LinearVariationalSolver(problem)
	solver.solve()
	"""

### post-processing
coor = mesh.coordinates()
u_nodal_values = u.vector()
u_array = u_nodal_values.array()	# to have numpy Python array
if mesh.num_vertices() == len(u_array):
	# if elements are Lagrangian with degree >1, if-test fails (nodes != vertices)
	for i in range(mesh.num_vertices()):
		print 'u(%8g,%8g) = %g' % (coor[i][0], coor[i][1], u_array[i])

print "------- Errors -------"
uex_array = interpolate(uex,V).vector().array()
h = (mesh.hmax()+mesh.hmin())/2.0
cellArea = h*h/2.0
errL1 = sum(abs(u_array-uex_array)*cellArea)
errL2 = vdot(abs(u_array-uex_array)*cellArea,abs(u_array-uex_array))
errLinf = abs(u_array-uex_array).max()
uexF = project(uex,V)
errH1_form = inner(nabla_grad(u-uexF),nabla_grad(u-uexF))*dx
errH1 = assemble(errH1_form)
print "   Error in L1 = %g" % errL1
print "   Error in L2 = %4g" % errL2
print "   Error in Linf = %4g" % errLinf
print "   Error in H1 = %4g" % errH1
errL1_form = abs(u-uexF)*dx
errL1 = assemble(errL1_form) # computes integral in errL1_form
print "   Control error in L1 = %g" % errL1
print "---- Errors - end ----"

u_grad = project(nabla_grad(u),VectorFunctionSpace(mesh,"CG",1))
u_dx, u_dy = u_grad.split(deepcopy=True) #extract components as Functions

### plots
plot(mesh,title="FE mesh")
file = File("myPoisson.pvd")
file << u
viz_u = plot(u,
			wireframe=False,
			title="Scaled solution",
			rescale=False,
			axes=True, # include axes
			basename="myPoisson", # default plotfile name
			)
viz_u.elevate(-65) # tilt camera -65 degrees (latitude dir)
viz_u.set_min_max(0, 0.5*u_array.max()) # color scale
viz_u.update(u) # bring settings above into action
viz_u.write_png("myPoisson.png")
viz_u.write_ps("myPoisson", format="eps")
#plot(u,wireframe=True,title="solution (wireframe=True)")

viz_u_grad = plot(u_grad,basename="myPoisson_grad")
viz_u_grad.write_png("myPoisson_grad.png")

#interactive() # to hold plots: to continue, you need to press 'q' in each plot window
raw_input("Press ENTER to continue")

### altre cose interessanti
#info(parameters,True)

print ""
print "Bien: On a fini Normalement! :)"