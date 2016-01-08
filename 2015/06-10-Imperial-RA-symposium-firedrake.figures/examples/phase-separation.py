import random
from firedrake import *
parameters["matnest"] = False
parameters["coffee"]["O2"] = True
parameters["assembly_cache"]["enabled"] = False

# Class representing the intial conditions
class InitialConditions(Expression):
    def __init__(self):
        # random.seed(2 + op2.MPI.comm.rank)
        super(InitialConditions, self).__init__()

    def eval(self, values, x):
        values[0] = 1 if random.random() < 0.5 else -1
    def value_shape(self):
        return ()

# Model parameters
dt = 5.0e-06  # time step

twoD = True
output = True
# Create mesh and define function spaces
if twoD:
    mesh = UnitSquareMesh(200, 200)
else:
    mesh = UnitCubeMesh(20, 20, 20)


V = FunctionSpace(mesh, "Lagrange", 1)
W = V*V


# homogeneous free energy:
# f = alpha*(c**2 - 1)**2
# dc/dt = div D(c) grad (mu)
# mu = df/dc - gamma div grad c
gamma = Constant(0.005)
D = Constant(10)

q, v  = TestFunctions(W)

u   = Function(W, name="u_(n+1)") # current solution
u0  = Function(W, name="u_n")     # solution from previous converged step

# Split mixed functions
c,  mu  = split(u)
c0, mu0 = split(u0)


c = variable(c)
alpha = Constant(10)
f = alpha*(c**2 - 1)**2
fprime = diff(f, c)

# mu_(n+theta)
# Crank Nicholson timestepping scheme
theta  = 0.5
mu_theta = (1.0-theta)*mu0 + theta*mu

# Weak statement of the equations
L0 = (c - c0)*q*dx + dt*D*dot(grad(mu_theta), grad(q))*dx
L1 = (mu - fprime)*v*dx - gamma*dot(grad(c), grad(v))*dx
F = L0 + L1

# Create nonlinear problem and Newton solver
problem = NonlinearVariationalProblem(F, u)
solver = NonlinearVariationalSolver(problem, options_prefix="ch_")

c, mu = u.split()
c.rename("order parameter")
c.interpolate(InitialConditions())
# Output file
t = 0.0

# This is conserved, since dc/dt == div(j(x)); j(x) == D grad(mu)
order_parameter = c*dx
# This decays to zero (possibly walking over humps)
free_energy = (f + gamma/2 * dot(grad(c), grad(c)))*dx

if output:
    if twoD:
        cout = File("ch_output/order-parameter-2D.pvd")
    else:
        cout = File("ch_output/order-parameter-3D.pvd")
    cout << (c, t)

# Step in time
nstep = 400
if op2.MPI.comm.rank == 0:
    print assemble(order_parameter), assemble(free_energy)
else:
    assemble(order_parameter), assemble(free_energy)
import time, sys
start = time.time()
for _ in range(nstep):
    t += dt
    u0.assign(u)
    try:
        solver.solve()
    except RuntimeError:
        pass
    if op2.MPI.comm.rank == 0:
        sys.stdout.flush()
        current = time.time()
        print t, 'estimated time to solution %.2g hours' % ((float(nstep)/(_+1) - 1) * (current - start) / (3600.0))
        print assemble(order_parameter), assemble(free_energy)
    else:
        assemble(order_parameter), assemble(free_energy)
    if output:
        cout << (c, t)
