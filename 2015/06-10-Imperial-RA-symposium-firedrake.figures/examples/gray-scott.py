# Gray-Scott coupled reaction diffusion
# Exploring the parameter space in space (ala Xmorphia)

from firedrake import *
import numpy as np
parameters["matnest"] = False
parameters["assembly_cache"]["enabled"] = False
parameters["coffee"]["O2"] = False
parameters["pyop2_options"]["lazy_evaluation"] = False
# Define the discretisation settings
n = 250
dt = Constant(0.3) # Timestep
d = 2.5

mesh = RectangleMesh(n, n, d, d, quadrilateral=True)

V = FunctionSpace(mesh, "CG", 1)

Z = MixedFunctionSpace([V, V])

u, v = TrialFunctions(Z)
z = Function(Z, name = "z")
z_old = Function(Z, name = "z_old")
u_old, v_old = split(z_old)

q, p = TestFunctions(Z)

z_vtu = File("results/z.pvd")


# Define the model constants
eps1 = Constant(2*10**-5) # Diffusion coefficient for variable u
eps2 = Constant(10**-5) # Diffusion coefficient for variable v

k = Function(V)
F = Function(V)

F.interpolate(Expression("0.08*(x[1]/d) + 0.006", d=d))
k.interpolate(Expression("0.04*(x[0]/d) + 0.03", d=d))
# k = Constant(0.053) # Rate constant for the second reaction
# F = Constant(0.022) # Feed rate


# Define the weak formulation of the Gray-Scott equations
a = inner(u, q)*dx + inner(v, p)*dx

L1 = inner(u_old, q)*dx + dt*(inner(F*(1 - u_old), q) - eps1 * inner(grad(u_old), grad(q)) - inner(u_old*v_old*v_old, q))*dx

L2 = inner(v_old, p)*dx + dt*(inner(u_old*v_old*v_old, p) - eps2 * inner(grad(v_old), grad(p)) - inner((k + F) * v_old, p))*dx

L = L1 + L2


solver = LinearSolver(assemble(a), solver_parameters={'ksp_type': 'cg',
                                                      'pc_type': 'icc',
                                                      'sub_pc_type': 'icc'})


u, v = z_old.split()

u.dat.data[:] = 0.5 * (1.0 +  0.01 * np.random.uniform(0, 1, u.dof_dset.size))
v.dat.data[:] = 0.25 * (1.0 + 0.01 * np.random.uniform(0, 1, v.dof_dset.size))

b = assemble(L)

# Run the timeloop
u, v = z.split()
import time
start = time.time()
nstep = 200000
for step in range(nstep):
    last = time.time()
    b = assemble(L, tensor=b)
    solver.solve(z, b)
    if step % 100 == 0:
        print 'estimated tts %.3g hours (step took %.2g seconds)' % (((float(nstep)/(step+1) - 1)*(time.time() - start)/3600.0), (time.time() - last))
        z_vtu << u
    z_old.assign(z)
    step += 1
z_vtu << u
