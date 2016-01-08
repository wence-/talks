from firedrake import *

mesh = Mesh("wave.msh")
V = FunctionSpace(mesh, 'CG', 1)

mass_lump = False
output = True

t = 0
dt = 0.00075
T = 8
omega = 4*2*pi
forcing = Constant(omega*t)

phi = Function(V)
p = Function(V)

bc = DirichletBC(V, forcing, 1)

v = TestFunction(V)

update = dt*dot(grad(v), grad(phi))*dx

if mass_lump:
    lumped = assemble(v*dx)
else:
    u = TrialFunction(V)
    A = assemble(u*v*dx, bcs=bc)
    solver = LinearSolver(A, solver_parameters={'ksp_type': 'cg',
                                                'pc_type': 'bjacobi',
                                                'sub_pc_type': 'ilu',
                                                'ksp_monitor': False,})
    update += v*p*dx

b = assemble(update)

if output:
    pout = File("/data/lmitche1/models/wave-2d/p.pvd")
    phiout = File("/data/lmitche1/models/wave-2d/phi.pvd")

step = 0
while t <= T:
    forcing.assign(cos(omega*t))
    phi -= (dt/2) * p

    if mass_lump:
        p += assemble(update, tensor=b)/lumped
        bc.apply(p)
    else:
        b = assemble(update, tensor=b)
        solver.solve(p, b)

    phi -= (dt/2) * p
    
    t += dt

    if output and (step % 10 == 0):
        pout << p
        phiout << phi
        if op2.MPI.comm.rank == 0:
            print t

    step += 1
