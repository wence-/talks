from firedrake import *
from firedrake.petsc import PETSc

# LES parametrisation breaks with assembly cache on
parameters["assembly_cache"]["enabled"] = False
mesh = Mesh("cylinder.msh")

V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)

u = TrialFunction(V)
p = TrialFunction(Q)
v = TestFunction(V)
q = TestFunction(Q)

dt = 0.0005
t = 0
T = 8

nu = Constant(2e-4)

Um = 1.5
H = 0.41

# Constant parabolic inflow
# Max inflow is (1.5, 0)
inflow_expr = Expression(("4*Um*x[1]*(H-x[1])/(H*H)", "0"),
                         t=t, Um=Um, H=H)

# No slip on side walls and cylinder
noslip = DirichletBC(V, (0, 0), 1)

# Inflow on left
inflow = DirichletBC(V, inflow_expr, 2)

# zero pressure outflow
outflow = DirichletBC(Q, 0, 3)

bcu = [inflow, noslip]
bcp = [outflow]

u0 = Function(V)
u1 = Function(V)

# Initial guess is just the inflow bc
u0.interpolate(inflow_expr)
for bc in bcu:
    bc.apply(u0)

p1 = Function(Q)

f = Constant((0, 0))

k = Constant(dt)

use_les = False
if use_les:
    # Doesn't work?
    S = 0.5*(grad(u0) + grad(u0).T)
    dim = len(u0)
    second_inv = 0
    for i in range(dim):
        for j in range(dim):
            second_inv += 2*S[i, j]**2

    second_inv = sqrt(second_inv)
    width = CellVolume(mesh)**(1.0/dim)

    smagorinsky_coefficient = Constant(0.1)

    L_eddy = q*second_inv*(smagorinsky_coefficient * width)**2*dx
    a_eddy = p*q*dx
    eddy_viscosity = Function(Q)
    A_eddy = assemble(a_eddy)
    b_eddy = assemble(L_eddy)
    eddy_solver = LinearSolver(A_eddy, solver_parameters={'ksp_type': 'cg'},
                               options_prefix="eddy_")
    viscosity = nu + eddy_viscosity
else:
    viscosity = nu


# Tentative velocity step
F1 = (1/k)*inner(u - u0, v)*dx + inner(grad(u0)*u0, v)*dx + \
     viscosity*inner(grad(u), grad(v))*dx - inner(f, v)*dx

a1 = lhs(F1)
L1 = rhs(F1)

# Pressure update
a2 = inner(grad(p), grad(q))*dx
L2 = -(1/k)*div(u1)*q*dx

# Velocity update

a3 = inner(u, v)*dx
L3 = inner(u1, v)*dx - k*inner(grad(p1), v)*dx

A1 = assemble(a1, bcs=bcu)
A2 = assemble(a2, bcs=bcp)
A3 = assemble(a3, bcs=bcu)

tentative = LinearSolver(A1, solver_parameters={'ksp_type': 'gmres',
                                                'ksp_monitor': False},
                         options_prefix="tentative_")

pressure = LinearSolver(A2, solver_parameters={'ksp_type': 'cg',
                                               'pc_type': 'hypre',
                                               'ksp_monitor': False},
                        options_prefix="pressure_")

velocity = LinearSolver(A3, solver_parameters={'ksp_type': 'cg',
                                               'pc_type': 'bjacobi',
                                               'sub_pc_type': 'ilu',
                                               'ksp_monitor': False},
                        options_prefix="velocity_")


b1 = assemble(L1)
b2 = assemble(L2)
b3 = assemble(L3)

output = True

if output:
    # pout = File("/data/lmitche1/models/ns-cylinder/pressure.pvd")
    vout = File("/data/lmitche1/models/ns-cylinder/velocity.pvd")
    vortout = File("/data/lmitche1/models/ns-cylinder/vorticity.pvd")
    Vout = VectorFunctionSpace(mesh, 'CG', 1)

    u = TrialFunction(Vout)
    v = TestFunction(Vout)

    proj = assemble(inner(u, v)*dx)

    projector = LinearSolver(proj, solver_parameters={'ksp_type': 'cg',
                                                      'pc_type': 'bjacobi',
                                                      'sub_pc_type': 'ilu'},
                             options_prefix="project_")

    projrhs = inner(u1, v)*dx

    prhs = assemble(projrhs)

    uout = Function(Vout)

    Vvort = FunctionSpace(mesh, 'CG', 1)

    u = TrialFunction(Vvort)
    v = TestFunction(Vvort)

    vort = Function(Vvort)

    vortp = assemble(u*v*dx)

    vortproj = LinearSolver(vortp, solver_parameters={'ksp_type': 'cg',
                                                      'pc_type': 'bjacobi',
                                                      'sub_pc_type': 'ilu'},
                            options_prefix="vorticity_")


    vortrhs = rot(u1)*v*dx

    vrhs = assemble(vortrhs)

step = 0

import time
start = time.time()

if use_les:
    les_stage = PETSc.Log.Stage("LES")
tentative_stage = PETSc.Log.Stage("tentative")
pressure_stage = PETSc.Log.Stage("pressure")
velocity_stage = PETSc.Log.Stage("velocity")
output_stage = PETSc.Log.Stage("output")
while t <= T:

    if use_les:
        with les_stage:
            b_eddy = assemble(L_eddy, tensor=b_eddy)
            eddy_solver.solve(eddy_viscosity, b_eddy)
    with tentative_stage:
        b1 = assemble(L1, tensor=b1)
        tentative.solve(u1, b1)

    with pressure_stage:
        b2 = assemble(L2, tensor=b2)
        pressure.solve(p1, b2)

    with velocity_stage:
        b3 = assemble(L3, tensor=b3)
        velocity.solve(u1, b3)
        u0.assign(u1)

    t += dt
    if output and (step % 100 == 0):
        # prhs = assemble(projrhs, tensor=prhs)

        # projector.solve(uout, prhs)

        with output_stage:
            vrhs = assemble(vortrhs, tensor=vrhs)

            vortproj.solve(vort, vrhs)
        
            vortout << vort
        # vout << uout
        # pout << p1
    if step % 10 == 0:
        if op2.MPI.comm.rank == 0:
            import sys
            sys.stdout.flush()
            current = time.time()
            print t, 'estimated time to solution %.2g hours' % ((T/t - 1) * (current - start) / (3600.0))
    step += 1

