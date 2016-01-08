from firedrake import *

m = RectangleMesh(20, 20, 1, 1)


mh = MeshHierarchy(m, 4)


for m_ in mh:
    m_.coordinates *= 2
    m_.coordinates -= 1
Vh = FunctionSpaceHierarchy(mh, 'CG', 1)


V = Vh[-1]

v = TestFunction(V)

u = FunctionHierarchy(Vh)

u_ = u[-1]

parameters["coffee"]["O2"] = False
p = Constant(4)
lmbda = Constant(0)
epsilon = Constant(1e-2)
source = Constant(1)
nu = (epsilon**2 + 0.5 * inner(grad(u_), grad(u_)))**((p - 2.0)/2.0)


initial_expr = Expression("(1 - x[0]*x[0]) * (1 - x[1]*x[1]) * x[0] * x[1]")

for _u in u:
    _u.interpolate(initial_expr)
u_.interpolate(initial_expr)
# u_.assign(0)
# p-Bratu
# F = inner(nu*grad(u_), grad(v))*dx(degree=4) - lmbda*exp(u_)*v*dx(degree=4) - inner(source, v)*dx

# p-Laplace
F = inner(nu*grad(u_), grad(v))*dx(degree=4) - v*dx
bcs = DirichletBC(V, Constant(0), (1, 2, 3, 4))
problem = NonlinearVariationalProblem(F, u_, bcs=bcs)
parameters = {
    "snes_type": "newtonls",
    "snes_monitor": True,
    # "snes_rtol": 1e-15,
    # "pc_type": "mg",
    # "ksp_max_it": 10,
    # "ksp_convergence_test": "skip",

    "npc_snes_type": "fas",
    "npc_snes_fas_type": "full",
    "npc_snes_monitor": True,
    "npc_fas_coarse_snes_type": "newtonls",
    "npc_fas_coarse_snes_max_it": 20,
    "npc_fas_coarse_ksp_type": "preonly",
    "npc_fas_coarse_pc_type": "lu",
    "npc_fas_coarse_snes_linesearch_type": "bt",
    "npc_fas_levels_snes_type": "ncg",
    "npc_fas_levels_snes_max_it": 2,
    "npc_fas_levels_snes_rtol": 1e-4,
    "npc_fas_levels_snes_monitor": True,
}
solver = NLVSHierarchy(problem, solver_parameters=parameters,
                       options_prefix="")

solver.solve()

File("u.vtu") << u_

print op2.MPI.comm.allreduce(V.dof_dset.size)
