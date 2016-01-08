from firedrake import *

m = UnitSquareMesh(10, 10)

nlevel = 4

mh = MeshHierarchy(m, nlevel)

V = FunctionSpaceHierarchy(mh, "CG", 1)

u = FunctionHierarchy(V)

V = V[-1]
u = u[-1]

v = TestFunction(V)

delta = Constant(1)

F = delta*inner(grad(v), grad(u))*dx + 1.0/delta * inner(v, u**3 - u)*dx

bcs = [DirichletBC(V, Constant(1), (1, 2)),
       DirichletBC(V, Constant(-1), (3, 4))]

problem = NonlinearVariationalProblem(F, u, bcs=bcs)

fas_parameters = {
    "snes_type": "newtonls",
    "pc_type": "mg",
    "pc_type": "lu",
    "ksp_type": "gmres",
    "ksp_max_it": 7,
    "mg_levels_ksp_max_it": 4,
    "mg_levels_ksp_type": "richardson",
    "mg_levels_ksp_richardson_self_scale": True,
    "mg_levels_pc_type": "sor",
    "mg_levels_ksp_chebyshev_estimate_eigenvalues": True,
    "mg_levels_ksp_chebyshev_estimate_eigenvalues_random": True,
    "snes_monitor": True,
    "ksp_convergence_test": "skip",
    "snes_max_it": 100,
    "snes_linesearch_type": "l2",
    "snes_converged_reason": True,
    "snes_fas_type": "full",
                  #"snes_linesearch_monitor": None,
                  #"fas_levels_snes_linesearch_monitor": None,
                  #"fas_coarse_snes_linesearch_monitor": None,
                 }

solver = NLVSHierarchy(problem, solver_parameters=fas_parameters)

solver.solve()

File("u.vtu") << u
