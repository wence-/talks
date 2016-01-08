from firedrake import *

mesh = UnitSquareMesh(10, 10)

nlevel = 4

mh = MeshHierarchy(mesh, nlevel)

V = FunctionSpaceHierarchy(mh, 'CG', 1)

u = FunctionHierarchy(V)

coeff_expr = Expression("""(((int)(floor(x[0]*H)) % 2 == 0 &&
                             (int)(floor(x[1]*H)) % 2 == 0) ? 20 :
                            ((int)(floor(x[0]*H)) % 2 == 1 &&
                             (int)(floor(x[1]*H)) % 2 == 0) ? 0.002 :
                            ((int)(floor(x[0]*H)) % 2 == 0 &&
                             (int)(floor(x[1]*H)) % 2 == 1) ? 0.2 :
                            ((int)(floor(x[0]*H)) % 2 == 1 &&
                             (int)(floor(x[1]*H)) % 2 == 1) ? 2000 : 10000)""",
                        H=10)


DG = FunctionSpaceHierarchy(mh, "DG", 0)

C = FunctionHierarchy(DG)
for c in C:
    c.interpolate(coeff_expr)
c = C[-1]


V = V[-1]
v = TestFunction(V)

u_ = u[-1]

F = c*dot(grad(u_), grad(v))*dx - Constant(1)*v*dx

bcs = DirichletBC(V, Constant(0), (1, 2, 3, 4))

problem = NonlinearVariationalProblem(F, u_, bcs=bcs)

fas_parms = {"snes_type": "fas",
             "snes_fas_type": "full",
             "fas_coarse_snes_type": "newtonls",
             "fas_coarse_ksp_type": "preonly",
             "fas_coarse_pc_type": "redundant",
             "fas_coarse_redundant_pc_type": "lu",
             "fas_coarse_snes_linesearch_type": "basic",
             "fas_levels_snes_type": "newtonls",
             "fas_levels_snes_linesearch_type": "basic",
             "fas_levels_snes_max_it": 1,
             "fas_levels_ksp_type": "chebyshev",
             "fas_levels_ksp_max_it": 2,
             "fas_levels_pc_type": "jacobi",
             "fas_levels_ksp_convergence_test": "skip",
             "snes_max_it": 1,
             "snes_monitor": True,
             "snes_convergence_test": "skip"}
mg_parms = {'snes_type': 'ksponly',
            'ksp_type': 'cg',
            'ksp_monitor': True,
            'mg_levels_ksp_type': 'richardson',
            'mg_levels_ksp_chebyshev_estimate_eigenvalues': True,
            'mg_levels_ksp_chebyshev_estimate_eigenvalues_random': True,
            'mg_levels_pc_type': 'ilu',
            'mg_levels_ksp_max_it': 2}

solver = NLVSHierarchy(problem, solver_parameters=mg_parms,
                       options_prefix="")

import time

start = time.time()
solver.solve()
print time.time() - start

File("coeff.vtu") << c
File("u.vtu") << u_


File("coeff0.vtu") << C[0]
