from firedrake import *

mesh = Mesh("doughnut.msh")

nlevel = 5
mh = MeshHierarchy(mesh, nlevel)

V = FunctionSpaceHierarchy(mh, "CG", 1)

u = FunctionHierarchy(V)

u = u[-1]
V = V[-1]

w = TestFunction(V)
x = SpatialCoordinate(V.mesh())

r = sqrt(x[0]**2 + x[1]**2)
rho = 1.0/r**3

F = (8*inner(grad(u), grad(w)) + rho * inner(u**5, w) + \
     Constant(-1.0/10.0)*inner(u, w))*dx(degree=5)

bcs = DirichletBC(V, Constant(1.0), (21, 22))

u.assign(1)

problem = NonlinearVariationalProblem(F, u, bcs=bcs)
fas_parameters = {
                  "snes_type": "fas",
                  "snes_monitor": None,
                  "snes_ngmres_monitor": None,
                  "snes_max_it": 100,
                  "snes_stol": 0.0,
                  "snes_rtol": 0.0,
                  "snes_atol": 1.0e-10,
                  "snes_linesearch_type": "basic",
                  "snes_converged_reason": None,
                  "fas_levels_snes_monitor": None,
                  "fas_coarse_snes_monitor": None,
                  "snes_fas_type": "full",
                  "fas_coarse_pc_type": "lu",
                  "fas_coarse_snes_max_it": 200,
                  "fas_coarse_snes_atol": 1.0e-11,
                  "fas_coarse_snes_rtol": 0.0,
                  "fas_coarse_snes_type": "newtonls",
                  "fas_coarse_snes_linesearch_type": "l2",
                  "fas_coarse_ksp_type": "gmres",
                  "fas_levels_snes_type": "ngs",
                  "fas_levels_snes_max_it": 5,
                  "fas_levels_ksp_type": "preonly",
                  "fas_levels_ksp_max_it": 1,
                  "fas_levels_pc_type": "lu",
                  "fas_levels_snes_linesearch_type": "basic",
                  "fas_levels_snes_atol": 1.0e-11,
                  "fas_levels_snes_rtol": 0.0,
                  "fas_levels_snes_norm_schedule": "always", # PETSc bug
                  "fas_levels_ksp_convergence_test": "skip", # PETSc bug
                  "fas_levels_snes_convergence_test": "skip",
                  "snes_view": None,
                  "mat_mumps_icntl_14": 100, # PETSc bug
                 }

newton_parameters = {"snes_type": "newtonls",
                     "snes_linesearch_type": "basic",
                     "pc_type": "mg",
                     "ksp_max_it": 10,
                     "ksp_type": "gmres",
                     "ksp_convergence_test": "skip",
                     "ksp_monitor": True,
                     "mg_coarse_pc_type": "redundant",
                     "mg_coarse_redundant_pc_type": "lu",
                     "mg_coarse_ksp_type": "preonly",
                     "mg_coarse_ksp_max_it": 30,
                     "mg_levels_ksp_type": "chebyshev",
                     "mg_levels_ksp_max_it": 3,
                     "mg_levels_pc_type": "bjacobi",
                     "mg_levels_sub_pc_type": "ilu",
                     "snes_monitor": True}
solver = NLVSHierarchy(problem,
                       solver_parameters=newton_parameters,
                       options_prefix="")

solver.solve()

File("solution.pvd") << u
                       
print V.dof_dset.size
