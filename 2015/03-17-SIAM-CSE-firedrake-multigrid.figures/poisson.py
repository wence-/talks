from firedrake import *

mesh = UnitSquareMesh(10, 10)

nlevel = 6

mh = MeshHierarchy(mesh, nlevel)

V = FunctionSpaceHierarchy(mh, 'CG', 1)

u = FunctionHierarchy(V)

exact_expr = Expression("pow(x[0]*x[0] + x[1]*x[1], 1.0/3) * sin(2.0/3.0 * (atan2(x[1], x[0]) < 0 ? atan2(x[1], x[0]) + 2*pi : atan2(x[1], x[0])))")

v = TestFunction(V[-1])

u_ = u[-1]

F = dot(grad(u_), grad(v))*dx

bcs = DirichletBC(V[-1], exact_expr, (1, 2, 3, 4, 7))

problem = NonlinearVariationalProblem(F, u_, bcs=bcs)

solver = NLVSHierarchy(problem, solver_parameters={'snes_type': 'ksponly',
                                                   'ksp_type': 'cg',
                                                   'pc_type': 'hypre',
                                                   'ksp_monitor': True,
                                                   'ksp_max_it': 3,
                                                   'mg_levels_ksp_type': 'chebyshev',
                                                   'mg_levels_pc_type': 'jacobi',
                                                   'mg_levels_ksp_max_it': 3,
                                                   'ksp_convergence_test': 'skip',})

import time

start = time.time()
solver.solve()
print time.time() - start

# u[-1].assign(0)
# start = time.time()
# solver.solve()
# print time.time() - start

# exact = Function(V[-1])

# exact.interpolate(exact_expr)

# # File("u.pvd") << u_
# # File("exact.pvd") << exact
# # File("diff.pvd") << assemble(exact - u_)
# print norm(assemble(exact - u_))

print V[-1].dof_dset.size
