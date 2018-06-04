from firedrake import *
from firedrake.petsc import *

N = 4 # base mesh size
L = 5  # number of levels employed in multigrid
p = 2  # element degree used (velocity degree for Taylor-Hood)


distribution_parameters={"partition": True, "overlap_type": (DistributedMeshOverlapType.VERTEX, 1)}

mesh = UnitSquareMesh(N, N, distribution_parameters=distribution_parameters, quadrilateral=False)
# Construct mesh hierarchy
mh = MeshHierarchy(mesh, L, distribution_parameters=distribution_parameters)
mesh = mh[-1] # mesh to solve on

# Construct function space
Ve = VectorElement("P", mesh.ufl_cell(), p)
Qe = FiniteElement("P", mesh.ufl_cell(), p-1) # high-order Taylor-Hood
Ze = MixedElement([Ve, Qe])
Q  = FunctionSpace(mesh, Qe)
Z  = FunctionSpace(mesh, Ze)

x, y = SpatialCoordinate(mesh)
x = x
y = y
Re = (1 + 1e2 * exp(-((x - 0.5)**2 + (y - 0.5)**2)*1000) +
      1e7*exp(-((x - 0.75)**2 + (y - 0.25)**2)*400) +
      1e1*exp(-((x - 0.25)**2 + (y - 0.75)**2)*250) +
      1e4*exp(-((x - 0.8)**2 + (y - 0.9)**2)*400) +
      1e5*exp(-((x - 0.1)**2 + (y - 0.8)**2)*350))
z = Function(Z)
w = TestFunction(Z)
(u, p) = split(z)
(v, q) = split(w)

# Define PDE
F = (
      1.0/Re * inner(grad(u), grad(v)) * dx
    - p * div(v) * dx
    - div(u) * q * dx + dot(Constant((3, 3)), v)*dx
    )

solution = as_vector([x**2 + y**2, 2*x**2 - 2*x*y])
# Dirichlet BCs for driven cavity, and nullspace to eliminate pressure scaling
nsp = MixedVectorSpaceBasis(Z, [Z.sub(0), VectorSpaceBasis(constant=True)])
bcs = [DirichletBC(Z.sub(0), solution, [1, 2, 3, 4])]

solver_parameters = {
    "mat_type": "matfree",
    "snes_type": "ksponly",
    "ksp_rtol": 1e-8,
    "ksp_type": "fgmres",
    "pc_type": "mg",
    "mg_levels": {
        "ksp_type": "gmres",
        "ksp_max_it": 5,
        "pc_type": "python",
        "pc_python_type": "ssc.PatchPC",
        "patch_pc_patch_save_operators": True,
        "patch_pc_patch_construction_type": "vanka",
        "patch_pc_patch_partition_of_unity": True,
        "patch_pc_patch_vanka_dim": 0,
        "patch_pc_patch_construction_dim": 0,
        "patch_pc_patch_exclude_subspace": 1,
        "patch_pc_patch_sub_mat_type": "seqaij",
        "patch_sub_ksp_type": "preonly",
        "patch_sub_pc_type": "lu",
        "patch_sub_pc_factor_shift_type": "nonzero"
    },
    "mg_coarse": {
        "ksp_type": "preonly",
        "pc_type": "python",
        "pc_python_type": "firedrake.AssembledPC",
        "assembled_pc_type": "svd",
    }
}

problem = NonlinearVariationalProblem(F, z, bcs)

solver = NonlinearVariationalSolver(problem, nullspace=nsp,
                                    solver_parameters=solver_parameters,
                                    options_prefix="", appctx=appctx)

pvd = File("output/output.pvd")

# Actually solve the problem
solver.solve()

u, p = z.split()
u.rename("Velocity")
p.rename("Pressure")

P = FunctionSpace(mesh, "DP", 0)
nu = Function(P, name="Viscosity")
nu.interpolate(1/Re)

pvd.write(u, p, nu)
