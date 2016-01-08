from firedrake import *
mesh = Mesh("omega.msh")
V = FunctionSpace(mesh, "CG", 2)
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(1)
bc = DirichletBC(V, 0, 1)
sol = Function(V)
a = dot(grad(u), grad(v))*dx
L = f*v*dx
solve(a == L, sol, bcs=bc)
File("omega.pvd") << sol
