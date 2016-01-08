# FIXME not working right, I think.
from firedrake import *
parameters["matnest"] = False
m = UnitSquareMesh(50, 50)

V = VectorFunctionSpace(m, 'CG', 2)
Q = FunctionSpace(m, 'CG', 1)
T = FunctionSpace(m, 'CG', 1)
W = V*Q*T

v, w, t = TestFunctions(W)
sol = Function(W)
pr = Constant(10)
gr = Constant(1000)
u, omega, temp = split(sol)

sol.project(Expression(("0", "0", "0", "x[0]")))


F = inner(grad(u), grad(v))*dx  + \
    dot(as_vector((-omega.dx(1), omega.dx(0))), v)*dx + \
    dot(grad(omega), grad(w))*dx + div(as_vector((u[0]*omega, u[1]*omega)))*w*dx - \
    gr*temp.dx(0)*w*dx + \
    dot(grad(temp), grad(t))*dx + pr*div(as_vector((u[0]*temp, u[1]*temp)))*t*dx


noslip = DirichletBC(W.sub(0), Constant((0, 0)), (1, 2, 3))

lid = DirichletBC(W.sub(0), Constant((1, 0)), 4)

omega_bc = DirichletBC(W.sub(1), Constant(0), (1, 2, 3, 4))

t_bc_1 = DirichletBC(W.sub(2), Constant(0), 1)
t_bc_2 = DirichletBC(W.sub(2), Constant(1), 2)



bcs = [lid, noslip, omega_bc, t_bc_1, t_bc_2]
solve(F == 0, sol, bcs=bcs,
      solver_parameters={'snes_monitor': True,
                         'pc_type': 'lu',
                         'ksp_monitor': True,})

u, omega, t = sol.split()

File("u.vtu") << u
File("omega.vtu") << omega
File("t.vtu") << t
