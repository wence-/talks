m = UnitIntervalMesh(5)

mesh = ExtrudedMesh(m, layers=5)

U0 = FiniteElement("CG", interval, 1)
U1 = FiniteElement("DG", interval, 0)
V0 = FiniteElement("CG", interval, 1)
V1 = FiniteElement("DG", interval, 0)

W0_elt = OuterProductElement(U0, V0)
W1_a = HDiv(OuterProductElement(U1, V0))
W1_b = HDiv(OuterProductElement(U0, V1))
W1_elt = W1_a + W1_b

W0 = FunctionSpace(mesh, W0_elt)
W1 = FunctionSpace(mesh, W1_elt)

W = W0*W1
sigma, u = TrialFunctions(W)
tau, v = TestFunctions(W)

L = assemble((sigma*tau - inner(rot(tau), u) + inner(rot(sigma), v) +
              div(u)*div(v))*dx)
