U0 = FiniteElement("CG", interval, 1)
U1 = FiniteElement("DG", interval, 0)
V0 = FiniteElement("CG", interval, 1)
V1 = FiniteElement("DG", interval, 0)

W1_a = OuterProductElement(U0, V1)
W1_b = OuterProductElement(U1, V0)

W1 = W1_a + W1_b

RT1 = HDiv(full)
N1 = HCurl(full)
