import paraFEM
import paraEigen
# ._____._____.
# |     |
# |     |
# ._____._____.

n1 = paraFEM.Node(-1, 0, 0)
n2 = paraFEM.Node(-1, -1, 0)
n3 = paraFEM.Node(-1, -2, 0)
n4 = paraFEM.Node(1, -2, 0)
n5 = paraFEM.Node(1, -1, 0)
n6 = paraFEM.Node(1, 0, 0)

n1.fixed = paraEigen.vector3(0, 0, 0)
n6.fixed = paraEigen.vector3(0, 0, 0)

n3.add_external_force(paraEigen.vector3(0, -100, 0))
n4.add_external_force(paraEigen.vector3(0, -0, 0))

membrane_mat = paraFEM.MembraneMaterial(100, 0.3)
membrane_mat.d_velocity = 0.0
membrane_mat.rho = 0.001

truss_mat = paraFEM.TrussMaterial(10000)
truss_mat.d_velocity = 0.000
truss_mat.d_structural = 0
truss_mat.rho = 0.001

m1 = paraFEM.Membrane3([n2, n3, n4], membrane_mat)
m2 = paraFEM.Membrane3([n2, n5, n4], membrane_mat)
t1 = paraFEM.Truss([n1, n2], truss_mat)
t2 = paraFEM.Truss([n6, n5], truss_mat)

case = paraFEM.Case([m1, m2, t1])
writer = paraFEM.vtkWriter("/tmp/test_mbr_truss/bla")

for i in range(50000):
    case.makeStep(h=0.000001)
    if i % 1000 == 0:
        writer.writeCase(case, 0.1)
