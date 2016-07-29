import paraFEM
import paraEigen as eigen

mat = paraFEM.TrussMaterial(1000)
mat.rho = 1
mat.d_structural = 0.1
mat.d_velocity = 5

node1 = paraFEM.Node(-1, 0, 0)
node2 = paraFEM.Node(0, 0, 0)
node3 = paraFEM.Node(1, 0, 0)

node1.fixed = eigen.vector3(0, 0, 0)
node3.fixed = eigen.vector3(0, 0, 0)
node2.add_external_force(eigen.vector3(0, 1, 0))

truss1 = paraFEM.Truss([node1, node2], mat)
truss2 = paraFEM.Truss([node2, node3], mat)

case = paraFEM.Case([truss1, truss2])

writer = paraFEM.vtkWriter("/tmp/paraFEM/test")

for i in range(1000):
    case.explicitStep(0.01)
    if (i % 10) == 0:
    	writer.writeCase(case, 0.3)
