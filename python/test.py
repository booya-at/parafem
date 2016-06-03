import paraFEM_py as fem
import paraEigen as eigen

mat = fem.TrussMaterial(1000)
mat.rho = 1
mat.d_structural = 0.1
mat.d_velocity = 5

node1 = fem.Node(-1, 0, 0)
node2 = fem.Node(0, 0, 0)
node3 = fem.Node(1, 0, 0)

node1.fixed = eigen.vector3(0, 0, 0)
node3.fixed = eigen.vector3(0, 0, 0)
node2.add_external_force(eigen.vector3(0, 1, 0))

truss1 = fem.Truss([node1, node2], mat)
truss2 = fem.Truss([node2, node3], mat)

case = fem.Case([truss1, truss2])

for i in range(1000):
    case.makeStep(0.01)
    print(node2.acceleration)
