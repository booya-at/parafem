import parafem
import numpy as np

mat = parafem.TrussMaterial(1000)
mat.rho = 1
mat.d_structural = 0.0
mat.d_velocity = 1

node1 = parafem.Node(-1, 0, 0)
node2 = parafem.Node(0, 0, 0)
node3 = parafem.Node(1, 0, 0)

node1.fixed = np.array([0, 0, 0])
node3.fixed = np.array([0, 0, 0])
node2.add_external_force(np.array([0, 1, 0]))

truss1 = parafem.Truss([node1, node2], mat)
truss2 = parafem.Truss([node2, node3], mat)

case = parafem.Case([truss1, truss2])

writer = parafem.vtkWriter("/tmp/parafem/truss1_py/output")

for i in range(10000):
    case.explicit_step(0.01)
    if (i % 10) == 0:
        writer.writeCase(case, 0.3)
