import parafem
import numpy
# ._____._____.
# |     |
# |     |
# ._____._____.

n1 = parafem.Node(-1, 0, 0)
n2 = parafem.Node(-1, -1, 0)
n3 = parafem.Node(-1, -2, 0)
n4 = parafem.Node(1, -2, 0)
n5 = parafem.Node(1, -1, 0)
n6 = parafem.Node(1, 0, 0)

n7 = parafem.Node(0, -3, 0)
#m3: n7 n3 n4

n1.fixed = numpy.array([0, 0, 0])
n6.fixed = numpy.array([0, 0, 0])

#n3.add_external_force(numpy.array([0, -100, 0]))
#n4.add_external_force(numpy.array([0, -0, 0]))
n7.add_external_force(numpy.array([0, -100, 0]))

membrane_mat = parafem.MembraneMaterial(100, 0.3)
membrane_mat.d_velocity = 0.0
membrane_mat.rho = 0.001

truss_mat = parafem.TrussMaterial(10000)
truss_mat.d_velocity = 0.000
truss_mat.d_structural = 0
truss_mat.rho = 0.001

m1 = parafem.Membrane3([n2, n3, n4], membrane_mat)
m2 = parafem.Membrane3([n2, n5, n4], membrane_mat)
m3 = parafem.Membrane3([n7, n3, n4], membrane_mat)
t1 = parafem.Truss([n1, n2], truss_mat)
t2 = parafem.Truss([n6, n5], truss_mat)

case = parafem.Case([m1, m2, m3, t1, t2])
writer = parafem.vtkWriter("/tmp/test_mbr_truss/bla")

for i in range(50000):
    case.explicit_step(h=0.00001)
    if i % 1000 == 0:
        print(i)
        writer.writeCase(case, 0.1)
