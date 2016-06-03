import paraFEM_py as paraFEM

mat = paraFEM.TrussMaterial(1000)
node1 = paraFEM.Node(0,0,0)
node2 = paraFEM.Node(0,0,1)
truss = paraFEM.Truss([node1, node2], mat)