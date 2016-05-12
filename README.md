# paraFEM

The goal of this project is to have a working fem code supporting membrane and truss elements. The theory is taken from the paper `EXPLICIT AlGORITHMS FOR THE NONLINEAR DYNAMICS OF SHELLS`

## tasks

 - a working c++ source
 - test framework
 - json inputs
 - python binding with pybind11
 - make faster (fixed size arrays, ...)

 
## python bindings:

```python
import paraFEM
nodes = []
elements = []

mat1 = paraFEM.material.ripstop
mat2 = paraFEM.material.lines

for point in mesh.points:
    node = paraFEM.Node(*point)
    node.setExternalForce(0, 100, 10)
    node.setFixed(0, 0, 0)
    nodes.append(node)

for node_nr in mesh.elements:
    element = paraFEM.Element(*node_nr)
    element.setPressure(10, update=False)
    element.material = mat1
    element.writeStress = "scalar"  # or vector
    element.groupTag = "groupa"
    element.tag = "10"
    elements.append(element)
    
for hinge_nodes in mesh.hinges:
    hinge = paraFEM.hinge(*hinge)
    elemnts.append(hinge)

    
case = paraFEM.case(elements)
case.store(100)
case.maxU(0.0001)
case.run(0, 10)

paraFEM.vtk.storeCase(case, "filename")
```