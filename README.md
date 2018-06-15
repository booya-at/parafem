# paraFEM

The goal of this project is to have a working fem code supporting membrane and truss elements. The theory is taken from the paper `EXPLICIT AlGORITHMS FOR THE NONLINEAR DYNAMICS OF SHELLS`

## tasks

 - a working c++ source
 - compareison of elements
 - test framework
 - json inputs
 - python binding with pybind11
 - make faster (fixed size arrays, ...)

## tests and use-cases

### pillow-test
compare [klatsch-stange](https://pinsundmehr.de/wp-content/uploads/2017/08/klatschstangen-pum-01-bertrandt-ingenieur-dienstleister-02.jpg)
<img src="./images/pillow.png" alt="result" width="400"/>

### deformation of glider due to aerodynamic forces acting on mebran-elements
<img src="./images/glider.png" alt="result" width="400"/>

### computing line-forces. (forces have impact on line-geometry and line-lengths)
<img src="./images/line-forces.png" alt="result" width="400"/>

## python bindings:

```python
import paraFEM
import numpy as np

mat = paraFEM.TrussMaterial(1000)
mat.rho = 1
mat.d_structural = 0.0
mat.d_velocity = 1

node1 = paraFEM.Node(-1, 0, 0)
node2 = paraFEM.Node(0, 0, 0)
node3 = paraFEM.Node(1, 0, 0)

node1.fixed = np.array([0, 0, 0])
node3.fixed = np.array([0, 0, 0])
node2.add_external_force(np.array([0, 1, 0]))

truss1 = paraFEM.Truss([node1, node2], mat)
truss2 = paraFEM.Truss([node2, node3], mat)

case = paraFEM.Case([truss1, truss2])

writer = paraFEM.vtkWriter("/tmp/paraFEM/truss1_py/output")

for i in range(10000):
    case.explicitStep(0.01)
    if (i % 10) == 0:
        writer.writeCase(case, 0.3)
```


## TODO:

- line-forces can be computed with a simplified geometry. Using a __implicit integration__ will make the computation faster.