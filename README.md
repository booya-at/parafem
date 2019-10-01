# paraFEM

The goal of this project is to have a working fem code supporting membrane and truss elements. The theory is taken from the paper [`EXPLICIT AlGORITHMS FOR THE NONLINEAR DYNAMICS OF SHELLS`][ResearchPaper]

[ResearchPaper]: https://s3.amazonaws.com/academia.edu.documents/46423439/0045-7825_2884_2990026-420160612-12328-se460b.pdf?response-content-disposition=inline%3B%20filename%3DExplicit_algorithms_for_the_nonlinear_dy.pdf&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAIWOWYYGZ2Y53UL3A%2F20190930%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20190930T182236Z&X-Amz-Expires=3600&X-Amz-SignedHeaders=host&X-Amz-Signature=04ebbc85c0c783c09f590ac0753ab36a420506205082f774879e8f72b27e3fb6

## Tasks/Todos

- [ ] A working `c++` source
- [ ] Comparison of elements
- [ ] Test framework
- [ ] JSON inputs
- [ ] Python binding with pybind11
- [ ] Make faster (fixed size arrays, ...)
- [ ] line-forces can be computed with a simplified geometry.  
  - Using a __implicit integration__ will make the computation faster.

## Tests and use-cases

### Pillow-test

Compare [klatsch-stange](https://pinsundmehr.de/wp-content/uploads/2017/08/klatschstangen-pum-01-bertrandt-ingenieur-dienstleister-02.jpg)

<img src="./images/pillow.png" alt="result" width="600"/>

### Deformation

Deformation of glider due to aerodynamic forces acting on membrane-elements

<img src="./images/glider.png" alt="result" width="600"/>

### Computing line-forces

Computing the forces having an impact on line-geometry and line-lengths

<img src="./images/line-forces.png" alt="result" width="600"/>

## Python bindings

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

## License
GPLv3 (see [`LICENSE`](LICENSE)]
