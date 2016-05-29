#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include "node.h"
#include "element.h"
#include "case.h"
#include "material.h"

namespace py = pybind11;

void init_paraFEM(py::module &m){
    
    py::class_<paraFEM::Node> (m, "Node")
        .def(py::init<double, double, double>());
    
    py::class_<paraFEM::TrussMaterial>(m, "TrussMaterial");
    
    py::class_<paraFEM::MembraneMaterial>(m, "MembraneMaterial");
        
    py::class_<paraFEM::Truss> (m, "Truss")
        .def(py::init<std::vector<paraFEM::NodePtr>, paraFEM::TrussMaterialPtr>());
    
    py::class_<paraFEM::Membrane3> (m, "Membrane3")
        .def(py::init<std::vector<paraFEM::NodePtr>, paraFEM::MembraneMaterialPtr>());
      
    py::class_<paraFEM::Membrane4> (m, "Membrane4")
        .def(py::init<std::vector<paraFEM::NodePtr>, paraFEM::MembraneMaterialPtr>());
    
}


PYBIND11_PLUGIN(paraFEM_py){
    py::module m("paraFEM_py");
    init_paraFEM(m);
    return m.ptr();
};