#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/eigen.h>

#include "node.h"
#include "element.h"
#include "case.h"
#include "material.h"

#ifdef BUILD_VTK_WRITER
    #include "vtkWriter.h"
#endif


// PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);


namespace py = pybind11;

PYBIND11_MODULE(_paraFEM, m){
    #ifdef BUILD_VTK_WRITER
        py::class_<paraFEM::VtkWriter, std::shared_ptr<paraFEM::VtkWriter>>(m, "vtkWriter")
            .def(py::init<const char*>())
            .def("writeCase", &paraFEM::VtkWriter::writeCase);
    #endif

    py::class_<paraFEM::Node, paraFEM::NodePtr> (m, "Node")
        .def(py::init<double, double, double>())
        .def_readonly("position", &paraFEM::Node::position)
        .def_readonly("velocity", &paraFEM::Node::velocity)
        .def_readonly("acceleration", &paraFEM::Node::acceleration)
        .def_readwrite("fixed", &paraFEM::Node::fixed)
        .def_readwrite("massInfluence", &paraFEM::Node::massInfluence)
        .def("add_external_force", &paraFEM::Node::add_external_force);

    py::class_<paraFEM::Material, paraFEM::MaterialPtr>(m, "Material")
        .def_readwrite("rho", &paraFEM::Material::rho)
        .def_readwrite("d_structural", &paraFEM::Material::d_structural)
        .def_readwrite("d_velocity", &paraFEM::Material::d_velocity)
        .def_readwrite("elasticity", &paraFEM::Material::elasticity);
    
    py::class_<paraFEM::TrussMaterial, paraFEM::TrussMaterialPtr, paraFEM::Material>(m, "TrussMaterial")
        .def(py::init<double>());
    
    py::class_<paraFEM::MembraneMaterial, paraFEM::MembraneMaterialPtr, paraFEM::Material>(m, "MembraneMaterial")
        .def(py::init<double, double>());

    py::class_<paraFEM::Element, paraFEM::ElementPtr>(m, "__Element")
        .def("getStress", &paraFEM::Element::getStress)
        .def_readonly("nodes", &paraFEM::Element::nodes);

    py::class_<paraFEM::Membrane, paraFEM::MembranePtr, paraFEM::Element>(m, "__Membrane")
        .def_readwrite("pressure", &paraFEM::Membrane::pressure);
        
    py::class_<paraFEM::Truss, paraFEM::TrussPtr, paraFEM::Element> (m, "Truss")
        .def(py::init<std::vector<paraFEM::NodePtr>, paraFEM::TrussMaterialPtr>());
    
    py::class_<paraFEM::Membrane3, paraFEM::Membrane3Ptr, paraFEM::Membrane> (m, "Membrane3")
        .def(py::init<std::vector<paraFEM::NodePtr>, paraFEM::MembraneMaterialPtr>());
      
    py::class_<paraFEM::Membrane4, paraFEM::Membrane4Ptr, paraFEM::Membrane> (m, "Membrane4")
        .def(py::init<std::vector<paraFEM::NodePtr>, paraFEM::MembraneMaterialPtr>())
        .def(py::init<std::vector<paraFEM::NodePtr>, paraFEM::MembraneMaterialPtr, bool>());

    py::class_<paraFEM::FemCase, paraFEM::FemCasePtr> (m, "Case")
        .def(py::init<std::vector<paraFEM::ElementPtr>>())
        .def("explicitStep", &paraFEM::FemCase::explicitStep, "make one explicit step",
            py::arg("h") = 0.0001, py::arg("externalFactor") = 1)
        .def("getExplicitMaxTimeStep", &paraFEM::FemCase::getExplicitMaxTimeStep, "cfl - value");
};