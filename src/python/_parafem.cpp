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

PYBIND11_MODULE(_parafem, m){
    #ifdef BUILD_VTK_WRITER
        py::class_<parafem::VtkWriter, std::shared_ptr<parafem::VtkWriter>>(m, "vtkWriter")
            .def(py::init<const char*>())
            .def("writeCase", &parafem::VtkWriter::writeCase);
    #endif

    py::class_<parafem::Node, parafem::NodePtr> (m, "Node")
        .def(py::init<double, double, double>())
        .def_readonly("position", &parafem::Node::position)
        .def_readonly("velocity", &parafem::Node::velocity)
        .def_readonly("acceleration", &parafem::Node::acceleration)
        .def_readwrite("fixed", &parafem::Node::fixed)
        .def_readwrite("mass_influence", &parafem::Node::mass_influence)
        .def("add_external_force", &parafem::Node::add_external_force);

    py::class_<parafem::Material, parafem::MaterialPtr>(m, "Material")
        .def_readwrite("rho", &parafem::Material::rho)
        .def_readwrite("d_structural", &parafem::Material::d_structural)
        .def_readwrite("d_velocity", &parafem::Material::d_velocity)
        .def_readwrite("elasticity", &parafem::Material::elasticity);
    
    py::class_<parafem::TrussMaterial, parafem::TrussMaterialPtr, parafem::Material>(m, "TrussMaterial")
        .def(py::init<double>());
    
    py::class_<parafem::MembraneMaterial, parafem::MembraneMaterialPtr, parafem::Material>(m, "MembraneMaterial")
        .def(py::init<double, double>());

    py::class_<parafem::Element, parafem::ElementPtr>(m, "__Element")
        .def("get_stress", &parafem::Element::get_stress)
        .def_readonly("nodes", &parafem::Element::nodes);

    py::class_<parafem::Membrane, parafem::MembranePtr, parafem::Element>(m, "__Membrane")
        .def_readwrite("pressure", &parafem::Membrane::pressure);
        
    py::class_<parafem::Truss, parafem::TrussPtr, parafem::Element> (m, "Truss")
        .def(py::init<std::vector<parafem::NodePtr>, parafem::TrussMaterialPtr>());
    
    py::class_<parafem::Membrane3, parafem::Membrane3Ptr, parafem::Membrane> (m, "Membrane3")
        .def(py::init<std::vector<parafem::NodePtr>, parafem::MembraneMaterialPtr>());
      
    py::class_<parafem::Membrane4, parafem::Membrane4Ptr, parafem::Membrane> (m, "Membrane4")
        .def(py::init<std::vector<parafem::NodePtr>, parafem::MembraneMaterialPtr>())
        .def(py::init<std::vector<parafem::NodePtr>, parafem::MembraneMaterialPtr, bool>());

    py::class_<parafem::FemCase, parafem::FemCasePtr> (m, "Case")
        .def(py::init<std::vector<parafem::ElementPtr>>())
        .def("explicit_step", &parafem::FemCase::explicit_step, "make one explicit step",
            py::arg("h") = 0.0001, py::arg("external_factor") = 1)
        .def("get_explicit_max_time_step", &parafem::FemCase::get_explicit_max_time_step, "cfl - value")
        .def("get_max_velocity", &parafem::FemCase::get_max_velocity, "get Max node velocity");
};