#include <iostream>
#include <vector>

#include "node.h"
#include "element.h"
#include "case.h"
#include "vtkWriter.h"
#include <memory>


int main(int argc, char **argv) {
    parafem::NodeVec nodes;
    parafem::ElementVec elements;
    parafem::TrussPtr truss;
    std::shared_ptr<parafem::TrussMaterial> mat (new parafem::TrussMaterial(5000));
    mat->d_structural = 10;
    mat->d_velocity = 0.1;
    for (int i = 0; i < 100; i++)
    {
        nodes.push_back(std::make_shared<parafem::Node>(double(i), 0, 0));
        nodes.back()->external_force.y() = -1;
        if (i != 0)
        {
            truss = make_shared<parafem::Truss>(parafem::NodeVec{nodes[i-1], nodes[i]}, mat);
            elements.push_back(truss);
        }
    }
    nodes[50]->external_force.z() = 100;
    nodes[0]->fixed = parafem::Vector3(0, 0, 0);
    nodes.back()->fixed = parafem::Vector3(0, 0, 0);
    parafem::FemCasePtr c (new parafem::FemCase(elements));

    parafem::VtkWriter writer = parafem::VtkWriter("/tmp/parafem/truss2/output");
    for (int i=0; i<2000; i++)
    {
        c->explicit_step(0.01);
        if (i % 20 == 0)
            writer.writeCase(c);
        
    }
}
