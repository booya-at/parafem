#include <iostream>
#include <vector>

#include "node.h"
#include "element.h"
#include "case.h"
#include "vtkWriter.h"
#include <memory>


int main(int argc, char **argv) {
    paraFEM::NodeVec nodes;
    paraFEM::ElementVec elements;
    paraFEM::TrussPtr truss;
    std::shared_ptr<paraFEM::TrussMaterial> mat (new paraFEM::TrussMaterial(5000));
    mat->d_structural = 10;
    mat->d_velocity = 0.1;
    for (int i = 0; i < 100; i++)
    {
        nodes.push_back(std::make_shared<paraFEM::Node>(double(i), 0, 0));
        nodes.back()->externalForce.y() = -1;
        if (i != 0)
        {
            truss = make_shared<paraFEM::Truss>(paraFEM::NodeVec{nodes[i-1], nodes[i]}, mat);
            elements.push_back(truss);
        }
    }
    nodes[50]->externalForce.z() = 100;
    nodes[0]->fixed = paraFEM::Vector3(0, 0, 0);
    nodes.back()->fixed = paraFEM::Vector3(0, 0, 0);
    paraFEM::FemCasePtr c (new paraFEM::FemCase(elements));

    paraFEM::VtkWriter writer = paraFEM::VtkWriter("/tmp/paraFEM/truss2_");
    for (int i=0; i<2000; i++)
    {
        c->makeStep(0.01);
        if (i % 20 == 0)
            writer.writeCase(c);
        
    }
}
