#include <iostream>
#include <vector>

#include "node.h"
#include "element.h"
#include "case.h"
#include "vtkWriter.h"


int main(int argc, char **argv) {
    paraFEM::NodePtr n1 (new paraFEM::Node(0, 0, 0));
    paraFEM::NodePtr n2 (new paraFEM::Node(1, 0, 0));
    paraFEM::NodePtr n3 (new paraFEM::Node(0, 1, 0));
    n2->fixed = paraFEM::Vector3(0,0,0);
    n3->fixed = paraFEM::Vector3(0,0,0);
    n1->externalForce.y() -= 1;
    std::shared_ptr<paraFEM::MembraneMaterial> mat (new paraFEM::MembraneMaterial(100, 0.3));
    mat->d_structural = 1;
    paraFEM::Membrane3Ptr m1 (new paraFEM::Membrane3(paraFEM::NodeVec{n1, n2, n3}, mat));
    paraFEM::FemCasePtr c1 (new paraFEM::FemCase(paraFEM::ElementVec{m1}));
    c1->d_velocity = 0.1;
    VtkWriter writer = VtkWriter("/tmp/paraFEM/membrane3_1/output");

    for (int i=0; i<10000; i++)
    {
        c1->makeStep(0.0005);
        if (i % 30 == 0)
            writer.writeCase(c1, 0.4);
    }
}
