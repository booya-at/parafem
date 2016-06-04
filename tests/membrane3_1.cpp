#include <iostream>
#include <vector>

#include "node.h"
#include "element.h"
#include "case.h"
#include "vtkWriter.h"


int main(int argc, char **argv) {
    paraFEM::NodePtr n1 (new paraFEM::Node(0, 0, 0));
    paraFEM::NodePtr n2 (new paraFEM::Node(1, 0, 0));
    paraFEM::NodePtr n3 (new paraFEM::Node(1, 1, 0));
    paraFEM::NodePtr n4 (new paraFEM::Node(0, 1, 0));
    n1->fixed = paraFEM::Vector3(0,0,0);
    n4->fixed = paraFEM::Vector3(0,1,0);
    n2->fixed = paraFEM::Vector3(0,1,0);
    n3->fixed = paraFEM::Vector3(0,1,0);
    std::shared_ptr<paraFEM::MembraneMaterial> mat (new paraFEM::MembraneMaterial(100, 0.3));
    mat->d_velocity = 0.1;
    paraFEM::Membrane3Ptr m1 (new paraFEM::Membrane3(paraFEM::NodeVec{n4, n1, n3}, mat));
    paraFEM::Membrane3Ptr m2 (new paraFEM::Membrane3(paraFEM::NodeVec{n2, n3, n1}, mat));
    paraFEM::FemCasePtr c1 (new paraFEM::FemCase(paraFEM::ElementVec{m1, m2}));
    paraFEM::VtkWriter writer = paraFEM::VtkWriter("/tmp/paraFEM/membrane3_1/output");

    for (int i=0; i<100; i++)
    {
        n2->position += paraFEM::Vector3(0.01, 0, 0);
        n3->position += paraFEM::Vector3(0.01, 0, 0);
        n2->velocity.x() = 0.01 / 0.001;
        n3->velocity.x() = 0.01 / 0.001;
        c1->makeStep(0.001);
        if (i % 10 == 0)
            writer.writeCase(c1, 0.4);
    }
}
