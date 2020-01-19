#include <iostream>
#include <vector>

#include "node.h"
#include "element.h"
#include "case.h"
#include "vtkWriter.h"


int main(int argc, char **argv) {
    parafem::NodePtr n1 (new parafem::Node(0, 0, 0));
    parafem::NodePtr n2 (new parafem::Node(1, 0, 0));
    parafem::NodePtr n3 (new parafem::Node(1, 1, 0));
    parafem::NodePtr n4 (new parafem::Node(0, 1, 0));
    n1->fixed = parafem::Vector3(0,0,0);
    n4->fixed = parafem::Vector3(0,1,0);
    n2->fixed = parafem::Vector3(0,1,0);
    n3->fixed = parafem::Vector3(0,1,0);
    std::shared_ptr<parafem::MembraneMaterial> mat (new parafem::MembraneMaterial(100, 0.3));
    mat->d_velocity = 0.1;
    parafem::Membrane3Ptr m1 (new parafem::Membrane3(parafem::NodeVec{n4, n1, n3}, mat));
    parafem::Membrane3Ptr m2 (new parafem::Membrane3(parafem::NodeVec{n2, n3, n1}, mat));
    parafem::FemCasePtr c1 (new parafem::FemCase(parafem::ElementVec{m1, m2}));
    parafem::VtkWriter writer = parafem::VtkWriter("/tmp/parafem/membrane3_1/output");

    for (int i=0; i<100; i++)
    {
        n2->position += parafem::Vector3(0.01, 0, 0);
        n3->position += parafem::Vector3(0.01, 0, 0);
        n2->velocity.x() = 0.01 / 0.001;
        n3->velocity.x() = 0.01 / 0.001;
        c1->explicit_step(0.001);
        if (i % 10 == 0)
            writer.writeCase(c1, 0.4);
    }
}
