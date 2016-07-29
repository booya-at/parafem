#include <iostream>
#include <vector>

#include "node.h"
#include "element.h"
#include "case.h"
#include "vtkWriter.h"


int main(int argc, char **argv) {
    paraFEM::NodePtr n1 (new paraFEM::Node(-1, 0, 0));
    paraFEM::NodePtr n2 (new paraFEM::Node(0, 0, 0));
    paraFEM::NodePtr n3 (new paraFEM::Node(1, 0, 0));
    n1->fixed = paraFEM::Vector3(0,0,0);
    n3->fixed = paraFEM::Vector3(0,0,0);
    n2->externalForce.y() -= 1;
    std::shared_ptr<paraFEM::TrussMaterial> mat (new paraFEM::TrussMaterial(1000));
    mat->d_structural = 0.00;
    paraFEM::TrussPtr t1 (new paraFEM::Truss(paraFEM::NodeVec{n1, n2}, mat));
    paraFEM::TrussPtr t2 (new paraFEM::Truss(paraFEM::NodeVec{n2, n3}, mat));
    paraFEM::FemCasePtr c1 (new paraFEM::FemCase(paraFEM::ElementVec{t1, t2}));
    paraFEM::VtkWriter writer = paraFEM::VtkWriter("/tmp/mytest/atest");
    for (int i=0; i<1000; i++)
    {
        c1->explicitStep(0.01);
        if (i % 10 == 0)
            writer.writeCase(c1);
        
    }
}
