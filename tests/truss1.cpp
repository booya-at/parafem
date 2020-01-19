#include <iostream>
#include <vector>

#include "node.h"
#include "element.h"
#include "case.h"
#include "vtkWriter.h"


int main(int argc, char **argv) {
    parafem::NodePtr n1 (new parafem::Node(-1, 0, 0));
    parafem::NodePtr n2 (new parafem::Node(0, 0, 0));
    parafem::NodePtr n3 (new parafem::Node(1, 0, 0));
    n1->fixed = parafem::Vector3(0,0,0);
    n3->fixed = parafem::Vector3(0,0,0);
    n2->external_force.y() -= 1;
    std::shared_ptr<parafem::TrussMaterial> mat (new parafem::TrussMaterial(1000));
    mat->d_structural = 0.00;
    parafem::TrussPtr t1 (new parafem::Truss(parafem::NodeVec{n1, n2}, mat));
    parafem::TrussPtr t2 (new parafem::Truss(parafem::NodeVec{n2, n3}, mat));
    parafem::FemCasePtr c1 (new parafem::FemCase(parafem::ElementVec{t1, t2}));
    parafem::VtkWriter writer = parafem::VtkWriter("/tmp/parafem/truss1/output");
    for (int i=0; i<1000; i++)
    {
        c1->explicit_step(0.01);
        if (i % 10 == 0)
            writer.writeCase(c1);
        
    }
}
