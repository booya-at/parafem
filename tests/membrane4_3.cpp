#include <iostream>
#include <vector>

#include "node.h"
#include "element.h"
#include "case.h"
#include "vtkWriter.h"


int main(int argc, char **argv) {
    
    
    ///////////////----INPUT----////////////////

    // SIM
    double stepsize = 0.0001;
    int interations = 100000;
    int num_export = 200;
    
    // MATERIAL
    double E = 100;
    double nue = 0.0;
    double structural_damping = 0.0;
    double velocity_damping = 0;
    
    // FORCE
    double force = 100;
    ////////////////////////////////////////////
    
    //INTEGRATION
    bool reduced_integration = true;

    paraFEM::NodePtr n1 (new paraFEM::Node(0, 0, 0));
    paraFEM::NodePtr n2 (new paraFEM::Node(1, 0, 0));
    paraFEM::NodePtr n3 (new paraFEM::Node(1, 1, 0));
    paraFEM::NodePtr n4 (new paraFEM::Node(0, 1, 0));

    n1->fixed << 0, 0, 0;
    n2->fixed << 0, 0, 0;
    n3->fixed << 1, 0, 0;
    n4->fixed << 1, 0, 0;
    
    n3->externalForce << force/2, 0, 0;
    n4->externalForce << force/2, 0, 0;
  
    std::shared_ptr<paraFEM::MembraneMaterial> mat (new paraFEM::MembraneMaterial(E, nue));
    mat->d_structural = structural_damping;

    paraFEM::Membrane4Ptr m1 (new paraFEM::Membrane4(paraFEM::NodeVec{n1, n2, n3, n4}, mat, reduced_integration));
    paraFEM::FemCasePtr c1 (new paraFEM::FemCase(paraFEM::ElementVec{m1}));

    c1->d_velocity = velocity_damping;
    VtkWriter writer = VtkWriter("/tmp/paraFEM/membrane4_3/output");

    for (int i=0; i<interations; i++)
    {
        c1->makeStep(stepsize);
        if (i % int(interations / num_export) == 0)
            writer.writeCase(c1, 0.3);
    }
    cout << "spannung in x richtung: " << m1->getStress().x() << endl;
    cout << "spannung in y richtung: " << m1->getStress().y() << endl;
    cout << "verschiebung in x richtung: " << n3->position.x()-1 << endl;
}
