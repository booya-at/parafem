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
    int iterations = 100000;
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

    parafem::NodePtr n1 (new parafem::Node(0, 0, 0));
    parafem::NodePtr n2 (new parafem::Node(1, 0, 0));
    parafem::NodePtr n3 (new parafem::Node(1, 1, 0));
    parafem::NodePtr n4 (new parafem::Node(0, 1, 0));

    n1->fixed << 0, 0, 0;
    n2->fixed << 0, 0, 0;
    n3->fixed << 1, 0, 0;
    n4->fixed << 1, 0, 0;
    
    n3->external_force << force/2, 0, 0;
    n4->external_force << force/2, 0, 0;
  
    std::shared_ptr<parafem::MembraneMaterial> mat (new parafem::MembraneMaterial(E, nue));
    mat->d_structural = structural_damping;
    mat->d_velocity = velocity_damping;

    parafem::Membrane4Ptr m1 (new parafem::Membrane4(parafem::NodeVec{n1, n2, n3, n4}, mat, reduced_integration));
    parafem::FemCasePtr c1 (new parafem::FemCase(parafem::ElementVec{m1}));

    parafem::VtkWriter writer = parafem::VtkWriter("/tmp/parafem/membrane4_3/output");

    for (int i=0; i<iterations; i++)
    {
        c1->explicit_step(stepsize);
        if (i % int(iterations / num_export) == 0)
            writer.writeCase(c1, 0.3);
    }
    cout << "spannung in x richtung: " << m1->get_stress().x() << endl;
    cout << "spannung in y richtung: " << m1->get_stress().y() << endl;
    cout << "verschiebung in x richtung: " << n3->position.x()-1 << endl;
}
