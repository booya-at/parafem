#include <iostream>
#include <vector>

#include "node.h"
#include "element.h"
#include "case.h"
#include "vtkWriter.h"


int main(int argc, char **argv) {

    ///////////////----INPUT----////////////////
    // MESH
    int num_nodes = 40; // num_nodes x num_nodes

    // SIM
    double stepsize = 0.001;
    int interations = 5000;
    int num_export = 100;
    
    // MATERIAL
    double E = 1000;
    double nue = 0.4;
    double structural_damping = 1. / E;
    double velocity_damping = 0.2;
    double rho = 0.1;
    
    // FORCE
    double pressure = 10;
    ////////////////////////////////////////////


    //INTEGRATION
    bool reduced_integration = 1;


    // DECLARATION
    paraFEM::NodeVec grid;
    paraFEM::ElementVec elements;
    paraFEM::Membrane4Ptr m1;

    int pos1;
    int pos2;
    int pos3;
    int pos4;

    paraFEM::MembraneMaterialPtr mat (new paraFEM::MembraneMaterial(E, nue));
    mat->d_structural = structural_damping;
    mat->rho = rho;

    // NODES
    for (int x=0; x < num_nodes; x++)
    {
        for (int y=0; y < num_nodes; y++)
        {
            grid.push_back(std::make_shared<paraFEM::Node>(x, y, 0));
            if (x==0 or x==num_nodes-1 or y==0 or y == num_nodes-1)
                grid.back()->fixed << 1, 1, 0;
        }
    }

    // ELEMENTS
    for (int ex=0; ex < (num_nodes -1); ex++)
    {
        for (int ey=0; ey < num_nodes -1; ey++)
        {
            pos1 = ex * num_nodes + ey;
            pos2 = pos1 + 1;
            pos3 = pos2 + num_nodes;
            pos4 = pos3 -1;
            m1 = std::make_shared<paraFEM::Membrane4>(
                    paraFEM::NodeVec{grid[pos1], grid[pos4], grid[pos3], grid[pos2]},
                    mat, reduced_integration);
            m1->coordSys = paraFEM::CoordSys(m1->coordSys.n, paraFEM::Vector3(1, 0, 0));
            m1->pressure = pressure;
            elements.push_back(m1);
        }
    }

    // CASE
    paraFEM::FemCasePtr c1 (new paraFEM::FemCase(elements));
    c1->d_velocity = velocity_damping;

    // WRITER
    VtkWriter writer = VtkWriter("/tmp/paraFEM/membrane4_2/int0output");

    // LOOP;
    for (int i=0; i<interations; i++)
    {
        c1->makeStep(stepsize);
        if (i % int(interations / num_export) == 0)
        {
            writer.writeCase(c1);
            cout << "time: "<< c1->time << " of " << interations * stepsize << endl;
        }
    }
}
