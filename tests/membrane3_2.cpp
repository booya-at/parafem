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
    double stepsize = 0.005;
    int iterations = 500;
    int num_export = 100;
    
    // MATERIAL
    double E = 1000;
    double nue = 0.3;
    double structural_damping = 1. / E;
    double velocity_damping = 0.5;
    double rho = 0.1;
    
    // FORCE
    double pressure = 10;
    ////////////////////////////////////////////



    // DECLARATION
    parafem::NodeVec grid;
    parafem::ElementVec elements;
    parafem::Membrane3Ptr m1;
    parafem::Membrane3Ptr m2;

    int pos1;
    int pos2;
    int pos3;
    int pos4;

    parafem::MembraneMaterialPtr mat (new parafem::MembraneMaterial(E, nue));
    mat->d_structural = structural_damping;
    mat->rho = rho;
    mat->d_velocity = velocity_damping;

    // NODES
    for (int x=0; x < num_nodes; x++)
    {
        for (int y=0; y < num_nodes; y++)
        {
            grid.push_back(std::make_shared<parafem::Node>(x, y, 0));
            if (x==0 or x==num_nodes-1 or y==0 or y == num_nodes-1)
                grid.back()->fixed << 0, 0, 0;
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
            m1 = std::make_shared<parafem::Membrane3>(
                    parafem::NodeVec{grid[pos4], grid[pos3], grid[pos1]}, mat);
            m2 = std::make_shared<parafem::Membrane3>(
                    parafem::NodeVec{grid[pos2], grid[pos1], grid[pos3]}, mat);
            m1->pressure = pressure;
            m2->pressure = pressure;
            elements.push_back(m1);
            elements.push_back(m2);
        }
    }

    // CASE
    parafem::FemCasePtr c1 (new parafem::FemCase(elements));

    // WRITER
    parafem::VtkWriter writer = parafem::VtkWriter("/tmp/parafem/membrane3_2/output");

    // LOOP
    for (int i=0; i<iterations; i++)
    {
        c1->explicit_step(stepsize);
        if (i % int(iterations / num_export) == 0)
        {
            writer.writeCase(c1);
            cout << "time: "<< c1->time << " of " << iterations * stepsize << endl;
        }
    }
}
