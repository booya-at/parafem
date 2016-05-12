#include "case.h"
#include <iostream>

namespace paraFEM {

FemCase::FemCase(ElementVec elements)
{
    for (auto element: elements)
    {
        for (auto node: element->nodes)
            this->nodes.insert(node);
        this->elements.push_back(element);
    }
    int count = 0;
    for (auto node: nodes)
        node->nr = count ++;
}

void FemCase::makeStep(double h)
{
    time += h;

    // do in parallel
    for (auto element: elements)
    {
        element->makeStep(h);  // computing the internal forces
    }
    //do in parallel
    for (auto node: nodes)
    {
        node->internalForce += node->velocity * d_velocity;
        node->solveEquilibrium(h);
        node->internalForce.setZero();
    }
    
}


}