#include "case.h"
#include <iostream>

#include "omp.h"

namespace paraFEM {

FemCase::FemCase(std::vector<ElementPtr> elements)
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

void FemCase::makeStep(double h, double externalFactor)
{
    time += h;
    auto node_cp = this->get_nodes();

// for with iterator
#pragma omp parallel for
    for(int i = 0; i < elements.size(); i++)
    {
        auto element = elements[i];
        element->makeStep(h);  // computing the internal forces
    }

// for with iterator
#pragma omp parallel for
    for(int j = 0; j < node_cp.size(); j++)
    {
        auto node = node_cp[j];
        node->solveEquilibrium(h, externalFactor);
        node->internalForce.setZero();
    }
    
}

NodeVec FemCase::get_nodes()
{
    return NodeVec (nodes.begin(), nodes.end());
}


}