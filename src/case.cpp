#include "case.h"
#include <iostream>

#include "omp.h"

namespace paraFEM {

FemCase::FemCase(std::vector<ElementPtr> elements, 
                 std::vector<PositionConstraintPtr> posConstraints,
                 std::vector<ForceConstraintPtr> forceConstraints): FemCase::FemCase(elements)
{
    this->posConstraints = posConstraints;
    this->forceConstraints = forceConstraints;
}

FemCase::FemCase(std::vector<ElementPtr> elements)
{
    for (auto element: elements)
    {
        for (auto node: element->nodes)
            this->nodes.insert(node);
        this->elements.push_back(element);
    }
    std::cout << "maxTimeStep= " << std::get<0>(this->getExplicitMaxTimeStep()) << std::endl;
    int count = 0;
    for (auto node: nodes)
        node->nr = count ++;
}

std::tuple<double, ElementPtr> FemCase::getExplicitMaxTimeStep()
{

    double maxTimeStep = 1;
    double elTimeStep;
    ElementPtr badElement;
    for (auto element: this->elements)
    {
        elTimeStep = element->characteristicLength / element->getMaterial()->waveSpeed();
        if (elTimeStep < maxTimeStep and elTimeStep != 0)
        {
            badElement = element;
            maxTimeStep = elTimeStep;
        }
    }
    return std::make_tuple(maxTimeStep, badElement);
}


void FemCase::explicitStep(double h, double externalFactor)
{
    time += h;
    auto node_cp = this->get_nodes();

// for with iterator
#pragma omp parallel for
    for(int i = 0; i < elements.size(); i++)
    {
        auto element = elements[i];
        element->explicitStep(h);  // computing the internal forces
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

// void FemCase::implicitStep(double h, double externalFactor)
// {
//     time += h;
//     auto node_cp = this->get_nodes();
//     Eigen::Triplet<double> Kt, Dt, Mt;
//     Eigen::MatrixXd K, D, M;

// #pragma omp parallel for
//     for(int j = 0; j < elements.size(); j++)
//     {
//         auto element = elements[j];
//         element->implicitStep(h, Kt, Dt, Mt);
//     }
//     // create global matrices from triplets

//     // solve equilibrium equation

//     // set new positions

//     // integrate force
    
// }

NodeVec FemCase::get_nodes()
{
    return NodeVec (nodes.begin(), nodes.end());
}


}