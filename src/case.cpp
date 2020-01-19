#include "case.h"
#include <iostream>
#include <algorithm>    // std::max
#include "omp.h"

namespace parafem {

FemCase::FemCase(std::vector<ElementPtr> elements, 
                 std::vector<PositionConstraintPtr> pos_constraints,
                 std::vector<ForceConstraintPtr> force_constraints): FemCase::FemCase(elements)
{
    this->pos_constraints = pos_constraints;
    this->force_constraints = force_constraints;
}

FemCase::FemCase(std::vector<ElementPtr> elements)
{
    for (auto element: elements)
    {
        for (auto node: element->nodes)
            this->nodes.insert(node);
        this->elements.push_back(element);
    }
    std::cout << "maxTimeStep= " << std::get<0>(this->get_explicit_max_time_step()) << std::endl;
    int count = 0;
    for (auto node: nodes)
        node->nr = count ++;
}

std::tuple<double, ElementPtr> FemCase::get_explicit_max_time_step()
{

    double maxTimeStep = 1;
    double elTimeStep;
    ElementPtr badElement;
    for (auto element: this->elements)
    {
        elTimeStep = element->characteristicLength / element->get_material()->waveSpeed();
        if (elTimeStep < maxTimeStep and elTimeStep != 0)
        {
            badElement = element;
            maxTimeStep = elTimeStep;
        }
    }
    return std::make_tuple(maxTimeStep, badElement);
}

double FemCase::get_max_velocity() {
    double velocity = 0;
    for (auto node: this->get_nodes()) {
        velocity = std::max(velocity, node->velocity.norm());
    }
    return velocity;
}


void FemCase::explicit_step(double h, double external_factor)
{
    time += h;
    auto node_cp = this->get_nodes();

    #pragma omp parallel for
    for (int i=0; i<this->elements.size(); i++) {
        auto element = elements[i];
        element->geometry_step();
    }

    for(int i = 0; i < elements.size(); i++)
    {
        auto element = elements[i];
        element->explicit_step(h);  // computing the internal forces
    }

    #pragma omp parallel for
    for(int j = 0; j < node_cp.size(); j++)
    {
        auto node = node_cp[j];
        node->solveEquilibrium(h, external_factor);
        node->internal_force.setZero();
    }
    
}

void FemCase::apply_forces()
{
    std::vector<trip> Kt;

    #pragma omp parallel for
    for(int j = 0; j < elements.size(); j++)
    {
        auto element = elements[j];
        element->implicit_step(Kt);
    }

    // applying forces
    // apply boundary conditions
    // solve system -> u, residuals
    // apply forces
    // add residuals to nodes.

}

NodeVec FemCase::get_nodes()
{
    return NodeVec (nodes.begin(), nodes.end());
}


}