#include "element.h"

#include <iostream>

namespace paraFEM {

Truss::Truss(std::vector<NodePtr> points, std::shared_ptr<TrussMaterial> mat)
{
    nodes = points;
    stress = 0;
    tangent = (nodes[1]->position - nodes[0]->position);
    length = tangent.norm();
    tangent.normalize();
    material = mat;
    nodes[0]->massInfluence += 0.5 * material->rho * length;
    nodes[1]->massInfluence += 0.5 * material->rho * length;
}


void Truss::makeStep(double h)
{
    // update the coordinate system
    tangent =  nodes[1]->position - nodes[0]->position;
    double new_length = tangent.norm();
    tangent.normalize();
    
    // strainrate
    double strain_rate =  2 / (length + new_length) * 
                          (nodes[0]->velocity - nodes[1]->velocity).dot(tangent);
    // stressrate
    double stress_rate = material->elasticity * strain_rate;
    length = new_length;
    // stress
    stress += h * stress_rate;

    // internal force with virtual power
    nodes[1]->internalForce -= tangent * (stress + strain_rate * material->d_structural);
    nodes[0]->internalForce += tangent * (stress + strain_rate * material->d_structural);
    for (auto node: nodes)
        node->internalForce += node->velocity * material->d_velocity * new_length / nodes.size();
}

Vector3 Truss::getStress()
{
    return stress * tangent;
};

}