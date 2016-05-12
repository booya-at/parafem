#include "element.h"

#include <iostream>

namespace paraFEM {

Truss::Truss(std::vector<NodePtr> points, std::shared_ptr<TrussMaterial> mat)
{
    nodes = points;
    stress = 0;
    length = (nodes[0]->position - nodes[1]->position).norm();
    material = mat;
    nodes[0]->massInfluence += 0.5 * material->rho * length;
    nodes[1]->massInfluence += 0.5 * material->rho * length;
}


void Truss::makeStep(double h)
{
    // update the coordinate system
    Vector3 p0 = nodes[0]->position;
    Vector3 p1 = nodes[1]->position;
    double new_length = (p1 - p0).norm();
    Vector3 tangent = p1 - p0;
    tangent.normalize();
    
    // velocity vectors
    Vector3 u0 = nodes[0]->velocity;
    Vector3 u1 = nodes[1]->velocity;
    // strainrate
    double strain_rate =  2 / (length + new_length) * (u1- u0).dot(tangent);
    // stressrate
    double stress_rate = material->elasticity * strain_rate;
    length = new_length;
    // stress
    stress += h * stress_rate;

    // internal force with virtual power
    nodes[1]->internalForce -= tangent * (stress + strain_rate * material->d_structural);
    nodes[0]->internalForce += tangent * (stress + strain_rate * material->d_structural);
}

}