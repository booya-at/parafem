#include "element.h"
#include "assert.h"

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
    dViscous = length * material->waveSpeed() * material->rho * material->d_minMode;
    this->characteristicLength = length;
}

MaterialPtr Truss::getMaterial()
{
    return std::dynamic_pointer_cast<Material>(this->material);
}


void Truss::explicitStep(double h)
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
    double absStress = stress + stress_rate * material->d_structural + dViscous * strain_rate;
    nodes[1]->internalForce -= tangent * absStress;
    nodes[0]->internalForce += tangent * absStress;
    for (auto node: nodes)
        node->internalForce += node->velocity * material->d_velocity * new_length / nodes.size();
}

Vector3 Truss::getStress()
{
    return stress * tangent;
};

LineJoint::LineJoint(const std::vector<NodePtr> nodes, std::shared_ptr<TrussMaterial> material)
{
    this->material = material;
    this->nodes = nodes;
    assert(nodes.size() > 3);
}

void LineJoint::explicitStep(double h)
{
    Vector3 t = nodes[1]->position - nodes[0]->position;
    t.normalize();
    for (int i=2; i < nodes.size(); i++)
    {
        Vector3 ti = nodes[i]->position - nodes[0]->position;
        Vector3 n_diff = nodes[0]->position + ti.dot(t) * t - nodes[i]->position;
        nodes[i]->internalForce += this->material->elasticity * n_diff;
    }
}

Vector3 LineJoint::getStress()
    {
        return Vector3(0, 0, 0);
    }

}