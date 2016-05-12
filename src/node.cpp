#include "node.h"


namespace paraFEM {

Node::Node(double x, double y, double z)
{
    position = Eigen::Vector3d(x, y, z);
    velocity.setZero();
    acceleration.setZero();
    massInfluence = 0;
    fixed.setOnes();
    externalForce.setZero();
    internalForce.setZero();
}

void Node::solveEquilibrium(double h)
{
    acceleration = 1 / massInfluence * (externalForce - internalForce);
    velocity += h * (acceleration.cwiseProduct(fixed));   // at time t + 0.5h
    position += velocity * h;       // at time t + 1
}


}