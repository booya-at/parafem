#include "node.h"
#include <iostream>
#include "utils.h"

namespace paraFEM {

Node::Node(double x, double y, double z)
{
    this->position = Eigen::Vector3d(x, y, z);
    this->velocity.setZero();
    this->acceleration.setZero();
    this->massInfluence = 1.;
    this->fixed.setOnes();
    this->externalForce.setZero();
    this->internalForce.setZero();
}

void Node::solveEquilibrium(double h, double externalFactor)
{
    this->acceleration = 1. / this->massInfluence * (this->externalForce * externalFactor - this->internalForce);
    this->velocity += h * (this->acceleration.cwiseProduct(fixed));   // at time t + 0.5h
    this->position += this->velocity * h;       // at time t + 1

    /*
    if (is_nan(this->acceleration) || is_nan(this->velocity) || is_nan(this->position)){
        std::cout << "ERR" << this->massInfluence << "/" << this->externalForce << "/" << this->internalForce << std::endl;
    }*/
    check_nan(this->acceleration, "node_solve_equilibrium:acceleration");
    check_nan(this->velocity, "node_solve_equilibrium:velocity");
    check_nan(this->position, "node_solve_equilibrium:position");

}

void Node::add_external_force(Vector3 force)
{
    externalForce += force;
}

}