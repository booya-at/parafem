#include "node.h"
#include <iostream>
#include "utils.h"

namespace parafem {

Node::Node(double x, double y, double z)
{
    this->position = Eigen::Vector3d(x, y, z);
    this->velocity.setZero();
    this->acceleration.setZero();
    this->mass_influence = 1.;
    this->fixed.setOnes();
    this->external_force.setZero();
    this->internal_force.setZero();
}

void Node::solveEquilibrium(double h, double external_factor)
{
    this->acceleration = 1. / this->mass_influence * (this->external_force * external_factor - this->internal_force);
    this->velocity += h * (this->acceleration.cwiseProduct(fixed));   // at time t + 0.5h
    this->position += this->velocity * h;       // at time t + 1

    /*
    if (is_nan(this->acceleration) || is_nan(this->velocity) || is_nan(this->position)){
        std::cout << "ERR" << this->mass_influence << "/" << this->external_force << "/" << this->internal_force << std::endl;
    }*/
    check_nan(this->acceleration, "node_solve_equilibrium:acceleration");
    check_nan(this->velocity, "node_solve_equilibrium:velocity");
    check_nan(this->position, "node_solve_equilibrium:position");

}

void Node::add_external_force(Vector3 force)
{
    external_force += force;
}

}