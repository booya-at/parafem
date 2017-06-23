#ifndef NODE_H
#define NODE_H

#include "base.h"
#include "Eigen/Geometry"
#include <memory>
#include <vector>


namespace paraFEM
{

typedef Eigen::Vector3d Vector3;
typedef Eigen::Vector2d Vector2;
    
struct Node: Base
{
public:
    Node(double x, double y, double z);
    Vector3 previous_position;
    Vector3 position;
    Vector3 velocity;
    Vector3 acceleration;
    
    Vector3 externalForce;
    Vector3 internalForce;

    void add_external_force(Vector3 force);
    
    double massInfluence = 1;

    void solveEquilibrium(double h, double externalFactor=1);
    Vector3 fixed;
    
    int nr;
};

typedef std::shared_ptr<Node> NodePtr;
typedef std::vector<NodePtr> NodeVec;

}
#endif