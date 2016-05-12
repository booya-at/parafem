#include "element.h"
#include <iostream>


namespace paraFEM {

void Membrane::setConstPressure(double p)
{
    for (auto node: nodes)
        node->externalForce += coordSys.n * area * p / nodes.size();
}
    
    
    
Membrane3::Membrane3(std::vector<NodePtr> points, std::shared_ptr<MembraneMaterial> mat)
{
    nodes = points;
    center.setZero();
    for (auto node: nodes)
        center += node->position;
    center /= nodes.size();
    Vector3 n = (nodes[1]->position - nodes[0]->position).cross(
                         nodes[2]->position - nodes[1]->position);
    area = n.norm() / 2;
    coordSys = CoordSys(n / area / 2.);
    for (auto node: nodes)
    {
        local_position.push_back(coordSys.toLocal(node->position - center));
        node->massInfluence += area / 3. * mat->rho;
    }
    material = mat;
    stress.setZero();
}

void Membrane3::makeStep(double h)
{
    // 0: properties
    
    std::vector<Vector2> new_local_position;
    center.setZero();
    
    //      find new center
    for (auto node: nodes)
        center += node->position;
    center /= nodes.size();
    
    // 1: new coordSys
    //      compute the new normal of the element and update the coordinate system
    Vector3 n = (nodes[1]->position - nodes[0]->position).cross(
                 nodes[2]->position - nodes[1]->position);
    area = n.norm() / 2.;
    coordSys.update(n / area / 2);

    //      get the local position for the updated coordinate system
    for (auto node: nodes)
        new_local_position.push_back(coordSys.toLocal(node->position - center));
    
    //      find rotation of the local coordinates and rotate the coordinate system
    coordSys.rotate(local_position, new_local_position);
    
    //      set the local positions by appling the transformation matrix on global positions
    local_position.resize(0);
    local_velocity.resize(0);
    for (auto node: nodes)
    {
        local_position.push_back(coordSys.toLocal(node->position - center));
        local_velocity.push_back(coordSys.toLocal(node->velocity));
    }
                     
    // 2: B-matrix
    Eigen::Matrix<double, 2, 3> B;
    B << local_position[1].y() - local_position[2].y(),
         local_position[2].y() - local_position[0].y(),
         local_position[0].y() - local_position[1].y(),
         local_position[2].x() - local_position[1].x(),
         local_position[0].x() - local_position[2].x(),
         local_position[1].x() - local_position[0].x();
    B /= area * 2;
    
    // 3: strain rate  //? local_velocity as matrix ?
    Vector3 strain_rate;
    strain_rate.setZero();
    for (int i = 0; i < nodes.size(); i++)
    {
        strain_rate += Vector3(B(0,i) * local_velocity[i].x(),
                               B(1,i) * local_velocity[i].y(),
                               B(0,i) * local_velocity[i].y() +
                               B(1,i) * local_velocity[i].x());
    }

    // 4: stress rate (hook) ->material,...
    Vector3 stress_rate;
    stress_rate = material->C * strain_rate;
    
    // 5: integrate stress_rate (time, area)
    stress += h * stress_rate;
    Vector3 local_force =  stress * area;
    
    // 6: nodal forces with virtual power
    // T.T * B.T * stress
    Vector2 local_node_force;
    for (int i = 0; i < nodes.size(); i++)
    {
        local_node_force.setZero();
        local_node_force.x() += B(0, i) * (local_force(0) + material->d_structural * strain_rate(0))
                             +  B(1, i) * (local_force(2) + material->d_structural * strain_rate(2));
        local_node_force.y() += B(1, i) * (local_force(1) + material->d_structural * strain_rate(1))
                             +  B(0, i) * (local_force(2) + material->d_structural * strain_rate(2));
        nodes[i]->internalForce += coordSys.toGlobal(local_node_force);
    }
}


Membrane4::Membrane4(std::vector<NodePtr> points, std::shared_ptr<MembraneMaterial> mat)
{
    nodes = points;
    for (auto node: nodes)
        center += node->position;
    center /= nodes.size();
    Vector3 n = (nodes[2]->position - nodes[0]->position).cross(
                 nodes[3]->position - nodes[1]->position);
    area = n.norm() / 2;
    coordSys = CoordSys(n / area / 2);
    for (auto node: nodes)
    {
        local_position.push_back(coordSys.toLocal(node->position - center));
        local_velocity.push_back(coordSys.toLocal(node->velocity));
        node->massInfluence += area / 4. * mat->rho;
    }
    material = mat;
    stress.setZero();
    
    initHG();
}

void Membrane4::initHG()
{
    hg_stress.setZero();
    Eigen::Vector4d h;
    h << 1, -1, 1, -1;
    Eigen::Vector4d x;
    x << nodes[0]->position.x(),
         nodes[1]->position.x(),
         nodes[2]->position.x(),
         nodes[3]->position.x();
    Eigen::Vector4d y;
    y << nodes[0]->position.y(),
         nodes[1]->position.y(),
         nodes[2]->position.y(),
         nodes[3]->position.y();
    Eigen::Vector4d B0;
    B0 << local_position[1].y() - local_position[3].y(),
          local_position[2].y() - local_position[0].y(),
          local_position[3].y() - local_position[1].y(),
          local_position[0].y() - local_position[2].y();
    Eigen::Vector4d B1;
    B1 << local_position[3].x() - local_position[1].x(),
          local_position[0].x() - local_position[2].x(),
          local_position[1].x() - local_position[3].x(),
          local_position[2].x() - local_position[0].x();
    B0 /= 2.*area;
    B1 /= 2.*area;
    hg_const = 1. / 8. * 0.03 * (material->C(0,0) * area) * (B0.dot(B0) + B1.dot(B1)) * 1000.;
    hg_gamma = (h - ((h.transpose() * x)(0,0) * B0 + (h.transpose() * y)(0,0) * B1));
}
    
void Membrane4::makeStep(double h)
{    // 0: properties
    
    std::vector<Vector2> new_local_position;
    center.setZero();
    
    //      find new center
    for (auto node: nodes)
        center += node->position;
    center /= nodes.size();
    
    // 1: new coordSys
    //      compute the new normal of the element and update the coordinate system
    Vector3 n = (nodes[2]->position - nodes[0]->position).cross(
                 nodes[3]->position - nodes[1]->position);
    area = n.norm() / 2.;
    coordSys.update(n / area / 2.);

    //      get the local position for the updated coordinate system
    for (auto node: nodes)
        new_local_position.push_back(coordSys.toLocal(node->position - center));
    
    //      find rotation of the local coordinates and rotate the coordinate system
    coordSys.rotate(local_position, new_local_position);
    
    //      set the local positions by appling the transformation matrix on global positions
    local_position.resize(0);
    local_velocity.resize(0);
    for (auto node: nodes)
    {
        local_position.push_back(coordSys.toLocal(node->position - center));
        local_velocity.push_back(coordSys.toLocal(node->velocity));
    }
    
    // 2: B-matrix
    Eigen::Matrix<double, 2, 4> B;
    B << local_position[1].y() - local_position[3].y(),
         local_position[2].y() - local_position[0].y(),
         local_position[3].y() - local_position[1].y(),
         local_position[0].y() - local_position[2].y(),
         local_position[3].x() - local_position[1].x(),
         local_position[0].x() - local_position[2].x(),
         local_position[1].x() - local_position[3].x(),
         local_position[2].x() - local_position[0].x();
    B /= area * 2;

    // 3: strain rate
    Vector3 strain_rate;
    strain_rate.setZero();
    for (int i = 0; i < nodes.size(); i++)
    {
        strain_rate += Vector3(B(0,i) * local_velocity[i].x(),
                               B(1,i) * local_velocity[i].y(),
                               B(0,i) * local_velocity[i].y() +
                               B(1,i) * local_velocity[i].x());
    }

    // 4: stress rate (hook)
    Vector3 stress_rate;
    stress_rate = material->C * strain_rate;

    // 5: integrate sigma
    stress += h * stress_rate;
    Vector3 local_force = stress * area;

    // 6: nodal forces due to internal work
    Vector2 local_node_force;
    for (int i = 0; i < nodes.size(); i++)
    {
        local_node_force.setZero();
        local_node_force.x() += B(0, i) * (local_force(0) + material->d_structural * strain_rate(0));
                             +  B(1, i) * (local_force(2) + material->d_structural * strain_rate(2));
        local_node_force.y() += B(1, i) * (local_force(1) + material->d_structural * strain_rate(1));
                             +  B(0, i) * (local_force(2) + material->d_structural * strain_rate(2));
        nodes[i]->internalForce += coordSys.toGlobal(local_node_force);
    }
    // 7: hourglasscontrol
    Eigen::Matrix<double, 2, 4> local_vel_mat;
    local_vel_mat << local_velocity[0].x(),
                     local_velocity[1].x(),
                     local_velocity[2].x(),
                     local_velocity[3].x(),
                     local_velocity[0].y(),
                     local_velocity[1].y(),
                     local_velocity[2].y(),
                     local_velocity[3].y();

    //      hourglass stressrate + integration
    hg_stress += hg_const * (local_vel_mat * hg_gamma) * h;
    for (int i = 0; i < nodes.size(); i++)
    {
        local_node_force.setZero();
        local_node_force.x() += hg_gamma[i] * hg_stress.x();
        local_node_force.y() += hg_gamma[i] * hg_stress.y();
        nodes[i]->internalForce += coordSys.toGlobal(local_node_force);
    }
}


}