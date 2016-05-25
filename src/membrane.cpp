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
    material = mat;
    center.setZero();
    for (auto node: nodes)
        center += node->position;
    center /= nodes.size();
    Vector3 n = (nodes[1]->position - nodes[0]->position).cross(
                         nodes[2]->position - nodes[1]->position);
    area = n.norm() / 2;
    coordSys = CoordSys(n / area / 2., nodes[1]->position - nodes[0]->position);
    int row_count = 0;
    for (auto node: nodes)
    {
        pos_mat.row(row_count) = coordSys.toLocal(node->position - center);
        vel_mat.row(row_count) = coordSys.toLocal(node->velocity);
        node->massInfluence += material->rho * area / nodes.size();
        row_count ++;
    }
    stress.setZero();
}

void Membrane3::makeStep(double h)
{
    // 0: properties
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
    Eigen::Matrix<double, 3, 2> new_pos_mat;
    int row_count = 0;
    for (auto node: nodes)
    {
        new_pos_mat.row(row_count) = coordSys.toLocal(node->position - center);
        row_count ++;
    }
    
    //      find rotation of the local coordinates and rotate the coordinate system
    coordSys.rotate(pos_mat, new_pos_mat);

    row_count = 0;
    for (auto node: nodes)
    {
        pos_mat.row(row_count) = coordSys.toLocal(node->position - center);
        vel_mat.row(row_count) = coordSys.toLocal(node->velocity);
        row_count ++;
    }

    // 2: B-matrix
    dN << -1, 1, 0,
          -1, 0, 1;
    B = (dN * pos_mat).inverse() * dN;

    // 3: strain rate  //? local_velocity as matrix ?
    Vector3 strain_rate;
    strain_rate.setZero();
    for (int i = 0; i < nodes.size(); i++)
    {
        strain_rate += Vector3(B(0,i) * vel_mat(i, 0),
                               B(1,i) * vel_mat(i, 1),
                               B(0,i) * vel_mat(i, 1) +
                               B(1,i) * vel_mat(i, 0));
    }

    // 4: stress rate (hook) ->material,...
    Vector3 stress_rate;
    stress_rate = material->C * strain_rate;
    
    // 5: integrate stress_rate (time, area)
    stress += h * stress_rate;
    Vector3 local_force =  stress * area;
    Vector3 structural_damping = material->d_structural * stress_rate * area;
    
    // 6: nodal forces with virtual power
    // T.T * B.T * stress
    Vector2 local_node_force;
    
    for (int i = 0; i < nodes.size(); i++)
    {
        local_node_force.setZero();
        local_node_force.x() += B(0, i) * (local_force(0) + structural_damping(0))
                             +  B(1, i) * (local_force(2) + structural_damping(2));
        local_node_force.y() += B(1, i) * (local_force(1) + structural_damping(1))
                             +  B(0, i) * (local_force(2) + structural_damping(2));
        nodes[i]->internalForce += coordSys.toGlobal(local_node_force);
    }
    for (auto node: nodes)
        node->internalForce -= pressure / nodes.size() * this->coordSys.n;
    
}

Vector3 Membrane3::getStress()
{
    return stress;
};


Membrane4::Membrane4(std::vector<NodePtr> points, 
                     std::shared_ptr<MembraneMaterial> mat,
                     bool reduced_integration)
{
    
    material = mat;
    nodes = points;
    for (auto node: nodes)
        center += node->position;
    center /= nodes.size();
    Vector3 n = (nodes[2]->position - nodes[0]->position).cross(
                 nodes[3]->position - nodes[1]->position);
    area = n.norm() / 2;
    coordSys = CoordSys(n / area / 2, nodes[1]->position - nodes[0]->position);
    int row_count = 0;
    for (auto node: nodes)
    {
        pos_mat.row(row_count) = coordSys.toLocal(node->position - center);
        vel_mat.row(row_count) = coordSys.toLocal(node->velocity);
        node->massInfluence += material->rho * area / nodes.size();
        row_count ++;
    }
    if (reduced_integration)
    {
        integration_points.push_back(IntegrationPoint(0, 0, 1));
        initHG();
    }
    else
    {
        std::array<double, 2> signs {-1, 1};
        double value =  0.57735026918962584;
        for (auto sign_eta: signs)
        {
            for (auto sign_zeta: signs) 
            {
                integration_points.push_back(IntegrationPoint(sign_eta * value, sign_zeta * value, 1./4.));
            }
        }
    }
}

void Membrane4::initHG()
{
    hg_stress.setZero();
    Eigen::Vector4d h;

    h << 1, -1, 1, -1;
    dN <<  1, 1, -1, -1,
           -1, 1, 1, -1;
    dN /= -4;
    B = (dN * pos_mat).inverse() * dN;

    Eigen::Vector4d B0 =  B.row(0);
    Eigen::Vector4d B1 = B.row(1);

    hg_const = 1. / 8. * 0.03 * (material->C(0,0) * area) * (B0.dot(B0) + B1.dot(B1));
    hg_gamma = h - (h.dot(pos_mat.col(0)) * B0 + h.dot(pos_mat.col(0)) * B1);
}
    
void Membrane4::makeStep(double h)
{    // 0: properties
    
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
    Eigen::Matrix<double, 4, 2> new_pos_mat;
    int row_count = 0;
    for (auto node: nodes)
    {
        new_pos_mat.row(row_count) = coordSys.toLocal(node->position - center);
        row_count ++;
    }
    
    //      find rotation of the local coordinates and rotate the coordinate system
    coordSys.rotate(pos_mat, new_pos_mat);
    
    //      set the local positions by appling the transformation matrix on global positions
    row_count = 0;
    for (auto node: nodes)
    {
        pos_mat.row(row_count) = coordSys.toLocal(node->position - center);
        vel_mat.row(row_count) = coordSys.toLocal(node->velocity);
        row_count ++;
    }
    
    // 2: B-matrix
    // 3: strain rate

    Vector3 strain_rate;
    Vector3 stress_rate;
    Vector2 local_node_force;

    for (IntegrationPoint& int_point: integration_points)
    {
        strain_rate.setZero();
        
        dN <<  1 + int_point.zeta, 1 - int_point.zeta, -1 + int_point.zeta, -1 - int_point.zeta, 
               1 + int_point.eta, -1 - int_point.eta, -1 + int_point.eta, 1 - int_point.eta;
        dN /= -4;
        B = (dN * pos_mat).inverse() * dN;

        for (int i = 0; i < nodes.size(); i++)
        {
            strain_rate += Vector3(B(0,i) * vel_mat(i, 0),
                                B(1,i) * vel_mat(i, 1),
                                B(0,i) * vel_mat(i, 1) +
                                B(1,i) * vel_mat(i, 0));
        }
        // 4: stress rate (hook)
        stress_rate = material->C * strain_rate;

        // 5: integrate sigma
        int_point.stress += h * stress_rate;
        Vector3 local_force = int_point.stress * area;
        Vector3 structural_damping = material->d_structural * stress_rate * area;

        // 6: nodal forces due to internal work
        Vector2 local_node_force;
        
        for (int i = 0; i < nodes.size(); i++)
            {
                local_node_force.setZero();
                local_node_force.x() += B(0, i) * (local_force(0) + structural_damping(0))
                                     +  B(1, i) * (local_force(2) + structural_damping(2));
                local_node_force.y() += B(1, i) * (local_force(1) + structural_damping(1))
                                     +  B(0, i) * (local_force(2) + structural_damping(2));
                nodes[i]->internalForce += int_point.weight * coordSys.toGlobal(local_node_force);
            }

    }
    
    for (auto node: nodes)
        node->internalForce -= pressure / nodes.size() * this->coordSys.n;
    // 7: hourglasscontrol
    if (integration_points.size() == 1)
    {
        //      hourglass stressrate + integration
        hg_stress += hg_const * (vel_mat.transpose() * hg_gamma) * h;
        for (int i = 0; i < nodes.size(); i++)
        {
            local_node_force.setZero();
            local_node_force.x() += hg_gamma[i] * hg_stress.x();
            local_node_force.y() += hg_gamma[i] * hg_stress.y();
            nodes[i]->internalForce += coordSys.toGlobal(local_node_force);
        }
    }
}

Vector3 Membrane4::getStress()
{
    Vector3 stress;
    stress.setZero();
    for (auto int_point: integration_points)
        stress += int_point.stress * int_point.weight;
    return stress;
};

}