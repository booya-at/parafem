#include "element.h"
#include <iostream>
#include "utils.h"

namespace paraFEM {

MaterialPtr Membrane::getMaterial()
{
    return std::dynamic_pointer_cast<Material>(this->material);
}

Vector3 Membrane::calculate_center() {
    this->center.setZero();

    //      find new center
    for (auto node: this->nodes) {
        this->center += node->position;
    }
    this->center /= this->nodes.size();

    return this->center;

}
    
    
Membrane3::Membrane3(std::vector<NodePtr> points, std::shared_ptr<MembraneMaterial> mat)
{
    this->nodes = points;
    this->material = mat;
    this->pressure = 0;
    this->calculate_center();

    Vector3 n = (nodes[1]->position - nodes[0]->position).cross(
                         nodes[2]->position - nodes[1]->position);
    this->area = n.norm() / 2;
    if (this->area == 0 || std::isnan(this->area))
    {

        std::cout << "warning area " << this->area << std::endl;
        for (auto node: nodes)
        {
            std::cout << node->position << std::endl;
        }
        throw FemException("Zero-Area membrane");
    }
    this->coordSys = CoordSys(n / area / 2., nodes[1]->position - nodes[0]->position);
    int row_count = 0;
    for (auto node: nodes)
    {
        pos_mat.row(row_count) = coordSys.toLocal(node->position - center);
        vel_mat.row(row_count) = coordSys.toLocal(node->velocity);
        node->massInfluence += material->rho * area / nodes.size();
        row_count ++;
    }
    stress.setZero();
    double smallestSide = (nodes.back()->position - nodes[0]->position).norm();
    for (int i=0; i < nodes.size() - 1; i++)
    {
        double sideLen = (nodes[i]->position - nodes[i+1]->position).norm();
        if (sideLen < smallestSide)
            smallestSide = sideLen;
    }
    dViscous = smallestSide * material->waveSpeed() * 
               material->rho * material->d_minMode;
    
    characteristicLength = smallestSide;
}

void Membrane3::geometryStep() {
    this->calculate_center();


    // 1: new coordSys
    //      compute the new normal of the element and update the coordinate system
    Vector3 n = (this->nodes[1]->position - this->nodes[0]->position).cross(
                 this->nodes[2]->position - this->nodes[1]->position);
    this->area = n.norm() / 2.;
    check_nan(n, "Membrane3::explicitStep:n");
    this->coordSys.update(n / area / 2, this->nodes[1]->position - this->nodes[0]->position);

    //      get the local position for the updated coordinate system
    Eigen::Matrix<double, 3, 2> new_pos_mat;
    int row_count = 0;
    for (auto node: this->nodes)
    {
        new_pos_mat.row(row_count) = this->coordSys.toLocal(node->position - this->center);
        row_count ++;
    }

    //      find rotation of the local coordinates and rotate the coordinate system
    // coordSys.rotate(pos_mat, new_pos_mat);

    row_count = 0;
    for (auto node: nodes)
    {
        this->pos_mat.row(row_count) = coordSys.toLocal(node->position - center);
        this->vel_mat.row(row_count) = coordSys.toLocal(node->velocity);
        row_count ++;
    }



    // 2: B-matrix
    this->dN << -1, 1, 0,
                -1, 0, 1;
    this->B = (this->dN * this->pos_mat).inverse() * this->dN;

}

void Membrane3::explicitStep(double h)
{
    double step_size = h;
    if (not this->is_valid)
        return;


    // 3: strain rate  //? local_velocity as matrix ?
    Vector3 strain_rate;
    strain_rate.setZero();
    for (int i = 0; i < nodes.size(); i++)
    {
        strain_rate += Vector3(this->B(0,i) * this->vel_mat(i, 0),
                               this->B(1,i) * this->vel_mat(i, 1),
                               this->B(0,i) * this->vel_mat(i, 1) +
                               this->B(1,i) * this->vel_mat(i, 0));
    }

    // 4: stress rate (hook) ->material,...
    Vector3 stress_rate;
    stress_rate = this->material->C * strain_rate;
    
    // 5: integrate stress_rate (time, area)
    stress += step_size * stress_rate;
    Vector3 local_force =  stress * this->area;
    Vector3 structural_damping = this->material->d_structural * stress_rate * this->area;

    Vector3 hydrostatic_damping = Vector3(1, 1, 0);
    hydrostatic_damping *= this->dViscous * (strain_rate.x() + strain_rate.y());
    
    // 6: nodal forces with virtual power
    // T.T * B.T * stress
    Vector2 local_node_force;
    
    for (int i = 0; i < this->nodes.size(); i++)
    {
        local_node_force.setZero();
        local_node_force.x() += this->B(0, i) * (local_force(0) + structural_damping(0) + hydrostatic_damping(0))
                             +  this->B(1, i) * (local_force(2) + structural_damping(2));
        local_node_force.y() += this->B(1, i) * (local_force(1) + structural_damping(1) + hydrostatic_damping(1))
                             +  this->B(0, i) * (local_force(2) + structural_damping(2));
        this->nodes[i]->internalForce += this->coordSys.toGlobal(local_node_force);
    }
    Vector3 force_cp = this->pressure / nodes.size() * this->coordSys.n;
    double damping = this->material->d_velocity * this->area / this->nodes.size();
    for (auto node: this->nodes)
    {
        node->internalForce -= force_cp;
        node->internalForce += node->velocity * damping;

        check_nan(node->internalForce, "membrane3:explicitstep:node_internal_force");
    }

    
}

Vector3 Membrane3::getStress()
{
    return stress;
};


Membrane4::Membrane4(std::vector<NodePtr> points, 
                     std::shared_ptr<MembraneMaterial> mat,
                     bool reduced_integration)
{
    
    this->material = mat;
    this->nodes = points;
    this->pressure = 0;
    this->calculate_center();
    Vector3 n = (nodes[2]->position - nodes[0]->position).cross(
                 nodes[3]->position - nodes[1]->position);
    area = n.norm() / 2;
    if (area == 0)
    {
        is_valid = false;
        return;
    }
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
    double smallestSide = (nodes.back()->position - nodes[0]->position).norm();
    for (int i=0; i < nodes.size() - 1; i++)
    {
        double sideLen = (nodes[i]->position - nodes[i+1]->position).norm();
        if (sideLen < smallestSide)
            smallestSide = sideLen;
    }
    dViscous = smallestSide * material->waveSpeed() * 
               material->rho * material->d_minMode;
    characteristicLength = smallestSide;

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

void Membrane4::geometryStep() {
    this->calculate_center();

    // 1: new coordSys
    //      compute the new normal of the element and update the coordinate system
    Vector3 n = (nodes[2]->position - nodes[0]->position).cross(
                 nodes[3]->position - nodes[1]->position);
    area = n.norm() / 2.;
    check_nan(n, "Membrane3::explicitStep:n");
    this->coordSys.update(n / area / 2.);

    //      get the local position for the updated coordinate system
    Eigen::Matrix<double, 4, 2> new_pos_mat;
    int row_count = 0;
    for (auto node: this->nodes)
    {
        new_pos_mat.row(row_count) = this->coordSys.toLocal(node->position - this->center);
        row_count ++;
    }

    //      find rotation of the local coordinates and rotate the coordinate system
    //coordSys.rotate(pos_mat, new_pos_mat); TODO: CHECK!!!

    //      set the local positions by appling the transformation matrix on global positions
    row_count = 0;
    for (auto node: this->nodes)
    {
        this->pos_mat.row(row_count) = this->coordSys.toLocal(node->position - this->center);
        this->vel_mat.row(row_count) = this->coordSys.toLocal(node->velocity);
        row_count ++;
    }

}
    
void Membrane4::explicitStep(double h)
{
    double step_size = h;
    // 0: properties
    if (not is_valid)
        return;

    // 2: B-matrix
    // 3: strain rate

    Vector3 strain_rate;
    Vector3 stress_rate;
    Vector2 local_node_force;
    for (auto node: this->nodes) {
        check_nan(node->internalForce, "membrane4:explicitstep:node_internal_force0");
    }

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
        int_point.stress += step_size * stress_rate;
        Vector3 local_force = int_point.stress * area;
        Vector3 structural_damping = material->d_structural * strain_rate * area;
        Vector3 hydrostatic_damping = Vector3(1, 1, 0);
        hydrostatic_damping *= dViscous * (strain_rate.x() + strain_rate.y());
        // 6: nodal forces due to internal work
        Vector2 local_node_force;
        
        for (int i = 0; i < nodes.size(); i++)
            {
                local_node_force.setZero();
                local_node_force.x() += B(0, i) * (local_force(0) + structural_damping(0) + hydrostatic_damping(0))
                                     +  B(1, i) * (local_force(2) + structural_damping(2));
                local_node_force.y() += B(1, i) * (local_force(1) + structural_damping(1) + hydrostatic_damping(1))
                                     +  B(0, i) * (local_force(2) + structural_damping(2));
                nodes[i]->internalForce += int_point.weight * coordSys.toGlobal(local_node_force);
            }

    }

    for (auto node: this->nodes) {
        check_nan(node->internalForce, "membrane4:explicitstep:node_internal_force1");
    }
    for (auto node: nodes)
    {
        node->internalForce -= pressure / nodes.size() * this->coordSys.n;
        node->internalForce += node->velocity * material->d_velocity * area / nodes.size();
    }

    for (auto node: this->nodes) {
        check_nan(node->internalForce, "membrane4:explicitstep:node_internal_force2");
    }
    // 7: hourglasscontrol
    if (integration_points.size() == 1)
    {
        //      hourglass stressrate + integration
        hg_stress += hg_const * (vel_mat.transpose() * hg_gamma) * step_size;
        for (int i = 0; i < nodes.size(); i++)
        {
            local_node_force.setZero();
            local_node_force.x() += hg_gamma[i] * hg_stress.x();
            local_node_force.y() += hg_gamma[i] * hg_stress.y();
            nodes[i]->internalForce += coordSys.toGlobal(local_node_force);
        }
    }

    for (auto node: this->nodes) {
        check_nan(node->internalForce, "membrane4:explicitstep:node_internal_force3");
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