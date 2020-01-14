#include "element.h"
#include <iostream>

namespace paraFEM {
    
CoordSys::CoordSys(Vector3 n_vec)
{
    check_nan(n_vec, "CoordSys:n");
    this->n = n_vec;
    this->n.normalize();
    this->t1 = Vector3(1, 0, 0);
    if (this->n.dot(this->t1) > 0.9){
        this->t1 = Vector3(0, 1, 0);
    }
    this->t2 = n.cross(t1);
    this->t1 = t2.cross(n);

    this->t1.normalize();
    this->t2.normalize();
    this->mat.row(0) = t1;
    this->mat.row(1) = t2;
    check_nan(this->t1, "CoordSys:t1");
    check_nan(this->t2, "CoordSys:t2");
}

CoordSys::CoordSys(Vector3 n_vec, Vector3 t1_vec)
{
    check_nan(n_vec, "CoordSys:n");
    this->n = n_vec;
    this->n.normalize();

    this->t2 = n.cross(t1_vec);
    this->t1 = t2.cross(n);
    this->t1.normalize();
    this->t2.normalize();
    this->mat.row(0) = t1;
    this->mat.row(1) = t2;
    check_nan(this->t1, "CoordSys:t1");
    check_nan(this->t2, "CoordSys:t2");
}

void CoordSys::update(Vector3 n_vec)
{
    this->update(n_vec, this->t1);
    check_nan(n_vec, "CoordSys:n:update");
    this->n = n_vec;
    this->n.normalize();
    this->t2 = n.cross(t1);
    this->t2.normalize();
    this->t1 = t2.cross(n);
    this->t1.normalize();
    this->mat.row(0) = t1;
    this->mat.row(1) = t2;
    check_nan(this->t1, "CoordSys:t1:update");
    check_nan(this->t2, "CoordSys:t2:update");
}

void CoordSys::update(Vector3 n_vec, Vector3 t1_vec)
{
    check_nan(n_vec, "CoordSys:n:update");
    this->n = n_vec;
    this->n.normalize();
    this->t1 = t1_vec;
    this->t1.normalize();
    this->t2 = n.cross(t1);
    this->t2.normalize();
    this->mat.row(0) = t1;
    this->mat.row(1) = t2;
    check_nan(this->t1, "CoordSys:t1:update");
    check_nan(this->t2, "CoordSys:t2:update");
}

void CoordSys::rotate(std::vector< Vector2 > first, std::vector< Vector2 > second)
{
   Eigen::Matrix2d rot = find_rot_mat(first, second);
   mat = rot.transpose() * mat;   //rot.T?
   t1 = mat.row(0);
   t2 = mat.row(1);
    check_nan(this->t1, "CoordSys:t1:rotate");
    check_nan(this->t2, "CoordSys:t2:rotate");
}

void CoordSys::rotate(Eigen::MatrixX2d first, Eigen::MatrixX2d second)
{
   Eigen::Matrix2d rot = find_rot_mat(first.transpose(), second.transpose());
   mat = rot.transpose() * mat;   //rot.T?
   t1 = mat.row(0);
   t2 = mat.row(1);
    check_nan(this->t1, "CoordSys:t1:rotate");
    check_nan(this->t2, "CoordSys:t2:rotate");
}



Vector2 CoordSys::toLocal(Vector3 vec)
{
    return mat * vec;
}

Vector3 CoordSys::toGlobal(Vector2 vec)
{
    Vector3 result = mat.transpose() * vec;
    if (is_nan(result)) {
        std::cout << "Error" << this->mat << "//" << vec << "->" << result << std::endl;
        std::cout << "N" << this->n << "//t1" << this->t1 << "//t2: " << this->t2 << std::endl;
    }
    check_nan(result, "CoordSys::toGlobal");

    return result;
}


std::vector< int > Element::getNr()
{
    std::vector< int > numbers;
    for (auto node: nodes)
        numbers.push_back(node->nr);
    return numbers;
}


IntegrationPoint::IntegrationPoint(double eta_in, double zeta_in, double weight_in)
{
    eta = eta_in;
    zeta = zeta_in;
    weight = weight_in;
    stress.setZero();
}


Eigen::Matrix2Xd to_matrix(std::vector<Vector2> in)
{
    Eigen::Matrix2Xd out;
    out.resize(2, in.size());
    for (int i = 0; i < in.size(); i++)
    {
        out.col(i) = in[i];
    }
    return out;
}

Eigen::Matrix2d find_rot_mat(Eigen::Matrix2Xd in, Eigen::Matrix2Xd out)
{
    Eigen::Affine2d A;
    A.linear() = Eigen::Matrix2d::Identity(2, 2);
    A.translation() = Eigen::Vector2d::Zero();

    if (in.cols() != out.cols())
        throw "Find3DAffineTransform(): input data mismatch";

    // First find the scale, by finding the ratio of sums of some distances,
    // then bring the datasets to the same scale.
    double dist_in = 0, dist_out = 0;
    for (int col = 0; col < in.cols()-1; col++) {
        dist_in  += (in.col(col+1) - in.col(col)).norm();
        dist_out += (out.col(col+1) - out.col(col)).norm();
    }
    if (dist_in <= 0 || dist_out <= 0)
        return A.linear();
    double scale = dist_out/dist_in;
    out /= scale;

    // Find the centroids then shift to the origin
    Vector2 in_ctr = Vector2::Zero();
    Vector2 out_ctr = Vector2::Zero();
    for (int col = 0; col < in.cols(); col++) {
        in_ctr  += in.col(col);
        out_ctr += out.col(col);
    }
    in_ctr /= in.cols();
    out_ctr /= out.cols();
    for (int col = 0; col < in.cols(); col++) {
        in.col(col)  -= in_ctr;
        out.col(col) -= out_ctr;
    }

    // SVD
    Eigen::MatrixXd Cov = in * out.transpose();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // Find the rotation
    double d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
    if (d > 0)
        d = 1.0;
    else
        d = -1.0;
    Eigen::Matrix2d I = Eigen::Matrix2d::Identity(2, 2);
    I(1, 1) = d;
    Eigen::Matrix2d R = svd.matrixV() * I * svd.matrixU().transpose();

    // The final transform
    A.linear() = scale * R;
    A.translation() = scale*(out_ctr - R*in_ctr);

    return A.linear();
}

Eigen::Matrix2d find_rot_mat(std::vector<Vector2> p, std::vector<Vector2> q)
{
    return find_rot_mat(to_matrix(p), to_matrix(q));
}

}