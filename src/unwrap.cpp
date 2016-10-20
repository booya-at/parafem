#include "unwrap.h"
#include <Eigen/SparseQR>
#include <Eigen/SVD>

namespace paraFEM{
typedef Eigen::Triplet<double> T;
typedef Eigen::SparseMatrix<double> spMat;

LscmRelax::LscmRelax(std::vector<Eigen::Vector3D> vertices, 
                     std::vector<std::array<3, int>> triangles,
                     std::vector<unsigned int> fixed_pins)
{
    this->vertices = vertices;
    this->triangles = triangles;
    this->flat_vertices.resize(this->vertices.size());
    this->fixed_pins = fixed_pins;

    // set the fixed pins of the flat-mesh:
    this->set_fixed_pins();

    
    unsigned int fixed_count = 0;
    for (unsigned int i=0; i < this->vertices.size(); i++)
    {
        if (i == fixed_pins[fixed_count])
            fixed_count ++;
        else
            this->new_vertex_order.push_back(i - fixed_count);
    }

    for (auto fixed_index: fixed_pins)
        this->new_vertex_order.push_back(fixed_index);

    // get the reversed map:
    old_vertex_order.resize(new_vertex_order.size());
    unsigned int j = 0;
    for (auto index: this->new_vertex_order)
    {
        old_vertex_order[index] = j;
        j++;
    }

}

//////////////////////////////////////////////////////////////////////////
/////////////////                 F.E.M                      /////////////
//////////////////////////////////////////////////////////////////////////
LscmRelax::relax(double step_size)
{
        
}


//////////////////////////////////////////////////////////////////////////
/////////////////                 L.S.C.M                    /////////////
//////////////////////////////////////////////////////////////////////////
LscmRelax::lscm()
{
    this->set_q_l_0();
    std::vector<T> triple_list;
    unsigned int i = 0;
    double x21, x31, y31. x32;

    // 1. create the triplet list (t * 2, v * 2)
    for(auto triange: this->triangles)
    {
        x21 = this->q_l_0(i, 0);
        x31 = this->q_l_0(i, 1):
        y31 = this->q_l_0(i, 2);
        x32 = x31 - x21;

        triple_list.push_back(T(2 * i, get->col(tri[0]) * 2, x32));
        triple_list.push_back(T(2 * i, get->col(tri[0]) * 2 + 1, -y31));
        triple_list.push_back(T(2 * i, get->col(tri[1]) * 2, -x31));
        triple_list.push_back(T(2 * i, get->col(tri[1]) * 2 + 1, y31));
        triple_list.push_back(T(2 * i, get->col(tri[2]) * 2, x21));

        triple_list.push_back(T(2 * i + 1, get->col(tri[0]) * 2, y31));
        triple_list.push_back(T(2 * i + 1, get->col(tri[0]) * 2 + 1, x32));
        triple_list.push_back(T(2 * i + 1, get->col(tri[1]) * 2, -y31));
        triple_list.push_back(T(2 * i + 1, get->col(tri[1]) * 2 + 1, -x31));
        triple_list.push_back(T(2 * i + 1, get->col(tri[2]) * 2 + 1, x21));

        i++;
    }

    // 2. divide the triplets in matrix(unknown part) and rhs(known part) and reset the position 
    std::vector<T> rhs_triplets;
    std::vector<T> mat_triplets;
    for (auto triplet: triple_list)
    {
        if (triplet.col() > (this->vertices.size() - this->fixed_pins.size()) * 2)
            rhs_triplets.push_back(triplet)
        else
            mat_triplets.push_back(triplet);
    }

    // 3. create a rhs_pos vector
    Eigen::VectorXd rhs_pos(this->vertices.size() * 2);
    rhs_pos.setZero();
    for (auto index: this->fixed_pins)
    {
        rhs_pos[index * 2] = this->flat_vertices[index].x();      //TODO: not yet set
        rhs_pos[index * 2 + 1] = this->flat_vertices[index].y();
    }

    // 4. fill a sparse matrix and calculdate the rhs
    VectorXd rhs(this->triangles.size() * 2); // maybe use a sparse vector
    spMat B(this->triangles.size() * 2, this->vertices.size() * 2);
    B.setFromTriplets(rhs_triplets);
    rhs = B * rhs_pos;

    // 5. create the lhs matrix
    spMat A(this->triangles.size() * 2, (this->vertices.size() - this->fixed_pins.size()) * 2);
    A.setFromTriplets(mat_triplets);

    // 6. solve the system and set the flatted coordinates
    Eigen::SparseQR<spMat> solver;
    VectorXd sol(this->vertices.size() * 2);
    sol = solver.compute(A).solve(rhs);

        // TODO: create function, is needed also in the fem step
    this->flat_vertices.resize(this->vertices.size());
    for (unsigned int i=0; i < this->vertices.size(); i++)
        this->flat_vertices[this->get_old_order(i)] = Vector2(sol[i * 2], sol[i * 2 + 1]);
    // 7. if size of fixed pins <= 2: scale the map to the same area as the 3d mesh

}

Eigen::MatrixXd LscmRelax::set_q_l_g()
{
    // get the coordinates of a triangle in local coordinates from the 3d mesh
    // x1, y1, y2 = 0
    // -> vector<x2, x3, y3>
    this->q_l_0.resize(this->triangle.size(), 3);
    for (auto triangle: this->triangles)
    {
        Vector3 r1 = this->vertices[tri[0]];
        Vector3 r2 = this->vertices[tri[1]];
        Vector3 r3 = this->vertices[tri[2]];
        Vector3 r21 = r2 - r1;
        Vector3 r31 = r3 - r1;
        this->q_l_0 << r21.norm();
        r21.normalize();
        this->q_l_0 << r31.dot(r21), r31.cross(r21).norm();
    }
}

Eigen::MatrixXd LscmRelax::set_q_l_m()
{
    // get the coordinates of a triangle in local coordinates from the 2d map
    // x1, y1, y2 = 0
    // -> vector<x2, x3, y3>
    this->q_l_1.resize(this->triangle.size(), 3);
    this->alpha.resize(this->triangle.size(), 2);
    for (auto triangle: this->triangles)
    {
        Vector2 r1 = this->flat_vertices[tri[0]];
        Vector2 r2 = this->flat_vertices[tri[1]];
        Vector2 r3 = this->flat_vertices[tri[2]];
        Vector2 r21 = r2 - r1;
        Vector2 r31 = r3 - r1;
        this->q_l_1 << r21.norm();
        r21.normalize();
        this->q_l_1 << r31.dot(r21), r31.x() * r21.y() - r31.y() * r21.x();
        this->alpha << r21.x(), r21.y();
    }
}

std::array<unsigned int> LscmRelax::set_fixed_pins()
{
    // TODO!!!!
    // if less then one fixed pin is set find two by an automated algorithm and align them to y = 0
    // if more then two pins are choosen find a leastsquare-plane and project the points on it
    // insert the points in the flat-vertices vector
    // TODO!!!!
    if (this->fixed_pins.size() == 0)
        this->fixed_pins.push_back(0);
    if (this->fixed_pins.size() == 1)
    {
        this->fixed_pins.push_back(get_max_distance(this->vertices[0], this->vertices));
        double dist = (this->vertices[this->fixed_pins[0]] - this->vertices[this->fixed_pins[1]]).norm();
        this->flat_vertices[this->fixed_pins[0]] = Vector2(0, 0);
        this->flat_vertices[this->fixed_pins[1]] = Vector2(dist, 0);
    }
    if (this->fixed_pins > 2)
    {
        std::vector<Vector3> fixed_3d_points;
        for (unsigned int index: this->fixed_pins)
            fixed_3d_points.push_back(this->vertices[index]);
        std::vector<Vector2> maped_points = map_to_2D(fixed_3d_points);
        unsigned int i = 0;
        for (unsigned int index: this->fixed_pins)
        {
            this->flat_vertices[i] = maped_points[i];
            i++;
        }
    }
}



unsigned int LscmRelax::get_new_order(unsigned int i)
{
    return this->new_vertex_order[i];
}

unsigned int LscmRelax::get_old_order(unsigned int i)
{
    return this->old_vertex_order[i];
}

std::vector<Vector2> map_to_2D(std::vector<Vector3> points)
{
    // debugging!!!
    Eigen::MatrixXd mat(points.size(), 4, 1);
    for(auto point: points)
    {
        mat << point;
        mat << 1;
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Vector3 n = svd.matrixV().row(2)
    n.normalize();
    Vector3 y(0, 1, 0);
    if ((n - y).norm() < 0.0001)
        y = Vector3(1, 0, 0);
    Vector3 x = n.cross(y);
    x.normalize();
    y = np.cross(n, x)
    Matrix<2, 3> transform;
    transform.col(0) = x;
    transform.col(1) = y;
    return mat_to_vector<Vector2>(mat_from_vector<Vector3>(points) * transform);
}


unsigned int get_max_distance(Vector3 point, std::vector<Vector3> vertices)
{
    double max_dist = 0;
    double dist;
    unsigned int max_dist_index;
    unsigned int j = 0;
    for (auto point_j: vertices)
    {
        dist = (point - point_j).norm();
        if (dist > max_dist)
        {
            max_dist = dist;
            max_dist_index = j;
        }
        j ++;
    }
    return max_dist_index;
}

}
