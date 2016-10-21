#include "unwrap.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SVD>
#include <iostream>
#include <algorithm>

// TODO:
// vertices, flat vertices = EIgen::MatrixXd (better python conversation, parallell)


namespace paraFEM{

typedef Eigen::Triplet<double> T;
typedef Eigen::SparseMatrix<double> spMat;



std::vector<Vector2> map_to_2D(std::vector<Vector3> points)
{
    // debugging!!!
    Eigen::MatrixXd mat(points.size(), 4);
    for(auto point: points)
    {
        mat << point.x(), point.y(), point.z(), 1;
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Vector3 n = svd.matrixV().row(2);
    n.normalize();
    Vector3 y(0, 1, 0);
    if ((n - y).norm() < 0.0001)
        y = Vector3(1, 0, 0);
    Vector3 x = n.cross(y);
    x.normalize();
    y = n.cross(x);
    Eigen::Matrix<double, 3, 2> transform;
    transform.col(0) = x;
    transform.col(1) = y;
    return mat2vec<Vector2>(vec2mat<Vector3>(points) * transform);
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

void LscmRelax::init(std::vector<Vector3> vertices, 
                     std::vector<std::array<int, 3>> triangles,
                     std::vector<int> fixed_pins)
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
        if (fixed_count < this->fixed_pins.size())
        {
            if (i == this->fixed_pins[fixed_count])
                fixed_count ++;
            else
                this->old_vertex_order.push_back(i);
        }
        else
            this->old_vertex_order.push_back(i);
    }

    for (auto fixed_index: this->fixed_pins)
        this->old_vertex_order.push_back(fixed_index);

    // get the reversed map:
    this->new_vertex_order.resize(this->old_vertex_order.size());
    unsigned int j = 0;
    for (auto index: this->old_vertex_order)
    {
        this->new_vertex_order[index] = j;
        j++;
    }
}

LscmRelax::LscmRelax(std::vector<Vector3> vertices, 
                     std::vector<std::array<int, 3>> triangles,
                     std::vector<int> fixed_pins)
{
    this->init(vertices, triangles, fixed_pins);
}

LscmRelax::LscmRelax(std::vector<std::array<double, 3>> vertices, 
                     std::vector<std::array<int, 3>> triangles,
                     std::vector<int> fixed_pins)
{
    std::vector<Vector3> v;
    for (auto arr: vertices)
        v.push_back(Vector3(arr[0], arr[1], arr[2]));
    this->init(v, triangles, fixed_pins);
}

//////////////////////////////////////////////////////////////////////////
/////////////////                 F.E.M                      /////////////
//////////////////////////////////////////////////////////////////////////
void LscmRelax::relax(double step_size)
{
        
}


//////////////////////////////////////////////////////////////////////////
/////////////////                 L.S.C.M                    /////////////
//////////////////////////////////////////////////////////////////////////
void LscmRelax::lscm()
{
    this->set_q_l_g();
    std::vector<T> triple_list;
    unsigned int i = 0;
    double x21, x31, y31, x32;

    // 1. create the triplet list (t * 2, v * 2)
    for(auto triangle: this->triangles)
    {
        x21 = this->q_l_g(i, 0);
        x31 = this->q_l_g(i, 1);
        y31 = this->q_l_g(i, 2);
        x32 = x31 - x21;

        triple_list.push_back(T(2 * i, this->get_new_order(triangle[0]) * 2, x32));
        triple_list.push_back(T(2 * i, this->get_new_order(triangle[0]) * 2 + 1, -y31));
        triple_list.push_back(T(2 * i, this->get_new_order(triangle[1]) * 2, -x31));
        triple_list.push_back(T(2 * i, this->get_new_order(triangle[1]) * 2 + 1, y31));
        triple_list.push_back(T(2 * i, this->get_new_order(triangle[2]) * 2, x21));

        triple_list.push_back(T(2 * i + 1, this->get_new_order(triangle[0]) * 2, y31));
        triple_list.push_back(T(2 * i + 1, this->get_new_order(triangle[0]) * 2 + 1, x32));
        triple_list.push_back(T(2 * i + 1, this->get_new_order(triangle[1]) * 2, -y31));
        triple_list.push_back(T(2 * i + 1, this->get_new_order(triangle[1]) * 2 + 1, -x31));
        triple_list.push_back(T(2 * i + 1, this->get_new_order(triangle[2]) * 2 + 1, x21));

        i++;
    }
    // 2. divide the triplets in matrix(unknown part) and rhs(known part) and reset the position
    std::vector<T> rhs_triplets;
    std::vector<T> mat_triplets;
    for (auto triplet: triple_list)
    {
        if (triplet.col() > (this->vertices.size() - this->fixed_pins.size()) * 2 - 1)
            rhs_triplets.push_back(triplet);
        else
            mat_triplets.push_back(triplet);
    }

    // 3. create a rhs_pos vector
    Eigen::VectorXd rhs_pos(this->vertices.size() * 2);
    rhs_pos.setZero();
    for (auto index: this->fixed_pins)
    {
        rhs_pos[this->get_new_order(index) * 2] = this->flat_vertices[index].x();      //TODO: not yet set
        rhs_pos[this->get_new_order(index) * 2 + 1] = this->flat_vertices[index].y();
    }

    // 4. fill a sparse matrix and calculdate the rhs
    Eigen::VectorXd rhs(this->triangles.size() * 2); // maybe use a sparse vector
    spMat B(this->triangles.size() * 2, this->vertices.size() * 2);
    B.setFromTriplets(rhs_triplets.begin(), rhs_triplets.end());
    rhs = B * rhs_pos;

    // 5. create the lhs matrix
    spMat A(this->triangles.size() * 2, (this->vertices.size() - this->fixed_pins.size()) * 2);
    A.setFromTriplets(mat_triplets.begin(), mat_triplets.end());

    // 6. solve the system and set the flatted coordinates
    // Eigen::SparseQR<spMat, Eigen::COLAMDOrdering<int> > solver;
    Eigen::LeastSquaresConjugateGradient<spMat > solver;
    Eigen::VectorXd sol(this->vertices.size() * 2);
    solver.compute(A);
    sol = solver.solve(-rhs);

        // TODO: create function, is needed also in the fem step
    this->flat_vertices.resize(this->vertices.size());
    for (unsigned int i=0; i < this->vertices.size(); i++)
    {
        if (sol.size() > i * 2 + 1)
            this->flat_vertices[this->get_old_order(i)] = Vector2(sol[i * 2], sol[i * 2 + 1]);
    }
    // 7. if size of fixed pins <= 2: scale the map to the same area as the 3d mesh

}

void LscmRelax::set_q_l_g()
{
    // get the coordinates of a triangle in local coordinates from the 3d mesh
    // x1, y1, y2 = 0
    // -> vector<x2, x3, y3>
    this->q_l_g.resize(this->triangles.size(), 3);
    unsigned int row_count = 0;
    for (auto triangle: this->triangles)
    {
        Vector3 r1 = this->vertices[triangle[0]];
        Vector3 r2 = this->vertices[triangle[1]];
        Vector3 r3 = this->vertices[triangle[2]];
        Vector3 r21 = r2 - r1;
        Vector3 r31 = r3 - r1;
        double r21_norm = r21.norm();
        r21.normalize();
        this->q_l_g.row(row_count) << r21_norm, r31.dot(r21), r31.cross(r21).norm();
        row_count ++;
    }
}

void LscmRelax::set_q_l_m()
{
    // get the coordinates of a triangle in local coordinates from the 2d map
    // x1, y1, y2 = 0
    // -> vector<x2, x3, y3>
    this->q_l_m.resize(this->triangles.size(), 3);
    unsigned int row_count = 0;
    for (auto triangle: this->triangles)
    {
        Vector2 r1 = this->flat_vertices[triangle[0]];
        Vector2 r2 = this->flat_vertices[triangle[1]];
        Vector2 r3 = this->flat_vertices[triangle[2]];
        Vector2 r21 = r2 - r1;
        Vector2 r31 = r3 - r1;
        double r21_norm = r21.norm();
        r21.normalize();
        this->q_l_g.row(row_count) << r21_norm, r31.dot(r21), r31.x() * r21.y() - r31.y() * r21.x();
        row_count ++;
    }
}

void LscmRelax::set_fixed_pins()
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
    std::sort(this->fixed_pins.begin(), this->fixed_pins.end());
    if (this->fixed_pins.size() > 2)
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

std::vector<Vector2> LscmRelax::get_flat_vertices()
{
    return this->flat_vertices;
}

std::vector<Vector3> LscmRelax::get_flat_vertices_3D()
{
    std::vector<Vector3> flat_3D;
    for (auto vertex: this->flat_vertices)
        flat_3D.push_back(Vector3(vertex.x(), vertex.y(), 0));
    return flat_3D;
}

std::vector<std::array<double, 2>> LscmRelax::get_flat_list()
{
    std::vector<std::array<double, 2>> flat;
    for (auto vertex: this->flat_vertices)
        flat.push_back(std::array<double, 2>{vertex.x(), vertex.y()});
    return flat;
}

std::vector<std::array<double, 3>> LscmRelax::get_flat_list_3D()
{
    std::vector<std::array<double, 3>> flat;
    for (auto vertex: this->flat_vertices)
        flat.push_back(std::array<double, 3>{vertex.x(), vertex.y(), 0});
    return flat;
}

}
