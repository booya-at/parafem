#include "unwrap.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SVD>
#include <iostream>
#include <algorithm>
#include <cmath>

// TODO:
// vertices, flat vertices = EIgen::MatrixXd (better python conversation, parallell)
// make it a own library
// area constrained (scale the final unwrapped mesh to the original area)
// FEM approach

namespace paraFEM{

typedef Eigen::Triplet<double> trip;
typedef Eigen::SparseMatrix<double> spMat;



ColMat<double, 2> map_to_2D(ColMat<double, 3> points)
{
    ColMat<double, 4> mat(points.size(), 4);
    mat << points, ColMat<double, 1>::Ones(points.size());
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
    return points * transform;
}


unsigned int get_max_distance(Vector3 point, RowMat<double, 3> vertices, double & max_dist)
{
    max_dist = 0;
    double dist;
    unsigned long max_dist_index;
    unsigned long j = 0;
    for (j=0; j < vertices.cols(); j++)
    {
        // debugging
        dist = (point - vertices.col(j)).norm();
        if (dist > max_dist)
        {
            max_dist = dist;
            max_dist_index = j;
        }
    }
    return max_dist_index;
}

void LscmRelax::init(
        RowMat<double, 3> vertices, 
        RowMat<long, 3> triangles,
        std::vector<long> fixed_pins)
{
    this->vertices = vertices;
    this->triangles = triangles;
    this->flat_vertices.resize(2, this->vertices.cols());
    this->fixed_pins = fixed_pins;

    // set the fixed pins of the flat-mesh:
    this->set_fixed_pins();

    
    int fixed_count = 0;
    for (long i=0; i < this->vertices.cols(); i++)
    {   
        if (fixed_count < this->fixed_pins.size())
        {
            if (i == this->fixed_pins[fixed_count])
                fixed_count ++;
            else
                this->old_order.push_back(i);
        }
        else
            this->old_order.push_back(i);
    }

    for (auto fixed_index: this->fixed_pins)
        this->old_order.push_back(fixed_index);

    // get the reversed map:
    this->new_order.resize(this->old_order.size());
    long j = 0;
    for (auto index: this->old_order)
    {
        this->new_order[index] = j;
        j++;
    }

    // set the C-mat
    double rho = (1 - this->nue) / (1 - 2 * this->nue);
    double rho1 = this->nue / (1 - 2 * this->nue);
    this->C << 1,   nue, 0,
                nue, 1,  0,
                0,   0, (1 - nue) / 2;
    this->C *= elasticity / (1 - nue * nue);
}

LscmRelax::LscmRelax(
        RowMat<double, 3> vertices, 
        RowMat<long, 3> triangles,
        std::vector<long> fixed_pins)
{
    this->init(vertices, triangles, fixed_pins);
}

LscmRelax::LscmRelax(std::vector<std::array<double, 3>> vertices, 
                     std::vector<std::array<long, 3>> triangles,
                     std::vector<long> fixed_pins)
{

    RowMat<double, 3> verts(3, vertices.size());
    RowMat<long, 3> tris(3, triangles.size());
    std::array<double, 3> vert_i_vals; 
    std::array<long, 3> tri_i_vals;
    long i, j; 
    for (i=0; i < vertices.size(); i++)
    {
        vert_i_vals = vertices[i];
        for (j=0; j < 3; j++)
            verts(j, i) = vert_i_vals[j];
    }
    for (i=0; i < triangles.size(); i++)
    {
        tri_i_vals = triangles[i];
        for (j=0; j < 3; j++)
            tris(j, i) = tri_i_vals[j];
    }
    this->init(verts, tris, fixed_pins);
}

//////////////////////////////////////////////////////////////////////////
/////////////////                 F.E.M                      /////////////
//////////////////////////////////////////////////////////////////////////
void LscmRelax::relax(double weight)
{
    this->set_q_l_m();
    ColMat<double, 3> d_q_l_g = this->q_l_g - this->q_l_m;
    // for every triangle
    Eigen::Matrix<double, 3, 6> B;
    Eigen::Matrix<double, 2, 2> T;
    Eigen::Matrix<double, 6, 6> K_m;
    Eigen::Matrix<double, 6, 1> u_m, rhs_m;
    Eigen::VectorXd rhs(this->vertices.cols() * 2);
    Eigen::VectorXd sol(this->vertices.cols() * 2);
    spMat K_g(this->vertices.cols() * 2, this->vertices.cols() * 2);
    std::vector<trip> K_g_triplets;
    Vector2 v0, v1, v2, v12, v23, v31;
    long row_pos, col_pos;
    double A;

    rhs.setZero();

    for (long i=0; i<this->triangles.cols(); i++)
    {
        // 1: construct B-mat in m-system
        v0 = this->flat_vertices.col(this->triangles(0, i));
        v1 = this->flat_vertices.col(this->triangles(1, i));
        v2 = this->flat_vertices.col(this->triangles(2, i));
        v12 = v1 - v0;
        v23 = v2 - v1;
        v31 = v0 - v1;
        B << -v23.y(),   0,        -v31.y(),   0,        -v12.y(),   0,
              0,         v23.x(),   0,         v31.x(),   0,         v12.x(),
              v23.x(),  -v23.y(),   v31.x(),  -v31.y(),   v12.x(),  -v12.y();
        T << v12.x(), -v12.y(),
             v12.y(), v12.x();
        T /= v12.norm();
        A = std::abs(this->q_l_m(i, 0) * this->q_l_m(i, 2) / 2);
        B /= A * 2; // (2*area)

        // 2: sigma due dqlg in m-system
        u_m << Vector2(0, 0), T * Vector2(d_q_l_g(i, 0), 0), T * Vector2(d_q_l_g(i, 1), d_q_l_g(i, 2));

        // 3: rhs_m = B.T * C * B * dqlg_m
        //    K_m = B.T * C * B
        rhs_m = B.transpose() * this->C * B * u_m * A;
        K_m = B.transpose() * this->C * B * A;

        // 5: add to rhs_g, K_g
        for (int j=0; j < 3; j++)
        {
            row_pos = this->triangles(j, i);

            rhs[row_pos * 2]     += rhs_m[j * 2];
            rhs[row_pos * 2 + 1] += rhs_m[j * 2 +1];
            for (int k=0; k < 3; k++)
            {
                col_pos = this->triangles(k, i);
                K_g_triplets.push_back(trip(row_pos * 2,     col_pos * 2,        K_m(j * 2,      k * 2)));
                K_g_triplets.push_back(trip(row_pos * 2 + 1, col_pos * 2,        K_m(j * 2 + 1,  k * 2)));
                K_g_triplets.push_back(trip(row_pos * 2 + 1, col_pos * 2 + 1,    K_m(j * 2 + 1,  k * 2 + 1)));
                K_g_triplets.push_back(trip(row_pos * 2,     col_pos * 2 + 1,    K_m(j * 2,      k * 2 + 1)));
                // we don't have to fill all because the matrix is symetric.
            }
        }
        K_g.setFromTriplets(K_g_triplets.begin(), K_g_triplets.end());
    }
    // solve linear system (privately store the value for guess in next step)
    Eigen::ConjugateGradient<spMat, Eigen::Lower> solver;
    solver.compute(K_g);
    sol = solver.solve(rhs);
    // std::cout << "qlm" << this->q_l_m << std::endl;
    // std::cout << "qlg" << this->q_l_g << std::endl;
    std::cout << "u_m  " << u_m << std::endl;
    // std::cout << "rhs_m  " << rhs_m << std::endl;
    std::cout << "B  " << B << std::endl;
    // std::cout << "sol  " << sol << std::endl;
    // reset coordinates of the flat vertices, recalc the qlm-mat
    this->set_shift(sol * weight);
    this->set_q_l_m();
    this->sol = sol;
    this->rhs = rhs;
    this->MATRIX = K_g;

}


//////////////////////////////////////////////////////////////////////////
/////////////////                 L.S.C.M                    /////////////
//////////////////////////////////////////////////////////////////////////
void LscmRelax::lscm()
{
    this->set_q_l_g();
    std::vector<trip> triple_list;
    long i;
    double x21, x31, y31, x32;

    // 1. create the triplet list (t * 2, v * 2)
    for(i=0; i<this->triangles.cols(); i++)
    {
        x21 = this->q_l_g(i, 0);
        x31 = this->q_l_g(i, 1);
        y31 = this->q_l_g(i, 2);
        x32 = x31 - x21;

        triple_list.push_back(trip(2 * i, this->new_order[this->triangles(0, i)] * 2, x32));
        triple_list.push_back(trip(2 * i, this->new_order[this->triangles(0, i)] * 2 + 1, -y31));
        triple_list.push_back(trip(2 * i, this->new_order[this->triangles(1, i)] * 2, -x31));
        triple_list.push_back(trip(2 * i, this->new_order[this->triangles(1, i)] * 2 + 1, y31));
        triple_list.push_back(trip(2 * i, this->new_order[this->triangles(2, i)] * 2, x21));

        triple_list.push_back(trip(2 * i + 1, this->new_order[this->triangles(0, i)] * 2, y31));
        triple_list.push_back(trip(2 * i + 1, this->new_order[this->triangles(0, i)] * 2 + 1, x32));
        triple_list.push_back(trip(2 * i + 1, this->new_order[this->triangles(1, i)] * 2, -y31));
        triple_list.push_back(trip(2 * i + 1, this->new_order[this->triangles(1, i)] * 2 + 1, -x31));
        triple_list.push_back(trip(2 * i + 1, this->new_order[this->triangles(2, i)] * 2 + 1, x21));

    }
    // 2. divide the triplets in matrix(unknown part) and rhs(known part) and reset the position
    std::vector<trip> rhs_triplets;
    std::vector<trip> mat_triplets;
    for (auto triplet: triple_list)
    {
        if (triplet.col() > (this->vertices.cols() - this->fixed_pins.size()) * 2 - 1)
            rhs_triplets.push_back(triplet);
        else
            mat_triplets.push_back(triplet);
    }

    // 3. create a rhs_pos vector
    Eigen::VectorXd rhs_pos(this->vertices.cols() * 2);
    rhs_pos.setZero();
    for (auto index: this->fixed_pins)
    {
        rhs_pos[this->new_order[index] * 2] = this->flat_vertices(0, index);      //TODO: not yet set
        rhs_pos[this->new_order[index] * 2 + 1] = this->flat_vertices(1, index);
    }

    // 4. fill a sparse matrix and calculdate the rhs
    Eigen::VectorXd rhs(this->triangles.cols() * 2); // maybe use a sparse vector
    spMat B(this->triangles.cols() * 2, this->vertices.cols() * 2);
    B.setFromTriplets(rhs_triplets.begin(), rhs_triplets.end());
    rhs = B * rhs_pos;

    // 5. create the lhs matrix
    spMat A(this->triangles.cols() * 2, (this->vertices.cols() - this->fixed_pins.size()) * 2);
    A.setFromTriplets(mat_triplets.begin(), mat_triplets.end());

    // 6. solve the system and set the flatted coordinates
    // Eigen::SparseQR<spMat, Eigen::COLAMDOrdering<int> > solver;
    Eigen::LeastSquaresConjugateGradient<spMat > solver;
    Eigen::VectorXd sol(this->vertices.size() * 2);
    solver.compute(A);
    sol = solver.solve(-rhs);

    // TODO: create function, is needed also in the fem step
    this->set_position(sol);
    // 7. if size of fixed pins <= 2: scale the map to the same area as the 3d mesh

}

void LscmRelax::set_q_l_g()
{
    // get the coordinates of a triangle in local coordinates from the 3d mesh
    // x1, y1, y2 = 0
    // -> vector<x2, x3, y3>
    this->q_l_g.resize(this->triangles.cols(), 3);
    for (long i = 0; i < this->triangles.cols(); i++)
    {
        Vector3 r1 = this->vertices.col(this->triangles(0, i));
        Vector3 r2 = this->vertices.col(this->triangles(1, i));
        Vector3 r3 = this->vertices.col(this->triangles(2, i));
        Vector3 r21 = r2 - r1;
        Vector3 r31 = r3 - r1;
        double r21_norm = r21.norm();
        r21.normalize();
        this->q_l_g.row(i) << r21_norm, r31.dot(r21), r31.cross(r21).norm();
    }
}

void LscmRelax::set_q_l_m()
{
    // get the coordinates of a triangle in local coordinates from the 2d map
    // x1, y1, y2 = 0
    // -> vector<x2, x3, y3>
    this->q_l_m.resize(this->triangles.cols(), 3);
    for (long i = 0; i < this->triangles.cols(); i++)
    {
        Vector2 r1 = this->flat_vertices.col(this->triangles(0, i));
        Vector2 r2 = this->flat_vertices.col(this->triangles(1, i));
        Vector2 r3 = this->flat_vertices.col(this->triangles(2, i));
        Vector2 r21 = r2 - r1;
        Vector2 r31 = r3 - r1;
        double r21_norm = r21.norm();
        r21.normalize();
        this->q_l_m.row(i) << r21_norm, r31.dot(r21), r31.x() * r21.y() - r31.y() * r21.x();
    }
}

void LscmRelax::set_fixed_pins()
{
    // if less then one fixed pin is set find two by an automated algorithm and align them to y = 0
    // if more then two pins are choosen find a leastsquare-plane and project the points on it
    // insert the points in the flat-vertices vector
    if (this->fixed_pins.size() == 0)
        this->fixed_pins.push_back(0);
    if (this->fixed_pins.size() == 1)
    {

        double dist;
        this->fixed_pins.push_back(get_max_distance(this->vertices.col(this->fixed_pins[0]), this->vertices, dist));
        this->flat_vertices.col(this->fixed_pins[0]) = Vector2(0, 0);
        this->flat_vertices.col(this->fixed_pins[1]) = Vector2(dist, 0);
    }
    std::sort(this->fixed_pins.begin(), this->fixed_pins.end());
    // not yet working
    // if (this->fixed_pins.size() > 2)
    // {
    //     std::vector<Vector3> fixed_3d_points;
    //     for (unsigned int index: this->fixed_pins)
    //         fixed_3d_points.push_back(this->vertices[index]);
    //     std::vector<Vector2> maped_points = map_to_2D(fixed_3d_points);
    //     unsigned int i = 0;
    //     for (unsigned int index: this->fixed_pins)
    //     {
    //         this->flat_vertices[i] = maped_points[i];
    //         i++;
    //     }
    // }
}


ColMat<double, 3> LscmRelax::get_flat_vertices_3D()
{
    ColMat<double, 2> mat = this->flat_vertices.transpose();
    ColMat<double, 3> mat3d(mat.rows(), 3);
    mat3d << mat, ColMat<double, 1>::Zero(mat.rows());
    return mat3d;
}

void LscmRelax::set_position(Eigen::VectorXd sol)
{
    for (long i=0; i < this->vertices.size(); i++)
    {
        if (sol.size() > i * 2 + 1)
            this->flat_vertices.col(this->old_order[i]) << sol[i * 2], sol[i * 2 + 1];
    }
}

void LscmRelax::set_shift(Eigen::VectorXd sol)
{
    for (long i=0; i < this->vertices.size(); i++)
    {
        if (sol.size() > i * 2 + 1)
            this->flat_vertices.col(i) += Vector2(sol[i * 2], sol[i * 2 + 1]);
    }
}
}
