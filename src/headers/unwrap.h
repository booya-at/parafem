// LeastSquareConformalMapping + fem relaxing
// ------------------------------------------
// 
#ifndef UNWRAP_H
#define UNWRAP_H

// 1: local coordinates 2d representation  q_l_0
// 2: least square conformal map -> flat_vertices_0
// 3: local coordinates of mapped mesh q_l_1
// 4: diff in local coordinates -> forces R.B^T.B.(x1-x0)
// 5: stiffnes mat K
// 6: K.u=forces ->u
// 7: x1, y1 += w * u

#include <vector>
#include <memory>
#include <tuple>

#include "material.h"
#include "Eigen/Geometry"

namespace paraFEM{

typedef Eigen::Vector3d Vector3;
typedef Eigen::Vector2d Vector2;

class LscmRelax{
private:
    Eigen::MatrixXd q_l_g;              // the position of the 3d triangles at there locale coord sys
    Eigen::MatrixXd q_l_m;              // the mapped position in local coord sys
    void set_q_l_g();
    void set_q_l_m();
    std::vector<Vector2> flat_vertices;
    std::vector<int> fixed_pins;

    void set_fixed_pins();

    std::vector<int> new_vertex_order;
    std::vector<int> old_vertex_order;

    unsigned int get_new_order(unsigned int i);
    unsigned int get_old_order(unsigned int i);

    void init(
        std::vector<Vector3> vertices, 
        std::vector<std::array<int, 3>> triangles,
        std::vector<int> fixed_pins);
public:
    // set fixed pins (index pos)
    // create new vertex order (fixed pins are send to the end)
    LscmRelax(
        std::vector<Vector3> vertices, 
        std::vector<std::array<int, 3>> triangles,
        std::vector<int> fixed_pins={});

    LscmRelax(
        std::vector<std::array<double, 3>> vertices, 
        std::vector<std::array<int, 3>> triangles,
        std::vector<int> fixed_pins={});

    std::vector<Vector3> vertices;
    std::vector<std::array<int, 3>> triangles;


    Eigen::Matrix<int, Eigen::Dynamic, 1> _fixed_pins;
    Eigen::Matrix<double, Eigen::Dynamic, 3> _vertices;
    Eigen::Matrix<long, Eigen::Dynamic, 3> _triangles;
    Eigen::Matrix<double, Eigen::Dynamic, 2> _flat_vertices;

//////////////////////////////////// TODO: remove this functions
    std::vector<Vector2> get_flat_vertices();
    std::vector<Vector3> get_flat_vertices_3D();
    std::vector<std::array<double, 2>> get_flat_list();
    std::vector<std::array<double, 3>> get_flat_list_3D();
//////////////////////////////////////////////////////////////

    MembraneMaterial mat = MembraneMaterial(1000, 0.5); // elasticity, nue

    void lscm();
    void relax(double step_size);
};

template <typename T>
Eigen::MatrixXd vec2mat(std::vector<T> vector_matrix)
{
    Eigen::MatrixXd mat(vector_matrix.size(), vector_matrix[0].size());
    for (auto vec: vector_matrix)
        mat << vec;
    return mat;
}

template <typename T>
std::vector<T> mat2vec(Eigen::MatrixXd mat)
{
    std::vector<T> vector_matrix;
    for (int i = 0; i < mat.rows(); i++)
        vector_matrix.push_back(T(mat.row(i)));
    return vector_matrix;
}

}


#endif