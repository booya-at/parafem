// LeastSquareConformalMapping + fem relaxing
// ------------------------------------------
// 
// 1: local coordinates 2d representation  q_l_0
// 2: least square conformal map -> flat_vertices_0
// 3: local coordinates of mapped mesh q_l_1
// 4: diff in local coordinates -> forces R.B^T.B.(x1-x0)
// 5: stiffnes mat K
// 6: K.u=forces ->u
// 7: x1, y1 += w * u



#ifndef UNWRAP_H
#define UNWRAP_H

#include "material.h"

#include <vector>
#include <memory>
#include <tuple>
#include "Eigen/Geometry"
namespace paraFEM{

typedef Eigen::Vector3d Vector3;
typedef Eigen::Vector2d Vector2;

class LscmRelax{
private:
    Eigen::MatrixXd q_l_0;              // the position of the 3d triangles at there locale coord sys
    Eigen::MatrixXd q_l_1;              // the mapped position in local coord sys
    void LscmRelax::set_q_l_0();
    void LscmRelax::set_q_l_1();
    std::vector<Vector2> flat_vertices;

    std::vector<unsigned int> fixed_pins;
    void set_fixed_pins();

    std::vector<int> new_vertex_order;
    std::vector<int> old_vertex_order;

    unsigned int get_new_order(unsigned int i);
    unsigned int get_old_order(unsigned int i);
public:
    // set fixed pins (index pos)
    // create new vertex order (fixed pins are send to the end)
    LscmRelax(
        std::vector<Eigen::Vector3D> vertices, 
        std::vector<std::array<3, int>> triangles,
        std::vector<unsigned int> fixed_pins_indices=NULL);
    
    std::vector<Eigen::Vector3D> vertices;
    std::vector<std::array<3, int>> triangles;
    std::vector<Vector2> get_flat_vertices();

    MembraneMaterial mat = MembraneMaterial(1000, 0.5); // elasticity, nue

    void lscm();
    void relax(double step_size);
};

template <typename T>
Eigen::MatrixXd mat_from_vector(std::vector<T> vector_matrix)
{
    Eigen::MatrixXd mat(vector_matrix.size(), vector_matrix[0].size());
    for (auto vec: vector_matrix)
        mat << vec_;
}

template <typename T>
std::vector<T> vector_from_matrix(Eigen::MatrixXd mat)
{
    std::vector<T> vector_matrix;
    for (int i == 0; i < mat.rows(); i++)
        vector_matrix.push_back(T(mat[i]))
}

}


#endif