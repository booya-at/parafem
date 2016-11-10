#include "nurbs.h"

namespace paraFEM{


// DE BOOR ALGORITHM FROM OPENGLIDER
std::function<double(double)> get_basis(int degree, int i, Eigen::VectorXi knots)
    // Return a basis_function for the given degree """
{
    if (degree == 0)
    {
        return [degree, i, knots](double t)
        {
            // The basis function for degree = 0 as per eq. 7
            int t_this = knots[i];
            int t_next = knots[i+1];
            return (t_next >= t > t_this);
        };
    }
    else
    {
        return [degree, i, knots](double t)
        {
            // The basis function for degree > 0 as per eq. 8
            if (i == t == 0)
                return 1.;
            double out = 0.;
            int t_this = knots[i];
            int t_next = knots[i + 1];
            int t_precog  = knots[i + degree];
            int t_horizon = knots[i + degree + 1];

            double top = (t - t_this);
            double bottom = (t_precog - t_this);

            if (bottom != 0)
                out = top / bottom * get_basis(degree - 1, i, knots)(t);

            top = (t_horizon - t);
            bottom = (t_horizon - t_next);
            if (bottom != 0)
                out += top / bottom * get_basis(degree-1, i + 1, knots)(t);
            if (bottom == top == 0)
                out = 0;
            return out;
        };
    }
};


std::function<double(double)> get_basis_derivative(int order, int degree, int i, Eigen::VectorXi knots)
    // Return the derivation of the basis function """
{
    if (order == 1)
    {
        return [degree, i, knots](double t)
        {
            double out = 0;
            out +=  get_basis(degree - 1, i, knots)(t) *
                    degree / (knots[i + degree] - knots[i]);
            out -=  get_basis(degree - 1, i + 1, knots)(t) *
                    degree / (knots[i + degree + 1] - knots[i + 1]);
            return out;
        };
    }
    else
    {   return [degree, i, knots, order](double t)
        {
            double out = 0;
            out +=  get_basis_derivative(order, degree - 1, i, knots)(t) *
                    degree / (knots[i + degree] - knots[i]);
            out -=  get_basis_derivative(order, degree - 1, i + 1, knots)(t) *
                    degree / (knots[i + degree + 1] - knots[i + 1]);
            return out;
        };
    }
}


NurbsBase::NurbsBase(Eigen::VectorXi u_knots, Eigen::VectorXi v_knots,
                     Eigen::VectorXd weights,
                     int degree_u, int degree_v)
{
    assert(weights.size() == u_knots.size() * v_knots.size());
    this->u_knots = u_knots;
    this->v_knots = v_knots;
    this->weights = weights;
    this->degree_u = degree_u;
    this->degree_v = degree_v;
    for (int u_i = 0; u_i < u_knots.size(); u_i ++)
        this->u_functions.push_back(get_basis(degree_u, u_i, u_knots));
    for (int v_i = 0; v_i < v_knots.size(); v_i ++)
        this->u_functions.push_back(get_basis(degree_v, v_i, v_knots));
}

Eigen::VectorXd NurbsBase::getInfluenceVector(Eigen::Vector2d u)
{
    double n_u, n_v;
    double sum_weights = 0;
    Eigen::VectorXd infl(this->u_functions.size() * this->v_functions.size());
    int i = 0;
    int u_i = 0;
    int v_i = 0;
    for (auto foo_u: this->u_functions)
    {
        n_u = foo_u(u.x());
        v_i = 0;
        for (auto foo_v: this->v_functions)
        {
            n_v = foo_v(u.y());
            sum_weights += weights[i] * n_u * n_v;
            infl[i] = weights[i] * n_u * n_v;
            i ++;
            v_i ++;
        }
        u_i ++;
    }
    return infl / sum_weights;
}

void add_triplets(Eigen::VectorXd values, double row, std::vector<trip> &triplets)
{
    for (int i=0; i < values.size(); i++)
    {
        if (values(i) != 0.)
            triplets.push_back(trip(row, i, values(i)));
    }
}

spMat NurbsBase::getInfluenceMatrix(Eigen::Matrix<double, Eigen::Dynamic, 2> U)
{
    std::vector<trip> triplets;
    for (int row_index; row_index < U.rows(); row_index++)
        add_triplets(this->getInfluenceVector(U.row(row_index)), row_index, triplets);
    spMat mat(U.rows(), this->u_functions.size() * this->v_functions.size());
    mat.setFromTriplets(triplets.begin(), triplets.end());
    return mat;
}

void NurbsBase::computeFirstDerivatives()
{
    for (int u_i = 0; u_i < u_knots.size(); u_i ++)
        this->Du_functions.push_back(get_basis_derivative(1, this->degree_u, u_i, this->u_knots));
    for (int v_i = 0; v_i < u_knots.size(); v_i ++)
        this->Dv_functions.push_back(get_basis_derivative(1, this->degree_v, v_i, this->v_knots));
}

// void NurbsBase::computeSecondDerivatives()
// {
//     for (int u_i = 0; u_i < u_knots.size(); u_i ++)
//         this->DDu_functions.push_back(get_basis_derivative(2, this->degree_u, u_i, this->u_knots));
//     for (int v_i = 0; v_i < u_knots.size(); v_i ++)
//         this->DDv_functions.push_back(get_basis_derivative(2, this->degree_v, v_i, this->v_knots));
// }

Eigen::VectorXd NurbsBase::getDuVector(Eigen::Vector2d u)
{
    Eigen::VectorXd A1, A2;
    double C1, C2;
    A1.resize(this->u_functions.size() * this->v_functions.size());
    A2.resize(this->u_functions.size() * this->v_functions.size());
    double A3 = 0;
    double A5 = 0;
    int i = 0;
    int u_i = 0;
    int v_i = 0;
    Eigen::VectorXd n_u, n_v, Dn_u;
    n_u.resize(this->u_functions.size());
    Dn_u.resize(this->v_functions.size());
    n_v.resize(this->v_functions.size());
    for (int u_i=0; u_i < this->u_functions.size(); u_i++)
    {
        n_u[u_i] = this->u_functions[u_i](u.x());
        Dn_u[u_i] = this->Du_functions[u_i](u.x());
    }
    for (int v_i=0; v_i < this->v_functions.size(); v_i++)
    {
        n_v[v_i] = this->v_functions[v_i](u.y());
    }

    for (int u_i=0; u_i < this->u_functions.size(); u_i++)
    {
        v_i = 0;
        for (int v_i=0; v_i < this->v_functions.size(); v_i++)
        {
            C1 = weights[i] * Dn_u[u_i] * n_v[v_i];
            C2 = weights[i] * n_u[u_i] * n_v[v_i];
            A1[i] = C1; A2[i] = C1;
            A3 += C2; A5 += C1;
            i ++;
            v_i ++;
        }
        u_i ++;
    }
    return A1 / A3 - A2 * A5 / A3;
}

Eigen::VectorXd NurbsBase::getDvVector(Eigen::Vector2d u)
{
    Eigen::VectorXd A1, A2;
    double C1, C2;
    A1.resize(this->u_functions.size() * this->v_functions.size());
    A2.resize(this->u_functions.size() * this->v_functions.size());
    double A3 = 0;
    double A5 = 0;
    int i = 0;
    int u_i = 0;
    int v_i = 0;
    Eigen::VectorXd n_u, n_v, Dn_v;
    n_u.resize(this->u_functions.size());
    Dn_v.resize(this->v_functions.size());
    n_v.resize(this->v_functions.size());
    for (int u_i=0; u_i < this->u_functions.size(); u_i++)
    {
        n_u[u_i] = this->u_functions[u_i](u.x());
    }
    for (int v_i=0; v_i < this->v_functions.size(); v_i++)
    {
        n_v[v_i] = this->v_functions[v_i](u.x());
        Dn_v[v_i] = this->Dv_functions[v_i](u.y());
    }

    for (int u_i=0; u_i < this->u_functions.size(); u_i++)
    {
        v_i = 0;
        for (int v_i=0; v_i < this->v_functions.size(); v_i++)
        {
            C1 = weights[i] * Dn_v[v_i] * n_u[u_i];
            C2 = weights[i] * n_v[v_i] * n_u[u_i];
            A1[i] = C1; A2[i] = C1;
            A3 += C2; A5 += C1;
            i ++;
            v_i ++;
        }
        u_i ++;
    }
    return A1 / A3 - A2 * A5 / A3;
}


spMat NurbsBase::getDuMatrix(Eigen::Matrix<double, Eigen::Dynamic, 2> U)
{
    std::vector<trip> triplets;
    for (int row_index; row_index < U.rows(); row_index++)
        add_triplets(this->getDuVector(U.row(row_index)), row_index, triplets);
    spMat mat(U.rows(), this->u_functions.size() * this->v_functions.size());
    mat.setFromTriplets(triplets.begin(), triplets.end());
    return mat;
}

spMat NurbsBase::getDvMatrix(Eigen::Matrix<double, Eigen::Dynamic, 2> U)
{
    std::vector<trip> triplets;
    for (int row_index; row_index < U.rows(); row_index++)
        add_triplets(this->getDvVector(U.row(row_index)), row_index, triplets);
    spMat mat(U.rows(), this->u_functions.size() * this->v_functions.size());
    mat.setFromTriplets(triplets.begin(), triplets.end());
    return mat;
}

}