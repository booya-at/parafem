#ifndef ELEMENT_H
#define ELEMENT_H

#include "base.h"
#include "node.h"
#include "utils.h"
#include "material.h"
#include "Eigen/Geometry"
#include <Eigen/IterativeLinearSolvers>
#include <memory>
#include <vector>

namespace paraFEM
{

typedef Eigen::Triplet<double> trip;

Eigen::Matrix2Xd toMatrix(std::vector<Vector2>);
Eigen::Matrix2d findRotMat(std::vector<Vector2>, std::vector<Vector2>);
Eigen::Matrix2d findRotMat(Eigen::Matrix2Xd in, Eigen::Matrix2Xd out);


struct CoordSys: Base
{
    CoordSys(){};
    CoordSys(Vector3);
    CoordSys(Vector3, Vector3);

    Eigen::Matrix<double, 2, 3> mat;            // matrix to transform global(3d) to local(2d)
    Vector3 n, t1, t2;
    void rotate(std::vector<Vector2>, std::vector<Vector2>);
    void rotate(Eigen::MatrixX2d, Eigen::MatrixX2d);
    void update(Vector3);
    
    void update(Vector3, Vector3);
    Vector2 toLocal(Vector3);
    Vector3 toGlobal(Vector2);
    
};

struct IntegrationPoint: Base
{
    IntegrationPoint(double eta, double zeta, double weight);
    double eta;
    double zeta;
    double weight;
    Vector3 stress;
};

struct Element: Base
{
    std::vector<NodePtr> nodes;
    virtual void geometryStep() = 0;
    virtual void explicitStep(double h) = 0;  // compute the internal forces acting on the nodes.
    // virtual void implicitStep(double h) = 0;
    std::vector<int> getNr();
    virtual Vector3 getStress()=0;
    bool is_valid = true;
    double dViscous;
    double characteristicLength;
    virtual MaterialPtr getMaterial() = 0;
};

struct Truss: public Element
{
    Vector3 pressure;
    Vector3 tangent;
    double length;
    Truss(const std::vector<NodePtr>, std::shared_ptr<TrussMaterial>);
    virtual Vector3 getStress();
    double stress;           // at timestep n
    virtual void geometryStep();
    virtual void explicitStep(double h);
    virtual void implicitStep(std::vector<trip> & Kt);
    std::shared_ptr<TrussMaterial> material;
    void addNodalPressure(Vector3);
    MaterialPtr getMaterial();
};

struct LineJoint: public Element
{
// elment that connects nodes that are not exactly matching
    LineJoint(const std::vector<NodePtr>, std::shared_ptr<TrussMaterial>);
    virtual void geometryStep();
    virtual void explicitStep(double h);
    virtual Vector3 getStress();
    std::shared_ptr<TrussMaterial> material;
};

struct Membrane: public Element
{
    Vector3 center;
    CoordSys coordSys;
    double area;
    double pressure;                 // pressure acting on internal forces

    Vector3 calculate_center();
    
    std::shared_ptr<MembraneMaterial> material;
    MaterialPtr getMaterial();
};

struct Membrane3: public Membrane
{
    Membrane3(const std::vector<NodePtr>, std::shared_ptr<MembraneMaterial>);
    virtual void geometryStep();
    virtual void explicitStep(double h);
    Vector3 stress;
    virtual Vector3 getStress();

    Eigen::Matrix<double, 3, 2> pos_mat, vel_mat;
    Eigen::Matrix<double, 2, 3> B, dN;
};

struct Membrane4: public  Membrane
{
    Membrane4(const std::vector<NodePtr>, std::shared_ptr<MembraneMaterial>, bool reduced_integration=true);
    std::vector<IntegrationPoint> integration_points;  //eta, zeta, weight, stress
    virtual void geometryStep();
    virtual void explicitStep(double h);
    
    // hourglass control
    void initHG();
    Eigen::Vector4d hg_gamma;
    double hg_const;
    Vector2 hg_stress;
    
    
    Eigen::Matrix<double, 4, 2> pos_mat, vel_mat;
    Eigen::Matrix<double, 2, 4> B, dN;
    virtual Vector3 getStress();
    
};

typedef std::shared_ptr<Truss> TrussPtr;
typedef std::shared_ptr<Membrane> MembranePtr;
typedef std::shared_ptr<Membrane3> Membrane3Ptr;
typedef std::shared_ptr<Membrane4> Membrane4Ptr;
typedef std::shared_ptr<Element> ElementPtr;
typedef std::vector<ElementPtr> ElementVec;
}

#endif