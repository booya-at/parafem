#ifndef ELEMENT_H
#define ELEMENT_H

#include "base.h"
#include "node.h"
#include "material.h"
#include "eigen3/Eigen/Geometry"
#include <memory>
#include <vector>

namespace paraFEM
{


Eigen::Matrix2Xd toMatrix(std::vector<Vector2>);
Eigen::Matrix2d findRotMat(std::vector<Vector2>, std::vector<Vector2>);
Eigen::Matrix2d findRotMat(Eigen::Matrix2Xd in, Eigen::Matrix2Xd out);


class CoordSys: Base
{
public:
    CoordSys(){};
    CoordSys(Vector3);
    CoordSys(Vector3, Vector3);

    Eigen::Matrix<double, 2, 3> mat;            // matrix to transform global(3d) to local(2d)
    Vector3 n, t1, t2;
    void rotate(std::vector<Vector2>, std::vector<Vector2>);
    void update(Vector3);
    Vector2 toLocal(Vector3);
    Vector3 toGlobal(Vector2);
};

class Element: Base
{
public:
    std::vector<NodePtr> nodes;
    virtual void makeStep(double h) = 0;  // compute the internal forces acting on the nodes.
    std::vector<int> getNr();
};

class Truss: public Element
{
public:
    Vector3 pressure;
    double length;
    Truss(std::vector<NodePtr>, std::shared_ptr<TrussMaterial>);
    double stress;           // at timestep n
    virtual void makeStep(double h);
    std::shared_ptr<TrussMaterial> material;
    void addNodalPressure(Vector3);
};

class Membrane: public Element
{
public:
    Vector3 center;
    CoordSys coordSys;
    double area;
    double pressure;                 // pressure acting on internal forces
    void setConstPressure(double);
    std::vector<Vector2> position;   // storing the local position of all points
    Vector3 stress;
    std::vector<Vector2> local_position;
    std::vector<Vector2> local_velocity;
    std::shared_ptr<MembraneMaterial> material;
    void addNodalPressure(double);
};

class Membrane3: public Membrane
{
public:
    Membrane3(std::vector<NodePtr>, std::shared_ptr<MembraneMaterial>);
    virtual void makeStep(double h);
};

class Membrane4: public  Membrane
{
public:
    Membrane4(std::vector<NodePtr>, std::shared_ptr<MembraneMaterial>);
    virtual void makeStep(double h);
    void initHG();
    Eigen::Vector4d hg_gamma;
    double hg_const;
    Vector2 hg_stress;
    void addNodalPressure(double);
    
};

typedef std::shared_ptr<Truss> TrussPtr;
typedef std::shared_ptr<Membrane> MembranePtr;
typedef std::shared_ptr<Membrane3> Membrane3Ptr;
typedef std::shared_ptr<Membrane4> Membrane4Ptr;
typedef std::vector<std::shared_ptr<Element>> ElementVec;
}
#endif