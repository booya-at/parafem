#ifndef MATERIAL_H
#define MATERIAL_H

#include "base.h"
#include "Eigen/Core"
#include <memory>

namespace paraFEM{

class Material: Base
{
public:
    double elasticity = 10000.;
    double rho = 1.;
    double d_structural = 0.0001;
    double d_velocity = 1.;
    double d_minMode = 0.000;
    virtual double waveSpeed() = 0;
};


class TrussMaterial: public Material
{
public:
    TrussMaterial(double elasticity);
    virtual double waveSpeed();
};


class MembraneMaterial: public Material
{
public:
    MembraneMaterial(double elasticity, double nue);
    double nue;
    Eigen::Matrix3d C;
    virtual double waveSpeed();
};

typedef std::shared_ptr<Material> MaterialPtr;
typedef std::shared_ptr<MembraneMaterial> MembraneMaterialPtr;
typedef std::shared_ptr<TrussMaterial> TrussMaterialPtr;

}


#endif