#include "base.h"
#include "eigen3/Eigen/Core"
#include <memory>

namespace paraFEM{

class Material: Base
{
public:
    double rho = 1;
    double d_structural = 1;
    double elasticity = 10000;
};


class TrussMaterial: public Material
{
public:
    TrussMaterial(double elasticity);
};


class MembraneMaterial: public Material
{
public:
    MembraneMaterial(double elasticity, double nue);
    double nue;
    Eigen::Matrix3d C;
};

typedef std::shared_ptr<MembraneMaterial> MembraneMaterialPtr;
typedef std::shared_ptr<TrussMaterial> TrussMaterialPtr;

}