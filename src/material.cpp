#include "material.h"


namespace paraFEM {
    
    
TrussMaterial::TrussMaterial(double elasticity)
{
    this->elasticity = elasticity;
}

MembraneMaterial::MembraneMaterial(double elasticity, double nue)
{
    this->elasticity = elasticity;
    this->nue = nue;
    this->C << 1,   nue, 0,
                nue, 1,  0,
                0,   0, (1 - nue) / 2;
    this->C *= elasticity / (1 - nue * nue);
}



}