#include "material.h"


namespace parafem {


double TrussMaterial::waveSpeed()
{
	return pow(this->elasticity / this->rho, 0.5);
}


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


double MembraneMaterial::waveSpeed()
{
	return pow(this->elasticity / this->rho / (1 - this->nue), 0.5);
}

}