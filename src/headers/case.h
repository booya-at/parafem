#ifndef CASE_H
#define CASE_H

#include <memory>
#include <vector>
#include <set>

#include "base.h"
#include "node.h"
#include "element.h"

namespace paraFEM
{

class FemCase: Base
{
public:
    double time=0;
    FemCase(std::vector<std::shared_ptr<Element>>);
    std::set<NodePtr> nodes;
    ElementVec elements;
    void makeStep(double h);
    
    double d_velocity = 0;
};

typedef std::shared_ptr<FemCase> FemCasePtr;

};
#endif