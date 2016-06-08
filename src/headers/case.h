#ifndef CASE_H
#define CASE_H

#include <memory>
#include <vector>
#include <set>

#include "base.h"
#include "node.h"
#include "element.h"

#include "eigen3/Eigen/Core"

namespace paraFEM
{

class FemCase: Base
{
public:
    double time=0;
    FemCase(std::vector<ElementPtr>);
    std::set<NodePtr> nodes;
    ElementVec elements;
    void makeStep(double h, double externalFactor=1);

    NodeVec get_nodes();

};

typedef std::shared_ptr<FemCase> FemCasePtr;

};
#endif