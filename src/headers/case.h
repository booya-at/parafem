#ifndef CASE_H
#define CASE_H

#include <memory>
#include <vector>
#include <set>

#include "base.h"
#include "node.h"
#include "element.h"
#include "constraint.h"

#include "Eigen/Core"

namespace paraFEM
{

class FemCase: Base
{
public:
    double time=0;
    FemCase(std::vector<ElementPtr>);
	FemCase(std::vector<ElementPtr> elements, 
            std::vector<PositionConstraintPtr>,
            std::vector<ForceConstraintPtr>);
    std::set<NodePtr> nodes;
    std::vector<ElementPtr> elements;
    std::vector<ForceConstraintPtr> forceConstraints;
    std::vector<PositionConstraintPtr> posConstraints;
    void explicitStep(double h, double externalFactor=1);
    void implicitStep(double h, double externalFactor=1);

    NodeVec get_nodes();

};

typedef std::shared_ptr<FemCase> FemCasePtr;

};
#endif