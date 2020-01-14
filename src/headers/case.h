#ifndef CASE_H
#define CASE_H

#include <memory>
#include <vector>
#include <set>
#include <tuple>

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
    std::vector<ForceConstraintPtr> force_constraints;
    std::vector<PositionConstraintPtr> pos_constraints;
    void explicit_step(double h, double external_factor=1);
    void implicit_step(double h, double external_factor=1);
    double getMaxVelocity();
    std::tuple<double, ElementPtr> get_explicit_max_time_step();
    NodeVec get_nodes();

};

typedef std::shared_ptr<FemCase> FemCasePtr;

};
#endif