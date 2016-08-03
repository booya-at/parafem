#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include "base.h"
#include "node.h"
#include "material.h"
#include "eigen3/Eigen/Geometry"
#include <memory>
#include <vector>
#include <tuple>

namespace paraFEM
{

struct ForceConstraint
{
	virtual void updateForce() = 0;
};

struct PositionConstraint
{
	virtual void updatePosition() = 0;
};

struct FixConstraint: public PositionConstraint
{
	std::vector<NodePtr> points;
	FixConstraint(std::vector<NodePtr>, Eigen::Vector3d fix);
	Eigen::Vector3d fix;
	virtual void updatePosition() = 0;
};

struct SymmetricPosition: public PositionConstraint
{
	// updates the position of the symmetric points.

	std::vector<std::array<NodePtr, 2>> symPoints;
	Eigen::Vector3d mirrorPoint, mirrorNormal;
	SymmetricPosition(std::vector<NodePtr> points, Eigen::Vector3d mirrorPoint, Eigen::Vector3d mirrorNormal);
	virtual void updatePosition() = 0;
};

struct SymmetricForce: public ForceConstraint
{
	// applying the mirrored internal force
	std::vector<NodePtr> symPoints; 
	Eigen::Vector3d mirrorPoint, mirrorNormal;
	SymmetricForce(std::vector<NodePtr> points, Eigen::Vector3d mirrorPoint, Eigen::Vector3d mirrorNormal);
	virtual void updateForce() = 0;
};


typedef std::shared_ptr<PositionConstraint> PositionConstraintPtr;
typedef std::shared_ptr<ForceConstraint> ForceConstraintPtr;

}

#endif