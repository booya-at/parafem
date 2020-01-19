#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include "base.h"
#include "node.h"
#include "material.h"
#include "Eigen/Geometry"
#include <memory>
#include <vector>
#include <tuple>

namespace parafem
{

struct ForceConstraint
{
	virtual void update_force() = 0;
};

struct PositionConstraint
{
	virtual void update_position() = 0;
};

struct FixConstraint: public PositionConstraint
{
	std::vector<NodePtr> points;
	FixConstraint(std::vector<NodePtr>, Eigen::Vector3d fix);
	Eigen::Vector3d fix;
	virtual void update_position() = 0;
};

struct SymmetricPosition: public PositionConstraint
{
	// updates the position of the symmetric points.

	std::vector<std::array<NodePtr, 2>> sym_points;
	Eigen::Vector3d mirror_point, mirror_normal;
	SymmetricPosition(std::vector<NodePtr> points, Eigen::Vector3d mirror_point, Eigen::Vector3d mirror_normal);
	virtual void update_position() = 0;
};

struct SymmetricForce: public ForceConstraint
{
	// applying the mirrored internal force
	std::vector<NodePtr> sym_points; 
	Eigen::Vector3d mirror_point, mirror_normal;
	SymmetricForce(std::vector<NodePtr> points, Eigen::Vector3d mirror_point, Eigen::Vector3d mirror_normal);
	virtual void update_force() = 0;
};


typedef std::shared_ptr<PositionConstraint> PositionConstraintPtr;
typedef std::shared_ptr<ForceConstraint> ForceConstraintPtr;

}

#endif