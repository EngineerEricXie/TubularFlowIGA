#ifndef IGA_NATIVE_TET_ALE_BOUNDARY_HPP
#define IGA_NATIVE_TET_ALE_BOUNDARY_HPP

#include "NativeTetFem.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <stdexcept>
#include <string>

namespace iga {

enum class NativeTetAleBoundaryMode
{
	PrescribedDisplacement,
	Fixed,
	SlidingNormalConstraint
};

struct NativeTetAleBoundaryRule
{
	int label = -1;
	NativeTetAleBoundaryMode mode = NativeTetAleBoundaryMode::Fixed;
};

struct NativeTetAleBoundaryConstraints
{
	std::map<std::uint32_t,std::array<double,3>> full_displacement_m;
	std::set<int> applied_labels;
};

// Convert label-level policy into explicit nodal Dirichlet data. Sliding is a
// declared mode, but fails closed until the harmonic solver has normal MPCs.
inline NativeTetAleBoundaryConstraints ResolveNativeTetAleBoundaryConstraints(
	const NativeTetMesh& mesh,const std::vector<NativeTetAleBoundaryRule>& rules,
	const std::map<std::uint32_t,std::array<double,3>>& prescribed_displacement_m,
	double conflict_tolerance_m = 1.0e-12)
{
	if (rules.empty() || !(conflict_tolerance_m >= 0.0)
		|| !std::isfinite(conflict_tolerance_m))
		throw std::invalid_argument("native ALE boundary policy is invalid");
	std::map<int,NativeTetAleBoundaryMode> modes;
	for (const auto& rule : rules) {
		if (rule.label < 0 || !modes.emplace(rule.label,rule.mode).second)
			throw std::invalid_argument("native ALE boundary labels must be unique and nonnegative");
	}
	NativeTetAleBoundaryConstraints result;
	for (const auto& triangle : mesh.boundary_triangles) {
		const auto found = modes.find(triangle.boundary_label);
		if (found == modes.end())
			throw std::runtime_error("native ALE boundary label has no explicit rule");
		result.applied_labels.insert(triangle.boundary_label);
		if (found->second == NativeTetAleBoundaryMode::SlidingNormalConstraint)
			throw std::runtime_error("native ALE sliding-normal MPC is declared but not implemented");
		for (const auto node : triangle.nodes) {
			if (node >= mesh.points.size())
				throw std::invalid_argument("native ALE boundary node is out of range");
			std::array<double,3> value{{0,0,0}};
			if (found->second == NativeTetAleBoundaryMode::PrescribedDisplacement) {
				const auto prescribed = prescribed_displacement_m.find(node);
				if (prescribed == prescribed_displacement_m.end())
					throw std::runtime_error("native ALE moving boundary node has no displacement");
				value = prescribed->second;
			}
			for (const auto component : value)
				if (!std::isfinite(component))
					throw std::invalid_argument("native ALE boundary displacement is nonfinite");
			const auto previous = result.full_displacement_m.find(node);
			if (previous == result.full_displacement_m.end())
				result.full_displacement_m.emplace(node,value);
			else
				for (int component = 0; component < 3; ++component)
					if (std::abs(previous->second[component]-value[component])
						> conflict_tolerance_m)
						throw std::runtime_error(
							"native ALE boundary rules conflict at a shared node");
		}
	}
	for (const auto& rule : rules)
		if (result.applied_labels.count(rule.label) == 0)
			throw std::runtime_error("native ALE boundary rule names an absent label");
	return result;
}

} // namespace iga

#endif
