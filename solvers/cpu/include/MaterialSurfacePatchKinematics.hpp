#ifndef IGA_MATERIAL_SURFACE_PATCH_KINEMATICS_HPP
#define IGA_MATERIAL_SURFACE_PATCH_KINEMATICS_HPP

// Stateless composition of a trial P1 patch publication into a whole closed
// material surface.  No FSI runtime is introduced here; this is the narrow
// bridge that makes a later runtime unable to invent geometry outside its map.
#include "MaterialSurfacePatchMap.hpp"
#include "FsiDomainRuntime.hpp"

#include <limits>

namespace iga {

class MaterialSurfacePatchKinematics {
public:
	struct Result {
		MaterialSurfaceKinematics target;
		std::string composition_identity_sha256;
	};

	static Result ComposeTarget(const MaterialSurfacePatchMap& map,
		const MaterialSurfaceKinematics& committed_full, const SurfaceKinematics& trial_patch,
		const FsiTrialContext& context)
	{
		map.FullReference().Validate(); committed_full.Validate(); ValidateFsiTrialContext(context);
		ValidateSurfaceKinematics(trial_patch, map.Layout());
		if (!(trial_patch.interface == map.Interface().id)
			|| trial_patch.stamp.time_s != context.EndTime() || trial_patch.stamp.step != context.step
			|| trial_patch.stamp.coupling_iteration != context.coupling_iteration)
			throw std::runtime_error("patch kinematics does not match its exact FSI trial context");
		if (committed_full.MaterialIdentitySha256() != map.FullReference().MaterialIdentitySha256()
			|| committed_full.TopologyIdentitySha256() != map.FullReference().TopologyIdentitySha256()
			|| committed_full.EvaluatedTimeS() != context.start_time_s
			|| committed_full.StepEndS() != context.start_time_s)
			throw std::runtime_error("committed full material state does not match map or context start");

		auto positions = map.FullReference().ReferenceMaterialVerticesM();
		std::vector<std::array<double, 3>> velocities(positions.size(), {{0.0, 0.0, 0.0}});
		for (std::size_t local = 0; local < map.Layout().owned_global_node_ids.size(); ++local) {
			const auto global = map.Layout().owned_global_node_ids[local];
			const auto source = map.SourceVertexForGlobalNode(global);
			const bool clamped = std::binary_search(map.ConfiguredClampedGlobalNodeIds().begin(),
				map.ConfiguredClampedGlobalNodeIds().end(), global);
			if (clamped && (!ExactlyZero(trial_patch.displacement_m[local]) || !ExactlyZero(trial_patch.velocity_m_per_s[local])))
				throw std::runtime_error("clamped material patch seam must publish exact zero displacement and velocity");
			for (std::size_t axis = 0; axis < 3; ++axis) {
				positions[source][axis] = map.FullReference().ReferenceMaterialVerticesM()[source][axis]
					+ trial_patch.displacement_m[local][axis];
				velocities[source][axis] = trial_patch.velocity_m_per_s[local][axis];
			}
		}
		for (std::size_t source = 0; source < positions.size(); ++source)
			for (std::size_t axis = 0; axis < 3; ++axis)
				if (!std::isfinite(positions[source][axis])) throw std::runtime_error("patch target position is nonfinite");
		ValidateBackwardEuler(committed_full.SourceVerticesM(), positions, velocities, context.dt_s);

		std::vector<RawSurfaceTriangle> source_triangles;
		source_triangles.reserve(committed_full.SourceTriangles().size());
		for (const auto& source : committed_full.SourceTriangles()) {
			RawSurfaceTriangle triangle;
			for (std::size_t corner = 0; corner < 3; ++corner) triangle.indices[corner] = source.source_vertex_indices[corner];
			triangle.boundary_id = source.boundary_id;
			source_triangles.push_back(triangle);
		}
		auto target = MaterialSurfaceKinematics::CreateFromSourceTopology(
			map.FullReference().ReferenceMaterialVerticesM(), std::move(positions), std::move(velocities),
			std::move(source_triangles), context.EndTime(), context.start_time_s, context.EndTime());
		if (target.MaterialIdentitySha256() != committed_full.MaterialIdentitySha256()
			|| target.TopologyIdentitySha256() != committed_full.TopologyIdentitySha256())
			throw std::runtime_error("patch composition changed immutable full material/topology identity");
		Result result{std::move(target), {}};
		result.composition_identity_sha256 = HashComposition(map, committed_full, trial_patch, context, result.target);
		return result;
	}

private:
	static bool ExactlyZero(const std::array<double, 3>& value)
	{ return value[0] == 0.0 && value[1] == 0.0 && value[2] == 0.0; }
	static void ValidateBackwardEuler(const std::vector<std::array<double, 3>>& committed,
		const std::vector<std::array<double, 3>>& target, const std::vector<std::array<double, 3>>& velocity, double dt)
	{
		if (committed.size() != target.size() || target.size() != velocity.size())
			throw std::runtime_error("patch backward-Euler fields have incompatible source sizes");
		double minimum[3] = {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
		double maximum[3] = {-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
		for (const auto& vertices : {&committed, &target}) for (const auto& vertex : *vertices)
			for (std::size_t axis = 0; axis < 3; ++axis) { minimum[axis] = std::min(minimum[axis], vertex[axis]); maximum[axis] = std::max(maximum[axis], vertex[axis]); }
		double geometry_length = 0.0;
		for (std::size_t axis = 0; axis < 3; ++axis) geometry_length = std::hypot(geometry_length, maximum[axis]-minimum[axis]);
		for (std::size_t vertex = 0; vertex < target.size(); ++vertex)
			for (std::size_t axis = 0; axis < 3; ++axis) {
				const double expected = dt*velocity[vertex][axis];
				const double actual = target[vertex][axis]-committed[vertex][axis];
				if (!std::isfinite(expected) || !std::isfinite(actual))
					throw std::runtime_error("patch target violates backward-Euler displacement/velocity consistency");
				if (actual == expected) continue;
				const double scale = std::max({geometry_length, std::fabs(target[vertex][axis]), std::fabs(committed[vertex][axis]), std::fabs(expected), std::fabs(actual)});
				const double tolerance = 2048.0*std::numeric_limits<double>::epsilon()*scale;
				if (std::fabs(actual-expected) > tolerance)
					throw std::runtime_error("patch target violates backward-Euler displacement/velocity consistency");
			}
	}
	static std::string HashComposition(const MaterialSurfacePatchMap& map,
		const MaterialSurfaceKinematics& committed, const SurfaceKinematics& patch,
		const FsiTrialContext& context, const MaterialSurfaceKinematics& target)
	{
		Sha256 hash;
		distributed_surface_detail::AppendString(hash, "MaterialSurfacePatchKinematics/composition/v1");
		distributed_surface_detail::AppendString(hash, map.IdentitySha256());
		distributed_surface_detail::AppendString(hash, committed.ContentIdentitySha256());
		distributed_surface_detail::AppendString(hash, BuildSurfaceKinematicsIdentitySha256(patch, map.Layout()));
		hash.AppendLittleEndian64(context.step); hash.AppendNormalizedDouble(context.start_time_s);
		hash.AppendNormalizedDouble(context.dt_s); hash.AppendLittleEndian64(context.coupling_iteration);
		distributed_surface_detail::AppendString(hash, target.ContentIdentitySha256());
		return hash.Hex();
	}
};

} // namespace iga

#endif
