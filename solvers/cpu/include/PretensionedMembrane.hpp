#ifndef IGA_PRETENSIONED_MEMBRANE_HPP
#define IGA_PRETENSIONED_MEMBRANE_HPP

// Bounded, rank-local P1 membrane kernel for the first FSI structural slice.
// It deliberately owns no runtime or MPI communication.  The scalar unknown
// is normal displacement on the immutable reference material surface.
#include "DistributedSurfaceInterface.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

class PretensionedMembraneFsiRuntime;

struct PretensionedMembraneMaterial {
	// rho_A [kg/m^2], c_A [kg/(m^2 s)], T0 [N/m], k_A [N/m^3].
	double areal_mass_kg_per_m2 = 0.0;
	double damping_kg_per_m2_s = 0.0;
	double pretension_n_per_m = 0.0;
	double foundation_n_per_m3 = 0.0;
};

struct PretensionedMembraneOptions {
	// Dense direct solve is intentionally limited to a benchmark-sized,
	// single-partition surface; this is not an MPI-scalable structural solver.
	std::size_t maximum_nodes = 4096;
};

struct PretensionedMembraneState {
	// Scalars aligned with DistributedSurfaceLayout::owned_global_node_ids.
	std::vector<double> displacement_m;
	std::vector<double> velocity_m_per_s;
};

struct PretensionedMembraneTrialContext {
	std::uint64_t step = 0;
	double start_time_s = 0.0;
	double dt_s = 0.0;
	std::uint64_t coupling_iteration = 0;
	// The fluid stamp is an exact expected input contract, not a template.
	SurfaceFieldStamp expected_traction_stamp;

	double EndTime() const { return start_time_s+dt_s; }
};

enum class PretensionedMembraneTrialPhase : std::uint8_t { Solved, Prepared };

struct PretensionedMembraneTrial {
	PretensionedMembraneState state;
	SurfaceKinematics kinematics;
	std::string base_committed_state_identity_sha256;
	std::string prepared_committed_state_identity_sha256;
	std::string trial_identity_sha256;
	bool prepared = false;

	private:
	friend class PretensionedMembrane;
	// These are opaque lifecycle capabilities issued by one membrane instance.
	// They are never caller-settable and are checked against owner-held state.
	std::uint64_t owner_token = 0;
	std::uint64_t generation = 0;
	PretensionedMembraneTrialPhase phase = PretensionedMembraneTrialPhase::Solved;
};

namespace pretensioned_membrane_detail {

inline std::array<double, 3> Subtract(const std::array<double, 3>& left,
	const std::array<double, 3>& right)
{ return {{left[0]-right[0], left[1]-right[1], left[2]-right[2]}}; }

inline std::array<double, 3> Cross(const std::array<double, 3>& left,
	const std::array<double, 3>& right)
{
	return {{left[1]*right[2]-left[2]*right[1], left[2]*right[0]-left[0]*right[2],
		left[0]*right[1]-left[1]*right[0]}};
}

inline double Dot(const std::array<double, 3>& left, const std::array<double, 3>& right)
{ return left[0]*right[0]+left[1]*right[1]+left[2]*right[2]; }

inline std::array<std::uint64_t, 3> CanonicalCyclicTriangle(std::array<std::uint64_t, 3> triangle)
{
	const auto first = triangle;
	const std::array<std::uint64_t, 3> second{{triangle[1], triangle[2], triangle[0]}};
	const std::array<std::uint64_t, 3> third{{triangle[2], triangle[0], triangle[1]}};
	return std::min(first, std::min(second, third));
}

inline void AppendState(Sha256& hash, const PretensionedMembraneState& state)
{
	hash.AppendLittleEndian64(static_cast<std::uint64_t>(state.displacement_m.size()));
	for (const double value : state.displacement_m) hash.AppendNormalizedDouble(value);
	hash.AppendLittleEndian64(static_cast<std::uint64_t>(state.velocity_m_per_s.size()));
	for (const double value : state.velocity_m_per_s) hash.AppendNormalizedDouble(value);
}

inline void RequireFinite(double value, const char* name)
{
	if (!std::isfinite(value)) throw std::runtime_error(std::string(name)+" must be finite");
}

} // namespace pretensioned_membrane_detail

class PretensionedMembrane {
public:
	PretensionedMembrane(const PretensionedMembrane&) = delete;
	PretensionedMembrane& operator=(const PretensionedMembrane&) = delete;
	PretensionedMembrane(PretensionedMembrane&&) = delete;
	PretensionedMembrane& operator=(PretensionedMembrane&&) = delete;

	PretensionedMembrane(DistributedSurfaceLayout layout,
		DistributedSurfaceInterface structure_interface,
		SurfaceInterfaceRef traction_interface,
		PretensionedMembraneMaterial material,
		std::vector<std::uint64_t> clamped_global_node_ids,
		PretensionedMembraneOptions options = PretensionedMembraneOptions{},
		PretensionedMembraneState initial_state = PretensionedMembraneState{})
		: layout_(std::move(layout)), structure_interface_(std::move(structure_interface)),
		  traction_interface_(std::move(traction_interface)), material_(material),
		  clamped_global_node_ids_(std::move(clamped_global_node_ids)), options_(options),
		  owner_token_(AllocateOwnerToken())
	{
		ValidateDistributedSurfaceLayout(layout_);
		if (layout_.partition_count != 1 || layout_.partition_rank != 0)
			throw std::runtime_error("pretensioned membrane supports one partition only");
		if (layout_.owned_global_node_ids.size() > options_.maximum_nodes)
			throw std::runtime_error("pretensioned membrane exceeds its dense node limit");
		if (options_.maximum_nodes == 0) throw std::runtime_error("pretensioned membrane node limit must be positive");
		ValidateDistributedSurfaceInterface(structure_interface_);
		ValidateSurfaceInterfaceRef(traction_interface_);
		if (structure_interface_.reference_mesh_identity_sha256 != layout_.reference_mesh_identity_sha256)
			throw std::runtime_error("pretensioned membrane structure interface and layout differ");
		const std::vector<SurfaceFieldQuantity> provides{SurfaceFieldQuantity::Displacement, SurfaceFieldQuantity::Velocity};
		const std::vector<SurfaceFieldQuantity> requires{SurfaceFieldQuantity::TractionOnStructure};
		if (structure_interface_.provides != provides || structure_interface_.requires != requires)
			throw std::runtime_error("pretensioned membrane structure interface must provide displacement/velocity and require traction");
		ValidateMaterial();
		if (!std::is_sorted(clamped_global_node_ids_.begin(), clamped_global_node_ids_.end())
			|| std::adjacent_find(clamped_global_node_ids_.begin(), clamped_global_node_ids_.end()) != clamped_global_node_ids_.end())
			throw std::runtime_error("pretensioned membrane clamped IDs must be canonical sorted and unique");
		BuildIndexAndGeometry();
		for (const auto id : clamped_global_node_ids_) {
			auto found = std::lower_bound(layout_.owned_global_node_ids.begin(), layout_.owned_global_node_ids.end(), id);
			if (found == layout_.owned_global_node_ids.end() || *found != id)
				throw std::runtime_error("pretensioned membrane clamps an unknown global node");
			clamped_[static_cast<std::size_t>(found-layout_.owned_global_node_ids.begin())] = true;
		}
		AssembleReferenceOperators();
		RequireStaticCoercivity();
		if (initial_state.displacement_m.empty() && initial_state.velocity_m_per_s.empty()) {
			initial_state.displacement_m.assign(NodeCount(), 0.0);
			initial_state.velocity_m_per_s.assign(NodeCount(), 0.0);
		}
		ValidateState(initial_state, "initial state");
		committed_state_ = std::move(initial_state);
		committed_state_identity_sha256_ = BuildStateIdentity(committed_state_);
		model_identity_sha256_ = BuildModelIdentity();
	}

	const DistributedSurfaceLayout& Layout() const noexcept { return layout_; }
	const std::vector<std::array<double, 3>>& ReferenceVertexNormals() const noexcept { return normals_; }
	const std::vector<bool>& ClampedNodeMask() const noexcept { return clamped_; }
	const PretensionedMembraneState& CommittedState() const noexcept { return committed_state_; }
	const std::string& ModelIdentitySha256() const noexcept { return model_identity_sha256_; }
	const std::string& CommittedStateIdentitySha256() const noexcept { return committed_state_identity_sha256_; }

	PretensionedMembraneTrial SolveTrial(const PretensionedMembraneTrialContext& context,
		const SurfaceTraction& traction) const
	{
		if (has_active_trial_)
			throw std::runtime_error("pretensioned membrane requires RejectTrial or AbortTrial before a new solve");
		ValidateContext(context);
		ValidateSurfaceTraction(traction, layout_);
		if (!(traction.interface == traction_interface_))
			throw std::runtime_error("pretensioned membrane traction belongs to the wrong interface");
		if (BuildSurfaceFieldStampIdentitySha256(traction.stamp, layout_)
			!= BuildSurfaceFieldStampIdentitySha256(context.expected_traction_stamp, layout_))
			throw std::runtime_error("pretensioned membrane traction stamp is stale or unexpected");
		const std::size_t count = NodeCount();
		std::vector<double> load(count, 0.0);
		for (std::size_t node = 0; node < count; ++node)
			load[node] = pretensioned_membrane_detail::Dot(traction.consistent_nodal_force_n[node], normals_[node]);
		const double inverse_dt = 1.0/context.dt_s;
		std::vector<double> system(count*count, 0.0), rhs(count, 0.0);
		for (std::size_t row = 0; row < count; ++row) {
			for (std::size_t column = 0; column < count; ++column) {
				const double mass = mass_[row*count+column];
				const double stiffness = stiffness_[row*count+column];
				system[row*count+column] = (material_.areal_mass_kg_per_m2*inverse_dt
					+material_.damping_kg_per_m2_s)*mass+context.dt_s*stiffness;
				rhs[row] += (material_.areal_mass_kg_per_m2*inverse_dt*mass)*committed_state_.velocity_m_per_s[column]
					-stiffness*committed_state_.displacement_m[column];
			}
			rhs[row] += load[row];
		}
		std::vector<double> velocity = SolveWithClamps(system, rhs);
		std::vector<double> displacement(count, 0.0);
		for (std::size_t node = 0; node < count; ++node)
			displacement[node] = clamped_[node] ? 0.0
				: committed_state_.displacement_m[node]+context.dt_s*velocity[node];
		PretensionedMembraneTrial result;
		result.state.displacement_m = std::move(displacement);
		result.state.velocity_m_per_s = std::move(velocity);
		result.base_committed_state_identity_sha256 = committed_state_identity_sha256_;
		if (next_generation_ == std::numeric_limits<std::uint64_t>::max())
			throw std::runtime_error("pretensioned membrane trial generation is exhausted");
		result.owner_token = owner_token_;
		result.generation = next_generation_+1;
		result.phase = PretensionedMembraneTrialPhase::Solved;
		result.prepared = false;
		result.kinematics.interface = structure_interface_.id;
		result.kinematics.displacement_m.resize(count);
		result.kinematics.velocity_m_per_s.resize(count);
		for (std::size_t node = 0; node < count; ++node) {
			for (int component = 0; component < 3; ++component) {
				result.kinematics.displacement_m[node][component] = result.state.displacement_m[node]*normals_[node][component];
				result.kinematics.velocity_m_per_s[node][component] = result.state.velocity_m_per_s[node]*normals_[node][component];
			}
		}
		result.kinematics.stamp = MakeOutputStamp(context, traction, result.state);
		result.trial_identity_sha256 = BuildTrialIdentity(result);
		next_generation_ = result.generation;
		active_generation_ = result.generation;
		active_phase_ = PretensionedMembraneTrialPhase::Solved;
		active_trial_identity_sha256_ = result.trial_identity_sha256;
		has_active_trial_ = true;
		return result;
	}

	PretensionedMembraneTrial PrepareTrial(PretensionedMembraneTrial trial) const
	{
		ValidateTrial(trial, false);
		// Complete all allocation/validation before changing the active
		// capability.  The prepared identity is what FinalizePreparedTrialNoexcept
		// exchanges with the committed identity.
		trial.prepared_committed_state_identity_sha256 = BuildStateIdentity(trial.state);
		trial.prepared = true;
		trial.phase = PretensionedMembraneTrialPhase::Prepared;
		trial.trial_identity_sha256 = BuildTrialIdentity(trial);
		// Copy while the active capability still describes the solved trial.  A
		// later swap and phase assignment are noexcept, so allocation failure
		// cannot strand this owner in Prepared with a stale identity.
		std::string prepared_active_identity = trial.trial_identity_sha256;
		active_trial_identity_sha256_.swap(prepared_active_identity);
		active_phase_ = PretensionedMembraneTrialPhase::Prepared;
		return trial;
	}

	void FinalizeTrial(const PretensionedMembraneTrial& trial)
	{
		RequireFinalizeAllowed(trial);
		// Preserve the standalone value-taking API.  The runtime uses the rvalue
		// variant below to avoid copying numerical state at commit.
		PretensionedMembraneTrial candidate = trial;
		FinalizePreparedTrialNoexcept(std::move(candidate));
	}

	void RequireFinalizeAllowed(const PretensionedMembraneTrial& trial) const
	{
		ValidateTrial(trial, true);
		if (trial.prepared_committed_state_identity_sha256.empty()
			|| trial.prepared_committed_state_identity_sha256 != BuildStateIdentity(trial.state))
			throw std::runtime_error("pretensioned membrane prepared trial state identity is invalid");
	}

	void RejectTrial(const PretensionedMembraneTrial& trial) const
	{
		ValidateTrial(trial, false);
		ClearActiveTrial();
	}
	void AbortTrial() const noexcept { ClearActiveTrial(); }

private:
	friend class PretensionedMembraneFsiRuntime;
	void FinalizePreparedTrialNoexcept(PretensionedMembraneTrial&& trial) noexcept
	{
		// RequireFinalizeAllowed() must have succeeded immediately beforehand.
		committed_state_.displacement_m.swap(trial.state.displacement_m);
		committed_state_.velocity_m_per_s.swap(trial.state.velocity_m_per_s);
		committed_state_identity_sha256_.swap(trial.prepared_committed_state_identity_sha256);
		ClearActiveTrial();
	}
	static std::uint64_t AllocateOwnerToken()
	{
		static std::atomic<std::uint64_t> next(0);
		const std::uint64_t token = next.fetch_add(1, std::memory_order_relaxed)+1;
		if (token == 0) throw std::runtime_error("pretensioned membrane owner token is exhausted");
		return token;
	}

	void ClearActiveTrial() const noexcept
	{
		has_active_trial_ = false;
		active_generation_ = 0;
		active_trial_identity_sha256_.clear();
	}

	std::size_t NodeCount() const noexcept { return layout_.owned_global_node_ids.size(); }

	void ValidateMaterial() const
	{
		pretensioned_membrane_detail::RequireFinite(material_.areal_mass_kg_per_m2, "membrane areal mass");
		pretensioned_membrane_detail::RequireFinite(material_.damping_kg_per_m2_s, "membrane damping");
		pretensioned_membrane_detail::RequireFinite(material_.pretension_n_per_m, "membrane pretension");
		pretensioned_membrane_detail::RequireFinite(material_.foundation_n_per_m3, "membrane foundation");
		if (!(material_.areal_mass_kg_per_m2 > 0.0) || material_.damping_kg_per_m2_s < 0.0
			|| material_.pretension_n_per_m < 0.0 || material_.foundation_n_per_m3 < 0.0)
			throw std::runtime_error("pretensioned membrane physical coefficients have invalid signs");
	}

	void BuildIndexAndGeometry()
	{
		const std::size_t count = NodeCount();
		positions_.resize(count);
		for (std::size_t node = 0; node < count; ++node) {
			const auto id = layout_.owned_global_node_ids[node];
			auto found = std::lower_bound(layout_.reference_positions.begin(), layout_.reference_positions.end(), id,
				[](const SurfaceReferencePosition& value, std::uint64_t needle) { return value.global_node_id < needle; });
			if (found == layout_.reference_positions.end() || found->global_node_id != id)
				throw std::runtime_error("pretensioned membrane layout lost an owned position");
			positions_[node] = found->position_m;
		}
		std::vector<std::array<std::uint64_t, 3>> triangles = layout_.reference_triangles;
		for (auto& triangle : triangles) triangle = pretensioned_membrane_detail::CanonicalCyclicTriangle(triangle);
		std::sort(triangles.begin(), triangles.end());
		for (std::size_t index = 1; index < triangles.size(); ++index)
			if (triangles[index] == triangles[index-1]) throw std::runtime_error("pretensioned membrane has duplicate triangles");
		struct Edge { std::uint64_t low, high; int direction; };
		std::vector<Edge> edges;
		edges.reserve(3*triangles.size());
		for (const auto& triangle : triangles) {
			for (int edge = 0; edge < 3; ++edge) {
				const auto first = triangle[static_cast<std::size_t>(edge)];
				const auto second = triangle[static_cast<std::size_t>((edge+1)%3)];
				edges.push_back({std::min(first, second), std::max(first, second), first < second ? 1 : -1});
			}
		}
		std::sort(edges.begin(), edges.end(), [](const Edge& left, const Edge& right) {
			return left.low != right.low ? left.low < right.low : left.high < right.high;
		});
		for (std::size_t first = 0; first < edges.size();) {
			std::size_t last = first+1;
			while (last < edges.size() && edges[last].low == edges[first].low && edges[last].high == edges[first].high) ++last;
			if (last-first > 2 || (last-first == 2 && edges[first].direction == edges[first+1].direction))
				throw std::runtime_error("pretensioned membrane triangle winding is nonorientable or inconsistent");
			first = last;
		}
		triangles_ = std::move(triangles);
		normals_.assign(count, {{0.0, 0.0, 0.0}});
		for (const auto& triangle : triangles_) {
			const auto indices = TriangleIndices(triangle);
			const auto cross = pretensioned_membrane_detail::Cross(
				pretensioned_membrane_detail::Subtract(positions_[indices[1]], positions_[indices[0]]),
				pretensioned_membrane_detail::Subtract(positions_[indices[2]], positions_[indices[0]]));
			const double squared = pretensioned_membrane_detail::Dot(cross, cross);
			if (!std::isfinite(squared) || !(squared > 0.0)) throw std::runtime_error("pretensioned membrane has degenerate triangle");
			for (const auto index : indices)
				for (int component = 0; component < 3; ++component) normals_[index][component] += cross[component];
		}
		for (auto& normal : normals_) {
			const double length = std::sqrt(pretensioned_membrane_detail::Dot(normal, normal));
			if (!std::isfinite(length) || !(length > 0.0)) throw std::runtime_error("pretensioned membrane has a zero reference vertex normal");
			for (double& component : normal) component /= length;
		}
		clamped_.assign(count, false);
	}

	std::array<std::size_t, 3> TriangleIndices(const std::array<std::uint64_t, 3>& triangle) const
	{
		std::array<std::size_t, 3> result{};
		for (int local = 0; local < 3; ++local) {
			auto found = std::lower_bound(layout_.owned_global_node_ids.begin(), layout_.owned_global_node_ids.end(), triangle[local]);
			if (found == layout_.owned_global_node_ids.end() || *found != triangle[local])
				throw std::runtime_error("pretensioned membrane triangle is not locally owned");
			result[local] = static_cast<std::size_t>(found-layout_.owned_global_node_ids.begin());
		}
		return result;
	}

	void AssembleReferenceOperators()
	{
		const std::size_t count = NodeCount();
		mass_.assign(count*count, 0.0); stiffness_.assign(count*count, 0.0);
		for (const auto& triangle : triangles_) {
			const auto indices = TriangleIndices(triangle);
			const auto first_edge = pretensioned_membrane_detail::Subtract(positions_[indices[1]], positions_[indices[0]]);
			const auto second_edge = pretensioned_membrane_detail::Subtract(positions_[indices[2]], positions_[indices[0]]);
			const auto cross = pretensioned_membrane_detail::Cross(first_edge, second_edge);
			const double twice_area = std::sqrt(pretensioned_membrane_detail::Dot(cross, cross));
			const double area = 0.5*twice_area;
			std::array<double, 3> unit_normal = cross;
			for (double& component : unit_normal) component /= twice_area;
			const std::array<double, 3> gradient1 = Scale(pretensioned_membrane_detail::Cross(second_edge, unit_normal), 1.0/twice_area);
			const std::array<double, 3> gradient2 = Scale(pretensioned_membrane_detail::Cross(unit_normal, first_edge), 1.0/twice_area);
			const std::array<double, 3> gradient0 = Scale(Add(gradient1, gradient2), -1.0);
			const std::array<std::array<double, 3>, 3> gradients{{gradient0, gradient1, gradient2}};
			for (int row = 0; row < 3; ++row) for (int column = 0; column < 3; ++column) {
				mass_[indices[row]*count+indices[column]] += area*(row == column ? 2.0 : 1.0)/12.0;
				stiffness_[indices[row]*count+indices[column]] += material_.pretension_n_per_m
					*area*pretensioned_membrane_detail::Dot(gradients[row], gradients[column])
					+material_.foundation_n_per_m3*area*(row == column ? 2.0 : 1.0)/12.0;
			}
		}
	}

	static std::array<double, 3> Add(const std::array<double, 3>& left, const std::array<double, 3>& right)
	{ return {{left[0]+right[0], left[1]+right[1], left[2]+right[2]}}; }
	static std::array<double, 3> Scale(std::array<double, 3> value, double scale)
	{ for (double& component : value) component *= scale; return value; }

	void RequireStaticCoercivity() const
	{
		std::vector<std::size_t> free_nodes;
		for (std::size_t node = 0; node < NodeCount(); ++node) if (!clamped_[node]) free_nodes.push_back(node);
		if (free_nodes.empty()) return;
		std::vector<double> matrix(free_nodes.size()*free_nodes.size());
		for (std::size_t row = 0; row < free_nodes.size(); ++row)
			for (std::size_t column = 0; column < free_nodes.size(); ++column)
				matrix[row*free_nodes.size()+column] = stiffness_[free_nodes[row]*NodeCount()+free_nodes[column]];
		if (!FactorSpdCholesky(matrix))
			throw std::runtime_error("pretensioned membrane static operator is not coercive after Dirichlet elimination");
	}

	static bool FactorSpdCholesky(std::vector<double>& matrix)
	{
		const std::size_t count = static_cast<std::size_t>(std::sqrt(static_cast<double>(matrix.size())));
		if (count == 0 || count*count != matrix.size()) return false;
		double scale = 0.0;
		for (const double value : matrix) {
			if (!std::isfinite(value)) return false;
			scale = std::max(scale, std::abs(value));
		}
		if (!(scale > 0.0)) return false;
		for (double& value : matrix) value /= scale;
		const double symmetry_tolerance = 512.0*std::numeric_limits<double>::epsilon();
		for (std::size_t row = 0; row < count; ++row) for (std::size_t column = 0; column < row; ++column) {
			const double pair_scale = std::max(std::abs(matrix[row*count+column]), std::abs(matrix[column*count+row]));
			if (std::abs(matrix[row*count+column]-matrix[column*count+row]) > symmetry_tolerance*pair_scale) return false;
		}
		const double pivot_tolerance = 128.0*std::numeric_limits<double>::epsilon();
		for (std::size_t row = 0; row < count; ++row) {
			for (std::size_t column = 0; column <= row; ++column) {
				double value = matrix[row*count+column];
				for (std::size_t inner = 0; inner < column; ++inner) value -= matrix[row*count+inner]*matrix[column*count+inner];
				if (row == column) {
					if (!std::isfinite(value) || !(value > pivot_tolerance)) return false;
					matrix[row*count+column] = std::sqrt(value);
				} else {
					if (!std::isfinite(value)) return false;
					matrix[row*count+column] = value/matrix[column*count+column];
				}
			}
		}
		return true;
	}

	std::vector<double> SolveWithClamps(const std::vector<double>& matrix, const std::vector<double>& rhs) const
	{
		std::vector<std::size_t> free_nodes;
		for (std::size_t node = 0; node < NodeCount(); ++node) if (!clamped_[node]) free_nodes.push_back(node);
		std::vector<double> result(NodeCount(), 0.0);
		if (free_nodes.empty()) return result;
		const std::size_t count = free_nodes.size();
		std::vector<double> scaled_matrix(count*count), value(count);
		for (std::size_t row = 0; row < count; ++row) {
			value[row] = rhs[free_nodes[row]];
			for (std::size_t column = 0; column < count; ++column)
				scaled_matrix[row*count+column] = matrix[free_nodes[row]*NodeCount()+free_nodes[column]];
		}
		double scale = 0.0;
		for (const double coefficient : scaled_matrix) {
			if (!std::isfinite(coefficient)) throw std::runtime_error("pretensioned membrane dynamic system is non-finite");
			scale = std::max(scale, std::abs(coefficient));
		}
		if (!(scale > 0.0)) throw std::runtime_error("pretensioned membrane dynamic system is singular");
		for (double& coefficient : scaled_matrix) coefficient /= scale;
		for (double& entry : value) {
			entry /= scale;
			if (!std::isfinite(entry)) throw std::runtime_error("pretensioned membrane dynamic right-hand side is non-finite");
		}
		const std::vector<double> scaled_rhs = value;
		std::vector<double> factor = scaled_matrix;
		if (!FactorSpdCholesky(factor)) throw std::runtime_error("pretensioned membrane dynamic system is not SPD");
		for (std::size_t row = 0; row < count; ++row) {
			for (std::size_t column = 0; column < row; ++column) value[row] -= factor[row*count+column]*value[column];
			value[row] /= factor[row*count+row];
		}
		for (std::size_t reverse = count; reverse-- > 0;) {
			for (std::size_t column = reverse+1; column < count; ++column) value[reverse] -= factor[column*count+reverse]*value[column];
			value[reverse] /= factor[reverse*count+reverse];
		}
		double matrix_norm = 0.0, rhs_norm = 0.0, solution_norm = 0.0, residual_norm = 0.0;
		for (std::size_t row = 0; row < count; ++row) {
			double row_sum = 0.0, residual = -scaled_rhs[row];
			for (std::size_t column = 0; column < count; ++column) {
				row_sum += std::abs(scaled_matrix[row*count+column]);
				residual += scaled_matrix[row*count+column]*value[column];
			}
			matrix_norm = std::max(matrix_norm, row_sum);
			rhs_norm = std::max(rhs_norm, std::abs(scaled_rhs[row]));
			solution_norm = std::max(solution_norm, std::abs(value[row]));
			residual_norm = std::max(residual_norm, std::abs(residual));
			if (!std::isfinite(value[row])) throw std::runtime_error("pretensioned membrane dynamic solution is non-finite");
			result[free_nodes[row]] = value[row];
		}
		const double residual_scale = matrix_norm*solution_norm+rhs_norm;
		if (!std::isfinite(matrix_norm) || !std::isfinite(rhs_norm) || !std::isfinite(solution_norm)
			|| !std::isfinite(residual_scale))
			throw std::runtime_error("pretensioned membrane dynamic residual scale is non-finite");
		const double residual_tolerance = 4096.0*std::numeric_limits<double>::epsilon()
			*std::max(1.0, static_cast<double>(count))*residual_scale;
		if (!std::isfinite(residual_norm) || residual_norm > residual_tolerance)
			throw std::runtime_error("pretensioned membrane dynamic solve residual is too large");
		return result;
	}

	void ValidateContext(const PretensionedMembraneTrialContext& context) const
	{
		if (!std::isfinite(context.start_time_s) || !std::isfinite(context.dt_s) || !(context.dt_s > 0.0)
			|| !std::isfinite(context.EndTime())) throw std::runtime_error("pretensioned membrane trial requires finite time and positive dt");
		ValidateSurfaceFieldStamp(context.expected_traction_stamp, layout_);
		if (context.expected_traction_stamp.time_s != context.EndTime() || context.expected_traction_stamp.step != context.step
			|| context.expected_traction_stamp.coupling_iteration != context.coupling_iteration)
			throw std::runtime_error("pretensioned membrane expected traction stamp does not match its trial context");
	}

	void ValidateState(const PretensionedMembraneState& state, const char* name) const
	{
		if (state.displacement_m.size() != NodeCount() || state.velocity_m_per_s.size() != NodeCount())
			throw std::runtime_error(std::string("pretensioned membrane ")+name+" does not align with owned nodes");
		for (std::size_t node = 0; node < NodeCount(); ++node) {
			pretensioned_membrane_detail::RequireFinite(state.displacement_m[node], name);
			pretensioned_membrane_detail::RequireFinite(state.velocity_m_per_s[node], name);
			if (clamped_[node] && (state.displacement_m[node] != 0.0 || state.velocity_m_per_s[node] != 0.0))
				throw std::runtime_error("pretensioned membrane clamped state must be exactly zero");
		}
	}

	std::string BuildStateIdentity(const PretensionedMembraneState& state) const
	{
		ValidateState(state, "state"); Sha256 hash;
		distributed_surface_detail::AppendString(hash, "PretensionedMembrane/committed-state/v1");
		distributed_surface_detail::AppendString(hash, layout_.layout_identity_sha256);
		pretensioned_membrane_detail::AppendState(hash, state); return hash.Hex();
	}

	std::string BuildModelIdentity() const
	{
		Sha256 hash; distributed_surface_detail::AppendString(hash, "PretensionedMembrane/model/v1");
		distributed_surface_detail::AppendString(hash, BuildDistributedSurfaceInterfaceIdentitySha256(structure_interface_));
		distributed_surface_detail::AppendString(hash, traction_interface_.domain_id);
		distributed_surface_detail::AppendString(hash, traction_interface_.subsystem_id);
		distributed_surface_detail::AppendString(hash, traction_interface_.interface_id);
		distributed_surface_detail::AppendString(hash, layout_.layout_identity_sha256);
		hash.AppendNormalizedDouble(material_.areal_mass_kg_per_m2); hash.AppendNormalizedDouble(material_.damping_kg_per_m2_s);
		hash.AppendNormalizedDouble(material_.pretension_n_per_m); hash.AppendNormalizedDouble(material_.foundation_n_per_m3);
		for (const auto id : clamped_global_node_ids_) hash.AppendLittleEndian64(id);
		return hash.Hex();
	}

	SurfaceFieldStamp MakeOutputStamp(const PretensionedMembraneTrialContext& context,
		const SurfaceTraction& traction, const PretensionedMembraneState& trial_state) const
	{
		Sha256 hash; distributed_surface_detail::AppendString(hash, "PretensionedMembrane/producer-state/v1");
		distributed_surface_detail::AppendString(hash, model_identity_sha256_);
		distributed_surface_detail::AppendString(hash, committed_state_identity_sha256_);
		// A kinematics publication is a statement about this computed trial, not
		// merely about the inputs that happened to determine it.
		distributed_surface_detail::AppendString(hash, BuildStateIdentity(trial_state));
		distributed_surface_detail::AppendString(hash, BuildSurfaceTractionIdentitySha256(traction, layout_));
		hash.AppendNormalizedDouble(context.start_time_s); hash.AppendNormalizedDouble(context.dt_s);
		hash.AppendLittleEndian64(context.step); hash.AppendLittleEndian64(context.coupling_iteration);
		SurfaceFieldStamp output; output.time_s = context.EndTime(); output.step = context.step;
		output.coupling_iteration = context.coupling_iteration; output.reference_mesh_identity_sha256 = layout_.reference_mesh_identity_sha256;
		output.layout_identity_sha256 = layout_.layout_identity_sha256; output.partition_identity_sha256 = BuildDistributedSurfacePartitionIdentitySha256(layout_);
		output.producer_state_identity_sha256 = hash.Hex(); return output;
	}

	std::string BuildTrialIdentity(const PretensionedMembraneTrial& trial) const
	{
		ValidateState(trial.state, "trial state"); ValidateSurfaceKinematics(trial.kinematics, layout_);
		Sha256 hash; distributed_surface_detail::AppendString(hash, "PretensionedMembrane/trial/v1");
		distributed_surface_detail::AppendString(hash, model_identity_sha256_);
		distributed_surface_detail::AppendString(hash, trial.base_committed_state_identity_sha256);
		distributed_surface_detail::AppendString(hash, trial.prepared_committed_state_identity_sha256);
		pretensioned_membrane_detail::AppendState(hash, trial.state);
		distributed_surface_detail::AppendString(hash, BuildSurfaceKinematicsIdentitySha256(trial.kinematics, layout_));
		hash.AppendLittleEndian64(trial.owner_token);
		hash.AppendLittleEndian64(trial.generation);
		hash.AppendLittleEndian32(static_cast<std::uint32_t>(trial.phase));
		hash.AppendLittleEndian32(trial.prepared ? 1U : 0U); return hash.Hex();
	}

	void ValidateTrial(const PretensionedMembraneTrial& trial, bool require_prepared) const
	{
		if (!has_active_trial_ || trial.owner_token != owner_token_ || trial.generation != active_generation_
			|| trial.phase != active_phase_ || trial.trial_identity_sha256 != active_trial_identity_sha256_)
			throw std::runtime_error("pretensioned membrane trial is not the current lifecycle capability");
		if (trial.base_committed_state_identity_sha256 != committed_state_identity_sha256_)
			throw std::runtime_error("pretensioned membrane trial is stale against committed state");
		if (trial.prepared != (trial.phase == PretensionedMembraneTrialPhase::Prepared))
			throw std::runtime_error("pretensioned membrane trial phase is inconsistent");
		if (require_prepared != trial.prepared)
			throw std::runtime_error(require_prepared ? "pretensioned membrane finalize requires a prepared trial"
				: "pretensioned membrane reject requires a solved trial");
		if (!(trial.kinematics.interface == structure_interface_.id) || trial.trial_identity_sha256 != BuildTrialIdentity(trial))
			throw std::runtime_error("pretensioned membrane trial identity is invalid");
	}

	DistributedSurfaceLayout layout_;
	DistributedSurfaceInterface structure_interface_;
	SurfaceInterfaceRef traction_interface_;
	PretensionedMembraneMaterial material_;
	std::vector<std::uint64_t> clamped_global_node_ids_;
	PretensionedMembraneOptions options_;
	std::vector<std::array<double, 3>> positions_, normals_;
	std::vector<std::array<std::uint64_t, 3>> triangles_;
	std::vector<bool> clamped_;
	std::vector<double> mass_, stiffness_;
	PretensionedMembraneState committed_state_;
	std::string model_identity_sha256_, committed_state_identity_sha256_;
	const std::uint64_t owner_token_;
	mutable std::uint64_t next_generation_ = 0, active_generation_ = 0;
	mutable bool has_active_trial_ = false;
	mutable PretensionedMembraneTrialPhase active_phase_ = PretensionedMembraneTrialPhase::Solved;
	mutable std::string active_trial_identity_sha256_;
};

} // namespace iga

#endif
