#include "DistributedSurfaceInterface.hpp"
#include "DynamicWeightedAitkenRelaxation.hpp"

#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

template <class Function>
void RequireRejected(Function&& function)
{
	bool rejected = false;
	try { function(); }
	catch (const std::runtime_error&) { rejected = true; }
	assert(rejected);
}

const std::string kMesh(64, 'a');
const std::string kState(64, 'b');
const std::string kProjection(64, 'c');

iga::DistributedSurfaceLayout MakeLayout()
{
	iga::DistributedSurfaceLayout layout;
	layout.reference_mesh_identity_sha256 = kMesh;
	layout.global_node_count = 3;
	layout.partition_count = 1;
	layout.partition_rank = 0;
	layout.owned_global_node_ids = {10, 20, 30};
	layout.reference_positions = {{10, {{0.0, 0.0, 0.0}}}, {20, {{1.0, 0.0, 0.0}}}, {30, {{0.0, 1.0, 0.0}}}};
	layout.reference_triangles = {{{10, 20, 30}}};
	layout.owned_reference_lumped_areas_m2 = {1.0, 2.0, 3.0};
	layout.layout_identity_sha256 = iga::BuildDistributedSurfaceLayoutIdentitySha256(layout);
	return layout;
}

iga::DistributedSurfaceLayout MakeFourNodeLayout(std::uint64_t rank,
	const std::vector<std::uint64_t>& ids, const std::vector<double>& areas)
{
	iga::DistributedSurfaceLayout layout;
	layout.reference_mesh_identity_sha256 = kMesh;
	layout.global_node_count = 4;
	layout.partition_count = 2;
	layout.partition_rank = rank;
	layout.owned_global_node_ids = ids;
	layout.reference_positions = {{10, {{0.0, 0.0, 0.0}}}, {20, {{1.0, 0.0, 0.0}}},
		{30, {{0.0, 1.0, 0.0}}}, {40, {{1.0, 1.0, 0.0}}}};
	layout.reference_triangles = {{{10, 20, 30}}, {{20, 40, 30}}};
	layout.owned_reference_lumped_areas_m2 = areas;
	layout.layout_identity_sha256 = iga::BuildDistributedSurfaceLayoutIdentitySha256(layout);
	return layout;
}

iga::SurfaceFieldStamp MakeStamp(const iga::DistributedSurfaceLayout& layout)
{
	iga::SurfaceFieldStamp stamp;
	stamp.time_s = 0.125;
	stamp.step = std::numeric_limits<std::uint64_t>::max();
	stamp.coupling_iteration = std::numeric_limits<std::uint64_t>::max()-1;
	stamp.reference_mesh_identity_sha256 = layout.reference_mesh_identity_sha256;
	stamp.layout_identity_sha256 = layout.layout_identity_sha256;
	stamp.partition_identity_sha256 = iga::BuildDistributedSurfacePartitionIdentitySha256(layout);
	stamp.producer_state_identity_sha256 = kState;
	return stamp;
}

iga::SurfaceKinematics MakeKinematics(const iga::DistributedSurfaceLayout& layout)
{
	iga::SurfaceKinematics value;
	value.interface = {"fluid", "fluid_backend", "membrane"};
	value.stamp = MakeStamp(layout);
	value.displacement_m = {{{0.0, 0.0, 0.01}}, {{0.0, 0.0, 0.02}}, {{0.0, 0.0, 0.03}}};
	value.velocity_m_per_s = {{{0.0, 0.0, 1.0}}, {{0.0, 0.0, 2.0}}, {{0.0, 0.0, 3.0}}};
	return value;
}

iga::SurfaceTraction MakeTraction(const iga::DistributedSurfaceLayout& layout)
{
	iga::SurfaceTraction value;
	value.interface = {"fluid", "fluid_backend", "membrane"};
	value.stamp = MakeStamp(layout);
	value.traction_on_structure_pa = {{{0.0, 0.0, -10.0}}, {{0.0, 0.0, -20.0}}, {{0.0, 0.0, -30.0}}};
	value.consistent_nodal_force_n = {{{0.0, 0.0, -1.0}}, {{0.0, 0.0, -2.0}}, {{0.0, 0.0, -3.0}}};
	value.projection_identity_sha256 = kProjection;
	return value;
}

iga::DistributedSurfaceInterface MakeInterface()
{
	iga::DistributedSurfaceInterface surface;
	surface.id = {"fluid", "fluid_backend", "membrane"};
	surface.subsystem_id = "fluid_backend";
	surface.boundary_labels = {0, 7};
	surface.reference_mesh_identity_sha256 = kMesh;
	surface.provides = {iga::SurfaceFieldQuantity::TractionOnStructure};
	surface.requires = {iga::SurfaceFieldQuantity::Displacement, iga::SurfaceFieldQuantity::Velocity};
	return surface;
}

} // namespace

int main()
{
	const auto layout = MakeLayout();
	iga::ValidateDistributedSurfaceLayout(layout);
	assert(iga::BuildDistributedSurfaceLayoutIdentitySha256(layout) == layout.layout_identity_sha256);
	const auto kinematics = MakeKinematics(layout);
	const auto traction = MakeTraction(layout);
	iga::ValidateSurfaceKinematics(kinematics, layout);
	iga::ValidateSurfaceTraction(traction, layout);
	assert(iga::BuildSurfaceKinematicsIdentitySha256(kinematics, layout)
		== iga::BuildSurfaceKinematicsIdentitySha256(kinematics, layout));
	assert(iga::BuildSurfaceTractionIdentitySha256(traction, layout)
		== iga::BuildSurfaceTractionIdentitySha256(traction, layout));
	assert(iga::MakeCartesianSurfaceLumpedAreaWeights(layout)
		== std::vector<double>({1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0}));

	// A surface communicator may include ranks with no publication authority.
	const auto empty_partition = MakeFourNodeLayout(0, {}, {});
	const auto other_empty_partition = MakeFourNodeLayout(1, {}, {});
	iga::ValidateDistributedSurfaceLayout(empty_partition);
	assert(empty_partition.layout_identity_sha256 == other_empty_partition.layout_identity_sha256);
	assert(iga::BuildDistributedSurfacePartitionIdentitySha256(empty_partition)
		!= iga::BuildDistributedSurfacePartitionIdentitySha256(other_empty_partition));
	assert(iga::MakeCartesianSurfaceLumpedAreaWeights(empty_partition).empty());
	auto empty_kinematics = MakeKinematics(empty_partition);
	empty_kinematics.displacement_m.clear(); empty_kinematics.velocity_m_per_s.clear();
	iga::ValidateSurfaceKinematics(empty_kinematics, empty_partition);
	auto empty_traction = MakeTraction(empty_partition);
	empty_traction.traction_on_structure_pa.clear(); empty_traction.consistent_nodal_force_n.clear();
	iga::ValidateSurfaceTraction(empty_traction, empty_partition);
	RequireRejected([&] { iga::ValidateSurfaceKinematics(empty_kinematics, other_empty_partition); });
	empty_traction.consistent_nodal_force_n.push_back({{1., 0., 0.}});
	RequireRejected([&] { iga::ValidateSurfaceTraction(empty_traction, empty_partition); });
	auto empty_single = layout;
	empty_single.owned_global_node_ids.clear(); empty_single.owned_reference_lumped_areas_m2.clear();
	RequireRejected([&] { iga::BuildDistributedSurfaceLayoutIdentitySha256(empty_single); });

	// Same global topology and local vector shape/weights, but a different
	// rank-owned ID slice. Equal-sized slices are not interchangeable.
	auto first_partition = MakeFourNodeLayout(0, {10, 20}, {1.0, 2.0});
	auto different_partition = MakeFourNodeLayout(1, {30, 40}, {1.0, 2.0});
	iga::ValidateDistributedSurfaceLayout(first_partition);
	iga::ValidateDistributedSurfaceLayout(different_partition);
	assert(first_partition.layout_identity_sha256 == different_partition.layout_identity_sha256);
	assert(first_partition.owned_global_node_ids.size() == different_partition.owned_global_node_ids.size());
	assert(first_partition.owned_reference_lumped_areas_m2 == different_partition.owned_reference_lumped_areas_m2);
	const auto first_stamp = MakeStamp(first_partition);
	const auto second_stamp = MakeStamp(different_partition);
	assert(first_stamp.partition_identity_sha256 != second_stamp.partition_identity_sha256);
	assert(iga::BuildSurfaceFieldStampIdentitySha256(first_stamp, first_partition)
		!= iga::BuildSurfaceFieldStampIdentitySha256(second_stamp, different_partition));
	RequireRejected([&] { iga::ValidateSurfaceFieldStamp(first_stamp, different_partition); });
	RequireRejected([&] { iga::BuildSurfaceFieldStampIdentitySha256(first_stamp, different_partition); });
	iga::SurfaceKinematics equal_shape_field;
	equal_shape_field.interface = {"fluid", "fluid_backend", "membrane"};
	equal_shape_field.stamp = first_stamp;
	equal_shape_field.displacement_m = {{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}};
	equal_shape_field.velocity_m_per_s = {{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}};
	RequireRejected([&] { iga::ValidateSurfaceKinematics(equal_shape_field, different_partition); });

	const auto surface = MakeInterface();
	iga::ValidateDistributedSurfaceInterface(surface);
	assert(surface.boundary_labels.front() == 0);
	const auto surface_identity = iga::BuildDistributedSurfaceInterfaceIdentitySha256(surface);
	assert(surface_identity == iga::BuildDistributedSurfaceInterfaceIdentitySha256(surface));
	auto interface_mutation = surface;
	interface_mutation.boundary_labels = {0, 8};
	assert(surface_identity != iga::BuildDistributedSurfaceInterfaceIdentitySha256(interface_mutation));
	auto subsystem_mutation = surface;
	subsystem_mutation.id.subsystem_id = "other_backend";
	subsystem_mutation.subsystem_id = "other_backend";
	assert(surface_identity != iga::BuildDistributedSurfaceInterfaceIdentitySha256(subsystem_mutation));
	// These had identical concatenated quantity bytes before collection counts.
	auto split_one = surface;
	split_one.provides = {iga::SurfaceFieldQuantity::Displacement};
	split_one.requires = {iga::SurfaceFieldQuantity::Velocity, iga::SurfaceFieldQuantity::TractionOnStructure};
	auto split_two = surface;
	split_two.provides = {iga::SurfaceFieldQuantity::Displacement, iga::SurfaceFieldQuantity::Velocity};
	split_two.requires = {iga::SurfaceFieldQuantity::TractionOnStructure};
	assert(iga::BuildDistributedSurfaceInterfaceIdentitySha256(split_one)
		!= iga::BuildDistributedSurfaceInterfaceIdentitySha256(split_two));

	RequireRejected([] { iga::ValidateSurfaceInterfaceRef({"", "fluid_backend", "membrane"}); });
	RequireRejected([] { iga::ValidateSurfaceInterfaceRef({"fluid", "", "membrane"}); });
	RequireRejected([&] { auto invalid = layout; invalid.reference_mesh_identity_sha256 = "not-a-hash"; iga::ValidateDistributedSurfaceLayout(invalid); });
	RequireRejected([&] { auto invalid = layout; invalid.reference_positions[0].position_m[0] = NAN; iga::ValidateDistributedSurfaceLayout(invalid); });
	RequireRejected([&] { auto invalid = layout; invalid.owned_global_node_ids = {10, 10, 30}; iga::ValidateDistributedSurfaceLayout(invalid); });
	RequireRejected([&] { auto invalid = layout; invalid.owned_global_node_ids = {20, 10, 30}; iga::ValidateDistributedSurfaceLayout(invalid); });
	RequireRejected([&] { auto invalid = layout; invalid.partition_count = 0; iga::ValidateDistributedSurfaceLayout(invalid); });
	RequireRejected([&] { auto invalid = layout; invalid.partition_rank = 1; iga::ValidateDistributedSurfaceLayout(invalid); });
	RequireRejected([&] { auto invalid = layout; invalid.owned_reference_lumped_areas_m2[0] = 0.0; iga::ValidateDistributedSurfaceLayout(invalid); });
	RequireRejected([&] { auto invalid = layout; invalid.reference_triangles[0][2] = 99; iga::ValidateDistributedSurfaceLayout(invalid); });
	RequireRejected([&] { auto invalid = kinematics; invalid.displacement_m.pop_back(); iga::ValidateSurfaceKinematics(invalid, layout); });
	RequireRejected([&] { auto invalid = MakeStamp(layout); invalid.partition_identity_sha256 = kState; iga::ValidateSurfaceFieldStamp(invalid, layout); });
	RequireRejected([&] { auto invalid = MakeStamp(layout); invalid.time_s = NAN; iga::ValidateSurfaceFieldStamp(invalid); });
	RequireRejected([&] { auto invalid = traction; invalid.projection_identity_sha256 = "bad"; iga::ValidateSurfaceTraction(invalid, layout); });
	RequireRejected([&] { auto invalid = surface; invalid.boundary_labels.clear(); iga::ValidateDistributedSurfaceInterface(invalid); });
	RequireRejected([&] { auto invalid = surface; invalid.boundary_labels = {-1}; iga::ValidateDistributedSurfaceInterface(invalid); });
	RequireRejected([&] { auto invalid = surface; invalid.boundary_labels = {7, 0}; iga::ValidateDistributedSurfaceInterface(invalid); });
	RequireRejected([&] { auto invalid = surface; invalid.provides = {iga::SurfaceFieldQuantity::Velocity, iga::SurfaceFieldQuantity::Displacement}; iga::ValidateDistributedSurfaceInterface(invalid); });
	RequireRejected([&] { auto invalid = surface; invalid.provides.push_back(iga::SurfaceFieldQuantity::TractionOnStructure); iga::ValidateDistributedSurfaceInterface(invalid); });
	RequireRejected([&] { auto invalid = surface; invalid.requires.push_back(iga::SurfaceFieldQuantity::TractionOnStructure); iga::ValidateDistributedSurfaceInterface(invalid); });
	RequireRejected([&] { auto invalid = surface; invalid.provides.clear(); iga::ValidateDistributedSurfaceInterface(invalid); });
	RequireRejected([&] { auto invalid = surface; invalid.requires.clear(); iga::ValidateDistributedSurfaceInterface(invalid); });

	// Manual two-rank reduction with unequal local areas, compared to one rank.
	const auto rank_zero = MakeFourNodeLayout(0, {10, 20}, {1.0, 4.0});
	const auto rank_one = MakeFourNodeLayout(1, {30, 40}, {2.0, 3.0});
	const std::string rank_zero_id = iga::BuildDistributedSurfacePartitionIdentitySha256(rank_zero);
	const std::string rank_one_id = iga::BuildDistributedSurfacePartitionIdentitySha256(rank_one);
	const double total_area = 10.0;
	iga::DynamicWeightedAitkenRelaxation distributed_zero(rank_zero_id, {1.0, 4.0}, total_area);
	iga::DynamicWeightedAitkenRelaxation distributed_one(rank_one_id, {2.0, 3.0}, total_area);
	const std::vector<double> previous_zero{2.0, -1.0};
	const std::vector<double> previous_one{3.0, 1.0};
	const std::vector<double> current_zero{1.0, 1.0};
	const std::vector<double> current_one{2.0, 4.0};
	assert(distributed_zero.ControlStateIdentitySha256() == distributed_one.ControlStateIdentitySha256());
	const auto initial_control_identity = distributed_zero.ControlStateIdentitySha256();
	const auto initial_zero = distributed_zero.ProposeWithGlobalReduction({0.0, 0.0}, previous_zero,
		1.0, 3.0, 0.0, 0.0, initial_control_identity, rank_zero_id);
	// A locally pending proposal and a proposal-ready rank advertise distinct
	// collective phases. Once every rank has the same global proposal facts,
	// their pending identities match even though their local vectors differ.
	assert(distributed_zero.ControlStateIdentitySha256() != distributed_one.ControlStateIdentitySha256());
	const auto initial_one = distributed_one.ProposeWithGlobalReduction({0.0, 0.0}, previous_one,
		1.0, 3.0, 0.0, 0.0, initial_control_identity, rank_one_id);
	assert(distributed_zero.ControlStateIdentitySha256() == distributed_one.ControlStateIdentitySha256());
	distributed_zero.AcceptApplied(initial_zero, previous_zero, initial_zero.relaxation_factor, rank_zero_id);
	distributed_one.AcceptApplied(initial_one, previous_one, initial_one.relaxation_factor, rank_one_id);
	assert(distributed_zero.ControlStateIdentitySha256() == distributed_one.ControlStateIdentitySha256());
	const auto local_zero = distributed_zero.LocalContributions(current_zero, rank_zero_id);
	const auto local_one = distributed_one.LocalContributions(current_one, rank_one_id);
	assert(std::abs(local_zero.numerator+10.0) < 1.0e-15 && std::abs(local_zero.denominator-17.0) < 1.0e-15);
	assert(std::abs(local_one.numerator-3.0) < 1.0e-15 && std::abs(local_one.denominator-29.0) < 1.0e-15);
	const double numerator = local_zero.numerator+local_one.numerator;
	const double denominator = local_zero.denominator+local_one.denominator;
	const double global_scale = 4.0;
	const auto reduced_control_identity = distributed_zero.ControlStateIdentitySha256();
	const auto reduced_zero = distributed_zero.ProposeWithGlobalReduction({0.0, 0.0}, current_zero,
		1.0, global_scale, numerator, denominator, reduced_control_identity, rank_zero_id);
	const auto reduced_one = distributed_one.ProposeWithGlobalReduction({0.0, 0.0}, current_one,
		1.0, global_scale, numerator, denominator, reduced_control_identity, rank_one_id);
	assert(reduced_zero.status == iga::AitkenRelaxationStatus::Dynamic);
	assert(std::abs(reduced_zero.relaxation_factor-7.0/92.0) < 1.0e-15);
	assert(reduced_zero.relaxation_factor == reduced_one.relaxation_factor);
	assert(distributed_zero.ControlStateIdentitySha256() == distributed_one.ControlStateIdentitySha256());
	RequireRejected([&] { distributed_zero.AcceptApplied(reduced_one, current_zero,
		reduced_one.relaxation_factor, rank_zero_id); });
	distributed_zero.AcceptApplied(reduced_zero, current_zero, reduced_zero.relaxation_factor, rank_zero_id);
	distributed_one.AcceptApplied(reduced_one, current_one, reduced_one.relaxation_factor, rank_one_id);
	assert(distributed_zero.ControlStateIdentitySha256() == distributed_one.ControlStateIdentitySha256());
	auto one_rank = MakeFourNodeLayout(0, {10, 20, 30, 40}, {1.0, 4.0, 2.0, 3.0});
	one_rank.partition_count = 1;
	one_rank.layout_identity_sha256 = iga::BuildDistributedSurfaceLayoutIdentitySha256(one_rank);
	const std::string one_rank_id = iga::BuildDistributedSurfacePartitionIdentitySha256(one_rank);
	iga::DynamicWeightedAitkenRelaxation single(one_rank_id, {1.0, 4.0, 2.0, 3.0}, total_area);
	const auto single_initial = single.Propose({0.0, 0.0, 0.0, 0.0}, {2.0, -1.0, 3.0, 1.0}, 1.0, one_rank_id);
	single.AcceptApplied(single_initial, {2.0, -1.0, 3.0, 1.0}, single_initial.relaxation_factor, one_rank_id);
	const auto one_proposal = single.Propose({0.0, 0.0, 0.0, 0.0}, {1.0, 1.0, 2.0, 4.0}, 1.0, one_rank_id);
	assert(one_proposal.relaxation_factor == reduced_zero.relaxation_factor);
	RequireRejected([&] { distributed_zero.Propose({0.0, 0.0}, current_zero, 1.0, rank_zero_id); });
	RequireRejected([&] { distributed_zero.ProposeWithGlobalReduction({0.0, 0.0}, current_zero, 1.0, 1.0, numerator, denominator, reduced_control_identity, rank_zero_id); });
	RequireRejected([&] { distributed_zero.LocalContributions(current_zero, rank_one_id); });

	// This two-component example proves that reference-area weights affect omega.
	iga::DynamicWeightedAitkenRelaxation weighted(one_rank_id, {1.0, 3.0}, 4.0);
	iga::DynamicWeightedAitkenRelaxation uniform(one_rank_id, {1.0, 1.0}, 2.0);
	const auto weighted_initial = weighted.Propose({0.0, 0.0}, {2.0, -1.0}, 1.0, one_rank_id);
	const auto uniform_initial = uniform.Propose({0.0, 0.0}, {2.0, -1.0}, 1.0, one_rank_id);
	weighted.AcceptApplied(weighted_initial, {2.0, -1.0}, weighted_initial.relaxation_factor, one_rank_id);
	uniform.AcceptApplied(uniform_initial, {2.0, -1.0}, uniform_initial.relaxation_factor, one_rank_id);
	const auto weighted_result = weighted.Propose({0.0, 0.0}, {1.0, 1.0}, 1.0, one_rank_id);
	const auto uniform_result = uniform.Propose({0.0, 0.0}, {1.0, 1.0}, 1.0, one_rank_id);
	assert(weighted_result.status == iga::AitkenRelaxationStatus::Dynamic);
	assert(uniform_result.status == iga::AitkenRelaxationStatus::Dynamic);
	assert(std::abs(weighted_result.relaxation_factor-uniform_result.relaxation_factor) > 1.0e-12);

	// A proposal is an owned transaction: an arbitrary scalar, a stale proposal,
	// or a reset-invalidated proposal cannot advance the residual history.
	iga::DynamicWeightedAitkenRelaxation ownership_guard(one_rank_id, {1.0, 1.0}, 2.0);
	const auto guard_initial = ownership_guard.Propose({0.0, 0.0}, {2.0, -1.0}, 1.0, one_rank_id);
	auto changed_next = guard_initial;
	changed_next.next[0] += 1.0;
	RequireRejected([&] { ownership_guard.AcceptApplied(changed_next, {2.0, -1.0}, changed_next.relaxation_factor, one_rank_id); });
	auto changed_status = guard_initial;
	changed_status.status = iga::AitkenRelaxationStatus::Dynamic;
	RequireRejected([&] { ownership_guard.AcceptApplied(changed_status, {2.0, -1.0}, changed_status.relaxation_factor, one_rank_id); });
	auto changed_diagnostic = guard_initial;
	changed_diagnostic.unclamped_relaxation_factor += 0.125;
	RequireRejected([&] { ownership_guard.AcceptApplied(changed_diagnostic, {2.0, -1.0}, changed_diagnostic.relaxation_factor, one_rank_id); });
	auto changed_global = guard_initial;
	changed_global.global_numerator += 1.0;
	RequireRejected([&] { ownership_guard.AcceptApplied(changed_global, {2.0, -1.0}, changed_global.relaxation_factor, one_rank_id); });
	RequireRejected([&] { ownership_guard.AcceptApplied(guard_initial, {2.0, -1.0}, 0.4, one_rank_id); });
	RequireRejected([&] { ownership_guard.AcceptApplied(guard_initial, {1.0, -1.0}, guard_initial.relaxation_factor, one_rank_id); });
	ownership_guard.AcceptApplied(guard_initial, {2.0, -1.0}, guard_initial.relaxation_factor, one_rank_id);
	RequireRejected([&] { ownership_guard.AcceptApplied(guard_initial, {2.0, -1.0}, guard_initial.relaxation_factor, one_rank_id); });
	const auto stale_proposal = ownership_guard.Propose({0.0, 0.0}, {1.0, 1.0}, 1.0, one_rank_id);
	ownership_guard.Reset();
	RequireRejected([&] { ownership_guard.AcceptApplied(stale_proposal, {1.0, 1.0}, stale_proposal.relaxation_factor, one_rank_id); });

	// The caller allgathers or broadcasts this identity before global reduction.
	// A reset (phase/generation) mismatch fails closed before any reduction is used.
	iga::DynamicWeightedAitkenRelaxation phase_zero(one_rank_id, {1.0, 1.0}, 2.0);
	iga::DynamicWeightedAitkenRelaxation phase_one(one_rank_id, {1.0, 1.0}, 2.0);
	const auto phase_identity = phase_zero.ControlStateIdentitySha256();
	phase_one.Reset();
	RequireRejected([&] { phase_one.ProposeWithGlobalReduction({0.0, 0.0}, {1.0, 1.0},
		1.0, 1.0, 0.0, 0.0, phase_identity, one_rank_id); });

	// Equal generations and accepted counts still reject if the accepted omega
	// differs; rank-local residual histories themselves need not match.
	iga::DynamicWeightedAitkenRelaxation prior_state_zero(one_rank_id, {1.0, 1.0}, 2.0);
	iga::DynamicWeightedAitkenRelaxation prior_state_one(one_rank_id, {1.0, 1.0}, 2.0);
	const auto previous_zero_initial = prior_state_zero.Propose({0.0, 0.0}, {2.0, -1.0}, 1.0, one_rank_id);
	const auto previous_one_initial = prior_state_one.Propose({0.0, 0.0}, {3.0, 1.0}, 1.0, one_rank_id);
	prior_state_zero.AcceptApplied(previous_zero_initial, {2.0, -1.0}, previous_zero_initial.relaxation_factor, one_rank_id);
	prior_state_one.AcceptApplied(previous_one_initial, {3.0, 1.0}, previous_one_initial.relaxation_factor, one_rank_id);
	const auto previous_zero_dynamic = prior_state_zero.Propose({0.0, 0.0}, {1.0, 1.0}, 1.0, one_rank_id);
	const auto previous_one_dynamic = prior_state_one.Propose({0.0, 0.0}, {0.0, 1.0}, 1.0, one_rank_id);
	prior_state_zero.AcceptApplied(previous_zero_dynamic, {1.0, 1.0}, previous_zero_dynamic.relaxation_factor, one_rank_id);
	prior_state_one.AcceptApplied(previous_one_dynamic, {0.0, 1.0}, previous_one_dynamic.relaxation_factor, one_rank_id);
	assert(previous_zero_dynamic.relaxation_factor != previous_one_dynamic.relaxation_factor);
	RequireRejected([&] { prior_state_one.ProposeWithGlobalReduction({0.0, 0.0}, {1.0, 1.0},
		1.0, 1.0, 0.0, 0.0, prior_state_zero.ControlStateIdentitySha256(), one_rank_id); });
}
