#include "PretensionedMembraneFsiRuntime.hpp"

#include <array>
#include <cassert>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

template <class Function>
void Reject(Function&& function)
{
	bool rejected = false;
	try { function(); }
	catch (const std::runtime_error&) { rejected = true; }
	assert(rejected);
}

std::string Hash(const char* text)
{
	iga::Sha256 hash;
	hash.Append(text, std::char_traits<char>::length(text));
	return hash.Hex();
}

iga::DistributedSurfaceLayout Layout()
{
	iga::DistributedSurfaceLayout value;
	value.reference_mesh_identity_sha256 = Hash("membrane-fsi-runtime-mesh");
	value.global_node_count = 3; value.partition_count = 1; value.partition_rank = 0;
	value.owned_global_node_ids = {0, 1, 2};
	value.reference_positions = {{0, {{0.0, 0.0, 0.0}}}, {1, {{1.0, 0.0, 0.0}}},
		{2, {{0.0, 1.0, 0.0}}}};
	value.reference_triangles = {{{0, 1, 2}}};
	value.owned_reference_lumped_areas_m2 = {1.0/6.0, 1.0/6.0, 1.0/6.0};
	value.layout_identity_sha256 = iga::BuildDistributedSurfaceLayoutIdentitySha256(value);
	return value;
}

iga::DistributedSurfaceInterface Structure(const std::string& mesh)
{
	iga::DistributedSurfaceInterface value;
	value.id = {"structure", "membrane", "patch"}; value.subsystem_id = "membrane";
	value.boundary_labels = {0}; value.reference_mesh_identity_sha256 = mesh;
	value.provides = {iga::SurfaceFieldQuantity::Displacement, iga::SurfaceFieldQuantity::Velocity};
	value.requires = {iga::SurfaceFieldQuantity::TractionOnStructure};
	return value;
}

iga::DistributedSurfaceInterface Fluid(const std::string& mesh)
{
	iga::DistributedSurfaceInterface value;
	value.id = {"fluid", "immersed", "patch"}; value.subsystem_id = "immersed";
	value.boundary_labels = {0}; value.reference_mesh_identity_sha256 = mesh;
	value.provides = {iga::SurfaceFieldQuantity::TractionOnStructure};
	value.requires = {iga::SurfaceFieldQuantity::Displacement, iga::SurfaceFieldQuantity::Velocity};
	return value;
}

iga::SurfaceFieldStamp Stamp(const iga::DistributedSurfaceLayout& layout, double time,
	std::uint64_t step, std::uint64_t iteration, const char* state)
{
	iga::SurfaceFieldStamp value;
	value.time_s = time; value.step = step; value.coupling_iteration = iteration;
	value.reference_mesh_identity_sha256 = layout.reference_mesh_identity_sha256;
	value.layout_identity_sha256 = layout.layout_identity_sha256;
	value.partition_identity_sha256 = iga::BuildDistributedSurfacePartitionIdentitySha256(layout);
	value.producer_state_identity_sha256 = Hash(state);
	return value;
}

iga::SurfaceFieldStampEnvelope Envelope(const iga::DistributedSurfaceLayout& layout, double time,
	std::uint64_t step, std::uint64_t iteration)
{
	iga::SurfaceFieldStampEnvelope value;
	value.time_s = time; value.step = step; value.coupling_iteration = iteration;
	value.reference_mesh_identity_sha256 = layout.reference_mesh_identity_sha256;
	value.layout_identity_sha256 = layout.layout_identity_sha256;
	value.partition_identity_sha256 = iga::BuildDistributedSurfacePartitionIdentitySha256(layout);
	return value;
}

iga::SurfaceTraction Traction(const iga::DistributedSurfaceLayout& layout,
	const iga::SurfaceInterfaceRef& fluid, const iga::SurfaceFieldStamp& stamp, double load = 2.0)
{
	iga::SurfaceTraction value;
	value.interface = fluid; value.stamp = stamp;
	value.traction_on_structure_pa.assign(layout.owned_global_node_ids.size(), {{0.0, 0.0, 0.0}});
	value.consistent_nodal_force_n.assign(layout.owned_global_node_ids.size(), {{0.0, 0.0, 0.0}});
	value.consistent_nodal_force_n[2][2] = load;
	value.projection_identity_sha256 = Hash("runtime-projection");
	return value;
}

iga::PretensionedMembraneFsiRuntime Runtime(const iga::DistributedSurfaceLayout& layout,
	const iga::DistributedSurfaceInterface& structure, const iga::DistributedSurfaceInterface& fluid,
	const iga::FsiCouplingEdge& edge)
{
	return iga::PretensionedMembraneFsiRuntime("structure", "membrane", edge, structure, fluid,
		layout, layout, {1.0, .1, 2.0, 1.0}, {0, 1});
}

void Solve(iga::PretensionedMembraneFsiRuntime& runtime, const iga::SurfaceTraction& traction,
	const iga::SurfaceFieldStampEnvelope& output_envelope, std::uint64_t iteration)
{
	runtime.BeginCouplingIteration(iteration, traction.stamp, output_envelope);
	runtime.SetSurfaceTraction("patch", traction);
	runtime.SolveMembraneTrial();
}

} // namespace

static_assert(!std::is_copy_constructible<iga::PretensionedMembraneFsiRuntime>::value,
	"runtime must have one numerical/lifecycle owner");
static_assert(!std::is_move_constructible<iga::PretensionedMembraneFsiRuntime>::value,
	"runtime must not move issued membrane capabilities");

int main()
{
	const auto layout = Layout(); const auto structure = Structure(layout.reference_mesh_identity_sha256);
	const auto fluid = Fluid(layout.reference_mesh_identity_sha256);
	const iga::FsiCouplingEdge edge("membrane-wall", fluid.id, structure.id,
		iga::FsiCouplingLaw::FluidStructureTractionKinematics);
	auto runtime = Runtime(layout, structure, fluid, edge);
	const iga::DomainStepContext step{4, .1, .05};
	auto wrong_structure = structure; wrong_structure.id.interface_id = "other";
	Reject([&] { Runtime(layout, wrong_structure, fluid, edge); });
	auto wrong_edge = edge; wrong_edge.structure = fluid.id;
	Reject([&] { Runtime(layout, structure, fluid, wrong_edge); });

	Reject([&] { runtime.GetSurfaceKinematics("patch"); });
	Reject([&] { runtime.BeginCouplingIteration(0, {}, {}); });
	runtime.BeginMacroStep(step);
	const auto input0 = Stamp(layout, .1+.05, 4, 0, "fluid-0");
	const auto traction0 = Traction(layout, fluid.id, input0);
	const auto envelope0 = Envelope(layout, .1+.05, 4, 0);
	auto invalid_envelope = envelope0; invalid_envelope.coupling_iteration = 9;
	Reject([&] { runtime.BeginCouplingIteration(0, input0, invalid_envelope); });
	Reject([&] { runtime.SetSurfaceTraction("other", traction0); });
	runtime.BeginCouplingIteration(0, input0, envelope0);
	Reject([&] { runtime.BeginCouplingIteration(0, input0, envelope0); });
	auto wrong_interface = traction0; wrong_interface.interface = structure.id;
	Reject([&] { runtime.SetSurfaceTraction("patch", wrong_interface); });
	auto wrong_stamp = traction0; wrong_stamp.stamp.coupling_iteration = 9;
	Reject([&] { runtime.SetSurfaceTraction("patch", wrong_stamp); });
	runtime.SetSurfaceTraction("patch", traction0);
	Reject([&] { runtime.SetSurfaceTraction("patch", traction0); });
	runtime.SolveTrial();
	auto trial_snapshot = runtime.GetSurfaceKinematics("patch");
	const auto output0 = trial_snapshot.stamp;
	iga::ValidateSurfaceFieldStampMatchesEnvelope(output0, envelope0, layout);
	trial_snapshot.displacement_m[2][2] = 123.0;
	assert(runtime.GetSurfaceKinematics("patch").displacement_m[2][2] != 123.0);
	Reject([&] { runtime.GetCommittedSurfaceKinematics("patch", output0); });
	runtime.RejectCouplingIteration();
	Reject([&] { runtime.GetSurfaceKinematics("patch"); });
	// A rejected iteration accepts the exact same fluid publication again; no
	// coordinator needs to predict or duplicate the membrane numerical solve.
	Solve(runtime, traction0, envelope0, 0);
	assert(runtime.GetSurfaceKinematics("patch").stamp.producer_state_identity_sha256
		== output0.producer_state_identity_sha256);
	runtime.RejectCouplingIteration();

	const auto input1 = Stamp(layout, .1+.05, 4, 1, "fluid-1");
	const auto traction1 = Traction(layout, fluid.id, input1, 3.0);
	const auto envelope1 = Envelope(layout, .1+.05, 4, 1);
	Solve(runtime, traction1, envelope1, 1);
	const auto output1 = runtime.GetSurfaceKinematics("patch").stamp;
	runtime.PrepareCommitStep();
	Reject([&] { runtime.PrepareCommitStep(); });
	Reject([&] { runtime.GetCommittedSurfaceKinematics("patch", output1); });
	runtime.RejectCouplingIteration();
	Reject([&] { runtime.GetSurfaceKinematics("patch"); });

	const auto input2 = Stamp(layout, .1+.05, 4, 2, "fluid-2");
	const auto traction2 = Traction(layout, fluid.id, input2, 4.0);
	const auto envelope2 = Envelope(layout, .1+.05, 4, 2);
	Solve(runtime, traction2, envelope2, 2);
	const auto output2 = runtime.GetSurfaceKinematics("patch").stamp;
	Reject([&] { runtime.FinalizeCommitStep(); });
	runtime.PrepareCommit();
	const auto before_finalize = runtime.GetSurfaceKinematics("patch");
	runtime.FinalizeCommit();
	Reject([&] { runtime.GetSurfaceKinematics("patch"); });
	const auto committed = runtime.GetCommittedSurfaceKinematics("patch", output2);
	assert(committed.stamp.producer_state_identity_sha256 == output2.producer_state_identity_sha256);
	assert(committed.displacement_m == before_finalize.displacement_m);
	Reject([&] { runtime.GetCommittedSurfaceKinematics("patch", output1); });

	runtime.BeginStep(step);
	const auto bad_input = Stamp(layout, .1+.05, 4, 3, "fluid-3");
	const auto bad_traction = Traction(layout, fluid.id, bad_input);
	auto incoherent_envelope = Envelope(layout, .1+.05, 4, 3);
	incoherent_envelope.reference_mesh_identity_sha256 = Hash("foreign-mesh");
	Reject([&] { runtime.BeginIteration(3, bad_input, incoherent_envelope); });
	const auto envelope3 = Envelope(layout, .1+.05, 4, 3);
	runtime.BeginIteration(3, bad_input, envelope3);
	runtime.SetSurfaceTraction("patch", bad_traction);
	runtime.SolveMembraneTrial();
	iga::ValidateSurfaceFieldStampMatchesEnvelope(runtime.GetSurfaceKinematics("patch").stamp,
		envelope3, layout);
	runtime.RollbackTrial();
	const auto input4 = Stamp(layout, .1+.05, 4, 4, "fluid-4");
	const auto traction4 = Traction(layout, fluid.id, input4);
	const auto envelope4 = Envelope(layout, .1+.05, 4, 4);
	const auto state_snapshot = runtime.CommittedStateSnapshot();
	assert(!state_snapshot.displacement_m.empty());
	Solve(runtime, traction4, envelope4, 4);
	runtime.PrepareCommitStep();
	assert(runtime.GetCommittedSurfaceKinematics("patch", output2).stamp.producer_state_identity_sha256
		== output2.producer_state_identity_sha256);
	runtime.AbortStep();
	Reject([&] { runtime.GetSurfaceKinematics("patch"); });
	assert(runtime.GetCommittedSurfaceKinematics("patch", output2).stamp.step == 4);

	return 0;
}
