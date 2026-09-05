#include "CoupledDomainRuntime.hpp"
#include "FsiCouplingEdge.hpp"
#include "FsiDomainRuntime.hpp"

#include <cassert>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
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
const std::string kOtherMesh(64, 'd');
const std::string kFluidState(64, 'b');
const std::string kStructureState(64, 'e');
const std::string kProjection(64, 'c');

iga::DistributedSurfaceLayout MakeLayout(const std::string& mesh = kMesh)
{
	iga::DistributedSurfaceLayout layout;
	layout.reference_mesh_identity_sha256 = mesh;
	layout.global_node_count = 3;
	layout.partition_count = 1;
	layout.partition_rank = 0;
	layout.owned_global_node_ids = {10, 20, 30};
	layout.reference_positions = {{10, {{0.0, 0.0, 0.0}}}, {20, {{1.0, 0.0, 0.0}}},
		{30, {{0.0, 1.0, 0.0}}}};
	layout.reference_triangles = {{{10, 20, 30}}};
	layout.owned_reference_lumped_areas_m2 = {1.0, 2.0, 3.0};
	layout.layout_identity_sha256 = iga::BuildDistributedSurfaceLayoutIdentitySha256(layout);
	return layout;
}

iga::SurfaceFieldStamp MakeStamp(const iga::DistributedSurfaceLayout& layout,
	const std::string& producer_state, std::uint64_t step = 7, std::uint64_t iteration = 2)
{
	iga::SurfaceFieldStamp stamp;
	stamp.time_s = 0.2;
	stamp.step = step;
	stamp.coupling_iteration = iteration;
	stamp.reference_mesh_identity_sha256 = layout.reference_mesh_identity_sha256;
	stamp.layout_identity_sha256 = layout.layout_identity_sha256;
	stamp.partition_identity_sha256 = iga::BuildDistributedSurfacePartitionIdentitySha256(layout);
	stamp.producer_state_identity_sha256 = producer_state;
	return stamp;
}

iga::DistributedSurfaceInterface MakeFluidSurface(const std::string& id = "wall",
	const std::string& mesh = kMesh)
{
	iga::DistributedSurfaceInterface surface;
	surface.id = {"fluid", "fluid_backend", id};
	surface.subsystem_id = "fluid_backend";
	surface.boundary_labels = {0, 7};
	surface.reference_mesh_identity_sha256 = mesh;
	surface.provides = {iga::SurfaceFieldQuantity::TractionOnStructure};
	surface.requires = {iga::SurfaceFieldQuantity::Displacement,
		iga::SurfaceFieldQuantity::Velocity};
	return surface;
}

iga::DistributedSurfaceInterface MakeStructureSurface(const std::string& id = "wall",
	const std::string& mesh = kMesh)
{
	iga::DistributedSurfaceInterface surface;
	surface.id = {"structure", "structure_backend", id};
	surface.subsystem_id = "structure_backend";
	surface.boundary_labels = {0, 7};
	surface.reference_mesh_identity_sha256 = mesh;
	surface.provides = {iga::SurfaceFieldQuantity::Displacement,
		iga::SurfaceFieldQuantity::Velocity};
	surface.requires = {iga::SurfaceFieldQuantity::TractionOnStructure};
	return surface;
}

iga::SurfaceKinematics MakeKinematics(const iga::SurfaceInterfaceRef& interface,
	const iga::DistributedSurfaceLayout&, const iga::SurfaceFieldStamp& stamp)
{
	iga::SurfaceKinematics value;
	value.interface = interface;
	value.stamp = stamp;
	value.displacement_m = {{{0.0, 0.0, 0.01}}, {{0.0, 0.0, 0.02}}, {{0.0, 0.0, 0.03}}};
	value.velocity_m_per_s = {{{0.0, 0.0, 1.0}}, {{0.0, 0.0, 2.0}}, {{0.0, 0.0, 3.0}}};
	return value;
}

iga::SurfaceTraction MakeTraction(const iga::SurfaceInterfaceRef& interface,
	const iga::DistributedSurfaceLayout&, const iga::SurfaceFieldStamp& stamp)
{
	iga::SurfaceTraction value;
	value.interface = interface;
	value.stamp = stamp;
	value.traction_on_structure_pa = {{{0.0, 0.0, -10.0}}, {{0.0, 0.0, -20.0}}, {{0.0, 0.0, -30.0}}};
	value.consistent_nodal_force_n = {{{0.0, 0.0, -1.0}}, {{0.0, 0.0, -2.0}}, {{0.0, 0.0, -3.0}}};
	value.projection_identity_sha256 = kProjection;
	return value;
}

class FluidMock final : public iga::CoupledDomainRuntime, public iga::FsiFluidDomainRuntime {
public:
	FluidMock(iga::FsiCouplingEdge edge, iga::DistributedSurfaceInterface surface,
		iga::DistributedSurfaceInterface peer, iga::DistributedSurfaceLayout layout,
		iga::SurfaceFieldStamp expected_input, iga::SurfaceFieldStamp expected_output)
		: edge_(std::move(edge)), catalog_{std::move(surface)}, layout_(std::move(layout)),
		  expected_input_(std::move(expected_input)), expected_output_(std::move(expected_output)),
		  lifecycle_("fluid", "fluid_backend", edge_, catalog_.front(), peer, layout_, layout_),
		  traction_(MakeTraction(edge_.fluid, layout_, expected_output_))
	{
		port_.id = "outlet";
		port_.subsystem_id = "fluid_backend";
		port_.locator_kind = "boundary";
		port_.locator = "outlet";
		port_.provides = {iga::PortQuantity::FlowRate};
		ports_ = {port_};
	}

	const std::string& DomainId() const noexcept override { return domain_id_; }
	iga::DomainKind Kind() const noexcept override { return iga::DomainKind::ThreeDBodyFittedFlow; }
	const std::vector<iga::CouplingPort>& Ports() const noexcept override { return ports_; }
	void BeginStep(const iga::DomainStepContext& step) override { lifecycle_.BeginStep(step); }
	void SetPortInput(const std::string& port_id, const iga::PortBoundaryData& input) override
	{
		if (port_id != port_.id) throw std::runtime_error("unknown scalar port");
		iga::ValidatePortBoundaryData(input);
		last_scalar_input_ = input;
	}
	void SolveTrial() override { lifecycle_.RequireSolveAllowed(); lifecycle_.MarkSolved(edge_.fluid, traction_.stamp); }
	iga::PortState GetPortState(const std::string& port_id) const override
	{
		if (port_id != port_.id) throw std::runtime_error("unknown scalar port");
		return {};
	}
	void RollbackTrial() override { lifecycle_.RollbackTrial(); }
	void AbortStep() override { lifecycle_.AbortStep(); }
	void PrepareCommitStep() override { lifecycle_.PrepareCommit(); }
	void FinalizeCommitStep() noexcept override { lifecycle_.FinalizeCommit(); }

	void BeginIteration(std::uint64_t iteration)
	{ lifecycle_.BeginIteration(iteration, expected_input_, iga::MakeSurfaceFieldStampEnvelope(expected_output_)); }
	const iga::FsiTrialLifecycle& Lifecycle() const { return lifecycle_; }
	const std::vector<iga::DistributedSurfaceInterface>& SurfaceInterfaces() const noexcept override { return catalog_; }
	void SetSurfaceKinematics(const std::string& interface_id, const iga::SurfaceKinematics& kinematics) override
	{
		if (interface_id != edge_.fluid.interface_id) throw std::runtime_error("unknown fluid surface interface");
		iga::ValidateFsiFluidKinematicsInput(edge_, kinematics, layout_, expected_input_);
		lifecycle_.MarkInput(kinematics.interface, kinematics.stamp);
		kinematics_ = kinematics;
	}
	iga::SurfaceTraction GetSurfaceTraction(const std::string& interface_id) const override
	{
		if (interface_id != edge_.fluid.interface_id) throw std::runtime_error("unknown fluid surface interface");
		lifecycle_.RequireTrialOutput(traction_.interface, traction_.stamp);
		return traction_;
	}
	bool HasScalarInput() const { return last_scalar_input_.has_value(); }

private:
	std::string domain_id_ = "fluid";
	iga::FsiCouplingEdge edge_;
	std::vector<iga::DistributedSurfaceInterface> catalog_;
	iga::DistributedSurfaceLayout layout_;
	iga::SurfaceFieldStamp expected_input_;
	iga::SurfaceFieldStamp expected_output_;
	iga::FsiTrialLifecycle lifecycle_;
	iga::SurfaceTraction traction_;
	std::optional<iga::SurfaceKinematics> kinematics_;
	iga::CouplingPort port_;
	std::vector<iga::CouplingPort> ports_;
	std::optional<iga::PortBoundaryData> last_scalar_input_;
};

class StructureMock final : public iga::FsiStructureDomainRuntime {
public:
	StructureMock(iga::FsiCouplingEdge edge, iga::DistributedSurfaceInterface surface,
		iga::DistributedSurfaceInterface peer, iga::DistributedSurfaceLayout layout,
		iga::SurfaceFieldStamp expected_input, iga::SurfaceFieldStamp expected_output)
		: edge_(std::move(edge)), catalog_{std::move(surface)}, layout_(std::move(layout)),
		  expected_input_(std::move(expected_input)), expected_output_(std::move(expected_output)),
		  lifecycle_("structure", "structure_backend", edge_, catalog_.front(), peer, layout_, layout_),
		  kinematics_(MakeKinematics(edge_.structure, layout_, expected_output_)) {}

	void BeginStep(const iga::DomainStepContext& step) { lifecycle_.BeginStep(step); }
	void BeginIteration(std::uint64_t iteration)
	{ lifecycle_.BeginIteration(iteration, expected_input_, iga::MakeSurfaceFieldStampEnvelope(expected_output_)); }
	void SolveTrial() { lifecycle_.RequireSolveAllowed(); lifecycle_.MarkSolved(edge_.structure, kinematics_.stamp); }
	void RollbackTrial() { lifecycle_.RollbackTrial(); }
	void AbortStep() { lifecycle_.AbortStep(); }
	void PrepareCommitStep() { lifecycle_.PrepareCommit(); }
	void FinalizeCommitStep() { lifecycle_.FinalizeCommit(); }
	const iga::FsiTrialLifecycle& Lifecycle() const { return lifecycle_; }
	const std::vector<iga::DistributedSurfaceInterface>& SurfaceInterfaces() const noexcept override { return catalog_; }
	void SetSurfaceTraction(const std::string& interface_id, const iga::SurfaceTraction& traction) override
	{
		if (interface_id != edge_.structure.interface_id) throw std::runtime_error("unknown structure surface interface");
		iga::ValidateFsiStructureTractionInput(edge_, traction, layout_, expected_input_);
		lifecycle_.MarkInput(traction.interface, traction.stamp);
		traction_ = traction;
	}
	iga::SurfaceKinematics GetSurfaceKinematics(const std::string& interface_id) const override
	{
		if (interface_id != edge_.structure.interface_id) throw std::runtime_error("unknown structure surface interface");
		lifecycle_.RequireTrialOutput(kinematics_.interface, kinematics_.stamp);
		return kinematics_;
	}

private:
	iga::FsiCouplingEdge edge_;
	std::vector<iga::DistributedSurfaceInterface> catalog_;
	iga::DistributedSurfaceLayout layout_;
	iga::SurfaceFieldStamp expected_input_;
	iga::SurfaceFieldStamp expected_output_;
	iga::FsiTrialLifecycle lifecycle_;
	iga::SurfaceKinematics kinematics_;
	std::optional<iga::SurfaceTraction> traction_;
};

} // namespace

int main()
{
	const auto layout = MakeLayout();
	const auto fluid_surface = MakeFluidSurface();
	const auto structure_surface = MakeStructureSurface();
	const iga::FsiCouplingEdge edge("fsi_wall", fluid_surface.id, structure_surface.id,
		iga::FsiCouplingLaw::FluidStructureTractionKinematics);
	iga::ValidateFsiCouplingEdgeLayouts(edge, fluid_surface, structure_surface, layout, layout);
	const std::string identity = iga::BuildFsiCouplingEdgeIdentitySha256(edge);
	assert(identity == iga::BuildFsiCouplingEdgeIdentitySha256(edge));
	auto subsystem_mutation = edge;
	subsystem_mutation.fluid.subsystem_id = "other_backend";
	assert(identity != iga::BuildFsiCouplingEdgeIdentitySha256(subsystem_mutation));
	auto swapped = edge;
	std::swap(swapped.fluid, swapped.structure);
	assert(identity != iga::BuildFsiCouplingEdgeIdentitySha256(swapped));
	iga::ValidateFsiCouplingEdgeBindings({edge}, {fluid_surface}, {structure_surface});

	RequireRejected([&] { auto invalid = fluid_surface; invalid.subsystem_id = "other";
		iga::ValidateFsiCouplingEdgeInterfaces(edge, invalid, structure_surface); });
	auto first_slice = layout;
	first_slice.partition_count = 2;
	first_slice.owned_global_node_ids = {10, 20};
	first_slice.owned_reference_lumped_areas_m2 = {1.0, 2.0};
	first_slice.layout_identity_sha256 = iga::BuildDistributedSurfaceLayoutIdentitySha256(first_slice);
	auto different_slice = first_slice;
	different_slice.owned_global_node_ids = {20, 30};
	different_slice.owned_reference_lumped_areas_m2 = {2.0, 3.0};
	assert(first_slice.layout_identity_sha256 == different_slice.layout_identity_sha256);
	assert(iga::BuildDistributedSurfacePartitionIdentitySha256(first_slice)
		!= iga::BuildDistributedSurfacePartitionIdentitySha256(different_slice));
	RequireRejected([&] { iga::ValidateFsiCouplingEdgeLayouts(edge, fluid_surface, structure_surface,
		first_slice, different_slice); });
	RequireRejected([&] { auto other_structure = MakeStructureSurface("other");
		iga::FsiCouplingEdge duplicate("another", fluid_surface.id, other_structure.id, edge.law);
		iga::ValidateFsiCouplingEdgeBindings({edge, duplicate}, {fluid_surface}, {structure_surface, other_structure}); });
	RequireRejected([&] { auto other_fluid = MakeFluidSurface("other_fluid");
		auto other_structure = MakeStructureSurface("other_structure");
		iga::FsiCouplingEdge duplicate_id(edge.id, other_fluid.id, other_structure.id, edge.law);
		iga::ValidateFsiCouplingEdgeBindings({edge, duplicate_id}, {fluid_surface, other_fluid},
			{structure_surface, other_structure}); });
	RequireRejected([&] { iga::ValidateFsiCouplingEdge({"bad", {"fluid", "fluid_backend", "wall"},
		{"fluid", "structure_backend", "wall"}, edge.law}); });
	RequireRejected([&] { auto other = layout; other.reference_positions[1].position_m[0] = 2.0;
		other.layout_identity_sha256 = iga::BuildDistributedSurfaceLayoutIdentitySha256(other);
		iga::ValidateFsiCouplingEdgeLayouts(edge, fluid_surface, structure_surface, layout, other); });

	const auto fluid_output = MakeStamp(layout, kFluidState);
	const auto structure_output = MakeStamp(layout, kStructureState);
	FluidMock fluid(edge, fluid_surface, structure_surface, layout, structure_output, fluid_output);
	StructureMock structure(edge, structure_surface, fluid_surface, layout, fluid_output, structure_output);
	assert(fluid.DomainId() == "fluid" && fluid.Ports().size() == 1);
	iga::PortBoundaryData scalar_input;
	scalar_input.time_s = 0.2;
	fluid.SetPortInput("outlet", scalar_input);
	assert(fluid.HasScalarInput());

	RequireRejected([&] { fluid.GetSurfaceTraction("wall"); });
	RequireRejected([&] { structure.GetSurfaceKinematics("wall"); });
	RequireRejected([&] { fluid.BeginIteration(2); });
	iga::DomainStepContext step{7, 0.125, 0.075};
	fluid.BeginStep(step);
	structure.BeginStep(step);
	fluid.BeginIteration(2);
	structure.BeginIteration(2);
	const std::string first_context = fluid.Lifecycle().ContextIdentitySha256();
	assert(first_context == fluid.Lifecycle().ContextIdentitySha256());
	RequireRejected([&] { fluid.SolveTrial(); });
	RequireRejected([&] { fluid.GetSurfaceTraction("wall"); });
	RequireRejected([&] { fluid.SetSurfaceKinematics("wall", MakeKinematics(fluid_surface.id, layout, structure_output)); });
	const auto kinematics = MakeKinematics(structure_surface.id, layout, structure_output);
	fluid.SetSurfaceKinematics("wall", kinematics);
	RequireRejected([&] { fluid.SetSurfaceKinematics("wall", kinematics); });
	fluid.SolveTrial();
	const auto traction = fluid.GetSurfaceTraction("wall");
	RequireRejected([&] { structure.GetSurfaceKinematics("wall"); });
	structure.SetSurfaceTraction("wall", traction);
	RequireRejected([&] { structure.SetSurfaceTraction("wall", traction); });
	structure.SolveTrial();
	assert(structure.GetSurfaceKinematics("wall").interface == structure_surface.id);
	fluid.PrepareCommitStep();
	structure.PrepareCommitStep();
	assert(!fluid.Lifecycle().HasCommittedOutput());
	RequireRejected([&] { fluid.Lifecycle().RequireCommittedOutput(fluid_surface.id, fluid_output); });
	fluid.FinalizeCommitStep();
	structure.FinalizeCommitStep();
	assert(fluid.Lifecycle().HasCommittedOutput());
	fluid.Lifecycle().RequireCommittedOutput(fluid_surface.id, fluid_output);
	RequireRejected([&] { fluid.GetSurfaceTraction("wall"); });

	fluid.BeginStep(step);
	fluid.BeginIteration(2);
	auto iteration_three_kinematics = kinematics;
	iteration_three_kinematics.stamp = MakeStamp(layout, kStructureState, 7, 3);
	RequireRejected([&] { fluid.SetSurfaceKinematics("wall", iteration_three_kinematics); });
	fluid.RollbackTrial();
	RequireRejected([&] { fluid.GetSurfaceTraction("wall"); });
	fluid.BeginIteration(2);
	fluid.SetSurfaceKinematics("wall", kinematics);
	fluid.AbortStep();
	assert(fluid.Lifecycle().Phase() == iga::FsiTrialPhase::Idle);
	RequireRejected([&] { fluid.GetSurfaceTraction("wall"); });
	RequireRejected([&] { fluid.Lifecycle().ContextIdentitySha256(); });

	return 0;
}
