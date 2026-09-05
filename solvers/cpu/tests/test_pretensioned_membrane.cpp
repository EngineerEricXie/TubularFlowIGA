#include "PretensionedMembrane.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <vector>

static_assert(!std::is_copy_constructible<iga::PretensionedMembrane>::value,
	"PretensionedMembrane must retain unique lifecycle ownership");
static_assert(!std::is_copy_assignable<iga::PretensionedMembrane>::value,
	"PretensionedMembrane must retain unique lifecycle ownership");
static_assert(!std::is_move_constructible<iga::PretensionedMembrane>::value,
	"PretensionedMembrane must retain unique lifecycle ownership");
static_assert(!std::is_move_assignable<iga::PretensionedMembrane>::value,
	"PretensionedMembrane must retain unique lifecycle ownership");

namespace {

bool Near(double left, double right, double tolerance = 2.0e-9)
{ return std::abs(left-right) <= tolerance*std::max({1.0, std::abs(left), std::abs(right)}); }

template <class Function> void Reject(Function&& function)
{ bool rejected = false; try { function(); } catch (const std::exception&) { rejected = true; } assert(rejected); }

std::string HashText(const char* text)
{ iga::Sha256 hash; hash.Append(text, std::char_traits<char>::length(text)); return hash.Hex(); }

iga::DistributedSurfaceInterface StructureInterface(const std::string& mesh)
{
	iga::DistributedSurfaceInterface value;
	value.id = {"structure", "membrane", "patch"}; value.subsystem_id = "membrane";
	value.boundary_labels = {0}; value.reference_mesh_identity_sha256 = mesh;
	value.provides = {iga::SurfaceFieldQuantity::Displacement, iga::SurfaceFieldQuantity::Velocity};
	value.requires = {iga::SurfaceFieldQuantity::TractionOnStructure}; return value;
}

iga::SurfaceInterfaceRef FluidInterface() { return {"fluid", "immersed", "patch"}; }

iga::DistributedSurfaceLayout SquareLayout(int intervals)
{
	iga::DistributedSurfaceLayout value;
	value.reference_mesh_identity_sha256 = HashText("pretensioned-membrane-square");
	value.global_node_count = static_cast<std::uint64_t>((intervals+1)*(intervals+1));
	value.partition_count = 1; value.partition_rank = 0;
	for (int y = 0; y <= intervals; ++y) for (int x = 0; x <= intervals; ++x) {
		const auto id = static_cast<std::uint64_t>(y*(intervals+1)+x);
		value.owned_global_node_ids.push_back(id);
		value.reference_positions.push_back({id, {{double(x)/intervals, double(y)/intervals, 0.0}}});
	}
	for (int y = 0; y < intervals; ++y) for (int x = 0; x < intervals; ++x) {
		const auto a = static_cast<std::uint64_t>(y*(intervals+1)+x), b = a+1;
		const auto c = static_cast<std::uint64_t>((y+1)*(intervals+1)+x), d = c+1;
		value.reference_triangles.push_back({{a,b,d}}); value.reference_triangles.push_back({{a,d,c}});
	}
	value.owned_reference_lumped_areas_m2.assign(value.owned_global_node_ids.size(), 0.0);
	for (const auto& triangle : value.reference_triangles) {
		const auto& a = value.reference_positions[triangle[0]].position_m;
		const auto& b = value.reference_positions[triangle[1]].position_m;
		const auto& c = value.reference_positions[triangle[2]].position_m;
		const double area = .5*std::abs((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]));
		for (const auto id : triangle) value.owned_reference_lumped_areas_m2[id] += area/3.0;
	}
	value.layout_identity_sha256 = iga::BuildDistributedSurfaceLayoutIdentitySha256(value); return value;
}

std::vector<std::uint64_t> SquareBoundary(int intervals)
{
	std::vector<std::uint64_t> result;
	for (int y = 0; y <= intervals; ++y) for (int x = 0; x <= intervals; ++x)
		if (x == 0 || y == 0 || x == intervals || y == intervals) result.push_back(static_cast<std::uint64_t>(y*(intervals+1)+x));
	return result;
}

iga::SurfaceFieldStamp Stamp(const iga::DistributedSurfaceLayout& layout, double time,
	std::uint64_t step, std::uint64_t iteration, const char* tag)
{
	iga::SurfaceFieldStamp value; value.time_s = time; value.step = step; value.coupling_iteration = iteration;
	value.reference_mesh_identity_sha256 = layout.reference_mesh_identity_sha256;
	value.layout_identity_sha256 = layout.layout_identity_sha256;
	value.partition_identity_sha256 = iga::BuildDistributedSurfacePartitionIdentitySha256(layout);
	value.producer_state_identity_sha256 = HashText(tag); return value;
}

iga::SurfaceTraction Traction(const iga::DistributedSurfaceLayout& layout, const iga::SurfaceFieldStamp& stamp,
	const std::vector<double>& scalar_force)
{
	(void)layout;
	iga::SurfaceTraction value; value.interface = FluidInterface(); value.stamp = stamp;
	value.traction_on_structure_pa.assign(scalar_force.size(), {{0.0, 0.0, 0.0}});
	value.consistent_nodal_force_n.resize(scalar_force.size());
	for (std::size_t node = 0; node < scalar_force.size(); ++node)
		value.consistent_nodal_force_n[node] = {{0.0, 0.0, scalar_force[node]}};
	value.projection_identity_sha256 = HashText("known-projection"); return value;
}

iga::PretensionedMembraneTrialContext Context(const iga::DistributedSurfaceLayout& layout,
	double start, double dt, std::uint64_t step, std::uint64_t iteration, const char* tag)
{
	iga::PretensionedMembraneTrialContext value;
	value.start_time_s = start; value.dt_s = dt; value.step = step; value.coupling_iteration = iteration;
	value.expected_traction_stamp = Stamp(layout, start+dt, step, iteration, tag); return value;
}

iga::PretensionedMembrane MakeSquare(int intervals, iga::PretensionedMembraneMaterial material,
	iga::PretensionedMembraneState state = {})
{
	auto layout = SquareLayout(intervals);
	return iga::PretensionedMembrane(layout, StructureInterface(layout.reference_mesh_identity_sha256),
		FluidInterface(), material, SquareBoundary(intervals), {4096}, std::move(state));
}

void ZeroLoadAndClamps()
{
	iga::PretensionedMembraneMaterial material{2.0, .3, 4.0, 1.0}; auto membrane = MakeSquare(2, material);
	auto context = Context(membrane.Layout(), 1.0, .05, 8, 2, "zero");
	auto trial = membrane.SolveTrial(context, Traction(membrane.Layout(), context.expected_traction_stamp,
		std::vector<double>(membrane.Layout().owned_global_node_ids.size(), 0.0)));
	for (const auto value : trial.state.displacement_m) assert(value == 0.0);
	for (const auto value : trial.state.velocity_m_per_s) assert(value == 0.0);
	for (const auto id : SquareBoundary(2)) {
		assert(trial.state.displacement_m[id] == 0.0 && trial.state.velocity_m_per_s[id] == 0.0);
		assert((trial.kinematics.displacement_m[id] == std::array<double,3>{{0.0,0.0,0.0}}));
	}
}

void LinearityAndRetry()
{
	iga::PretensionedMembraneMaterial material{1.0, .1, 3.0, 2.0}; auto first = MakeSquare(2, material);
	auto second = MakeSquare(2, material), combined = MakeSquare(2, material);
	auto context = Context(first.Layout(), 0.0, .1, 1, 0, "linearity");
	std::vector<double> a(first.Layout().global_node_count, 0.0), b(a); a[4] = 1.25; b[4] = -.35;
	auto trial_a = first.SolveTrial(context, Traction(first.Layout(), context.expected_traction_stamp, a));
	auto trial_b = second.SolveTrial(context, Traction(second.Layout(), context.expected_traction_stamp, b));
	for (double& value : a) value += b[&value-&a[0]];
	auto trial_ab = combined.SolveTrial(context, Traction(combined.Layout(), context.expected_traction_stamp, a));
	for (std::size_t node = 0; node < a.size(); ++node) assert(Near(trial_ab.state.displacement_m[node],
		trial_a.state.displacement_m[node]+trial_b.state.displacement_m[node]));
	const auto original = trial_a.trial_identity_sha256; first.RejectTrial(trial_a);
	auto retry = first.SolveTrial(context, Traction(first.Layout(), context.expected_traction_stamp, std::vector<double>(a.size(), 0.0)));
	auto clean = MakeSquare(2, material).SolveTrial(context, Traction(first.Layout(), context.expected_traction_stamp, std::vector<double>(a.size(), 0.0)));
	assert(retry.trial_identity_sha256 != clean.trial_identity_sha256 && original != retry.trial_identity_sha256);
	assert(retry.state.displacement_m == clean.state.displacement_m && retry.state.velocity_m_per_s == clean.state.velocity_m_per_s);
	Reject([&] { first.SolveTrial(context, Traction(first.Layout(), context.expected_traction_stamp, std::vector<double>(a.size(), 0.0))); });
	first.AbortTrial();
	auto accepted = first.SolveTrial(context, Traction(first.Layout(), context.expected_traction_stamp, a));
	auto prepared = first.PrepareTrial(accepted); first.FinalizeTrial(prepared);
	assert(first.CommittedStateIdentitySha256() != prepared.base_committed_state_identity_sha256);
	Reject([&] { first.FinalizeTrial(prepared); });
}

double SpatialError(int intervals)
{
	iga::PretensionedMembraneMaterial material{1.0, 0.0, 1.0, 1.0}; auto membrane = MakeSquare(intervals, material);
	const auto& layout = membrane.Layout(); std::vector<double> load(layout.global_node_count);
	// Assemble q_i = int N_i f_h dA: a consistent P1 nodal force, with
	// f = (2 pi^2 + 1) sin(pi x) sin(pi y) for (-Delta+1)d=f.
	for (const auto& triangle : layout.reference_triangles) {
		std::array<double, 3> source{};
		for (int local = 0; local < 3; ++local) {
			const auto& x = layout.reference_positions[triangle[local]].position_m;
			source[local] = (2.0*M_PI*M_PI+1.0)*std::sin(M_PI*x[0])*std::sin(M_PI*x[1]);
		}
		const auto& a = layout.reference_positions[triangle[0]].position_m;
		const auto& b = layout.reference_positions[triangle[1]].position_m;
		const auto& c = layout.reference_positions[triangle[2]].position_m;
		const double area = .5*std::abs((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]));
		for (int row = 0; row < 3; ++row) for (int column = 0; column < 3; ++column)
			load[triangle[row]] += area*(row == column ? 2.0 : 1.0)*source[column]/12.0;
	}
	auto context = Context(layout, 0.0, 1.0e6, 0, 0, "spatial");
	auto trial = membrane.SolveTrial(context, Traction(layout, context.expected_traction_stamp, load));
	double squared = 0.0, weight = 0.0;
	for (std::size_t node = 0; node < load.size(); ++node) {
		const auto& x = layout.reference_positions[node].position_m;
		const double error = trial.state.displacement_m[node]-std::sin(M_PI*x[0])*std::sin(M_PI*x[1]);
		squared += layout.owned_reference_lumped_areas_m2[node]*error*error; weight += layout.owned_reference_lumped_areas_m2[node];
	}
	return std::sqrt(squared/weight);
}

void SpatialRefinement()
{
	const double coarse = SpatialError(6), middle = SpatialError(10), fine = SpatialError(14);
	const double rate_coarse = std::log(coarse/middle)/std::log(10.0/6.0);
	const double rate_fine = std::log(middle/fine)/std::log(14.0/10.0);
	std::cout << "spatial errors " << coarse << " " << middle << " " << fine << "; rates " << rate_coarse << " " << rate_fine << "\n";
	assert(middle < coarse && fine < middle && rate_coarse > 1.25 && rate_fine > 1.25);
}

iga::PretensionedMembrane MakeModalMembrane()
{
	iga::DistributedSurfaceLayout layout; layout.reference_mesh_identity_sha256 = HashText("modal");
	layout.global_node_count = 5; layout.partition_count = 1; layout.partition_rank = 0;
	layout.owned_global_node_ids = {0,1,2,3,4};
	layout.reference_positions = {{0,{{0,0,0}}},{1,{{1,0,0}}},{2,{{1,1,0}}},{3,{{0,1,0}}},{4,{{.5,.5,0}}}};
	layout.reference_triangles = {{{{0,1,4}},{{1,2,4}},{{2,3,4}},{{3,0,4}}}};
	layout.owned_reference_lumped_areas_m2 = {1.0/6,1.0/6,1.0/6,1.0/6,1.0/3};
	layout.layout_identity_sha256 = iga::BuildDistributedSurfaceLayoutIdentitySha256(layout);
	return iga::PretensionedMembrane(layout, StructureInterface(layout.reference_mesh_identity_sha256), FluidInterface(),
		{2.0, 20.0, 1.0, 3.0}, {0,1,2,3});
}

double ModalExact(double time)
{
	const double mass = 2.0/6.0, damping = 20.0/6.0, stiffness = 4.0+3.0/6.0, force = 1.0;
	const double discriminant = std::sqrt(damping*damping-4.0*mass*stiffness);
	const double first = (-damping+discriminant)/(2.0*mass), second = (-damping-discriminant)/(2.0*mass), steady = force/stiffness;
	const double a = steady*second/(first-second), b = -steady-a;
	return steady+a*std::exp(first*time)+b*std::exp(second*time);
}

double TemporalError(int steps)
{
	auto membrane = MakeModalMembrane(); const double dt = 1.0/steps;
	for (int step = 0; step < steps; ++step) {
		auto context = Context(membrane.Layout(), step*dt, dt, step, 0, "temporal");
		std::vector<double> force(5, 0.0); force[4] = 1.0;
		auto trial = membrane.SolveTrial(context, Traction(membrane.Layout(), context.expected_traction_stamp, force));
		membrane.FinalizeTrial(membrane.PrepareTrial(trial));
	}
	return std::abs(membrane.CommittedState().displacement_m[4]-ModalExact(1.0));
}

void TemporalRefinement()
{
	const double coarse = TemporalError(20), middle = TemporalError(40), fine = TemporalError(80);
	const double first = std::log(coarse/middle)/std::log(2.0), second = std::log(middle/fine)/std::log(2.0);
	std::cout << "temporal errors " << coarse << " " << middle << " " << fine << "; rates " << first << " " << second << "\n";
	assert(middle < coarse && fine < middle && first > .75 && second > .75);
}

void ContractRejections()
{
	auto layout = SquareLayout(2); auto interface = StructureInterface(layout.reference_mesh_identity_sha256);
	iga::PretensionedMembraneMaterial material{1.0, 0.0, 1.0, 1.0};
	auto bad_capability = interface; bad_capability.requires.clear();
	Reject([&] { iga::PretensionedMembrane(layout, bad_capability, FluidInterface(), material, SquareBoundary(2)); });
	auto distributed = layout; distributed.partition_count = 2; distributed.layout_identity_sha256 = iga::BuildDistributedSurfaceLayoutIdentitySha256(distributed);
	Reject([&] { iga::PretensionedMembrane(distributed, StructureInterface(distributed.reference_mesh_identity_sha256), FluidInterface(), material, SquareBoundary(2)); });
	Reject([&] { iga::PretensionedMembrane(layout, interface, FluidInterface(), {0.0,0.0,1.0,0.0}, SquareBoundary(2)); });
	Reject([&] { MakeSquare(2, {1.0,0.0,0.0,0.0}); });
	MakeSquare(2, {1.0,0.0,1.0e-200,0.0});
	MakeSquare(2, {1.0,0.0,0.0,1.0e-200});
	auto membrane = iga::PretensionedMembrane(layout, interface, FluidInterface(), material, SquareBoundary(2));
	auto context = Context(layout, 0.0, .1, 2, 3, "contract"); std::vector<double> force(layout.global_node_count, 0.0);
	auto stale = Traction(layout, context.expected_traction_stamp, force); stale.stamp.step += 1;
	Reject([&] { membrane.SolveTrial(context, stale); });
	// A failed solve leaves no lifecycle capability behind.
	auto valid_after_failure = membrane.SolveTrial(context, Traction(layout, context.expected_traction_stamp, force));
	membrane.RejectTrial(valid_after_failure);
	auto wrong_interface = Traction(layout, context.expected_traction_stamp, force); wrong_interface.interface.domain_id = "other";
	Reject([&] { membrane.SolveTrial(context, wrong_interface); });
	auto first = membrane.SolveTrial(context, Traction(layout, context.expected_traction_stamp, force));
	Reject([&] { membrane.SolveTrial(context, Traction(layout, context.expected_traction_stamp, force)); });
	membrane.RejectTrial(first);
	Reject([&] { membrane.FinalizeTrial(first); });
	auto retry_context = Context(layout, 0.0, .1, 2, 4, "contract-retry");
	auto coupling_retry = membrane.SolveTrial(retry_context, Traction(layout, retry_context.expected_traction_stamp, force));
	assert(coupling_retry.base_committed_state_identity_sha256 == first.base_committed_state_identity_sha256);
	assert(coupling_retry.state.displacement_m == first.state.displacement_m
		&& coupling_retry.state.velocity_m_per_s == first.state.velocity_m_per_s);
	auto forged = coupling_retry; forged.state.displacement_m[4] = 1.0;
	Reject([&] { membrane.PrepareTrial(forged); });
	auto foreign_owner = MakeSquare(2, material);
	auto foreign = foreign_owner.SolveTrial(retry_context, Traction(layout, retry_context.expected_traction_stamp, force));
	Reject([&] { membrane.PrepareTrial(foreign); });
	auto prepared = membrane.PrepareTrial(coupling_retry);
	Reject([&] { membrane.PrepareTrial(coupling_retry); });
	Reject([&] { membrane.FinalizeTrial(coupling_retry); });
	Reject([&] { membrane.RejectTrial(prepared); });
	membrane.AbortTrial();
	Reject([&] { membrane.FinalizeTrial(prepared); });
	auto replacement = membrane.SolveTrial(context, Traction(layout, context.expected_traction_stamp, force));
	membrane.RejectTrial(replacement);
}

} // namespace

int main()
{
	ZeroLoadAndClamps(); LinearityAndRetry(); SpatialRefinement(); TemporalRefinement(); ContractRejections();
	std::cout << "pretensioned membrane tests passed\n";
}
