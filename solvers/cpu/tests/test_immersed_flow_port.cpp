#include "ImmersedStaticFlowRuntime.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

iga::RawSurfaceTriangle Face(std::int64_t a, std::int64_t b, std::int64_t c, std::uint32_t label)
{ iga::RawSurfaceTriangle result; result.indices = {{a,b,c}}; result.boundary_id = label; return result; }

iga::RawSurfaceSoup PortCube()
{
	iga::RawSurfaceSoup result;
	result.vertices = {{{{.1,.1,.1}},{{.9,.1,.1}},{{.9,.9,.1}},{{.1,.9,.1}},
		{{.1,.1,.9}},{{.9,.1,.9}},{{.9,.9,.9}},{{.1,.9,.9}}}};
	// z-low and z-high are open caps; every other closed-surface triangle is a wall.
	result.triangles = {Face(0,2,1,1),Face(0,3,2,1),Face(4,5,6,2),Face(4,6,7,2),
		Face(0,1,5,0),Face(0,5,4,0),Face(1,2,6,0),Face(1,6,5,0),
		Face(2,3,7,0),Face(2,7,6,0),Face(3,0,4,0),Face(3,4,7,0)};
	return result;
}

iga::CartesianDomainClassification Domain()
{
	const iga::CubicCartesianGridSpec spec{{{0,0,0}},{{1,1,1}},{{3,3,3}}};
	return iga::CartesianDomainClassification(iga::CubicCartesianBackground(spec),
		iga::SurfaceSpatialIndex(iga::ClosedTriangulatedSurface::Build(PortCube())));
}

iga::RawSurfaceSoup ZeroPortCube()
{
	auto result = PortCube();
	// Label zero is an open z-low port here; the side wall deliberately uses a
	// different label so this exercises the valid non-overlapping case.
	result.triangles[0].boundary_id = 0;
	result.triangles[1].boundary_id = 0;
	for (std::size_t i = 4; i < result.triangles.size(); ++i) result.triangles[i].boundary_id = 3;
	return result;
}

iga::CartesianDomainClassification ZeroPortDomain()
{
	const iga::CubicCartesianGridSpec spec{{{0,0,0}},{{1,1,1}},{{3,3,3}}};
	return iga::CartesianDomainClassification(iga::CubicCartesianBackground(spec),
		iga::SurfaceSpatialIndex(iga::ClosedTriangulatedSurface::Build(ZeroPortCube())));
}

template <class Function> void Reject(Function&& function)
{ bool rejected = false; try { function(); } catch (const std::exception&) { rejected = true; } assert(rejected); }

const iga::ImmersedStaticFlowDiagnostics::Port& Port(const iga::ImmersedStaticFlowRuntime& runtime, const char* id)
{
	const auto& ports = runtime.Diagnostics().ports;
	const auto found = std::find_if(ports.begin(), ports.end(), [id](const auto& port) { return port.id == id; });
	assert(found != ports.end()); return *found;
}

bool HasLabel(const iga::SurfaceQuadratureRule& rule, int label)
{
	return std::any_of(rule.Points().begin(), rule.Points().end(), [label](const iga::SurfaceQuadraturePoint& point) { return point.boundary_id == label; });
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int status = 0;
	try {
		// The topology count seam exercises PETSc-index boundaries without
		// materializing a correspondingly enormous Cartesian background.
		const std::size_t petsc_max = static_cast<std::size_t>(std::numeric_limits<PetscInt>::max());
		assert(iga::CheckedImmersedStaticFlowRowCount(1, 2, true) == 7);
		assert(iga::CheckedImmersedStaticFlowRowCount(petsc_max/4, petsc_max%4, false) == petsc_max);
		Reject([&] { (void)iga::CheckedImmersedStaticFlowRowCount(petsc_max/4, petsc_max%4, true); });
		Reject([&] { (void)iga::CheckedImmersedStaticFlowRowCount(std::numeric_limits<std::size_t>::max()/4+1, 0, false); });

		const auto domain = Domain();
		// Port algebra is insensitive to octree refinement; the dedicated aneurysm
		// regressions cover the deeper physical-refinement path.
		const iga::CutCellVolumeQuadratureCatalog expanded(domain, {2,500000,500000,3000000});
		const iga::CutCellVolumeQuadratureCatalog compact(domain, {2,500000,500000,3000000}, iga::CutCellVolumeQuadratureStorageMode::Compact);
		const iga::ImmersedSurfaceQuadratureCatalog surface(domain);
		const iga::CutCellGhostPenaltyCatalog expanded_ghost(domain, expanded);
		const iga::CutCellGhostPenaltyCatalog compact_ghost(domain, compact);
		iga::ImmersedStaticFlowOptions base; base.parameters = {1.0, 1.0, 0.0}; base.wall_labels = {0};

		// Zero is a valid port label when it does not overlap a wall.  Assemble a
		// pressure load and measure the same selected cap through the full runtime.
		const auto zero_domain = ZeroPortDomain();
		const iga::CutCellVolumeQuadratureCatalog zero_volume(zero_domain, {2,500000,500000,3000000});
		const iga::ImmersedSurfaceQuadratureCatalog zero_surface(zero_domain);
		const iga::CutCellGhostPenaltyCatalog zero_ghost(zero_domain, zero_volume);
		iga::ImmersedStaticFlowOptions zero_port_options;
		zero_port_options.parameters = {1.0, 1.0, 0.0}; zero_port_options.wall_labels = {3};
		zero_port_options.ports = {{"zero",0,iga::ImmersedFlowPortControlMode::Pressure,5.0},
			{"other",2,iga::ImmersedFlowPortControlMode::Pressure,7.0}};
		zero_port_options.assemble_volume = false; zero_port_options.assemble_wall = false;
		zero_port_options.assemble_ghost = false;
		iga::ImmersedStaticFlowRuntime zero_port_runtime(zero_domain, zero_volume, zero_surface, zero_ghost, zero_port_options);
		zero_port_runtime.Assemble();
		const auto& zero_port = Port(zero_port_runtime, "zero");
		assert(zero_port.assembled_surface_points > 0 && zero_port.area_m2 > 0.0);
		assert(std::abs(zero_port.measurement.area_m2-zero_port.area_m2) < 2e-12);
		assert(std::abs(zero_port.measurement.mean_pressure_pa) < 2e-12);
		assert(std::any_of(zero_port_runtime.AssembledNegativeResidual().begin(), zero_port_runtime.AssembledNegativeResidual().end(),
			[](PetscScalar value) { return PetscRealPart(value) != 0.0; }));

		// A cap-only cut cell is legal: its open patch has no selected Nitsche
		// points, but it still retains required ghost coverage.
		bool cap_only = false;
		for (std::uint64_t id = 0; id < domain.Cells().size(); ++id) {
			if (domain.Cells()[static_cast<std::size_t>(id)].classification != iga::CellClassification::Cut || !expanded.Cell(id).usable) continue;
			const auto& rule = surface.UsableRule(domain, id);
			cap_only = cap_only || ((HasLabel(rule, 1) || HasLabel(rule, 2)) && !HasLabel(rule, 0));
			if (expanded.Cell(id).diagnostics.estimated_reference_volume > 0.0) assert(expanded_ghost.Covered(id));
		}
		assert(cap_only);

		// No-port construction restores the Phase-5 selected-wall preflight;
		// configured ports retain the cap-only cut-cell exception below.
		Reject([&] { iga::ImmersedStaticFlowRuntime rejected_no_port(domain, expanded, surface, expanded_ghost, base); });

		// No port with all surface labels selected preserves the Phase-5 gauge row.
		iga::ImmersedStaticFlowOptions closed = base; closed.wall_labels = {0,1,2};
		iga::ImmersedStaticFlowRuntime no_port(domain, expanded, surface, expanded_ghost, closed);
		assert(no_port.HasGauge() && no_port.Diagnostics().total_dofs == 4*no_port.Diagnostics().active_nodes+1);
		assert(no_port.ScalarDiagonalStructureVerified());

		iga::ImmersedStaticFlowOptions flow = base;
		flow.ports = {{"in",1,iga::ImmersedFlowPortControlMode::FlowRate,0.0}, {"out",2,iga::ImmersedFlowPortControlMode::FlowRate,0.0}};
		// Isolate the controller block: wall, volume, ghost, and gauge terms do
		// not participate in these residual/Jacobian assertions.
		flow.assemble_volume = false; flow.assemble_wall = false; flow.assemble_ghost = false; flow.assemble_gauge = false;
		iga::ImmersedStaticFlowRuntime flow_runtime(domain, expanded, surface, expanded_ghost, flow);
		assert(flow_runtime.HasGauge());
		assert(flow_runtime.Diagnostics().total_dofs == 4*flow_runtime.Diagnostics().active_nodes+3);
		assert(flow_runtime.ScalarDiagonalStructureVerified());
		const auto scalar_seed = flow_runtime.AssembledJacobianDense();
		const std::size_t scalar_n = flow_runtime.Diagnostics().total_dofs;
		for (const auto* id : {"in", "out"}) {
			const auto row = static_cast<std::size_t>(flow_runtime.PortMultiplierDof(id));
			assert(PetscRealPart(scalar_seed[row*scalar_n+row]) == 0.0);
		}
		const auto gauge_row = static_cast<std::size_t>(flow_runtime.GaugeDof());
		assert(PetscRealPart(scalar_seed[gauge_row*scalar_n+gauge_row]) == 0.0);
		const double area = Port(flow_runtime, "in").area_m2;
		assert(std::abs(area-.64) < 2e-12 && std::abs(Port(flow_runtime, "out").area_m2-area) < 2e-12);
		assert(std::abs(area-surface.Diagnostics().source_area_by_boundary_id.at(1))
			<= 1024.0*std::numeric_limits<double>::epsilon()*area);
		flow_runtime.SetPortControlValue("in", -2.0*area); flow_runtime.SetPortControlValue("out", 2.0*area);
		std::vector<PetscScalar> state(flow_runtime.Diagnostics().total_dofs, 0.0);
		for (std::size_t a = 0; a < flow_runtime.ActiveNodes().size(); ++a) state[4*a+2] = 2.0;
		state[static_cast<std::size_t>(flow_runtime.PortMultiplierDof("in"))] = 3.0;
		flow_runtime.SetCommittedState(state); flow_runtime.Assemble();
		const auto flow_rhs = flow_runtime.AssembledNegativeResidual();
		const auto& in = Port(flow_runtime, "in"); const auto& out = Port(flow_runtime, "out");
		assert(std::abs(in.measurement.outward_flow_m3_s+2.0*area) <= 1e-12*area);
		assert(std::abs(out.measurement.outward_flow_m3_s-2.0*area) <= 1e-12*area);
		assert(std::abs(in.constraint_residual) <= 1e-12*area && std::abs(out.constraint_residual) <= 1e-12*area);
		assert(std::abs(PetscRealPart(flow_rhs[static_cast<std::size_t>(in.multiplier_row)])) <= 1e-12*area);
		// R_u += lambda*c, rhs=-R, and both Jacobian off-diagonal blocks are +c.
		std::vector<PetscScalar> direction(state.size(), 0.0); direction[static_cast<std::size_t>(in.multiplier_row)] = 1.0;
		const auto action = flow_runtime.AssembledJacobianAction(direction);
		bool nonzero_c = false;
		double c_squared = 0.0;
		for (std::size_t row = 0; row < action.size(); ++row) if (row < flow_runtime.Diagnostics().physical_dofs) {
			nonzero_c = nonzero_c || std::abs(PetscRealPart(action[row])) > 0.0;
			assert(std::abs(PetscRealPart(flow_rhs[row]+3.0*action[row])) < 3e-11);
			c_squared += PetscRealPart(action[row])*PetscRealPart(action[row]);
		}
		assert(nonzero_c);
		std::vector<PetscScalar> velocity_direction = action;
		velocity_direction[static_cast<std::size_t>(in.multiplier_row)] = 0.0;
		velocity_direction[static_cast<std::size_t>(out.multiplier_row)] = 0.0;
		velocity_direction[static_cast<std::size_t>(flow_runtime.GaugeDof())] = 0.0;
		const auto row_action = flow_runtime.AssembledJacobianAction(velocity_direction);
		assert(std::abs(PetscRealPart(row_action[static_cast<std::size_t>(in.multiplier_row)])-c_squared) < 3e-11);

		// Constant p and u are analytic polynomial traces: Q, pbar, and tbar_n
		// are evaluated on exactly the retained rule used by the load/controller.
		assert(std::abs(out.measurement.mean_pressure_pa) < 1e-13 && std::abs(out.measurement.mean_normal_traction_pa) < 1e-13);
		for (std::size_t a = 0; a < flow_runtime.ActiveNodes().size(); ++a) state[4*a+3] = 7.0;
		flow_runtime.SetCommittedState(state); flow_runtime.Assemble();
		assert(std::abs(Port(flow_runtime, "out").measurement.mean_pressure_pa-7.0) < 2e-12);
		assert(std::abs(Port(flow_runtime, "out").measurement.mean_normal_traction_pa+7.0) < 2e-12);
		// Affine u_z=4z has symmetric-gradient contribution 2*mu*4*n_z^2.
		// Integrate that analytic traction independently over the retained outlet.
		state.assign(state.size(), 0.0);
		const auto grid = domain.Background().Spec();
		const std::uint32_t nx = grid.cells[0]+3, ny = grid.cells[1]+3;
		for (const auto node : flow_runtime.ActiveNodes()) {
			const std::uint32_t i = static_cast<std::uint32_t>(node)%nx;
			const std::uint32_t j = (static_cast<std::uint32_t>(node)/nx)%ny;
			const std::uint32_t k = static_cast<std::uint32_t>(node)/(nx*ny);
			const auto x = domain.Background().Greville(i, j, k);
			state[static_cast<std::size_t>(flow_runtime.Dof(node, 2))] = 4.0*x[2];
			state[static_cast<std::size_t>(flow_runtime.Dof(node, 3))] = 3.0;
		}
		flow_runtime.SetCommittedState(state); flow_runtime.Assemble();
		double outlet_area = 0.0, analytic_traction_integral = 0.0;
		for (std::uint64_t id = 0; id < domain.Cells().size(); ++id)
			for (const auto& point : surface.UsableRule(domain, id).Points()) if (point.boundary_id == 2) {
				outlet_area += point.weight;
				analytic_traction_integral += (-3.0+8.0*point.normal[2]*point.normal[2])*point.weight;
			}
		assert(std::abs(Port(flow_runtime, "out").measurement.mean_normal_traction_pa-analytic_traction_integral/outlet_area) < 3e-11);

		// A compatible tiny SI target does not accept the zero-flow state.
		iga::ImmersedStaticFlowOptions tiny = flow;
		tiny.ports[0].value = -5e-12; tiny.ports[1].value = 5e-12;
		iga::ImmersedStaticFlowRuntime tiny_runtime(domain, expanded, surface, expanded_ghost, tiny);
		tiny_runtime.Assemble();
		assert(Port(tiny_runtime, "in").absolute_flow_residual_m3_s == 5e-12);
		assert(Port(tiny_runtime, "in").absolute_flow_residual_m3_s > Port(tiny_runtime, "in").flow_tolerance_m3_s);
		assert(std::isfinite(Port(tiny_runtime, "in").normalized_flow_residual) && Port(tiny_runtime, "in").normalized_flow_residual > 1.0);
		iga::ImmersedStaticFlowOptions tiny_solve = base;
		tiny_solve.ports = tiny.ports;
		tiny_solve.lu_pivot_shift = 1e-12;
		tiny_solve.flow_controller_absolute_tolerance_m3_s = 0.0;
		tiny_solve.flow_controller_relative_tolerance = 1e-10;
		tiny_solve.flow_controller_reference_flow_m3_s = 5e-12;
		iga::ImmersedStaticFlowRuntime tiny_solve_runtime(domain, expanded, surface, expanded_ghost, tiny_solve);
		assert(tiny_solve_runtime.SolveTrial());
		for (const auto* id : {"in", "out"}) {
			const auto& port = Port(tiny_solve_runtime, id);
			assert(port.absolute_flow_residual_m3_s <= 1e-10*std::abs(port.target));
			assert(port.absolute_flow_residual_m3_s <= port.flow_tolerance_m3_s && std::isfinite(port.normalized_flow_residual) && port.normalized_flow_residual <= 1.0);
		}
		// A zero absolute tolerance cannot permit a relative tolerance/reference
		// pair that underflows to zero before residual normalization.
		iga::ImmersedStaticFlowOptions underflow = flow;
		underflow.flow_controller_absolute_tolerance_m3_s = 0.0;
		underflow.flow_controller_relative_tolerance = std::numeric_limits<double>::denorm_min();
		underflow.flow_controller_reference_flow_m3_s = std::numeric_limits<double>::denorm_min();
		iga::ImmersedStaticFlowRuntime underflow_runtime(domain, expanded, surface, expanded_ghost, underflow);
		Reject([&] { underflow_runtime.Assemble(); });

		iga::ImmersedStaticFlowOptions pressure = base;
		// Both caps remain open ports; label 0 is the only wall label.
		pressure.ports = {{"pressure1",1,iga::ImmersedFlowPortControlMode::Pressure,5.0}, {"pressure2",2,iga::ImmersedFlowPortControlMode::Pressure,7.0}};
		pressure.assemble_volume = false; pressure.assemble_wall = false; pressure.assemble_ghost = false;
		iga::ImmersedStaticFlowRuntime pressure_runtime(domain, expanded, surface, expanded_ghost, pressure);
		assert(!pressure_runtime.HasGauge() && pressure_runtime.Diagnostics().total_dofs == pressure_runtime.Diagnostics().physical_dofs);
		pressure_runtime.Assemble(); const auto pressure_rhs = pressure_runtime.AssembledNegativeResidual();
		iga::ImmersedStaticFlowOptions traction = base;
		traction.ports = {{"traction1",1,iga::ImmersedFlowPortControlMode::MeanNormalTraction,5.0}, {"traction2",2,iga::ImmersedFlowPortControlMode::MeanNormalTraction,7.0}};
		traction.assemble_volume = false; traction.assemble_wall = false; traction.assemble_ghost = false;
		iga::ImmersedStaticFlowRuntime traction_runtime(domain, expanded, surface, expanded_ghost, traction);
		traction_runtime.Assemble(); const auto traction_rhs = traction_runtime.AssembledNegativeResidual();
		// Independent oracle over *all retained cap rules* (no copied production
		// active-cell predicate).  Each cap label is visited once per retained
		// quadrature point, its eligible area is audited against source area, and
		// the combined body-fitted load is compared below.
		std::vector<PetscScalar> independent_pressure(pressure_rhs.size(), 0.0);
		std::array<double, 3> eligible_area{{0.0,0.0,0.0}};
		std::array<std::size_t, 3> visits{{0,0,0}};
		for (std::uint64_t id = 0; id < domain.Cells().size(); ++id) {
			const auto element = domain.Background().MaterializeElement(id);
			for (const auto& point : surface.UsableRule(domain, id).Points()) {
				if (point.boundary_id != 1 && point.boundary_id != 2) continue;
				const int cap = point.boundary_id; const double value = cap == 1 ? 5.0 : 7.0;
				eligible_area[static_cast<std::size_t>(cap)] += point.weight; ++visits[static_cast<std::size_t>(cap)];
				const auto basis = iga::EvaluateBasis(element, point.parametric[0], point.parametric[1], point.parametric[2], false);
				for (std::size_t a = 0; a < element.connectivity.size(); ++a)
					for (int component = 0; component < 3; ++component)
						independent_pressure[static_cast<std::size_t>(pressure_runtime.Dof(element.connectivity[a], component))]
							-= value*basis.value[a]*point.normal[component]*point.weight;
			}
		}
		assert(visits[1] > 0 && visits[2] > 0);
		assert(Port(pressure_runtime, "pressure1").assembled_surface_points == visits[1]);
		assert(Port(pressure_runtime, "pressure2").assembled_surface_points == visits[2]);
		assert(std::abs(eligible_area[1]-surface.Diagnostics().source_area_by_boundary_id.at(1)) < 2e-12);
		assert(std::abs(eligible_area[2]-surface.Diagnostics().source_area_by_boundary_id.at(2)) < 2e-12);
		bool pressure_sign = false;
		for (std::size_t row = 0; row < pressure_rhs.size(); ++row) {
			pressure_sign = pressure_sign || std::abs(PetscRealPart(pressure_rhs[row])) > 0.0;
			assert(std::abs(PetscRealPart(pressure_rhs[row]-independent_pressure[row])) < 3e-11);
			assert(std::abs(PetscRealPart(pressure_rhs[row]+traction_rhs[row])) < 3e-11);
		}
		assert(pressure_sign); // t_n maps p_b=-t_n and therefore reverses the body-fitted load sign.

		// The conservative resolved mixed form is completed on every retained
		// physical label.  Summing physical pressure residuals therefore yields
		// the open velocity trace plus the prescribed wall-normal trace, not the
		// independently sampled volume-divergence-minus-wall-flow diagnostic.
		iga::ImmersedStaticFlowOptions mixed_trace = base;
		mixed_trace.ports = { {"pressure1",1,iga::ImmersedFlowPortControlMode::Pressure,0.0},
			{"pressure2",2,iga::ImmersedFlowPortControlMode::Pressure,0.0} };
		mixed_trace.assemble_ghost = false;
		mixed_trace.wall_velocity = [](const std::array<double, 3>& x, int) {
			return std::array<double, 3>{{x[0],x[1],0.0}};
		};
		iga::ImmersedStaticFlowRuntime mixed_trace_runtime(domain, expanded, surface, expanded_ghost, mixed_trace);
		std::vector<PetscScalar> mixed_trace_state(mixed_trace_runtime.Diagnostics().total_dofs, 0.0);
		const auto trace_grid = domain.Background().Spec();
		const std::uint32_t trace_nx = trace_grid.cells[0]+3, trace_ny = trace_grid.cells[1]+3;
		for (const auto node : mixed_trace_runtime.ActiveNodes()) {
			const std::uint32_t i = static_cast<std::uint32_t>(node)%trace_nx;
			const std::uint32_t j = (static_cast<std::uint32_t>(node)/trace_nx)%trace_ny;
			const std::uint32_t k = static_cast<std::uint32_t>(node)/(trace_nx*trace_ny);
			const auto x = domain.Background().Greville(i, j, k);
			mixed_trace_state[static_cast<std::size_t>(mixed_trace_runtime.Dof(node, 2))] = 2.0*x[2];
		}
		mixed_trace_runtime.SetCommittedState(mixed_trace_state); mixed_trace_runtime.Assemble();
		const auto mixed_trace_rhs = mixed_trace_runtime.AssembledNegativeResidual();
		double pressure_residual_sum = 0.0;
		for (std::size_t node = 0; node < mixed_trace_runtime.Diagnostics().active_nodes; ++node)
			pressure_residual_sum -= PetscRealPart(mixed_trace_rhs[4*node+3]); // Convert stored -R to R.
		const auto mixed_trace_conservation = mixed_trace_runtime.ConservationDiagnostics();
		double label_sum = 0.0, open_label_sum = 0.0, wall_label_sum = 0.0, prescribed_wall_trace = 0.0;
		for (const auto& entry : mixed_trace_conservation.surface_flow_by_boundary_label_m3_s) {
			label_sum += entry.second;
			if (entry.first == 0) wall_label_sum += entry.second;
			else open_label_sum += entry.second;
		}
		for (std::uint64_t id = 0; id < domain.Cells().size(); ++id) {
			const auto element = domain.Background().MaterializeElement(id);
			for (const auto& point : surface.UsableRule(domain, id).Points()) if (point.boundary_id == 0)
				prescribed_wall_trace += (point.physical[0]*point.normal[0]+point.physical[1]*point.normal[1])*point.weight;
		}
		assert(std::abs(label_sum-mixed_trace_conservation.total_surface_outward_flow_m3_s) < 2e-12);
		assert(std::abs(open_label_sum-mixed_trace_conservation.open_port_outward_flow_m3_s) < 2e-12);
		assert(std::abs(wall_label_sum-mixed_trace_conservation.wall_outward_flow_m3_s) < 2e-12);
		assert(std::abs(pressure_residual_sum-(mixed_trace_conservation.open_port_outward_flow_m3_s+prescribed_wall_trace)) < 3e-11);
		const double old_volume_minus_wall = mixed_trace_conservation.volume_divergence_integral_m3_s
			-mixed_trace_conservation.wall_outward_flow_m3_s;
		assert(std::abs(pressure_residual_sum-old_volume_minus_wall) > 1e-4);

		iga::ImmersedStaticFlowOptions mixed = flow; mixed.ports[1] = {"pressure",2,iga::ImmersedFlowPortControlMode::Pressure,0.0};
		iga::ImmersedStaticFlowRuntime mixed_runtime(domain, expanded, surface, expanded_ghost, mixed);
		assert(!mixed_runtime.HasGauge() && mixed_runtime.Diagnostics().total_dofs == mixed_runtime.Diagnostics().physical_dofs+1);
		iga::ImmersedStaticFlowOptions two_pressure = base; two_pressure.ports = {{"pressure1",1,iga::ImmersedFlowPortControlMode::Pressure,0.0}, {"pressure2",2,iga::ImmersedFlowPortControlMode::Pressure,0.0}};
		iga::ImmersedStaticFlowRuntime two_pressure_runtime(domain, expanded, surface, expanded_ghost, two_pressure);
		assert(!two_pressure_runtime.HasGauge());

		Reject([&] { auto bad = flow; bad.ports[0].value = 1.0; iga::ImmersedStaticFlowRuntime x(domain, expanded, surface, expanded_ghost, bad); });
		Reject([&] { auto bad = pressure; bad.ports.push_back({"duplicate",2,iga::ImmersedFlowPortControlMode::Pressure,0.0}); iga::ImmersedStaticFlowRuntime x(domain, expanded, surface, expanded_ghost, bad); });
		Reject([&] { auto bad = pressure; bad.ports[0].boundary_label = 9; iga::ImmersedStaticFlowRuntime x(domain, expanded, surface, expanded_ghost, bad); });
		Reject([&] { auto bad = pressure; bad.wall_labels.clear(); iga::ImmersedStaticFlowRuntime x(domain, expanded, surface, expanded_ghost, bad); });
		Reject([&] { auto bad = pressure; bad.ports[0].control_mode = iga::ImmersedFlowPortControlMode::TotalPressure; iga::ImmersedStaticFlowRuntime x(domain, expanded, surface, expanded_ghost, bad); });
		Reject([&] { auto bad = pressure; bad.ports[0].boundary_label = -1; iga::ImmersedStaticFlowRuntime x(domain, expanded, surface, expanded_ghost, bad); });
		// Port arithmetic fails deterministically before any non-finite PETSc insertion.
		Reject([&] { (void)iga::CheckedImmersedFlowPortProduct(std::numeric_limits<double>::max(), 2.0, "test"); });
		for (std::uint64_t id = 0; id < domain.Cells().size(); ++id) {
			const auto& rule = surface.UsableRule(domain, id);
			const auto found = std::find_if(rule.Points().begin(), rule.Points().end(), [](const auto& point) { return point.boundary_id == 1; });
			if (found == rule.Points().end()) continue;
			auto absurd = *found; absurd.weight = std::numeric_limits<double>::max();
			const iga::SurfaceQuadratureRule overflow_rule({absurd});
			const auto element = domain.Background().MaterializeElement(id);
			std::vector<std::array<double, 4>> zero_nodal(element.connectivity.size());
			Reject([&] { (void)iga::BuildImmersedFlowPortElement(element, rule, -1, iga::ImmersedFlowPortControlMode::Pressure, 0.0, zero_nodal); });
			Reject([&] { (void)iga::MeasureImmersedFlowPortElement(element, rule, -1, zero_nodal, 1.0); });
			Reject([&] { (void)iga::BuildImmersedFlowPortElement(element, overflow_rule, 1, iga::ImmersedFlowPortControlMode::Pressure, std::numeric_limits<double>::max(), zero_nodal); });
			break;
		}

		// Compact rules keep the same fixed port topology, retained areas, and
		// zero-state load result without expanding logical volume points.
		iga::ImmersedStaticFlowRuntime compact_runtime(domain, compact, surface, compact_ghost, pressure);
		compact_runtime.Assemble();
		for (const auto* id : {"pressure1", "pressure2"}) {
			const auto& expanded_port = Port(pressure_runtime, id);
			const auto& compact_port = Port(compact_runtime, id);
			assert(std::abs(compact_port.area_m2-expanded_port.area_m2) < 2e-12);
			assert(compact_port.assembled_surface_points == expanded_port.assembled_surface_points);
			assert(std::abs(compact_port.measurement.area_m2-expanded_port.measurement.area_m2) < 2e-12);
			assert(std::abs(compact_port.measurement.outward_flow_m3_s-expanded_port.measurement.outward_flow_m3_s) < 2e-12);
			assert(std::abs(compact_port.measurement.mean_pressure_pa-expanded_port.measurement.mean_pressure_pa) < 2e-12);
			assert(std::abs(compact_port.measurement.mean_normal_traction_pa-expanded_port.measurement.mean_normal_traction_pa) < 2e-12);
			assert(std::abs(compact_port.measurement.mean_velocity_squared_m2_s2-expanded_port.measurement.mean_velocity_squared_m2_s2) < 2e-12);
		}
		assert(compact_runtime.AssembledNegativeResidual() == pressure_rhs);
	} catch (const std::exception& error) {
		std::cerr << "immersed_flow_port_test: " << error.what() << '\n'; status = 1;
	}
	PetscFinalize();
	return status;
}
