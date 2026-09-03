#include "ImmersedStaticFlowRuntime.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

iga::RawSurfaceTriangle Face(std::int64_t a, std::int64_t b, std::int64_t c)
{
	iga::RawSurfaceTriangle result; result.indices = {{a,b,c}}; result.boundary_id = 7; return result;
}
iga::RawSurfaceSoup Cube(double lower, double upper)
{
	iga::RawSurfaceSoup result;
	result.vertices = {{{{lower,lower,lower}},{{upper,lower,lower}},{{upper,upper,lower}},{{lower,upper,lower}},
		{{lower,lower,upper}},{{upper,lower,upper}},{{upper,upper,upper}},{{lower,upper,upper}}}};
	result.triangles = {Face(0,2,1),Face(0,3,2),Face(4,5,6),Face(4,6,7),Face(0,1,5),Face(0,5,4),
		Face(1,2,6),Face(1,6,5),Face(2,3,7),Face(2,7,6),Face(3,0,4),Face(3,4,7)};
	return result;
}
iga::CartesianDomainClassification Domain()
{
	const iga::CubicCartesianGridSpec spec{{{0,0,0}},{{1,1,1}},{{3,3,3}}};
	return iga::CartesianDomainClassification(iga::CubicCartesianBackground(spec),
		iga::SurfaceSpatialIndex(iga::ClosedTriangulatedSurface::Build(Cube(.1,.9))));
}
struct Depth2CompactFixture {
	iga::CartesianDomainClassification domain;
	iga::CutCellVolumeQuadratureCatalog volume;
	iga::ImmersedSurfaceQuadratureCatalog surface;
	iga::CutCellGhostPenaltyCatalog ghost;

	Depth2CompactFixture()
		: domain(Domain()), volume(domain, {2,500000,500000,3000000}, iga::CutCellVolumeQuadratureStorageMode::Compact),
		surface(domain), ghost(domain, volume) {}
};
void AssertSameTopology(const iga::CartesianDomainClassification& expected,
	const iga::CartesianDomainClassification& actual)
{
	assert(expected.SurfaceCanonicalHash() == actual.SurfaceCanonicalHash());
	assert(expected.Background().Spec().lower_m == actual.Background().Spec().lower_m);
	assert(expected.Background().Spec().upper_m == actual.Background().Spec().upper_m);
	assert(expected.Background().Spec().cells == actual.Background().Spec().cells);
	assert(expected.Cells().size() == actual.Cells().size());
	for (std::size_t i = 0; i < expected.Cells().size(); ++i) {
		const auto& left = expected.Cells()[i]; const auto& right = actual.Cells()[i];
		assert(left.id == right.id && left.index == right.index && left.classification == right.classification
			&& left.triangle_ids == right.triangle_ids && left.boundary_only_contact == right.boundary_only_contact
			&& left.ambiguous == right.ambiguous);
	}
}
template <class Function> void Reject(Function&& function)
{
	bool rejected = false; try { function(); } catch (const std::exception&) { rejected = true; } assert(rejected);
}
double DirectMeasure(const iga::CartesianDomainClassification& domain, const iga::CutCellVolumeQuadratureCatalog& volume)
{
	double total = 0.0;
	for (std::uint64_t id = 0; id < domain.Cells().size(); ++id) {
		const auto& cell = volume.Cell(id);
		const auto classification = domain.Cells()[static_cast<std::size_t>(id)].classification;
		if ((classification != iga::CellClassification::Inside && classification != iga::CellClassification::Cut)
			|| !cell.usable || !(cell.diagnostics.estimated_reference_volume > 0.0)) continue;
		const auto element = domain.Background().MaterializeElement(id);
		const auto add = [&](const iga::VolumeQuadraturePoint& point) { total += point.weight*iga::EvaluateBasis(element, point.parametric[0], point.parametric[1], point.parametric[2]).raw_determinant; };
		if (volume.StorageMode() == iga::CutCellVolumeQuadratureStorageMode::Expanded)
			for (const auto& point : volume.UsableRule(domain, id).Points()) add(point);
		else iga::ForEachVolumePoint(volume.UsableCompactRule(domain, id), add);
	}
	return total;
}

double VectorNorm(const std::vector<PetscScalar>& values)
{
	double squared = 0.0;
	for (const auto value : values) squared += PetscRealPart(value)*PetscRealPart(value);
	return std::sqrt(squared);
}

std::array<double, 3> ResidualBlockNorms(const iga::ImmersedStaticFlowRuntime& runtime,
	const std::vector<PetscScalar>& residual)
{
	assert(residual.size() == runtime.Diagnostics().total_dofs);
	std::array<double, 3> result{{0.0,0.0,0.0}};
	for (std::size_t row = 0; row < residual.size(); ++row) {
		const std::size_t block = row >= runtime.Diagnostics().physical_dofs ? 2 : ((row%4) == 3 ? 1 : 0);
		result[block] += PetscRealPart(residual[row])*PetscRealPart(residual[row]);
	}
	for (auto& value : result) value = std::sqrt(value);
	return result;
}

double HydrostaticPressure(const std::array<double, 3>& x)
{
	return .3*x[0]-.2*x[1]+.1*x[2]-.1;
}

std::vector<PetscScalar> HydrostaticState(const iga::CartesianDomainClassification& domain,
	const iga::ImmersedStaticFlowRuntime& runtime)
{
	std::vector<PetscScalar> result(runtime.Diagnostics().total_dofs, 0.0);
	const auto& background = domain.Background();
	const auto& spec = background.Spec();
	const std::uint64_t nx = static_cast<std::uint64_t>(spec.cells[0])+3;
	const std::uint64_t ny = static_cast<std::uint64_t>(spec.cells[1])+3;
	for (const auto node : runtime.ActiveNodes()) {
		const auto raw = static_cast<std::uint64_t>(node);
		const std::uint32_t i = static_cast<std::uint32_t>(raw%nx);
		const std::uint32_t j = static_cast<std::uint32_t>((raw/nx)%ny);
		const std::uint32_t k = static_cast<std::uint32_t>(raw/(nx*ny));
		result[static_cast<std::size_t>(runtime.Dof(node, 3))] = HydrostaticPressure(background.Greville(i, j, k));
	}
	return result;
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int status = 0;
	try {
	{
	const auto domain = Domain();
	const iga::CutCellVolumeQuadratureCatalog volume(domain, {4,500000,500000,3000000});
	const iga::ImmersedSurfaceQuadratureCatalog surface(domain);
	const iga::CutCellGhostPenaltyCatalog ghost(domain, volume);
	iga::ImmersedStaticFlowOptions options; options.parameters = {1.0, 1.0, 0.0}; options.wall_labels = {7};
	iga::ImmersedStaticFlowRuntime zero(domain, volume, surface, ghost, options);
	zero.Assemble(); const auto first = zero.AssembledNegativeResidual(); zero.Assemble();
	assert(first == zero.AssembledNegativeResidual());
	// A pivot shift is strictly an LU-preconditioner option: it must not alter
	// the assembled operator or residual presented to GMRES.
	iga::ImmersedStaticFlowOptions shifted_assembly_options = options;
	shifted_assembly_options.lu_pivot_shift = 1e-12;
	iga::ImmersedStaticFlowRuntime shifted_assembly(domain, volume, surface, ghost, shifted_assembly_options);
	shifted_assembly.Assemble();
	assert(first == shifted_assembly.AssembledNegativeResidual());
	assert(zero.AssembledJacobianDense() == shifted_assembly.AssembledJacobianDense());
	std::vector<std::int32_t> expected_active;
	for (std::uint64_t id = 0; id < domain.Cells().size(); ++id) {
		const auto classification = domain.Cells()[static_cast<std::size_t>(id)].classification;
		const auto& cell = volume.Cell(id);
		if ((classification != iga::CellClassification::Inside && classification != iga::CellClassification::Cut) || !cell.usable
			|| !(cell.diagnostics.estimated_reference_volume > 0.0)) continue;
		const auto element = domain.Background().MaterializeElement(id);
		expected_active.insert(expected_active.end(), element.connectivity.begin(), element.connectivity.end());
	}
	std::sort(expected_active.begin(), expected_active.end()); expected_active.erase(std::unique(expected_active.begin(), expected_active.end()), expected_active.end());
	assert(expected_active == zero.ActiveNodes());
	assert(zero.Diagnostics().total_dofs == 4*zero.Diagnostics().active_nodes+1);
	assert(zero.Diagnostics().volume_cells > 0 && zero.Diagnostics().surface_cells > 0 && zero.Diagnostics().ghost_faces == ghost.Faces().size());
	assert(std::abs(zero.Diagnostics().pressure_measure-DirectMeasure(domain, volume)) < 2e-12);
	assert(std::isfinite(zero.Diagnostics().constant_pressure_defect));
	for (const auto node : zero.ActiveNodes()) assert(node >= 0);
	for (std::uint64_t raw_node = 0; raw_node < domain.Background().NodeCount(); ++raw_node)
		if (!std::binary_search(zero.ActiveNodes().begin(), zero.ActiveNodes().end(), static_cast<std::int32_t>(raw_node))) {
			Reject([&] { zero.Dof(static_cast<std::int32_t>(raw_node), 0); }); break;
		}
	Reject([&] { zero.Dof(-1, 0); });
	// R=[F+g lambda;g^T p], while rhs stores -R.  Test both signs through
	// independent lambda and pressure states, plus the assembled matrix action.
	std::vector<PetscScalar> gauge_state(zero.Diagnostics().total_dofs, 0.0);
	gauge_state[static_cast<std::size_t>(zero.GaugeDof())] = 2.0;
	zero.SetCommittedState(gauge_state); zero.Assemble(); const auto lambda_rhs = zero.AssembledNegativeResidual();
	std::vector<PetscScalar> lambda_direction(gauge_state.size(), 0.0); lambda_direction[static_cast<std::size_t>(zero.GaugeDof())] = 1.0;
	const auto lambda_action = zero.AssembledJacobianAction(lambda_direction);
	for (std::size_t a = 0; a < zero.ActiveNodes().size(); ++a) {
		const auto p = 4*a+3; assert(std::abs(PetscRealPart(lambda_rhs[p]+2.0*lambda_action[p])) < 2e-11);
	}
	gauge_state.assign(gauge_state.size(), 0.0); for (std::size_t a = 0; a < zero.ActiveNodes().size(); ++a) gauge_state[4*a+3] = 1.0;
	zero.SetCommittedState(gauge_state); zero.Assemble();
	assert(std::abs(PetscRealPart(zero.AssembledNegativeResidual()[static_cast<std::size_t>(zero.GaugeDof())])+zero.Diagnostics().pressure_measure) < 2e-11);

	const std::array<double, 3> force_density{{.3,-.2,.1}};
	options.body_force = [force_density](const std::array<double, 3>&) { return force_density; };
	iga::ImmersedStaticFlowOptions force_options = options;
	force_options.parameters = {2.5, 1.0, 0.0}; // Non-unit rho makes the PSPG rho scaling observable.
	force_options.body_force = [force_density](const std::array<double, 3>&) { return force_density; };
	iga::ImmersedStaticFlowOptions no_force_options = force_options;
	no_force_options.body_force = [](const std::array<double, 3>&) { return std::array<double, 3>{{0.0,0.0,0.0}}; };
	iga::ImmersedStaticFlowRuntime force_free(domain, volume, surface, ghost, no_force_options); force_free.Assemble();
	iga::ImmersedStaticFlowRuntime forced(domain, volume, surface, ghost, force_options); forced.Assemble();
	bool differs = false, pressure_differs = false;
	const auto force_free_rhs = force_free.AssembledNegativeResidual(); const auto forced_rhs = forced.AssembledNegativeResidual();
	for (std::size_t i = 0; i < first.size(); ++i) differs = differs || std::abs(PetscRealPart(forced_rhs[i]-force_free_rhs[i])) > 1e-12;
	for (std::size_t a = 0; a < forced.ActiveNodes().size(); ++a) pressure_differs = pressure_differs || std::abs(PetscRealPart(forced_rhs[4*a+3]-force_free_rhs[4*a+3])) > 1e-12;
	assert(differs); // weak force and its strong VMS residual path both enter the element operator.
	assert(pressure_differs); // PSPG responds in an explicit pressure row.
	// At zero velocity and pressure, fine_velocity=tau_m*f/rho, so the force
	// difference in a pressure rhs entry (-R) is
	//   + weight*|J|*tau_m/rho * grad(N_a).f .
	// Evaluate that formula directly from each cut-cell quadrature point; this
	// deliberately does not reuse the element assembly implementation.
	std::vector<PetscScalar> expected_force_difference(forced.Diagnostics().total_dofs, 0.0);
	for (std::uint64_t id = 0; id < domain.Cells().size(); ++id) {
		const auto& cell = volume.Cell(id);
		const auto classification = domain.Cells()[static_cast<std::size_t>(id)].classification;
		if ((classification != iga::CellClassification::Inside && classification != iga::CellClassification::Cut)
			|| !cell.usable || !(cell.diagnostics.estimated_reference_volume > 0.0)) continue;
		const auto element = domain.Background().MaterializeElement(id);
		for (const auto& point : volume.UsableRule(domain, id).Points()) {
			const auto basis = iga::EvaluateBasis(element, point.parametric[0], point.parametric[1], point.parametric[2], true);
			double metric[3][3]{};
			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j)
					for (int k = 0; k < 3; ++k) metric[i][j] += basis.inverse_jacobian[k][i]*basis.inverse_jacobian[k][j];
			double metric_norm = 0.0;
			for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) metric_norm += metric[i][j]*metric[i][j];
			const double nu = force_options.parameters.dynamic_viscosity/force_options.parameters.density;
			const double tau_m = 1.0/std::sqrt((nu*nu/12.0)*metric_norm);
			const double factor = point.weight*basis.raw_determinant*tau_m/force_options.parameters.density;
			for (std::size_t a = 0; a < element.connectivity.size(); ++a) {
				double gradient_force = 0.0;
				for (int direction = 0; direction < 3; ++direction) gradient_force += basis.gradient[a][direction]*force_density[direction];
				expected_force_difference[static_cast<std::size_t>(forced.Dof(element.connectivity[a], 3))] += factor*gradient_force;
			}
		}
	}
	for (std::size_t a = 0; a < forced.ActiveNodes().size(); ++a) {
		const double actual = PetscRealPart(forced_rhs[4*a+3]-force_free_rhs[4*a+3]);
		const double expected = PetscRealPart(expected_force_difference[4*a+3]);
		assert(std::abs(actual-expected) < 1e-11*std::max({1.0, std::abs(actual), std::abs(expected)}));
	}
	const auto other_domain = Domain();
	const iga::CutCellVolumeQuadratureCatalog other_volume(other_domain, {4,500000,500000,3000000});
	const iga::CutCellGhostPenaltyCatalog other_ghost(other_domain, other_volume);
	Reject([&] { iga::ImmersedStaticFlowRuntime invalid(domain, volume, surface, other_ghost, options); });
	// Fault injection is performed after the ghost catalog is bound.  Runtime
	// construction must reject the classified cell before it can allocate or
	// mutate a PETSc trial/assembly state.
	iga::CutCellVolumeQuadratureCatalog fault_volume(domain, {4,500000,500000,3000000});
	const iga::CutCellGhostPenaltyCatalog fault_ghost(domain, fault_volume);
	auto& fault_cells = const_cast<std::vector<iga::CutCellVolumeQuadratureCell>&>(fault_volume.Cells());
	auto fault = std::find_if(fault_cells.begin(), fault_cells.end(), [](const iga::CutCellVolumeQuadratureCell& cell) {
		return cell.classification == iga::CellClassification::Inside || cell.classification == iga::CellClassification::Cut;
	});
	assert(fault != fault_cells.end()); fault->usable = false;
	Reject([&] { iga::ImmersedStaticFlowRuntime invalid(domain, fault_volume, surface, fault_ghost, options); });

	// Keep depth-4 compact-rule assembly and measure coverage, but isolate all
	// repeated nonlinear work below.  This verifies compact storage against the
	// production-depth catalog without turning its expanded quadrature into a
	// line-search multiplier.
	const iga::CutCellVolumeQuadratureCatalog compact(domain, {4,500000,500000,3000000}, iga::CutCellVolumeQuadratureStorageMode::Compact);
	const iga::CutCellGhostPenaltyCatalog compact_ghost(domain, compact);
	iga::ImmersedStaticFlowRuntime compact_runtime(domain, compact, surface, compact_ghost, options);
	compact_runtime.Assemble();
	assert(std::abs(compact_runtime.Diagnostics().pressure_measure-DirectMeasure(domain, compact)) < 2e-12);

	// The hydrostatic root and every nonzero closed-wall solve use an otherwise
	// identical independent depth-2 compact fixture.  Keep the physical
	// topology, wall label, body force, and nonlinear options explicit here so
	// this is a runtime-cost isolation rather than a weaker regression.
	Depth2CompactFixture depth2;
	AssertSameTopology(domain, depth2.domain);
	assert(depth2.volume.StorageMode() == iga::CutCellVolumeQuadratureStorageMode::Compact);
	assert(depth2.volume.Options().max_depth == 2);
	assert(depth2.surface.Cells().size() == surface.Cells().size());
	assert(depth2.ghost.Faces().size() == ghost.Faces().size());
	iga::ImmersedStaticFlowOptions solve_options = options; solve_options.lu_pivot_shift = 1e-12;
	assert(solve_options.wall_labels == options.wall_labels && solve_options.parameters.density == options.parameters.density
		&& solve_options.parameters.dynamic_viscosity == options.parameters.dynamic_viscosity && solve_options.parameters.dt == options.parameters.dt);
	assert(solve_options.body_force({{.25,.5,.75}}) == force_density);

	// The conservative resolved mixed form, closed-wall Nitsche topology, and
	// pressure gauge must retain the exactly representable hydrostatic root.
	// The pressure has zero mean on this centered cube, so no gauge shift is
	// needed: u=0 and grad(p)=(.3,-.2,.1) balances the supplied body force.
	iga::ImmersedStaticFlowOptions hydrostatic_options = solve_options;
	hydrostatic_options.parameters = {1.0, 1.0, 0.0};
	hydrostatic_options.body_force = [](const std::array<double, 3>&) {
		return std::array<double, 3>{{.3,-.2,.1}};
	};
	assert(hydrostatic_options.wall_labels == options.wall_labels && hydrostatic_options.body_force({{.25,.5,.75}}) == force_density);
	iga::ImmersedStaticFlowRuntime hydrostatic(depth2.domain, depth2.volume, depth2.surface, depth2.ghost, hydrostatic_options);
	assert(hydrostatic.HasGauge() && hydrostatic.PortDefinitions().empty() && hydrostatic.ActiveNodes() == zero.ActiveNodes());
	hydrostatic.SetCommittedState(HydrostaticState(depth2.domain, hydrostatic));
	hydrostatic.Assemble();
	const auto hydrostatic_residual = hydrostatic.AssembledNegativeResidual();
	assert(VectorNorm(hydrostatic_residual) < 1e-14);
	for (const auto block : ResidualBlockNorms(hydrostatic, hydrostatic_residual)) assert(block < 1e-14);

	iga::ImmersedStaticFlowRuntime solve(depth2.domain, depth2.volume, depth2.surface, depth2.ghost, solve_options);
	assert(solve.HasGauge() && solve.PortDefinitions().empty() && solve.ActiveNodes() == zero.ActiveNodes());
	solve.Assemble(); const auto initial_solve_residual = solve.AssembledNegativeResidual();
	const double initial_solve_norm = VectorNorm(initial_solve_residual);
	const double solve_tolerance = std::max(solve_options.nonlinear_absolute_tolerance,
		solve_options.nonlinear_relative_tolerance*initial_solve_norm);
	assert(solve.SolveTrial());
	const auto final_solve_residual = solve.AssembledNegativeResidual();
	const auto final_blocks = ResidualBlockNorms(solve, final_solve_residual);
	assert(solve.Diagnostics().converged && solve.Diagnostics().nonlinear_iterations >= 2
		&& solve.Diagnostics().newton_steps.size() >= 2);
	assert(std::any_of(solve.Diagnostics().newton_steps.begin(), solve.Diagnostics().newton_steps.end(),
		[](const iga::ImmersedStaticFlowNewtonStep& step) {
			return step.damping == 1.0 && step.candidate_residual_norm < step.residual_norm;
		}));
	assert(solve.Diagnostics().residual_norm <= solve_tolerance && VectorNorm(final_solve_residual) <= solve_tolerance);
	for (const auto block : final_blocks) assert(block <= solve_tolerance);
	Reject([&] { solve.SetCommittedState(solve.CommittedState()); }); solve.Commit();
	for (const auto value : solve.CommittedState()) assert(std::isfinite(PetscRealPart(value)));
	assert(solve.Diagnostics().commit_count == 1 && solve.Diagnostics().committed);
	Reject([&] { solve.Commit(); });
	// Count a whole assembly first, then inject at the first body-force point of
	// the candidate reassembly: one initial assembly worth of calls lies between
	// the current count and the injected call.  A matching control trial proves
	// that this nonzero candidate update is accepted when assembly is allowed.
	std::size_t callback_calls = 0, calls_per_assembly = 0, throw_on_call = 0;
	bool candidate_state_was_mutated = false;
	std::vector<PetscScalar> accepted_candidate;
	iga::ImmersedStaticFlowRuntime* accepted_pointer = nullptr;
	std::size_t capture_candidate_call = 0;
	iga::ImmersedStaticFlowRuntime* reusable_pointer = nullptr;
	iga::ImmersedStaticFlowOptions reusable_options = solve_options;
	reusable_options.body_force = [&](const std::array<double, 3>&) {
		++callback_calls;
		if (accepted_pointer != nullptr && callback_calls == capture_candidate_call)
			accepted_candidate = accepted_pointer->TrialState();
		if (throw_on_call != 0 && callback_calls == throw_on_call) {
			assert(reusable_pointer != nullptr);
			candidate_state_was_mutated = reusable_pointer->TrialState() == accepted_candidate;
			throw std::runtime_error("injected candidate body-force callback failure");
		}
		return force_density;
	};
	assert(reusable_options.wall_labels == solve_options.wall_labels);
	iga::ImmersedStaticFlowRuntime accepted_control(depth2.domain, depth2.volume, depth2.surface, depth2.ghost, reusable_options);
	accepted_pointer = &accepted_control; const std::size_t calls_before_control_measure = callback_calls;
	accepted_control.Assemble(); const std::size_t control_calls_per_assembly = callback_calls-calls_before_control_measure;
	assert(control_calls_per_assembly > 0); capture_candidate_call = callback_calls+control_calls_per_assembly+1;
	assert(accepted_control.SolveTrial()); accepted_pointer = nullptr;
	assert(!accepted_candidate.empty());
	assert(!accepted_control.Diagnostics().newton_steps.empty());
	const auto accepted_first_step = accepted_control.Diagnostics().newton_steps.front();
	assert(accepted_candidate != accepted_control.CommittedState()); accepted_control.Rollback();
	iga::ImmersedStaticFlowRuntime reusable(depth2.domain, depth2.volume, depth2.surface, depth2.ghost, reusable_options);
	reusable_pointer = &reusable; const std::size_t calls_before_measure = callback_calls;
	reusable.Assemble(); calls_per_assembly = callback_calls-calls_before_measure;
	assert(calls_per_assembly > 0);
	const auto reusable_before = reusable.CommittedState();
	throw_on_call = callback_calls+calls_per_assembly+1;
	Reject([&] { reusable.SolveTrial(); });
	assert(candidate_state_was_mutated);
	assert(reusable.CommittedState() == reusable_before && reusable.TrialState() == reusable_before);
	assert(reusable.Diagnostics().rollback_count == 1 && reusable.Diagnostics().ksp_iterations == accepted_first_step.ksp_iterations
		&& reusable.Diagnostics().ksp_reason == accepted_first_step.ksp_reason && reusable.Diagnostics().damping == 0.0);
	throw_on_call = 0; reusable.Assemble(); // Reuse after the throwing callback proves PETSc arrays were restored.
	assert(reusable.SolveTrial()); reusable.Commit();
	assert(reusable.Diagnostics().commit_count == 1);
	}
	} catch (const std::exception& error) {
		// The inner scope has already destroyed every PETSc object.  This keeps a
		// failed solve from terminating the process before PetscFinalize().
		std::cerr << "immersed_static_flow_test: " << error.what() << '\n';
		status = 1;
	}
	PetscFinalize();
	return status;
}
