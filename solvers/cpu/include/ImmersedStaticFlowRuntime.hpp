#ifndef IGA_IMMERSED_STATIC_FLOW_RUNTIME_HPP
#define IGA_IMMERSED_STATIC_FLOW_RUNTIME_HPP

// Serial-only Phase 5 global assembly.  This is intentionally a small
// orchestration layer: geometry catalogs retain rules, while this class owns
// active-row compaction, PETSc insertion, the pressure gauge, and the
// transactional nonlinear trial state.
#include "CutCellGhostPenalty.hpp"
#include "ImmersedFlowPort.hpp"
#include "ImmersedNitscheWall.hpp"

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

struct ImmersedStaticFlowOptions {
	NavierStokesParameters parameters{1.0, 1.0, 0.0};
	std::vector<int> wall_labels;
	// Definitions fix scalar-row topology.  Values are finite constants for a
	// static solve and may be refreshed through SetPortControlValue().
	std::vector<ImmersedFlowPortDefinition> ports;
	double wall_gamma0 = 2.0;
	// The supplied value is evaluated only at selected immersed-wall points.
	// Keeping it here makes the global runtime able to represent nonzero exact
	// wall data without changing the catalog-bound Nitsche formulation.
	ImmersedWallVelocityEvaluator wall_velocity = [](const std::array<double, 3>&, int) {
		return std::array<double, 3>{{0.0,0.0,0.0}};
	};
	NavierStokesBodyForceEvaluator body_force = [](const std::array<double, 3>&) {
		return std::array<double, 3>{{0.0,0.0,0.0}};
	};
	PetscInt nonlinear_maximum_iterations = 12;
	PetscInt ksp_maximum_iterations = 2000;
	double ksp_relative_tolerance = 1e-10;
	double nonlinear_relative_tolerance = 1e-9;
	double nonlinear_absolute_tolerance = 1e-11;
	// SI flow-controller convergence gate.  The floor prevents a tiny target
	// from being silently accepted as zero by a dimensionless residual test.
	double flow_controller_relative_tolerance = 1e-10;
	double flow_controller_absolute_tolerance_m3_s = 1e-15;
	double flow_controller_reference_flow_m3_s = 1e-12;
	// Optional closure-test guard: after the ordinary global rtol/atol is met,
	// continue until every initially nonzero field block has this reduction.
	// Zero-initial gauge blocks use the scale-aware global stopping threshold.
	double nonlinear_block_reduction = 0.0;
	double minimum_damping = 1.0/128.0;
	// These switches are diagnostic-only.  They retain the same active map and
	// gauge, while allowing a closure test to isolate assembled contributions.
	bool assemble_volume = true;
	bool assemble_wall = true;
	bool assemble_ghost = true;
	bool assemble_gauge = true;
	// Applied only inside the LU preconditioner (zero leaves PETSc's ordinary
	// factorization policy in effect).  The assembled operator and all
	// reported residual/gauge checks remain exactly unshifted.
	double lu_pivot_shift = 0.0;
};

struct ImmersedStaticFlowNewtonStep {
	PetscInt iteration = 0, ksp_iterations = 0;
	KSPConvergedReason ksp_reason = KSP_CONVERGED_ITERATING;
	double residual_norm = 0.0, update_norm = 0.0, ksp_residual_norm = 0.0, linear_residual_norm = 0.0, linear_relative_residual = 0.0;
	double damping = 0.0, candidate_residual_norm = 0.0;
};

struct ImmersedStaticFlowDiagnostics {
	std::size_t active_nodes = 0, physical_dofs = 0, total_dofs = 0;
	std::size_t volume_cells = 0, surface_cells = 0, ghost_faces = 0;
	double pressure_measure = 0.0, constant_pressure_defect = 0.0;
	PetscInt nonlinear_iterations = 0, ksp_iterations = 0;
	KSPConvergedReason ksp_reason = KSP_CONVERGED_ITERATING;
	double residual_norm = 0.0, damping = 0.0;
	bool committed = true, trial_active = false, converged = false;
	bool prepared = false;
	std::size_t commit_count = 0, rollback_count = 0, prepare_count = 0,
		finalize_count = 0;
	struct Port {
		std::string id;
		int boundary_label = -1;
		ImmersedFlowPortControlMode control_mode = ImmersedFlowPortControlMode::Pressure;
		double area_m2 = 0.0, target = 0.0, multiplier = 0.0;
		PetscInt multiplier_row = -1;
		ImmersedFlowPortMeasurement measurement{};
		double constraint_residual = 0.0;
		double absolute_flow_residual_m3_s = 0.0;
		double normalized_flow_residual = 0.0;
		double flow_tolerance_m3_s = 0.0;
		std::size_t assembled_surface_points = 0;
	};
	std::vector<Port> ports;
	bool gauge_present = true;
	PetscInt gauge_row = -1;
	// Every appended controller/gauge row has an explicit zero diagonal in the
	// AIJ sparsity pattern.  This is structural only, not a regularization.
	bool scalar_diagonal_structure_verified = false;
	std::size_t wall_selected_points = 0;
	std::map<int, std::size_t> wall_selected_points_by_label;
	std::vector<ImmersedStaticFlowNewtonStep> newton_steps;
};

// Keep the scalar topology in size_t until every contribution has been
// checked.  This lets construction reject an unrepresentable PETSc layout
// before assigning or incrementing any PetscInt row index.
inline std::size_t CheckedImmersedStaticFlowRowCount(std::size_t active_nodes,
	std::size_t controller_count, bool include_gauge)
{
	if (active_nodes > std::numeric_limits<std::size_t>::max()/4)
		throw std::overflow_error("immersed physical dof count overflows");
	const std::size_t physical_dofs = 4*active_nodes;
	if (controller_count > std::numeric_limits<std::size_t>::max()-physical_dofs)
		throw std::overflow_error("immersed controller dof count overflows");
	std::size_t total_dofs = physical_dofs+controller_count;
	if (include_gauge) {
		if (total_dofs == std::numeric_limits<std::size_t>::max())
			throw std::overflow_error("immersed gauge dof count overflows");
		++total_dofs;
	}
	if (total_dofs > static_cast<std::size_t>(std::numeric_limits<PetscInt>::max()))
		throw std::overflow_error("PETSc row count overflows");
	return total_dofs;
}

class ImmersedStaticFlowRuntime {
public:
	ImmersedStaticFlowRuntime(const CartesianDomainClassification& domain,
		const CutCellVolumeQuadratureCatalog& volume,
		const ImmersedSurfaceQuadratureCatalog& surface,
		const CutCellGhostPenaltyCatalog& ghost, ImmersedStaticFlowOptions options = {})
		: domain_(domain), volume_(volume), surface_(surface), ghost_(ghost), options_(std::move(options))
	{
		if (!options_.wall_velocity || !options_.body_force || !std::isfinite(options_.parameters.density) || !(options_.parameters.density > 0.0)
			|| !std::isfinite(options_.parameters.dynamic_viscosity) || !(options_.parameters.dynamic_viscosity > 0.0)
			|| !std::isfinite(options_.parameters.dt) || options_.parameters.dt != 0.0 || !std::isfinite(options_.wall_gamma0)
			|| !(options_.wall_gamma0 > 0.0) || options_.nonlinear_maximum_iterations <= 0
			|| options_.ksp_maximum_iterations <= 0 || !std::isfinite(options_.ksp_relative_tolerance) || !(options_.ksp_relative_tolerance > 0.0)
			|| !std::isfinite(options_.nonlinear_relative_tolerance) || !(options_.nonlinear_relative_tolerance > 0.0)
			|| !std::isfinite(options_.nonlinear_absolute_tolerance) || !(options_.nonlinear_absolute_tolerance > 0.0)
			|| !std::isfinite(options_.flow_controller_relative_tolerance) || !(options_.flow_controller_relative_tolerance > 0.0)
			|| !std::isfinite(options_.flow_controller_absolute_tolerance_m3_s) || !(options_.flow_controller_absolute_tolerance_m3_s >= 0.0)
			|| !std::isfinite(options_.flow_controller_reference_flow_m3_s) || !(options_.flow_controller_reference_flow_m3_s > 0.0)
			|| !std::isfinite(options_.nonlinear_block_reduction) || options_.nonlinear_block_reduction < 0.0 || !std::isfinite(options_.minimum_damping) || !(options_.minimum_damping > 0.0)
			|| options_.minimum_damping > 1.0 || !std::isfinite(options_.lu_pivot_shift) || options_.lu_pivot_shift < 0.0)
			throw std::invalid_argument("immersed static-flow options are invalid");
		ghost_.ValidateBinding(domain_, volume_);
		if (!surface_.Usable() || !ghost_.Usable()) throw std::invalid_argument("immersed static-flow catalogs are unusable");
		ValidateImmersedNitscheWallLabels(options_.wall_labels);
		if (surface_.GridSpec().lower_m != volume_.GridSpec().lower_m || surface_.GridSpec().upper_m != volume_.GridSpec().upper_m
			|| surface_.GridSpec().cells != volume_.GridSpec().cells || surface_.SurfaceCanonicalHash() != volume_.SurfaceCanonicalHash()
			|| surface_.Cells().size() != volume_.Cells().size()) throw std::invalid_argument("immersed surface and volume catalogs do not have an exact common binding");
		ConfigurePorts();
		if (options_.wall_labels.empty() && options_.ports.empty()) throw std::invalid_argument("immersed static-flow requires a wall or open port label");
		PreflightCatalogs();
		BuildActiveMap();
		// Catalogs and active compaction are immutable, so stream this once.
		// This preserves compact-rule laziness while avoiding an extra full
		// volume traversal for every Newton assembly.
		BuildGaugeWeights();
		const PetscInt n = CheckedPetscCount(diagnostics_.total_dofs);
		try {
			Check(MatCreateSeqAIJ(PETSC_COMM_SELF, n, n, 300, nullptr, &jacobian_), "MatCreateSeqAIJ");
			Check(MatSetOption(jacobian_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE), "MatSetOption");
			// All appended scalar rows are saddle-point rows.  Seed their exact-zero
			// diagonals before zeros are ignored, retaining the intended AIJ pattern
			// without altering the mathematical operator.
			Check(MatSetOption(jacobian_, MAT_IGNORE_ZERO_ENTRIES, PETSC_FALSE), "MatSetOption retain scalar diagonals");
			for (const PetscInt row : AppendedScalarRows())
				Check(MatSetValue(jacobian_, row, row, 0.0, INSERT_VALUES), "MatSetValue scalar structural diagonal");
			Check(MatAssemblyBegin(jacobian_, MAT_FINAL_ASSEMBLY), "MatAssemblyBegin scalar structural diagonals"); Check(MatAssemblyEnd(jacobian_, MAT_FINAL_ASSEMBLY), "MatAssemblyEnd scalar structural diagonals");
			AuditAppendedScalarDiagonals();
			Check(MatSetOption(jacobian_, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE), "MatSetOption ignore zero entries");
			Check(VecCreateSeq(PETSC_COMM_SELF, n, &state_), "VecCreateSeq");
			Check(VecDuplicate(state_, &committed_), "VecDuplicate committed"); Check(VecDuplicate(state_, &prepared_), "VecDuplicate prepared"); Check(VecDuplicate(state_, &rhs_), "VecDuplicate rhs"); Check(VecDuplicate(state_, &update_), "VecDuplicate update");
			Check(VecDuplicate(state_, &action_input_), "VecDuplicate Jacobian input"); Check(VecDuplicate(state_, &action_output_), "VecDuplicate Jacobian output");
			Check(VecDuplicate(state_, &constant_pressure_), "VecDuplicate pressure constant"); Check(VecDuplicate(state_, &pressure_defect_), "VecDuplicate pressure defect");
			Check(VecSet(state_, 0.0), "VecSet state"); Check(VecSet(committed_, 0.0), "VecSet committed");
			Check(VecSet(constant_pressure_, 0.0), "VecSet pressure constant");
			for (std::size_t a = 0; a < active_nodes_.size(); ++a) Check(VecSetValue(constant_pressure_, 4*CheckedPetscCount(a)+3, 1.0, INSERT_VALUES), "VecSetValue pressure constant");
			Check(VecAssemblyBegin(constant_pressure_), "VecAssemblyBegin pressure constant"); Check(VecAssemblyEnd(constant_pressure_), "VecAssemblyEnd pressure constant");
			Check(KSPCreate(PETSC_COMM_SELF, &ksp_), "KSPCreate"); Check(KSPSetType(ksp_, KSPGMRES), "KSPSetType");
			Check(KSPSetTolerances(ksp_, options_.ksp_relative_tolerance, PETSC_DEFAULT, PETSC_DEFAULT, options_.ksp_maximum_iterations), "KSPSetTolerances");
			PC pc = nullptr; Check(KSPGetPC(ksp_, &pc), "KSPGetPC"); Check(PCSetType(pc, PCLU), "PCSetType");
			// A requested pivot shift is confined to the LU preconditioner; GMRES
			// continues to apply the unmodified assembled operator.
			if (options_.lu_pivot_shift > 0.0) {
				Check(PCFactorSetShiftType(pc, MAT_SHIFT_NONZERO), "PCFactorSetShiftType");
				Check(PCFactorSetShiftAmount(pc, options_.lu_pivot_shift), "PCFactorSetShiftAmount");
			}
		} catch (...) { Destroy(); throw; }
	}

	~ImmersedStaticFlowRuntime()
	{
		Destroy();
	}
	ImmersedStaticFlowRuntime(const ImmersedStaticFlowRuntime&) = delete;
	ImmersedStaticFlowRuntime& operator=(const ImmersedStaticFlowRuntime&) = delete;

	const std::vector<std::int32_t>& ActiveNodes() const noexcept { return active_nodes_; }
	const ImmersedStaticFlowDiagnostics& Diagnostics() const noexcept { return diagnostics_; }
	const std::vector<ImmersedFlowPortDefinition>& PortDefinitions() const noexcept { return options_.ports; }
	// Focused transactional tests can prove that an adapter leaves its published
	// image untouched when backend preparation fails before publication.
	void FailNextPrepareForTesting() noexcept { fail_next_prepare_for_testing_ = true; }
	bool ScalarDiagonalStructureVerified() const noexcept { return diagnostics_.scalar_diagonal_structure_verified; }
	bool HasGauge() const noexcept { return diagnostics_.gauge_present; }
	PetscInt GaugeDof() const noexcept { return diagnostics_.gauge_row; }
	PetscInt PortMultiplierDof(const std::string& id) const
	{
		for (const auto& port : diagnostics_.ports) if (port.id == id) return port.multiplier_row;
		throw std::out_of_range("immersed flow port id is not configured");
	}
	void SetPortControlValue(const std::string& id, double value)
	{
		if (diagnostics_.trial_active) throw std::logic_error("cannot change immersed port control during an active trial");
		if (!std::isfinite(value)) throw std::invalid_argument("immersed flow port value must be finite");
		for (std::size_t i = 0; i < options_.ports.size(); ++i) if (options_.ports[i].id == id) {
			options_.ports[i].value = value; diagnostics_.ports[i].target = value;
			return;
		}
		throw std::out_of_range("immersed flow port id is not configured");
	}
	PetscInt Dof(std::int32_t node, int field) const
	{
		if (field < 0 || field > 3 || node < 0 || static_cast<std::size_t>(node) >= node_to_active_.size()
			|| node_to_active_[static_cast<std::size_t>(node)] < 0) throw std::out_of_range("background node is not active");
		return 4*node_to_active_[static_cast<std::size_t>(node)]+field;
	}
	std::vector<PetscScalar> CommittedState() const { return CopyVector(committed_); }
	std::vector<PetscScalar> TrialState() const { return CopyVector(state_); }
	std::vector<PetscScalar> AssembledNegativeResidual() const { return CopyVector(rhs_); }
	std::vector<PetscScalar> AssembledJacobianAction(const std::vector<PetscScalar>& values) const
	{
		if (values.size() != diagnostics_.total_dofs) throw std::invalid_argument("Jacobian action vector size is invalid");
		PetscScalar* data = nullptr; Check(VecGetArray(action_input_, &data), "VecGetArray Jacobian input");
		for (std::size_t i = 0; i < values.size(); ++i) data[i] = values[i];
		Check(VecRestoreArray(action_input_, &data), "VecRestoreArray Jacobian input");
		Check(MatMult(jacobian_, action_input_, action_output_), "MatMult Jacobian action");
		return CopyVector(action_output_);
	}
	std::vector<PetscScalar> AssembledJacobianDense() const
	{
		const PetscInt n = CheckedPetscCount(diagnostics_.total_dofs);
		std::vector<PetscInt> rows(static_cast<std::size_t>(n)); std::iota(rows.begin(), rows.end(), PetscInt{0});
		std::vector<PetscScalar> result(static_cast<std::size_t>(n)*static_cast<std::size_t>(n));
		Check(MatGetValues(jacobian_, n, rows.data(), n, rows.data(), result.data()), "MatGetValues dense Jacobian");
		return result;
	}

	void SetCommittedState(const std::vector<PetscScalar>& values)
	{
		if (diagnostics_.trial_active) throw std::logic_error("cannot replace committed state during an active trial");
		if (values.size() != diagnostics_.total_dofs) throw std::invalid_argument("committed state size is invalid");
		for (const auto value : values) if (!std::isfinite(PetscRealPart(value))) throw std::invalid_argument("committed state is not finite");
		PetscScalar* target = nullptr; Check(VecGetArray(committed_, &target), "VecGetArray committed");
		for (std::size_t i = 0; i < values.size(); ++i) target[i] = values[i];
		Check(VecRestoreArray(committed_, &target), "VecRestoreArray committed"); Check(VecCopy(committed_, state_), "VecCopy committed state");
	}

	// No global dense temporary is used: the only dense objects are existing
	// 64-node element and <=80-node face blocks supplied by the catalog APIs.
	void Assemble()
	{
		ValidateAllFlowCompatibility();
		Check(MatZeroEntries(jacobian_), "MatZeroEntries"); Check(VecSet(rhs_, 0.0), "VecSet rhs");
		diagnostics_.volume_cells = diagnostics_.surface_cells = diagnostics_.ghost_faces = 0;
		diagnostics_.wall_selected_points = 0;
		diagnostics_.wall_selected_points_by_label.clear();
		for (auto& port : diagnostics_.ports) port.assembled_surface_points = 0;
		for (std::uint64_t id = 0; id < domain_.Cells().size(); ++id) {
			if (!UsablePositive(id)) continue;
			const auto element = domain_.Background().MaterializeElement(id);
			const auto nodal = Gather(element);
			NavierStokesSystem volume_system = BuildVolumeSystem(element, nodal, id);
			if (options_.assemble_volume) { ScatterElement(element.connectivity, volume_system); ++diagnostics_.volume_cells; }
			if (domain_.Cells()[static_cast<std::size_t>(id)].classification == CellClassification::Cut) {
				if (!ghost_.Covered(id)) throw std::runtime_error("covered immersed wall and ghost policy mismatch");
				const auto& rule = surface_.UsableRule(domain_, id);
				// With no ports, retain the Phase-5 selected-wall preflight exactly:
				// every positive cut cell must produce a selected wall contribution.
				if (options_.ports.empty() || HasSelectedWallPoint(rule)) {
					const auto wall = BuildImmersedNitscheWallElementFromVolumeSystem(domain_, volume_, surface_, id, nodal, {},
						options_.parameters, options_.wall_labels, volume_system, ghost_, options_.wall_gamma0, options_.wall_velocity);
					const std::size_t wall_points = WallContributionCount(wall.diagnostics);
					if (wall_points == 0) throw std::runtime_error("selected immersed wall labels produced no surface contribution");
					// The wall builder returns volume plus wall terms.  Always form the
					// wall delta, even for a wall-only diagnostic assembly: otherwise
					// disabling the volume switch silently double-counts cut-cell volume.
					NavierStokesSystem wall_delta = wall.system;
					for (std::size_t i = 0; i < wall_delta.jacobian.size(); ++i) wall_delta.jacobian[i] -= volume_system.jacobian[i];
					for (std::size_t i = 0; i < wall_delta.negative_residual.size(); ++i) wall_delta.negative_residual[i] -= volume_system.negative_residual[i];
					if (options_.assemble_wall) {
						ScatterElement(element.connectivity, wall_delta); ++diagnostics_.surface_cells; diagnostics_.wall_selected_points += wall_points;
						for (const auto& entry : wall.diagnostics.by_boundary_id) diagnostics_.wall_selected_points_by_label[entry.first] += entry.second.selected_points;
					}
				}
				for (std::size_t port = 0; port < options_.ports.size(); ++port)
					if (RuleHasLabel(rule, options_.ports[port].boundary_label))
						ScatterPortElement(element, nodal, rule, port);
			}
		}
		for (std::size_t face = 0; options_.assemble_ghost && face < ghost_.Faces().size(); ++face) {
			VecReadArray values(state_, "VecGetArrayRead ghost face");
			const auto block = ghost_.AssembleFaceLocal(face, domain_, volume_, [this, data = values.Data()](std::int32_t node, int field) { return PetscRealPart(data[Dof(node, field)]); }, options_.parameters.dynamic_viscosity);
			values.Close("VecRestoreArrayRead ghost face");
			ScatterBlock(block.connectivity, block.jacobian, block.negative_residual); ++diagnostics_.ghost_faces;
		}
		Check(MatAssemblyBegin(jacobian_, MAT_FINAL_ASSEMBLY), "MatAssemblyBegin physical"); Check(MatAssemblyEnd(jacobian_, MAT_FINAL_ASSEMBLY), "MatAssemblyEnd physical");
		Check(VecAssemblyBegin(rhs_), "VecAssemblyBegin physical"); Check(VecAssemblyEnd(rhs_), "VecAssemblyEnd physical");
		MeasurePorts();
		MeasureConstantPressureDefect();
		if (options_.assemble_gauge && HasGauge()) InsertGauge();
		Check(MatAssemblyBegin(jacobian_, MAT_FINAL_ASSEMBLY), "MatAssemblyBegin gauge"); Check(MatAssemblyEnd(jacobian_, MAT_FINAL_ASSEMBLY), "MatAssemblyEnd gauge");
		Check(VecAssemblyBegin(rhs_), "VecAssemblyBegin gauge"); Check(VecAssemblyEnd(rhs_), "VecAssemblyEnd gauge");
	}

	// Start each Newton attempt from an exact copy of the committed vector.
	bool SolveTrial()
	{
		Check(VecCopy(committed_, state_), "VecCopy begin trial"); diagnostics_.trial_active = true; diagnostics_.committed = false;
		diagnostics_.converged = false; diagnostics_.nonlinear_iterations = diagnostics_.ksp_iterations = 0;
		diagnostics_.ksp_reason = KSP_CONVERGED_ITERATING; diagnostics_.damping = 0.0; diagnostics_.residual_norm = 0.0;
		diagnostics_.newton_steps.clear();
		try {
			double initial = -1.0; std::array<double, 3> initial_blocks{{0.0,0.0,0.0}};
			for (PetscInt iteration = 0; iteration < options_.nonlinear_maximum_iterations; ++iteration) {
				Assemble(); PetscReal residual = 0.0; Check(VecNorm(rhs_, NORM_2, &residual), "VecNorm residual");
				if (!std::isfinite(residual)) throw std::runtime_error("nonlinear residual is not finite");
				if (initial < 0.0) { initial = residual; initial_blocks = ResidualBlockNorms(CopyVector(rhs_)); }
				diagnostics_.residual_norm = residual; diagnostics_.nonlinear_iterations = iteration;
				if (residual <= std::max(options_.nonlinear_absolute_tolerance, options_.nonlinear_relative_tolerance*initial)
					&& BlockReductionSatisfied(initial_blocks, CopyVector(rhs_), initial) && ControllerConvergenceSatisfied()) {
					diagnostics_.converged = true; return true;
				}
				ImmersedStaticFlowNewtonStep step; step.iteration = iteration; step.residual_norm = residual;
				// Preserve the assembled, unpreconditioned right-hand side for the
				// independently reported true linear residual.  KSP implementations
				// are permitted to use their input work vector internally.
				Check(VecCopy(rhs_, action_input_), "VecCopy linear right hand side");
				Check(KSPSetOperators(ksp_, jacobian_, jacobian_), "KSPSetOperators"); Check(KSPSolve(ksp_, rhs_, update_), "KSPSolve");
				Check(KSPGetConvergedReason(ksp_, &diagnostics_.ksp_reason), "KSPGetConvergedReason"); PetscInt work = 0; Check(KSPGetIterationNumber(ksp_, &work), "KSPGetIterationNumber"); diagnostics_.ksp_iterations += work;
				PetscReal ksp_residual = 0.0; Check(KSPGetResidualNorm(ksp_, &ksp_residual), "KSPGetResidualNorm");
				step.ksp_iterations = work; step.ksp_reason = diagnostics_.ksp_reason; step.ksp_residual_norm = ksp_residual;
				if (diagnostics_.ksp_reason <= 0) throw std::runtime_error("static immersed-flow KSP failed with reason "+std::to_string(static_cast<int>(diagnostics_.ksp_reason)));
				EnsureFinite(update_);
				PetscReal update_norm = 0.0, linear_residual = 0.0; Check(VecNorm(update_, NORM_2, &update_norm), "VecNorm update");
				Check(MatMult(jacobian_, update_, action_output_), "MatMult linear residual"); Check(VecAXPY(action_output_, -1.0, action_input_), "VecAXPY linear residual"); Check(VecNorm(action_output_, NORM_2, &linear_residual), "VecNorm linear residual");
				step.update_norm = update_norm; step.linear_residual_norm = linear_residual; step.linear_relative_residual = linear_residual/residual;
				double damping = 1.0, old = residual; bool accepted = false;
				while (damping >= options_.minimum_damping) {
					Check(VecAXPY(state_, damping, update_), "VecAXPY trial update"); Assemble(); PetscReal candidate = 0.0; Check(VecNorm(rhs_, NORM_2, &candidate), "VecNorm candidate");
				if (std::isfinite(candidate) && candidate < old) {
					accepted = true; step.damping = damping; step.candidate_residual_norm = candidate; diagnostics_.damping = damping; diagnostics_.residual_norm = candidate;
					++diagnostics_.nonlinear_iterations;
					diagnostics_.newton_steps.push_back(step);
					if (candidate <= std::max(options_.nonlinear_absolute_tolerance, options_.nonlinear_relative_tolerance*initial)
						&& BlockReductionSatisfied(initial_blocks, CopyVector(rhs_), initial) && ControllerConvergenceSatisfied()) {
						diagnostics_.converged = true; return true;
					}
					break;
				}
					Check(VecAXPY(state_, -damping, update_), "VecAXPY trial undo"); damping *= .5;
				}
				if (!accepted) { diagnostics_.newton_steps.push_back(step); throw std::runtime_error("static immersed-flow backtracking failed"); }
			}
			throw std::runtime_error("static immersed-flow nonlinear iteration cap reached");
		} catch (...) { Rollback(); throw; }
	}
	void Commit()
	{
		PrepareCommit(); FinalizeCommit();
	}
	// All PETSc work is deliberately completed here.  FinalizeCommit only swaps
	// already-owned handles and updates scalar bookkeeping, so callers can use
	// it as the nonthrowing publication half of a coupled transaction.
	void PrepareCommit()
	{
		if (!diagnostics_.trial_active || !diagnostics_.converged || diagnostics_.prepared)
			throw std::logic_error("cannot prepare an unconverged or already prepared static-flow trial");
		if (fail_next_prepare_for_testing_) {
			fail_next_prepare_for_testing_ = false;
			throw std::runtime_error("injected immersed static-flow prepare failure");
		}
		Check(VecCopy(state_, prepared_), "VecCopy prepare commit");
		diagnostics_.prepared = true; ++diagnostics_.prepare_count;
	}
	void FinalizeCommit() noexcept
	{
		if (!diagnostics_.prepared) return;
		std::swap(committed_, prepared_);
		diagnostics_.prepared = false; diagnostics_.trial_active = false;
		diagnostics_.committed = true; ++diagnostics_.commit_count;
		++diagnostics_.finalize_count;
	}
	void AbortPrepared() noexcept
	{
		diagnostics_.prepared = false;
	}
	void Rollback()
	{
		if (!diagnostics_.trial_active) return;
		AbortPrepared();
		if (state_ && committed_) Check(VecCopy(committed_, state_), "VecCopy rollback");
		diagnostics_.trial_active = false; diagnostics_.committed = true; diagnostics_.converged = false; ++diagnostics_.rollback_count;
	}

private:
	static void Check(PetscErrorCode code, const char* operation)
	{
		if (code != 0) throw std::runtime_error(std::string("PETSc ")+operation+" failed with error "+std::to_string(static_cast<long long>(code)));
	}
	class VecReadArray {
	public:
		VecReadArray(Vec vector, const char* operation) : vector_(vector)
		{
			Check(VecGetArrayRead(vector_, &values_), operation); open_ = true;
		}
		~VecReadArray() noexcept
		{
			if (open_) { const PetscScalar* values = values_; VecRestoreArrayRead(vector_, &values); }
		}
		const PetscScalar* Data() const noexcept { return values_; }
		void Close(const char* operation)
		{
			if (!open_) return;
			const PetscScalar* values = values_;
			const PetscErrorCode code = VecRestoreArrayRead(vector_, &values);
			if (code == 0) { values_ = nullptr; open_ = false; return; }
			Check(code, operation);
		}
	private:
		Vec vector_ = nullptr;
		const PetscScalar* values_ = nullptr;
		bool open_ = false;
	};
	void Destroy() noexcept
	{
		if (ksp_) KSPDestroy(&ksp_);
		if (update_) VecDestroy(&update_);
		if (pressure_defect_) VecDestroy(&pressure_defect_);
		if (constant_pressure_) VecDestroy(&constant_pressure_);
		if (action_output_) VecDestroy(&action_output_);
		if (action_input_) VecDestroy(&action_input_);
		if (rhs_) VecDestroy(&rhs_);
		if (prepared_) VecDestroy(&prepared_);
		if (committed_) VecDestroy(&committed_);
		if (state_) VecDestroy(&state_);
		if (jacobian_) MatDestroy(&jacobian_);
	}
	static PetscInt CheckedPetscCount(std::size_t value)
	{
		if (value > static_cast<std::size_t>(std::numeric_limits<PetscInt>::max())) throw std::overflow_error("PETSc row count overflows");
		return static_cast<PetscInt>(value);
	}
	std::vector<PetscInt> AppendedScalarRows() const
	{
		std::vector<PetscInt> rows;
		for (const auto& port : diagnostics_.ports) if (port.multiplier_row >= 0) rows.push_back(port.multiplier_row);
		if (HasGauge()) rows.push_back(GaugeDof());
		return rows;
	}
	void AuditAppendedScalarDiagonals()
	{
		for (const PetscInt row : AppendedScalarRows()) {
			PetscInt count = 0;
			const PetscInt* columns = nullptr;
			const PetscScalar* values = nullptr;
			Check(MatGetRow(jacobian_, row, &count, &columns, &values), "MatGetRow scalar structural diagonal");
			bool found = false, exact_zero = false;
			for (PetscInt entry = 0; entry < count; ++entry) if (columns[entry] == row) {
				found = true; exact_zero = PetscRealPart(values[entry]) == 0.0; break;
			}
			Check(MatRestoreRow(jacobian_, row, &count, &columns, &values), "MatRestoreRow scalar structural diagonal");
			if (!found || !exact_zero) throw std::runtime_error("immersed scalar row lacks an exact-zero structural diagonal");
		}
		diagnostics_.scalar_diagonal_structure_verified = true;
	}
	bool UsablePositive(std::uint64_t id) const
	{
		const auto& cell = volume_.Cell(id);
		const auto classification = domain_.Cells()[static_cast<std::size_t>(id)].classification;
		return cell.usable && !CertifiedEmpty(id)
			&& (classification == CellClassification::Inside || classification == CellClassification::Cut);
	}
	bool CertifiedEmpty(std::uint64_t id) const
	{
		const auto& cell = volume_.Cell(id);
		if (volume_.StorageMode() == CutCellVolumeQuadratureStorageMode::Expanded)
			return cell.rule.Points().empty();
		return CompactCutCellVolumeLogicalPointCount(cell.compact_rule) == 0;
	}
	static bool RuleHasLabel(const SurfaceQuadratureRule& rule, int label)
	{
		return std::any_of(rule.Points().begin(), rule.Points().end(), [label](const SurfaceQuadraturePoint& point) { return point.boundary_id == label; });
	}
	bool HasSelectedWallPoint(const SurfaceQuadratureRule& rule) const
	{
		return std::any_of(rule.Points().begin(), rule.Points().end(), [this](const SurfaceQuadraturePoint& point) {
			return std::binary_search(options_.wall_labels.begin(), options_.wall_labels.end(), point.boundary_id);
		});
	}
	void ValidateAllFlowCompatibility() const
	{
		if (options_.ports.empty() || std::any_of(options_.ports.begin(), options_.ports.end(), [](const ImmersedFlowPortDefinition& port) { return IsImmersedFlowPressureLike(port.control_mode); })) return;
		long double sum = 0.0L, scale = 0.0L;
		for (const auto& port : options_.ports) { sum += port.value; scale += std::abs(static_cast<long double>(port.value)); }
		const long double tolerance = static_cast<long double>(options_.flow_controller_absolute_tolerance_m3_s)
			+ static_cast<long double>(options_.flow_controller_relative_tolerance)*std::max(scale, static_cast<long double>(options_.flow_controller_reference_flow_m3_s));
		if (std::abs(sum) > tolerance) throw std::invalid_argument("all flow-controlled immersed ports require compatible net outward flow");
	}
	void ConfigurePorts()
	{
		// Preserve the Phase-5 closed-surface contract exactly: there is no
		// new label-completeness or source-area policy when no ports are present.
		if (options_.ports.empty()) { diagnostics_.ports.clear(); diagnostics_.gauge_present = true; return; }
		std::vector<std::string> ids; std::vector<int> labels;
		for (const auto& port : options_.ports) { ValidateImmersedFlowPortDefinition(port); ids.push_back(port.id); labels.push_back(port.boundary_label); }
		std::sort(ids.begin(), ids.end()); if (std::adjacent_find(ids.begin(), ids.end()) != ids.end()) throw std::invalid_argument("immersed flow port ids must be unique");
		std::sort(labels.begin(), labels.end()); if (std::adjacent_find(labels.begin(), labels.end()) != labels.end()) throw std::invalid_argument("immersed flow port labels must be unique");
		for (const int label : labels) if (std::binary_search(options_.wall_labels.begin(), options_.wall_labels.end(), label)) throw std::invalid_argument("immersed flow port labels must be disjoint from wall labels");
		const auto& catalog = surface_.Diagnostics();
		for (const int label : options_.wall_labels) {
			const auto area = catalog.source_area_by_boundary_id.find(static_cast<std::uint32_t>(label));
			if (area == catalog.source_area_by_boundary_id.end() || !(area->second > 0.0) || !std::isfinite(area->second)) throw std::invalid_argument("immersed wall label has no positive source area");
		}
		for (const auto& entry : catalog.source_area_by_boundary_id) {
			const int label = static_cast<int>(entry.first);
			if (!std::binary_search(options_.wall_labels.begin(), options_.wall_labels.end(), label)
				&& !std::binary_search(labels.begin(), labels.end(), label))
				throw std::invalid_argument("immersed surface label is neither a Nitsche wall nor an open port");
		}
		diagnostics_.ports.clear();
		for (const auto& port : options_.ports) {
			const auto source = catalog.source_area_by_boundary_id.find(static_cast<std::uint32_t>(port.boundary_label));
			if (source == catalog.source_area_by_boundary_id.end() || !(source->second > 0.0) || !std::isfinite(source->second)) throw std::invalid_argument("immersed flow port label has no positive source area");
			long double retained = 0.0L;
			for (std::uint64_t id = 0; id < surface_.Cells().size(); ++id)
				for (const auto& point : surface_.UsableRule(domain_, id).Points()) if (point.boundary_id == port.boundary_label) retained += point.weight;
			const long double scale = std::max(std::abs(retained), std::abs(static_cast<long double>(source->second)));
			if (!std::isfinite(static_cast<double>(retained)) || std::abs(retained-static_cast<long double>(source->second)) > 1024.0L*std::numeric_limits<double>::epsilon()*std::max(scale, std::numeric_limits<long double>::denorm_min())) throw std::runtime_error("immersed flow port retained area fails the surface catalog audit tolerance");
			ImmersedStaticFlowDiagnostics::Port diagnostic; diagnostic.id = port.id; diagnostic.boundary_label = port.boundary_label; diagnostic.control_mode = port.control_mode; diagnostic.area_m2 = static_cast<double>(retained); diagnostic.target = port.value;
			diagnostics_.ports.push_back(std::move(diagnostic));
		}
		ValidateAllFlowCompatibility();
		diagnostics_.gauge_present = std::none_of(options_.ports.begin(), options_.ports.end(), [](const ImmersedFlowPortDefinition& port) { return IsImmersedFlowPressureLike(port.control_mode); });
	}
	void PreflightCatalogs() const
	{
		for (std::uint64_t id = 0; id < domain_.Cells().size(); ++id) {
			const auto classification = domain_.Cells()[static_cast<std::size_t>(id)].classification;
			if (classification != CellClassification::Inside && classification != CellClassification::Cut) continue;
			if (!volume_.Cell(id).usable)
				throw std::runtime_error("classified immersed-flow cell has an unusable volume rule");
			if (volume_.StorageMode() == CutCellVolumeQuadratureStorageMode::Expanded) volume_.ValidateUsableRule(domain_, id);
			else volume_.ValidateUsableCompactRule(domain_, id);
			if (CertifiedEmpty(id)) {
				if (classification == CellClassification::Inside)
					throw std::runtime_error("inside immersed-flow cell has a certified-empty volume rule");
				continue;
			}
			if (!UsablePositive(id)) continue;
			if (classification != CellClassification::Cut) continue;
			if (!ghost_.Covered(id)) throw std::runtime_error("covered immersed Nitsche policy requires every positive cut cell to be ghost-covered");
			surface_.ValidateUsableRule(domain_, id);
			const auto& rule = surface_.UsableRule(domain_, id);
			// Preserve the Phase-5 no-port preflight at construction time.  Open
			// ports deliberately relax this only for cap-only cells; ghost coverage
			// and ConfigurePorts() label completeness remain mandatory in both modes.
			if (options_.ports.empty() && !HasSelectedWallPoint(rule))
				throw std::runtime_error("positive cut cell has no selected immersed wall surface contribution");
		}
	}
	void BuildActiveMap()
	{
		node_to_active_.assign(static_cast<std::size_t>(domain_.Background().NodeCount()), -1);
		for (std::uint64_t id = 0; id < domain_.Cells().size(); ++id) if (UsablePositive(id)) {
			const auto e = domain_.Background().MaterializeElement(id);
			active_nodes_.insert(active_nodes_.end(), e.connectivity.begin(), e.connectivity.end());
		}
		std::sort(active_nodes_.begin(), active_nodes_.end()); active_nodes_.erase(std::unique(active_nodes_.begin(), active_nodes_.end()), active_nodes_.end());
		if (active_nodes_.empty()) throw std::runtime_error("active-node compaction is empty");
		const std::size_t controllers = static_cast<std::size_t>(std::count_if(options_.ports.begin(), options_.ports.end(), [](const ImmersedFlowPortDefinition& port) { return port.control_mode == ImmersedFlowPortControlMode::FlowRate; }));
		const std::size_t total_dofs = CheckedImmersedStaticFlowRowCount(active_nodes_.size(), controllers, HasGauge());
		diagnostics_.active_nodes = active_nodes_.size();
		diagnostics_.physical_dofs = 4*active_nodes_.size();
		diagnostics_.total_dofs = total_dofs;
		for (std::size_t i = 0; i < active_nodes_.size(); ++i) {
			const PetscInt index = CheckedPetscCount(i);
			node_to_active_[static_cast<std::size_t>(active_nodes_[i])] = index;
		}
		std::size_t next_scalar_row = diagnostics_.physical_dofs;
		for (std::size_t i = 0; i < diagnostics_.ports.size(); ++i) if (options_.ports[i].control_mode == ImmersedFlowPortControlMode::FlowRate) {
			if (next_scalar_row >= total_dofs) throw std::logic_error("immersed controller row topology is inconsistent");
			diagnostics_.ports[i].multiplier_row = CheckedPetscCount(next_scalar_row);
			++next_scalar_row;
		}
		if (HasGauge()) {
			if (next_scalar_row >= total_dofs) throw std::logic_error("immersed gauge row topology is inconsistent");
			diagnostics_.gauge_row = CheckedPetscCount(next_scalar_row);
			++next_scalar_row;
		} else diagnostics_.gauge_row = -1;
		if (next_scalar_row != total_dofs) throw std::logic_error("immersed scalar row count is inconsistent");
	}
	std::vector<std::array<double, 4>> Gather(const Element& element) const
	{
		VecReadArray values(state_, "VecGetArrayRead gather"); std::vector<std::array<double, 4>> result(element.connectivity.size());
		for (std::size_t a = 0; a < result.size(); ++a) for (int f = 0; f < 4; ++f) result[a][f] = PetscRealPart(values.Data()[Dof(element.connectivity[a], f)]);
		values.Close("VecRestoreArrayRead gather"); return result;
	}
	NavierStokesSystem BuildVolumeSystem(const Element& element, const std::vector<std::array<double, 4>>& nodal, std::uint64_t id) const
	{
		if (volume_.StorageMode() == CutCellVolumeQuadratureStorageMode::Expanded)
			return BuildNavierStokesElement(element, nodal, {}, options_.parameters, volume_.UsableRule(domain_, id), options_.body_force);
		const auto& rule = volume_.UsableCompactRule(domain_, id);
		return BuildNavierStokesElementFromPoints(element, nodal, {}, options_.parameters,
			[&rule](const auto& consume) { ForEachVolumePoint(rule, consume); }, options_.body_force);
	}
	void ScatterElement(const std::vector<std::int32_t>& nodes, const NavierStokesSystem& system)
	{
		std::vector<PetscInt> rows; rows.reserve(4*nodes.size()); for (const auto node : nodes) for (int f = 0; f < 4; ++f) rows.push_back(Dof(node, f));
		for (const auto value : system.jacobian) if (!std::isfinite(PetscRealPart(value))) throw std::overflow_error("immersed static-flow element matrix value is not finite");
		for (const auto value : system.negative_residual) if (!std::isfinite(PetscRealPart(value))) throw std::overflow_error("immersed static-flow element residual value is not finite");
		Check(MatSetValues(jacobian_, CheckedPetscCount(rows.size()), rows.data(), CheckedPetscCount(rows.size()), rows.data(), system.jacobian.data(), ADD_VALUES), "MatSetValues element");
		Check(VecSetValues(rhs_, CheckedPetscCount(rows.size()), rows.data(), system.negative_residual.data(), ADD_VALUES), "VecSetValues element");
	}
	void ScatterBlock(const std::vector<std::int32_t>& nodes, const std::vector<PetscScalar>& matrix, const std::vector<PetscScalar>& residual)
	{
		std::vector<PetscInt> rows; rows.reserve(4*nodes.size()); for (const auto node : nodes) for (int f = 0; f < 4; ++f) rows.push_back(Dof(node, f));
		for (const auto value : matrix) if (!std::isfinite(PetscRealPart(value))) throw std::overflow_error("immersed static-flow block matrix value is not finite");
		for (const auto value : residual) if (!std::isfinite(PetscRealPart(value))) throw std::overflow_error("immersed static-flow block residual value is not finite");
		Check(MatSetValues(jacobian_, CheckedPetscCount(rows.size()), rows.data(), CheckedPetscCount(rows.size()), rows.data(), matrix.data(), ADD_VALUES), "MatSetValues block"); Check(VecSetValues(rhs_, CheckedPetscCount(rows.size()), rows.data(), residual.data(), ADD_VALUES), "VecSetValues block");
	}
	void ScatterPortElement(const Element& element, const std::vector<std::array<double, 4>>& nodal,
		const SurfaceQuadratureRule& rule, std::size_t port_index)
	{
		const auto& definition = options_.ports.at(port_index);
		const auto local = BuildImmersedFlowPortElement(element, rule, definition.boundary_label,
			definition.control_mode, definition.value, nodal);
		for (const auto& point : rule.Points()) if (point.boundary_id == definition.boundary_label)
			++diagnostics_.ports.at(port_index).assembled_surface_points;
		std::vector<PetscInt> rows; rows.reserve(4*element.connectivity.size());
		for (const auto node : element.connectivity) for (int field = 0; field < 4; ++field) rows.push_back(Dof(node, field));
		if (IsImmersedFlowPressureLike(definition.control_mode)) {
			Check(VecSetValues(rhs_, CheckedPetscCount(rows.size()), rows.data(), local.negative_residual.data(), ADD_VALUES), "VecSetValues port pressure traction");
			return;
		}
		if (definition.control_mode != ImmersedFlowPortControlMode::FlowRate) throw std::logic_error("invalid immersed port control mode");
		const PetscInt multiplier = diagnostics_.ports.at(port_index).multiplier_row;
		if (multiplier < 0) throw std::logic_error("flow-controlled immersed port has no multiplier row");
		double lambda = 0.0; { VecReadArray state(state_, "VecGetArrayRead port multiplier"); lambda = PetscRealPart(state.Data()[multiplier]); state.Close("VecRestoreArrayRead port multiplier"); }
		for (std::size_t row = 0; row < rows.size(); ++row) {
			const double coefficient = PetscRealPart(local.flow_coefficient[row]);
			if (coefficient == 0.0) continue;
			// R_u += lambda*c and J_{u,lambda}=c.  rhs=-R.
			Check(MatSetValue(jacobian_, rows[row], multiplier, coefficient, ADD_VALUES), "MatSetValue flow controller column");
			Check(MatSetValue(jacobian_, multiplier, rows[row], coefficient, ADD_VALUES), "MatSetValue flow controller row");
				const double momentum_residual = -CheckedImmersedFlowPortProduct(lambda, coefficient, "controller lambda coefficient");
				Check(VecSetValue(rhs_, rows[row], momentum_residual, ADD_VALUES), "VecSetValue flow controller momentum residual");
		}
	}
	void MeasurePorts()
	{
		for (std::size_t port = 0; port < options_.ports.size(); ++port) {
			ImmersedFlowPortMeasurement total{};
			for (std::uint64_t id = 0; id < surface_.Cells().size(); ++id) {
				const auto& rule = surface_.UsableRule(domain_, id); if (!RuleHasLabel(rule, options_.ports[port].boundary_label)) continue;
				const auto local = MeasureImmersedFlowPortElement(domain_.Background().MaterializeElement(id), rule,
					options_.ports[port].boundary_label, Gather(domain_.Background().MaterializeElement(id)), options_.parameters.dynamic_viscosity);
				AddImmersedFlowPortFinite(total.area_m2, local.area_m2, "runtime measurement area");
				AddImmersedFlowPortFinite(total.outward_flow_m3_s, local.outward_flow_m3_s, "runtime measurement flow");
				AddImmersedFlowPortFinite(total.mean_pressure_pa, CheckedImmersedFlowPortProduct(local.mean_pressure_pa, local.area_m2, "runtime measurement pressure"), "runtime measurement pressure accumulation");
				AddImmersedFlowPortFinite(total.mean_normal_traction_pa, CheckedImmersedFlowPortProduct(local.mean_normal_traction_pa, local.area_m2, "runtime measurement traction"), "runtime measurement traction accumulation");
				AddImmersedFlowPortFinite(total.mean_velocity_squared_m2_s2, CheckedImmersedFlowPortProduct(local.mean_velocity_squared_m2_s2, local.area_m2, "runtime measurement kinetic"), "runtime measurement kinetic accumulation");
			}
			if (!(total.area_m2 > 0.0) || !std::isfinite(total.area_m2)) throw std::runtime_error("immersed flow port measurement area is invalid");
			total.mean_pressure_pa /= total.area_m2; total.mean_normal_traction_pa /= total.area_m2; total.mean_velocity_squared_m2_s2 /= total.area_m2;
			if (!std::isfinite(total.outward_flow_m3_s) || !std::isfinite(total.mean_pressure_pa)
				|| !std::isfinite(total.mean_normal_traction_pa) || !std::isfinite(total.mean_velocity_squared_m2_s2)) throw std::overflow_error("immersed flow port runtime measurement is not finite");
			auto& diagnostic = diagnostics_.ports[port]; diagnostic.measurement = total;
			if (options_.ports[port].control_mode == ImmersedFlowPortControlMode::FlowRate) {
				diagnostic.constraint_residual = options_.ports[port].value-total.outward_flow_m3_s;
				if (!std::isfinite(diagnostic.constraint_residual)) throw std::overflow_error("immersed flow controller residual overflows");
				diagnostic.absolute_flow_residual_m3_s = std::abs(diagnostic.constraint_residual);
				diagnostic.flow_tolerance_m3_s = options_.flow_controller_absolute_tolerance_m3_s;
				AddImmersedFlowPortFinite(diagnostic.flow_tolerance_m3_s, CheckedImmersedFlowPortProduct(options_.flow_controller_relative_tolerance, std::max(std::abs(options_.ports[port].value), options_.flow_controller_reference_flow_m3_s), "controller tolerance"), "controller tolerance accumulation");
				if (!std::isfinite(diagnostic.constraint_residual) || !std::isfinite(diagnostic.flow_tolerance_m3_s)
					|| !(diagnostic.flow_tolerance_m3_s > 0.0)) throw std::overflow_error("immersed flow controller tolerance is invalid");
				diagnostic.normalized_flow_residual = diagnostic.absolute_flow_residual_m3_s/diagnostic.flow_tolerance_m3_s;
				if (!std::isfinite(diagnostic.normalized_flow_residual)) throw std::overflow_error("immersed normalized flow controller residual is not finite");
				VecReadArray state(state_, "VecGetArrayRead port multiplier diagnostic"); diagnostic.multiplier = PetscRealPart(state.Data()[diagnostic.multiplier_row]); state.Close("VecRestoreArrayRead port multiplier diagnostic");
				Check(VecSetValue(rhs_, diagnostic.multiplier_row, diagnostic.constraint_residual, ADD_VALUES), "VecSetValue flow controller residual");
			}
		}
	}
	bool ControllerConvergenceSatisfied() const
	{
		for (const auto& port : diagnostics_.ports) if (port.control_mode == ImmersedFlowPortControlMode::FlowRate)
			if (!std::isfinite(port.absolute_flow_residual_m3_s) || !std::isfinite(port.flow_tolerance_m3_s)
				|| port.absolute_flow_residual_m3_s > port.flow_tolerance_m3_s) return false;
		return true;
	}
	void BuildGaugeWeights()
	{
		gauge_weights_.assign(diagnostics_.physical_dofs/4, 0.0); diagnostics_.pressure_measure = 0.0;
		for (std::uint64_t id = 0; id < domain_.Cells().size(); ++id) if (UsablePositive(id)) {
			const auto element = domain_.Background().MaterializeElement(id);
			const auto add = [&](const VolumeQuadraturePoint& point) { const auto basis = EvaluateBasis(element, point.parametric[0], point.parametric[1], point.parametric[2], false); const double measure = point.weight*basis.raw_determinant; if (!std::isfinite(measure) || !(measure > 0.0)) throw std::runtime_error("gauge volume measure is invalid"); diagnostics_.pressure_measure += measure; for (std::size_t a = 0; a < element.connectivity.size(); ++a) gauge_weights_[static_cast<std::size_t>(node_to_active_[static_cast<std::size_t>(element.connectivity[a])])] += basis.value[a]*measure; };
			if (volume_.StorageMode() == CutCellVolumeQuadratureStorageMode::Expanded) for (const auto& point : volume_.UsableRule(domain_, id).Points()) add(point);
			else ForEachVolumePoint(volume_.UsableCompactRule(domain_, id), add);
		}
	}
	void InsertGauge()
	{
		if (!HasGauge()) throw std::logic_error("cannot insert an absent pressure gauge");
		if (!std::isfinite(diagnostics_.pressure_measure) || !(diagnostics_.pressure_measure > 0.0)) throw std::runtime_error("pressure gauge has zero measure");
		double pressure_integral = 0.0, lambda = 0.0;
		{
			VecReadArray state(state_, "VecGetArrayRead gauge");
			const PetscInt q = GaugeDof(); lambda = PetscRealPart(state.Data()[q]);
			for (std::size_t node = 0; node < gauge_weights_.size(); ++node) {
				const double weight = gauge_weights_[node]; if (!std::isfinite(weight)) throw std::runtime_error("pressure gauge entry is not finite");
				pressure_integral += weight*PetscRealPart(state.Data()[4*CheckedPetscCount(node)+3]);
			}
			state.Close("VecRestoreArrayRead gauge");
		}
		if (!std::isfinite(pressure_integral)) throw std::runtime_error("pressure gauge residual is not finite");
		const PetscInt q = GaugeDof();
		// Constructor seeding reserved q,q structurally.  Do not insert a
		// numerical zero here: the true gauge diagonal remains exactly zero.
		for (std::size_t node = 0; node < gauge_weights_.size(); ++node) {
			const double weight = gauge_weights_[node]; const PetscInt p = 4*CheckedPetscCount(node)+3;
			Check(MatSetValue(jacobian_, p, q, weight, ADD_VALUES), "MatSetValue gauge column"); Check(MatSetValue(jacobian_, q, p, weight, ADD_VALUES), "MatSetValue gauge row");
				Check(VecSetValue(rhs_, p, -CheckedImmersedFlowPortProduct(weight, lambda, "gauge multiplier weight"), ADD_VALUES), "VecSetValue gauge pressure residual");
		}
		// R_g = g^T p and rhs stores -R.
		Check(VecSetValue(rhs_, q, -pressure_integral, ADD_VALUES), "VecSetValue gauge residual");
	}
	void MeasureConstantPressureDefect()
	{
		Check(MatMult(jacobian_, constant_pressure_, pressure_defect_), "MatMult pressure defect");
		PetscReal norm = 0.0; Check(VecNorm(pressure_defect_, NORM_INFINITY, &norm), "VecNorm pressure defect"); diagnostics_.constant_pressure_defect = norm;
	}
	std::array<double, 3> ResidualBlockNorms(const std::vector<PetscScalar>& residual) const
	{
		if (residual.size() != diagnostics_.total_dofs) throw std::invalid_argument("residual block vector size is invalid");
		std::array<double, 3> norms{{0.0,0.0,0.0}};
		for (std::size_t row = 0; row < residual.size(); ++row) {
			const double value = PetscRealPart(residual[row]);
			const std::size_t block = row >= diagnostics_.physical_dofs ? 2 : ((row%4) == 3 ? 1 : 0);
			norms[block] += value*value;
		}
		for (auto& value : norms) value = std::sqrt(value);
		return norms;
	}
	bool BlockReductionSatisfied(const std::array<double, 3>& initial, const std::vector<PetscScalar>& residual, double global_initial) const
	{
		if (options_.nonlinear_block_reduction == 0.0) return true;
		const auto current = ResidualBlockNorms(residual);
		const double zero_block_tolerance = options_.nonlinear_absolute_tolerance+options_.nonlinear_relative_tolerance*global_initial;
		for (std::size_t block = 0; block < current.size(); ++block)
			if (initial[block] > 0.0 ? current[block] > initial[block]/options_.nonlinear_block_reduction : current[block] > zero_block_tolerance) return false;
		return true;
	}
	void EnsureFinite(Vec vector) const { const auto values = CopyVector(vector); for (const auto value : values) if (!std::isfinite(PetscRealPart(value))) throw std::runtime_error("PETSc linear update is not finite"); }
	static std::size_t WallContributionCount(const ImmersedNitscheWallDiagnostics& diagnostics) { std::size_t total = 0; for (const auto& entry : diagnostics.by_boundary_id) total += entry.second.selected_points; return total; }
	std::vector<PetscScalar> CopyVector(Vec vector) const { PetscInt n = 0; Check(VecGetSize(vector, &n), "VecGetSize"); VecReadArray data(vector, "VecGetArrayRead copy"); std::vector<PetscScalar> r(data.Data(), data.Data()+n); data.Close("VecRestoreArrayRead copy"); return r; }

	const CartesianDomainClassification& domain_; const CutCellVolumeQuadratureCatalog& volume_; const ImmersedSurfaceQuadratureCatalog& surface_; const CutCellGhostPenaltyCatalog& ghost_; ImmersedStaticFlowOptions options_;
	std::vector<std::int32_t> active_nodes_; std::vector<PetscInt> node_to_active_; std::vector<double> gauge_weights_;
	Mat jacobian_ = nullptr; Vec state_ = nullptr, committed_ = nullptr, prepared_ = nullptr, rhs_ = nullptr, update_ = nullptr;
	mutable Vec action_input_ = nullptr, action_output_ = nullptr;
	Vec constant_pressure_ = nullptr, pressure_defect_ = nullptr;
	KSP ksp_ = nullptr; ImmersedStaticFlowDiagnostics diagnostics_{};
	bool fail_next_prepare_for_testing_ = false;
};

} // namespace iga

#endif
