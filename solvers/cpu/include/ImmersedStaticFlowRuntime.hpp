#ifndef IGA_IMMERSED_STATIC_FLOW_RUNTIME_HPP
#define IGA_IMMERSED_STATIC_FLOW_RUNTIME_HPP

// Serial-only Phase 5 global assembly.  This is intentionally a small
// orchestration layer: geometry catalogs retain rules, while this class owns
// active-row compaction, PETSc insertion, the pressure gauge, and the
// transactional nonlinear trial state.
#include "CutCellGhostPenalty.hpp"
#include "ImmersedNitscheWall.hpp"

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

struct ImmersedStaticFlowOptions {
	NavierStokesParameters parameters{1.0, 1.0, 0.0};
	std::vector<int> wall_labels;
	double wall_gamma0 = 2.0;
	NavierStokesBodyForceEvaluator body_force = [](const std::array<double, 3>&) {
		return std::array<double, 3>{{0.0,0.0,0.0}};
	};
	PetscInt nonlinear_maximum_iterations = 12;
	PetscInt ksp_maximum_iterations = 2000;
	double nonlinear_relative_tolerance = 1e-9;
	double nonlinear_absolute_tolerance = 1e-11;
	double minimum_damping = 1.0/128.0;
};

struct ImmersedStaticFlowDiagnostics {
	std::size_t active_nodes = 0, physical_dofs = 0, total_dofs = 0;
	std::size_t volume_cells = 0, surface_cells = 0, ghost_faces = 0;
	double pressure_measure = 0.0, constant_pressure_defect = 0.0;
	PetscInt nonlinear_iterations = 0, ksp_iterations = 0;
	KSPConvergedReason ksp_reason = KSP_CONVERGED_ITERATING;
	double residual_norm = 0.0, damping = 0.0;
	bool committed = true, trial_active = false, converged = false;
	std::size_t commit_count = 0, rollback_count = 0;
};

class ImmersedStaticFlowRuntime {
public:
	ImmersedStaticFlowRuntime(const CartesianDomainClassification& domain,
		const CutCellVolumeQuadratureCatalog& volume,
		const ImmersedSurfaceQuadratureCatalog& surface,
		const CutCellGhostPenaltyCatalog& ghost, ImmersedStaticFlowOptions options = {})
		: domain_(domain), volume_(volume), surface_(surface), ghost_(ghost), options_(std::move(options))
	{
		if (!options_.body_force || !std::isfinite(options_.parameters.density) || !(options_.parameters.density > 0.0)
			|| !std::isfinite(options_.parameters.dynamic_viscosity) || !(options_.parameters.dynamic_viscosity > 0.0)
			|| !std::isfinite(options_.parameters.dt) || options_.parameters.dt != 0.0 || !std::isfinite(options_.wall_gamma0)
			|| !(options_.wall_gamma0 > 0.0) || options_.nonlinear_maximum_iterations <= 0
			|| options_.ksp_maximum_iterations <= 0 || !std::isfinite(options_.nonlinear_relative_tolerance) || !(options_.nonlinear_relative_tolerance > 0.0)
			|| !std::isfinite(options_.nonlinear_absolute_tolerance) || !(options_.nonlinear_absolute_tolerance > 0.0) || !std::isfinite(options_.minimum_damping) || !(options_.minimum_damping > 0.0)
			|| options_.minimum_damping > 1.0)
			throw std::invalid_argument("immersed static-flow options are invalid");
		ghost_.ValidateBinding(domain_, volume_);
		if (!surface_.Usable() || !ghost_.Usable()) throw std::invalid_argument("immersed static-flow catalogs are unusable");
		ValidateImmersedNitscheWallLabels(options_.wall_labels);
		if (options_.wall_labels.empty()) throw std::invalid_argument("immersed static-flow requires at least one wall label");
		if (surface_.GridSpec().lower_m != volume_.GridSpec().lower_m || surface_.GridSpec().upper_m != volume_.GridSpec().upper_m
			|| surface_.GridSpec().cells != volume_.GridSpec().cells || surface_.SurfaceCanonicalHash() != volume_.SurfaceCanonicalHash()
			|| surface_.Cells().size() != volume_.Cells().size()) throw std::invalid_argument("immersed surface and volume catalogs do not have an exact common binding");
		PreflightCatalogs();
		BuildActiveMap();
		const PetscInt n = CheckedPetscCount(diagnostics_.total_dofs);
		try {
			Check(MatCreateSeqAIJ(PETSC_COMM_SELF, n, n, 300, nullptr, &jacobian_), "MatCreateSeqAIJ");
			Check(MatSetOption(jacobian_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE), "MatSetOption");
			Check(VecCreateSeq(PETSC_COMM_SELF, n, &state_), "VecCreateSeq");
			Check(VecDuplicate(state_, &committed_), "VecDuplicate committed"); Check(VecDuplicate(state_, &rhs_), "VecDuplicate rhs"); Check(VecDuplicate(state_, &update_), "VecDuplicate update");
			Check(VecSet(state_, 0.0), "VecSet state"); Check(VecSet(committed_, 0.0), "VecSet committed");
			Check(KSPCreate(PETSC_COMM_SELF, &ksp_), "KSPCreate"); Check(KSPSetType(ksp_, KSPGMRES), "KSPSetType");
			Check(KSPSetTolerances(ksp_, 1e-10, PETSC_DEFAULT, PETSC_DEFAULT, options_.ksp_maximum_iterations), "KSPSetTolerances");
			PC pc = nullptr; Check(KSPGetPC(ksp_, &pc), "KSPGetPC"); Check(PCSetType(pc, PCLU), "PCSetType");
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
	PetscInt GaugeDof() const noexcept { return static_cast<PetscInt>(diagnostics_.physical_dofs); }
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
		Vec input = nullptr, output = nullptr;
		try { Check(VecDuplicate(state_, &input), "VecDuplicate Jacobian input"); Check(VecDuplicate(state_, &output), "VecDuplicate Jacobian output");
			PetscScalar* data = nullptr; Check(VecGetArray(input, &data), "VecGetArray Jacobian input"); for (std::size_t i = 0; i < values.size(); ++i) data[i] = values[i]; Check(VecRestoreArray(input, &data), "VecRestoreArray Jacobian input"); Check(MatMult(jacobian_, input, output), "MatMult Jacobian action"); auto result = CopyVector(output); Check(VecDestroy(&output), "VecDestroy Jacobian output"); Check(VecDestroy(&input), "VecDestroy Jacobian input"); return result;
		} catch (...) { if (output) VecDestroy(&output); if (input) VecDestroy(&input); throw; }
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
		Check(MatZeroEntries(jacobian_), "MatZeroEntries"); Check(VecSet(rhs_, 0.0), "VecSet rhs");
		diagnostics_.volume_cells = diagnostics_.surface_cells = diagnostics_.ghost_faces = 0;
		diagnostics_.pressure_measure = 0.0;
		for (std::uint64_t id = 0; id < domain_.Cells().size(); ++id) {
			if (!UsablePositive(id)) continue;
			const auto element = domain_.Background().MaterializeElement(id);
			const auto nodal = Gather(element);
			NavierStokesSystem local = BuildVolumeSystem(element, nodal, id);
			++diagnostics_.volume_cells;
			if (domain_.Cells()[static_cast<std::size_t>(id)].classification == CellClassification::Cut) {
				if (!ghost_.Covered(id)) throw std::runtime_error("covered immersed wall and ghost policy mismatch");
				const auto wall = BuildImmersedNitscheWallElementFromVolumeSystem(domain_, volume_, surface_, id, nodal, {},
					options_.parameters, options_.wall_labels, local, ghost_, options_.wall_gamma0);
				local = wall.system;
				if (WallContributionCount(wall.diagnostics) == 0) throw std::runtime_error("selected immersed wall labels produced no surface contribution");
				++diagnostics_.surface_cells;
			}
			ScatterElement(element.connectivity, local);
		}
		for (std::size_t face = 0; face < ghost_.Faces().size(); ++face) {
			VecReadArray values(state_, "VecGetArrayRead ghost face");
			const auto block = ghost_.AssembleFaceLocal(face, domain_, volume_, [this, data = values.Data()](std::int32_t node, int field) { return PetscRealPart(data[Dof(node, field)]); }, options_.parameters.dynamic_viscosity);
			values.Close("VecRestoreArrayRead ghost face");
			ScatterBlock(block.connectivity, block.jacobian, block.negative_residual); ++diagnostics_.ghost_faces;
		}
		Check(MatAssemblyBegin(jacobian_, MAT_FINAL_ASSEMBLY), "MatAssemblyBegin physical"); Check(MatAssemblyEnd(jacobian_, MAT_FINAL_ASSEMBLY), "MatAssemblyEnd physical");
		Check(VecAssemblyBegin(rhs_), "VecAssemblyBegin physical"); Check(VecAssemblyEnd(rhs_), "VecAssemblyEnd physical");
		MeasureConstantPressureDefect();
		BuildGaugeWeights();
		InsertGauge();
		Check(MatAssemblyBegin(jacobian_, MAT_FINAL_ASSEMBLY), "MatAssemblyBegin gauge"); Check(MatAssemblyEnd(jacobian_, MAT_FINAL_ASSEMBLY), "MatAssemblyEnd gauge");
		Check(VecAssemblyBegin(rhs_), "VecAssemblyBegin gauge"); Check(VecAssemblyEnd(rhs_), "VecAssemblyEnd gauge");
	}

	// Start each Newton attempt from an exact copy of the committed vector.
	bool SolveTrial()
	{
		Check(VecCopy(committed_, state_), "VecCopy begin trial"); diagnostics_.trial_active = true; diagnostics_.committed = false;
		diagnostics_.converged = false; diagnostics_.nonlinear_iterations = diagnostics_.ksp_iterations = 0;
		diagnostics_.ksp_reason = KSP_CONVERGED_ITERATING; diagnostics_.damping = 0.0; diagnostics_.residual_norm = 0.0;
		try {
			double initial = -1.0;
			for (PetscInt iteration = 0; iteration < options_.nonlinear_maximum_iterations; ++iteration) {
				Assemble(); PetscReal residual = 0.0; Check(VecNorm(rhs_, NORM_2, &residual), "VecNorm residual");
				if (!std::isfinite(residual)) throw std::runtime_error("nonlinear residual is not finite");
				if (initial < 0.0) initial = residual;
				diagnostics_.residual_norm = residual; diagnostics_.nonlinear_iterations = iteration;
				if (residual <= std::max(options_.nonlinear_absolute_tolerance, options_.nonlinear_relative_tolerance*initial)) {
					diagnostics_.converged = true; return true;
				}
				Check(KSPSetOperators(ksp_, jacobian_, jacobian_), "KSPSetOperators"); Check(KSPSolve(ksp_, rhs_, update_), "KSPSolve");
				Check(KSPGetConvergedReason(ksp_, &diagnostics_.ksp_reason), "KSPGetConvergedReason"); PetscInt work = 0; Check(KSPGetIterationNumber(ksp_, &work), "KSPGetIterationNumber"); diagnostics_.ksp_iterations += work;
				if (diagnostics_.ksp_reason <= 0) throw std::runtime_error("static immersed-flow KSP failed");
				EnsureFinite(update_);
				double damping = 1.0, old = residual; bool accepted = false;
				while (damping >= options_.minimum_damping) {
					Check(VecAXPY(state_, damping, update_), "VecAXPY trial update"); Assemble(); PetscReal candidate = 0.0; Check(VecNorm(rhs_, NORM_2, &candidate), "VecNorm candidate");
				if (std::isfinite(candidate) && candidate < old) {
					accepted = true; diagnostics_.damping = damping; diagnostics_.residual_norm = candidate;
					++diagnostics_.nonlinear_iterations;
					if (candidate <= std::max(options_.nonlinear_absolute_tolerance, options_.nonlinear_relative_tolerance*initial)) {
						diagnostics_.converged = true; return true;
					}
					break;
				}
					Check(VecAXPY(state_, -damping, update_), "VecAXPY trial undo"); damping *= .5;
				}
				if (!accepted) throw std::runtime_error("static immersed-flow backtracking failed");
			}
			throw std::runtime_error("static immersed-flow nonlinear iteration cap reached");
		} catch (...) { Rollback(); throw; }
	}
	void Commit()
	{
		if (!diagnostics_.trial_active || !diagnostics_.converged) throw std::logic_error("cannot commit an unconverged static-flow trial");
		Check(VecCopy(state_, committed_), "VecCopy commit"); diagnostics_.trial_active = false; diagnostics_.committed = true; ++diagnostics_.commit_count;
	}
	void Rollback()
	{
		if (!diagnostics_.trial_active) return;
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
		if (rhs_) VecDestroy(&rhs_);
		if (committed_) VecDestroy(&committed_);
		if (state_) VecDestroy(&state_);
		if (jacobian_) MatDestroy(&jacobian_);
	}
	static PetscInt CheckedPetscCount(std::size_t value)
	{
		if (value > static_cast<std::size_t>(std::numeric_limits<PetscInt>::max())) throw std::overflow_error("PETSc row count overflows");
		return static_cast<PetscInt>(value);
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
			bool selected = false;
			for (const auto& point : surface_.UsableRule(domain_, id).Points())
				selected = selected || std::binary_search(options_.wall_labels.begin(), options_.wall_labels.end(), point.boundary_id);
			if (!selected) throw std::runtime_error("positive cut cell has no selected immersed wall surface contribution");
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
		if (active_nodes_.empty() || active_nodes_.size() > std::numeric_limits<std::size_t>::max()/4) throw std::runtime_error("active-node compaction is empty or overflows");
		for (std::size_t i = 0; i < active_nodes_.size(); ++i) { const PetscInt index = CheckedPetscCount(i); node_to_active_[static_cast<std::size_t>(active_nodes_[i])] = index; }
		diagnostics_.active_nodes = active_nodes_.size(); diagnostics_.physical_dofs = 4*active_nodes_.size();
		if (diagnostics_.physical_dofs == std::numeric_limits<std::size_t>::max()) throw std::overflow_error("gauge dof overflows");
		diagnostics_.total_dofs = diagnostics_.physical_dofs+1;
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
		Check(MatSetValues(jacobian_, CheckedPetscCount(rows.size()), rows.data(), CheckedPetscCount(rows.size()), rows.data(), system.jacobian.data(), ADD_VALUES), "MatSetValues element");
		Check(VecSetValues(rhs_, CheckedPetscCount(rows.size()), rows.data(), system.negative_residual.data(), ADD_VALUES), "VecSetValues element");
	}
	void ScatterBlock(const std::vector<std::int32_t>& nodes, const std::vector<PetscScalar>& matrix, const std::vector<PetscScalar>& residual)
	{
		std::vector<PetscInt> rows; rows.reserve(4*nodes.size()); for (const auto node : nodes) for (int f = 0; f < 4; ++f) rows.push_back(Dof(node, f));
		Check(MatSetValues(jacobian_, CheckedPetscCount(rows.size()), rows.data(), CheckedPetscCount(rows.size()), rows.data(), matrix.data(), ADD_VALUES), "MatSetValues block"); Check(VecSetValues(rhs_, CheckedPetscCount(rows.size()), rows.data(), residual.data(), ADD_VALUES), "VecSetValues block");
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
		for (std::size_t node = 0; node < gauge_weights_.size(); ++node) {
			const double weight = gauge_weights_[node]; const PetscInt p = 4*CheckedPetscCount(node)+3;
			Check(MatSetValue(jacobian_, p, q, weight, ADD_VALUES), "MatSetValue gauge column"); Check(MatSetValue(jacobian_, q, p, weight, ADD_VALUES), "MatSetValue gauge row");
			Check(VecSetValue(rhs_, p, -weight*lambda, ADD_VALUES), "VecSetValue gauge pressure residual");
		}
		// R_g = g^T p and rhs stores -R.
		Check(VecSetValue(rhs_, q, -pressure_integral, ADD_VALUES), "VecSetValue gauge residual");
	}
	void MeasureConstantPressureDefect()
	{
		Vec constant = nullptr, defect = nullptr; try { Check(VecDuplicate(state_, &constant), "VecDuplicate pressure constant"); Check(VecDuplicate(state_, &defect), "VecDuplicate pressure defect"); Check(VecSet(constant, 0.0), "VecSet pressure constant");
			for (std::size_t a = 0; a < active_nodes_.size(); ++a) Check(VecSetValue(constant, 4*CheckedPetscCount(a)+3, 1.0, INSERT_VALUES), "VecSetValue pressure constant");
			Check(VecAssemblyBegin(constant), "VecAssemblyBegin pressure constant"); Check(VecAssemblyEnd(constant), "VecAssemblyEnd pressure constant"); Check(MatMult(jacobian_, constant, defect), "MatMult pressure defect"); PetscReal norm = 0.0; Check(VecNorm(defect, NORM_INFINITY, &norm), "VecNorm pressure defect"); diagnostics_.constant_pressure_defect = norm; Check(VecDestroy(&defect), "VecDestroy pressure defect"); Check(VecDestroy(&constant), "VecDestroy pressure constant");
		} catch (...) { if (defect) VecDestroy(&defect); if (constant) VecDestroy(&constant); throw; }
	}
	void EnsureFinite(Vec vector) const { const auto values = CopyVector(vector); for (const auto value : values) if (!std::isfinite(PetscRealPart(value))) throw std::runtime_error("PETSc linear update is not finite"); }
	static std::size_t WallContributionCount(const ImmersedNitscheWallDiagnostics& diagnostics) { std::size_t total = 0; for (const auto& entry : diagnostics.by_boundary_id) total += entry.second.selected_points; return total; }
	std::vector<PetscScalar> CopyVector(Vec vector) const { PetscInt n = 0; Check(VecGetSize(vector, &n), "VecGetSize"); VecReadArray data(vector, "VecGetArrayRead copy"); std::vector<PetscScalar> r(data.Data(), data.Data()+n); data.Close("VecRestoreArrayRead copy"); return r; }

	const CartesianDomainClassification& domain_; const CutCellVolumeQuadratureCatalog& volume_; const ImmersedSurfaceQuadratureCatalog& surface_; const CutCellGhostPenaltyCatalog& ghost_; ImmersedStaticFlowOptions options_;
	std::vector<std::int32_t> active_nodes_; std::vector<PetscInt> node_to_active_; std::vector<double> gauge_weights_;
	Mat jacobian_ = nullptr; Vec state_ = nullptr, committed_ = nullptr, rhs_ = nullptr, update_ = nullptr; KSP ksp_ = nullptr; ImmersedStaticFlowDiagnostics diagnostics_{};
};

} // namespace iga

#endif
