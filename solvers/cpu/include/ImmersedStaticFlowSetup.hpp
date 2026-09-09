#ifndef IGA_IMMERSED_STATIC_FLOW_SETUP_HPP
#define IGA_IMMERSED_STATIC_FLOW_SETUP_HPP

#include "CutCellGhostPenalty.hpp"
#include "ImmersedFlowPort.hpp"
#include "ImmersedNitscheWall.hpp"
#include "PhaseProfile.hpp"
#include <petscksp.h>
#include <algorithm>
#include <array>
#include <limits>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <map>
#include <string>
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
	// Wall-clock observations are diagnostics only.  Aggregates cover this
	// runtime object's lifetime and intentionally do not participate in any
	// nonlinear, controller, or acceptance criterion.
	double last_assembly_seconds = 0.0, aggregate_assembly_seconds = 0.0;
	double last_linear_solve_seconds = 0.0, aggregate_linear_solve_seconds = 0.0;
	double last_line_search_candidate_assembly_seconds = 0.0,
		aggregate_line_search_candidate_assembly_seconds = 0.0;
	std::vector<ImmersedStaticFlowNewtonStep> newton_steps;
};

// These are trace and volume diagnostics only: they neither enter the
// controller constraint nor alter the assembled continuity equations.  They
// deliberately stream the exact catalog rules used by production assembly so
// a cut-domain divergence-theorem discrepancy is observable without mixing in
// full-cell or independently reconstructed quadrature.
struct ImmersedStaticFlowConservationDiagnostics {
	std::map<int, double> surface_flow_by_boundary_label_m3_s;
	double open_port_outward_flow_m3_s = 0.0;
	double wall_outward_flow_m3_s = 0.0;
	double total_surface_outward_flow_m3_s = 0.0;
	double volume_divergence_integral_m3_s = 0.0;
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

// Shared serial/MPI geometry preflight and stable active/scalar row numbering.
// Construction is local: this object owns no distributed PETSc resource.
class ImmersedStaticFlowSetup {
public:
	ImmersedStaticFlowSetup(const CartesianDomainClassification& domain,
		const CutCellVolumeQuadratureCatalog& volume,
		const ImmersedSurfaceQuadratureCatalog& surface,
		const CutCellGhostPenaltyCatalog& ghost, ImmersedStaticFlowOptions options)
		: domain_(domain), volume_(volume), surface_(surface), ghost_(ghost), options_(std::move(options))
	{
		PhaseScope geometry_phase(ProfilePhase::Geometry);
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
	}
	const ImmersedStaticFlowOptions& Options() const noexcept { return options_; }
	const ImmersedStaticFlowDiagnostics& Diagnostics() const noexcept { return diagnostics_; }
	const std::vector<std::int32_t>& ActiveNodes() const noexcept { return active_nodes_; }
	bool HasGauge() const noexcept { return diagnostics_.gauge_present; }
	PetscInt GaugeDof() const noexcept { return diagnostics_.gauge_row; }
	PetscInt Dof(std::int32_t node, int field) const
	{
		if (field < 0 || field > 3 || node < 0 || static_cast<std::size_t>(node) >= node_to_active_.size()
			|| node_to_active_[static_cast<std::size_t>(node)] < 0) throw std::out_of_range("background node is not active");
		return 4*node_to_active_[static_cast<std::size_t>(node)]+field;
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
protected:
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

protected:
	static PetscInt CheckedPetscCount(std::size_t value)
	{
		if (value > static_cast<std::size_t>(std::numeric_limits<PetscInt>::max())) throw std::overflow_error("PETSc row count overflows");
		return static_cast<PetscInt>(value);
	}
	const CartesianDomainClassification& domain_;
	const CutCellVolumeQuadratureCatalog& volume_;
	const ImmersedSurfaceQuadratureCatalog& surface_;
	const CutCellGhostPenaltyCatalog& ghost_;
	ImmersedStaticFlowOptions options_;
	std::vector<std::int32_t> active_nodes_;
	std::vector<PetscInt> node_to_active_;
	ImmersedStaticFlowDiagnostics diagnostics_{};
};

} // namespace iga
#endif
