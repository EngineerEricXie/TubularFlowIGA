#ifndef IGA_IMMERSED_NITSCHE_WALL_HPP
#define IGA_IMMERSED_NITSCHE_WALL_HPP

// Catalog-bound static immersed-wall assembly.  Surface weights are physical
// dA values supplied by ImmersedSurfaceQuadratureCatalog; this deliberately
// does not use body-fitted face utilities or a surface Jacobian.
#include "CutCellVolumeQuadrature.hpp"
#include "CutCellGhostPenalty.hpp"
#include "ImmersedSurfaceQuadrature.hpp"
#include "NavierStokesElement.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace iga {

using ImmersedWallVelocityEvaluator = std::function<std::array<double, 3>(
	const std::array<double, 3>& physical, int boundary_id)>;
using ImmersedMaterialWallVelocityEvaluator = std::function<std::array<double, 3>(
	const SurfaceQuadraturePoint& point, const ImmersedSurfaceQuadraturePointProvenance& provenance)>;

inline ImmersedMaterialWallVelocityEvaluator AdaptImmersedWallVelocityEvaluator(
	const ImmersedWallVelocityEvaluator& evaluator)
{
	if (!evaluator) return {};
	return [&evaluator](const SurfaceQuadraturePoint& point,
		const ImmersedSurfaceQuadraturePointProvenance&) { return evaluator(point.physical, point.boundary_id); };
}

struct ImmersedNitscheWallLabelDiagnostics {
	std::size_t selected_points = 0;
	std::size_t skipped_points = 0;
	double selected_area_m2 = 0.0;
	double skipped_area_m2 = 0.0;
};

struct ImmersedNitscheWallDiagnostics {
	std::map<int, ImmersedNitscheWallLabelDiagnostics> by_boundary_id;
	double fraction_lower = 0.0;
	double fraction_estimate = 0.0;
	double fraction_upper = 0.0;
	double minimum_h_n_m = 0.0;
	double maximum_h_n_m = 0.0;
	double minimum_eta = 0.0;
	double maximum_eta = 0.0;
	double maximum_eta_h_n_over_mu = 0.0;
	double maximum_gap_norm = 0.0;
	bool ghost_covered_policy = false;
};

struct ImmersedNitscheWallAssembly {
	NavierStokesSystem system;
	ImmersedNitscheWallDiagnostics diagnostics;
};

inline void ValidateImmersedNitscheWallLabels(const std::vector<int>& labels)
{
	for (std::size_t i = 0; i < labels.size(); ++i) {
		if (labels[i] < 0) throw std::invalid_argument("immersed Nitsche wall boundary labels must be nonnegative");
		if (i && labels[i-1] >= labels[i])
			throw std::invalid_argument("immersed Nitsche wall boundary labels must be sorted and unique");
	}
}

inline void ValidateImmersedNitscheWallProvenance(const CartesianDomainClassification& domain,
	const ImmersedSurfaceQuadraturePointProvenance& provenance)
{
	if (provenance.canonical_triangle >= domain.SurfaceIndex().Surface().Triangles().size())
		throw std::invalid_argument("immersed Nitsche wall provenance is not bound to this domain surface");
	const double tolerance = 64.0*std::numeric_limits<double>::epsilon();
	double barycentric_sum = 0.0;
	for (const double value : provenance.canonical_barycentric) {
		if (!std::isfinite(value) || value < -tolerance || value > 1.0+tolerance)
			throw std::invalid_argument("immersed Nitsche wall provenance is outside the simplex");
		barycentric_sum += value;
	}
	if (!std::isfinite(barycentric_sum) || std::abs(barycentric_sum-1.0) > tolerance)
		throw std::invalid_argument("immersed Nitsche wall provenance does not sum to one");
}

inline double ImmersedNitscheWallFraction(const CutCellVolumeQuadratureCatalog& volume_catalog,
	std::uint64_t cell_id, ImmersedNitscheWallDiagnostics& diagnostics)
{
	const auto& volume = volume_catalog.Cell(cell_id).diagnostics;
	diagnostics.fraction_lower = volume.lower_reference_volume;
	diagnostics.fraction_estimate = volume.estimated_reference_volume;
	diagnostics.fraction_upper = volume.upper_reference_volume;
	const auto finite_fraction = [](double value) { return std::isfinite(value) && value >= 0.0 && value <= 1.0; };
	const double tolerance = 2.0e-12*std::max({1.0, std::abs(diagnostics.fraction_lower),
		std::abs(diagnostics.fraction_estimate), std::abs(diagnostics.fraction_upper)});
	if (!finite_fraction(diagnostics.fraction_lower) || !finite_fraction(diagnostics.fraction_estimate)
		|| !finite_fraction(diagnostics.fraction_upper)
		|| diagnostics.fraction_lower > diagnostics.fraction_estimate+tolerance
		|| diagnostics.fraction_estimate > diagnostics.fraction_upper+tolerance)
		throw std::runtime_error("cut-cell reference fraction is not a finite certified fraction in [0,1]");
	if (!(diagnostics.fraction_estimate > 0.0))
		throw std::runtime_error("immersed Nitsche wall requires a positive cut-cell reference fraction");
	return diagnostics.fraction_estimate;
}

inline void RequireFiniteImmersedNitscheWallValue(double value, const char* description)
{
	if (!std::isfinite(value)) throw std::overflow_error(description);
}

inline double ImmersedNitscheWallDot(const std::array<double, 3>& left,
	const std::array<double, 3>& right, const char* description)
{
	const double value = left[0]*right[0]+left[1]*right[1]+left[2]*right[2];
	RequireFiniteImmersedNitscheWallValue(value, description);
	return value;
}

inline void AccumulateImmersedNitscheWallValue(PetscScalar& entry, double contribution,
	const char* description)
{
	RequireFiniteImmersedNitscheWallValue(contribution, description);
	if (!std::isfinite(PetscRealPart(entry)))
		throw std::overflow_error("immersed Nitsche wall cannot add to a nonfinite system entry");
	const PetscScalar updated = entry+contribution;
	if (!std::isfinite(PetscRealPart(updated))) throw std::overflow_error(description);
	entry = updated;
}

inline void AccumulateImmersedNitscheWallDiagnostic(double& entry, double contribution,
	const char* description)
{
	RequireFiniteImmersedNitscheWallValue(entry, description);
	RequireFiniteImmersedNitscheWallValue(contribution, description);
	const double updated = entry+contribution;
	RequireFiniteImmersedNitscheWallValue(updated, description);
	entry = updated;
}

inline void ValidateImmersedNitscheWallPreassembledVolumeSystem(const Element& element,
	const std::vector<std::array<double, 4>>& nodal_state,
	const std::vector<std::array<double, 4>>& previous_nodal_state,
	const NavierStokesParameters& parameters, const NavierStokesSystem& system)
{
	if (!std::isfinite(parameters.density) || !(parameters.density > 0.0)
		|| !std::isfinite(parameters.dynamic_viscosity) || !(parameters.dynamic_viscosity > 0.0)
		|| !std::isfinite(parameters.dt) || parameters.dt < 0.0)
		throw std::invalid_argument("immersed Nitsche wall Navier-Stokes parameters are invalid");
	const std::size_t nen = element.connectivity.size();
	if (nen > std::numeric_limits<std::size_t>::max()/4)
		throw std::overflow_error("immersed Nitsche wall local degree-of-freedom count overflows");
	const std::size_t ndof = 4*nen;
	if (ndof != 0 && ndof > std::numeric_limits<std::size_t>::max()/ndof)
		throw std::overflow_error("immersed Nitsche wall local Jacobian count overflows");
	if (nodal_state.size() != nen || (parameters.dt > 0.0 && previous_nodal_state.size() != nen))
		throw std::invalid_argument("immersed Nitsche wall preassembled state size does not match Navier-Stokes element");
	if (system.negative_residual.size() != ndof || system.jacobian.size() != ndof*ndof)
		throw std::invalid_argument("immersed Nitsche wall preassembled volume-system block size is invalid");
	for (const auto value : system.negative_residual)
		if (!std::isfinite(PetscRealPart(value))) throw std::invalid_argument("immersed Nitsche wall preassembled volume residual is not finite");
	for (const auto value : system.jacobian)
		if (!std::isfinite(PetscRealPart(value))) throw std::invalid_argument("immersed Nitsche wall preassembled volume Jacobian is not finite");
}

inline ImmersedNitscheWallAssembly BuildImmersedNitscheWallElementImpl(
	const CartesianDomainClassification& domain,
	const CutCellVolumeQuadratureCatalog& volume_catalog,
	const ImmersedSurfaceQuadratureCatalog& surface_catalog, std::uint64_t cell_id,
	const std::vector<std::array<double, 4>>& nodal_state,
	const std::vector<std::array<double, 4>>& previous_nodal_state,
	const NavierStokesParameters& parameters, const std::vector<int>& selected_boundary_labels,
	double gamma0, const ImmersedMaterialWallVelocityEvaluator& wall_velocity,
	const CutCellGhostPenaltyCatalog* ghost_catalog, const NavierStokesSystem* volume_system = nullptr)
{
	if (!std::isfinite(gamma0) || !(gamma0 > 0.0))
		throw std::invalid_argument("immersed Nitsche wall gamma0 must be finite and positive");
	if (!std::isfinite(parameters.density) || !(parameters.density > 0.0)
		|| !std::isfinite(parameters.dynamic_viscosity) || !(parameters.dynamic_viscosity > 0.0)
		|| !std::isfinite(parameters.dt) || parameters.dt < 0.0)
		throw std::invalid_argument("immersed Nitsche wall Navier-Stokes parameters are invalid");
	if (!wall_velocity) throw std::invalid_argument("immersed Nitsche wall velocity evaluator is required");
	ValidateImmersedNitscheWallLabels(selected_boundary_labels);
	if (ghost_catalog) {
		ghost_catalog->ValidateBinding(domain, volume_catalog);
		if (!ghost_catalog->Covered(cell_id))
			throw std::runtime_error("covered immersed Nitsche policy requires ghost-covered cut cell");
	}

	// Validate exact domain/grid/surface bindings.  A caller that streamed a
	// compact volume rule may supply its already assembled volume base; this
	// adapter then adds only the surface terms.
	if (volume_catalog.StorageMode() == CutCellVolumeQuadratureStorageMode::Expanded)
		volume_catalog.ValidateUsableRule(domain, cell_id);
	else
		volume_catalog.ValidateUsableCompactRule(domain, cell_id);
	const auto& surface_rule = surface_catalog.UsableRule(domain, cell_id);
	const auto& surface_provenance = surface_catalog.UsableProvenance(domain, cell_id);
	if (surface_provenance.size() != surface_rule.Points().size())
		throw std::runtime_error("immersed Nitsche wall provenance must be one-to-one with surface quadrature");
	for (std::size_t point_index = 0; point_index < surface_rule.Points().size(); ++point_index) {
		const auto& point = surface_rule.Points()[point_index];
		if (!std::isfinite(point.weight) || !(point.weight > 0.0))
			throw std::invalid_argument("immersed Nitsche wall surface quadrature weight is invalid");
		ValidateImmersedNitscheWallProvenance(domain, surface_provenance[point_index]);
	}
	const auto element = domain.Background().MaterializeElement(cell_id);
	if (!selected_boundary_labels.empty()) {
		if (nodal_state.size() != element.connectivity.size())
			throw std::runtime_error("immersed Nitsche wall nodal-state size does not match element connectivity");
		for (const auto& node : nodal_state)
			for (double value : node)
				if (!std::isfinite(value)) throw std::runtime_error("immersed Nitsche wall nodal state is not finite");
	}
	ImmersedNitscheWallAssembly result;
	if (volume_system) {
		ValidateImmersedNitscheWallPreassembledVolumeSystem(element, nodal_state, previous_nodal_state,
			parameters, *volume_system);
		result.system = *volume_system;
	}
	else {
		if (volume_catalog.StorageMode() != CutCellVolumeQuadratureStorageMode::Expanded)
			throw std::invalid_argument("compact immersed Nitsche wall assembly requires a preassembled volume system");
		result.system = BuildNavierStokesElement(element, nodal_state, previous_nodal_state,
			parameters, volume_catalog.UsableRule(domain, cell_id));
	}
	if (selected_boundary_labels.empty()) return result;
	for (const auto value : result.system.jacobian)
		if (!std::isfinite(PetscRealPart(value)))
			throw std::overflow_error("immersed Nitsche wall volume Jacobian is not finite");
	for (const auto value : result.system.negative_residual)
		if (!std::isfinite(PetscRealPart(value)))
			throw std::overflow_error("immersed Nitsche wall volume residual is not finite");
	const auto alpha = ImmersedNitscheWallFraction(volume_catalog, cell_id, result.diagnostics);
	result.diagnostics.ghost_covered_policy = ghost_catalog != nullptr;

	const auto& surface_diagnostics = surface_catalog.Diagnostics();
	for (int label : selected_boundary_labels) {
		const auto found = surface_diagnostics.area_by_boundary_id.find(static_cast<std::uint32_t>(label));
		if (found == surface_diagnostics.area_by_boundary_id.end() || !std::isfinite(found->second)
			|| !(found->second > 0.0))
			throw std::runtime_error("selected immersed Nitsche wall label has no positive globally configured area");
		result.diagnostics.by_boundary_id.emplace(label, ImmersedNitscheWallLabelDiagnostics{});
	}
	if (surface_rule.Points().empty()) return result;
	ValidateSurfaceQuadratureRule(element, surface_rule);

	const auto selected = [&selected_boundary_labels](int label) {
		return std::binary_search(selected_boundary_labels.begin(), selected_boundary_labels.end(), label);
	};
	const auto nen = element.connectivity.size();
	const auto ndof = 4*nen;
	bool selected_any = false;
	for (std::size_t point_index = 0; point_index < surface_rule.Points().size(); ++point_index) {
		const auto& point = surface_rule.Points()[point_index];
		const auto& provenance = surface_provenance[point_index];
		auto& label_diagnostics = result.diagnostics.by_boundary_id[point.boundary_id];
		if (!selected(point.boundary_id)) {
			++label_diagnostics.skipped_points;
			AccumulateImmersedNitscheWallDiagnostic(label_diagnostics.skipped_area_m2, point.weight,
				"immersed Nitsche wall skipped-area diagnostic is not finite");
			continue;
		}
		++label_diagnostics.selected_points;
		AccumulateImmersedNitscheWallDiagnostic(label_diagnostics.selected_area_m2, point.weight,
			"immersed Nitsche wall selected-area diagnostic is not finite");
		auto basis = EvaluateBasis(element, point.parametric[0], point.parametric[1], point.parametric[2], false);
		const auto inverse_normal = std::array<double, 3>{{
			basis.inverse_jacobian[0][0]*point.normal[0] + basis.inverse_jacobian[0][1]*point.normal[1] + basis.inverse_jacobian[0][2]*point.normal[2],
			basis.inverse_jacobian[1][0]*point.normal[0] + basis.inverse_jacobian[1][1]*point.normal[1] + basis.inverse_jacobian[1][2]*point.normal[2],
			basis.inverse_jacobian[2][0]*point.normal[0] + basis.inverse_jacobian[2][1]*point.normal[1] + basis.inverse_jacobian[2][2]*point.normal[2]}};
		const double inverse_normal_norm = std::hypot(inverse_normal[0], inverse_normal[1], inverse_normal[2]);
		if (!std::isfinite(inverse_normal_norm) || !(inverse_normal_norm > 0.0))
			throw std::runtime_error("immersed Nitsche wall normal length is not finite and positive");
		const double h_n = 1.0/inverse_normal_norm;
		const double eta = ghost_catalog ? 16.0*gamma0*parameters.dynamic_viscosity/h_n
			: 16.0*gamma0*parameters.dynamic_viscosity/(alpha*h_n);
		if (!std::isfinite(h_n) || !(h_n > 0.0) || !std::isfinite(eta) || !(eta > 0.0))
			throw std::overflow_error("immersed Nitsche wall penalty is not finite and positive");
		const double eta_h_n_over_mu = eta*h_n/parameters.dynamic_viscosity;
		if (!std::isfinite(eta_h_n_over_mu) || !(eta_h_n_over_mu > 0.0))
			throw std::overflow_error("immersed Nitsche wall penalty diagnostic overflows");
		const auto g = wall_velocity(point, provenance);
		if (!QuadratureFinite(g)) throw std::runtime_error("immersed Nitsche wall velocity is not finite");
		std::array<double, 4> state{};
		double gradient[3][3]{};
		for (std::size_t a = 0; a < nen; ++a) {
			for (int i = 0; i < 4; ++i) state[i] += nodal_state[a][i]*basis.value[a];
			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j) gradient[i][j] += nodal_state[a][i]*basis.gradient[a][j];
		}
		std::array<double, 3> gap{{state[0]-g[0], state[1]-g[1], state[2]-g[2]}};
		for (double value : state) RequireFiniteImmersedNitscheWallValue(value,
			"immersed Nitsche wall interpolated state is not finite");
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j) RequireFiniteImmersedNitscheWallValue(gradient[i][j],
				"immersed Nitsche wall interpolated gradient is not finite");
		for (double value : gap) RequireFiniteImmersedNitscheWallValue(value,
			"immersed Nitsche wall gap is not finite");
		const double gap_norm = std::hypot(gap[0], gap[1], gap[2]);
		RequireFiniteImmersedNitscheWallValue(gap_norm, "immersed Nitsche wall gap diagnostic is not finite");
		result.diagnostics.maximum_gap_norm = std::max(result.diagnostics.maximum_gap_norm, gap_norm);
		result.diagnostics.minimum_h_n_m = selected_any ? std::min(result.diagnostics.minimum_h_n_m, h_n) : h_n;
		result.diagnostics.maximum_h_n_m = std::max(result.diagnostics.maximum_h_n_m, h_n);
		result.diagnostics.minimum_eta = selected_any ? std::min(result.diagnostics.minimum_eta, eta) : eta;
		result.diagnostics.maximum_eta = std::max(result.diagnostics.maximum_eta, eta);
		result.diagnostics.maximum_eta_h_n_over_mu = std::max(result.diagnostics.maximum_eta_h_n_over_mu,
			eta_h_n_over_mu);
		selected_any = true;

		std::array<double, 3> traction{};
		for (int i = 0; i < 3; ++i) {
			traction[i] = -state[3]*point.normal[i];
			for (int j = 0; j < 3; ++j)
				traction[i] += parameters.dynamic_viscosity*(gradient[i][j]+gradient[j][i])*point.normal[j];
			RequireFiniteImmersedNitscheWallValue(traction[i], "immersed Nitsche wall traction is not finite");
		}
		for (std::size_t a = 0; a < nen; ++a) {
			const auto na = basis.value[a];
			const auto& ga = basis.gradient[a];
			const auto ga_n = ImmersedNitscheWallDot(ga, point.normal, "immersed Nitsche wall test-normal derivative is not finite");
			const auto ga_gap = ImmersedNitscheWallDot(ga, gap, "immersed Nitsche wall test-gap derivative is not finite");
			for (int i = 0; i < 3; ++i) {
				const double residual = -na*traction[i]
					-parameters.dynamic_viscosity*(ga_n*gap[i]+point.normal[i]*ga_gap)+eta*na*gap[i];
				RequireFiniteImmersedNitscheWallValue(residual, "immersed Nitsche wall velocity residual is not finite");
				AccumulateImmersedNitscheWallValue(result.system.negative_residual[4*a+i], -residual*point.weight,
					"immersed Nitsche wall velocity residual contribution is not finite");
			}
			const double pressure_residual = -na*ImmersedNitscheWallDot(point.normal, gap,
				"immersed Nitsche wall pressure residual is not finite");
			RequireFiniteImmersedNitscheWallValue(pressure_residual, "immersed Nitsche wall pressure residual is not finite");
			AccumulateImmersedNitscheWallValue(result.system.negative_residual[4*a+3], -pressure_residual*point.weight,
				"immersed Nitsche wall pressure residual contribution is not finite");
			for (std::size_t b = 0; b < nen; ++b) {
				const auto nb = basis.value[b];
				const auto& gb = basis.gradient[b];
				const auto gb_n = ImmersedNitscheWallDot(gb, point.normal, "immersed Nitsche wall trial-normal derivative is not finite");
				for (int i = 0; i < 3; ++i) {
					for (int j = 0; j < 3; ++j) {
						const double tangent = -parameters.dynamic_viscosity*na*((i == j ? gb_n : 0.0)+gb[i]*point.normal[j])
							-parameters.dynamic_viscosity*nb*((i == j ? ga_n : 0.0)+ga[j]*point.normal[i])
							+(i == j ? eta*na*nb : 0.0);
						RequireFiniteImmersedNitscheWallValue(tangent, "immersed Nitsche wall velocity tangent is not finite");
						AccumulateImmersedNitscheWallValue(result.system.jacobian[(4*a+i)*ndof+4*b+j], tangent*point.weight,
							"immersed Nitsche wall velocity tangent contribution is not finite");
					}
					const double velocity_pressure = na*nb*point.normal[i]*point.weight;
					AccumulateImmersedNitscheWallValue(result.system.jacobian[(4*a+i)*ndof+4*b+3], velocity_pressure,
						"immersed Nitsche wall velocity-pressure contribution is not finite");
					AccumulateImmersedNitscheWallValue(result.system.jacobian[(4*a+3)*ndof+4*b+i], -velocity_pressure,
						"immersed Nitsche wall pressure-velocity contribution is not finite");
				}
			}
		}
	}
	return result;
}

inline ImmersedNitscheWallAssembly BuildImmersedNitscheWallElement(
	const CartesianDomainClassification& domain,
	const CutCellVolumeQuadratureCatalog& volume_catalog,
	const ImmersedSurfaceQuadratureCatalog& surface_catalog, std::uint64_t cell_id,
	const std::vector<std::array<double, 4>>& nodal_state,
	const std::vector<std::array<double, 4>>& previous_nodal_state,
	const NavierStokesParameters& parameters, const std::vector<int>& selected_boundary_labels,
	double gamma0 = 2.0)
{
	return BuildImmersedNitscheWallElementImpl(domain, volume_catalog, surface_catalog, cell_id,
		nodal_state, previous_nodal_state, parameters, selected_boundary_labels, gamma0,
		AdaptImmersedWallVelocityEvaluator([](const std::array<double, 3>&, int) { return std::array<double, 3>{{0.0, 0.0, 0.0}}; }), nullptr);
}

// Existing velocity-evaluator overload, retained bit-for-bit in its legacy
// path.  The ghost-covered overload below is explicit so callers cannot
// accidentally switch the penalty rule.
inline ImmersedNitscheWallAssembly BuildImmersedNitscheWallElement(
	const CartesianDomainClassification& domain,
	const CutCellVolumeQuadratureCatalog& volume_catalog,
	const ImmersedSurfaceQuadratureCatalog& surface_catalog, std::uint64_t cell_id,
	const std::vector<std::array<double, 4>>& nodal_state,
	const std::vector<std::array<double, 4>>& previous_nodal_state,
	const NavierStokesParameters& parameters, const std::vector<int>& selected_boundary_labels,
	double gamma0, const ImmersedWallVelocityEvaluator& wall_velocity)
{
	return BuildImmersedNitscheWallElementImpl(domain, volume_catalog, surface_catalog, cell_id,
		nodal_state, previous_nodal_state, parameters, selected_boundary_labels, gamma0, AdaptImmersedWallVelocityEvaluator(wall_velocity), nullptr);
}

inline ImmersedNitscheWallAssembly BuildImmersedNitscheWallElementMaterialAware(
	const CartesianDomainClassification& domain,
	const CutCellVolumeQuadratureCatalog& volume_catalog,
	const ImmersedSurfaceQuadratureCatalog& surface_catalog, std::uint64_t cell_id,
	const std::vector<std::array<double, 4>>& nodal_state,
	const std::vector<std::array<double, 4>>& previous_nodal_state,
	const NavierStokesParameters& parameters, const std::vector<int>& selected_boundary_labels,
	double gamma0, const ImmersedMaterialWallVelocityEvaluator& wall_velocity)
{
	return BuildImmersedNitscheWallElementImpl(domain, volume_catalog, surface_catalog, cell_id,
		nodal_state, previous_nodal_state, parameters, selected_boundary_labels, gamma0, wall_velocity, nullptr);
}

inline ImmersedNitscheWallAssembly BuildImmersedNitscheWallElement(
	const CartesianDomainClassification& domain,
	const CutCellVolumeQuadratureCatalog& volume_catalog,
	const ImmersedSurfaceQuadratureCatalog& surface_catalog, std::uint64_t cell_id,
	const std::vector<std::array<double, 4>>& nodal_state,
	const std::vector<std::array<double, 4>>& previous_nodal_state,
	const NavierStokesParameters& parameters, const std::vector<int>& selected_boundary_labels,
	const CutCellGhostPenaltyCatalog& ghost_catalog, double gamma0 = 2.0,
	const ImmersedWallVelocityEvaluator& wall_velocity = [](const std::array<double, 3>&, int) { return std::array<double, 3>{{0.0,0.0,0.0}}; })
{
	return BuildImmersedNitscheWallElementImpl(domain, volume_catalog, surface_catalog, cell_id,
		nodal_state, previous_nodal_state, parameters, selected_boundary_labels, gamma0, AdaptImmersedWallVelocityEvaluator(wall_velocity), &ghost_catalog);
}

// Streaming callers retain their pointwise volume system and hand it to the
// wall adapter.  This avoids an expanded-only catalog lookup while preserving
// the established overloads above.
inline ImmersedNitscheWallAssembly BuildImmersedNitscheWallElementFromVolumeSystem(
	const CartesianDomainClassification& domain,
	const CutCellVolumeQuadratureCatalog& volume_catalog,
	const ImmersedSurfaceQuadratureCatalog& surface_catalog, std::uint64_t cell_id,
	const std::vector<std::array<double, 4>>& nodal_state,
	const std::vector<std::array<double, 4>>& previous_nodal_state,
	const NavierStokesParameters& parameters, const std::vector<int>& selected_boundary_labels,
	const NavierStokesSystem& volume_system, double gamma0 = 2.0,
	const ImmersedWallVelocityEvaluator& wall_velocity = [](const std::array<double, 3>&, int) { return std::array<double, 3>{{0.0,0.0,0.0}}; })
{
	return BuildImmersedNitscheWallElementImpl(domain, volume_catalog, surface_catalog, cell_id,
		nodal_state, previous_nodal_state, parameters, selected_boundary_labels, gamma0, AdaptImmersedWallVelocityEvaluator(wall_velocity),
		nullptr, &volume_system);
}

inline ImmersedNitscheWallAssembly BuildImmersedNitscheWallElementFromVolumeSystem(
	const CartesianDomainClassification& domain,
	const CutCellVolumeQuadratureCatalog& volume_catalog,
	const ImmersedSurfaceQuadratureCatalog& surface_catalog, std::uint64_t cell_id,
	const std::vector<std::array<double, 4>>& nodal_state,
	const std::vector<std::array<double, 4>>& previous_nodal_state,
	const NavierStokesParameters& parameters, const std::vector<int>& selected_boundary_labels,
	const NavierStokesSystem& volume_system, const CutCellGhostPenaltyCatalog& ghost_catalog,
	double gamma0 = 2.0,
	const ImmersedWallVelocityEvaluator& wall_velocity = [](const std::array<double, 3>&, int) { return std::array<double, 3>{{0.0,0.0,0.0}}; })
{
	return BuildImmersedNitscheWallElementImpl(domain, volume_catalog, surface_catalog, cell_id,
		nodal_state, previous_nodal_state, parameters, selected_boundary_labels, gamma0, AdaptImmersedWallVelocityEvaluator(wall_velocity),
		&ghost_catalog, &volume_system);
}

} // namespace iga

#endif
