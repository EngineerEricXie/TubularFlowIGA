#ifndef IGA_IMMERSED_FLOW_PORT_HPP
#define IGA_IMMERSED_FLOW_PORT_HPP

// Static open-patch flow ports for the immersed Cartesian runtime.  The
// catalog owns surface geometry and supplies physical dA/outward normals; this
// adapter only evaluates field traces and emits bounded element-local terms.
#include "ImmersedSurfaceQuadrature.hpp"
#include "NavierStokesElement.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

enum class ImmersedFlowPortControlMode {
	Pressure,
	MeanNormalTraction,
	FlowRate,
	TotalPressure // Deliberately rejected: Phase 6 has no averaging contract.
};

struct ImmersedFlowPortDefinition {
	std::string id;
	int boundary_label = -1;
	ImmersedFlowPortControlMode control_mode = ImmersedFlowPortControlMode::Pressure;
	double value = 0.0;
};

struct ImmersedFlowPortMeasurement {
	double area_m2 = 0.0;
	double outward_flow_m3_s = 0.0;
	double mean_pressure_pa = 0.0;
	double mean_normal_traction_pa = 0.0;
	// This is a trace diagnostic (m^2/s^2), not a total-pressure control.
	// Keeping it finite is useful when callers form a total-pressure diagnostic.
	double mean_velocity_squared_m2_s2 = 0.0;
};

struct ImmersedFlowPortElementAssembly {
	// c = integral N^T n dA, stored in the velocity positions of a 4*nen
	// local vector.  It is used for both flow-controller off-diagonal blocks.
	std::vector<PetscScalar> flow_coefficient;
	std::vector<PetscScalar> negative_residual;
};

inline bool IsImmersedFlowPressureLike(ImmersedFlowPortControlMode mode)
{
	return mode == ImmersedFlowPortControlMode::Pressure
		|| mode == ImmersedFlowPortControlMode::MeanNormalTraction;
}

inline void ValidateImmersedFlowPortDefinition(const ImmersedFlowPortDefinition& port)
{
	if (port.id.empty()) throw std::invalid_argument("immersed flow port id must be nonempty");
	if (port.boundary_label <= 0) throw std::invalid_argument("immersed flow port boundary label must be positive");
	if (port.control_mode == ImmersedFlowPortControlMode::TotalPressure)
		throw std::invalid_argument("total-pressure immersed port control is unsupported in Phase 6");
	if (!std::isfinite(port.value)) throw std::invalid_argument("immersed flow port value must be finite");
}

inline double CheckedImmersedFlowPortProduct(double left, double right, const char* what)
{
	const double result = left*right;
	if (!std::isfinite(result)) throw std::overflow_error(std::string("immersed flow port ")+what+" overflows or is not finite");
	return result;
}

inline void AddImmersedFlowPortFinite(double& total, double value, const char* what)
{
	if (!std::isfinite(value) || !std::isfinite(total)) throw std::overflow_error(std::string("immersed flow port ")+what+" is not finite");
	total += value;
	if (!std::isfinite(total)) throw std::overflow_error(std::string("immersed flow port ")+what+" overflows");
}

// Physical all-label trace that completes the conservative resolved mixed
// volume form.  It intentionally has no port-control/load policy: every
// retained surface point contributes exactly once, independently of whether
// it is later used by Nitsche or a port controller.
//
// R_u,p = -int_Gamma N_a p n,  R_p,u = +int_Gamma N_a (u.n).
// The stored vector is -R, hence the +N_a p n and -N_a(u.n) additions below.
inline NavierStokesSystem BuildImmersedConservativeMixedTraceElement(
	const Element& element, const SurfaceQuadratureRule& rule,
	const std::vector<std::array<double, 4>>& nodal_state)
{
	if (nodal_state.size() != element.connectivity.size())
		throw std::invalid_argument("immersed conservative mixed trace nodal-state size is invalid");
	ValidateSurfaceQuadratureRule(element, rule);
	const std::size_t ndof = 4*element.connectivity.size();
	NavierStokesSystem result{std::vector<PetscScalar>(ndof*ndof, 0.0),
		std::vector<PetscScalar>(ndof, 0.0)};
	for (const auto& point : rule.Points()) {
		const auto basis = EvaluateBasis(element, point.parametric[0], point.parametric[1], point.parametric[2], false);
		double pressure = 0.0;
		std::array<double, 3> velocity{};
		for (std::size_t b = 0; b < element.connectivity.size(); ++b) {
			AddImmersedFlowPortFinite(pressure, CheckedImmersedFlowPortProduct(nodal_state[b][3], basis.value[b], "mixed trace pressure"), "mixed trace pressure accumulation");
			for (int component = 0; component < 3; ++component)
				AddImmersedFlowPortFinite(velocity[component], CheckedImmersedFlowPortProduct(nodal_state[b][component], basis.value[b], "mixed trace velocity"), "mixed trace velocity accumulation");
		}
		double normal_velocity = 0.0;
		for (int component = 0; component < 3; ++component)
			AddImmersedFlowPortFinite(normal_velocity, CheckedImmersedFlowPortProduct(velocity[component], point.normal[component], "mixed trace normal velocity"), "mixed trace normal velocity accumulation");
		for (std::size_t a = 0; a < element.connectivity.size(); ++a) {
			const double test_weight = CheckedImmersedFlowPortProduct(basis.value[a], point.weight, "mixed trace test weight");
			const std::size_t pressure_row = 4*a+3;
			double continuity_rhs = PetscRealPart(result.negative_residual[pressure_row]);
			AddImmersedFlowPortFinite(continuity_rhs, -CheckedImmersedFlowPortProduct(test_weight, normal_velocity, "mixed trace continuity residual"), "mixed trace continuity residual accumulation");
			result.negative_residual[pressure_row] = continuity_rhs;
			for (int component = 0; component < 3; ++component) {
				const std::size_t velocity_row = 4*a+component;
				double momentum_rhs = PetscRealPart(result.negative_residual[velocity_row]);
				AddImmersedFlowPortFinite(momentum_rhs, CheckedImmersedFlowPortProduct(CheckedImmersedFlowPortProduct(test_weight, pressure, "mixed trace momentum pressure"), point.normal[component], "mixed trace momentum normal"), "mixed trace momentum residual accumulation");
				result.negative_residual[velocity_row] = momentum_rhs;
			}
			for (std::size_t b = 0; b < element.connectivity.size(); ++b)
				for (int component = 0; component < 3; ++component) {
					const double coefficient = CheckedImmersedFlowPortProduct(CheckedImmersedFlowPortProduct(test_weight, basis.value[b], "mixed trace Jacobian basis"), point.normal[component], "mixed trace Jacobian normal");
					const std::size_t up = (4*a+component)*ndof+4*b+3;
					const std::size_t pu = pressure_row*ndof+4*b+component;
					double up_total = PetscRealPart(result.jacobian[up]);
					AddImmersedFlowPortFinite(up_total, -coefficient, "mixed trace velocity-pressure Jacobian accumulation");
					result.jacobian[up] = up_total;
					double pu_total = PetscRealPart(result.jacobian[pu]);
					AddImmersedFlowPortFinite(pu_total, coefficient, "mixed trace pressure-velocity Jacobian accumulation");
					result.jacobian[pu] = pu_total;
				}
		}
	}
	return result;
}

inline ImmersedFlowPortElementAssembly BuildImmersedFlowPortElement(
	const Element& element, const SurfaceQuadratureRule& rule, int boundary_label,
	ImmersedFlowPortControlMode mode, double value,
	const std::vector<std::array<double, 4>>& nodal_state)
{
	if (boundary_label <= 0 || !std::isfinite(value) || nodal_state.size() != element.connectivity.size())
		throw std::invalid_argument("immersed flow port element inputs are invalid");
	if (mode == ImmersedFlowPortControlMode::TotalPressure)
		throw std::invalid_argument("total-pressure immersed port control is unsupported in Phase 6");
	ValidateSurfaceQuadratureRule(element, rule);
	const std::size_t ndof = 4*element.connectivity.size();
	ImmersedFlowPortElementAssembly result{std::vector<PetscScalar>(ndof, 0.0),
		std::vector<PetscScalar>(ndof, 0.0)};
	const double pressure = mode == ImmersedFlowPortControlMode::MeanNormalTraction ? -value : value;
	for (const auto& point : rule.Points()) {
		if (point.boundary_id != boundary_label) continue;
		const auto basis = EvaluateBasis(element, point.parametric[0], point.parametric[1], point.parametric[2], false);
		for (std::size_t a = 0; a < element.connectivity.size(); ++a) {
			for (int component = 0; component < 3; ++component) {
				const double c = CheckedImmersedFlowPortProduct(CheckedImmersedFlowPortProduct(basis.value[a], point.normal[component], "coefficient"), point.weight, "coefficient");
				const std::size_t row = 4*a+component;
				double coefficient_total = PetscRealPart(result.flow_coefficient[row]);
				AddImmersedFlowPortFinite(coefficient_total, c, "coefficient accumulation");
				result.flow_coefficient[row] = coefficient_total;
				// Body-fitted convention: sigma*n=-p_b*n contributes -p_b*N*n
				// to the stored negative residual.
				if (IsImmersedFlowPressureLike(mode)) {
					double residual_total = PetscRealPart(result.negative_residual[row]);
					AddImmersedFlowPortFinite(residual_total, -CheckedImmersedFlowPortProduct(pressure, c, "pressure load"), "pressure load accumulation");
					result.negative_residual[row] = residual_total;
				}
			}
		}
	}
	return result;
}

inline ImmersedFlowPortMeasurement MeasureImmersedFlowPortElement(
	const Element& element, const SurfaceQuadratureRule& rule, int boundary_label,
	const std::vector<std::array<double, 4>>& nodal_state, double viscosity)
{
	if (boundary_label <= 0 || nodal_state.size() != element.connectivity.size()
		|| !std::isfinite(viscosity) || !(viscosity > 0.0))
		throw std::invalid_argument("immersed flow port measurement inputs are invalid");
	ValidateSurfaceQuadratureRule(element, rule);
	ImmersedFlowPortMeasurement result;
	for (const auto& point : rule.Points()) {
		if (point.boundary_id != boundary_label) continue;
		const auto basis = EvaluateBasis(element, point.parametric[0], point.parametric[1], point.parametric[2], false);
		std::array<double, 4> state{}; double gradient[3][3]{};
		for (std::size_t a = 0; a < element.connectivity.size(); ++a)
			for (int field = 0; field < 4; ++field) {
				AddImmersedFlowPortFinite(state[field], CheckedImmersedFlowPortProduct(nodal_state[a][field], basis.value[a], "state trace"), "state trace accumulation");
				if (field < 3) for (int direction = 0; direction < 3; ++direction)
					AddImmersedFlowPortFinite(gradient[field][direction], CheckedImmersedFlowPortProduct(nodal_state[a][field], basis.gradient[a][direction], "gradient trace"), "gradient trace accumulation");
			}
		double normal_velocity = 0.0, normal_traction = -state[3];
		double velocity_squared = 0.0;
		for (int i = 0; i < 3; ++i) {
			AddImmersedFlowPortFinite(normal_velocity, CheckedImmersedFlowPortProduct(state[i], point.normal[i], "normal velocity"), "normal velocity accumulation");
			AddImmersedFlowPortFinite(velocity_squared, CheckedImmersedFlowPortProduct(state[i], state[i], "kinetic trace"), "kinetic trace accumulation");
			for (int j = 0; j < 3; ++j) {
				const double symmetric_gradient = gradient[i][j]+gradient[j][i];
				if (!std::isfinite(symmetric_gradient)) throw std::overflow_error("immersed flow port traction gradient is not finite");
				const double traction_term = CheckedImmersedFlowPortProduct(
					CheckedImmersedFlowPortProduct(CheckedImmersedFlowPortProduct(viscosity, point.normal[i], "traction"), symmetric_gradient, "traction"),
					point.normal[j], "traction");
				AddImmersedFlowPortFinite(normal_traction, traction_term, "traction accumulation");
			}
		}
		if (!std::isfinite(normal_velocity) || !std::isfinite(normal_traction) || !std::isfinite(velocity_squared))
			throw std::overflow_error("immersed flow port measurement is not finite");
		AddImmersedFlowPortFinite(result.area_m2, point.weight, "measurement area");
		AddImmersedFlowPortFinite(result.outward_flow_m3_s, CheckedImmersedFlowPortProduct(normal_velocity, point.weight, "measurement flow"), "measurement flow accumulation");
		AddImmersedFlowPortFinite(result.mean_pressure_pa, CheckedImmersedFlowPortProduct(state[3], point.weight, "measurement pressure"), "measurement pressure accumulation");
		AddImmersedFlowPortFinite(result.mean_normal_traction_pa, CheckedImmersedFlowPortProduct(normal_traction, point.weight, "measurement traction"), "measurement traction accumulation");
		AddImmersedFlowPortFinite(result.mean_velocity_squared_m2_s2, CheckedImmersedFlowPortProduct(velocity_squared, point.weight, "measurement kinetic"), "measurement kinetic accumulation");
	}
	if (!std::isfinite(result.area_m2) || !(result.area_m2 > 0.0))
		throw std::runtime_error("immersed flow port has zero retained area");
	result.mean_pressure_pa /= result.area_m2;
	result.mean_normal_traction_pa /= result.area_m2;
	result.mean_velocity_squared_m2_s2 /= result.area_m2;
	if (!std::isfinite(result.mean_pressure_pa) || !std::isfinite(result.mean_normal_traction_pa)
		|| !std::isfinite(result.mean_velocity_squared_m2_s2)) throw std::overflow_error("immersed flow port measurement average is not finite");
	return result;
}

} // namespace iga

#endif
