#ifndef NAVIER_STOKES_ELEMENT_HPP
#define NAVIER_STOKES_ELEMENT_HPP

#include "IgaDatabase.hpp"
#include "Quadrature.hpp"
#include "TransportElement.hpp"

#include <array>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <vector>

namespace iga {

struct NavierStokesSystem {
	std::vector<PetscScalar> jacobian;
	std::vector<PetscScalar> negative_residual;
};

struct NavierStokesParameters {
	double density = 1.0;
	double dynamic_viscosity = 0.0;
	double dt = 0.0;
};

// The established body-fitted weak form integrates the resolved mixed pair by
// parts in the volume.  Immersed cut cells can instead retain the resolved
// gradients in the volume and supply the cancelling physical trace explicitly.
// VMS/PSPG remains identical in both forms.
enum class NavierStokesResolvedMixedForm {
	LegacyBodyFitted,
	Conservative
};

// Physical body-force density per unit volume (N/m^3).  Keeping it as a
// point evaluator makes manufactured loads and spatially varying gravity
// unambiguous, while the zero evaluator preserves the established API.
using NavierStokesBodyForceEvaluator = std::function<std::array<double, 3>(
	const std::array<double, 3>& physical)>;

inline void Stabilization(const std::array<std::array<double, 3>, 3>& inverse_jacobian,
	const std::array<double, 4>& state, double kinematic_viscosity, double dt,
	double& tau_m, double& tau_c)
{
	double metric[3][3]{};
	double direction[3]{};
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k) metric[i][j] += inverse_jacobian[k][i] * inverse_jacobian[k][j];
			direction[i] += inverse_jacobian[j][i];
		}
	double metric_norm = 0.0, direction_norm = 0.0, velocity_metric = 0.0;
	for (int i = 0; i < 3; ++i) {
		direction_norm += direction[i] * direction[i];
			for (int j = 0; j < 3; ++j) {
				metric_norm += metric[i][j] * metric[i][j];
				velocity_metric += state[i] * metric[i][j] * state[j];
			}
	}
	const auto temporal_scale = dt > 0.0 ? 4.0/(dt*dt) : 0.0;
	tau_m = 1.0 / std::sqrt(temporal_scale + velocity_metric
		+ (1.0/12.0) * kinematic_viscosity * kinematic_viscosity * metric_norm);
	tau_c = 1.0 / (tau_m * direction_norm);
}

// The point visitor lets cut-cell callers stream compact quadrature without
// materializing its logical points.  The callback receives one point at a
// time and must invoke the supplied consumer in canonical rule order.
template <class PointVisitor> inline NavierStokesSystem BuildNavierStokesElementFromPoints(const Element& element,
	const std::vector<std::array<double, 4>>& nodal_state,
	const std::vector<std::array<double, 4>>& previous_nodal_state,
	const NavierStokesParameters& parameters, PointVisitor&& visit_points,
	const NavierStokesBodyForceEvaluator& body_force,
	NavierStokesResolvedMixedForm resolved_mixed_form = NavierStokesResolvedMixedForm::LegacyBodyFitted)
{
	if (!(parameters.density > 0.0) || !(parameters.dynamic_viscosity > 0.0) || parameters.dt < 0.0)
		throw std::runtime_error("invalid Navier-Stokes density, dynamic viscosity, or time step");
	if (nodal_state.size() != element.connectivity.size())
		throw std::runtime_error("Navier-Stokes nodal-state size does not match element connectivity");
	if (parameters.dt > 0.0 && previous_nodal_state.size() != element.connectivity.size())
		throw std::runtime_error("transient Navier-Stokes requires a matching previous nodal state");
	if (!body_force) throw std::invalid_argument("Navier-Stokes body-force evaluator is required");
	const auto density = parameters.density;
	const auto viscosity = parameters.dynamic_viscosity;
	const auto kinematic_viscosity = viscosity/density;
	const auto nen = element.connectivity.size();
	const auto ndof = 4 * nen;
	NavierStokesSystem system{std::vector<PetscScalar>(ndof*ndof, 0.0), std::vector<PetscScalar>(ndof, 0.0)};
	visit_points([&](const VolumeQuadraturePoint& point) {
				auto basis = EvaluateBasis(element, point.parametric[0], point.parametric[1],
					point.parametric[2], true);
				const auto measure = point.weight*basis.raw_determinant;
				std::array<double, 4> state{};
				const auto force = body_force(EvaluateElementGeometry(element, point.parametric).physical);
				for (const double value : force)
					if (!std::isfinite(value)) throw std::runtime_error("Navier-Stokes body force is not finite");
				std::array<double, 3> previous_velocity{};
				double gradient[4][3]{};
				double hessian[4][3][3]{};
				for (std::size_t a = 0; a < nen; ++a)
					for (int field = 0; field < 4; ++field) {
						state[field] += nodal_state[a][field] * basis.value[a];
						for (int i = 0; i < 3; ++i) {
							gradient[field][i] += nodal_state[a][field] * basis.gradient[a][i];
							for (int j = 0; j < 3; ++j) hessian[field][i][j] += nodal_state[a][field] * basis.hessian[a][i][j];
						}
					}
				if (parameters.dt > 0.0)
					for (std::size_t a = 0; a < nen; ++a)
						for (int component = 0; component < 3; ++component)
							previous_velocity[component] += previous_nodal_state[a][component] * basis.value[a];
				std::array<double, 3> time_derivative{};
				if (parameters.dt > 0.0)
					for (int component = 0; component < 3; ++component)
						time_derivative[component] = (state[component]-previous_velocity[component])/parameters.dt;
				double tau_m = 0.0, tau_c = 0.0;
				Stabilization(basis.inverse_jacobian, state, kinematic_viscosity,
					parameters.dt, tau_m, tau_c);
				std::array<double, 3> fine_velocity{};
				for (int component = 0; component < 3; ++component) {
					const auto convection = state[0]*gradient[component][0] + state[1]*gradient[component][1] + state[2]*gradient[component][2];
					const auto pressure_gradient = gradient[3][component];
					const auto laplacian = hessian[component][0][0] + hessian[component][1][1] + hessian[component][2][2];
					fine_velocity[component] = -tau_m * (time_derivative[component] + convection
						+ pressure_gradient/density - kinematic_viscosity*laplacian-force[component]/density);
				}
				const auto fine_pressure = -density*tau_c
					* (gradient[0][0] + gradient[1][1] + gradient[2][2]);
				double metric[3][3]{};
				for (int i = 0; i < 3; ++i)
					for (int j = 0; j < 3; ++j)
						for (int k = 0; k < 3; ++k)
							metric[i][j] += basis.inverse_jacobian[k][i]*basis.inverse_jacobian[k][j];

				for (std::size_t a = 0; a < nen; ++a) {
					const auto na = basis.value[a];
					const auto& ga = basis.gradient[a];
					std::array<double, 4> residual{};
					for (int component = 0; component < 3; ++component) {
						residual[component] = density*na*time_derivative[component]-na*force[component]
							- ga[component]*fine_pressure;
						if (resolved_mixed_form == NavierStokesResolvedMixedForm::Conservative)
							residual[component] += na*gradient[3][component];
						else residual[component] -= ga[component]*state[3];
						for (int direction = 0; direction < 3; ++direction) {
							residual[component] += viscosity * ga[direction] * (gradient[component][direction] + gradient[direction][component]);
							residual[component] += density*na * (state[direction]+fine_velocity[direction]) * gradient[component][direction];
							residual[component] -= density*ga[direction] * fine_velocity[component] * (state[direction]+fine_velocity[direction]);
						}
					}
					residual[3] = -ga[0]*fine_velocity[0] - ga[1]*fine_velocity[1] - ga[2]*fine_velocity[2];
					if (resolved_mixed_form == NavierStokesResolvedMixedForm::Conservative)
						for (int direction = 0; direction < 3; ++direction)
							residual[3] -= ga[direction]*state[direction];
					else residual[3] = na*(gradient[0][0]+gradient[1][1]+gradient[2][2])
						- ga[0]*fine_velocity[0] - ga[1]*fine_velocity[1] - ga[2]*fine_velocity[2];
					for (int field = 0; field < 4; ++field)
						system.negative_residual[4*a+field] -= residual[field]*measure;

					// Differentiate the residual above exactly.  In particular, the
					// continuity row contains -grad(N_a).fine_velocity: its velocity
					// derivative includes the strong viscous Laplacian and both
					// stabilization parameters depend on the resolved velocity.
					for (std::size_t b = 0; b < nen; ++b) {
						const auto nb = basis.value[b];
						const auto& gb = basis.gradient[b];
						const double laplacian_b = basis.hessian[b][0][0]+basis.hessian[b][1][1]+basis.hessian[b][2][2];
						for (int field = 0; field < 4; ++field) {
							std::array<double, 3> delta_velocity{};
							std::array<double, 3> delta_laplacian{};
							std::array<double, 3> delta_pressure_gradient{};
							const double delta_pressure = field == 3 ? nb : 0.0;
							if (field < 3) {
								delta_velocity[field] = nb;
								delta_laplacian[field] = laplacian_b;
							} else delta_pressure_gradient = gb;
							double velocity_metric_direction = 0.0;
							for (int i = 0; i < 3; ++i)
								for (int j = 0; j < 3; ++j)
									velocity_metric_direction += delta_velocity[i]*metric[i][j]*state[j];
							const double delta_tau_m = -tau_m*tau_m*tau_m*velocity_metric_direction;
							const double delta_tau_c = tau_c*tau_m*tau_m*velocity_metric_direction;
							const double delta_divergence = field < 3 ? gb[field] : 0.0;
							std::array<double, 3> strong_residual_derivative{};
							std::array<double, 3> fine_velocity_derivative{};
							for (int component = 0; component < 3; ++component) {
								for (int direction = 0; direction < 3; ++direction)
									strong_residual_derivative[component] += delta_velocity[direction]*gradient[component][direction]
										+ state[direction]*(field == component ? gb[direction] : 0.0);
								if (parameters.dt > 0.0 && field == component)
									strong_residual_derivative[component] += nb/parameters.dt;
								strong_residual_derivative[component] += delta_pressure_gradient[component]/density
									- kinematic_viscosity*delta_laplacian[component];
								const double strong_residual = time_derivative[component]+state[0]*gradient[component][0]
									+ state[1]*gradient[component][1]+state[2]*gradient[component][2]
									+ gradient[3][component]/density-kinematic_viscosity
									*(hessian[component][0][0]+hessian[component][1][1]+hessian[component][2][2])
									- force[component]/density;
								fine_velocity_derivative[component] = -delta_tau_m*strong_residual
									- tau_m*strong_residual_derivative[component];
							}
							const double delta_fine_pressure = -density*(delta_tau_c
								*(gradient[0][0]+gradient[1][1]+gradient[2][2])+tau_c*delta_divergence);
							for (int i = 0; i < 3; ++i) {
								double tangent = density*na*(parameters.dt > 0.0 && field == i ? nb/parameters.dt : 0.0)
									- ga[i]*delta_fine_pressure;
								if (resolved_mixed_form == NavierStokesResolvedMixedForm::Conservative)
									tangent += na*delta_pressure_gradient[i];
								else tangent -= ga[i]*delta_pressure;
								for (int direction = 0; direction < 3; ++direction) {
									tangent += viscosity*ga[direction]*((field == i ? gb[direction] : 0.0)
										+ (field == direction ? gb[i] : 0.0));
									tangent += density*na*((delta_velocity[direction]+fine_velocity_derivative[direction])
										*gradient[i][direction]+(state[direction]+fine_velocity[direction])
										*(field == i ? gb[direction] : 0.0));
									tangent -= density*ga[direction]*(fine_velocity_derivative[i]
										*(state[direction]+fine_velocity[direction])+fine_velocity[i]
										*(delta_velocity[direction]+fine_velocity_derivative[direction]));
								}
								system.jacobian[(4*a+i)*ndof+4*b+field] += tangent*measure;
							}
							double continuity_tangent = resolved_mixed_form == NavierStokesResolvedMixedForm::Conservative
								? -(ga[0]*delta_velocity[0]+ga[1]*delta_velocity[1]+ga[2]*delta_velocity[2])
								: na*delta_divergence;
							for (int direction = 0; direction < 3; ++direction)
								continuity_tangent -= ga[direction]*fine_velocity_derivative[direction];
							system.jacobian[(4*a+3)*ndof+4*b+field] += continuity_tangent*measure;
						}
					}
				}
			});
	return system;
}

inline NavierStokesSystem BuildNavierStokesElement(const Element& element,
	const std::vector<std::array<double, 4>>& nodal_state,
	const std::vector<std::array<double, 4>>& previous_nodal_state,
	const NavierStokesParameters& parameters, const VolumeQuadratureRule& quadrature,
	const NavierStokesBodyForceEvaluator& body_force)
{
	ValidateVolumeQuadratureRule(element, quadrature);
	return BuildNavierStokesElementFromPoints(element, nodal_state, previous_nodal_state, parameters,
		[&quadrature](const auto& consume) { for (const auto& point : quadrature.Points()) consume(point); }, body_force);
}

inline NavierStokesSystem BuildNavierStokesElement(const Element& element,
	const std::vector<std::array<double, 4>>& nodal_state,
	const std::vector<std::array<double, 4>>& previous_nodal_state,
	const NavierStokesParameters& parameters, const VolumeQuadratureRule& quadrature)
{
	return BuildNavierStokesElement(element, nodal_state, previous_nodal_state, parameters, quadrature,
		[](const std::array<double, 3>&) { return std::array<double, 3>{{0.0,0.0,0.0}}; });
}

inline NavierStokesSystem BuildNavierStokesElement(const Element& element,
	const std::vector<std::array<double, 4>>& nodal_state,
	const std::vector<std::array<double, 4>>& previous_nodal_state,
	const NavierStokesParameters& parameters)
{
	FullCell4x4x4VolumeQuadratureProvider quadrature(element);
	return BuildNavierStokesElement(element, nodal_state, previous_nodal_state,
		parameters, quadrature.Rule());
}

inline NavierStokesSystem BuildNavierStokesElement(const Element& element,
	const std::vector<std::array<double, 4>>& nodal_state, double viscosity)
{
	return BuildNavierStokesElement(element, nodal_state, {}, {1.0, viscosity, 0.0});
}

} // namespace iga

#endif
