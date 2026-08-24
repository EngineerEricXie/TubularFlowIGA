#ifndef NAVIER_STOKES_ELEMENT_HPP
#define NAVIER_STOKES_ELEMENT_HPP

#include "IgaDatabase.hpp"
#include "TransportElement.hpp"

#include <array>
#include <cmath>
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

inline NavierStokesSystem BuildNavierStokesElement(const Element& element,
	const std::vector<std::array<double, 4>>& nodal_state,
	const std::vector<std::array<double, 4>>& previous_nodal_state,
	const NavierStokesParameters& parameters)
{
	if (!(parameters.density > 0.0) || !(parameters.dynamic_viscosity > 0.0) || parameters.dt < 0.0)
		throw std::runtime_error("invalid Navier-Stokes density, dynamic viscosity, or time step");
	if (nodal_state.size() != element.connectivity.size())
		throw std::runtime_error("Navier-Stokes nodal-state size does not match element connectivity");
	if (parameters.dt > 0.0 && previous_nodal_state.size() != element.connectivity.size())
		throw std::runtime_error("transient Navier-Stokes requires a matching previous nodal state");
	const auto density = parameters.density;
	const auto viscosity = parameters.dynamic_viscosity;
	const auto kinematic_viscosity = viscosity/density;
	constexpr std::array<double, 4> points{{0.06943184420297371, 0.33000947820757187, 0.6699905217924281, 0.9305681557970262}};
	constexpr std::array<double, 4> weights{{0.3478548451374539, 0.6521451548625461, 0.6521451548625461, 0.3478548451374539}};
	const auto nen = element.connectivity.size();
	const auto ndof = 4 * nen;
	NavierStokesSystem system{std::vector<PetscScalar>(ndof*ndof, 0.0), std::vector<PetscScalar>(ndof, 0.0)};
	for (std::size_t qz = 0; qz < 4; ++qz)
		for (std::size_t qy = 0; qy < 4; ++qy)
			for (std::size_t qx = 0; qx < 4; ++qx) {
				auto basis = EvaluateBasis(element, points[qx], points[qy], points[qz], true);
				const auto measure = weights[qx] * weights[qy] * weights[qz] * basis.determinant;
				std::array<double, 4> state{};
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
						+ pressure_gradient/density - kinematic_viscosity*laplacian);
				}
				const auto fine_pressure = -density*tau_c
					* (gradient[0][0] + gradient[1][1] + gradient[2][2]);

				for (std::size_t a = 0; a < nen; ++a) {
					const auto na = basis.value[a];
					const auto& ga = basis.gradient[a];
					std::array<double, 4> residual{};
					for (int component = 0; component < 3; ++component) {
						residual[component] = density*na*time_derivative[component]
							- ga[component]*state[3] - ga[component]*fine_pressure;
						for (int direction = 0; direction < 3; ++direction) {
							residual[component] += viscosity * ga[direction] * (gradient[component][direction] + gradient[direction][component]);
							residual[component] += density*na * (state[direction]+fine_velocity[direction]) * gradient[component][direction];
							residual[component] -= density*ga[direction] * fine_velocity[component] * (state[direction]+fine_velocity[direction]);
						}
					}
					residual[3] = na*(gradient[0][0]+gradient[1][1]+gradient[2][2])
						- ga[0]*fine_velocity[0] - ga[1]*fine_velocity[1] - ga[2]*fine_velocity[2];
					for (int field = 0; field < 4; ++field)
						system.negative_residual[4*a+field] -= residual[field]*measure;

					for (std::size_t b = 0; b < nen; ++b) {
						const auto nb = basis.value[b];
						const auto& gb = basis.gradient[b];
						const auto convection_b = state[0]*gb[0] + state[1]*gb[1] + state[2]*gb[2];
						const auto streamline_a = state[0]*ga[0] + state[1]*ga[1] + state[2]*ga[2];
						double tangent[4][4]{};
						const auto mass_b = parameters.dt > 0.0 ? nb/parameters.dt : 0.0;
						const auto diagonal = density*na*(mass_b+convection_b)
							+ viscosity*(ga[0]*gb[0]+ga[1]*gb[1]+ga[2]*gb[2])
							+ density*tau_m*streamline_a*(mass_b+convection_b);
						for (int i = 0; i < 3; ++i)
							for (int j = 0; j < 3; ++j)
								tangent[i][j] = viscosity*ga[j]*gb[i] + density*tau_c*ga[i]*gb[j];
						for (int i = 0; i < 3; ++i) tangent[i][i] += diagonal;
						for (int i = 0; i < 3; ++i) {
							tangent[i][3] = -ga[i]*nb + tau_m*streamline_a*gb[i];
							tangent[3][i] = na*gb[i] + tau_m*ga[i]*(mass_b+convection_b);
						}
						tangent[3][3] = (tau_m/density)*(ga[0]*gb[0]+ga[1]*gb[1]+ga[2]*gb[2]);
						for (int i = 0; i < 4; ++i)
							for (int j = 0; j < 4; ++j)
								system.jacobian[(4*a+i)*ndof + 4*b+j] += tangent[i][j]*measure;
					}
				}
			}
	return system;
}

inline NavierStokesSystem BuildNavierStokesElement(const Element& element,
	const std::vector<std::array<double, 4>>& nodal_state, double viscosity)
{
	return BuildNavierStokesElement(element, nodal_state, {}, {1.0, viscosity, 0.0});
}

} // namespace iga

#endif
