#include "BoundaryFlow.hpp"
#include "ImmersedFlowPort.hpp"
#include "NavierStokesElement.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

iga::Element UnitElement()
{
	iga::Element element;
	element.connectivity.resize(64);
	element.extraction.resize(64);
	std::size_t p = 0;
	for (int k = 0; k < 4; ++k)
		for (int j = 0; j < 4; ++j)
			for (int i = 0; i < 4; ++i, ++p) {
				element.connectivity[p] = static_cast<std::int32_t>(p);
				element.extraction[p].fill(0.0);
				element.extraction[p][p] = 1.0;
				element.bezier_points[p] = {i/3.0, j/3.0, k/3.0};
			}
	return element;
}

iga::SurfaceQuadratureRule UnitSurfaceQuadrature()
{
	std::vector<iga::SurfaceQuadraturePoint> points;
	points.reserve(96);
	for (int face = 0; face < 6; ++face)
		for (std::size_t q1 = 0; q1 < 4; ++q1)
			for (std::size_t q2 = 0; q2 < 4; ++q2) {
				std::array<double, 3> x{};
				std::array<double, 3> n{};
				switch (face) {
				case 0: x = {{iga::kGaussFourPoints[q1], iga::kGaussFourPoints[q2], 0.0}}; n = {{0.0,0.0,-1.0}}; break;
				case 1: x = {{iga::kGaussFourPoints[q1], 0.0, iga::kGaussFourPoints[q2]}}; n = {{0.0,-1.0,0.0}}; break;
				case 2: x = {{1.0, iga::kGaussFourPoints[q1], iga::kGaussFourPoints[q2]}}; n = {{1.0,0.0,0.0}}; break;
				case 3: x = {{iga::kGaussFourPoints[q1], 1.0, iga::kGaussFourPoints[q2]}}; n = {{0.0,1.0,0.0}}; break;
				case 4: x = {{0.0, iga::kGaussFourPoints[q1], iga::kGaussFourPoints[q2]}}; n = {{-1.0,0.0,0.0}}; break;
				default: x = {{iga::kGaussFourPoints[q1], iga::kGaussFourPoints[q2], 1.0}}; n = {{0.0,0.0,1.0}}; break;
				}
				points.push_back({x, x, n, iga::kGaussFourWeights[q1]*iga::kGaussFourWeights[q2]/4.0, face});
			}
	return iga::SurfaceQuadratureRule(std::move(points));
}

void AddSystem(iga::NavierStokesSystem& left, const iga::NavierStokesSystem& right)
{
	assert(left.jacobian.size() == right.jacobian.size());
	assert(left.negative_residual.size() == right.negative_residual.size());
	for (std::size_t i = 0; i < left.jacobian.size(); ++i) left.jacobian[i] += right.jacobian[i];
	for (std::size_t i = 0; i < left.negative_residual.size(); ++i) left.negative_residual[i] += right.negative_residual[i];
}

void LegacyStabilization(const std::array<std::array<double, 3>, 3>& inverse_jacobian,
	const std::array<double, 4>& state, double kinematic_viscosity, double dt,
	double& tau_m, double& tau_c)
{
	double metric[3][3]{};
	double direction[3]{};
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k)
				metric[i][j] += inverse_jacobian[k][i]*inverse_jacobian[k][j];
			direction[i] += inverse_jacobian[j][i];
		}
	double metric_norm = 0.0;
	double direction_norm = 0.0;
	double velocity_metric = 0.0;
	for (int i = 0; i < 3; ++i) {
		direction_norm += direction[i]*direction[i];
		for (int j = 0; j < 3; ++j) {
			metric_norm += metric[i][j]*metric[i][j];
			velocity_metric += state[i]*metric[i][j]*state[j];
		}
	}
	const auto temporal_scale = dt > 0.0 ? 4.0/(dt*dt) : 0.0;
	tau_m = 1.0/std::sqrt(temporal_scale + velocity_metric
		+ (1.0/12.0)*kinematic_viscosity*kinematic_viscosity*metric_norm);
	tau_c = 1.0/(tau_m*direction_norm);
}

iga::NavierStokesSystem LegacyTensorNavierStokesElement(const iga::Element& element,
	const std::vector<std::array<double, 4>>& nodal_state,
	const std::vector<std::array<double, 4>>& previous_nodal_state,
	const iga::NavierStokesParameters& parameters)
{
	const auto density = parameters.density;
	const auto viscosity = parameters.dynamic_viscosity;
	const auto kinematic_viscosity = viscosity/density;
	constexpr std::array<double, 4> points{{0.06943184420297371, 0.33000947820757187,
		0.6699905217924281, 0.9305681557970262}};
	constexpr std::array<double, 4> weights{{0.3478548451374539, 0.6521451548625461,
		0.6521451548625461, 0.3478548451374539}};
	const auto nen = element.connectivity.size();
	const auto ndof = 4*nen;
	iga::NavierStokesSystem system{std::vector<PetscScalar>(ndof*ndof, 0.0),
		std::vector<PetscScalar>(ndof, 0.0)};
	for (std::size_t qz = 0; qz < 4; ++qz)
		for (std::size_t qy = 0; qy < 4; ++qy)
			for (std::size_t qx = 0; qx < 4; ++qx) {
				auto basis = iga::EvaluateBasis(element, points[qx], points[qy], points[qz], true);
				const auto measure = weights[qx]*weights[qy]*weights[qz]*basis.determinant;
				std::array<double, 4> state{};
				std::array<double, 3> previous_velocity{};
				double gradient[4][3]{};
				double hessian[4][3][3]{};
				for (std::size_t a = 0; a < nen; ++a)
					for (int field = 0; field < 4; ++field) {
						state[field] += nodal_state[a][field]*basis.value[a];
						for (int i = 0; i < 3; ++i) {
							gradient[field][i] += nodal_state[a][field]*basis.gradient[a][i];
							for (int j = 0; j < 3; ++j)
								hessian[field][i][j] += nodal_state[a][field]*basis.hessian[a][i][j];
						}
					}
				if (parameters.dt > 0.0)
					for (std::size_t a = 0; a < nen; ++a)
						for (int component = 0; component < 3; ++component)
							previous_velocity[component] += previous_nodal_state[a][component]*basis.value[a];
				std::array<double, 3> time_derivative{};
				if (parameters.dt > 0.0)
					for (int component = 0; component < 3; ++component)
						time_derivative[component] = (state[component]-previous_velocity[component])
							/parameters.dt;
				double tau_m = 0.0;
				double tau_c = 0.0;
				LegacyStabilization(basis.inverse_jacobian, state, kinematic_viscosity,
					parameters.dt, tau_m, tau_c);
				std::array<double, 3> fine_velocity{};
				for (int component = 0; component < 3; ++component) {
					const auto convection = state[0]*gradient[component][0]
						+ state[1]*gradient[component][1] + state[2]*gradient[component][2];
					const auto pressure_gradient = gradient[3][component];
					const auto laplacian = hessian[component][0][0] + hessian[component][1][1]
						+ hessian[component][2][2];
					fine_velocity[component] = -tau_m*(time_derivative[component] + convection
						+ pressure_gradient/density - kinematic_viscosity*laplacian);
				}
				const auto fine_pressure = -density*tau_c
					*(gradient[0][0] + gradient[1][1] + gradient[2][2]);
				for (std::size_t a = 0; a < nen; ++a) {
					const auto na = basis.value[a];
					const auto& ga = basis.gradient[a];
					std::array<double, 4> residual{};
					for (int component = 0; component < 3; ++component) {
						residual[component] = density*na*time_derivative[component]
							- ga[component]*state[3] - ga[component]*fine_pressure;
						for (int direction = 0; direction < 3; ++direction) {
							residual[component] += viscosity*ga[direction]
								*(gradient[component][direction] + gradient[direction][component]);
							residual[component] += density*na*(state[direction] + fine_velocity[direction])
								*gradient[component][direction];
							residual[component] -= density*ga[direction]*fine_velocity[component]
								*(state[direction] + fine_velocity[direction]);
						}
					}
					residual[3] = na*(gradient[0][0] + gradient[1][1] + gradient[2][2])
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
						const auto diagonal = density*na*(mass_b + convection_b)
							+ viscosity*(ga[0]*gb[0] + ga[1]*gb[1] + ga[2]*gb[2])
							+ density*tau_m*streamline_a*(mass_b + convection_b);
						for (int i = 0; i < 3; ++i)
							for (int j = 0; j < 3; ++j)
								tangent[i][j] = viscosity*ga[j]*gb[i] + density*tau_c*ga[i]*gb[j];
						for (int i = 0; i < 3; ++i) tangent[i][i] += diagonal;
						for (int i = 0; i < 3; ++i) {
							tangent[i][3] = -ga[i]*nb + tau_m*streamline_a*gb[i];
							tangent[3][i] = na*gb[i] + tau_m*ga[i]*(mass_b + convection_b);
						}
						tangent[3][3] = (tau_m/density)
							*(ga[0]*gb[0] + ga[1]*gb[1] + ga[2]*gb[2]);
						for (int i = 0; i < 4; ++i)
							for (int j = 0; j < 4; ++j)
								system.jacobian[(4*a+i)*ndof + 4*b+j] += tangent[i][j]*measure;
					}
				}
			}
	return system;
}

void RequireEqual(const std::vector<PetscScalar>& expected,
	const std::vector<PetscScalar>& actual, const char* name)
{
	assert(expected.size() == actual.size());
	double largest = 0.0;
	for (std::size_t i = 0; i < expected.size(); ++i)
		largest = std::max(largest, std::abs(PetscRealPart(expected[i]-actual[i])));
	if (largest > 2e-14) {
		std::cerr << name << " maximum difference " << largest << '\n';
		std::abort();
	}
}

} // namespace

int main()
{
	const auto element = UnitElement();
	std::vector<std::array<double, 4>> state(64);
	for (std::size_t a = 0; a < state.size(); ++a)
		state[a] = {0.1 + 0.001*a, -0.02, 0.03, 0.04 - 0.0002*a};
	const auto legacy = iga::BuildNavierStokesElement(element, state, 0.17);
	iga::FullCell4x4x4VolumeQuadratureProvider unit_quadrature(element);
	const auto explicit_steady = iga::BuildNavierStokesElement(
		element, state, {}, {1.0, 0.17, 0.0}, unit_quadrature.Rule());
	const auto frozen_steady = LegacyTensorNavierStokesElement(element, state, {},
		{1.0, 0.17, 0.0});
	RequireEqual(frozen_steady.negative_residual, explicit_steady.negative_residual,
		"frozen steady negative residual");
	RequireEqual(legacy.jacobian, explicit_steady.jacobian, "steady wrapper Jacobian");
	RequireEqual(legacy.negative_residual, frozen_steady.negative_residual,
		"steady wrapper negative residual");
	// The conservative resolved pair is explicitly opt-in.  Its all-label
	// physical trace recovers the legacy body-fitted pair exactly on this
	// polynomial cube, while the resolved J_pu and J_up blocks remain skew.
	const auto conservative = iga::BuildNavierStokesElementFromPoints(element, state, {},
		{1.0, 0.17, 0.0}, [&unit_quadrature](const auto& consume) {
			for (const auto& point : unit_quadrature.Rule().Points()) consume(point);
		}, [](const std::array<double, 3>&) { return std::array<double, 3>{{0.0,0.0,0.0}}; },
		iga::NavierStokesResolvedMixedForm::Conservative);
	const auto trace = iga::BuildImmersedConservativeMixedTraceElement(element,
		UnitSurfaceQuadrature(), state);
	auto completed_conservative = conservative;
	AddSystem(completed_conservative, trace);
	RequireEqual(explicit_steady.jacobian, completed_conservative.jacobian,
		"conservative mixed trace Jacobian completion");
	RequireEqual(explicit_steady.negative_residual, completed_conservative.negative_residual,
		"conservative mixed trace residual completion");
	constexpr std::size_t ndof = 256;
	for (std::size_t a = 0; a < 64; ++a)
		for (std::size_t b = 0; b < 64; ++b)
			for (int component = 0; component < 3; ++component) {
				const auto pressure_velocity = PetscRealPart(trace.jacobian[(4*a+3)*ndof+4*b+component]);
				const auto velocity_pressure = PetscRealPart(trace.jacobian[(4*b+component)*ndof+4*a+3]);
				assert(std::abs(pressure_velocity+velocity_pressure) < 2e-12);
			}
	auto curved = element;
	for (auto& point : curved.bezier_points) {
		point[0] += 0.08*point[0]*point[1];
		point[1] += 0.05*point[1]*point[2];
	}
	iga::FullCell4x4x4VolumeQuadratureProvider curved_quadrature(curved);
	const auto curved_legacy = iga::BuildNavierStokesElement(
		curved, state, {}, {1.0, 0.17, 0.0});
	const auto curved_explicit = iga::BuildNavierStokesElement(
		curved, state, {}, {1.0, 0.17, 0.0}, curved_quadrature.Rule());
	const auto frozen_curved_steady = LegacyTensorNavierStokesElement(curved, state, {},
		{1.0, 0.17, 0.0});
	RequireEqual(frozen_curved_steady.negative_residual, curved_explicit.negative_residual,
		"frozen curved steady negative residual");
	RequireEqual(curved_legacy.jacobian, curved_explicit.jacobian,
		"curved steady wrapper Jacobian");
	RequireEqual(curved_legacy.negative_residual, frozen_curved_steady.negative_residual,
		"curved steady wrapper negative residual");
	std::vector<std::array<double, 4>> nonconstant_previous(64);
	for (std::size_t a = 0; a < nonconstant_previous.size(); ++a)
		nonconstant_previous[a] = {0.04-0.0007*a, -0.01+0.0003*a,
			0.015-0.0002*a, 0.0};
	const iga::NavierStokesParameters nonconstant_transient_parameters{2.5, 0.17, 0.2};
	const auto nonconstant_transient = iga::BuildNavierStokesElement(element, state,
		nonconstant_previous, nonconstant_transient_parameters, unit_quadrature.Rule());
	const auto frozen_nonconstant_transient = LegacyTensorNavierStokesElement(element,
		state, nonconstant_previous, nonconstant_transient_parameters);
	RequireEqual(frozen_nonconstant_transient.negative_residual,
		nonconstant_transient.negative_residual, "frozen nonconstant transient negative residual");
	const auto curved_nonconstant_transient = iga::BuildNavierStokesElement(curved, state,
		nonconstant_previous, nonconstant_transient_parameters, curved_quadrature.Rule());
	const auto frozen_curved_nonconstant_transient = LegacyTensorNavierStokesElement(curved,
		state, nonconstant_previous, nonconstant_transient_parameters);
	RequireEqual(frozen_curved_nonconstant_transient.negative_residual,
		curved_nonconstant_transient.negative_residual,
		"frozen curved nonconstant transient negative residual");
	// The VMS residual retains the strong viscous term in fine_velocity.  Its
	// tangent must therefore include that term, as well as the velocity
	// dependence of tau_m and tau_c; verify the complete nonlinear action.
	std::vector<std::array<double, 4>> fd_plus = state, fd_minus = state;
	std::vector<double> fd_direction(256);
	constexpr double fd_epsilon = 1e-4;
	for (std::size_t a = 0; a < state.size(); ++a)
		for (int field = 0; field < 4; ++field) {
			const double direction = 1e-3*(1.0+static_cast<double>((7*a+3*field)%13));
			fd_direction[4*a+field] = direction;
			fd_plus[a][field] += fd_epsilon*direction;
			fd_minus[a][field] -= fd_epsilon*direction;
		}
	const auto fd_plus_system = iga::BuildNavierStokesElement(element, fd_plus,
		nonconstant_previous, nonconstant_transient_parameters, unit_quadrature.Rule());
	const auto fd_minus_system = iga::BuildNavierStokesElement(element, fd_minus,
		nonconstant_previous, nonconstant_transient_parameters, unit_quadrature.Rule());
	double action_squared = 0.0, defect_squared = 0.0;
	for (std::size_t row = 0; row < 256; ++row) {
		double action = 0.0;
		for (std::size_t column = 0; column < 256; ++column)
			action += PetscRealPart(nonconstant_transient.jacobian[row*256+column])*fd_direction[column];
		const double defect = action+(PetscRealPart(fd_plus_system.negative_residual[row])
			-PetscRealPart(fd_minus_system.negative_residual[row]))/(2.0*fd_epsilon);
		action_squared += action*action;
		defect_squared += defect*defect;
	}
	const double nonlinear_fd_defect = std::sqrt(defect_squared/action_squared);
	std::cerr << "nonlinear Navier-Stokes Jacobian relative defect " << nonlinear_fd_defect << '\n';
	assert(nonlinear_fd_defect <= 1e-8);

	constexpr double density = 2.5;
	constexpr double dt = 0.2;
	const std::array<double, 3> current{{1.0, -0.5, 0.25}};
	const std::array<double, 3> previous{{0.6, -0.1, 0.05}};
	for (auto& node : state) node = {current[0], current[1], current[2], 0.0};
	std::vector<std::array<double, 4>> previous_state(64);
	for (auto& node : previous_state) node = {previous[0], previous[1], previous[2], 0.0};
	const auto transient = iga::BuildNavierStokesElement(
		element, state, previous_state, {density, 0.17, dt});
	const auto transient_explicit = iga::BuildNavierStokesElement(
		element, state, previous_state, {density, 0.17, dt}, unit_quadrature.Rule());
	const auto frozen_transient = LegacyTensorNavierStokesElement(element, state,
		previous_state, {density, 0.17, dt});
	RequireEqual(frozen_transient.negative_residual, transient_explicit.negative_residual,
		"frozen transient negative residual");
	RequireEqual(transient.jacobian, transient_explicit.jacobian, "transient wrapper Jacobian");
	RequireEqual(transient.negative_residual, frozen_transient.negative_residual,
		"transient wrapper negative residual");
	for (int component = 0; component < 3; ++component) {
		double negative_residual_sum = 0.0;
		double mass_block_sum = 0.0;
		for (std::size_t a = 0; a < 64; ++a) {
			negative_residual_sum += PetscRealPart(transient.negative_residual[4*a+component]);
			for (std::size_t b = 0; b < 64; ++b)
				mass_block_sum += PetscRealPart(
					transient.jacobian[(4*a+component)*256+4*b+component]);
		}
		const auto expected_residual = -density*(current[component]-previous[component])/dt;
		assert(std::abs(negative_residual_sum-expected_residual) < 2e-11);
		assert(std::abs(mass_block_sum-density/dt) < 2e-11);
	}
	iga::VolumeQuadratureRule single_point({{{{0.5, 0.5, 0.5}}, 0.25}});
	const auto reduced = iga::BuildNavierStokesElement(
		element, state, previous_state, {density, 0.17, dt}, single_point);
	double reduced_mass_block_sum = 0.0;
	for (std::size_t a = 0; a < 64; ++a)
		for (std::size_t b = 0; b < 64; ++b)
			reduced_mass_block_sum += PetscRealPart(reduced.jacobian[4*a*256+4*b]);
	assert(std::abs(reduced_mass_block_sum-0.25*density/dt) < 2e-11);
	assert(std::abs(reduced_mass_block_sum-density/dt) > 1.0);
	const std::array<double, 3> off_center{{0.23, 0.61, 0.37}};
	const iga::VolumeQuadratureRule off_center_rule(
		std::vector<iga::VolumeQuadraturePoint>{{off_center, 0.41}});
	const auto off_center_basis = iga::EvaluateBasis(element, off_center[0], off_center[1],
		off_center[2], true);
	constexpr double off_center_velocity = 0.15;
	std::array<double, 4> off_center_state{{off_center_velocity, 0.0, 0.0, 0.0}};
	double tau_m = 0.0;
	double tau_c = 0.0;
	iga::Stabilization(off_center_basis.inverse_jacobian, off_center_state, 0.17/density,
		dt, tau_m, tau_c);
	const auto off_center_derivative = off_center_velocity/tau_m;
	std::vector<std::array<double, 4>> off_center_current(64, off_center_state);
	std::vector<std::array<double, 4>> off_center_previous(64,
		{off_center_velocity-dt*off_center_derivative, 0.0, 0.0, 0.0});
	const auto off_center_system = iga::BuildNavierStokesElement(element,
		off_center_current, off_center_previous, {density, 0.17, dt}, off_center_rule);
	for (std::size_t a = 0; a < 64; ++a) {
		const auto expected = -density*off_center_basis.value[a]*off_center_derivative
			*0.41*off_center_basis.raw_determinant;
		assert(std::abs(PetscRealPart(off_center_system.negative_residual[4*a])-expected)
			< 2e-13);
	}
	for (auto& node : state) node = {2.0, 0.0, 0.0, 0.0};
	assert(std::abs(iga::IntegrateBoundaryFlow(element, 2, state)-2.0) < 2e-13);
	assert(std::abs(iga::IntegrateBoundaryFlow(element, 4, state)+2.0) < 2e-13);
	std::cout << "Navier-Stokes steady compatibility and backward-Euler mass tests passed\n";
}
