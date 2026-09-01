#include "BoundaryFlow.hpp"
#include "OutletFlow.hpp"
#include "PressureTractionFlow.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

int main()
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
	std::vector<std::array<double, 4>> state(64);
	for (auto& node : state) node = {2.0, 0.0, 0.0, 0.0};
	const auto interior_basis = iga::EvaluateBoundaryBasis(element, 0.2, 0.3, 0.4);
	assert(std::abs(interior_basis.physical_coordinate[0]-0.2) < 2e-13);
	assert(std::abs(interior_basis.physical_coordinate[1]-0.3) < 2e-13);
	assert(std::abs(interior_basis.physical_coordinate[2]-0.4) < 2e-13);
	assert(std::abs(iga::IntegrateBoundaryFlow(element, 2, state)-2.0) < 2e-13);
	assert(std::abs(iga::IntegrateBoundaryFlow(element, 4, state)+2.0) < 2e-13);
	for (auto& node : state) node[3] = 5.0;
	const auto pressure_area = iga::IntegrateBoundaryScalarAndArea(element, 2, state);
	assert(std::abs(pressure_area[0]-5.0) < 2e-13);
	assert(std::abs(pressure_area[1]-1.0) < 2e-13);
	for (auto& node : state) node[3] = 0.0;
	assert(std::abs(iga::IntegrateVolumeDivergence(element, state)) < 2e-13);
	assert(std::abs(iga::IntegrateBoundarySpeciesFlux(element, 2, state,
		std::vector<double>(64, 3.0))-6.0) < 2e-13);
	assert(std::abs(iga::IntegrateBoundarySpeciesFlux(element, 4, state,
		std::vector<double>(64, 3.0))+6.0) < 2e-13);
	std::vector<std::vector<double>> linear_species(64, std::vector<double>(1));
	for (std::size_t node = 0; node < linear_species.size(); ++node)
		linear_species[node][0] = element.bezier_points[node][0];
	const auto right_species = iga::IntegrateBoundaryTransportFlux(
		element, 2, state, linear_species, 0, {1.0}, {0.5});
	const auto left_species = iga::IntegrateBoundaryTransportFlux(
		element, 4, state, linear_species, 0, {1.0}, {0.5});
	assert(std::abs(right_species.concentration_integral-1.0) < 2e-13);
	assert(std::abs(right_species.total_outward_flux-1.5) < 2e-13);
	assert(std::abs(left_species.concentration_integral) < 2e-13);
	assert(std::abs(left_species.total_outward_flux-0.5) < 2e-13);
	auto scaled = element;
	for (auto& point : scaled.bezier_points) {
		point[0] *= 2.0;
		point[1] *= 3.0;
		point[2] *= 4.0;
	}
	std::vector<std::array<double, 4>> scaled_state(64, {2.0, 0.0, 0.0, 5.0});
	const auto scaled_pressure_area = iga::IntegrateBoundaryScalarAndArea(
		scaled, 2, scaled_state);
	assert(std::abs(scaled_pressure_area[0]-60.0) < 2e-12);
	assert(std::abs(scaled_pressure_area[1]-12.0) < 2e-12);
	std::vector<std::vector<double>> scaled_species(64, std::vector<double>(1));
	for (std::size_t node = 0; node < scaled_species.size(); ++node)
		scaled_species[node][0] = scaled.bezier_points[node][0];
	const auto scaled_right = iga::IntegrateBoundaryTransportFlux(
		scaled, 2, scaled_state, scaled_species, 0, {1.0}, {0.5});
	const auto scaled_left = iga::IntegrateBoundaryTransportFlux(
		scaled, 4, scaled_state, scaled_species, 0, {1.0}, {0.5});
	assert(std::abs(scaled_right.concentration_integral-24.0) < 2e-12);
	assert(std::abs(scaled_right.total_outward_flux-42.0) < 2e-12);
	assert(std::abs(scaled_left.concentration_integral) < 2e-12);
	assert(std::abs(scaled_left.total_outward_flux-6.0) < 2e-12);
	for (std::size_t node = 0; node < state.size(); ++node)
		state[node] = {element.bezier_points[node][0], 0.0, 0.0, 0.0};
	const auto linear_surface = iga::IntegrateBoundaryFlow(element, 2, state)
		+iga::IntegrateBoundaryFlow(element, 4, state);
	const auto linear_volume = iga::IntegrateVolumeDivergence(element, state);
	assert(std::abs(linear_surface-1.0) < 2e-13);
	assert(std::abs(linear_volume-1.0) < 2e-13);
	assert(std::abs(linear_surface-linear_volume) < 2e-13);
	const auto right_traction = iga::IntegrateBoundaryPressureTraction(element, 2, 2.0);
	const auto left_traction = iga::IntegrateBoundaryPressureTraction(element, 4, 2.0);
	double right_x = 0.0, left_x = 0.0;
	for (std::size_t node = 0; node < 64; ++node) {
		right_x += right_traction[4*node];
		left_x += left_traction[4*node];
		for (int component = 1; component < 4; ++component) {
			assert(std::abs(right_traction[4*node+component]) < 2e-13);
			assert(std::abs(left_traction[4*node+component]) < 2e-13);
		}
	}
	assert(std::abs(right_x+2.0) < 2e-13);
	assert(std::abs(left_x-2.0) < 2e-13);
	element.boundary_labels[2] = 2;
	element.boundary_labels[4] = 3;
	iga::BodyFittedSurface4x4QuadratureProvider surface_quadrature(element);
	iga::FullCell4x4x4VolumeQuadratureProvider volume_quadrature(element);
	for (auto& node : state) node = {2.0, 0.0, 0.0, 5.0};
	assert(std::abs(iga::IntegrateBoundaryFlow(element, state, surface_quadrature.Rule(), 2)
		-iga::IntegrateBoundaryFlow(element, 2, state)) < 2e-13);
	assert(std::abs(iga::IntegrateBoundaryFlow(element, state, surface_quadrature.Rule(), 3)
		-iga::IntegrateBoundaryFlow(element, 4, state)) < 2e-13);
	const auto explicit_pressure_area = iga::IntegrateBoundaryScalarAndArea(
		element, state, surface_quadrature.Rule(), 2);
	const auto legacy_pressure_area = iga::IntegrateBoundaryScalarAndArea(element, 2, state);
	assert(std::abs(explicit_pressure_area[0]-legacy_pressure_area[0]) < 2e-13
		&& std::abs(explicit_pressure_area[1]-legacy_pressure_area[1]) < 2e-13);
	const auto explicit_traction = iga::IntegrateBoundaryPressureTraction(
		element, surface_quadrature.Rule(), 2, 2.0);
	for (std::size_t row = 0; row < right_traction.size(); ++row)
		assert(std::abs(explicit_traction[row]-right_traction[row]) < 2e-13);
	std::vector<std::array<double, 4>> divergence_state(64);
	for (std::size_t node = 0; node < divergence_state.size(); ++node)
		divergence_state[node] = {element.bezier_points[node][0], 0.0, 0.0, 0.0};
	assert(std::abs(iga::IntegrateVolumeDivergence(element, divergence_state,
		volume_quadrature.Rule())-iga::IntegrateVolumeDivergence(element, divergence_state))
		< 2e-13);
	auto permuted_points = surface_quadrature.Rule().Points();
	std::reverse(permuted_points.begin(), permuted_points.end());
	const iga::SurfaceQuadratureRule permuted(std::move(permuted_points));
	assert(std::abs(iga::IntegrateBoundaryFlow(element, state, permuted, 2)
		-iga::IntegrateBoundaryFlow(element, state, surface_quadrature.Rule(), 2)) < 2e-13);
	const std::array<double, 3> interior_coordinate{{0.27, 0.44, 0.63}};
	const auto synthetic_basis = iga::EvaluateBoundaryBasis(element, interior_coordinate[0],
		interior_coordinate[1], interior_coordinate[2]);
	iga::SurfaceQuadraturePoint interior_point;
	interior_point.parametric = interior_coordinate;
	interior_point.physical = synthetic_basis.physical_coordinate;
	interior_point.normal = {{0.6, 0.8, 0.0}};
	interior_point.weight = 0.37;
	interior_point.boundary_id = 71;
	const iga::SurfaceQuadratureRule interior_surface({interior_point});
	assert(std::abs(iga::IntegrateBoundaryFlow(element, state, interior_surface, 71)
		-2.0*0.6*0.37) < 2e-13);
	const auto interior_scalar_area = iga::IntegrateBoundaryScalarAndArea(element, state,
		interior_surface, 71);
	assert(std::abs(interior_scalar_area[0]-5.0*0.37) < 2e-13);
	assert(std::abs(interior_scalar_area[1]-0.37) < 2e-13);
	const auto interior_traction = iga::IntegrateBoundaryPressureTraction(element,
		interior_surface, 71, 7.0);
	for (std::size_t a = 0; a < 64; ++a) {
		assert(std::abs(interior_traction[4*a]
			+7.0*0.37*0.6*synthetic_basis.value[a]) < 2e-13);
		assert(std::abs(interior_traction[4*a+1]
			+7.0*0.37*0.8*synthetic_basis.value[a]) < 2e-13);
		assert(std::abs(interior_traction[4*a+2]) < 2e-13);
	}
	std::vector<double> interior_species(64);
	std::vector<std::vector<double>> interior_species_fields(64, std::vector<double>(1));
	for (std::size_t a = 0; a < 64; ++a) {
		interior_species[a] = element.bezier_points[a][0];
		interior_species_fields[a][0] = interior_species[a];
	}
	assert(std::abs(iga::IntegrateBoundarySpeciesFlux(element, state, interior_species,
		interior_surface, 71)-0.37*2.0*0.6*0.27) < 2e-13);
	const auto interior_transport = iga::IntegrateBoundaryTransportFlux(element, state,
		interior_species_fields, 0, {1.0}, {0.5}, interior_surface, 71);
	assert(std::abs(interior_transport.concentration_integral-0.37*0.27) < 2e-13);
	assert(std::abs(interior_transport.total_outward_flux-0.37*0.6*(2.0*0.27-0.5))
		< 2e-13);
	auto repeated = element;
	repeated.boundary_labels[4] = 2;
	iga::BodyFittedSurface4x4QuadratureProvider repeated_quadrature(repeated);
	assert(std::abs(iga::IntegrateBoundaryFlow(repeated, state, repeated_quadrature.Rule(), 2)
		-(iga::IntegrateBoundaryFlow(repeated, 2, state)
			+iga::IntegrateBoundaryFlow(repeated, 4, state))) < 2e-13);
	auto curved = element;
	for (auto& point : curved.bezier_points) {
		point[0] += 0.08*point[0]*point[1] + 0.03*point[2];
		point[1] += 0.05*point[1]*point[2];
		point[2] += 0.04*point[0]*point[2];
	}
	iga::BodyFittedSurface4x4QuadratureProvider curved_surface_quadrature(curved);
	iga::FullCell4x4x4VolumeQuadratureProvider curved_volume_quadrature(curved);
	std::vector<std::array<double, 4>> curved_state(64, {2.0, 0.0, 0.0, 5.0});
	assert(std::abs(iga::IntegrateBoundaryFlow(curved, curved_state,
		curved_surface_quadrature.Rule(), 2)-iga::IntegrateBoundaryFlow(curved, 2,
		curved_state)) < 2e-12);
	const auto curved_scalar_area = iga::IntegrateBoundaryScalarAndArea(curved,
		curved_state, curved_surface_quadrature.Rule(), 2);
	const auto curved_legacy_scalar_area = iga::IntegrateBoundaryScalarAndArea(
		curved, 2, curved_state);
	assert(std::abs(curved_scalar_area[0]-curved_legacy_scalar_area[0]) < 2e-12
		&& std::abs(curved_scalar_area[1]-curved_legacy_scalar_area[1]) < 2e-12);
	const auto curved_traction = iga::IntegrateBoundaryPressureTraction(curved,
		curved_surface_quadrature.Rule(), 2, 2.0);
	const auto curved_legacy_traction = iga::IntegrateBoundaryPressureTraction(curved, 2, 2.0);
	for (std::size_t row = 0; row < curved_traction.size(); ++row)
		assert(std::abs(curved_traction[row]-curved_legacy_traction[row]) < 2e-12);
	std::vector<std::array<double, 4>> curved_divergence_state(64);
	for (std::size_t node = 0; node < curved_divergence_state.size(); ++node)
		curved_divergence_state[node] = {curved.bezier_points[node][0], 0.0, 0.0, 0.0};
	assert(std::abs(iga::IntegrateVolumeDivergence(curved, curved_divergence_state,
		curved_volume_quadrature.Rule())-iga::IntegrateVolumeDivergence(curved,
		curved_divergence_state)) < 2e-12);
	const auto global_traction = iga::IntegratePressureTractionForces(
		{{2, 2.0}, {3, 2.0}}, {element}, 64);
	double global_x = 0.0;
	for (std::size_t node = 0; node < 64; ++node)
		global_x += global_traction[4*node];
	assert(std::abs(global_x) < 2e-13);
	for (auto& node : state) node = {2.0, 0.0, 0.0, 0.0};
	iga::OutletModelState right, left;
	right.label = 2;
	left.label = 3;
	std::vector<double> flat_state(64*4);
	for (std::size_t node = 0; node < 64; ++node)
		for (int field = 0; field < 4; ++field)
			flat_state[4*node+static_cast<std::size_t>(field)] = state[node][field];
	const auto flows = iga::IntegrateOutletModelFlows(
		{right, left}, {element}, flat_state);
	assert(std::abs(flows[0]-2.0) < 2e-13);
	assert(std::abs(flows[1]+2.0) < 2e-13);
	std::cout << "dependency-free boundary flow and pressure traction tests passed\n";
}
