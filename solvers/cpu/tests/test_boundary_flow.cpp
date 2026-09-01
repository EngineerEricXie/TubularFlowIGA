#include "BoundaryFlow.hpp"
#include "OutletFlow.hpp"
#include "PressureTractionFlow.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
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
