#include "BoundaryFlow.hpp"
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
	const auto explicit_steady = iga::BuildNavierStokesElement(
		element, state, {}, {1.0, 0.17, 0.0});
	RequireEqual(legacy.jacobian, explicit_steady.jacobian, "steady Jacobian");
	RequireEqual(legacy.negative_residual, explicit_steady.negative_residual,
		"steady negative residual");

	constexpr double density = 2.5;
	constexpr double dt = 0.2;
	const std::array<double, 3> current{{1.0, -0.5, 0.25}};
	const std::array<double, 3> previous{{0.6, -0.1, 0.05}};
	for (auto& node : state) node = {current[0], current[1], current[2], 0.0};
	std::vector<std::array<double, 4>> previous_state(64);
	for (auto& node : previous_state) node = {previous[0], previous[1], previous[2], 0.0};
	const auto transient = iga::BuildNavierStokesElement(
		element, state, previous_state, {density, 0.17, dt});
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
	for (auto& node : state) node = {2.0, 0.0, 0.0, 0.0};
	assert(std::abs(iga::IntegrateBoundaryFlow(element, 2, state)-2.0) < 2e-13);
	assert(std::abs(iga::IntegrateBoundaryFlow(element, 4, state)+2.0) < 2e-13);
	std::cout << "Navier-Stokes steady compatibility and backward-Euler mass tests passed\n";
}
