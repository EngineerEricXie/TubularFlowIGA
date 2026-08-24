#include "CaseInput.hpp"
#include "GenericCaseInput.hpp"
#include "GenericTransportElement.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
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
	iga::TransportParameters parameters;
	parameters.diffusion = 0.17;
	parameters.vplus = 1.3;
	parameters.kplus = 0.4;
	parameters.kminus = 0.2;
	parameters.detach_plus = 0.07;
	parameters.dt = 0.015;
	parameters.steps = 3;
	parameters.artificial_diffusion = 0.09;
	const auto element = UnitElement();
	std::vector<std::array<double, 3>> velocity(64);
	for (std::size_t i = 0; i < velocity.size(); ++i)
		velocity[i] = {0.2 + 0.001*i, -0.03, 0.05};
	const auto legacy = iga::BuildTransportElement(element, velocity, parameters);
	const auto configuration = iga::ConvertLegacyNeuronTransport(parameters);
	const auto system = iga::CompileLinearSystem(configuration, "neuron_transport");
	const auto generic = iga::BuildGenericTransportElement(element, velocity, system);
	RequireEqual(legacy.left, generic.left, "left matrix");
	RequireEqual(legacy.previous, generic.previous, "previous matrix");
	std::cout << "generic neuron transport regression passed\n";
	iga::SimulationConfiguration custom;
	custom.fields = {{"oxygen", iga::FieldKind::Scalar, 0.0},
		{"drug", iga::FieldKind::Scalar, 0.0}, {"metabolite", iga::FieldKind::Scalar, 0.0}};
	custom.time = {0.02, 1};
	iga::EquationSystemDefinition custom_definition;
	custom_definition.name = "custom_three_field";
	custom_definition.kind = iga::EquationKind::LinearTransport;
	custom_definition.unknowns = {"oxygen", "drug", "metabolite"};
	custom_definition.terms = {
		{iga::TermKind::TimeDerivative, "oxygen", "oxygen", 1.0, ""},
		{iga::TermKind::TimeDerivative, "drug", "drug", 1.0, ""},
		{iga::TermKind::TimeDerivative, "metabolite", "metabolite", 1.0, ""},
		{iga::TermKind::LinearCoupling, "oxygen", "metabolite", -0.5, ""}
	};
	custom.equation_systems.push_back(custom_definition);
	const auto custom_system = iga::CompileLinearSystem(custom, "custom_three_field");
	const auto custom_matrices = iga::BuildGenericTransportElement(element, velocity, custom_system);
	assert(custom_matrices.left.size() == 192*192);
	assert(std::abs(PetscRealPart(custom_matrices.left[2])) > 0.0);
	std::cout << "custom three-field assembly passed\n";
	custom.boundaries = {{0, "wall", {{"oxygen", iga::FieldBoundaryKind::Flux,
		{1.0}, 0.0, 0.0, "", 1.0}}}};
	bool rejected_surface = false;
	try {
		iga::ResolveScalarBoundaries(custom, custom_system, std::vector<int>(64, 0));
	} catch (const std::runtime_error&) {
		rejected_surface = true;
	}
	assert(rejected_surface);
	std::cout << "unsupported scalar surface condition rejected\n";
}
