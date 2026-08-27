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

std::vector<PetscScalar> Expand(const iga::FieldBlockElementMatrix& blocks)
{
	const auto fields = blocks.pattern().fields();
	const auto nodes = blocks.nodes();
	const auto dofs = fields*nodes;
	std::vector<PetscScalar> result(dofs*dofs, 0.0);
	for (const auto& pair : blocks.pattern().active_pairs())
		for (std::size_t a = 0; a < nodes; ++a)
			for (std::size_t b = 0; b < nodes; ++b)
				result[(a*fields+pair.first)*dofs+b*fields+pair.second]
					= blocks.At(pair.first, pair.second, a, b);
	return result;
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
	assert(!system.stabilization.empty());
	const auto generic = iga::BuildGenericTransportElement(element, velocity, system, configuration);
	RequireEqual(legacy.left, Expand(generic.left), "left matrix");
	RequireEqual(legacy.previous, Expand(generic.previous), "previous matrix");
	assert(generic.left.pattern().pairs() == 4);
	assert(generic.previous.pattern().pairs() == 2);
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
	auto diagonal_only = custom;
	diagonal_only.equation_systems[0].terms.pop_back();
	const auto diagonal_system = iga::CompileLinearSystem(diagonal_only, "custom_three_field");
	const auto diagonal_patterns = iga::BuildTransportCouplingPatterns(
		diagonal_system, diagonal_only);
	assert(diagonal_patterns.left.pairs() == 3);
	assert(diagonal_patterns.previous.pairs() == 3);
	const auto custom_system = iga::CompileLinearSystem(custom, "custom_three_field");
	const auto custom_matrices = iga::BuildGenericTransportElement(element, velocity, custom_system, custom);
	assert(custom_matrices.left.pattern().pairs() == 4);
	assert(custom_matrices.previous.pattern().pairs() == 3);
	assert(custom_matrices.left.values().size() == 4*64*64);
	assert(custom_matrices.previous.values().size() == 3*64*64);
	assert(custom_matrices.left.pattern().Active(0, 2));
	assert(!custom_matrices.left.pattern().Active(2, 0));
	assert(std::abs(PetscRealPart(custom_matrices.left.At(0, 2, 0, 0))) > 0.0);
	auto two_way = custom;
	two_way.equation_systems[0].terms.push_back(
		{iga::TermKind::LinearCoupling, "metabolite", "oxygen", 0.25, ""});
	const auto two_way_system = iga::CompileLinearSystem(two_way, "custom_three_field");
	const auto two_way_patterns = iga::BuildTransportCouplingPatterns(two_way_system, two_way);
	assert(two_way_patterns.left.pairs() == 5);
	assert(two_way_patterns.left.Active(0, 2));
	assert(two_way_patterns.left.Active(2, 0));
	std::cout << "custom three-field assembly passed\n";
	custom.boundaries = {{0, "wall", {
		{"oxygen", iga::FieldBoundaryKind::Flux, {1.0}, 0.0, 0.0, "", 1.0, ""},
		{"drug", iga::FieldBoundaryKind::Robin, {}, 2.0, 3.0, "", 1.0, ""}}}};
	auto robin_only = custom;
	auto& robin_terms = robin_only.equation_systems[0].terms;
	robin_terms.erase(std::remove_if(robin_terms.begin(), robin_terms.end(),
		[](const auto& term) { return term.equation == "drug"; }), robin_terms.end());
	const auto robin_system = iga::CompileLinearSystem(robin_only, "custom_three_field");
	const auto robin_patterns = iga::BuildTransportCouplingPatterns(robin_system, robin_only);
	assert(robin_patterns.left.Active(1, 1));
	auto dirichlet_only = robin_only;
	dirichlet_only.boundaries[0].conditions[0].kind = iga::FieldBoundaryKind::Dirichlet;
	const auto dirichlet_system = iga::CompileLinearSystem(dirichlet_only, "custom_three_field");
	const auto dirichlet_patterns = iga::BuildTransportCouplingPatterns(
		dirichlet_system, dirichlet_only);
	assert(dirichlet_patterns.left.Active(0, 0));
	auto surface_element = element;
	surface_element.boundary_labels[0] = 0;
	const auto surface_matrices = iga::BuildGenericTransportElement(
		surface_element, velocity, custom_system, custom);
	double oxygen_source = 0.0;
	double drug_source = 0.0;
	double drug_left = 0.0;
	for (std::size_t a = 0; a < 64; ++a) {
		oxygen_source += PetscRealPart(surface_matrices.source[3*a]
			- custom_matrices.source[3*a]);
		drug_source += PetscRealPart(surface_matrices.source[3*a+1]
			- custom_matrices.source[3*a+1]);
		for (std::size_t b = 0; b < 64; ++b)
			drug_left += PetscRealPart(surface_matrices.left.At(1, 1, a, b)
				- custom_matrices.left.At(1, 1, a, b));
	}
	assert(std::abs(oxygen_source - 0.02) < 2e-13);
	assert(std::abs(drug_source - 0.12) < 2e-13);
	assert(std::abs(drug_left - 0.04) < 2e-13);
	const auto resolved = iga::ResolveScalarBoundaries(custom, custom_system, std::vector<int>(64, 0));
	assert(resolved.constrained_dofs == 0);
	std::cout << "scalar flux and Robin surface assembly passed\n";
}
