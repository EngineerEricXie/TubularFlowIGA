#include "TransportBudget.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {

int checks = 0;

void Require(bool condition, const char* message)
{
	++checks;
	if (!condition) throw std::runtime_error(message);
}

template <class Action>
void Reject(Action&& action)
{
	bool rejected = false;
	try { action(); } catch (const std::exception&) { rejected = true; }
	Require(rejected, "invalid budget operation was accepted");
}

void Near(long double value, long double expected, long double tolerance = 1e-12L)
{
	Require(std::abs(value-expected) <= tolerance, "analytic transport budget differs");
}

iga::Element Cube(double x = 1.0, double y = 1.0, double z = 1.0)
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
				element.bezier_points[p] = {x*i/3.0, y*j/3.0, z*k/3.0};
			}
	return element;
}

struct Fixture {
	iga::CompiledLinearSystem system;
	iga::SimulationConfiguration configuration;
	iga::ResolvedScalarBoundaries boundaries;
	std::vector<double> previous = std::vector<double>(64, 1.0);
	std::vector<double> current = std::vector<double>(64, 2.0);
	std::vector<std::array<double, 3>> velocity = std::vector<std::array<double, 3>>(64);
	iga::Element element = Cube();

	Fixture()
	{
		system.name = "budget-test";
		system.fields = {"c"};
		system.field_index = {{"c", 0}};
		system.dt = 0.25;
		system.steps = 1;
		system.terms = {{iga::TermKind::TimeDerivative, 0, 0, 1.0, ""},
			{iga::TermKind::VolumeSource, 0, 0, 4.0, ""}};
		boundaries.constrained.assign(64, 0);
		boundaries.value.assign(64, 0.0);
	}

	iga::TransportBudgetAccumulator Accumulator()
	{
		return {64, 1, system, configuration, boundaries, previous, current, velocity};
	}

	iga::TransportFieldBudget Evaluate()
	{
		auto accumulator = Accumulator();
		accumulator.Add(element);
		return accumulator.Finish().at(0);
	}
};

void Passed(const iga::TransportFieldBudget& budget)
{
	Require(iga::TransportBudgetWithin(budget.volume_defect, budget.volume_scale, 1e-6, 1e-12), "volume budget failed");
	Require(iga::TransportBudgetWithin(budget.surface_defect, budget.surface_scale, 1e-6, 1e-12), "surface budget failed");
	Require(iga::TransportBudgetWithin(budget.free_residual_l2, budget.free_component_scale, 1e-6, 1e-12), "free residual failed");
	Require(iga::TransportBudgetWithin(budget.essential_error_l2, budget.essential_reference_l2, 1e-6, 1e-12), "essential values failed");
}

void SourceAndScaling()
{
	Fixture fixture;
	auto budget = fixture.Evaluate();
	Passed(budget);
	Near(budget.previous_integral, 1.0L);
	Near(budget.current_integral, 2.0L);
	Near(budget.capacity_rate, 4.0L);
	Near(budget.source, 4.0L);
	Require(budget.free_rows == 64, "unexpected free-row count");
	fixture.element = Cube(2.0, 3.0, 4.0);
	budget = fixture.Evaluate();
	Passed(budget);
	Near(budget.previous_integral, 24.0L);
	Near(budget.current_integral, 48.0L);
	Near(budget.capacity_rate, 96.0L);
	Near(budget.source, 96.0L);
	// Opposite coefficient errors have zero total integral; a total-only gate
	// would miss them, while the independently integrated free rows must fail.
	fixture.current.front() += 0.1;
	fixture.current.back() -= 0.1;
	budget = fixture.Evaluate();
	Near(budget.volume_defect, 0.0L);
	Require(!iga::TransportBudgetWithin(budget.free_residual_l2, budget.free_component_scale, 1e-6, 1e-12),
		"equal-and-opposite field corruption was hidden by the global budget");
}

void CompressibleAdvection()
{
	Fixture fixture;
	fixture.system.terms = {{iga::TermKind::TimeDerivative, 0, 0, 1.0, ""},
		{iga::TermKind::Advection, 0, 0, 1.0, "prescribed"},
		{iga::TermKind::LinearCoupling, 0, 0, -1.0, ""},
		{iga::TermKind::VolumeSource, 0, 0, -1.0, ""}};
	fixture.system.stabilization = {{"c", "supg", "prescribed"}};
	fixture.element.boundary_labels.fill(0);
	iga::NamedBoundaryDefinition boundary;
	boundary.label = 0;
	boundary.name = "surface";
	fixture.configuration.boundaries.push_back(boundary);
	for (std::size_t node = 0; node < 64; ++node) {
		const auto x = fixture.element.bezier_points[node][0];
		fixture.previous[node] = fixture.current[node] = 1.0+2.0*x;
		fixture.velocity[node] = {x, 0.0, 0.0};
	}
	const auto budget = fixture.Evaluate();
	Passed(budget);
	Near(budget.advection_volume, 1.0L);
	Near(budget.advection_outward, 3.0L);
	Near(budget.velocity_divergence, 2.0L);
	Near(budget.reaction, -2.0L);
	Near(budget.source, -1.0L);
	Require(std::abs(budget.surface_defect+budget.velocity_divergence) > 1.0L,
		"analytic fixture would not detect omitted velocity-divergence correction");
}

void DiffusiveFlux()
{
	Fixture fixture;
	fixture.system.terms = {{iga::TermKind::TimeDerivative, 0, 0, 1.0, ""},
		{iga::TermKind::Diffusion, 0, 0, 0.5, ""}};
	const std::array<double, 6> flux{{-0.5, -0.5, 0.5, 0.5, -0.5, 0.5}};
	for (std::size_t face = 0; face < 6; ++face) {
		fixture.element.boundary_labels[face] = static_cast<std::int32_t>(face);
		iga::NamedBoundaryDefinition boundary;
		boundary.label = static_cast<int>(face);
		boundary.name = "face"+std::to_string(face);
		iga::FieldBoundaryCondition condition;
		condition.field = "c";
		condition.kind = iga::FieldBoundaryKind::Flux;
		condition.value = {flux[face]};
		boundary.conditions.push_back(condition);
		fixture.configuration.boundaries.push_back(boundary);
	}
	for (std::size_t node = 0; node < 64; ++node) {
		const auto& p = fixture.element.bezier_points[node];
		fixture.previous[node] = fixture.current[node] = p[0]+p[1]+p[2];
	}
	Passed(fixture.Evaluate());
	fixture.configuration.boundaries[2].conditions[0].value[0] *= -1;
	const auto bad = fixture.Evaluate();
	Require(!iga::TransportBudgetWithin(bad.surface_defect, bad.surface_scale, 1e-6, 1e-12), "wrong flux sign accepted");
	Require(!iga::TransportBudgetWithin(bad.free_residual_l2, bad.free_component_scale, 1e-6, 1e-12), "wrong flux residual accepted");
}

void EssentialAndRejection()
{
	Fixture fixture;
	fixture.boundaries.constrained[0] = 1;
	fixture.boundaries.value[0] = 2.0;
	fixture.previous[0] = fixture.current[0] = 2.0;
	Near(fixture.Evaluate().essential_error_l2, 0.0L);
	fixture.previous[0] = 1.0;
	Require(fixture.Evaluate().essential_error_l2 > 0.9L, "previous essential value not checked");
	fixture.current[0] = 3.0;
	Require(fixture.Evaluate().essential_error_l2 > 1.4L, "current essential value not checked");
	auto missing = fixture.Accumulator();
	Reject([&] { missing.Finish(); });
	auto duplicate = fixture.Accumulator();
	duplicate.Add(fixture.element);
	Reject([&] { duplicate.Add(fixture.element); });
	Reject([&] { duplicate.Finish(); });
	auto bad_element = fixture.element;
	bad_element.bezier_points.fill({0.0, 0.0, 0.0});
	auto poisoned = fixture.Accumulator();
	Reject([&] { poisoned.Add(bad_element); });
	Reject([&] { poisoned.Finish(); });
	fixture.system.stabilization = {{"c", "supg", "prescribed"}};
	fixture.system.terms.push_back({iga::TermKind::Diffusion, 0, 0, 1.0, ""});
	Reject([&] { fixture.Accumulator(); });
	fixture.system.stabilization.clear();
	fixture.current[0] = std::numeric_limits<double>::quiet_NaN();
	Reject([&] { fixture.Accumulator(); });
	Require(iga::TransportBudgetWithin(0.0L, 0.0L, 1e-6, 1e-12), "zero budget rejected");
	Require(!iga::TransportBudgetWithin(1e-11L, 0.0L, 1e-6, 1e-12), "zero-reference absolute gate ignored");
	Reject([] { iga::TransportBudgetWithin(0.0L, 1.0L, -1.0, 0.0); });
}

void CoupledSpecies()
{
	Fixture fixture;
	fixture.system.fields = {"a", "b"};
	fixture.system.field_index = {{"a", 0}, {"b", 1}};
	fixture.system.terms = {{iga::TermKind::TimeDerivative, 0, 0, 1.0, ""},
		{iga::TermKind::TimeDerivative, 1, 1, 1.0, ""},
		{iga::TermKind::LinearCoupling, 0, 0, 4.0, ""},
		{iga::TermKind::LinearCoupling, 1, 0, -4.0, ""}};
	fixture.previous.resize(128);
	fixture.current.resize(128);
	fixture.boundaries.constrained.assign(128, 0);
	fixture.boundaries.value.assign(128, 0.0);
	for (std::size_t node = 0; node < 64; ++node) {
		fixture.previous[2*node] = 2.0;
		fixture.previous[2*node+1] = 1.0;
		fixture.current[2*node] = 1.0;
		fixture.current[2*node+1] = 2.0;
	}
	auto accumulator = fixture.Accumulator();
	accumulator.Add(fixture.element);
	const auto budgets = accumulator.Finish();
	Passed(budgets[0]);
	Passed(budgets[1]);
	Near(budgets[0].capacity_rate, -4.0L);
	Near(budgets[1].capacity_rate, 4.0L);
	Near(budgets[0].reaction, 4.0L);
	Near(budgets[1].reaction, -4.0L);
	Near(budgets[0].reaction+budgets[1].reaction, 0.0L);
	fixture.system.terms.back().trial = 1;
	auto wrong = fixture.Accumulator();
	wrong.Add(fixture.element);
	const auto bad = wrong.Finish();
	Require(!iga::TransportBudgetWithin(bad[1].free_residual_l2, bad[1].free_component_scale, 1e-6, 1e-12),
		"wrong interspecies trial field accepted");
}

void RobinExchange()
{
	Fixture fixture;
	fixture.system.terms = {{iga::TermKind::TimeDerivative, 0, 0, 1.0, ""},
		{iga::TermKind::Diffusion, 0, 0, 0.5, ""}};
	for (std::size_t node = 0; node < 64; ++node)
		fixture.previous[node] = fixture.current[node] = fixture.element.bezier_points[node][0];
	for (int face = 0; face < 6; ++face) {
		fixture.element.boundary_labels[face] = face;
		iga::NamedBoundaryDefinition boundary;
		boundary.label = face;
		boundary.name = "face"+std::to_string(face);
		iga::FieldBoundaryCondition condition;
		condition.field = "c";
		condition.kind = face == 2 ? iga::FieldBoundaryKind::Robin : iga::FieldBoundaryKind::Flux;
		condition.value = {face == 4 ? -0.5 : 0.0};
		condition.coefficient = 0.5;
		condition.exterior_value = 2.0;
		boundary.conditions.push_back(condition);
		fixture.configuration.boundaries.push_back(boundary);
	}
	const auto budget = fixture.Evaluate();
	Passed(budget);
	Near(budget.natural_outward, 0.0L);
	Near(budget.natural_absolute_flux, 1.0L);
	fixture.configuration.boundaries[2].conditions[0].exterior_value = 0.0;
	const auto bad = fixture.Evaluate();
	Require(!iga::TransportBudgetWithin(bad.surface_defect, bad.surface_scale, 1e-6, 1e-12),
		"incorrect Robin exterior concentration accepted");
}

} // namespace

int main()
{
	try {
		SourceAndScaling();
		CompressibleAdvection();
		DiffusiveFlux();
		EssentialAndRejection();
		CoupledSpecies();
		RobinExchange();
		std::cout << "transport_budget checks=" << checks << " passed\n";
		return 0;
	} catch (const std::exception& error) {
		std::cerr << "transport_budget: " << error.what() << '\n';
		return 1;
	}
}
