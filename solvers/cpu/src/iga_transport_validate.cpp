#include "TransportBudget.hpp"
#include "CheckedText.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

namespace {

std::vector<double> ReadField(const std::filesystem::path& path, std::size_t nodes, std::size_t fields)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot read budget field: "+path.string());
	std::vector<double> values(nodes*fields);
	for (std::size_t node = 0; node < nodes; ++node) {
		std::string line;
		if (!std::getline(input, line)) throw std::runtime_error("budget field is truncated");
		std::istringstream row(line);
		std::string id_text;
		if (!(row >> id_text) || id_text.empty() || id_text.front() == '-')
			throw std::runtime_error("budget node IDs must be exactly 0..N-1");
		std::size_t used = 0;
		const auto id = std::stoull(id_text, &used);
		if (used != id_text.size() || id != node)
			throw std::runtime_error("budget node IDs must be exactly 0..N-1");
		for (std::size_t field = 0; field < fields; ++field)
			if (!(row >> values[node*fields+field]) || !std::isfinite(values[node*fields+field]))
				throw std::runtime_error("invalid or nonfinite budget field value");
		std::string extra;
		if (row >> extra) throw std::runtime_error("budget field has extra columns");
	}
	std::string extra;
	if (std::getline(input, extra)) throw std::runtime_error("budget field has extra rows");
	if (input.bad() || !input.eof()) throw std::runtime_error("cannot finish reading budget field");
	return values;
}

double Tolerance(const std::string& text)
{
	std::size_t used = 0;
	const auto value = std::stod(text, &used);
	if (used != text.size() || !std::isfinite(value) || value < 0.0)
		throw std::invalid_argument("budget tolerances must be finite and nonnegative");
	return value;
}

std::string JsonString(const std::string& value)
{
	std::ostringstream output;
	output.exceptions(std::ios::badbit | std::ios::failbit);
	output << '"';
	for (const unsigned char c : value) {
		if (c == '"' || c == '\\') output << '\\' << static_cast<char>(c);
		else if (c < 0x20) output << "\\u" << std::hex << std::setw(4) << std::setfill('0') << unsigned(c) << std::dec;
		else output << static_cast<char>(c);
	}
	output << '"';
	return output.str();
}

void Metric(std::ostream& output, const char* name, long double value)
{
	if (!std::isfinite(value) || std::abs(value) > std::numeric_limits<double>::max())
		throw std::runtime_error("nonfinite or unrepresentable budget metric");
	output << ",\n      " << JsonString(name) << ": " << static_cast<double>(value);
}

} // namespace

int main(int argc, char** argv)
{
	try {
		if (argc != 10 || std::string(argv[6]) != "--rtol" || std::string(argv[8]) != "--zero-atol")
			throw std::invalid_argument("usage: iga_transport_validate DATABASE CASE_DIR SYSTEM PREVIOUS.txt CURRENT.txt --rtol R --zero-atol A");
		const auto relative = Tolerance(argv[7]);
		const auto absolute = Tolerance(argv[9]);
		iga::Database database(argv[1]);
		if (database.header().nodes > std::numeric_limits<std::size_t>::max()
			|| database.header().elements > std::numeric_limits<std::size_t>::max())
			throw std::overflow_error("budget database dimensions exceed this platform");
		const auto nodes = static_cast<std::size_t>(database.header().nodes);
		const std::filesystem::path directory(argv[2]);
		const auto configuration = iga::ReadSimulationConfiguration((directory/"simulation_config.json").string());
		const auto system = iga::CompileLinearSystem(configuration, argv[3]);
		if (system.fields.empty() || nodes > std::numeric_limits<std::size_t>::max()/system.fields.size())
			throw std::overflow_error("budget field dimensions overflow");
		const auto labels = iga::ReadPointLabels((directory/"controlmesh.vtk").string(), nodes);
		const auto boundaries = iga::ResolveScalarBoundaries(configuration, system, labels);
		const auto previous = ReadField(argv[4], nodes, system.fields.size());
		const auto current = ReadField(argv[5], nodes, system.fields.size());
		const auto velocity = iga::ReadVelocity((directory/"initial_velocityfield.txt").string(), nodes);
		iga::TransportBudgetAccumulator accumulator(nodes, static_cast<std::size_t>(database.header().elements),
			system, configuration, boundaries, previous, current, velocity);
		for (std::uint64_t index = 0; index < database.header().elements; ++index) {
			const auto element = database.Load(index);
			if (element.id != index) throw std::runtime_error("budget database element ID differs from index");
			accumulator.Add(element);
		}
		const auto budgets = accumulator.Finish();
		std::ostringstream output;
		output.exceptions(std::ios::badbit | std::ios::failbit);
		output << std::setprecision(17) << "{\n  \"schema_version\": 1,\n  \"kind\": \"transport_species_budget\",\n"
			<< "  \"system\": " << JsonString(system.name) << ",\n  \"dt\": " << system.dt
			<< ",\n  \"relative_tolerance\": " << relative << ",\n  \"zero_reference_absolute_tolerance\": " << absolute
			<< ",\n  \"elements_integrated_once\": " << database.header().elements << ",\n  \"fields\": [";
		bool accepted = true;
		long double net_reaction = 0.0L;
		for (std::size_t field = 0; field < budgets.size(); ++field) {
			const auto& budget = budgets[field];
			const bool volume = iga::TransportBudgetWithin(budget.volume_defect, budget.volume_scale, relative, absolute);
			const bool surface = iga::TransportBudgetWithin(budget.surface_defect, budget.surface_scale, relative, absolute);
			const bool residual = iga::TransportBudgetWithin(budget.free_residual_l2, budget.free_component_scale, relative, absolute);
			const bool boundary = iga::TransportBudgetWithin(budget.essential_error_l2, budget.essential_reference_l2, relative, absolute);
			const bool passed = volume && surface && residual && boundary;
			accepted = accepted && passed;
			net_reaction += budget.reaction;
			output << (field ? "," : "") << "\n    {\n      \"name\": " << JsonString(system.fields[field])
				<< ",\n      \"passed\": " << (passed ? "true" : "false")
				<< ",\n      \"free_rows\": " << budget.free_rows
				<< ",\n      \"gates\": {\"volume\": " << (volume ? "true" : "false")
				<< ", \"surface\": " << (surface ? "true" : "false")
				<< ", \"free_residual\": " << (residual ? "true" : "false")
				<< ", \"essential_values\": " << (boundary ? "true" : "false") << '}';
			Metric(output, "previous_integral", budget.previous_integral);
			Metric(output, "current_integral", budget.current_integral);
			Metric(output, "capacity_rate", budget.capacity_rate);
			Metric(output, "current_capacity_rate_scale", budget.current_capacity_rate_scale);
			Metric(output, "previous_capacity_rate_scale", budget.previous_capacity_rate_scale);
			Metric(output, "advection_volume", budget.advection_volume);
			Metric(output, "advection_outward", budget.advection_outward);
			Metric(output, "advection_absolute_flux", budget.advection_absolute_flux);
			Metric(output, "concentration_velocity_divergence", budget.velocity_divergence);
			Metric(output, "reaction", budget.reaction);
			Metric(output, "volume_source", budget.source);
			Metric(output, "natural_outward", budget.natural_outward);
			Metric(output, "natural_absolute_flux", budget.natural_absolute_flux);
			Metric(output, "essential_inward", budget.essential_inward);
			Metric(output, "essential_absolute_exchange", budget.essential_absolute_exchange);
			Metric(output, "volume_defect", budget.volume_defect);
			Metric(output, "surface_defect", budget.surface_defect);
			Metric(output, "volume_scale", budget.volume_scale);
			Metric(output, "surface_scale", budget.surface_scale);
			Metric(output, "free_residual_l2", budget.free_residual_l2);
			Metric(output, "free_component_scale", budget.free_component_scale);
			Metric(output, "essential_error_l2", budget.essential_error_l2);
			Metric(output, "essential_reference_l2", budget.essential_reference_l2);
			output << "\n    }";
		}
		output << "\n  ]";
		Metric(output, "net_linear_reaction", net_reaction);
		output << ",\n  \"accepted\": " << (accepted ? "true" : "false")
			<< ",\n  \"scope\": \"Static prescribed-velocity backward Euler; weak essential-boundary exchange; scalar integrals are concentration times physical volume, not automatically mass in kg\"\n}\n";
		std::cout << output.str();
		iga::FlushCheckedText(std::cout);
		return accepted ? 0 : 2;
	} catch (const std::exception& error) {
		std::cerr << "iga_transport_validate: " << error.what() << '\n';
		return 1;
	}
}
