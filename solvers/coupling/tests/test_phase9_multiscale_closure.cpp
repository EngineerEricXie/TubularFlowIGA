#include "CouplingFixture.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr double kDensity = 1.0e-6;
constexpr double kViscosity = 1.0;
constexpr double kFlow = 1.0e-3;
constexpr double kCapacitance = 5.0e-4;
constexpr double kSourceResistance = 20.0;
constexpr double kTerminalAProximalResistance = 10.0;
constexpr double kTerminalADistalResistance = 100.0;
constexpr double kTerminalACapacitance = 5.0e-4;
constexpr double kTerminalAInitialPressure = 0.0;
constexpr double kTerminalBProximalResistance = 20.0;
constexpr double kTerminalBDistalResistance = 80.0;
constexpr double kTerminalBCapacitance = 6.25e-4;
constexpr double kTerminalBInitialPressure = 0.0;
constexpr double kRadius = 0.56418958354775628; // pi*r^2 = 1 m^2.
constexpr int kStrongMaximumIterations = 50;
constexpr double kStrongRelaxation = 0.5;
constexpr double kPressureTolerance = 1.0e-6;
constexpr double kFlowTolerance = 1.0e-10;

using Row = std::map<std::string, std::string>;

void Require(bool condition, const std::string& message)
{
	if (!condition) throw std::runtime_error(message);
}

std::string Number(double value)
{
	std::ostringstream output;
	output << std::setprecision(17) << value;
	return output.str();
}

std::string Quote(const fs::path& path)
{
	return "'"+path.string()+"'";
}

std::string Read(const fs::path& path)
{
	std::ifstream input(path);
	Require(static_cast<bool>(input), "cannot read "+path.string());
	return std::string((std::istreambuf_iterator<char>(input)),
		std::istreambuf_iterator<char>());
}

std::vector<Row> ReadCsv(const fs::path& path)
{
	std::ifstream input(path);
	Require(static_cast<bool>(input), "cannot read CSV "+path.string());
	std::string line;
	Require(static_cast<bool>(std::getline(input, line)), "CSV lacks a header");
	std::vector<std::string> names;
	std::istringstream header(line);
	for (std::string name; std::getline(header, name, ',');) names.push_back(name);
	std::vector<Row> rows;
	while (std::getline(input, line)) {
		if (line.empty()) continue;
		std::istringstream fields(line);
		Row row;
		for (const auto& name : names) {
			std::string field;
			Require(static_cast<bool>(std::getline(fields, field, ',')), "short CSV row");
			row.emplace(name, field);
		}
		rows.push_back(std::move(row));
	}
	return rows;
}

double Value(const Row& row, const char* key)
{
	const auto found = row.find(key);
	Require(found != row.end(), std::string("CSV column is absent: ")+key);
	std::size_t used = 0;
	double value = 0.0;
	try { value = std::stod(found->second, &used); }
	catch (const std::exception&) { throw std::runtime_error(std::string("invalid CSV value: ")+key); }
	Require(used == found->second.size() && std::isfinite(value),
		std::string("nonfinite CSV value: ")+key);
	return value;
}

void WriteThreeDCase(const fs::path& directory, int transverse, int axial,
	double dt, int steps)
{
	iga::test::WriteC2SquareDuctThreeDCase(directory, transverse, dt, steps,
		kDensity, kViscosity, 1.0, axial);
}

void WriteTreeCase(const fs::path& directory, double dt, int steps)
{
	std::ofstream network(directory/"tree.swc");
	network << std::setprecision(17)
		<< "1 2 0 0 0 " << kRadius << " -1\n"
		<< "2 2 1 0.35 0 " << kRadius << " 1\n"
		<< "3 2 1 -0.35 0 " << kRadius << " 1\n";
	std::ofstream config(directory/"simulation_config.json");
	config << "{\"schema_version\":3,\"dimension\":\"1d\","
		"\"simulation_scope\":{\"mode\":\"flow_only\"},"
		"\"geometry\":{\"kind\":\"swc_network\",\"file\":\"tree.swc\",\"length_scale_to_m\":1},"
		"\"fields\":[{\"name\":\"area\",\"kind\":\"scalar\"},{\"name\":\"flow_rate\",\"kind\":\"scalar\"},{\"name\":\"pressure\",\"kind\":\"pressure\"}],"
		"\"time\":{\"dt\":" << Number(dt) << ",\"steps\":" << steps << ",\"output_every\":1},"
		"\"temporal_functions\":[{\"name\":\"seed\",\"kind\":\"constant\",\"units\":\"m3/s\",\"value\":" << Number(kFlow) << "}],"
		"\"equation_systems\":[{\"name\":\"flow\",\"kind\":\"network_flow_1d\",\"unknowns\":[\"area\",\"flow_rate\",\"pressure\"],\"model\":\"rigid\",\"scheme\":\"steady_poiseuille\",\"dynamic_viscosity\":" << Number(kViscosity) << ",\"density\":" << Number(kDensity) << ",\"discretization\":{\"cells_per_segment\":1}}],"
		"\"boundaries\":[{\"name\":\"inlet\",\"role\":\"inlet\",\"node_ids\":[1],\"conditions\":[{\"field\":\"flow_rate\",\"type\":\"dirichlet\",\"quantity\":\"flow_rate\",\"waveform\":\"seed\"}]},"
		"{\"name\":\"outlet_a\",\"role\":\"outlet\",\"node_ids\":[2],\"conditions\":[{\"field\":\"pressure\",\"type\":\"pressure\",\"value\":0}]},"
		"{\"name\":\"outlet_b\",\"role\":\"outlet\",\"node_ids\":[3],\"conditions\":[{\"field\":\"pressure\",\"type\":\"pressure\",\"value\":0}]}]}\n";
	Require(static_cast<bool>(network) && static_cast<bool>(config), "cannot write bifurcating 1D case");
}

void WriteZeroDModel(const fs::path& directory, const std::string& role)
{
	std::ofstream output(directory/"zero_d_model.json");
	if (role == "source")
		output << "{\"role\":\"source_reservoir\",\"capacitance_m3_pa\":" << Number(kCapacitance)
			<< ",\"resistance_pa_s_m3\":" << Number(kSourceResistance)
			<< ",\"initial_pressure_pa\":0,\"prescribed_flow_m3_s\":" << Number(kFlow) << "}";
	else {
		const bool terminal_a = role == "terminal_a";
		Require(terminal_a || role == "terminal_b", "unknown terminal RCR model");
		const double proximal = terminal_a ? kTerminalAProximalResistance : kTerminalBProximalResistance;
		const double distal = terminal_a ? kTerminalADistalResistance : kTerminalBDistalResistance;
		const double capacitance = terminal_a ? kTerminalACapacitance : kTerminalBCapacitance;
		const double initial_pressure = terminal_a ? kTerminalAInitialPressure : kTerminalBInitialPressure;
		output << "{\"role\":\"terminal_rcr\",\"proximal_resistance_pa_s_m3\":"
			<< Number(proximal) << ",\"distal_resistance_pa_s_m3\":" << Number(distal)
			<< ",\"capacitance_m3_pa\":" << Number(capacitance)
			<< ",\"distal_pressure_pa\":0,\"initial_pressure_pa\":" << Number(initial_pressure) << "}";
	}
	Require(static_cast<bool>(output), "cannot write 0D model");
}

void WriteGraph(const fs::path& root, const std::string& database, double dt, int steps)
{
	std::ofstream output(root/"simulation_config.json");
	output << "{\n\"schema_version\":5,\"time\":{\"dt\":" << Number(dt)
		<< ",\"steps\":" << steps << "},\"start_domain\":\"source\","
		"\"execution\":{\"kind\":\"fixed\",\"maximum_iterations\":" << kStrongMaximumIterations
		<< ",\"pressure_relative_tolerance\":" << Number(kPressureTolerance)
		<< ",\"pressure_reference_pa\":1,\"flow_relative_tolerance\":" << Number(kFlowTolerance)
		<< ",\"relaxation_factor\":" << Number(kStrongRelaxation)
		<< ",\"minimum_relaxation\":" << Number(kStrongRelaxation)
		<< ",\"maximum_relaxation\":" << Number(kStrongRelaxation) << "},\n\"domains\":["
		"{\"id\":\"source\",\"dimension\":\"0d\",\"kind\":\"zero_d_flow\",\"case\":\"source\",\"zero_d_model\":\"zero_d_model.json\",\"ports\":[{\"id\":\"port\",\"locator_kind\":\"zero_d_port\",\"locator\":\"port\",\"provides\":[\"mean_pressure\",\"flow_rate\"],\"requires\":[\"mean_pressure\"]}]},"
		"{\"id\":\"root3d\",\"dimension\":\"3d\",\"kind\":\"body_fitted_iga_flow\",\"case\":\"root3d\",\"database\":\"" << database << "\",\"ports\":["
		"{\"id\":\"inlet\",\"locator_kind\":\"boundary_label\",\"locator\":\"1\",\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],\"requires\":[\"flow_rate\"]},"
		"{\"id\":\"outlet\",\"locator_kind\":\"boundary_label\",\"locator\":\"2\",\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],\"requires\":[\"mean_pressure\"]},"
		"{\"id\":\"wall\",\"locator_kind\":\"boundary_label\",\"locator\":\"0\",\"provides\":[\"flow_rate\"],\"requires\":[]}]},"
		"{\"id\":\"tree\",\"dimension\":\"1d\",\"kind\":\"network_flow\",\"case\":\"tree\",\"inlet_policy\":\"coupled_root\",\"ports\":["
		"{\"id\":\"root\",\"locator_kind\":\"runtime_port\",\"locator\":\"root\",\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],\"requires\":[\"flow_rate\"]},"
		"{\"id\":\"terminal_a\",\"locator_kind\":\"runtime_port\",\"locator\":\"outlet:2\",\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],\"requires\":[\"mean_pressure\"]},"
		"{\"id\":\"terminal_b\",\"locator_kind\":\"runtime_port\",\"locator\":\"outlet:3\",\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],\"requires\":[\"mean_pressure\"]}]},";
	for (const char suffix : std::string("ab"))
		output << "{\"id\":\"rcr_" << suffix << "\",\"dimension\":\"0d\",\"kind\":\"zero_d_flow\",\"case\":\"rcr_" << suffix << "\",\"zero_d_model\":\"zero_d_model.json\",\"ports\":[{\"id\":\"port\",\"locator_kind\":\"zero_d_port\",\"locator\":\"port\",\"provides\":[\"mean_pressure\",\"flow_rate\"],\"requires\":[\"flow_rate\"]}]}" << (suffix == 'a' ? "," : "");
	output << "],\"couplings\":["
		"{\"id\":\"e_source_root\",\"a\":{\"domain\":\"source\",\"port\":\"port\"},\"b\":{\"domain\":\"root3d\",\"port\":\"inlet\"},\"mode\":\"pressure_flow\",\"initial_pressure_pa\":0},"
		"{\"id\":\"e_root_tree\",\"a\":{\"domain\":\"root3d\",\"port\":\"outlet\"},\"b\":{\"domain\":\"tree\",\"port\":\"root\"},\"mode\":\"pressure_flow\",\"initial_pressure_pa\":0},"
		"{\"id\":\"e_tree_a\",\"a\":{\"domain\":\"tree\",\"port\":\"terminal_a\"},\"b\":{\"domain\":\"rcr_a\",\"port\":\"port\"},\"mode\":\"pressure_flow\",\"initial_pressure_pa\":0},"
		"{\"id\":\"e_tree_b\",\"a\":{\"domain\":\"tree\",\"port\":\"terminal_b\"},\"b\":{\"domain\":\"rcr_b\",\"port\":\"port\"},\"mode\":\"pressure_flow\",\"initial_pressure_pa\":0}]}\n";
	Require(static_cast<bool>(output), "cannot write five-domain graph");
}

void PrepareCase(const fs::path& root, double dt, int steps, std::uint32_t ranks,
	const std::string& database)
{
	for (const char* directory : {"source", "root3d", "tree", "rcr_a", "rcr_b"})
		fs::create_directories(root/directory);
	WriteZeroDModel(root/"source", "source");
	WriteZeroDModel(root/"rcr_a", "terminal_a");
	WriteZeroDModel(root/"rcr_b", "terminal_b");
	// Two axial C2 body-fitted elements make two-rank ownership real while
	// preserving a separately auditable discrete duct resistance.
	iga::test::WriteC2SquareDuctDatabase(root/database, 1, ranks, 2);
	WriteThreeDCase(root/"root3d", 1, 2, dt, steps);
	WriteTreeCase(root/"tree", dt, steps); // one 1D macro step per 3D step.
	WriteGraph(root, database, dt, steps);
}

int Run(const fs::path& executable, const fs::path& root, const fs::path& output,
	const std::string& launcher = {}, bool inject = false)
{
	const fs::path log = root/(output.filename().string()+".log");
	const std::string solver = launcher.empty()
		? " -ksp_type preonly -pc_type lu"
		: " -ksp_type gmres -pc_type bjacobi -sub_pc_type lu -ksp_rtol 1e-12";
	return std::system((std::string(inject ? "TUBULARFLOWIGA_INJECT_BIFURCATION_FAILURE_STEP=1 " : "")
		+launcher+Quote(executable)+" --graph-case "+Quote(root)+" --output-dir "+Quote(output)
		+solver+" >"+Quote(log)+" 2>&1").c_str());
}

struct Metrics {
	double max_pressure_residual = 0.0;
	double max_flow_residual = 0.0;
	double max_zero_d_residual = 0.0;
	double max_three_d_balance = 0.0;
	int max_iterations = 0;
	double source_pressure = 0.0;
	double rcr_a_pressure = 0.0;
	double rcr_b_pressure = 0.0;
	double rcr_a_outward_flow = 0.0;
	double rcr_b_outward_flow = 0.0;
	double rcr_a_port_pressure = 0.0;
	double rcr_b_port_pressure = 0.0;
};

int Step(const Row& row)
{
	const double value = Value(row, "step");
	Require(value >= 1.0 && std::floor(value) == value && value <= 2147483647.0,
		"CSV step is not a positive integer");
	return static_cast<int>(value);
}

int Iteration(const Row& row)
{
	const double value = Value(row, "iteration");
	Require(value >= 1.0 && std::floor(value) == value && value <= 2147483647.0,
		"CSV iteration is not a positive integer");
	return static_cast<int>(value);
}

bool Close(double first, double second, double relative = 1.0e-12,
	double absolute = 1.0e-15)
{
	return std::abs(first-second) <= absolute+relative*std::max(std::abs(first), std::abs(second));
}

const std::set<std::string>& ExpectedEdges()
{
	static const std::set<std::string> edges{"e_source_root", "e_root_tree",
		"e_tree_a", "e_tree_b"};
	return edges;
}

const std::set<std::string>& ExpectedZeroDDomains()
{
	static const std::set<std::string> domains{"source", "rcr_a", "rcr_b"};
	return domains;
}

std::size_t JsonClosingBracket(const std::string& text, std::size_t open)
{
	Require(open < text.size() && (text[open] == '[' || text[open] == '{'), "invalid JSON bracket");
	const char close = text[open] == '[' ? ']' : '}';
	int depth = 0;
	bool quoted = false;
	for (std::size_t index = open; index < text.size(); ++index) {
		if (quoted) {
			if (text[index] == '\\') ++index;
			else if (text[index] == '"') quoted = false;
			continue;
		}
		if (text[index] == '"') quoted = true;
		else if (text[index] == text[open]) ++depth;
		else if (text[index] == close && --depth == 0) return index;
	}
	throw std::runtime_error("unterminated generated graph JSON");
}

std::vector<std::string> JsonObjects(const std::string& text, std::size_t begin, std::size_t end)
{
	std::vector<std::string> objects;
	for (std::size_t index = begin; index < end;) {
		if (text[index] != '{') { ++index; continue; }
		const std::size_t close = JsonClosingBracket(text, index);
		Require(close < end, "generated graph JSON object exceeds array");
		objects.push_back(text.substr(index, close-index+1));
		index = close+1;
	}
	return objects;
}

std::string JsonString(const std::string& object, const std::string& key)
{
	const std::string prefix = "\""+key+"\":\"";
	const std::size_t begin = object.find(prefix);
	Require(begin != std::string::npos, "generated graph JSON lacks string field '"+key+"'");
	const std::size_t value = begin+prefix.size();
	const std::size_t end = object.find('"', value);
	Require(end != std::string::npos, "generated graph JSON has unterminated string field '"+key+"'");
	return object.substr(value, end-value);
}

std::set<std::string> ExpectedPorts(const fs::path& graph_config)
{
	const std::string graph = Read(graph_config);
	const std::size_t domains_key = graph.find("\"domains\":[");
	Require(domains_key != std::string::npos, "generated graph config lacks domains");
	const std::size_t domains_open = graph.find('[', domains_key);
	const std::size_t domains_close = JsonClosingBracket(graph, domains_open);
	std::set<std::string> ports;
	for (const auto& domain : JsonObjects(graph, domains_open+1, domains_close)) {
		const std::string domain_id = JsonString(domain, "id");
		const std::size_t ports_key = domain.find("\"ports\":[");
		Require(ports_key != std::string::npos, "generated graph domain lacks ports");
		const std::size_t ports_open = domain.find('[', ports_key);
		const std::size_t ports_close = JsonClosingBracket(domain, ports_open);
		for (const auto& port : JsonObjects(domain, ports_open+1, ports_close))
			Require(ports.insert(domain_id+"."+JsonString(port, "id")).second,
				"generated graph has a duplicate logical port");
	}
	Require(ports.size() == 9, "generated graph does not define nine logical ports");
	return ports;
}

std::map<int, double> ExpectedClock(int steps, double configured_dt)
{
	Require(std::isfinite(configured_dt) && configured_dt > 0.0,
		"configured closure dt is invalid");
	std::map<int, double> clock;
	double start_time = 0.0;
	for (int step = 1; step <= steps; ++step) {
		const double end_time = start_time+configured_dt;
		Require(clock.emplace(step, end_time).second, "clock has a duplicate accepted step");
		start_time = end_time;
	}
	return clock;
}

void RequireTime(const Row& row, const std::map<int, double>& expected_times,
	const char* name)
{
	const int step = Step(row);
	const auto expected = expected_times.find(step);
	Require(expected != expected_times.end() && Close(Value(row, "time_s"), expected->second),
		std::string(name)+" time does not equal the repeated EndTime clock");
}

std::map<std::string, std::string> ZeroDModelIdentities(const fs::path& manifest)
{
	const std::string text = Read(manifest);
	const std::size_t domains_key = text.find("\"domains\"");
	Require(domains_key != std::string::npos, "completion marker lacks domains");
	const std::size_t domains_open = text.find('[', domains_key);
	Require(domains_open != std::string::npos, "completion marker has malformed domains");
	const std::size_t domains_close = JsonClosingBracket(text, domains_open);
	std::map<std::string, std::string> identities;
	for (const auto& domain : JsonObjects(text, domains_open+1, domains_close)) {
		if (domain.find("\"kind\":\"zero_d_flow\"") == std::string::npos) continue;
		const std::string id = JsonString(domain, "id");
		Require(identities.emplace(id, JsonString(domain, "model_identity_sha256")).second,
			"completion marker repeats a 0D model identity");
	}
	Require(identities.size() == ExpectedZeroDDomains().size(),
		"completion marker has incomplete 0D model identities");
	return identities;
}

void RequireDistinctTerminalBindings(const fs::path& output)
{
	const std::string graph = Read(output.parent_path()/"simulation_config.json");
	for (const auto& coupling : {
		std::string("{\"id\":\"e_tree_a\",\"a\":{\"domain\":\"tree\",\"port\":\"terminal_a\"},\"b\":{\"domain\":\"rcr_a\",\"port\":\"port\"}"),
		std::string("{\"id\":\"e_tree_b\",\"a\":{\"domain\":\"tree\",\"port\":\"terminal_b\"},\"b\":{\"domain\":\"rcr_b\",\"port\":\"port\"}")})
		Require(graph.find(coupling) != std::string::npos,
			"terminal edge is not bound to its declared model and port");
	const auto identities = ZeroDModelIdentities(output/"graph_binding_manifest.json");
	Require(identities.at("rcr_a") != identities.at("rcr_b"),
		"terminal models have indistinguishable runtime identities");
}

Metrics Validate(const fs::path& output, int steps, double configured_dt)
{
	Require(fs::is_regular_file(output/"graph_binding_manifest.json"), "completion marker is missing");
	RequireDistinctTerminalBindings(output);
	const auto expected_times = ExpectedClock(steps, configured_dt);
	const auto summary = ReadCsv(output/"pressure_flow_steps.csv");
	const auto edges = ReadCsv(output/"pressure_flow_edges.csv");
	const auto iterations = ReadCsv(output/"pressure_flow_iterations.csv");
	const auto ports = ReadCsv(output/"pressure_flow_ports.csv");
	const auto history = ReadCsv(output/"zero_d_flow_history.csv");
	const auto balances = ReadCsv(output/"pressure_flow_domain_balances.csv");
	const auto expected_ports = ExpectedPorts(output.parent_path()/"simulation_config.json");
	Metrics metrics;
	std::map<int, double> accepted_times;
	std::map<int, const Row*> summary_rows;
	for (const auto& row : summary) {
		const int step = Step(row);
		RequireTime(row, expected_times, "step-summary CSV");
		Require(summary_rows.emplace(step, &row).second, "step-summary CSV repeats an accepted step");
	}
	Require(static_cast<int>(summary_rows.size()) == steps, "step-summary CSV misses an accepted step");
	std::map<int, std::map<std::string, const Row*>> edge_rows;
	for (const auto& edge : edges) {
		const int step = Step(edge);
		Require(step <= steps, "edge CSV contains an unaccepted step");
		RequireTime(edge, expected_times, "edge CSV");
		const auto inserted_time = accepted_times.emplace(step, Value(edge, "time_s"));
		if (!inserted_time.second) Require(Close(inserted_time.first->second, Value(edge, "time_s")),
			"edge CSV time coverage is inconsistent");
		const auto inserted = edge_rows[step].emplace(edge.at("edge_id"), &edge);
		Require(inserted.second, "edge CSV repeats an edge in one accepted step");
		metrics.max_pressure_residual = std::max(metrics.max_pressure_residual,
			std::abs(Value(edge, "normalized_pressure_residual")));
		metrics.max_flow_residual = std::max(metrics.max_flow_residual,
			std::abs(Value(edge, "normalized_flow_residual")));
		// Independent cancellation: endpoint signs are outward from their own domains.
		Require(std::abs(Value(edge, "first_outward_flow_m3_s")
			+Value(edge, "second_outward_flow_m3_s")) <= 1.0e-12,
			"independent outward-flow cancellation failed");
		Require(std::abs(Value(edge, "measured_pressure_pa")-Value(edge, "applied_pressure_pa"))
			<= kPressureTolerance*std::max({1.0, std::abs(Value(edge, "measured_pressure_pa")),
				std::abs(Value(edge, "applied_pressure_pa"))}), "independent pressure continuity failed");
	}
	Require(static_cast<int>(edge_rows.size()) == steps, "edge CSV misses an accepted step");
	for (int step = 1; step <= steps; ++step) {
		const auto found = edge_rows.find(step);
		Require(found != edge_rows.end() && found->second.size() == ExpectedEdges().size()
			&& [&] { std::set<std::string> keys; for (const auto& item : found->second) keys.insert(item.first);
				return keys == ExpectedEdges(); }(), "edge CSV key coverage is incomplete");
	}

	std::map<int, std::set<int>> iteration_blocks;
	int previous_step = 0;
	int previous_iteration = 0;
	std::set<std::string> block_edges;
	auto finish_iteration_block = [&] {
		if (previous_step == 0) return;
		Require(block_edges == ExpectedEdges(), "iteration CSV has missing or duplicate edge diagnostics");
	};
	for (const auto& row : iterations) {
		const int step = Step(row);
		const int iteration = Iteration(row);
		Require(step <= steps, "iteration CSV contains an unaccepted step");
		RequireTime(row, expected_times, "iteration CSV");
		Require(std::abs(Value(row, "applied_relaxation")-kStrongRelaxation) <= 1.0e-12,
			"strong-fixed iteration relaxation differs from the configured value");
		if (step != previous_step || iteration != previous_iteration) {
			finish_iteration_block();
			Require(step > previous_step || (step == previous_step && iteration == previous_iteration+1),
				"iteration CSV blocks are not contiguous");
			previous_step = step;
			previous_iteration = iteration;
			block_edges.clear();
		}
		Require(block_edges.insert(row.at("edge_id")).second,
			"iteration CSV repeats an edge within one iteration");
		iteration_blocks[step].insert(iteration);
	}
	finish_iteration_block();
	Require(static_cast<int>(iteration_blocks.size()) == steps, "iteration CSV misses an accepted step");
	for (int step = 1; step <= steps; ++step) {
		const auto found = iteration_blocks.find(step);
		Require(found != iteration_blocks.end() && !found->second.empty()
			&& *found->second.begin() == 1 && *found->second.rbegin() == static_cast<int>(found->second.size()),
			"iteration CSV does not have consecutive iterations");
		metrics.max_iterations = std::max(metrics.max_iterations, *found->second.rbegin());
	}
	Require(metrics.max_iterations > 1, "strong-fixed closure never exercised more than one iteration");
	for (int step = 1; step <= steps; ++step) {
		const int final_iteration = *iteration_blocks.at(step).rbegin();
		int final_edges = 0;
		for (const auto& row : iterations) if (Step(row) == step && Iteration(row) == final_iteration) {
			++final_edges;
			Require(Value(row, "converged") == 1.0, "final strong-fixed iteration is not converged on every edge");
			Require(std::abs(Value(row, "normalized_pressure_residual")) <= kPressureTolerance
				&& std::abs(Value(row, "normalized_flow_residual")) <= kFlowTolerance,
				"final strong-fixed iteration misses residual gates");
		}
		Require(final_edges == static_cast<int>(ExpectedEdges().size()), "final iteration edge coverage is incomplete");
	}

	std::map<int, const Row*> balance_rows;
	for (const auto& row : balances) {
		const int step = Step(row);
		Require(step <= steps && row.at("domain_id") == "root3d",
			"3D balance CSV has an unexpected domain or step");
		RequireTime(row, expected_times, "3D balance CSV");
		Require(balance_rows.emplace(step, &row).second, "3D balance CSV repeats an accepted step");
		metrics.max_three_d_balance = std::max(metrics.max_three_d_balance,
			std::abs(Value(row, "normalized_mass_imbalance")));
		Require(std::abs(Value(row, "normalized_mass_imbalance")) <= 1.0e-8,
			"existing 3D mass gate failed");
	}
	Require(static_cast<int>(balance_rows.size()) == steps, "3D balance CSV misses an accepted step");

	std::map<int, std::map<std::string, const Row*>> history_rows;
	for (const auto& row : history) {
		const int step = Step(row);
		Require(step <= steps && ExpectedZeroDDomains().count(row.at("domain_id")) == 1,
			"0D history has an unexpected domain or step");
		RequireTime(row, expected_times, "0D history CSV");
		Require(history_rows[step].emplace(row.at("domain_id"), &row).second,
			"0D history repeats a domain in one accepted step");
		const double initial = Value(row, "initial_stored_volume_m3");
		const double final = Value(row, "final_stored_volume_m3");
		const double source = Value(row, "prescribed_source_amount_m3");
		const double sink = Value(row, "distal_sink_amount_m3");
		const double port = Value(row, "outward_graph_port_amount_m3");
		const double residual = Value(row, "residual_m3");
		const double recomputed = final-initial-source+sink+port;
		Require(std::abs(recomputed-residual) <= 128.0*std::numeric_limits<double>::epsilon()
			*std::max({1.0, std::abs(initial), std::abs(final), std::abs(source), std::abs(sink), std::abs(port)}),
			"0D BE accounting identity changed");
		metrics.max_zero_d_residual = std::max(metrics.max_zero_d_residual, std::abs(residual));
		Require(std::abs(residual) <= 1.0e-15, "0D BE storage residual is not roundoff scale");
		if (step == steps) {
			if (row.at("domain_id") == "source") metrics.source_pressure = Value(row, "stored_pressure_pa");
			if (row.at("domain_id") == "rcr_a") metrics.rcr_a_pressure = Value(row, "stored_pressure_pa");
			if (row.at("domain_id") == "rcr_b") metrics.rcr_b_pressure = Value(row, "stored_pressure_pa");
		}
	}
	Require(static_cast<int>(history_rows.size()) == steps, "0D history misses an accepted step");
	for (int step = 1; step <= steps; ++step)
		Require(history_rows.at(step).size() == ExpectedZeroDDomains().size(),
			"0D history key coverage is incomplete (both RCR states are required)");

	std::map<int, std::map<std::string, const Row*>> port_rows;
	for (const auto& row : ports) {
		const int step = Step(row);
		Require(step <= steps, "port CSV contains an unaccepted step");
		RequireTime(row, expected_times, "port CSV");
		const std::string key = row.at("domain_id")+"."+row.at("port_id");
		Require(expected_ports.count(key) == 1 && port_rows[step].emplace(key, &row).second,
			"port CSV has an unexpected or repeated logical port");
	}
	Require(static_cast<int>(port_rows.size()) == steps, "port CSV misses an accepted step");
	for (int step = 1; step <= steps; ++step) {
		Require(port_rows.at(step).size() == expected_ports.size(), "port CSV logical-port coverage is incomplete");
		const auto& root = port_rows.at(step);
		double root_sum = 0.0, root_absolute = 0.0;
		for (const char* port : {"root3d.inlet", "root3d.outlet", "root3d.wall"}) {
			const double flow = Value(*root.at(port), "outward_flow_m3_s");
			root_sum += flow;
			root_absolute += std::abs(flow);
		}
		const Row& balance = *balance_rows.at(step);
		Require(Close(root_sum, Value(balance, "boundary_flow_sum_m3_s"))
			&& Close(root_absolute, Value(balance, "boundary_flow_absolute_sum_m3_s")),
			"root3d boundary-balance identity disagrees with accepted ports");
		const double normalized = root_absolute == 0.0 ? 0.0 : 2.0*std::abs(root_sum)/root_absolute;
		Require(Close(normalized, Value(balance, "normalized_mass_imbalance")),
			"root3d normalized-balance identity changed");
		for (const auto& terminal : {
			std::make_pair("a", "rcr_a.port"), std::make_pair("b", "rcr_b.port")}) {
			const std::string edge_id = std::string("e_tree_")+terminal.first;
			const std::string tree_port = std::string("tree.terminal_")+terminal.first;
			const Row& edge = *edge_rows.at(step).at(edge_id);
			Require(Close(Value(*root.at(tree_port), "outward_flow_m3_s"),
				Value(edge, "first_outward_flow_m3_s"))
				&& Close(Value(*root.at(terminal.second), "outward_flow_m3_s"),
					Value(edge, "second_outward_flow_m3_s")),
				"terminal edge is not bound to its declared tree and RCR ports");
		}

		// Recompute each edge amount from the accepted endpoint flows.  The two
		// endpoints must cancel independently before the component balance uses it.
		const double previous_time = step == 1 ? 0.0 : expected_times.at(step-1);
		const double dt = expected_times.at(step)-previous_time;
		Require(Close(dt, configured_dt), "repeated EndTime clock changed its configured dt");
		double canceled_edge_amounts = 0.0;
		for (const auto& edge : edge_rows.at(step)) {
			const double amount = (Value(*edge.second, "first_outward_flow_m3_s")
				+Value(*edge.second, "second_outward_flow_m3_s"))*dt;
			Require(std::abs(amount) <= 1.0e-15, "accepted edge-port amount does not cancel");
			canceled_edge_amounts += amount;
		}
		const auto& zero = history_rows.at(step);
		Require(Close(Value(*zero.at("source"), "outward_graph_port_amount_m3"),
			Value(*edge_rows.at(step).at("e_source_root"), "first_outward_flow_m3_s")*dt)
			&& Close(Value(*zero.at("rcr_a"), "outward_graph_port_amount_m3"),
			Value(*edge_rows.at(step).at("e_tree_a"), "second_outward_flow_m3_s")*dt)
			&& Close(Value(*zero.at("rcr_b"), "outward_graph_port_amount_m3"),
			Value(*edge_rows.at(step).at("e_tree_b"), "second_outward_flow_m3_s")*dt),
			"0D graph-port amounts disagree with independently accepted edge flows");
		const double storage = Value(*zero.at("source"), "final_stored_volume_m3")
			-Value(*zero.at("source"), "initial_stored_volume_m3")
			+Value(*zero.at("rcr_a"), "final_stored_volume_m3")
			-Value(*zero.at("rcr_a"), "initial_stored_volume_m3")
			+Value(*zero.at("rcr_b"), "final_stored_volume_m3")
			-Value(*zero.at("rcr_b"), "initial_stored_volume_m3");
		const double pump = Value(*zero.at("source"), "prescribed_source_amount_m3");
		const double sinks = Value(*zero.at("rcr_a"), "distal_sink_amount_m3")
			+Value(*zero.at("rcr_b"), "distal_sink_amount_m3");
		Require(std::abs(storage-pump+sinks+canceled_edge_amounts) <= 3.0e-15,
			"per-step component-global volume balance failed");
		if (step == steps) {
			metrics.rcr_a_outward_flow = Value(*root.at("rcr_a.port"), "outward_flow_m3_s");
			metrics.rcr_b_outward_flow = Value(*root.at("rcr_b.port"), "outward_flow_m3_s");
			metrics.rcr_a_port_pressure = Value(*root.at("rcr_a.port"), "mean_pressure_pa");
			metrics.rcr_b_port_pressure = Value(*root.at("rcr_b.port"), "mean_pressure_pa");
		}
	}
	Require(metrics.max_pressure_residual <= kPressureTolerance
		&& metrics.max_flow_residual <= kFlowTolerance, "strong coupling residual gate failed");
	return metrics;
}

void RequireParity(const fs::path& one, const fs::path& two)
{
	auto require_identity = [](const std::vector<Row>& first, const std::vector<Row>& second,
		const std::vector<const char*>& keys, const std::string& name) {
		Require(first.size() == second.size(), "one/two-rank row count differs: "+name);
		for (std::size_t row = 0; row < first.size(); ++row)
			for (const auto* key : keys) {
				Require(first[row].at(key) == second[row].at(key),
					"one/two-rank identity differs: "+name+"."+key);
			}
	};
	auto require_close = [](double first, double second, const std::string& name) {
		Require(std::abs(first-second) <= 1.0e-10*std::max({1.0, std::abs(first), std::abs(second)}),
			"one/two-rank numerical parity differs: "+name);
	};
	auto require_numeric = [&](const std::vector<Row>& first, const std::vector<Row>& second,
		const std::vector<const char*>& keys, const std::string& name) {
		for (std::size_t row = 0; row < first.size(); ++row)
			for (const auto* key : keys)
				require_close(Value(first[row], key), Value(second[row], key), name+"."+key);
	};
	auto require_optional_numeric = [&](const std::vector<Row>& first, const std::vector<Row>& second,
		const std::vector<const char*>& keys, const std::string& name) {
		for (std::size_t row = 0; row < first.size(); ++row)
			for (const auto* key : keys) {
				const auto& left = first[row].at(key);
				const auto& right = second[row].at(key);
				Require(left.empty() == right.empty(), "one/two-rank optional-state presence differs: "+name+"."+key);
				if (!left.empty()) require_close(Value(first[row], key), Value(second[row], key), name+"."+key);
			}
	};
	const auto one_edges = ReadCsv(one/"pressure_flow_edges.csv");
	const auto two_edges = ReadCsv(two/"pressure_flow_edges.csv");
	require_identity(one_edges, two_edges, {"step", "edge_id"}, "edges");
	require_numeric(one_edges, two_edges, {"time_s", "applied_pressure_pa", "measured_pressure_pa",
		"pressure_residual_pa", "normalized_pressure_residual", "first_outward_flow_m3_s",
		"second_outward_flow_m3_s", "flow_residual_m3_s", "normalized_flow_residual"}, "edges");
	const auto one_ports = ReadCsv(one/"pressure_flow_ports.csv");
	const auto two_ports = ReadCsv(two/"pressure_flow_ports.csv");
	require_identity(one_ports, two_ports, {"step", "domain_id", "port_id"}, "ports");
	require_numeric(one_ports, two_ports, {"time_s", "outward_flow_m3_s"}, "ports");
	require_optional_numeric(one_ports, two_ports, {"area_m2", "mean_pressure_pa"}, "ports");
	const auto one_zero = ReadCsv(one/"zero_d_flow_history.csv");
	const auto two_zero = ReadCsv(two/"zero_d_flow_history.csv");
	require_identity(one_zero, two_zero, {"step", "domain_id"}, "zero_d_history");
	require_numeric(one_zero, two_zero, {"time_s", "stored_pressure_pa", "initial_stored_volume_m3",
		"final_stored_volume_m3", "prescribed_source_amount_m3", "distal_sink_amount_m3",
		"outward_graph_port_amount_m3", "residual_m3"}, "zero_d_history");
	const auto one_balances = ReadCsv(one/"pressure_flow_domain_balances.csv");
	const auto two_balances = ReadCsv(two/"pressure_flow_domain_balances.csv");
	require_identity(one_balances, two_balances, {"step", "domain_id"}, "domain_balances");
	require_numeric(one_balances, two_balances, {"time_s", "boundary_flow_sum_m3_s",
		"boundary_flow_absolute_sum_m3_s", "normalized_mass_imbalance"}, "domain_balances");
	const auto one_iterations = ReadCsv(one/"pressure_flow_iterations.csv");
	const auto two_iterations = ReadCsv(two/"pressure_flow_iterations.csv");
	require_identity(one_iterations, two_iterations, {"step", "iteration", "edge_id", "converged"}, "iterations");
	require_numeric(one_iterations, two_iterations, {"time_s", "applied_relaxation", "applied_pressure_pa",
		"measured_pressure_pa", "pressure_residual_pa", "normalized_pressure_residual",
		"first_outward_flow_m3_s", "second_outward_flow_m3_s", "flow_residual_m3_s",
		"normalized_flow_residual"}, "iterations");
}

void ValidateAnalyticLimit(const Metrics& metrics)
{
	const double branch_length = std::sqrt(1.0+0.35*0.35);
	const double branch_resistance = 8.0*kViscosity*branch_length
		/(iga::test::kPi*std::pow(kRadius, 4));
	const auto duct_profile = iga::test::BuildC2SquareDuctProfile(1);
	const double duct_resistance = kViscosity
		*iga::test::AuditC2SquareDuctRitzProfile(duct_profile).ritz_resistance;
	const double path_a = branch_resistance+kTerminalAProximalResistance+kTerminalADistalResistance;
	const double path_b = branch_resistance+kTerminalBProximalResistance+kTerminalBDistalResistance;
	const double parallel = path_a*path_b/(path_a+path_b);
	const double expected_source = kFlow*(kSourceResistance+duct_resistance+parallel);
	const double expected_a_flow = kFlow*path_b/(path_a+path_b);
	const double expected_b_flow = kFlow*path_a/(path_a+path_b);
	const double expected_a_capacitor = expected_a_flow*kTerminalADistalResistance;
	const double expected_b_capacitor = expected_b_flow*kTerminalBDistalResistance;
	const double expected_a_port = expected_a_flow
		*(kTerminalAProximalResistance+kTerminalADistalResistance);
	const double expected_b_port = expected_b_flow
		*(kTerminalBProximalResistance+kTerminalBDistalResistance);
	auto require_settled = [](double actual, double expected, const std::string& name) {
		Require(std::abs(actual-expected) <= 0.08*std::max(1.0e-12, std::abs(expected)),
			"settled "+name+" misses its independently assigned RCR value");
	};
	// The actual C2 body-fitted duct has one transverse element and two axial
	// elements. Its pressure drop is compared to the no-slip square-duct series resistance;
	// The slowest terminal has eight time constants and a separate C2 Ritz resistance audit
	// justify this 8% transient/discretization gate; it is not fitted to a run.
	Require(std::abs(metrics.source_pressure-expected_source)/expected_source <= 0.08,
		"settled source pressure misses the duct/tree/RCR analytic resistance network");
	require_settled(metrics.rcr_a_pressure, expected_a_capacitor, "rcr_a capacitor pressure");
	require_settled(metrics.rcr_b_pressure, expected_b_capacitor, "rcr_b capacitor pressure");
	require_settled(-metrics.rcr_a_outward_flow, expected_a_flow, "rcr_a port flow");
	require_settled(-metrics.rcr_b_outward_flow, expected_b_flow, "rcr_b port flow");
	require_settled(metrics.rcr_a_port_pressure, expected_a_port, "rcr_a port pressure");
	require_settled(metrics.rcr_b_port_pressure, expected_b_port, "rcr_b port pressure");
	Require(metrics.rcr_a_pressure > metrics.rcr_b_pressure+4.0e-3
		&& metrics.rcr_a_outward_flow > metrics.rcr_b_outward_flow+2.0e-5
		&& metrics.rcr_a_port_pressure > metrics.rcr_b_port_pressure+5.0e-4,
		"terminal RCR states, pressures, and flows do not retain their expected asymmetric ordering: "
		"a=("+Number(metrics.rcr_a_pressure)+","+Number(metrics.rcr_a_port_pressure)+","+Number(metrics.rcr_a_outward_flow)
		+"), b=("+Number(metrics.rcr_b_pressure)+","+Number(metrics.rcr_b_port_pressure)+","+Number(metrics.rcr_b_outward_flow)+")");
}

struct TemporalRatios {
	double source = 0.0;
	double rcr_a = 0.0;
	double rcr_b = 0.0;
};

TemporalRatios ValidateTemporal(const Metrics& coarse, const Metrics& fine,
	const Metrics& reference)
{
	const double coarse_source = std::abs(coarse.source_pressure-reference.source_pressure);
	const double fine_source = std::abs(fine.source_pressure-reference.source_pressure);
	const double coarse_rcr_a = std::abs(coarse.rcr_a_pressure-reference.rcr_a_pressure);
	const double fine_rcr_a = std::abs(fine.rcr_a_pressure-reference.rcr_a_pressure);
	const double coarse_rcr_b = std::abs(coarse.rcr_b_pressure-reference.rcr_b_pressure);
	const double fine_rcr_b = std::abs(fine.rcr_b_pressure-reference.rcr_b_pressure);
	Require(fine_source < coarse_source && fine_rcr_a < coarse_rcr_a && fine_rcr_b < coarse_rcr_b,
		"dt halving did not converge source/RCR pressures toward reference");
	const TemporalRatios ratios{coarse_source/fine_source, coarse_rcr_a/fine_rcr_a,
		coarse_rcr_b/fine_rcr_b};
	// With the h/4 run as reference, a first-order error has ratio
	// (4h-h)/(2h-h)=3 rather than 2.
	Require(ratios.source > 2.2 && ratios.source < 4.0
		&& ratios.rcr_a > 2.2 && ratios.rcr_a < 4.0
		&& ratios.rcr_b > 2.2 && ratios.rcr_b < 4.0,
		"dt-halving source/RCR pressure rate is not first order");
	return ratios;
}

} // namespace

int main(int argc, char** argv)
{
	int status = 0;
	fs::path root;
	std::error_code ignored;
	try {
		Require(argc == 2, "usage: phase9_multiscale_closure_test <iga_multidomain_flow>");
		const fs::path executable = fs::absolute(argv[1]);
		Require(fs::is_regular_file(executable), "production multidomain runner is missing");
		root = fs::temp_directory_path()/"tubularflowiga-phase9-multiscale";
		fs::remove_all(root);
		fs::create_directories(root);
		if (std::getenv("TUBULARFLOWIGA_PHASE9_TEMPORAL_ONLY")) {
			PrepareCase(root, 0.04, 1, 1, "temporal_coarse.ntiga");
			if (Run(executable, root, root/"temporal_coarse") != 0)
				throw std::runtime_error("coarse temporal closure failed");
			const Metrics coarse = Validate(root/"temporal_coarse", 1, 0.04);
			PrepareCase(root, 0.02, 2, 1, "temporal_fine.ntiga");
			if (Run(executable, root, root/"temporal_fine") != 0)
				throw std::runtime_error("fine temporal closure failed");
			const Metrics fine = Validate(root/"temporal_fine", 2, 0.02);
			PrepareCase(root, 0.01, 4, 1, "temporal_reference.ntiga");
			if (Run(executable, root, root/"temporal_reference") != 0)
				throw std::runtime_error("reference temporal closure failed");
			const Metrics reference = Validate(root/"temporal_reference", 4, 0.01);
			ValidateTemporal(coarse, fine, reference);
			std::cout << "phase9 temporal dt-halving closure passed\n";
			fs::remove_all(root, ignored);
			return 0;
		}
		if (std::getenv("TUBULARFLOWIGA_PHASE9_PARITY_ONLY")) {
			PrepareCase(root, 0.05, 8, 1, "root_1.ntiga");
			if (Run(executable, root, root/"serial") != 0)
				throw std::runtime_error("serial parity closure failed");
			Validate(root/"serial", 8, 0.05);
			PrepareCase(root, 0.05, 8, 2, "root_2.ntiga");
			if (Run(executable, root, root/"mpi2", "mpiexec -np 2 ") != 0)
				throw std::runtime_error("two-rank parity closure failed");
			Validate(root/"mpi2", 8, 0.05);
			RequireParity(root/"serial", root/"mpi2");
			std::cout << "phase9 serial/two-rank parity closure passed\n";
			fs::remove_all(root, ignored);
			return 0;
		}
		if (std::getenv("TUBULARFLOWIGA_PHASE9_RETRY_ONLY")) {
			PrepareCase(root, 0.05, 1, 1, "root_1.ntiga");
			if (Run(executable, root, root/"failed", {}, true) == 0
				|| fs::exists(root/"failed/graph_binding_manifest.json"))
				throw std::runtime_error("precommit failure published a completion marker");
			if (Run(executable, root, root/"retry") != 0
				|| Run(executable, root, root/"replay") != 0)
				throw std::runtime_error("clean retry failed: "+Read(root/"retry.log"));
			Require(Read(root/"retry/pressure_flow_edges.csv")
				== Read(root/"replay/pressure_flow_edges.csv")
				&& Read(root/"retry/zero_d_flow_history.csv")
				== Read(root/"replay/zero_d_flow_history.csv"),
				"clean retry did not replay deterministically");
			std::cout << "phase9 precommit failure/retry closure passed\n";
			fs::remove_all(root, ignored);
			return 0;
		}

		// The initial production gate is the smallest valid serial case.
		PrepareCase(root, 0.05, 8, 1, "root_1.ntiga");
		if (Run(executable, root, root/"serial") != 0)
			throw std::runtime_error("serial production closure failed: "+Read(root/"serial.log"));
		const Metrics serial = Validate(root/"serial", 8, 0.05);
		ValidateAnalyticLimit(serial);

		// Failure is deliberately before commit; it must not publish a final marker.
		if (Run(executable, root, root/"failed", {}, true) == 0
			|| fs::exists(root/"failed/graph_binding_manifest.json"))
			throw std::runtime_error("precommit failure published a completion marker");
		if (Run(executable, root, root/"retry") != 0) throw std::runtime_error("clean deterministic retry failed");
		Require(Read(root/"serial/pressure_flow_edges.csv") == Read(root/"retry/pressure_flow_edges.csv")
			&& Read(root/"serial/zero_d_flow_history.csv") == Read(root/"retry/zero_d_flow_history.csv"),
			"clean retry did not replay the accepted closure");

		PrepareCase(root, 0.05, 8, 2, "root_2.ntiga");
		if (Run(executable, root, root/"mpi2", "mpiexec -np 2 ") != 0)
			throw std::runtime_error("two-rank production closure failed");
		Validate(root/"mpi2", 8, 0.05);
		RequireParity(root/"serial", root/"mpi2");

		PrepareCase(root, 0.04, 1, 1, "temporal_coarse.ntiga");
		if (Run(executable, root, root/"temporal_coarse") != 0) throw std::runtime_error("coarse temporal closure failed");
		const Metrics coarse = Validate(root/"temporal_coarse", 1, 0.04);
		PrepareCase(root, 0.02, 2, 1, "temporal_fine.ntiga");
		if (Run(executable, root, root/"temporal_fine") != 0) throw std::runtime_error("fine temporal closure failed");
		const Metrics fine = Validate(root/"temporal_fine", 2, 0.02);
		PrepareCase(root, 0.01, 4, 1, "temporal_reference.ntiga");
		if (Run(executable, root, root/"temporal_reference") != 0) throw std::runtime_error("reference temporal closure failed");
		const Metrics reference = Validate(root/"temporal_reference", 4, 0.01);
		const TemporalRatios temporal_ratios = ValidateTemporal(coarse, fine, reference);

		std::cout << "{\"phase9_multiscale_closure\":{\"domains\":5,\"edges\":4,"
			"\"execution\":\"strong_fixed\",\"one_d_macro_substeps\":1,"
			"\"max_pressure_residual\":" << serial.max_pressure_residual
			<< ",\"max_flow_residual\":" << serial.max_flow_residual
			<< ",\"max_0d_residual_m3\":" << serial.max_zero_d_residual
			<< ",\"max_3d_normalized_balance\":" << serial.max_three_d_balance
			<< ",\"max_iterations\":" << serial.max_iterations
			<< ",\"rcr_a\":{\"stored_pressure_pa\":" << serial.rcr_a_pressure
			<< ",\"port_pressure_pa\":" << serial.rcr_a_port_pressure
			<< ",\"outward_flow_m3_s\":" << serial.rcr_a_outward_flow << "}"
			<< ",\"rcr_b\":{\"stored_pressure_pa\":" << serial.rcr_b_pressure
			<< ",\"port_pressure_pa\":" << serial.rcr_b_port_pressure
			<< ",\"outward_flow_m3_s\":" << serial.rcr_b_outward_flow << "}"
			<< ",\"temporal_error_ratios\":{\"source\":" << temporal_ratios.source
			<< ",\"rcr_a\":" << temporal_ratios.rcr_a
			<< ",\"rcr_b\":" << temporal_ratios.rcr_b << "}"
			<< ",\"mpi_ranks\":[1,2],\"precommit_no_marker\":true,\"deterministic_retry\":true}}\n";
	} catch (const std::exception& error) {
		std::cerr << "phase9 multiscale closure failure: " << error.what() << '\n';
		status = 1;
	}
	if (!root.empty()) fs::remove_all(root, ignored);
	return status;
}
