#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

void Require(bool condition, const std::string& message)
{
	if (!condition) throw std::runtime_error(message);
}

std::string Read(const fs::path& path)
{
	std::ifstream input(path);
	Require(static_cast<bool>(input), "cannot read "+path.string());
	return std::string((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
}

std::vector<std::string> Split(const std::string& line)
{
	std::vector<std::string> fields;
	std::stringstream stream(line);
	std::string field;
	while (std::getline(stream, field, ',')) fields.push_back(field);
	return fields;
}

std::vector<std::map<std::string, std::string>> Csv(const fs::path& path)
{
	std::ifstream input(path);
	Require(static_cast<bool>(input), "cannot read CSV "+path.string());
	std::string line;
	Require(static_cast<bool>(std::getline(input, line)), "CSV lacks header");
	const auto header = Split(line);
	std::vector<std::map<std::string, std::string>> rows;
	while (std::getline(input, line)) {
		if (line.empty()) continue;
		const auto fields = Split(line);
		Require(fields.size() == header.size(), "unexpected CSV shape");
		std::map<std::string, std::string> row;
		for (std::size_t i = 0; i < header.size(); ++i) row.emplace(header[i], fields[i]);
		rows.push_back(std::move(row));
	}
	return rows;
}

double FixtureFlowRelativeTolerance(const fs::path& root)
{
	// Explicit execution performs one sweep, so the closure evidence retains
	// the fixture's configured flow gate instead of inferring convergence.
	constexpr double maximum_closure_gate = 1.0e-10;
	const std::string configuration = Read(root/"simulation_config.json");
	const std::string key = "\"flow_relative_tolerance\":";
	const auto position = configuration.find(key);
	Require(position != std::string::npos, "aneurysm fixture lacks flow relative tolerance");
	std::size_t used = 0;
	double tolerance = 0.0;
	try { tolerance = std::stod(configuration.substr(position+key.size()), &used); }
	catch (const std::exception&) { throw std::runtime_error("aneurysm fixture flow relative tolerance is invalid"); }
	Require(used > 0 && std::isfinite(tolerance) && tolerance > 0.0,
		"aneurysm fixture flow relative tolerance is invalid");
	Require(tolerance <= maximum_closure_gate,
		"aneurysm fixture flow relative tolerance weakens the Phase 6 closure gate");
	return tolerance;
}

double Number(const std::map<std::string, std::string>& row, const char* key)
{
	const auto found = row.find(key);
	Require(found != row.end(), std::string("CSV key is absent: ")+key);
	std::size_t used = 0;
	double value = 0.0;
	try { value = std::stod(found->second, &used); }
	catch (const std::exception&) { throw std::runtime_error(std::string("invalid CSV number: ")+key); }
	Require(used == found->second.size() && std::isfinite(value), std::string("nonfinite CSV number: ")+key);
	return value;
}

fs::path CopyCase(const fs::path& source, const std::string& suffix)
{
	const auto nonce = std::chrono::steady_clock::now().time_since_epoch().count();
	const fs::path root = fs::temp_directory_path()/("tubularflowiga-phase6-"+suffix+"-"+std::to_string(nonce));
	fs::create_directory(root);
	fs::copy_file(source/"simulation_config.json", root/"simulation_config.json");
	for (const char* directory : {"source", "sink", "immersed"})
		fs::copy(source/directory, root/directory, fs::copy_options::recursive);
	return root;
}

void WriteFixedConfiguration(const fs::path& root)
{
	std::string text = Read(root/"simulation_config.json");
	const std::string before = "\"kind\": \"explicit\", \"maximum_iterations\": 1";
	const std::string after = "\"kind\": \"fixed\", \"maximum_iterations\": 50";
	const auto position = text.find(before);
	Require(position != std::string::npos, "explicit aneurysm fixture execution is malformed");
	text.replace(position, before.size(), after);
	const std::string relaxation_before = "\"relaxation_factor\": 0.5";
	const auto relaxation_position = text.find(relaxation_before);
	Require(relaxation_position != std::string::npos,
		"explicit aneurysm fixture relaxation is malformed");
	// This isolated strong-fixed run is a linear, steady low-Re pressure map.
	// A near-full update removes artificial under-relaxation while retaining the
	// configured [0.05, .999] stability bounds and all numerical gates.
	text.replace(relaxation_position, relaxation_before.size(),
		"\"relaxation_factor\": 0.999");
	std::ofstream output(root/"simulation_config.json", std::ios::trunc);
	output << text;
	Require(static_cast<bool>(output), "cannot write fixed aneurysm fixture");
}

int Run(const fs::path& executable, const fs::path& root, const fs::path& output, bool inject)
{
	const fs::path log = root/(output.filename().string()+".log");
	const std::string command = std::string(inject ? "TUBULARFLOWIGA_INJECT_EXPLICIT_COUPLING_FAILURE_STEP=1 " : "")
		+"\""+executable.string()+"\" --graph-case \""+root.string()+"\" --output-dir \""
		+output.string()+"\" -immersed_static_ksp_type preonly -immersed_static_pc_type lu >\""
		+log.string()+"\" 2>&1";
	return std::system(command.c_str());
}

void ValidateManifest(const fs::path& output, const char* execution)
{
	const std::string manifest = Read(output/"graph_binding_manifest.json");
	for (const std::string& token : std::vector<std::string>{"\"execution\": \""+std::string(execution)+"\"",
		"\"kind\":\"three_d_immersed_flow\"", "\"surface_hash\":\"",
		"\"grid\":{", "\"catalog_audit\":{", "\"volume_points\":",
		"\"surface_points\":", "\"ghost_faces\":"})
		Require(manifest.find(token) != std::string::npos, "immersed manifest audit is incomplete");
}

struct RunMaxima {
	double edge_normalized_flow = 0.0;
	double strong_normalized_pressure = 0.0;
	double strong_iterations = 0.0;
	double immersed_normalized_balance = 0.0;
	double net_external_normalized_balance = 0.0;
};

RunMaxima ValidateRun(const fs::path& output, bool fixed, double flow_relative_tolerance)
{
	Require(fs::is_regular_file(output/"graph_binding_manifest.json"), "completion marker is missing");
	ValidateManifest(output, fixed ? "fixed" : "explicit");
	const auto edges = Csv(output/"pressure_flow_edges.csv");
	const auto balances = Csv(output/"pressure_flow_domain_balances.csv");
	const auto summary = Csv(output/"pressure_flow_steps.csv");
	Require(edges.size() == 6 && balances.size() == 3 && summary.size() == 3,
		"three quasi-static closure samples were not emitted");
	RunMaxima maxima;
	for (const auto& edge : edges) {
		const double flow = std::abs(Number(edge, "normalized_flow_residual"));
		maxima.edge_normalized_flow = std::max(maxima.edge_normalized_flow, flow);
		Require(flow <= flow_relative_tolerance,
			"edge normalized-flow closure exceeds the fixture flow-relative tolerance");
		if (fixed) {
			const double pressure = std::abs(Number(edge, "normalized_pressure_residual"));
			maxima.strong_normalized_pressure = std::max(maxima.strong_normalized_pressure, pressure);
			Require(pressure <= 1e-6,
				"strong fixed pressure residual exceeds 1e-6");
		}
	}
	for (const auto& row : balances) {
		const double balance = std::abs(Number(row, "normalized_mass_imbalance"));
		maxima.immersed_normalized_balance = std::max(maxima.immersed_normalized_balance, balance);
		Require(balance <= 1e-3,
			"immersed open-port normalized balance exceeds 1e-3");
	}
	for (const auto& row : summary) {
		const double iterations = Number(row, "iterations");
		if (fixed) maxima.strong_iterations = std::max(maxima.strong_iterations, iterations);
		Require(iterations <= 50.0, "strong fixed coupling exceeded 50 iterations");
		const double external = std::abs(Number(row, "external_outward_flow_m3_s"))/1e-4;
		maxima.net_external_normalized_balance = std::max(maxima.net_external_normalized_balance, external);
		Require(external <= 2e-3,
			"net external 1D balance exceeds 2e-3");
	}
	return maxima;
}

void RequireExactReplay(const fs::path& first, const fs::path& second)
{
	for (const char* name : {"pressure_flow_steps.csv", "pressure_flow_edges.csv",
		"pressure_flow_iterations.csv", "pressure_flow_ports.csv",
		"pressure_flow_domain_balances.csv", "one_d_initialization.csv"})
		Require(Read(first/name) == Read(second/name), std::string("retry did not exactly replay ")+name);
}

} // namespace

int main(int argc, char** argv)
{
	int status = 0;
	std::error_code ignored;
	fs::path explicit_root, fixed_root;
	try {
		Require(argc == 2, "usage: phase6_aneurysm_closure_test <iga_multidomain_flow>");
		const fs::path executable = fs::absolute(argv[1]);
		Require(fs::is_regular_file(executable), "multidomain runner is missing");
		const fs::path source = "../../examples/vascular_flow/immersed_aneurysm_chain";
		explicit_root = CopyCase(source, "explicit");
		const double flow_relative_tolerance = FixtureFlowRelativeTolerance(explicit_root);
		const fs::path explicit_output = explicit_root/"explicit";
		Require(Run(executable, explicit_root, explicit_output, false) == 0, "production explicit aneurysm run failed");
		const RunMaxima explicit_maxima = ValidateRun(explicit_output, false, flow_relative_tolerance);
		const fs::path failed_output = explicit_root/"injected_failure";
		Require(Run(executable, explicit_root, failed_output, true) != 0,
			"injected precommit failure unexpectedly succeeded");
		Require(!fs::exists(failed_output/"graph_binding_manifest.json")
			&& !fs::exists(failed_output/"pressure_flow_steps.csv"),
			"injected precommit failure published completion output");
		const fs::path retry_output = explicit_root/"retry";
		Require(Run(executable, explicit_root, retry_output, false) == 0, "explicit retry failed");
		ValidateRun(retry_output, false, flow_relative_tolerance);
		RequireExactReplay(explicit_output, retry_output);
		fixed_root = CopyCase(source, "fixed");
		WriteFixedConfiguration(fixed_root);
		const fs::path fixed_output = fixed_root/"fixed";
		Require(Run(executable, fixed_root, fixed_output, false) == 0, "production strong-fixed aneurysm run failed");
		const RunMaxima fixed_maxima = ValidateRun(fixed_output, true, flow_relative_tolerance);
		std::cout << "{\"phase6_coupled_closure\":{\"explicit_samples\":3,\"strong_fixed_samples\":3,"
			"\"edge_normalized_flow_gate\":" << flow_relative_tolerance << ",\"strong_pressure_gate\":1e-6,"
			"\"strong_iteration_gate\":50,\"immersed_normalized_balance_gate\":1e-3,"
			"\"net_external_normalized_balance_gate\":2e-3,"
			"\"max_edge_normalized_flow\":"
			<< std::max(explicit_maxima.edge_normalized_flow, fixed_maxima.edge_normalized_flow)
			<< ",\"max_strong_normalized_pressure\":" << fixed_maxima.strong_normalized_pressure
			<< ",\"max_strong_iterations\":" << fixed_maxima.strong_iterations
			<< ",\"max_immersed_3d_normalized_balance\":"
			<< std::max(explicit_maxima.immersed_normalized_balance, fixed_maxima.immersed_normalized_balance)
			<< ",\"max_net_external_1d_normalized_balance\":"
			<< std::max(explicit_maxima.net_external_normalized_balance, fixed_maxima.net_external_normalized_balance)
			<< ",\"explicit_manifest_audited\":true,\"strong_fixed_manifest_audited\":true,"
			"\"failure_no_completion_marker\":true,\"retry_exact_replay\":true,"
			"\"manifest_immersed_kind_hash_grid_catalog\":true}}\n";
	} catch (const std::exception& error) {
		std::cerr << "phase6 aneurysm coupled closure failure: " << error.what() << '\n';
		status = 1;
	}
	if (!explicit_root.empty()) fs::remove_all(explicit_root, ignored);
	if (!fixed_root.empty()) fs::remove_all(fixed_root, ignored);
	return status;
}
