#define main BifurcationSmokeFixtureMain
#include "test_bifurcation_coupling_smoke.cpp"
#undef main

namespace {

std::string ReadFile(const fs::path& path)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot read multi-island configuration");
	return std::string((std::istreambuf_iterator<char>(input)),
		std::istreambuf_iterator<char>());
}

void WriteOneDDomain(std::ostream& output, const std::string& id,
	const std::string& case_directory, bool coupled_root, bool coupled_terminal,
	bool trailing_comma)
{
	output << "    {\"id\":\"" << id << "\",\"dimension\":\"1d\","
		"\"kind\":\"network_flow\",\"case\":\"" << case_directory << "\","
		"\"inlet_policy\":\"" << (coupled_root ? "coupled_root" : "configured_open_loop")
		<< "\",\"ports\":[{\"id\":\"root\",\"locator_kind\":\"runtime_port\","
		"\"locator\":\"root\",\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],"
		"\"requires\":[" << (coupled_root ? "\"flow_rate\"" : "") << "]},"
		"{\"id\":\"terminal\",\"locator_kind\":\"runtime_port\","
		"\"locator\":\"outlet:2\",\"provides\":[\"area\",\"flow_rate\","
		"\"mean_pressure\"],\"requires\":["
		<< (coupled_terminal ? "\"mean_pressure\"" : "") << "]}]}"
		<< (trailing_comma ? ",\n" : "\n");
}

void WriteThreeDDomain(std::ostream& output, const std::string& id,
	const std::string& case_directory, const std::string& database,
	bool trailing_comma)
{
	output << "    {\"id\":\"" << id << "\",\"dimension\":\"3d\","
		"\"kind\":\"body_fitted_iga_flow\",\"case\":\"" << case_directory
		<< "\",\"database\":\"" << database << "\",\"ports\":["
		"{\"id\":\"inlet\",\"locator_kind\":\"boundary_label\",\"locator\":\"1\","
		"\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],"
		"\"requires\":[\"flow_rate\"]},"
		"{\"id\":\"outlet_a\",\"locator_kind\":\"boundary_label\",\"locator\":\"2\","
		"\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],"
		"\"requires\":[\"mean_pressure\"]},"
		"{\"id\":\"outlet_b\",\"locator_kind\":\"boundary_label\",\"locator\":\"3\","
		"\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],"
		"\"requires\":[\"mean_pressure\"]},"
		"{\"id\":\"wall\",\"locator_kind\":\"boundary_label\",\"locator\":\"0\","
		"\"provides\":[\"flow_rate\"],\"requires\":[]}]}"
		<< (trailing_comma ? ",\n" : "\n");
}

void WriteEdge(std::ostream& output, const std::string& id,
	const std::string& first_domain, const std::string& first_port,
	const std::string& second_domain, const std::string& second_port,
	double initial_pressure, bool trailing_comma)
{
	output << "    {\"id\":\"" << id << "\",\"a\":{\"domain\":\""
		<< first_domain << "\",\"port\":\"" << first_port << "\"},"
		"\"b\":{\"domain\":\"" << second_domain << "\",\"port\":\""
		<< second_port << "\"},\"mode\":\"pressure_flow\","
		"\"initial_pressure_pa\":" << Number(initial_pressure) << "}"
		<< (trailing_comma ? ",\n" : "\n");
}

void WriteMultiIslandGraph(const fs::path& root, const std::string& execution,
	bool permuted = false)
{
	const bool explicit_mode = execution == "explicit";
	std::ofstream output(root/"simulation_config.json");
	output << "{\n  \"schema_version\":5,\"time\":{\"dt\":" << Number(kDt)
		<< ",\"steps\":" << kSteps << "},\"start_domain\":\"source\",\n"
		<< "  \"execution\":{\"kind\":\"" << execution << "\","
		"\"maximum_iterations\":" << (explicit_mode ? 1 : 100)
		<< ",\"pressure_relative_tolerance\":1e-6,\"pressure_reference_pa\":1,"
		"\"flow_relative_tolerance\":1e-10,\"relaxation_factor\":0.5,"
		"\"minimum_relaxation\":0.05,\"maximum_relaxation\":0.999},\n"
		<< "  \"domains\":[\n";
	if (!permuted) {
		WriteOneDDomain(output, "source", "source", false, true, true);
		WriteThreeDDomain(output, "island_a", "three_a", "a.ntiga", true);
		WriteOneDDomain(output, "bridge", "bridge", true, true, true);
		WriteOneDDomain(output, "leaf_a", "leaf_a", true, false, true);
		WriteThreeDDomain(output, "island_b", "three_b", "b.ntiga", true);
		WriteOneDDomain(output, "leaf_b", "leaf_b", true, false, true);
		WriteOneDDomain(output, "leaf_c", "leaf_c", true, false, false);
	} else {
		WriteOneDDomain(output, "leaf_c", "leaf_c", true, false, true);
		WriteOneDDomain(output, "leaf_b", "leaf_b", true, false, true);
		WriteThreeDDomain(output, "island_b", "three_b", "b.ntiga", true);
		WriteOneDDomain(output, "leaf_a", "leaf_a", true, false, true);
		WriteOneDDomain(output, "bridge", "bridge", true, true, true);
		WriteThreeDDomain(output, "island_a", "three_a", "a.ntiga", true);
		WriteOneDDomain(output, "source", "source", false, true, false);
	}
	output << "  ],\n  \"couplings\":[\n";
	if (!permuted) {
		WriteEdge(output, "e0_source_a", "source", "terminal", "island_a", "inlet", 0.0, true);
		WriteEdge(output, "e1_a_bridge", "island_a", "outlet_a", "bridge", "root", 0.02, true);
		WriteEdge(output, "e2_a_leaf", "island_a", "outlet_b", "leaf_a", "root", 0.01, true);
		WriteEdge(output, "e3_bridge_b", "bridge", "terminal", "island_b", "inlet", 0.0, true);
		WriteEdge(output, "e4_b_leaf", "island_b", "outlet_a", "leaf_b", "root", 0.015, true);
		WriteEdge(output, "e5_b_leaf", "island_b", "outlet_b", "leaf_c", "root", 0.005, false);
	} else {
		WriteEdge(output, "e5_b_leaf", "leaf_c", "root", "island_b", "outlet_b", 0.005, true);
		WriteEdge(output, "e4_b_leaf", "leaf_b", "root", "island_b", "outlet_a", 0.015, true);
		WriteEdge(output, "e3_bridge_b", "island_b", "inlet", "bridge", "terminal", 0.0, true);
		WriteEdge(output, "e2_a_leaf", "leaf_a", "root", "island_a", "outlet_b", 0.01, true);
		WriteEdge(output, "e1_a_bridge", "bridge", "root", "island_a", "outlet_a", 0.02, true);
		WriteEdge(output, "e0_source_a", "island_a", "inlet", "source", "terminal", 0.0, false);
	}
	output << "  ]\n}\n";
	if (!output) throw std::runtime_error("cannot write multi-island graph");
}

int RunMultidomain(const fs::path& root, const fs::path& output,
	const std::string& launcher = {}, const std::string& prefix = {})
{
	const auto log = root/(output.filename().string()+".log");
	const std::string linear_solver = launcher.empty()
		? " -ksp_type preonly -pc_type lu"
		: " -ksp_type gmres -pc_type bjacobi -sub_pc_type lu -ksp_rtol 1e-12";
	return std::system((prefix+launcher+"./iga_multidomain_flow --graph-case "+Quote(root)
		+" --output-dir "+Quote(output)+linear_solver+">"
		+Quote(log)+" 2>&1").c_str());
}

int RunBifurcation(const fs::path& root, const fs::path& output)
{
	const auto log = root/(output.filename().string()+".log");
	return std::system(("./iga_1d_3d_bifurcation --graph-case "+Quote(root)
		+" --output-dir "+Quote(output)+" -ksp_type preonly -pc_type lu >"
		+Quote(log)+" 2>&1").c_str());
}

void ValidateMultidomain(const fs::path& output, const std::string& execution)
{
	if (!fs::is_regular_file(output/"graph_binding_manifest.json"))
		throw std::runtime_error("multi-island completion marker is missing");
	const auto steps = ReadCsv(output/"pressure_flow_steps.csv");
	const auto edges = ReadCsv(output/"pressure_flow_edges.csv");
	const auto ports = ReadCsv(output/"pressure_flow_ports.csv");
	const auto balances = ReadCsv(output/"pressure_flow_domain_balances.csv");
	const auto initialization = ReadCsv(output/"one_d_initialization.csv");
	if (steps.size() != kSteps || edges.size() != 6*kSteps
		|| ports.size() != 18*kSteps || balances.size() != 2*kSteps
		|| initialization.size() != 5)
		throw std::runtime_error("multi-island long-form row count is invalid");
	std::map<std::string, int> edge_counts;
	for (const auto& edge : edges) {
		++edge_counts[edge.at("edge_id")];
		if (std::abs(Value(edge, "flow_residual_m3_s")) > 1.0e-13
			|| Value(edge, "normalized_flow_residual") > 1.0e-10)
			throw std::runtime_error("multi-island edge conservation failed");
	}
	if (edge_counts.size() != 6)
		throw std::runtime_error("multi-island edge coverage is incomplete");
	std::map<std::string, int> balance_counts;
	for (const auto& row : balances) {
		++balance_counts[row.at("domain_id")];
		if (Value(row, "normalized_mass_imbalance") > 1.0e-3)
			throw std::runtime_error("multi-island 3D mass balance failed");
	}
	if (balance_counts != std::map<std::string, int>{{"island_a", kSteps},
		{"island_b", kSteps}})
		throw std::runtime_error("multi-island domain-balance coverage is invalid");
	for (const auto& row : steps) {
		const double external = std::abs(Value(row, "external_outward_flow_m3_s"));
		double scale = 0.0;
		const int step = static_cast<int>(Value(row, "step"));
		for (const auto& port : ports)
			if (static_cast<int>(Value(port, "step")) == step
				&& ((port.at("domain_id") == "source" && port.at("port_id") == "root")
					|| (port.at("domain_id").find("leaf_") == 0
						&& port.at("port_id") == "terminal")))
				scale += std::abs(Value(port, "outward_flow_m3_s"));
		if (external/std::max(1.0e-14, scale) > 1.0e-10)
			throw std::runtime_error("multi-island global external conservation failed");
	}
	if (execution == "aitken") {
		const auto iterations = ReadCsv(output/"pressure_flow_iterations.csv");
		if (iterations.size() <= edges.size())
			throw std::runtime_error("multi-island Aitken run did not exercise a retry");
	}
}

} // namespace

int main()
{
	try {
		const auto root = fs::temp_directory_path()
			/("tubularflowiga_multidomain_"+std::to_string(static_cast<long long>(getpid())));
		for (const auto& directory : {"three_a", "three_b", "source", "bridge",
			"leaf_a", "leaf_b", "leaf_c"})
			fs::create_directories(root/directory);
		WriteThreeDCase(root/"three_a");
		WriteThreeDCase(root/"three_b");
		WriteOneDCase(root/"source", 1.0e-3);
		WriteOneDCase(root/"bridge", 5.5e-4);
		WriteOneDCase(root/"leaf_a", 4.5e-4);
		WriteOneDCase(root/"leaf_b", 3.0e-4);
		WriteOneDCase(root/"leaf_c", 2.5e-4);
		WriteDatabase(root/"a.ntiga", 1);
		WriteDatabase(root/"b.ntiga", 1);
		WriteMultiIslandGraph(root, "explicit");
		if (RunBifurcation(root, root/"bifurcation_rejected") == 0
			|| fs::exists(root/"bifurcation_rejected/graph_binding_manifest.json"))
			throw std::runtime_error("legacy bifurcation entry accepted a multi-island graph");
		if (RunMultidomain(root, root/"explicit_one") != 0)
			throw std::runtime_error("one-rank explicit multi-island run failed");
		ValidateMultidomain(root/"explicit_one", "explicit");
		WriteMultiIslandGraph(root, "explicit", true);
		if (RunMultidomain(root, root/"explicit_permuted") != 0)
			throw std::runtime_error("permuted explicit multi-island run failed");
		ValidateMultidomain(root/"explicit_permuted", "explicit");
		const auto canonical_edges = ReadCsv(root/"explicit_one/pressure_flow_edges.csv");
		const auto permuted_edges = ReadCsv(root/"explicit_permuted/pressure_flow_edges.csv");
		if (canonical_edges.size() != permuted_edges.size())
			throw std::runtime_error("permuted multi-island row count changed");
		for (std::size_t row = 0; row < canonical_edges.size(); ++row)
			if (canonical_edges[row].at("edge_id") != permuted_edges[row].at("edge_id")
				|| std::abs(Value(canonical_edges[row], "measured_pressure_pa")
					-Value(permuted_edges[row], "measured_pressure_pa")) > 1.0e-12)
				throw std::runtime_error("permuted multi-island execution changed");
		WriteMultiIslandGraph(root, "explicit");
		if (RunMultidomain(root, root/"failed", {},
			"TUBULARFLOWIGA_INJECT_BIFURCATION_FAILURE_STEP=1 ") == 0
			|| fs::exists(root/"failed/graph_binding_manifest.json"))
			throw std::runtime_error("failed multi-island run published a completion marker");
		WriteMultiIslandGraph(root, "aitken");
		if (RunMultidomain(root, root/"aitken_one") != 0)
			throw std::runtime_error("one-rank Aitken multi-island run failed");
		ValidateMultidomain(root/"aitken_one", "aitken");
		WriteDatabase(root/"a_two.ntiga", 2);
		WriteDatabase(root/"b_two.ntiga", 2);
		WriteMultiIslandGraph(root, "explicit");
		{
			std::string config = ReadFile(root/"simulation_config.json");
			const auto first = config.find("a.ntiga");
			const auto second = config.find("b.ntiga");
			if (first == std::string::npos || second == std::string::npos)
				throw std::runtime_error("multi-island database paths are missing");
			config.replace(first, 7, "a_two.ntiga");
			config.replace(config.find("b.ntiga"), 7, "b_two.ntiga");
			std::ofstream rewritten(root/"simulation_config.json", std::ios::trunc);
			rewritten << config;
		}
		if (RunMultidomain(root, root/"explicit_two", "mpiexec -np 2 ") != 0)
			throw std::runtime_error("two-rank explicit multi-island run failed");
		ValidateMultidomain(root/"explicit_two", "explicit");
		const auto one = ReadCsv(root/"explicit_one/pressure_flow_edges.csv");
		const auto two = ReadCsv(root/"explicit_two/pressure_flow_edges.csv");
		if (one.size() != two.size()) throw std::runtime_error("multi-island MPI row mismatch");
		for (std::size_t row = 0; row < one.size(); ++row) {
			if (one[row].at("edge_id") != two[row].at("edge_id"))
				throw std::runtime_error("multi-island MPI ordering mismatch");
			for (const auto& name : {"measured_pressure_pa", "first_outward_flow_m3_s",
				"second_outward_flow_m3_s"})
				if (std::abs(Value(one[row], name)-Value(two[row], name))
					> 1.0e-11*std::max({1.0, std::abs(Value(one[row], name)),
						std::abs(Value(two[row], name))}))
					throw std::runtime_error("multi-island MPI numerical mismatch");
		}
		fs::remove_all(root);
		std::cout << "multidomain flow smoke test passed\n";
		return 0;
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n';
		return 1;
	}
}
