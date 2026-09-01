#include "IgaDatabase.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr double kDensity = 1.0e-6;
constexpr double kViscosity = 1.0;
constexpr double kDt = 1.0e-2;
constexpr int kSteps = 2;
constexpr double kRadius = 0.56418958354775628;

std::string Number(double value)
{
	std::ostringstream output;
	output << std::setprecision(17) << value;
	return output.str();
}

void WriteDatabase(const fs::path& path, std::uint32_t ranks)
{
	std::ofstream output(path, std::ios::binary | std::ios::trunc);
	if (!output) throw std::runtime_error("cannot create bifurcation database");
	constexpr std::uint64_t header_size = 88;
	const std::uint64_t element_offset = header_size+2*sizeof(std::uint64_t)+sizeof(std::int32_t);
	output.write(iga::kMagic.data(), iga::kMagic.size());
	iga::Write(output, iga::kVersion);
	iga::Write(output, ranks);
	iga::Write(output, std::uint64_t{1});
	iga::Write(output, std::uint64_t{64});
	iga::Write(output, iga::kBezierPointCount);
	iga::Write(output, std::uint32_t{0});
	const auto rank_index_position = output.tellp();
	iga::Write(output, std::uint64_t{0});
	for (int axis = 0; axis < 3; ++axis) iga::Write(output, 0.0);
	iga::Write(output, 1.0);
	iga::Write(output, 1.0);
	iga::Write(output, element_offset);
	iga::Write(output, std::uint64_t{0});
	iga::Write(output, std::int32_t{0});
	iga::Write(output, std::uint64_t{0});
	iga::Write(output, std::int32_t{0});
	iga::Write(output, std::int32_t{0});
	iga::Write(output, std::uint32_t{64});
	const std::array<std::int32_t, 6> labels{{0, 0, 2, 3, 1, 0}};
	output.write(reinterpret_cast<const char*>(labels.data()),
		static_cast<std::streamsize>(sizeof(labels)));
	for (std::int32_t node = 0; node < 64; ++node) iga::Write(output, node);
	for (std::uint8_t row = 0; row < 64; ++row) {
		iga::Write(output, std::uint8_t{1});
		iga::Write(output, row);
		iga::Write(output, 1.0);
	}
	for (int k = 0; k < 4; ++k)
		for (int j = 0; j < 4; ++j)
			for (int i = 0; i < 4; ++i) {
				const std::array<double, 3> point{{i/3.0, j/3.0, k/3.0}};
				output.write(reinterpret_cast<const char*>(point.data()),
					static_cast<std::streamsize>(sizeof(point)));
			}
	const auto rank_index_offset = static_cast<std::uint64_t>(output.tellp());
	for (std::uint32_t rank = 0; rank <= ranks; ++rank)
		iga::Write(output, static_cast<std::uint64_t>(rank));
	for (std::uint32_t rank = 0; rank < ranks; ++rank) iga::Write(output, std::uint64_t{0});
	output.seekp(header_size+sizeof(std::uint64_t));
	iga::Write(output, rank_index_offset);
	output.seekp(rank_index_position);
	iga::Write(output, rank_index_offset);
	if (!output) throw std::runtime_error("cannot finalize bifurcation database");
}

void WriteThreeDCase(const fs::path& directory)
{
	std::ofstream mesh(directory/"controlmesh.vtk");
	mesh << "# vtk DataFile Version 3.0\nbifurcation smoke\nASCII\n"
		<< "DATASET UNSTRUCTURED_GRID\nPOINTS 64 double\n";
	for (int k = 0; k < 4; ++k)
		for (int j = 0; j < 4; ++j)
			for (int i = 0; i < 4; ++i)
				mesh << i/3.0 << ' ' << j/3.0 << ' ' << k/3.0 << '\n';
	mesh << "CELLS 1 9\n8 0 3 15 12 48 51 63 60\nCELL_TYPES 1\n12\n"
		<< "POINT_DATA 64\nSCALARS boundary_label int 1\nLOOKUP_TABLE default\n";
	for (int k = 0; k < 4; ++k)
		for (int j = 0; j < 4; ++j)
			for (int i = 0; i < 4; ++i) {
				int label = -1;
				if (j == 0 || k == 0 || k == 3) label = 0;
				else if (i == 0) label = 1;
				else if (i == 3) label = 2;
				else if (j == 3) label = 3;
				mesh << label << '\n';
			}
	std::ofstream velocity(directory/"initial_velocityfield.txt");
	for (int k = 0; k < 4; ++k)
		for (int j = 0; j < 4; ++j)
			for (int i = 0; i < 4; ++i) {
				const double y = j/3.0;
				const double z = k/3.0;
				velocity << std::setprecision(17)
					<< (i == 0 ? 16.0*y*(1.0-y)*z*(1.0-z) : 0.0) << " 0 0\n";
			}
	std::ofstream config(directory/"simulation_config.json");
	config << "{\n  \"schema_version\":3,\"dimension\":\"3d\",\n"
		<< "  \"simulation_scope\":{\"mode\":\"flow_only\"},\n"
		<< "  \"coupling\":{\"scheme\":\"explicit_staggered\",\"flow_epsilon_m3_s\":1e-14,"
			"\"three_d_ports\":{\"inlet_label\":1,\"outlet_labels\":[2,3]}},\n"
		<< "  \"fields\":[{\"name\":\"velocity\",\"kind\":\"vector3\"},"
			"{\"name\":\"pressure\",\"kind\":\"pressure\"}],\n"
		<< "  \"time\":{\"dt\":" << Number(kDt) << ",\"steps\":" << kSteps << "},\n"
		<< "  \"equation_systems\":[{\"name\":\"flow\",\"kind\":\"navier_stokes\","
			"\"unknowns\":[\"velocity\",\"pressure\"],\"viscosity\":" << Number(kViscosity)
		<< ",\"density\":" << Number(kDensity)
		<< ",\"time_integration\":\"backward_euler\"}],\n"
		<< "  \"boundaries\":[\n"
		<< "    {\"label\":0,\"name\":\"wall\",\"conditions\":[{\"field\":\"velocity\","
			"\"type\":\"dirichlet\",\"value\":[0,0,0]}]},\n"
		<< "    {\"label\":1,\"name\":\"inlet\",\"conditions\":[{\"field\":\"velocity\","
			"\"type\":\"dirichlet\",\"profile\":\"initial_velocityfield.txt\",\"scale\":1}]},\n"
		<< "    {\"label\":2,\"name\":\"outlet_a\",\"conditions\":[{\"field\":\"pressure\","
			"\"type\":\"pressure_traction\",\"value\":0}]},\n"
		<< "    {\"label\":3,\"name\":\"outlet_b\",\"conditions\":[{\"field\":\"pressure\","
			"\"type\":\"pressure_traction\",\"value\":0}]}]\n}\n";
	if (!mesh || !velocity || !config)
		throw std::runtime_error("cannot write bifurcation 3D case");
}

void WriteOneDCase(const fs::path& directory, double inlet_flow)
{
	std::ofstream network(directory/"tree.swc");
	network << std::setprecision(17) << "1 2 0 0 0 " << kRadius << " -1\n"
		<< "2 2 1 0 0 " << kRadius << " 1\n";
	std::ofstream config(directory/"simulation_config.json");
	config << "{\n  \"schema_version\":3,\"dimension\":\"1d\","
		"\"simulation_scope\":{\"mode\":\"flow_only\"},\n"
		<< "  \"geometry\":{\"kind\":\"swc_network\",\"file\":\"tree.swc\","
			"\"length_scale_to_m\":1},\n"
		<< "  \"fields\":[{\"name\":\"area\",\"kind\":\"scalar\"},"
			"{\"name\":\"flow_rate\",\"kind\":\"scalar\"},"
			"{\"name\":\"pressure\",\"kind\":\"pressure\"}],\n"
		<< "  \"time\":{\"dt\":" << Number(kDt) << ",\"steps\":" << kSteps
		<< ",\"output_every\":1},\n"
		<< "  \"temporal_functions\":[{\"name\":\"inlet_flow\",\"kind\":\"constant\","
			"\"units\":\"m3/s\",\"value\":" << Number(inlet_flow) << "}],\n"
		<< "  \"equation_systems\":[{\"name\":\"flow\",\"kind\":\"network_flow_1d\","
			"\"unknowns\":[\"area\",\"flow_rate\",\"pressure\"],\"model\":\"rigid\","
			"\"scheme\":\"steady_poiseuille\",\"dynamic_viscosity\":" << Number(kViscosity)
		<< ",\"density\":" << Number(kDensity)
		<< ",\"discretization\":{\"cells_per_segment\":1}}],\n"
		<< "  \"boundaries\":[{\"name\":\"inlet\",\"role\":\"inlet\",\"node_ids\":[1],"
			"\"conditions\":[{\"field\":\"flow_rate\",\"type\":\"dirichlet\","
			"\"quantity\":\"flow_rate\",\"waveform\":\"inlet_flow\"}]},"
			"{\"name\":\"outlet\",\"role\":\"outlet\",\"node_ids\":[2],"
			"\"conditions\":[{\"field\":\"pressure\",\"type\":\"pressure\",\"value\":0}]}]\n}\n";
	if (!network || !config) throw std::runtime_error("cannot write bifurcation 1D case");
}

void WriteGraph(const fs::path& root, const std::string& database,
	const std::string& execution)
{
	const bool explicit_mode = execution == "explicit";
	std::ofstream output(root/"simulation_config.json");
	output << "{\n  \"schema_version\":5,\"time\":{\"dt\":" << Number(kDt)
		<< ",\"steps\":" << kSteps << "},\"start_domain\":\"source\",\n"
		<< "  \"execution\":{\"kind\":\"" << execution << "\",\"maximum_iterations\":"
		<< (explicit_mode ? 1 : 60) << ",\"pressure_relative_tolerance\":1e-6,"
			"\"pressure_reference_pa\":1,\"flow_relative_tolerance\":1e-10,"
			"\"relaxation_factor\":0.5,\"minimum_relaxation\":0.05,"
			"\"maximum_relaxation\":0.999},\n"
		<< "  \"domains\":[\n"
		<< "    {\"id\":\"source\",\"dimension\":\"1d\",\"kind\":\"network_flow\","
			"\"case\":\"source\",\"inlet_policy\":\"configured_open_loop\",\"ports\":["
			"{\"id\":\"root_obs\",\"locator_kind\":\"runtime_port\",\"locator\":\"root\","
			"\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],\"requires\":[]},"
			"{\"id\":\"terminal\",\"locator_kind\":\"runtime_port\",\"locator\":\"outlet:2\","
			"\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],"
			"\"requires\":[\"mean_pressure\"]}]},\n"
		<< "    {\"id\":\"junction\",\"dimension\":\"3d\","
			"\"kind\":\"body_fitted_iga_flow\",\"case\":\"three_d\",\"database\":\""
		<< database << "\",\"ports\":["
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
			"\"provides\":[\"flow_rate\"],\"requires\":[]}]},\n";
	for (const char branch : std::string("ab"))
		output << "    {\"id\":\"branch_" << branch << "\",\"dimension\":\"1d\","
			"\"kind\":\"network_flow\",\"case\":\"branch_" << branch << "\","
			"\"inlet_policy\":\"coupled_root\",\"ports\":["
			"{\"id\":\"root\",\"locator_kind\":\"runtime_port\",\"locator\":\"root\","
			"\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],"
			"\"requires\":[\"flow_rate\"]},"
			"{\"id\":\"terminal\",\"locator_kind\":\"runtime_port\","
			"\"locator\":\"outlet:2\",\"provides\":[\"area\",\"flow_rate\","
			"\"mean_pressure\"],\"requires\":[]}]}" << (branch == 'a' ? ",\n" : "\n");
	output << "  ],\n  \"couplings\":["
		<< "{\"id\":\"edge_in\",\"a\":{\"domain\":\"source\",\"port\":\"terminal\"},"
			"\"b\":{\"domain\":\"junction\",\"port\":\"inlet\"},"
			"\"mode\":\"pressure_flow\",\"initial_pressure_pa\":0},"
		<< "{\"id\":\"edge_a\",\"a\":{\"domain\":\"junction\",\"port\":\"outlet_a\"},"
			"\"b\":{\"domain\":\"branch_a\",\"port\":\"root\"},"
			"\"mode\":\"pressure_flow\",\"initial_pressure_pa\":0.02},"
		<< "{\"id\":\"edge_b\",\"a\":{\"domain\":\"junction\",\"port\":\"outlet_b\"},"
			"\"b\":{\"domain\":\"branch_b\",\"port\":\"root\"},"
			"\"mode\":\"pressure_flow\",\"initial_pressure_pa\":0.01}]}\n";
	if (!output) throw std::runtime_error("cannot write bifurcation graph");
}

std::string Quote(const fs::path& path) { return "'"+path.string()+"'"; }

int Run(const fs::path& root, const fs::path& output, const std::string& launcher = {},
	const std::string& prefix = {})
{
	const auto log = root/(output.filename().string()+".log");
	const std::string linear_solver = launcher.empty()
		? " -ksp_type preonly -pc_type lu"
		: " -ksp_type gmres -pc_type bjacobi -sub_pc_type lu -ksp_rtol 1e-12";
	return std::system((prefix+launcher+"./iga_1d_3d_bifurcation --graph-case "+Quote(root)
		+" --output-dir "+Quote(output)+linear_solver+">"
		+Quote(log)+" 2>&1").c_str());
}

using CsvRow = std::map<std::string, std::string>;

std::vector<CsvRow> ReadCsv(const fs::path& path)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("missing bifurcation CSV: "+path.string());
	std::string line;
	std::getline(input, line);
	std::vector<std::string> names;
	std::istringstream header(line);
	for (std::string name; std::getline(header, name, ',');) names.push_back(name);
	std::vector<CsvRow> rows;
	while (std::getline(input, line)) {
		if (line.empty()) continue;
		std::istringstream values(line);
		CsvRow row;
		for (const auto& name : names) {
			std::string value;
			if (!std::getline(values, value, ','))
				throw std::runtime_error("short bifurcation CSV row");
			row.emplace(name, value);
		}
		rows.push_back(std::move(row));
	}
	return rows;
}

double Value(const CsvRow& row, const std::string& name)
{
	return std::stod(row.at(name));
}

void Validate(const fs::path& output)
{
	if (!fs::is_regular_file(output/"graph_binding_manifest.json"))
		throw std::runtime_error("bifurcation completion marker is missing");
	std::ifstream marker_input(output/"graph_binding_manifest.json");
	const std::string marker((std::istreambuf_iterator<char>(marker_input)),
		std::istreambuf_iterator<char>());
	if (marker.find("\"benchmark\": \"one_d_three_d_bifurcation\"")
		== std::string::npos
		|| marker.find("\"three_d_domain\": \"junction\"") == std::string::npos
		|| marker.find("\"branch_count\": 2") == std::string::npos
		|| marker.find("\"domains\"") != std::string::npos)
		throw std::runtime_error("legacy bifurcation completion marker changed");
	const auto steps = ReadCsv(output/"pressure_flow_steps.csv");
	const auto edges = ReadCsv(output/"pressure_flow_edges.csv");
	const auto iterations = ReadCsv(output/"pressure_flow_iterations.csv");
	const auto ports = ReadCsv(output/"pressure_flow_ports.csv");
	const auto initialization = ReadCsv(output/"one_d_initialization.csv");
	if (steps.size() != kSteps || edges.size() != 3*kSteps || iterations.size() < edges.size()
		|| ports.size() != 10*kSteps)
		throw std::runtime_error("bifurcation long-form row count is invalid");
	if (initialization.size() != 3)
		throw std::runtime_error("bifurcation 1D initialization record is invalid");
	std::map<std::string, double> seeds;
	for (const auto& row : initialization)
		seeds[row.at("domain_id")] = Value(row, "native_inlet_flow_m3_s");
	if (std::abs(seeds.at("source")-1.0e-3) > 1.0e-15
		|| std::abs(seeds.at("branch_a")-6.0e-4) > 1.0e-15
		|| std::abs(seeds.at("branch_b")-4.0e-4) > 1.0e-15)
		throw std::runtime_error("downstream branches did not preserve independent native seeds");
	std::map<std::string, int> edge_counts;
	std::map<int, std::map<std::string, double>> branch_flows;
	for (const auto& edge : edges) {
		++edge_counts[edge.at("edge_id")];
		if (edge.at("edge_id") == "edge_a" || edge.at("edge_id") == "edge_b")
			branch_flows[static_cast<int>(Value(edge, "step"))][edge.at("edge_id")]
				= Value(edge, "first_outward_flow_m3_s");
		if (std::abs(Value(edge, "flow_residual_m3_s")) > 1.0e-13
			|| Value(edge, "normalized_flow_residual") > 1.0e-10)
			throw std::runtime_error("bifurcation interface conservation failed");
	}
	if (edge_counts != std::map<std::string, int>{{"edge_a", kSteps},
		{"edge_b", kSteps}, {"edge_in", kSteps}})
		throw std::runtime_error("bifurcation edge ordering/content is invalid");
	for (const auto& step : branch_flows)
		if (std::abs(step.second.at("edge_a")-step.second.at("edge_b")) < 1.0e-8)
			throw std::runtime_error("bifurcation did not route distinct sibling outlet flows");
	std::map<int, double> junction_sum;
	std::map<int, double> junction_scale;
	std::map<int, double> external_sum;
	std::map<int, double> external_scale;
	for (const auto& port : ports) {
		const int step = static_cast<int>(Value(port, "step"));
		const double flow = Value(port, "outward_flow_m3_s");
		const auto& domain = port.at("domain_id");
		const auto& port_id = port.at("port_id");
		if (domain == "junction") {
			junction_sum[step] += flow;
			junction_scale[step] += std::abs(flow);
		}
		if ((domain == "source" && port_id == "root_obs")
			|| ((domain == "branch_a" || domain == "branch_b")
				&& port_id == "terminal")) {
			external_sum[step] += flow;
			external_scale[step] += std::abs(flow);
		}
	}
	for (const auto& step_row : steps) {
		const int step = static_cast<int>(Value(step_row, "step"));
		const double reported_mass = Value(step_row, "three_d_mass_imbalance_m3_s");
		const double reported_external = Value(step_row, "external_outward_flow_m3_s");
		if (!std::isfinite(reported_mass) || !std::isfinite(reported_external)
			|| std::abs(reported_mass-junction_sum.at(step)) > 1.0e-14
			|| std::abs(reported_external-external_sum.at(step)) > 1.0e-14)
			throw std::runtime_error("bifurcation conservation diagnostic is invalid");
		const double junction_normalized = std::abs(reported_mass)
			/std::max(1.0e-14, junction_scale.at(step));
		const double external_normalized = std::abs(reported_external)
			/std::max(1.0e-14, external_scale.at(step));
		if (junction_normalized > 1.0e-3 || external_normalized > 1.0e-10)
			throw std::runtime_error("bifurcation junction/global conservation gate failed");
	}
}

} // namespace

int main()
{
	try {
		const auto root = fs::temp_directory_path()
			/("tubularflowiga_bifurcation_"+std::to_string(static_cast<long long>(getpid())));
		fs::create_directories(root/"three_d");
		fs::create_directories(root/"source");
		fs::create_directories(root/"branch_a");
		fs::create_directories(root/"branch_b");
		WriteThreeDCase(root/"three_d");
		WriteOneDCase(root/"source", 1.0e-3);
		WriteOneDCase(root/"branch_a", 6.0e-4);
		WriteOneDCase(root/"branch_b", 4.0e-4);
		WriteDatabase(root/"one.ntiga", 1);
		WriteGraph(root, "one.ntiga", "explicit");
		if (Run(root, root/"explicit_one") != 0)
			throw std::runtime_error("one-rank explicit bifurcation run failed");
		Validate(root/"explicit_one");
		if (Run(root, root/"failed", {},
			"TUBULARFLOWIGA_INJECT_BIFURCATION_FAILURE_STEP=1 ") == 0
			|| fs::exists(root/"failed/graph_binding_manifest.json"))
			throw std::runtime_error("failed bifurcation run published a completion marker");
		WriteGraph(root, "one.ntiga", "aitken");
		if (Run(root, root/"aitken_one") != 0)
			throw std::runtime_error("one-rank Aitken bifurcation run failed");
		Validate(root/"aitken_one");
		WriteDatabase(root/"two.ntiga", 2);
		WriteGraph(root, "two.ntiga", "explicit");
		if (Run(root, root/"explicit_two", "mpiexec -np 2 ") != 0)
			throw std::runtime_error("two-rank explicit bifurcation run failed");
		Validate(root/"explicit_two");
		const auto one = ReadCsv(root/"explicit_one/pressure_flow_edges.csv");
		const auto two = ReadCsv(root/"explicit_two/pressure_flow_edges.csv");
		if (one.size() != two.size()) throw std::runtime_error("MPI bifurcation row mismatch");
		for (std::size_t row = 0; row < one.size(); ++row) {
			if (one[row].at("edge_id") != two[row].at("edge_id"))
				throw std::runtime_error("MPI bifurcation edge order mismatch");
			for (const auto& name : {"measured_pressure_pa", "first_outward_flow_m3_s",
				"second_outward_flow_m3_s"})
				if (std::abs(Value(one[row], name)-Value(two[row], name))
					> 1.0e-11*std::max({1.0, std::abs(Value(one[row], name)),
						std::abs(Value(two[row], name))}))
					throw std::runtime_error("MPI bifurcation numerical mismatch");
		}
		fs::remove_all(root);
		std::cout << "bifurcation coupling smoke test passed\n";
		return 0;
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n';
		return 1;
	}
}
