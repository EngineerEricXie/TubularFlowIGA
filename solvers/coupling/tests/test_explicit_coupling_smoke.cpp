#include "IgaDatabase.hpp"

#include <algorithm>
#include <array>
#include <cassert>
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

constexpr double kFlowM3S = 1.0e-3;
constexpr double kDensityKgM3 = 1.0e-6;
constexpr double kViscosityPaS = 1.0;
constexpr double kDtS = 1.0e-2;
constexpr int kSteps = 3;
constexpr double kInterfaceAreaM2 = 1.0;
constexpr double kRadiusM = 0.56418958354775628; // sqrt(1/pi): circular area = 1 m^2.
constexpr double kSquareDuctResistancePaSM3 = 28.454153* kViscosityPaS;
constexpr double kCircularSegmentResistancePaSM3 = 8.0*kViscosityPaS
	/(3.14159265358979323846*std::pow(kRadiusM, 4));
constexpr double kPressureReferencePa = (2.0*kCircularSegmentResistancePaSM3
	+kSquareDuctResistancePaSM3)*kFlowM3S;
// Predeclared coarse one-element pressure gate; it is not a convergence claim.
constexpr double kPressureRelativeTolerance = 0.5;

std::string JsonNumber(double value)
{
	std::ostringstream output;
	output << std::setprecision(17) << value;
	return output.str();
}

void WriteUnitDatabase(const fs::path& path, std::uint32_t ranks)
{
	std::ofstream output(path, std::ios::binary | std::ios::trunc);
	if (!output) throw std::runtime_error("cannot create explicit coupling smoke database");
	constexpr std::uint64_t header_size = 88;
	const std::uint64_t element_offset = header_size + 2*sizeof(std::uint64_t)+sizeof(std::int32_t);
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
	iga::Write(output, 1.0); // source units per normalized unit
	iga::Write(output, 1.0); // source length scale to metres
	iga::Write(output, element_offset);
	iga::Write(output, std::uint64_t{0});
	iga::Write(output, std::int32_t{0});
	iga::Write(output, std::uint64_t{0});
	iga::Write(output, std::int32_t{0});
	iga::Write(output, std::int32_t{0});
	iga::Write(output, std::uint32_t{64});
	const std::array<std::int32_t, 6> boundary_labels{{0, 0, 2, 0, 1, 0}};
	output.write(reinterpret_cast<const char*>(boundary_labels.data()),
		static_cast<std::streamsize>(sizeof(boundary_labels)));
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
	for (std::uint32_t rank = 0; rank <= ranks; ++rank) iga::Write(output, static_cast<std::uint64_t>(rank));
	for (std::uint32_t rank = 0; rank < ranks; ++rank) iga::Write(output, std::uint64_t{0});
	output.seekp(header_size+sizeof(std::uint64_t));
	iga::Write(output, rank_index_offset);
	output.seekp(rank_index_position);
	iga::Write(output, rank_index_offset);
	if (!output) throw std::runtime_error("cannot finalize explicit coupling smoke database");
}

void WriteThreeDCase(const fs::path& directory)
{
	std::ofstream mesh(directory/"controlmesh.vtk");
	if (!mesh) throw std::runtime_error("cannot create smoke control mesh");
	mesh << "# vtk DataFile Version 3.0\nexplicit coupling smoke\nASCII\n"
		<< "DATASET UNSTRUCTURED_GRID\nPOINTS 64 double\n";
	for (int k = 0; k < 4; ++k)
		for (int j = 0; j < 4; ++j)
			for (int i = 0; i < 4; ++i) mesh << i/3.0 << ' ' << j/3.0 << ' ' << k/3.0 << '\n';
	mesh << "CELLS 1 9\n8 0 3 15 12 48 51 63 60\nCELL_TYPES 1\n12\n"
		<< "POINT_DATA 64\nSCALARS boundary_label int 1\nLOOKUP_TABLE default\n";
	for (int k = 0; k < 4; ++k)
		for (int j = 0; j < 4; ++j)
			for (int i = 0; i < 4; ++i) {
				const bool wall = j == 0 || j == 3 || k == 0 || k == 3;
				mesh << (wall ? 0 : (i == 0 ? 1 : (i == 3 ? 2 : -1))) << '\n';
			}
	std::ofstream velocity(directory/"initial_velocityfield.txt");
	if (!velocity) throw std::runtime_error("cannot create smoke velocity profile");
	for (int k = 0; k < 4; ++k)
		for (int j = 0; j < 4; ++j)
			for (int i = 0; i < 4; ++i) {
				const double y = j/3.0;
				const double z = k/3.0;
				const double axial_profile = i == 0 ? 16.0*y*(1.0-y)*z*(1.0-z) : 0.0;
				velocity << std::setprecision(17) << axial_profile << " 0 0\n";
			}
	std::ofstream config(directory/"simulation_config.json");
	if (!config) throw std::runtime_error("cannot create smoke 3D configuration");
	config << "{\n"
		<< "  \"schema_version\": 3,\n  \"dimension\": \"3d\",\n"
		<< "  \"simulation_scope\": {\"mode\": \"flow_only\"},\n"
		<< "  \"coupling\": {\"scheme\": \"explicit_staggered\", \"flow_epsilon_m3_s\": 1e-14,\n"
		<< "    \"three_d_ports\": {\"inlet_label\": 1, \"outlet_labels\": [2]}},\n"
		<< "  \"fields\": [{\"name\": \"velocity\", \"kind\": \"vector3\"}, {\"name\": \"pressure\", \"kind\": \"pressure\"}],\n"
		<< "  \"time\": {\"dt\": " << JsonNumber(kDtS) << ", \"steps\": " << kSteps << "},\n"
		<< "  \"equation_systems\": [{\"name\": \"flow\", \"kind\": \"navier_stokes\",\n"
		<< "    \"unknowns\": [\"velocity\", \"pressure\"], \"viscosity\": " << JsonNumber(kViscosityPaS)
		<< ", \"density\": " << JsonNumber(kDensityKgM3) << ", \"time_integration\": \"backward_euler\"}],\n"
		<< "  \"boundaries\": [\n"
		<< "    {\"label\": 0, \"name\": \"no_slip_wall\", \"conditions\": [{\"field\": \"velocity\", \"type\": \"dirichlet\", \"value\": [0.0, 0.0, 0.0]}]},\n"
		<< "    {\"label\": 1, \"name\": \"inlet\", \"conditions\": [{\"field\": \"velocity\", \"type\": \"dirichlet\", \"profile\": \"initial_velocityfield.txt\", \"scale\": 1.0}]},\n"
		<< "    {\"label\": 2, \"name\": \"outlet\", \"conditions\": [{\"field\": \"pressure\", \"type\": \"pressure_traction\", \"value\": 0.0}]}\n"
		<< "  ]\n}\n";
}

void WriteOneDCase(const fs::path& directory, double dt_s = kDtS, int steps = kSteps)
{
	std::ofstream network(directory/"tree.swc");
	if (!network) throw std::runtime_error("cannot create smoke 1D network");
	network << std::setprecision(17) << "1 2 0 0 0 " << kRadiusM << " -1\n"
		<< "2 2 1 0 0 " << kRadiusM << " 1\n";
	std::ofstream config(directory/"simulation_config.json");
	if (!config) throw std::runtime_error("cannot create smoke 1D configuration");
	config << "{\n  \"schema_version\": 3, \"dimension\": \"1d\",\n"
		<< "  \"simulation_scope\": {\"mode\": \"flow_only\"},\n"
		<< "  \"geometry\": {\"kind\": \"swc_network\", \"file\": \"tree.swc\", \"length_scale_to_m\": 1.0},\n"
		<< "  \"fields\": [{\"name\": \"area\", \"kind\": \"scalar\"}, {\"name\": \"flow_rate\", \"kind\": \"scalar\"}, {\"name\": \"pressure\", \"kind\": \"pressure\"}],\n"
		<< "  \"time\": {\"dt\": " << JsonNumber(dt_s) << ", \"steps\": " << steps << ", \"output_every\": 1},\n"
		<< "  \"temporal_functions\": [{\"name\": \"inlet_flow\", \"kind\": \"constant\", \"units\": \"m3/s\", \"value\": " << JsonNumber(kFlowM3S) << "}],\n"
		<< "  \"equation_systems\": [{\"name\": \"flow\", \"kind\": \"network_flow_1d\", \"unknowns\": [\"area\", \"flow_rate\", \"pressure\"],\n"
		<< "    \"model\": \"rigid\", \"scheme\": \"steady_poiseuille\", \"dynamic_viscosity\": " << JsonNumber(kViscosityPaS)
		<< ", \"density\": " << JsonNumber(kDensityKgM3) << ", \"discretization\": {\"cells_per_segment\": 1}}],\n"
		<< "  \"boundaries\": [\n"
		<< "    {\"name\": \"inlet\", \"role\": \"inlet\", \"node_ids\": [1], \"conditions\": [{\"field\": \"flow_rate\", \"type\": \"dirichlet\", \"quantity\": \"flow_rate\", \"waveform\": \"inlet_flow\"}]},\n"
		<< "    {\"name\": \"outlet\", \"role\": \"outlet\", \"node_ids\": [2], \"conditions\": [{\"field\": \"pressure\", \"type\": \"pressure\", \"value\": 0.0}]}\n"
		<< "  ]\n}\n";
}

void WriteGraphCase(const fs::path& root, const std::string& execution_kind = "explicit",
	int inlet_native_to_outward_sign = 1, const std::string& database = "one.ntiga")
{
	std::ofstream config(root/"simulation_config.json");
	if (!config) throw std::runtime_error("cannot create multidomain smoke configuration");
	config << "{\n"
		<< "  \"schema_version\": 5,\n"
		<< "  \"time\": {\"dt\": " << JsonNumber(kDtS) << ", \"steps\": "
		<< kSteps << "},\n"
		<< "  \"start_domain\": \"arterial_source\",\n"
		<< "  \"execution\": {\"kind\": \"" << execution_kind << "\", "
		<< "\"maximum_iterations\": " << (execution_kind == "explicit" ? 1 : 50) << ", "
		<< "\"pressure_relative_tolerance\": 1e-6, "
		<< "\"pressure_reference_pa\": " << JsonNumber(kPressureReferencePa) << ", "
		<< "\"flow_relative_tolerance\": 1e-10, \"relaxation_factor\": 0.5, "
		<< "\"minimum_relaxation\": 0.05, \"maximum_relaxation\": 0.999},\n"
		<< "  \"domains\": [\n"
		<< "    {\"id\":\"arterial_source\",\"dimension\":\"1d\",\"kind\":\"network_flow\","
		<< "\"case\":\"upstream\",\"inlet_policy\":\"configured_open_loop\",\"ports\":["
		<< "{\"id\":\"to_roi\",\"locator_kind\":\"runtime_port\",\"locator\":\"outlet:2\","
		<< "\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],\"requires\":[\"mean_pressure\"]},"
		<< "{\"id\":\"source_root\",\"locator_kind\":\"runtime_port\",\"locator\":\"root\","
		<< "\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],\"requires\":[]}]},\n"
		<< "    {\"id\":\"local_roi\",\"dimension\":\"3d\",\"kind\":\"body_fitted_iga_flow\","
		<< "\"case\":\"three_d\",\"database\":\"" << database << "\",\"ports\":["
		<< "{\"id\":\"roi_in\",\"locator_kind\":\"boundary_label\",\"locator\":\"1\","
		<< "\"native_to_outward_sign\":" << inlet_native_to_outward_sign << ","
		<< "\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],\"requires\":[\"flow_rate\"]},"
		<< "{\"id\":\"roi_out\",\"locator_kind\":\"boundary_label\",\"locator\":\"2\","
		<< "\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],\"requires\":[\"mean_pressure\"]},"
		<< "{\"id\":\"wall_flux\",\"locator_kind\":\"boundary_label\",\"locator\":\"0\","
		<< "\"provides\":[\"flow_rate\"],\"requires\":[]}]},\n"
		<< "    {\"id\":\"venous_sink\",\"dimension\":\"1d\",\"kind\":\"network_flow\","
		<< "\"case\":\"downstream\",\"inlet_policy\":\"coupled_root\",\"ports\":["
		<< "{\"id\":\"sink_root\",\"locator_kind\":\"runtime_port\",\"locator\":\"root\","
		<< "\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],\"requires\":[\"flow_rate\"]},"
		<< "{\"id\":\"sink_terminal\",\"locator_kind\":\"runtime_port\",\"locator\":\"outlet:2\","
		<< "\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],\"requires\":[]}] }\n"
		<< "  ],\n"
		<< "  \"couplings\": [\n"
		<< "    {\"id\":\"edge_left\",\"a\":{\"domain\":\"arterial_source\",\"port\":\"to_roi\"},"
		<< "\"b\":{\"domain\":\"local_roi\",\"port\":\"roi_in\"},\"mode\":\"pressure_flow\","
		<< "\"initial_pressure_pa\":0},\n"
		<< "    {\"id\":\"edge_right\",\"a\":{\"domain\":\"local_roi\",\"port\":\"roi_out\"},"
		<< "\"b\":{\"domain\":\"venous_sink\",\"port\":\"sink_root\"},\"mode\":\"pressure_flow\","
		<< "\"initial_pressure_pa\":"
		<< JsonNumber(kCircularSegmentResistancePaSM3*kFlowM3S) << "}\n"
		<< "  ]\n}\n";
}

std::string Quote(const fs::path& path)
{
	return "'"+path.string()+"'";
}

std::string ReadText(const fs::path& path)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot read smoke diagnostic: "+path.string());
	return std::string((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
}

int Run(const fs::path& database, const fs::path& three_d, const fs::path& upstream,
	const fs::path& downstream, const fs::path& output, const std::string& launcher,
	const std::string& prefix = {}, const std::string& arguments = {}, const fs::path& stderr_path = {})
{
	std::string command = prefix+launcher+"./iga_1d_3d_explicit "+Quote(database)+" "+Quote(three_d)
		+" "+Quote(upstream)+" "+Quote(downstream)+" --upstream-terminal-node 2 --output-dir "
		+Quote(output)+" -ksp_type preonly -pc_type lu"+arguments;
	if (!stderr_path.empty()) command += " 2>"+Quote(stderr_path);
	return std::system(command.c_str());
}

int RunGraph(const fs::path& graph_case, const fs::path& output,
	const std::string& launcher = {}, const std::string& arguments = {})
{
	const std::string command = launcher+"./iga_1d_3d_explicit --graph-case "
		+Quote(graph_case)+" --output-dir "+Quote(output)
		+" -ksp_type preonly -pc_type lu"+arguments;
	return std::system(command.c_str());
}

using CsvRow = std::map<std::string, double>;

std::vector<CsvRow> ReadHistory(const fs::path& path)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("missing explicit coupling history");
	std::string header;
	if (!std::getline(input, header)) throw std::runtime_error("empty explicit coupling history");
	std::vector<std::string> names;
	std::istringstream headings(header);
	for (std::string name; std::getline(headings, name, ',');) names.push_back(name);
	if (names.empty()) throw std::runtime_error("invalid explicit coupling history header");
	std::vector<CsvRow> rows;
	for (std::string line; std::getline(input, line);) {
		if (line.empty()) continue;
		std::istringstream values(line);
		CsvRow row;
		for (std::size_t column = 0; column < names.size(); ++column) {
			std::string value;
			if (!std::getline(values, value, ',')) throw std::runtime_error("truncated explicit coupling history row");
			row.emplace(names[column], std::stod(value));
		}
		if (values.rdbuf()->in_avail() != 0) throw std::runtime_error("extra explicit coupling history value");
		rows.push_back(std::move(row));
	}
	return rows;
}

double Value(const CsvRow& row, const char* name)
{
	const auto found = row.find(name);
	if (found == row.end() || !std::isfinite(found->second))
		throw std::runtime_error(std::string("missing or non-finite history value: ")+name);
	return found->second;
}

bool Close(double first, double second, double tolerance = 1.0e-10)
{
	return std::abs(first-second) <= tolerance*std::max({1.0, std::abs(first), std::abs(second)});
}

void RequireManifest(const fs::path& path)
{
	std::ifstream input(path);
	const std::string text((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
	if (!input || text.find("\"scheme\": \"explicit_staggered\"") == std::string::npos
		|| text.find("\"geometry_transform_product_m\": 1") == std::string::npos
		|| text.find("\"initial_lagged_pressures_pa\"") == std::string::npos)
		throw std::runtime_error("explicit coupling manifest is incomplete");
}

void RequireGraphBindingManifest(const fs::path& path, const std::string& execution)
{
	const auto text = ReadText(path);
	if (text.find("\"schema_version\": 5") == std::string::npos
		|| text.find("\"start_domain\": \"arterial_source\"") == std::string::npos
		|| text.find("\"execution\": \""+execution+"\"") == std::string::npos
		|| text.find("\"id\":\"local_roi\"") == std::string::npos
		|| text.find("\"id\":\"edge_left\"") == std::string::npos
		|| text.find("\"domain\":\"arterial_source\",\"port\":\"to_roi\"")
			== std::string::npos
		|| text.find("\"database\":") == std::string::npos)
		throw std::runtime_error("graph binding manifest is incomplete");
}

void ValidateRows(const std::vector<CsvRow>& rows, bool explicit_mode = true)
{
	if (rows.size() != static_cast<std::size_t>(kSteps)) throw std::runtime_error("unexpected smoke history row count");
	double maximum_cap_relative_imbalance = 0.0;
	double maximum_wall_relative_flow = 0.0;
	double maximum_total_three_d_relative_imbalance = 0.0;
	double maximum_external_relative_imbalance = 0.0;
	double maximum_trial_linear_iterations = 0.0;
	for (std::size_t row_index = 0; row_index < rows.size(); ++row_index) {
		const auto& row = rows[row_index];
		if (std::abs(Value(row, "time_s")-(row_index+1)*kDtS) > 1.0e-14)
			throw std::runtime_error("explicit coupling history time does not match the committed step");
		if (!(Value(row, "upstream_root_outward_flow_m3_s") < 0.0)
			|| !(Value(row, "upstream_terminal_outward_flow_m3_s") > 0.0)
			|| !(Value(row, "three_d_inlet_outward_flow_m3_s") < 0.0)
			|| !(Value(row, "three_d_outlet_outward_flow_m3_s") > 0.0)
			|| !(Value(row, "downstream_root_outward_flow_m3_s") < 0.0)
			|| !(Value(row, "downstream_terminal_outward_flow_m3_s") > 0.0))
			throw std::runtime_error("explicit coupling smoke flow sign violation");
		if (std::abs(Value(row, "upstream_three_d_normalized_residual")) > 1.0e-10
			|| std::abs(Value(row, "three_d_downstream_normalized_residual")) > 1.0e-10)
			throw std::runtime_error("explicit coupling interface residual is not near roundoff");
		const double inlet_flow = std::abs(Value(row, "three_d_inlet_outward_flow_m3_s"));
		if (!(inlet_flow > 0.0)) throw std::runtime_error("explicit coupling 3D inlet flow is zero");
		const double cap_relative_imbalance = std::abs(
			Value(row, "three_d_inlet_outward_flow_m3_s")
			+Value(row, "three_d_outlet_outward_flow_m3_s"))/inlet_flow;
		const double wall_relative_flow = std::abs(Value(row, "three_d_wall_outward_flow_m3_s"))/inlet_flow;
		const double total_three_d_relative_imbalance = std::abs(
			Value(row, "three_d_mass_imbalance_m3_s"))/inlet_flow;
		const double external_relative_imbalance = std::abs(
			Value(row, "net_external_outward_flow_m3_s"))/std::abs(
				Value(row, "upstream_root_outward_flow_m3_s"));
		maximum_cap_relative_imbalance = std::max(maximum_cap_relative_imbalance,
			cap_relative_imbalance);
		maximum_wall_relative_flow = std::max(maximum_wall_relative_flow, wall_relative_flow);
		maximum_total_three_d_relative_imbalance = std::max(maximum_total_three_d_relative_imbalance,
			total_three_d_relative_imbalance);
		maximum_external_relative_imbalance = std::max(maximum_external_relative_imbalance,
			external_relative_imbalance);
		if (wall_relative_flow > 1.0e-10)
			throw std::runtime_error("explicit coupling wall flow exceeds no-slip tolerance");
		if (cap_relative_imbalance >= 1.0e-3)
			throw std::runtime_error("explicit coupling cap flow imbalance exceeds tolerance");
		if (total_three_d_relative_imbalance >= 1.0e-3)
			throw std::runtime_error("explicit coupling total 3D mass imbalance exceeds tolerance");
		if (external_relative_imbalance > 1.0e-10)
			throw std::runtime_error("explicit coupling external rigid-1D mass imbalance exceeds tolerance");
		for (const char* area : {"upstream_terminal_area_m2", "three_d_inlet_area_m2",
			"three_d_outlet_area_m2", "downstream_root_area_m2"})
			if (!(Value(row, area) > 0.0) || std::abs(Value(row, area)-kInterfaceAreaM2) > 1.0e-10)
				throw std::runtime_error("explicit coupling interface area mismatch");
		const double trial_linear_iterations = Value(row, "three_d_trial_linear_iterations");
		maximum_trial_linear_iterations = std::max(maximum_trial_linear_iterations,
			trial_linear_iterations);
		if ((explicit_mode && (Value(row, "iteration_count") != 1.0 || Value(row, "relaxation_factor") != 1.0))
			|| (!explicit_mode && (!(Value(row, "iteration_count") >= 1.0)
				|| !(Value(row, "relaxation_factor") > 0.0) || !(Value(row, "relaxation_factor") <= 1.0)))
			|| trial_linear_iterations < 0.0)
			throw std::runtime_error("coupling iteration diagnostics are invalid");
		if (!(Value(row, "external_pressure_drop_pa") > 0.0)
			|| !(Value(row, "upstream_root_pressure_pa") > Value(row, "downstream_terminal_pressure_pa")))
			throw std::runtime_error("explicit coupling pressure does not decrease along the full path");
	}
	const auto& final = rows.back();
	if (!(maximum_trial_linear_iterations > 0.0))
		throw std::runtime_error("explicit coupling smoke recorded no provisional linear iterations");
	const double pressure_drop = Value(final, "external_pressure_drop_pa");
	if (std::abs(pressure_drop-kPressureReferencePa) > kPressureRelativeTolerance*kPressureReferencePa)
		throw std::runtime_error("explicit coupling pressure drop disagrees with the predeclared square-duct reference");
	if (rows.size() > 1 && std::abs(
		Value(rows[rows.size()-2], "external_pressure_drop_pa")-pressure_drop)
		> 0.25*std::max(std::abs(pressure_drop), 1.0e-12))
		throw std::runtime_error("explicit coupling one-step pressure lag does not settle for constant flow");
	std::cout << std::setprecision(12) << "explicit smoke final pressure_drop_pa=" << pressure_drop
		<< " square_duct_plus_two_circular_reference_pa=" << kPressureReferencePa
		<< " circular_three_segment_reference_pa=" << 3.0*kCircularSegmentResistancePaSM3*kFlowM3S
		<< " tolerance=" << kPressureRelativeTolerance
		<< " max_cap_relative_imbalance=" << maximum_cap_relative_imbalance
		<< " max_wall_relative_flow=" << maximum_wall_relative_flow
		<< " max_total_3d_relative_imbalance=" << maximum_total_three_d_relative_imbalance
		<< " max_external_relative_imbalance=" << maximum_external_relative_imbalance << '\n';
}

void RequireSameHistory(const std::vector<CsvRow>& first, const std::vector<CsvRow>& second)
{
	if (first.size() != second.size()) throw std::runtime_error("one-rank/two-rank history lengths differ");
	for (std::size_t row = 0; row < first.size(); ++row) {
		if (first[row].size() != second[row].size()) throw std::runtime_error("one-rank/two-rank history columns differ");
		for (const auto& value : first[row]) {
			const auto peer = second[row].find(value.first);
			if (peer == second[row].end() || !Close(value.second, peer->second))
				throw std::runtime_error("one-rank/two-rank histories differ");
		}
	}
}

void RequireSubcycledWork(const std::vector<CsvRow>& rows, bool strong)
{
	for (const auto& row : rows) {
		for (const char* prefix : {"upstream", "downstream"}) {
			const std::string base(prefix);
			if (Value(row, (base+"_configured_substeps_attempted").c_str()) != 4.0
				|| Value(row, (base+"_explicit_cfl_substeps").c_str()) != 0.0
				|| Value(row, (base+"_accepted_configured_substeps").c_str()) != 4.0
				|| Value(row, (base+"_accepted_explicit_cfl_substeps").c_str()) != 0.0)
				throw std::runtime_error("subcycled 1D attempt diagnostics are invalid");
			const double iterations = Value(row, "iteration_count");
			const double expected_all = strong ? 4.0*iterations : 4.0;
			if (Value(row, (base+"_all_configured_substeps").c_str()) != expected_all
				|| Value(row, (base+"_rejected_configured_substeps").c_str()) != expected_all-4.0
				|| Value(row, (base+"_all_explicit_cfl_substeps").c_str()) != 0.0
				|| Value(row, (base+"_rejected_explicit_cfl_substeps").c_str()) != 0.0)
				throw std::runtime_error("subcycled 1D work accounting is invalid");
		}
	}
}

struct PhysicalHistoryDifference {
	double maximum_pressure_difference_pa = 0.0;
	double maximum_nonpressure_difference = 0.0;
	std::string maximum_pressure_field;
	std::string maximum_nonpressure_field;
};

PhysicalHistoryDifference RequireEquivalentPhysicalHistory(const std::vector<CsvRow>& fixed, const std::vector<CsvRow>& aitken)
{
	if (fixed.size() != aitken.size()) throw std::runtime_error("fixed/Aitken physical step counts differ");
	const double pressure_tolerance = 4.0*kPressureReferencePa*1.0e-6+1.0e-12;
	PhysicalHistoryDifference difference;
	for (std::size_t i = 0; i < fixed.size(); ++i)
		for (const auto& value : fixed[i]) {
			if (value.first.find("iteration") != std::string::npos || value.first.find("ksp") != std::string::npos
				|| value.first.find("relaxation") != std::string::npos || value.first.find("residual") != std::string::npos
				|| value.first.find("substeps") != std::string::npos)
				continue;
			const double peer = Value(aitken[i], value.first.c_str());
			const double absolute_difference = std::abs(value.second-peer);
			const bool pressure = value.first.find("pressure") != std::string::npos || value.first.find("drop") != std::string::npos;
			if (pressure && absolute_difference > difference.maximum_pressure_difference_pa) {
				difference.maximum_pressure_difference_pa = absolute_difference;
				difference.maximum_pressure_field = value.first;
			}
			if (!pressure && absolute_difference > difference.maximum_nonpressure_difference) {
				difference.maximum_nonpressure_difference = absolute_difference;
				difference.maximum_nonpressure_field = value.first;
			}
			const double tolerance = pressure
				? pressure_tolerance : 1.0e-10*std::max({1.0, std::abs(value.second), std::abs(peer)});
			if (absolute_difference > tolerance)
				throw std::runtime_error("fixed/Aitken physical history differs: "+value.first);
		}
	std::cout << std::setprecision(17) << "fixed/Aitken physical equivalence max_pressure_difference_pa="
		<< difference.maximum_pressure_difference_pa << " field=" << difference.maximum_pressure_field
		<< " tolerance_pa=" << pressure_tolerance << " max_nonpressure_difference="
		<< difference.maximum_nonpressure_difference << " field=" << difference.maximum_nonpressure_field
		<< " tolerance_relative=1e-10\n";
	return difference;
}

double ManifestNumber(const fs::path& path, const std::string& key)
{
	std::ifstream input(path);
	const std::string text((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
	const auto position = text.find("\""+key+"\"");
	if (!input || position == std::string::npos) throw std::runtime_error("missing strong manifest number: "+key);
	const auto colon = text.find(':', position);
	const auto end = text.find_first_of(",}\n", colon+1);
	if (colon == std::string::npos || end == std::string::npos)
		throw std::runtime_error("malformed strong manifest number: "+key);
	return std::stod(text.substr(colon+1, end-colon-1));
}

std::string ManifestText(const fs::path& path)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot read coupling manifest");
	return std::string((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
}

double ManifestNumberAfter(const std::string& text, std::size_t start, const std::string& key)
{
	const auto position = text.find("\""+key+"\"", start);
	if (position == std::string::npos) throw std::runtime_error("missing manifest key: "+key);
	const auto colon = text.find(':', position);
	const auto end = text.find_first_of(",}", colon+1);
	if (colon == std::string::npos || end == std::string::npos) throw std::runtime_error("malformed manifest key: "+key);
	return std::stod(text.substr(colon+1, end-colon-1));
}

void RequireManifestSubcycling(const fs::path& path, const std::vector<CsvRow>& rows,
	double configured_dt_s, int configured_steps, int substeps)
{
	const std::string text = ManifestText(path);
	const auto plan = text.find("\"one_d_subcycling\"");
	const auto work = text.find("\"one_d_work_totals\"");
	const auto upstream = text.find("\"upstream\"", plan);
	const auto downstream = text.find("\"downstream\"", upstream);
	const auto upstream_work = text.find("\"upstream\"", work);
	const auto downstream_work = text.find("\"downstream\"", upstream_work);
	if (plan == std::string::npos || work == std::string::npos || upstream == std::string::npos
		|| downstream == std::string::npos || upstream_work == std::string::npos || downstream_work == std::string::npos
		|| text.find("zero-order held interface data") == std::string::npos
		|| text.find("configured open-loop sampling at substep endpoints") == std::string::npos
		|| text.find("endpoint residuals") == std::string::npos
		|| text.find("internal_substeps is explicit CFL only") == std::string::npos)
		throw std::runtime_error("subcycling manifest semantics are incomplete");
	if (!Close(ManifestNumberAfter(text, plan, "macro_dt_s"), kDtS)
		|| ManifestNumberAfter(text, plan, "macro_steps") != kSteps
		|| !Close(ManifestNumberAfter(text, plan, "ratio_relative_tolerance"), 1.0e-12))
		throw std::runtime_error("subcycling manifest macro timing is invalid");
	for (const auto section : {upstream, downstream})
		if (!Close(ManifestNumberAfter(text, section, "configured_dt_s"), configured_dt_s)
			|| ManifestNumberAfter(text, section, "configured_steps") != configured_steps
			|| ManifestNumberAfter(text, section, "N") != substeps)
			throw std::runtime_error("subcycling manifest domain plan is invalid");
	for (const auto& domain : std::array<std::pair<std::string, std::size_t>, 2>{{{"upstream", upstream_work}, {"downstream", downstream_work}}}) {
		long long accepted_configured = 0, all_configured = 0, rejected_configured = 0;
		long long accepted_cfl = 0, all_cfl = 0, rejected_cfl = 0;
		for (const auto& row : rows) {
			accepted_configured += static_cast<long long>(Value(row, (domain.first+"_accepted_configured_substeps").c_str()));
			all_configured += static_cast<long long>(Value(row, (domain.first+"_all_configured_substeps").c_str()));
			rejected_configured += static_cast<long long>(Value(row, (domain.first+"_rejected_configured_substeps").c_str()));
			accepted_cfl += static_cast<long long>(Value(row, (domain.first+"_accepted_explicit_cfl_substeps").c_str()));
			all_cfl += static_cast<long long>(Value(row, (domain.first+"_all_explicit_cfl_substeps").c_str()));
			rejected_cfl += static_cast<long long>(Value(row, (domain.first+"_rejected_explicit_cfl_substeps").c_str()));
		}
		if (ManifestNumberAfter(text, domain.second, "accepted_configured_substeps") != accepted_configured
			|| ManifestNumberAfter(text, domain.second, "all_configured_substeps") != all_configured
			|| ManifestNumberAfter(text, domain.second, "rejected_configured_substeps") != rejected_configured
			|| ManifestNumberAfter(text, domain.second, "accepted_explicit_cfl_substeps") != accepted_cfl
			|| ManifestNumberAfter(text, domain.second, "all_explicit_cfl_substeps") != all_cfl
			|| ManifestNumberAfter(text, domain.second, "rejected_explicit_cfl_substeps") != rejected_cfl)
			throw std::runtime_error("subcycling manifest work total differs from step CSV");
	}
}

void ValidateStrongRun(const std::vector<CsvRow>& history, const std::vector<CsvRow>& iterations,
	const fs::path& manifest, bool fixed_mode = true, int first_step_veto_iteration = 0)
{
	if (history.size() != static_cast<std::size_t>(kSteps) || iterations.size() <= history.size())
		throw std::runtime_error("strong coupling smoke did not perform more than one sweep per step");
	const double pressure_relative_tolerance = ManifestNumber(manifest, "pressure_relative_tolerance");
	const double flow_relative_tolerance = ManifestNumber(manifest, "flow_relative_tolerance");
	const double pressure_reference_pa = ManifestNumber(manifest, "pressure_reference_pa");
	if (!Close(pressure_relative_tolerance, 1.0e-6) || !Close(flow_relative_tolerance, 1.0e-10)
		|| !Close(pressure_reference_pa, kPressureReferencePa))
		throw std::runtime_error("strong manifest controls do not match the smoke configuration");
	ValidateRows(history, false);
	std::map<int, std::vector<CsvRow>> by_step;
	for (const auto& row : iterations) by_step[static_cast<int>(Value(row, "physical_step"))].push_back(row);
	if (by_step.size() != history.size()) throw std::runtime_error("strong iteration history has wrong physical-step groups");
	long long sum_all_ksp = 0;
	long long sum_accepted_ksp = 0;
	bool observed_veto = false;
	for (const auto& entry : by_step) {
		const auto& rows = entry.second;
		if (rows.empty() || rows.size() > 50) throw std::runtime_error("strong iteration count is invalid");
		long long attempt_sum = 0;
		long long upstream_configured_sum = 0, upstream_cfl_sum = 0;
		long long downstream_configured_sum = 0, downstream_cfl_sum = 0;
		for (std::size_t index = 0; index < rows.size(); ++index) {
			if (fixed_mode && (!Close(Value(rows[index], "relaxation_factor_for_next_guess"), 0.5)
				|| !Close(Value(rows[index], "unclamped_relaxation_factor"), 0.5)
				|| Value(rows[index], "aitken_status_code") != -1.0
				|| Value(rows[index], "aitken_has_previous_residual") != 0.0
				|| Value(rows[index], "relaxation_update_applied") != (index+1 == rows.size() ? 0.0 : 1.0)))
				throw std::runtime_error("fixed strong proposal diagnostics are inconsistent");
			const double upstream_q = Value(rows[index], "upstream_terminal_outward_flow_m3_s");
			const double inlet_q = Value(rows[index], "three_d_inlet_outward_flow_m3_s");
			const double outlet_q = Value(rows[index], "three_d_outlet_outward_flow_m3_s");
			const double downstream_q = Value(rows[index], "downstream_root_outward_flow_m3_s");
			const double upstream_flow_residual = upstream_q+inlet_q;
			const double downstream_flow_residual = outlet_q+downstream_q;
			const double normalized_upstream_flow_residual = upstream_flow_residual/std::max({std::abs(upstream_q), std::abs(inlet_q), 1.0e-30});
			const double normalized_downstream_flow_residual = downstream_flow_residual/std::max({std::abs(outlet_q), std::abs(downstream_q), 1.0e-30});
			const double upstream_pressure_raw = Value(rows[index], "measured_three_d_inlet_pressure_pa")
				-Value(rows[index], "applied_upstream_terminal_pressure_pa");
			const double downstream_pressure_raw = Value(rows[index], "measured_downstream_root_pressure_pa")
				-Value(rows[index], "applied_three_d_outlet_traction_pressure_pa");
			const double normalized_upstream_pressure = std::abs(upstream_pressure_raw)/std::max({pressure_reference_pa,
				std::abs(Value(rows[index], "measured_three_d_inlet_pressure_pa")),
				std::abs(Value(rows[index], "applied_upstream_terminal_pressure_pa"))});
			const double normalized_downstream_pressure = std::abs(downstream_pressure_raw)/std::max({pressure_reference_pa,
				std::abs(Value(rows[index], "measured_downstream_root_pressure_pa")),
				std::abs(Value(rows[index], "applied_three_d_outlet_traction_pressure_pa"))});
			if (Value(rows[index], "iteration") != static_cast<double>(index+1)
				|| Value(rows[index], "normalized_upstream_pressure_residual") < 0.0
				|| Value(rows[index], "normalized_downstream_pressure_residual") < 0.0
				|| std::abs(Value(rows[index], "normalized_upstream_three_d_flow_residual")) > 1.0e-10
				|| std::abs(Value(rows[index], "normalized_three_d_downstream_flow_residual")) > 1.0e-10
				|| !Close(Value(rows[index], "signed_upstream_pressure_residual_pa"), upstream_pressure_raw)
				|| !Close(Value(rows[index], "signed_downstream_pressure_residual_pa"), downstream_pressure_raw)
				|| !Close(Value(rows[index], "normalized_upstream_pressure_residual"), normalized_upstream_pressure)
				|| !Close(Value(rows[index], "normalized_downstream_pressure_residual"), normalized_downstream_pressure)
				|| !Close(Value(rows[index], "upstream_three_d_flow_residual_m3_s"), upstream_flow_residual)
				|| !Close(Value(rows[index], "three_d_downstream_flow_residual_m3_s"), downstream_flow_residual)
				|| !Close(Value(rows[index], "normalized_upstream_three_d_flow_residual"), normalized_upstream_flow_residual)
				|| !Close(Value(rows[index], "normalized_three_d_downstream_flow_residual"), normalized_downstream_flow_residual))
				throw std::runtime_error("strong iteration diagnostics are invalid");
			const bool expected_converged = normalized_upstream_pressure <= pressure_relative_tolerance
				&& normalized_downstream_pressure <= pressure_relative_tolerance
				&& std::abs(normalized_upstream_flow_residual) <= flow_relative_tolerance
				&& std::abs(normalized_downstream_flow_residual) <= flow_relative_tolerance;
			// A test may veto exactly one otherwise converged trial on another rank.
			// Its serialized flag must reflect the global decision, with another
			// iteration following it. All numerical checks above remain unchanged.
			const bool veto = first_step_veto_iteration > 0 && entry.first == 1
				&& static_cast<int>(index+1) == first_step_veto_iteration;
			if (veto) {
				if (!expected_converged || index+1 == rows.size())
					throw std::runtime_error("strong injected veto did not reject a converged trial");
				observed_veto = true;
			}
			if (Value(rows[index], "converged") != (expected_converged && !veto ? 1.0 : 0.0))
				throw std::runtime_error("strong serialized convergence flag is incorrect");
			if (fixed_mode && (!Close(Value(rows[index], "next_upstream_terminal_pressure_pa"),
				Value(rows[index], "applied_upstream_terminal_pressure_pa")+0.5*(Value(rows[index], "measured_three_d_inlet_pressure_pa")-Value(rows[index], "applied_upstream_terminal_pressure_pa")))
				|| !Close(Value(rows[index], "next_three_d_outlet_traction_pressure_pa"),
				Value(rows[index], "applied_three_d_outlet_traction_pressure_pa")+0.5*(Value(rows[index], "measured_downstream_root_pressure_pa")-Value(rows[index], "applied_three_d_outlet_traction_pressure_pa")))))
				throw std::runtime_error("strong fixed-relaxation next guess is not exact");
			if (index > 0 && (!Close(Value(rows[index], "applied_upstream_terminal_pressure_pa"),
				Value(rows[index-1], "next_upstream_terminal_pressure_pa"))
				|| !Close(Value(rows[index], "applied_three_d_outlet_traction_pressure_pa"),
				Value(rows[index-1], "next_three_d_outlet_traction_pressure_pa"))))
				throw std::runtime_error("strong next guess was not reapplied after rollback");
			attempt_sum += static_cast<long long>(Value(rows[index], "three_d_attempt_linear_iterations"));
			if (static_cast<long long>(Value(rows[index], "three_d_cumulative_step_linear_iterations")) != attempt_sum)
				throw std::runtime_error("strong cumulative KSP work is invalid");
			const long long upstream_attempt_configured = static_cast<long long>(Value(rows[index], "upstream_1d_attempt_configured_substeps"));
			const long long upstream_attempt_cfl = static_cast<long long>(Value(rows[index], "upstream_1d_attempt_explicit_cfl_substeps"));
			const long long downstream_attempt_configured = static_cast<long long>(Value(rows[index], "downstream_1d_attempt_configured_substeps"));
			const long long downstream_attempt_cfl = static_cast<long long>(Value(rows[index], "downstream_1d_attempt_explicit_cfl_substeps"));
			if (upstream_attempt_configured < 1 || upstream_attempt_cfl < 0
				|| downstream_attempt_configured < 1 || downstream_attempt_cfl < 0)
				throw std::runtime_error("strong iteration has invalid 1D attempt work");
			upstream_configured_sum += upstream_attempt_configured;
			upstream_cfl_sum += upstream_attempt_cfl;
			downstream_configured_sum += downstream_attempt_configured;
			downstream_cfl_sum += downstream_attempt_cfl;
			if (static_cast<long long>(Value(rows[index], "upstream_1d_cumulative_configured_substeps")) != upstream_configured_sum
				|| static_cast<long long>(Value(rows[index], "upstream_1d_cumulative_explicit_cfl_substeps")) != upstream_cfl_sum
				|| static_cast<long long>(Value(rows[index], "downstream_1d_cumulative_configured_substeps")) != downstream_configured_sum
				|| static_cast<long long>(Value(rows[index], "downstream_1d_cumulative_explicit_cfl_substeps")) != downstream_cfl_sum)
				throw std::runtime_error("strong cumulative 1D work is invalid");
		}
		const auto& final = rows.back();
		if (Value(final, "normalized_upstream_pressure_residual") > pressure_relative_tolerance
			|| Value(final, "normalized_downstream_pressure_residual") > pressure_relative_tolerance
			|| Value(final, "converged") != 1.0)
			throw std::runtime_error("strong coupling did not meet pressure tolerance");
		if (entry.first < 1 || entry.first > static_cast<int>(history.size()))
			throw std::runtime_error("strong iteration physical step is invalid");
		const auto& step = history[static_cast<std::size_t>(entry.first-1)];
		if (Value(step, "iteration_count") != static_cast<double>(rows.size())
			|| static_cast<long long>(Value(step, "all_three_d_ksp_iterations")) != attempt_sum
			|| static_cast<long long>(Value(step, "accepted_three_d_ksp_iterations"))
				!= static_cast<long long>(Value(final, "three_d_attempt_linear_iterations"))
			|| static_cast<long long>(Value(step, "rejected_three_d_ksp_iterations"))
			!= attempt_sum-static_cast<long long>(Value(final, "three_d_attempt_linear_iterations")))
			throw std::runtime_error("strong step work accounting does not equal iteration attempts");
		const auto RequireStepOneDWork = [&](const char* domain, long long attempt_configured,
			long long all_configured, long long attempt_cfl, long long all_cfl) {
			const long long accepted_configured = static_cast<long long>(Value(step, (std::string(domain)+"_accepted_configured_substeps").c_str()));
			const long long serialized_all_configured = static_cast<long long>(Value(step, (std::string(domain)+"_all_configured_substeps").c_str()));
			const long long rejected_configured = static_cast<long long>(Value(step, (std::string(domain)+"_rejected_configured_substeps").c_str()));
			const long long accepted_cfl = static_cast<long long>(Value(step, (std::string(domain)+"_accepted_explicit_cfl_substeps").c_str()));
			const long long serialized_all_cfl = static_cast<long long>(Value(step, (std::string(domain)+"_all_explicit_cfl_substeps").c_str()));
			const long long rejected_cfl = static_cast<long long>(Value(step, (std::string(domain)+"_rejected_explicit_cfl_substeps").c_str()));
			if (accepted_configured != attempt_configured || serialized_all_configured != all_configured
				|| rejected_configured != all_configured-attempt_configured || accepted_cfl != attempt_cfl
				|| serialized_all_cfl != all_cfl || rejected_cfl != all_cfl-attempt_cfl)
				throw std::runtime_error(std::string("strong ")+domain+" 1D step work does not equal iteration attempts");
		};
		RequireStepOneDWork("upstream", static_cast<long long>(Value(final, "upstream_1d_attempt_configured_substeps")),
			upstream_configured_sum, static_cast<long long>(Value(final, "upstream_1d_attempt_explicit_cfl_substeps")), upstream_cfl_sum);
		RequireStepOneDWork("downstream", static_cast<long long>(Value(final, "downstream_1d_attempt_configured_substeps")),
			downstream_configured_sum, static_cast<long long>(Value(final, "downstream_1d_attempt_explicit_cfl_substeps")), downstream_cfl_sum);
		sum_all_ksp += attempt_sum;
		sum_accepted_ksp += static_cast<long long>(Value(final, "three_d_attempt_linear_iterations"));
	}
	bool observed_multiple_iterations = false;
	for (const auto& row : history) {
		observed_multiple_iterations = observed_multiple_iterations || Value(row, "iteration_count") > 1.0;
		if (!(Value(row, "iteration_count") >= 1.0) || Value(row, "iteration_count") > 50.0
			|| (fixed_mode && !Close(Value(row, "relaxation_factor"), 0.5))
			|| !Close(Value(row, "pressure_reference_pa"), kPressureReferencePa)
			|| Value(row, "accepted_three_d_ksp_iterations") < 0.0
			|| Value(row, "all_three_d_ksp_iterations") < Value(row, "accepted_three_d_ksp_iterations")
			|| !Close(Value(row, "rejected_three_d_ksp_iterations"), Value(row, "all_three_d_ksp_iterations")-Value(row, "accepted_three_d_ksp_iterations")))
			throw std::runtime_error("strong coupling work accounting is invalid");
		const double updates = Value(row, "relaxation_updates_applied");
		if (fixed_mode && (updates != Value(row, "iteration_count")-1.0
			|| !Close(Value(row, "final_proposed_relaxation_factor"), 0.5)
			|| Value(row, "final_proposal_applied") != 0.0
			|| (updates == 0.0 && Value(row, "last_applied_relaxation_factor") != 0.0)
			|| (updates > 0.0 && !Close(Value(row, "last_applied_relaxation_factor"), 0.5))))
			throw std::runtime_error("strong step relaxation summary is invalid");
	}
	if (first_step_veto_iteration > 0 && !observed_veto)
		throw std::runtime_error("strong injected veto is missing from iteration history");
	if (!observed_multiple_iterations)
		throw std::runtime_error("strong coupling smoke never required a fixed-point correction");
	if (static_cast<long long>(ManifestNumber(manifest, "total_coupling_iterations"))
		!= static_cast<long long>(iterations.size())
		|| static_cast<long long>(ManifestNumber(manifest, "accepted")) != sum_accepted_ksp
		|| static_cast<long long>(ManifestNumber(manifest, "all_attempts")) != sum_all_ksp)
		throw std::runtime_error("strong manifest work totals do not match history");
}

void RequireStrongManifest(const fs::path& path)
{
	std::ifstream input(path);
	const std::string text((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
	if (!input || text.find("\"scheme\": \"strong_fixed\"") == std::string::npos
		|| text.find("pressure_fixed_point") == std::string::npos
		|| text.find("pressure-traction parameter") == std::string::npos
		|| text.find("geometry_transform_product_m") == std::string::npos
		|| text.find("newton_controls") == std::string::npos
		|| text.find("\"restart\": \"unsupported\"") == std::string::npos
		|| text.find("aitken_controls") != std::string::npos
		|| text.find("aitken_recurrence") != std::string::npos)
		throw std::runtime_error("strong coupling manifest is incomplete");
}

void RequireAitkenManifest(const fs::path& path)
{
	std::ifstream input(path);
	const std::string text((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
	if (!input || text.find("\"scheme\": \"strong_aitken\"") == std::string::npos
		|| text.find("\"aitken_recurrence\"") == std::string::npos
		|| text.find("omega_hat=-omega_previous") == std::string::npos
		|| text.find("delta=(r_current/scale)-(r_previous/scale)") == std::string::npos
		|| text.find("retains prior clamped omega") == std::string::npos
		|| text.find("reset each macro-step") == std::string::npos
		|| text.find("final converged proposal is unused") == std::string::npos)
		throw std::runtime_error("Aitken coupling manifest is incomplete");
}

void ValidateAitkenRun(const std::vector<CsvRow>& history, const std::vector<CsvRow>& iterations,
	const fs::path& manifest, std::size_t fixed_first_step_iterations)
{
	ValidateStrongRun(history, iterations, manifest, false);
	ValidateRows(history, false);
	const double reference = ManifestNumber(manifest, "pressure_reference_pa");
	const double initial = ManifestNumber(manifest, "initial_relaxation");
	const double minimum = ManifestNumber(manifest, "minimum_relaxation");
	const double maximum = ManifestNumber(manifest, "maximum_relaxation");
	const double threshold = ManifestNumber(manifest, "scaled_difference_threshold");
	if (!Close(initial, 0.5) || !Close(minimum, 0.05) || !Close(maximum, 0.999)
		|| !(threshold > 0.0))
		throw std::runtime_error("Aitken manifest controls are invalid");
	std::map<int, std::vector<CsvRow>> by_step;
	for (const auto& row : iterations) by_step[static_cast<int>(Value(row, "physical_step"))].push_back(row);
	bool dynamic = false;
	std::array<long long, 6> status_counts{{0, 0, 0, 0, 0, 0}};
	for (const auto& entry : by_step) {
		const auto& rows = entry.second;
		std::array<double, 2> previous_residual{{0.0, 0.0}};
		double previous_omega = initial;
		bool has_previous = false;
		for (std::size_t i = 0; i < rows.size(); ++i) {
			const auto& row = rows[i];
			const std::array<double, 2> x{{Value(row, "applied_upstream_terminal_pressure_pa"),
				Value(row, "applied_three_d_outlet_traction_pressure_pa")}};
			const std::array<double, 2> residual{{Value(row, "signed_upstream_pressure_residual_pa"),
				Value(row, "signed_downstream_pressure_residual_pa")}};
			double omega = previous_omega, unclamped = previous_omega, numerator = 0.0, denominator = 0.0;
			int status = 0;
			if (has_previous) {
				double scale = reference;
				for (int c = 0; c < 2; ++c) scale = std::max(scale, std::max(std::abs(residual[c]), std::abs(previous_residual[c])));
				for (int c = 0; c < 2; ++c) {
					const double prior = previous_residual[c]/scale;
					const double difference = residual[c]/scale-prior;
					numerator += prior*difference;
					denominator += difference*difference;
				}
				if (denominator <= threshold) status = 4;
				else {
					unclamped = -previous_omega*numerator/denominator;
					if (!std::isfinite(unclamped)) status = 5;
					else if (unclamped < minimum) { omega = minimum; status = 2; }
					else if (unclamped > maximum) { omega = maximum; status = 3; }
					else { omega = std::max(minimum, std::min(maximum, unclamped)); status = 1; }
				}
			}
			if (Value(row, "aitken_has_previous_residual") != (has_previous ? 1.0 : 0.0)
				|| !Close(Value(row, "aitken_scaled_numerator"), numerator)
				|| !Close(Value(row, "aitken_scaled_denominator"), denominator)
				|| !Close(Value(row, "unclamped_relaxation_factor"), unclamped)
				|| !Close(Value(row, "relaxation_factor_for_next_guess"), omega)
				|| Value(row, "aitken_status_code") != static_cast<double>(status)
				|| !Close(Value(row, "next_upstream_terminal_pressure_pa"), x[0]+omega*residual[0])
				|| !Close(Value(row, "next_three_d_outlet_traction_pressure_pa"), x[1]+omega*residual[1]))
				throw std::runtime_error("Aitken CSV proposal is not independently reproducible");
			dynamic = dynamic || std::abs(Value(row, "relaxation_factor_for_next_guess")-0.5) > 1.0e-12;
			if (Value(row, "relaxation_update_applied") != (i+1 == rows.size() ? 0.0 : 1.0))
				throw std::runtime_error("Aitken applied-update flag is invalid");
			++status_counts[static_cast<std::size_t>(status)];
			if (Value(row, "relaxation_update_applied") != 0.0) {
				previous_residual = residual;
				previous_omega = omega;
				has_previous = true;
			}
		}
		const auto& step = history.at(static_cast<std::size_t>(entry.first-1));
		const double updates = Value(step, "relaxation_updates_applied");
		if (updates != static_cast<double>(rows.size()-1)
			|| Value(step, "final_proposal_applied") != 0.0
			|| !Close(Value(step, "final_proposed_relaxation_factor"),
				Value(rows.back(), "relaxation_factor_for_next_guess")))
			throw std::runtime_error("Aitken step relaxation summary is invalid");
		if (rows.size() == 1) {
			if (Value(step, "last_applied_relaxation_factor") != 0.0
				|| !Close(Value(step, "relaxation_factor"), initial))
				throw std::runtime_error("one-sweep Aitken step relaxation summary is invalid");
		} else if (!Close(Value(step, "last_applied_relaxation_factor"),
			Value(rows[rows.size()-2], "relaxation_factor_for_next_guess"))
			|| !Close(Value(step, "relaxation_factor"), Value(step, "last_applied_relaxation_factor")))
			throw std::runtime_error("Aitken last-applied relaxation summary is invalid");
	}
	if (static_cast<long long>(ManifestNumber(manifest, "initial")) != status_counts[0]
		|| static_cast<long long>(ManifestNumber(manifest, "dynamic")) != status_counts[1]
		|| static_cast<long long>(ManifestNumber(manifest, "clamped_minimum")) != status_counts[2]
		|| static_cast<long long>(ManifestNumber(manifest, "clamped_maximum")) != status_counts[3]
		|| static_cast<long long>(ManifestNumber(manifest, "tiny_difference_fallback")) != status_counts[4]
		|| static_cast<long long>(ManifestNumber(manifest, "nonfinite_candidate_fallback")) != status_counts[5])
		throw std::runtime_error("Aitken manifest status totals do not match CSV rows");
	if (!dynamic || by_step.at(1).size() >= fixed_first_step_iterations)
		throw std::runtime_error("Aitken did not dynamically reduce first-step strong sweeps");
}

} // namespace

int main()
{
	const auto root = fs::temp_directory_path()/(
		"tubularflowiga-explicit-coupling-smoke-"+std::to_string(static_cast<long long>(getpid())));
	const bool keep_output = std::getenv("TUBULARFLOWIGA_KEEP_TEST_OUTPUT") != nullptr;
	std::error_code error;
	fs::remove_all(root, error);
	try {
		fs::create_directories(root/"three_d");
		fs::create_directories(root/"upstream");
		fs::create_directories(root/"downstream");
		WriteThreeDCase(root/"three_d");
		WriteOneDCase(root/"upstream");
		WriteOneDCase(root/"downstream");
		fs::create_directories(root/"upstream_subcycled");
		fs::create_directories(root/"downstream_subcycled");
		WriteOneDCase(root/"upstream_subcycled", 0.0025, 12);
		WriteOneDCase(root/"downstream_subcycled", 0.0025, 12);
		fs::create_directories(root/"upstream_noninteger");
		fs::create_directories(root/"upstream_bad_horizon");
		WriteOneDCase(root/"upstream_noninteger", 0.003, 10);
		WriteOneDCase(root/"upstream_bad_horizon", 0.0025, 11);
		WriteUnitDatabase(root/"one.ntiga", 1);
		WriteGraphCase(root);
		fs::create_directories(root/"stale_graph_output");
		std::ofstream(root/"stale_graph_output/graph_binding_manifest.json") << "stale-marker\n";
		if (RunGraph(root, root/"stale_graph_output") == 0
			|| ReadText(root/"stale_graph_output/graph_binding_manifest.json") != "stale-marker\n")
			throw std::runtime_error("graph rerun accepted or modified a stale completion marker");
		if (RunGraph(root, root/"invalid_graph_cli_controls", "",
			" --coupling-mode explicit") == 0
			|| fs::exists(root/"invalid_graph_cli_controls/explicit_coupling_history.csv"))
			throw std::runtime_error("graph case accepted conflicting CLI coupling controls");
		if (Run(root/"one.ntiga", root/"three_d", root/"upstream", root/"downstream", root/"invalid_explicit_aitken_min", "", "",
			" --coupling-mode explicit --strong-aitken-min-relaxation 0.05") == 0
			|| fs::exists(root/"invalid_explicit_aitken_min/explicit_coupling_history.csv")
			|| fs::exists(root/"invalid_explicit_aitken_min/strong_coupling_history.csv"))
			throw std::runtime_error("explicit mode accepted an Aitken minimum option or wrote output");
		if (Run(root/"one.ntiga", root/"three_d", root/"upstream", root/"downstream", root/"invalid_fixed_aitken_max", "", "",
			" --coupling-mode strong-fixed --strong-pressure-reference-pa "+JsonNumber(kPressureReferencePa)
			+" --strong-aitken-max-relaxation 0.9") == 0
			|| fs::exists(root/"invalid_fixed_aitken_max/explicit_coupling_history.csv")
			|| fs::exists(root/"invalid_fixed_aitken_max/strong_coupling_history.csv"))
			throw std::runtime_error("strong-fixed mode accepted an Aitken maximum option or wrote output");
		const std::array<std::string, 3> invalid_newton_options{{
			" --three-d-max-newton", " --three-d-max-newton 0", " --three-d-max-newton nope"}};
		for (std::size_t invalid_index = 0; invalid_index < invalid_newton_options.size(); ++invalid_index) {
			const auto output = root/("invalid_newton_"+std::to_string(invalid_index));
			if (Run(root/"one.ntiga", root/"three_d", root/"upstream", root/"downstream", output, "", "",
				invalid_newton_options[invalid_index]) == 0
				|| fs::exists(output/"explicit_coupling_history.csv")
				|| fs::exists(output/"explicit_coupling_manifest.json")
				|| fs::exists(output/"strong_coupling_history.csv")
				|| fs::exists(output/"strong_coupling_iterations.csv")
				|| fs::exists(output/"strong_coupling_manifest.json"))
				throw std::runtime_error("invalid --three-d-max-newton unexpectedly succeeded or wrote output");
		}
		if (Run(root/"one.ntiga", root/"three_d", root/"upstream_noninteger", root/"downstream_subcycled", root/"invalid_noninteger_subcycling", "") == 0
			|| fs::exists(root/"invalid_noninteger_subcycling/explicit_coupling_history.csv"))
			throw std::runtime_error("noninteger 1D subcycling ratio unexpectedly succeeded");
		if (Run(root/"one.ntiga", root/"three_d", root/"upstream_bad_horizon", root/"downstream_subcycled", root/"invalid_horizon_subcycling", "") == 0
			|| fs::exists(root/"invalid_horizon_subcycling/explicit_coupling_history.csv"))
			throw std::runtime_error("mismatched 1D subcycling horizon unexpectedly succeeded");
		if (Run(root/"one.ntiga", root/"three_d", root/"upstream", root/"downstream", root/"one", "", "", " --three-d-max-newton 31") != 0)
			throw std::runtime_error("one-rank explicit coupling smoke run failed");
		const auto one = ReadHistory(root/"one/explicit_coupling_history.csv");
		RequireManifest(root/"one/explicit_coupling_manifest.json");
		if (ManifestNumber(root/"one/explicit_coupling_manifest.json", "maximum_iterations") != 31.0)
			throw std::runtime_error("explicit manifest did not record custom 3D Newton budget");
		RequireManifestSubcycling(root/"one/explicit_coupling_manifest.json", one, kDtS, kSteps, 1);
		ValidateRows(one);
		if (RunGraph(root, root/"graph_one") != 0)
			throw std::runtime_error("schema-v5 graph-authored coupling smoke run failed");
		const auto graph_one = ReadHistory(root/"graph_one/explicit_coupling_history.csv");
		RequireManifest(root/"graph_one/explicit_coupling_manifest.json");
		RequireGraphBindingManifest(root/"graph_one/graph_binding_manifest.json", "explicit");
		ValidateRows(graph_one);
		RequireSameHistory(one, graph_one);
		WriteGraphCase(root, "explicit", -1);
		if (RunGraph(root, root/"graph_reversed_native_orientation") != 0)
			throw std::runtime_error("schema-v5 reversed native orientation smoke run failed");
		const auto graph_reversed = ReadHistory(
			root/"graph_reversed_native_orientation/explicit_coupling_history.csv");
		RequireGraphBindingManifest(
			root/"graph_reversed_native_orientation/graph_binding_manifest.json", "explicit");
		double maximum_reversed_residual = 0.0;
		for (const auto& row : graph_reversed)
			maximum_reversed_residual = std::max(maximum_reversed_residual,
				std::abs(Value(row, "upstream_three_d_flow_residual_m3_s")));
		if (maximum_reversed_residual > 1.0e-12)
			throw std::runtime_error("reversed native orientation violated upstream flow conservation");
		WriteGraphCase(root);
		WriteUnitDatabase(root/"two.ntiga", 2);
		if (Run(root/"two.ntiga", root/"three_d", root/"upstream", root/"downstream", root/"two", "mpiexec -np 2 ") != 0)
			throw std::runtime_error("two-rank explicit coupling smoke run failed");
		const auto two = ReadHistory(root/"two/explicit_coupling_history.csv");
		RequireManifest(root/"two/explicit_coupling_manifest.json");
		if (ManifestNumber(root/"two/explicit_coupling_manifest.json", "maximum_iterations") != 30.0)
			throw std::runtime_error("default explicit manifest did not retain the 3D Newton budget of 30");
		ValidateRows(two);
		RequireSameHistory(one, two);
		WriteGraphCase(root, "explicit", 1, "two.ntiga");
		if (RunGraph(root, root/"graph_two", "mpiexec -np 2 ") != 0)
			throw std::runtime_error("two-rank schema-v5 graph coupling smoke run failed");
		const auto graph_two = ReadHistory(root/"graph_two/explicit_coupling_history.csv");
		RequireGraphBindingManifest(root/"graph_two/graph_binding_manifest.json", "explicit");
		ValidateRows(graph_two);
		RequireSameHistory(two, graph_two);
		const std::string strong_arguments = " --coupling-mode strong-fixed --strong-max-iterations 50"
			" --strong-pressure-relative-tol 1e-6 --strong-pressure-reference-pa "+JsonNumber(kPressureReferencePa)
			+" --strong-flow-relative-tol 1e-10 --strong-relaxation 0.5 --three-d-max-newton 32";
		if (Run(root/"one.ntiga", root/"three_d", root/"upstream", root/"downstream", root/"strong_one", "", "", strong_arguments) != 0)
			throw std::runtime_error("one-rank strong coupling smoke run failed");
		const auto strong_one = ReadHistory(root/"strong_one/strong_coupling_history.csv");
		const auto strong_one_iterations = ReadHistory(root/"strong_one/strong_coupling_iterations.csv");
		RequireStrongManifest(root/"strong_one/strong_coupling_manifest.json");
		if (ManifestNumber(root/"strong_one/strong_coupling_manifest.json", "maximum_iterations") != 32.0)
			throw std::runtime_error("strong manifest did not record custom 3D Newton budget");
		RequireManifestSubcycling(root/"strong_one/strong_coupling_manifest.json", strong_one, kDtS, kSteps, 1);
		ValidateStrongRun(strong_one, strong_one_iterations, root/"strong_one/strong_coupling_manifest.json");
		WriteGraphCase(root, "fixed");
		if (RunGraph(root, root/"graph_fixed") != 0)
			throw std::runtime_error("schema-v5 fixed graph-authored coupling smoke run failed");
		const auto graph_fixed = ReadHistory(root/"graph_fixed/strong_coupling_history.csv");
		const auto graph_fixed_iterations = ReadHistory(
			root/"graph_fixed/strong_coupling_iterations.csv");
		ValidateStrongRun(graph_fixed, graph_fixed_iterations,
			root/"graph_fixed/strong_coupling_manifest.json");
		RequireGraphBindingManifest(root/"graph_fixed/graph_binding_manifest.json", "fixed");
		RequireSameHistory(strong_one, graph_fixed);
		RequireSameHistory(strong_one_iterations, graph_fixed_iterations);
		if (Run(root/"two.ntiga", root/"three_d", root/"upstream", root/"downstream", root/"strong_two", "mpiexec -np 2 ", "", strong_arguments) != 0)
			throw std::runtime_error("two-rank strong coupling smoke run failed");
		const auto strong_two = ReadHistory(root/"strong_two/strong_coupling_history.csv");
		const auto strong_two_iterations = ReadHistory(root/"strong_two/strong_coupling_iterations.csv");
		RequireStrongManifest(root/"strong_two/strong_coupling_manifest.json");
		ValidateStrongRun(strong_two, strong_two_iterations, root/"strong_two/strong_coupling_manifest.json");
		RequireSameHistory(strong_one, strong_two);
		RequireSameHistory(strong_one_iterations, strong_two_iterations);
		const std::string aitken_arguments = " --coupling-mode strong-aitken --strong-max-iterations 50"
			" --strong-pressure-relative-tol 1e-6 --strong-pressure-reference-pa "+JsonNumber(kPressureReferencePa)
			+" --strong-flow-relative-tol 1e-10 --strong-relaxation 0.5 --strong-aitken-min-relaxation 0.05 --strong-aitken-max-relaxation 0.999";
		if (Run(root/"one.ntiga", root/"three_d", root/"upstream", root/"downstream", root/"aitken_one", "", "", aitken_arguments) != 0)
			throw std::runtime_error("one-rank Aitken coupling smoke run failed");
		const auto aitken_one = ReadHistory(root/"aitken_one/strong_coupling_history.csv");
		const auto aitken_one_iterations = ReadHistory(root/"aitken_one/strong_coupling_iterations.csv");
		RequireAitkenManifest(root/"aitken_one/strong_coupling_manifest.json");
		RequireManifestSubcycling(root/"aitken_one/strong_coupling_manifest.json", aitken_one, kDtS, kSteps, 1);
		ValidateAitkenRun(aitken_one, aitken_one_iterations, root/"aitken_one/strong_coupling_manifest.json",
			static_cast<std::size_t>(Value(strong_one.front(), "iteration_count")));
		RequireEquivalentPhysicalHistory(strong_one, aitken_one);
		WriteGraphCase(root, "aitken");
		if (RunGraph(root, root/"graph_aitken") != 0)
			throw std::runtime_error("schema-v5 Aitken graph-authored coupling smoke run failed");
		const auto graph_aitken = ReadHistory(root/"graph_aitken/strong_coupling_history.csv");
		const auto graph_aitken_iterations = ReadHistory(
			root/"graph_aitken/strong_coupling_iterations.csv");
		ValidateAitkenRun(graph_aitken, graph_aitken_iterations,
			root/"graph_aitken/strong_coupling_manifest.json",
			static_cast<std::size_t>(Value(graph_fixed.front(), "iteration_count")));
		RequireGraphBindingManifest(root/"graph_aitken/graph_binding_manifest.json", "aitken");
		RequireSameHistory(aitken_one, graph_aitken);
		RequireSameHistory(aitken_one_iterations, graph_aitken_iterations);
		if (Run(root/"two.ntiga", root/"three_d", root/"upstream", root/"downstream", root/"aitken_two", "mpiexec -np 2 ", "", aitken_arguments) != 0)
			throw std::runtime_error("two-rank Aitken coupling smoke run failed");
		const auto aitken_two = ReadHistory(root/"aitken_two/strong_coupling_history.csv");
		const auto aitken_two_iterations = ReadHistory(root/"aitken_two/strong_coupling_iterations.csv");
		RequireAitkenManifest(root/"aitken_two/strong_coupling_manifest.json");
		ValidateAitkenRun(aitken_two, aitken_two_iterations, root/"aitken_two/strong_coupling_manifest.json",
			static_cast<std::size_t>(Value(strong_two.front(), "iteration_count")));
		RequireSameHistory(aitken_one, aitken_two);
		RequireSameHistory(aitken_one_iterations, aitken_two_iterations);
		if (Run(root/"one.ntiga", root/"three_d", root/"upstream_subcycled", root/"downstream_subcycled", root/"sub_explicit_one", "") != 0)
			throw std::runtime_error("one-rank subcycled explicit coupling smoke run failed");
		const auto sub_explicit_one = ReadHistory(root/"sub_explicit_one/explicit_coupling_history.csv");
		RequireManifestSubcycling(root/"sub_explicit_one/explicit_coupling_manifest.json", sub_explicit_one, 0.0025, 12, 4);
		ValidateRows(sub_explicit_one);
		RequireSubcycledWork(sub_explicit_one, false);
		RequireEquivalentPhysicalHistory(one, sub_explicit_one);
		if (Run(root/"two.ntiga", root/"three_d", root/"upstream_subcycled", root/"downstream_subcycled", root/"sub_explicit_two", "mpiexec -np 2 ") != 0)
			throw std::runtime_error("two-rank subcycled explicit coupling smoke run failed");
		const auto sub_explicit_two = ReadHistory(root/"sub_explicit_two/explicit_coupling_history.csv");
		RequireSubcycledWork(sub_explicit_two, false);
		RequireSameHistory(sub_explicit_one, sub_explicit_two);
		if (Run(root/"one.ntiga", root/"three_d", root/"upstream_subcycled", root/"downstream_subcycled", root/"sub_strong_one", "", "", strong_arguments) != 0)
			throw std::runtime_error("one-rank subcycled fixed strong coupling smoke run failed");
		const auto sub_strong_one = ReadHistory(root/"sub_strong_one/strong_coupling_history.csv");
		const auto sub_strong_one_iterations = ReadHistory(root/"sub_strong_one/strong_coupling_iterations.csv");
		ValidateStrongRun(sub_strong_one, sub_strong_one_iterations, root/"sub_strong_one/strong_coupling_manifest.json");
		RequireManifestSubcycling(root/"sub_strong_one/strong_coupling_manifest.json", sub_strong_one, 0.0025, 12, 4);
		RequireSubcycledWork(sub_strong_one, true);
		RequireEquivalentPhysicalHistory(strong_one, sub_strong_one);
		if (Run(root/"two.ntiga", root/"three_d", root/"upstream_subcycled", root/"downstream_subcycled", root/"sub_strong_two", "mpiexec -np 2 ", "", strong_arguments) != 0)
			throw std::runtime_error("two-rank subcycled fixed strong coupling smoke run failed");
		const auto sub_strong_two = ReadHistory(root/"sub_strong_two/strong_coupling_history.csv");
		RequireSubcycledWork(sub_strong_two, true);
		RequireSameHistory(sub_strong_one, sub_strong_two);
		if (Run(root/"one.ntiga", root/"three_d", root/"upstream_subcycled", root/"downstream_subcycled", root/"sub_aitken_one", "", "", aitken_arguments) != 0)
			throw std::runtime_error("one-rank subcycled Aitken coupling smoke run failed");
		const auto sub_aitken_one = ReadHistory(root/"sub_aitken_one/strong_coupling_history.csv");
		const auto sub_aitken_one_iterations = ReadHistory(root/"sub_aitken_one/strong_coupling_iterations.csv");
		ValidateAitkenRun(sub_aitken_one, sub_aitken_one_iterations, root/"sub_aitken_one/strong_coupling_manifest.json",
			static_cast<std::size_t>(Value(sub_strong_one.front(), "iteration_count")));
		RequireManifestSubcycling(root/"sub_aitken_one/strong_coupling_manifest.json", sub_aitken_one, 0.0025, 12, 4);
		RequireSubcycledWork(sub_aitken_one, true);
		RequireEquivalentPhysicalHistory(sub_strong_one, sub_aitken_one);
		if (Run(root/"two.ntiga", root/"three_d", root/"upstream_subcycled", root/"downstream_subcycled", root/"sub_aitken_two", "mpiexec -np 2 ", "", aitken_arguments) != 0)
			throw std::runtime_error("two-rank subcycled Aitken coupling smoke run failed");
		const auto sub_aitken_two = ReadHistory(root/"sub_aitken_two/strong_coupling_history.csv");
		RequireSubcycledWork(sub_aitken_two, true);
		RequireSameHistory(sub_aitken_one, sub_aitken_two);
		const auto require_subcycled_failure = [&root](const std::string& name, const std::string& prefix,
			const std::string& arguments, const std::string& expected_final_tuple) {
			const auto output = root/name;
			const auto log = root/(name+".log");
			if (Run(root/"one.ntiga", root/"three_d", root/"upstream_subcycled", root/"downstream_subcycled", output, "", prefix, arguments, log) == 0)
				throw std::runtime_error("subcycled strong failure unexpectedly succeeded: "+name);
			if (fs::exists(output/"strong_coupling_history.csv") || fs::exists(output/"strong_coupling_iterations.csv")
				|| fs::exists(output/"strong_coupling_manifest.json"))
				throw std::runtime_error("subcycled strong failure wrote persistent output: "+name);
			const std::string text = ReadText(log);
			if (text.find("iteration=") == std::string::npos || text.find("G=(") == std::string::npos
				|| text.find("normalized_q=(") == std::string::npos
				|| text.find("upstream_1d_work=(") == std::string::npos
				|| text.find("downstream_1d_work=(") == std::string::npos
				|| text.find(expected_final_tuple) == std::string::npos)
				throw std::runtime_error("subcycled strong failure lacks complete final work diagnostics: "+name);
		};
		require_subcycled_failure("sub_fixed_max_one", "",
			" --coupling-mode strong-fixed --strong-max-iterations 1 --strong-pressure-relative-tol 1e-20"
			" --strong-pressure-reference-pa "+JsonNumber(kPressureReferencePa)+" --strong-flow-relative-tol 1e-10 --strong-relaxation 0.5",
			"upstream_1d_work=(4,0,4,0)");
		require_subcycled_failure("sub_aitken_max_one", "",
			" --coupling-mode strong-aitken --strong-max-iterations 1 --strong-pressure-relative-tol 1e-20"
			" --strong-pressure-reference-pa "+JsonNumber(kPressureReferencePa)+" --strong-flow-relative-tol 1e-10 --strong-relaxation 0.5",
			"upstream_1d_work=(4,0,4,0)");
		require_subcycled_failure("sub_fixed_injected", "TUBULARFLOWIGA_INJECT_EXPLICIT_COUPLING_FAILURE_STEP=1 ", strong_arguments,
			"upstream_1d_work=(4,0,84,0)");
		require_subcycled_failure("sub_aitken_injected", "TUBULARFLOWIGA_INJECT_EXPLICIT_COUPLING_FAILURE_STEP=1 ", aitken_arguments,
			"upstream_1d_work=(4,0,12,0)");
		const auto aitken_nonconverged_log = root/"aitken_nonconverged.log";
		if (Run(root/"one.ntiga", root/"three_d", root/"upstream", root/"downstream", root/"aitken_nonconverged", "", "",
			" --coupling-mode strong-aitken --strong-max-iterations 1 --strong-pressure-relative-tol 1e-20"
			" --strong-pressure-reference-pa "+JsonNumber(kPressureReferencePa)+" --strong-flow-relative-tol 1e-10 --strong-relaxation 0.5", aitken_nonconverged_log) == 0)
			throw std::runtime_error("Aitken max-iteration failure unexpectedly succeeded");
		if (fs::exists(root/"aitken_nonconverged/strong_coupling_history.csv")
			|| fs::exists(root/"aitken_nonconverged/strong_coupling_iterations.csv")
			|| fs::exists(root/"aitken_nonconverged/strong_coupling_manifest.json"))
			throw std::runtime_error("nonconverged Aitken coupling wrote persistent output");
		if (ReadText(aitken_nonconverged_log).find("iteration=1") == std::string::npos
			|| ReadText(aitken_nonconverged_log).find("G=(") == std::string::npos
			|| ReadText(aitken_nonconverged_log).find("normalized_q=(") == std::string::npos
			|| ReadText(aitken_nonconverged_log).find("omega=") == std::string::npos)
			throw std::runtime_error("Aitken max-iteration failure did not emit a complete iteration diagnostic");
		const auto aitken_rejected_log = root/"aitken_rejected.log";
		if (Run(root/"one.ntiga", root/"three_d", root/"upstream", root/"downstream", root/"aitken_rejected", "",
			"TUBULARFLOWIGA_INJECT_EXPLICIT_COUPLING_FAILURE_STEP=1 ", aitken_arguments, aitken_rejected_log) == 0)
			throw std::runtime_error("injected Aitken coupling failure unexpectedly succeeded");
		if (fs::exists(root/"aitken_rejected/strong_coupling_history.csv")
			|| fs::exists(root/"aitken_rejected/strong_coupling_iterations.csv")
			|| fs::exists(root/"aitken_rejected/strong_coupling_manifest.json"))
			throw std::runtime_error("rejected Aitken coupling step wrote persistent output");
		if (ReadText(aitken_rejected_log).find("iteration=3") == std::string::npos
			|| ReadText(aitken_rejected_log).find("G=(") == std::string::npos
			|| ReadText(aitken_rejected_log).find("normalized_p=(") == std::string::npos
			|| ReadText(aitken_rejected_log).find("ksp=(") == std::string::npos)
			throw std::runtime_error("injected Aitken failure did not emit its complete final iteration diagnostic");
		if (Run(root/"one.ntiga", root/"three_d", root/"upstream", root/"downstream", root/"strong_nonconverged", "", "",
			" --coupling-mode strong-fixed --strong-max-iterations 1 --strong-pressure-relative-tol 1e-20"
			" --strong-pressure-reference-pa "+JsonNumber(kPressureReferencePa)+" --strong-flow-relative-tol 1e-10 --strong-relaxation 0.5") == 0)
			throw std::runtime_error("strong max-iteration failure unexpectedly succeeded");
		if (fs::exists(root/"strong_nonconverged/strong_coupling_history.csv")
			|| fs::exists(root/"strong_nonconverged/strong_coupling_iterations.csv")
			|| fs::exists(root/"strong_nonconverged/strong_coupling_manifest.json"))
			throw std::runtime_error("nonconverged strong coupling wrote persistent output");
		if (Run(root/"one.ntiga", root/"three_d", root/"upstream", root/"downstream", root/"strong_rejected", "",
			"TUBULARFLOWIGA_INJECT_EXPLICIT_COUPLING_FAILURE_STEP=1 ", strong_arguments) == 0)
			throw std::runtime_error("injected strong coupling failure unexpectedly succeeded");
		if (fs::exists(root/"strong_rejected/strong_coupling_history.csv")
			|| fs::exists(root/"strong_rejected/strong_coupling_iterations.csv")
			|| fs::exists(root/"strong_rejected/strong_coupling_manifest.json"))
			throw std::runtime_error("rejected strong coupling step wrote persistent output");
		if (Run(root/"one.ntiga", root/"three_d", root/"upstream", root/"downstream", root/"rejected", "",
			"TUBULARFLOWIGA_INJECT_EXPLICIT_COUPLING_FAILURE_STEP=1 ") == 0)
			throw std::runtime_error("injected explicit coupling failure unexpectedly succeeded");
		if (fs::exists(root/"rejected/explicit_coupling_history.csv")
			|| fs::exists(root/"rejected/explicit_coupling_manifest.json"))
			throw std::runtime_error("rejected explicit coupling step wrote persistent output");
		if (!keep_output) fs::remove_all(root);
		std::cout << "explicit 1D--3D--1D MPI smoke test passed\n";
		return 0;
	} catch (const std::exception& exception) {
		std::cerr << "explicit coupling smoke test failed: " << exception.what() << '\n';
		if (!keep_output) fs::remove_all(root, error);
		return 1;
	}
}
