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

void WriteOneDCase(const fs::path& directory)
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
		<< "  \"time\": {\"dt\": " << JsonNumber(kDtS) << ", \"steps\": " << kSteps << ", \"output_every\": 1},\n"
		<< "  \"temporal_functions\": [{\"name\": \"inlet_flow\", \"kind\": \"constant\", \"units\": \"m3/s\", \"value\": " << JsonNumber(kFlowM3S) << "}],\n"
		<< "  \"equation_systems\": [{\"name\": \"flow\", \"kind\": \"network_flow_1d\", \"unknowns\": [\"area\", \"flow_rate\", \"pressure\"],\n"
		<< "    \"model\": \"rigid\", \"scheme\": \"steady_poiseuille\", \"dynamic_viscosity\": " << JsonNumber(kViscosityPaS)
		<< ", \"density\": " << JsonNumber(kDensityKgM3) << ", \"discretization\": {\"cells_per_segment\": 1}}],\n"
		<< "  \"boundaries\": [\n"
		<< "    {\"name\": \"inlet\", \"role\": \"inlet\", \"node_ids\": [1], \"conditions\": [{\"field\": \"flow_rate\", \"type\": \"dirichlet\", \"quantity\": \"flow_rate\", \"waveform\": \"inlet_flow\"}]},\n"
		<< "    {\"name\": \"outlet\", \"role\": \"outlet\", \"node_ids\": [2], \"conditions\": [{\"field\": \"pressure\", \"type\": \"pressure\", \"value\": 0.0}]}\n"
		<< "  ]\n}\n";
}

std::string Quote(const fs::path& path)
{
	return "'"+path.string()+"'";
}

int Run(const fs::path& database, const fs::path& three_d, const fs::path& upstream,
	const fs::path& downstream, const fs::path& output, const std::string& launcher,
	const std::string& prefix = {}, const std::string& arguments = {})
{
	const std::string command = prefix+launcher+"./iga_1d_3d_explicit "+Quote(database)+" "+Quote(three_d)
		+" "+Quote(upstream)+" "+Quote(downstream)+" --upstream-terminal-node 2 --output-dir "
		+Quote(output)+" -ksp_type preonly -pc_type lu"+arguments;
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

void ValidateStrongRun(const std::vector<CsvRow>& history, const std::vector<CsvRow>& iterations,
	const fs::path& manifest)
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
	for (const auto& entry : by_step) {
		const auto& rows = entry.second;
		if (rows.empty() || rows.size() > 50) throw std::runtime_error("strong iteration count is invalid");
		long long attempt_sum = 0;
		for (std::size_t index = 0; index < rows.size(); ++index) {
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
			if (Value(rows[index], "converged") != (expected_converged ? 1.0 : 0.0))
				throw std::runtime_error("strong serialized convergence flag is incorrect");
			if (!Close(Value(rows[index], "next_upstream_terminal_pressure_pa"),
				Value(rows[index], "applied_upstream_terminal_pressure_pa")+0.5*(Value(rows[index], "measured_three_d_inlet_pressure_pa")-Value(rows[index], "applied_upstream_terminal_pressure_pa")))
				|| !Close(Value(rows[index], "next_three_d_outlet_traction_pressure_pa"),
				Value(rows[index], "applied_three_d_outlet_traction_pressure_pa")+0.5*(Value(rows[index], "measured_downstream_root_pressure_pa")-Value(rows[index], "applied_three_d_outlet_traction_pressure_pa"))))
				throw std::runtime_error("strong fixed-relaxation next guess is not exact");
			if (index > 0 && (!Close(Value(rows[index], "applied_upstream_terminal_pressure_pa"),
				Value(rows[index-1], "next_upstream_terminal_pressure_pa"))
				|| !Close(Value(rows[index], "applied_three_d_outlet_traction_pressure_pa"),
				Value(rows[index-1], "next_three_d_outlet_traction_pressure_pa"))))
				throw std::runtime_error("strong next guess was not reapplied after rollback");
			attempt_sum += static_cast<long long>(Value(rows[index], "three_d_attempt_linear_iterations"));
			if (static_cast<long long>(Value(rows[index], "three_d_cumulative_step_linear_iterations")) != attempt_sum)
				throw std::runtime_error("strong cumulative KSP work is invalid");
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
		sum_all_ksp += attempt_sum;
		sum_accepted_ksp += static_cast<long long>(Value(final, "three_d_attempt_linear_iterations"));
	}
	bool observed_multiple_iterations = false;
	for (const auto& row : history) {
		observed_multiple_iterations = observed_multiple_iterations || Value(row, "iteration_count") > 1.0;
		if (!(Value(row, "iteration_count") >= 1.0) || Value(row, "iteration_count") > 50.0
			|| !Close(Value(row, "relaxation_factor"), 0.5)
			|| !Close(Value(row, "pressure_reference_pa"), kPressureReferencePa)
			|| Value(row, "accepted_three_d_ksp_iterations") < 0.0
			|| Value(row, "all_three_d_ksp_iterations") < Value(row, "accepted_three_d_ksp_iterations")
			|| !Close(Value(row, "rejected_three_d_ksp_iterations"), Value(row, "all_three_d_ksp_iterations")-Value(row, "accepted_three_d_ksp_iterations")))
			throw std::runtime_error("strong coupling work accounting is invalid");
	}
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
		|| text.find("\"restart\": \"unsupported\"") == std::string::npos)
		throw std::runtime_error("strong coupling manifest is incomplete");
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
		WriteUnitDatabase(root/"one.ntiga", 1);
		if (Run(root/"one.ntiga", root/"three_d", root/"upstream", root/"downstream", root/"one", "") != 0)
			throw std::runtime_error("one-rank explicit coupling smoke run failed");
		const auto one = ReadHistory(root/"one/explicit_coupling_history.csv");
		RequireManifest(root/"one/explicit_coupling_manifest.json");
		ValidateRows(one);
		WriteUnitDatabase(root/"two.ntiga", 2);
		if (Run(root/"two.ntiga", root/"three_d", root/"upstream", root/"downstream", root/"two", "mpiexec -np 2 ") != 0)
			throw std::runtime_error("two-rank explicit coupling smoke run failed");
		const auto two = ReadHistory(root/"two/explicit_coupling_history.csv");
		RequireManifest(root/"two/explicit_coupling_manifest.json");
		ValidateRows(two);
		RequireSameHistory(one, two);
		const std::string strong_arguments = " --coupling-mode strong-fixed --strong-max-iterations 50"
			" --strong-pressure-relative-tol 1e-6 --strong-pressure-reference-pa "+JsonNumber(kPressureReferencePa)
			+" --strong-flow-relative-tol 1e-10 --strong-relaxation 0.5";
		if (Run(root/"one.ntiga", root/"three_d", root/"upstream", root/"downstream", root/"strong_one", "", "", strong_arguments) != 0)
			throw std::runtime_error("one-rank strong coupling smoke run failed");
		const auto strong_one = ReadHistory(root/"strong_one/strong_coupling_history.csv");
		const auto strong_one_iterations = ReadHistory(root/"strong_one/strong_coupling_iterations.csv");
		RequireStrongManifest(root/"strong_one/strong_coupling_manifest.json");
		ValidateStrongRun(strong_one, strong_one_iterations, root/"strong_one/strong_coupling_manifest.json");
		if (Run(root/"two.ntiga", root/"three_d", root/"upstream", root/"downstream", root/"strong_two", "mpiexec -np 2 ", "", strong_arguments) != 0)
			throw std::runtime_error("two-rank strong coupling smoke run failed");
		const auto strong_two = ReadHistory(root/"strong_two/strong_coupling_history.csv");
		const auto strong_two_iterations = ReadHistory(root/"strong_two/strong_coupling_iterations.csv");
		RequireStrongManifest(root/"strong_two/strong_coupling_manifest.json");
		ValidateStrongRun(strong_two, strong_two_iterations, root/"strong_two/strong_coupling_manifest.json");
		RequireSameHistory(strong_one, strong_two);
		RequireSameHistory(strong_one_iterations, strong_two_iterations);
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
