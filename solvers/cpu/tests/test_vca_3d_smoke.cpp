#include "IgaDatabase.hpp"
#include "SimulationConfig.hpp"
#include "VcaCheckpoint.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

namespace {

void WriteUnitDatabase(const fs::path& path, std::uint32_t ranks = 1)
{
	std::ofstream output(path, std::ios::binary | std::ios::trunc);
	if (!output) throw std::runtime_error("cannot create VCA smoke-test database");
	constexpr std::uint64_t header_size = 48;
	const std::uint64_t element_offset = header_size + 2*sizeof(std::uint64_t)
		+ sizeof(std::int32_t);
	output.write(iga::kMagic.data(), iga::kMagic.size());
	iga::Write(output, iga::kVersion);
	iga::Write(output, ranks);
	iga::Write(output, std::uint64_t{1});
	iga::Write(output, std::uint64_t{64});
	iga::Write(output, iga::kBezierPointCount);
	iga::Write(output, std::uint32_t{0});
	const auto rank_index_position = output.tellp();
	iga::Write(output, std::uint64_t{0});
	iga::Write(output, element_offset);
	iga::Write(output, std::uint64_t{0});
	iga::Write(output, std::int32_t{0});
	iga::Write(output, std::uint64_t{0});
	iga::Write(output, std::int32_t{0});
	iga::Write(output, std::int32_t{0});
	iga::Write(output, std::uint32_t{64});
	const std::array<std::int32_t, 6> boundary_labels{{3, 3, 2, 3, 1, 3}};
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
	for (std::uint32_t rank = 0; rank <= ranks; ++rank)
		iga::Write(output, static_cast<std::uint64_t>(rank));
	for (std::uint32_t rank = 0; rank < ranks; ++rank)
		iga::Write(output, std::uint64_t{0});
	output.seekp(header_size + sizeof(std::uint64_t));
	iga::Write(output, rank_index_offset);
	output.seekp(rank_index_position);
	iga::Write(output, rank_index_offset);
	if (!output) throw std::runtime_error("failed to finalize VCA smoke-test database");
}

void WriteCase(const fs::path& directory)
{
	std::ofstream mesh(directory/"controlmesh.vtk");
	if (!mesh) throw std::runtime_error("cannot create VCA smoke-test mesh");
	mesh << "# vtk DataFile Version 3.0\nVCA smoke fixture\nASCII\n"
		<< "DATASET UNSTRUCTURED_GRID\nPOINTS 64 double\n";
	for (int k = 0; k < 4; ++k)
		for (int j = 0; j < 4; ++j)
			for (int i = 0; i < 4; ++i) mesh << i/3.0 << ' ' << j/3.0 << ' ' << k/3.0 << '\n';
	mesh << "CELLS 1 65\n64";
	for (int node = 0; node < 64; ++node) mesh << ' ' << node;
	mesh << "\nCELL_TYPES 1\n72\n";
	mesh << "POINT_DATA 64\nSCALARS boundary_label int 1\nLOOKUP_TABLE default\n";
	for (int k = 0; k < 4; ++k)
		for (int j = 0; j < 4; ++j)
			for (int i = 0; i < 4; ++i)
				mesh << (i == 0 ? 1 : (i == 3 ? 2 : 3)) << '\n';
	std::ofstream velocity(directory/"initial_velocityfield.txt");
	if (!velocity) throw std::runtime_error("cannot create VCA smoke-test velocity field");
	for (int node = 0; node < 64; ++node) velocity << "1 0 0\n";
	std::ofstream config(directory/"simulation_config.json");
	if (!config) throw std::runtime_error("cannot create VCA smoke-test configuration");
	config << R"json({
  "schema_version": 3,
  "dimension": "3d",
  "simulation_scope": {"mode": "vca_closed_loop"},
  "coupling": {
    "scheme": "explicit_staggered",
    "flow_epsilon_m3_s": 1e-14,
    "three_d_ports": {"inlet_label": 1, "outlet_labels": [2, 3]}
  },
  "external_circuit": {
    "reservoir": {"volume_m3": 1.0, "temperature_c": 37.0,
      "species": {"oxygen": 0.2}},
    "pump": {"mode": "flow_control", "flow_m3_s": 0.001}
  },
  "fields": [
    {"name": "velocity", "kind": "vector3"},
    {"name": "pressure", "kind": "pressure"},
    {"name": "oxygen", "kind": "scalar", "initial_value": 0.2}
  ],
  "time": {"dt": 0.1, "steps": 2},
  "equation_systems": [
    {"name": "flow", "kind": "navier_stokes",
      "unknowns": ["velocity", "pressure"], "viscosity": 1.0,
      "density": 1.0, "time_integration": "backward_euler"},
    {"name": "oxygen_transport", "kind": "linear_transport",
      "unknowns": ["oxygen"],
      "terms": [
        {"operator": "time_derivative", "equation": "oxygen"},
        {"operator": "advection", "equation": "oxygen", "velocity": "prescribed"},
        {"operator": "diffusion", "equation": "oxygen", "coefficient": 1e-6}
      ]}
  ],
  "boundaries": [
    {"label": 1, "name": "arterial_inlet", "conditions": [
      {"field": "velocity", "type": "dirichlet",
       "profile": "initial_velocityfield.txt", "scale": 0.001},
      {"field": "oxygen", "type": "dirichlet", "value": 0.2}
    ]},
    {"label": 2, "name": "venous_outlet_a", "conditions": [
      {"field": "pressure", "type": "pressure_traction", "value": 0.0},
      {"field": "oxygen", "type": "advective_outflow"}
    ]},
    {"label": 3, "name": "venous_outlet_b", "conditions": [
      {"field": "pressure", "type": "pressure_traction", "value": 0.0},
      {"field": "oxygen", "type": "advective_outflow"}
    ]}
  ]
})json";
}

std::string Quote(const fs::path& path)
{
	return "'"+path.string()+"'";
}

void Run(const fs::path& database, const fs::path& directory, const fs::path& checkpoint,
	const std::string& extra, const std::string& launcher = {})
{
	fs::create_directories(checkpoint.parent_path());
	const std::string command = launcher+"./iga_navier_stokes "+Quote(database)+" "+Quote(directory)
		+" --max-newton 12 --checkpoint "+Quote(checkpoint)
		+" --output "+Quote(checkpoint.parent_path()/"flow.txt")+extra;
	if (std::system(command.c_str()) != 0)
		throw std::runtime_error("iga_navier_stokes VCA smoke run failed");
}

bool SameFile(const fs::path& first, const fs::path& second)
{
	std::ifstream left(first, std::ios::binary);
	std::ifstream right(second, std::ios::binary);
	if (!left || !right) return false;
	left.seekg(0, std::ios::end);
	right.seekg(0, std::ios::end);
	if (left.tellg() != right.tellg()) return false;
	left.seekg(0);
	right.seekg(0);
	return std::equal(std::istreambuf_iterator<char>(left), std::istreambuf_iterator<char>(),
		std::istreambuf_iterator<char>(right));
}

void RequireCloseFiles(const fs::path& first, const fs::path& second,
	double relative_tolerance, const char* description)
{
	std::ifstream left(first);
	std::ifstream right(second);
	double left_value = 0.0, right_value = 0.0;
	while (left >> left_value) {
		if (!(right >> right_value)
			|| std::abs(left_value-right_value) > relative_tolerance
				*std::max({1.0, std::abs(left_value), std::abs(right_value)}))
			throw std::runtime_error(std::string("single-rank and two-rank ")+description+" differ");
	}
	if (right >> right_value)
		throw std::runtime_error(std::string("single-rank and two-rank ")+description+" sizes differ");
}

void RequireContains(const fs::path& path, const std::string& expected)
{
	std::ifstream input(path);
	std::string text((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
	if (text.find(expected) == std::string::npos)
		throw std::runtime_error("VCA manifest is missing "+expected);
}

void RequireSameReservoir(const fs::path& first, const fs::path& second)
{
	const auto left = iga::ReadVcaCheckpointMetadata(first).reservoir;
	const auto right = iga::ReadVcaCheckpointMetadata(second).reservoir;
	auto close = [](double first, double second) {
		return std::abs(first-second) <= 1e-10*std::max({1.0, std::abs(first), std::abs(second)});
	};
	if (!close(left.volume_m3, right.volume_m3)
		|| !close(left.temperature_c, right.temperature_c)
		|| !close(left.hematocrit_percent, right.hematocrit_percent)
		|| left.species.size() != right.species.size())
		throw std::runtime_error("single-rank and two-rank reservoir states differ");
	for (const auto& species : left.species) {
		const auto found = right.species.find(species.first);
		if (found == right.species.end() || !close(species.second, found->second))
			throw std::runtime_error("single-rank and two-rank reservoir states differ");
	}
}

} // namespace

int main()
{
	const auto root = fs::temp_directory_path()/"tubularflowiga-vca-3d-smoke";
	std::error_code error;
	fs::remove_all(root, error);
	try {
		fs::create_directories(root);
		const auto database = root/"fixture.ntiga";
		WriteUnitDatabase(database);
		WriteCase(root);
		const auto configuration = iga::ReadSimulationConfiguration(
			(root/"simulation_config.json").string());
		if (configuration.coupling.external_circuit.reservoir.species.count("oxygen") != 1)
			throw std::runtime_error("VCA smoke-test configuration lost reservoir oxygen");
		Run(database, root, root/"full/checkpoint", "");
		const auto manifest = root/"full/coupling_manifest.json";
		RequireContains(manifest, "\"backend\": \"cpu_3d_navier_stokes\"");
		RequireContains(manifest, "\"outlet_ids\": [2, 3]");
		RequireContains(manifest, "\"oxygen\"");
		RequireContains(manifest, "\"vca_balance_history\"");
		Run(database, root, root/"split/checkpoint", " --stop-after-step 1");
		Run(database, root, root/"resumed/checkpoint",
			" --restart "+Quote(root/"split/checkpoint"));
		if (!SameFile(root/"full/checkpoint.state", root/"resumed/checkpoint.state"))
			throw std::runtime_error("flow checkpoint/restart state differs from uninterrupted VCA run");
		if (!SameFile(root/"full/checkpoint.vca_transport.state",
			root/"resumed/checkpoint.vca_transport.state"))
			throw std::runtime_error("transport checkpoint/restart state differs from uninterrupted VCA run");
		const auto two_rank_database = root/"fixture-2.ntiga";
		WriteUnitDatabase(two_rank_database, 2);
		Run(two_rank_database, root, root/"two-rank/checkpoint", "", "mpiexec -np 2 ");
		RequireCloseFiles(root/"full/flow.txt", root/"two-rank/flow.txt", 1e-10,
			"velocity field");
		RequireCloseFiles(root/"full/flow.txt.pressure", root/"two-rank/flow.txt.pressure",
			1e-10, "pressure field");
		RequireSameReservoir(root/"full/checkpoint", root/"two-rank/checkpoint");
		fs::remove_all(root);
		std::cout << "VCA 3D flow, transport, checkpoint, and restart smoke test passed\n";
		return 0;
	} catch (const std::exception& exception) {
		std::cerr << "VCA 3D smoke test failed: " << exception.what() << '\n';
		fs::remove_all(root, error);
		return 1;
	}
}
