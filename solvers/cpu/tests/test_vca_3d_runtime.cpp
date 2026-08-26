#include "IgaDatabase.hpp"
#include "TransientTransportRuntime.hpp"

#include <petscksp.h>

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void WriteUnitDatabase(const std::filesystem::path& path)
{
	std::ofstream output(path, std::ios::binary | std::ios::trunc);
	if (!output) throw std::runtime_error("cannot create VCA runtime test database");
	constexpr std::uint64_t header_size = 48;
	const std::uint64_t element_offset = header_size + 2*sizeof(std::uint64_t)
		+ sizeof(std::int32_t);
	std::uint64_t element_end = element_offset;
	const std::uint32_t ranks = 1;
	const std::uint64_t elements = 1;
	const std::uint64_t nodes = 64;
	output.write(iga::kMagic.data(), iga::kMagic.size());
	iga::Write(output, iga::kVersion);
	iga::Write(output, ranks);
	iga::Write(output, elements);
	iga::Write(output, nodes);
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
	const std::array<std::int32_t, 6> boundary_labels{{1, -1, -1, -1, -1, -1}};
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
	element_end = static_cast<std::uint64_t>(output.tellp());
	const auto rank_index_offset = element_end;
	iga::Write(output, std::uint64_t{0});
	iga::Write(output, std::uint64_t{1});
	iga::Write(output, std::uint64_t{0});
	output.seekp(header_size + sizeof(std::uint64_t));
	iga::Write(output, element_end);
	output.seekp(rank_index_position);
	iga::Write(output, rank_index_offset);
	if (!output) throw std::runtime_error("failed to finalize VCA runtime test database");
}

iga::SimulationConfiguration MakeConfiguration(double boundary_value)
{
	iga::SimulationConfiguration configuration;
	configuration.fields = {{"oxygen", iga::FieldKind::Scalar, 0.0}};
	configuration.time = {0.1, 1};
	iga::EquationSystemDefinition system;
	system.name = "oxygen_transport";
	system.kind = iga::EquationKind::LinearTransport;
	system.unknowns = {"oxygen"};
	system.terms = {{iga::TermKind::TimeDerivative, "oxygen", "oxygen", 1.0, ""}};
	configuration.equation_systems.push_back(system);
	iga::FieldBoundaryCondition oxygen;
	oxygen.field = "oxygen";
	oxygen.kind = iga::FieldBoundaryKind::Dirichlet;
	oxygen.value = {boundary_value};
	configuration.boundaries.push_back({1, "vca_inlet", {oxygen}});
	configuration.velocity_sources.push_back({"in_memory_flow", "prescribed", "", "", "error"});
	return configuration;
}

void RequireNear(double expected, double actual, double tolerance, const char* message)
{
	if (std::abs(expected-actual) > tolerance) {
		std::cerr << message << ": expected " << expected << ", got " << actual << '\n';
		std::abort();
	}
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	try {
		const auto database_path = std::filesystem::temp_directory_path()
			/ "tubularflowiga-vca-runtime-test.ntiga";
		const auto state_path = std::filesystem::temp_directory_path()
			/ "tubularflowiga-vca-runtime-test.state";
		WriteUnitDatabase(database_path);
		iga::Database database(database_path.string());
		auto initial_configuration = MakeConfiguration(0.0);
		auto system = iga::CompileLinearSystem(initial_configuration, "oxygen_transport");
		system.velocity_source = "prescribed";
		iga::TransientTransportRuntime runtime(database, PETSC_COMM_WORLD,
			initial_configuration, system, std::vector<int>(64, 1));
		auto step_configuration = MakeConfiguration(2.0);
		runtime.Advance(step_configuration, std::vector<std::array<double, 3>>(64, {0.0, 0.0, 0.0}));
		const auto state = runtime.GatherState();
		assert(state.size() == 64);
		for (const auto value : state) RequireNear(2.0, value, 1e-11, "Dirichlet transport state");
		const auto mass = runtime.TotalMass(state);
		RequireNear(2.0, mass.at("oxygen"), 1e-10, "unit-cube oxygen mass");
		RequireNear(0.0, runtime.SourceIntegrals().at("oxygen"), 1e-12, "oxygen source integral");
		runtime.WriteState(state_path);
		iga::TransientTransportRuntime restored(database, PETSC_COMM_WORLD,
			initial_configuration, system, std::vector<int>(64, 1));
		restored.ReadState(state_path);
		const auto restored_state = restored.GatherState();
		for (const auto value : restored_state)
			RequireNear(2.0, value, 1e-11, "restored transport state");
		std::filesystem::remove(state_path);
		std::filesystem::remove(database_path);
		std::cout << "VCA 3D transient transport runtime passed\n";
	} catch (const std::exception& error) {
		std::cerr << "VCA 3D transient transport runtime failed: " << error.what() << '\n';
		PetscFinalize();
		return 1;
	}
	PetscFinalize();
	return 0;
}
