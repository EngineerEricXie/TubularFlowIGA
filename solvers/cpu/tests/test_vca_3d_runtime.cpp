#include "CouplingPort.hpp"
#include "IgaDatabase.hpp"
#include "TransientFlowRuntime.hpp"
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
	constexpr std::uint64_t header_size = 88;
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

template <class Function>
void RequireRejected(Function&& function, const char* message)
{
	bool rejected = false;
	try {
		function();
	} catch (const std::runtime_error&) {
		rejected = true;
	}
	if (!rejected) throw std::runtime_error(std::string("expected rejection: ")+message);
}

iga::ResolvedBoundaryConditions MakeFlowBoundaries(std::size_t nodes)
{
	iga::ResolvedBoundaryConditions result;
	result.velocity_constrained.assign(nodes, 0);
	result.pressure_constrained.assign(nodes, 0);
	result.transport_constrained.assign(nodes, 0);
	result.velocity.assign(nodes, {0.0, 0.0, 0.0});
	result.pressure.assign(nodes, 0.0);
	result.n0.assign(nodes, 0.0);
	result.nplus.assign(nodes, 0.0);
	return result;
}

iga::CouplingPort MakeBoundaryPort(const std::string& id, const std::string& locator)
{
	iga::CouplingPort result;
	result.id = id;
	result.subsystem_id = "three_d";
	result.locator_kind = "boundary_label";
	result.locator = locator;
	result.provides = {iga::PortQuantity::Area, iga::PortQuantity::FlowRate,
		iga::PortQuantity::MeanPressure, iga::PortQuantity::SpeciesFlux};
	return result;
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
		const auto required_indices = database.RequiredElementIndices(0);
		for (std::uint64_t element = 0; element < database.header().elements; ++element)
			if (database.owners().at(static_cast<std::size_t>(element)) == 0)
				assert(std::find(required_indices.begin(), required_indices.end(), element)
					!= required_indices.end());
		iga::OwnedRowAssembler sparse_assembler(database, PETSC_COMM_WORLD, 6);
		iga::FieldCouplingPattern diagonal(6);
		for (std::size_t field = 0; field < 6; ++field) diagonal.Add(field, field);
		Mat sparse_matrix = sparse_assembler.CreateMatrix(diagonal);
		iga::FieldBlockElementMatrix sparse_blocks(diagonal);
		for (const auto& element : sparse_assembler.elements()) {
			sparse_blocks.Reset(element.connectivity.size());
			for (std::size_t field = 0; field < 6; ++field)
				for (std::size_t a = 0; a < element.connectivity.size(); ++a)
					for (std::size_t b = 0; b < element.connectivity.size(); ++b)
						sparse_blocks.At(field, field, a, b) = 1.0;
			sparse_assembler.AddElementMatrix(sparse_matrix, element, sparse_blocks);
		}
		iga::OwnedRowAssembler::Assemble(sparse_matrix);
		MatInfo sparse_info{};
		MatGetInfo(sparse_matrix, MAT_GLOBAL_SUM, &sparse_info);
		assert(sparse_info.mallocs == 0.0);
		assert(sparse_info.nz_used == 6.0*64.0*64.0);
		assert(sparse_info.nz_allocated == sparse_info.nz_used);
		MatDestroy(&sparse_matrix);
		iga::TransientFlowRuntime flow(database, PETSC_COMM_WORLD, false, false,
			{1.0, 1.0, 0.0}, MakeFlowBoundaries(64), std::vector<int>(64, 1),
			std::vector<std::array<double, 3>>(64, {0.0, 0.0, 0.0}), {}, {});
		std::vector<PetscInt> flow_rows;
		std::vector<PetscScalar> flow_values;
		flow_rows.reserve(4*64);
		flow_values.reserve(4*64);
		for (PetscInt node = 0; node < 64; ++node) {
			for (PetscInt field = 0; field < 4; ++field) {
				flow_rows.push_back(4*node+field);
				flow_values.push_back(field == 2 ? 2.0 : (field == 3 ? 7.0 : 0.0));
			}
		}
		VecSetValues(flow.State(), static_cast<PetscInt>(flow_rows.size()), flow_rows.data(),
			flow_values.data(), INSERT_VALUES);
		iga::OwnedRowAssembler::Assemble(flow.State());
		const std::vector<std::string> species_fields{"oxygen"};
		const std::vector<double> species_state(64, 2.0);
		const auto generic_ports = std::vector<iga::CouplingPort>{
			MakeBoundaryPort("outlet", "1")};
		const auto generic = flow.MeasurePorts(generic_ports, 1.25, species_fields, species_state);
		const auto& generic_state = generic.at("outlet");
		assert(generic_state.area_m2.has_value());
		assert(generic_state.outward_flow_m3_s.has_value());
		assert(generic_state.mean_pressure_pa.has_value());
		RequireNear(1.0, *generic_state.area_m2, 1e-12, "generic port area");
		RequireNear(-2.0, *generic_state.outward_flow_m3_s, 1e-12, "generic port flow");
		RequireNear(7.0, *generic_state.mean_pressure_pa, 1e-12, "generic port pressure");
		RequireNear(-4.0, generic_state.outward_species_flux.at("oxygen"), 1e-12,
			"generic port species flux");
		RequireNear(1.25, generic_state.time_s, 1e-12, "generic port time");
		auto reversed_port = MakeBoundaryPort("reversed_outlet", "1");
		reversed_port.orientation.native_to_outward_sign = -1;
		const auto reversed = flow.MeasurePorts({reversed_port}, 1.25,
			species_fields, species_state);
		RequireNear(2.0, *reversed.at("reversed_outlet").outward_flow_m3_s, 1e-12,
			"generic port orientation flow conversion");
		RequireNear(4.0, reversed.at("reversed_outlet").outward_species_flux.at("oxygen"),
			1e-12, "generic port orientation species flux conversion");
		iga::ThreeDVascularPortDefinition vca_ports;
		vca_ports.outlet_labels = {1};
		const auto vca = flow.MeasurePorts(vca_ports, species_fields, species_state);
		RequireNear(*generic_state.outward_flow_m3_s, vca.flows.at(1), 1e-12,
			"generic and VCA flow parity");
		RequireNear(*generic_state.mean_pressure_pa, vca.pressures.at(1), 1e-12,
			"generic and VCA pressure parity");
		RequireNear(generic_state.outward_species_flux.at("oxygen"),
			vca.species_fluxes.at(1).at("oxygen"), 1e-12,
			"generic and VCA species flux parity");
		RequireRejected([&flow, &species_fields, &species_state] {
			auto invalid = MakeBoundaryPort("invalid_kind", "1");
			invalid.locator_kind = "network_node";
			flow.MeasurePorts({invalid}, 0.0, species_fields, species_state);
		}, "unsupported generic port locator kind");
		RequireRejected([&flow, &species_fields, &species_state] {
			auto invalid = MakeBoundaryPort("invalid_label", "1x");
			flow.MeasurePorts({invalid}, 0.0, species_fields, species_state);
		}, "invalid generic boundary label locator");
		RequireRejected([&flow, &generic_ports, &species_fields, &species_state] {
			auto duplicate = MakeBoundaryPort("duplicate_label", "1");
			flow.MeasurePorts({generic_ports.front(), duplicate}, 0.0, species_fields,
				species_state);
		}, "duplicate generic boundary label locator");
		auto initial_configuration = MakeConfiguration(0.0);
		auto system = iga::CompileLinearSystem(initial_configuration, "oxygen_transport");
		system.velocity_source = "prescribed";
		iga::TransientTransportRuntime runtime(database, PETSC_COMM_WORLD,
			initial_configuration, system, std::vector<int>(64, 1));
		auto step_configuration = MakeConfiguration(2.0);
		assert(runtime.RequiredNodes().size() == 64);
		runtime.Advance(step_configuration, runtime.RequiredNodes(),
			std::vector<std::array<double, 3>>(64, {0.0, 0.0, 0.0}));
		const auto state = runtime.GatherState();
		const auto required_state = runtime.GatherRequiredState();
		assert(state.size() == 64);
		assert(required_state == state);
		for (const auto value : state) RequireNear(2.0, value, 1e-11, "Dirichlet transport state");
		const auto mass = runtime.TotalMass(state);
		RequireNear(2.0, mass.at("oxygen"), 1e-10, "unit-cube oxygen mass");
		RequireNear(2.0, runtime.TotalMass().at("oxygen"), 1e-10,
			"required-node unit-cube oxygen mass");
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
