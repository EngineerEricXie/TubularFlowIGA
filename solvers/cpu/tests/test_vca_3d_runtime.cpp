#include "CouplingPort.hpp"
#include "IgaDatabase.hpp"
#include "ThreeDBodyFittedFlowDomainAdapter.hpp"
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

std::vector<int> MakeLifecycleLabels(std::size_t nodes)
{
	std::vector<int> labels(nodes, 1);
	labels.front() = 2;
	return labels;
}

iga::ResolvedBoundaryConditions MakeLifecycleFlowBoundaries(const std::vector<int>& labels)
{
	auto result = MakeFlowBoundaries(labels.size());
	for (std::size_t node = 0; node < labels.size(); ++node) {
		if (labels[node] == 1) result.velocity_constrained[node] = 1;
		if (labels[node] == 2) result.pressure_constrained[node] = 1;
	}
	result.velocity_nodes = static_cast<std::size_t>(std::count(
		result.velocity_constrained.begin(), result.velocity_constrained.end(), 1));
	result.pressure_nodes = static_cast<std::size_t>(std::count(
		result.pressure_constrained.begin(), result.pressure_constrained.end(), 1));
	return result;
}

iga::SimulationConfiguration MakeFlowConfiguration(double velocity_x, double pressure)
{
	iga::SimulationConfiguration configuration;
	configuration.fields = {{"velocity", iga::FieldKind::Vector3, 0.0},
		{"pressure", iga::FieldKind::Pressure, 0.0}};
	configuration.time = {0.1, 1};
	iga::EquationSystemDefinition system;
	system.name = "flow";
	system.kind = iga::EquationKind::NavierStokes;
	system.unknowns = {"velocity", "pressure"};
	system.viscosity = 1.0;
	system.density = 1.0;
	system.time_integration = "backward_euler";
	configuration.equation_systems.push_back(system);
	iga::FieldBoundaryCondition velocity;
	velocity.field = "velocity";
	velocity.kind = iga::FieldBoundaryKind::Dirichlet;
	velocity.value = {velocity_x, 0.0, 0.0};
	iga::FieldBoundaryCondition pressure_condition;
	pressure_condition.field = "pressure";
	pressure_condition.kind = iga::FieldBoundaryKind::PressureTraction;
	pressure_condition.value = {pressure};
	configuration.boundaries.push_back({1, "constrained", {velocity, pressure_condition}});
	iga::FieldBoundaryCondition pressure_gauge;
	pressure_gauge.field = "pressure";
	pressure_gauge.kind = iga::FieldBoundaryKind::Dirichlet;
	pressure_gauge.value = {0.0};
	configuration.boundaries.push_back({2, "pressure_gauge", {pressure_gauge}});
	return configuration;
}

iga::SimulationConfiguration MakeOutletFlowConfiguration()
{
	auto configuration = MakeFlowConfiguration(0.125, 0.0);
	configuration.boundaries.front().conditions.pop_back();
	iga::FieldBoundaryCondition outlet;
	outlet.field = "pressure";
	outlet.kind = iga::FieldBoundaryKind::WindkesselRC;
	outlet.resistance = 1.0;
	outlet.capacitance = 1e-8;
	outlet.reference_pressure = 1.0;
	outlet.initial_pressure = 3.0;
	configuration.boundaries.front().conditions.push_back(outlet);
	return configuration;
}

iga::CouplingPort MakePressureInputPort()
{
	iga::CouplingPort result;
	result.id = "pressure_input";
	result.subsystem_id = "three_d";
	result.locator_kind = "boundary_label";
	result.locator = "1";
	result.provides = {iga::PortQuantity::Area, iga::PortQuantity::FlowRate};
	result.requires = {iga::PortQuantity::MeanPressure};
	return result;
}

std::vector<double> CopyVector(Vec vector)
{
	PetscInt local_size = 0;
	VecGetLocalSize(vector, &local_size);
	const PetscScalar* values = nullptr;
	VecGetArrayRead(vector, &values);
	std::vector<double> result(static_cast<std::size_t>(local_size));
	for (PetscInt i = 0; i < local_size; ++i)
		result[static_cast<std::size_t>(i)] = PetscRealPart(values[i]);
	VecRestoreArrayRead(vector, &values);
	return result;
}

double OutletCapacitorPressure(const iga::TransientFlowRuntime& runtime)
{
	return runtime.OutletModels().front().capacitor_pressure;
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

		const auto lifecycle_labels = MakeLifecycleLabels(64);
		iga::TransientFlowRuntime lifecycle(database, PETSC_COMM_WORLD, true, true,
			{1.0, 1.0, 0.1}, MakeLifecycleFlowBoundaries(lifecycle_labels),
			lifecycle_labels,
			std::vector<std::array<double, 3>>(64, {0.0, 0.0, 0.0}), {}, {});
		const auto initial_flow_configuration = MakeFlowConfiguration(0.0, 0.0);
		const auto trial_flow_configuration = MakeFlowConfiguration(0.125, 2.0);
		lifecycle.InitializeState(initial_flow_configuration);
		const auto committed_state = CopyVector(lifecycle.State());
		RequireRejected([&lifecycle] { lifecycle.SolveTrial(); },
			"SolveTrial before BeginStep");
		RequireRejected([&lifecycle] { lifecycle.CommitStep(); },
			"CommitStep before a solved trial");
		lifecycle.BeginStep(0, 0.1, 12, 1e-8, 1e-8, 3.0);
		RequireRejected([&lifecycle] {
			lifecycle.BeginStep(0, 0.1, 12, 1e-8, 1e-8, 3.0);
		}, "nested BeginStep");
		iga::PortBoundaryData pressure_input;
		pressure_input.time_s = 0.1;
		pressure_input.mean_pressure_pa = 4.0;
		const auto pressure_port = MakePressureInputPort();
		RequireRejected([&lifecycle, &pressure_port, &pressure_input] {
			lifecycle.SetPortInput(pressure_port, pressure_input);
		}, "port input before configured trial boundaries");
		lifecycle.SetTrialBoundaryConfiguration(trial_flow_configuration);
		auto no_velocity_constraint_configuration = trial_flow_configuration;
		no_velocity_constraint_configuration.boundaries.front().conditions.erase(
			no_velocity_constraint_configuration.boundaries.front().conditions.begin());
		RequireRejected([&lifecycle, &no_velocity_constraint_configuration] {
			lifecycle.SetTrialBoundaryConfiguration(no_velocity_constraint_configuration);
		}, "removing a velocity constraint changes fixed PETSc boundary rows");
		auto added_velocity_constraint_configuration = trial_flow_configuration;
		added_velocity_constraint_configuration.boundaries.back().conditions.push_back(
			trial_flow_configuration.boundaries.front().conditions.front());
		RequireRejected([&lifecycle, &added_velocity_constraint_configuration] {
			lifecycle.SetTrialBoundaryConfiguration(added_velocity_constraint_configuration);
		}, "adding a velocity constraint changes fixed PETSc boundary rows");
		auto pressure_dirichlet_configuration = trial_flow_configuration;
		pressure_dirichlet_configuration.boundaries.front().conditions.back().kind
			= iga::FieldBoundaryKind::Dirichlet;
		RequireRejected([&lifecycle, &pressure_dirichlet_configuration] {
			lifecycle.SetTrialBoundaryConfiguration(pressure_dirichlet_configuration);
		}, "adding a pressure constraint changes fixed PETSc boundary rows");
		auto no_pressure_constraint_configuration = trial_flow_configuration;
		no_pressure_constraint_configuration.boundaries.back().conditions.clear();
		RequireRejected([&lifecycle, &no_pressure_constraint_configuration] {
			lifecycle.SetTrialBoundaryConfiguration(no_pressure_constraint_configuration);
		}, "removing a pressure constraint changes fixed PETSc boundary rows");
		lifecycle.SetPortInput(pressure_port, pressure_input);
		assert(lifecycle.PressureTractionValue(1).has_value());
		RequireNear(4.0, *lifecycle.PressureTractionValue(1), 0.0,
			"trial pressure input");
		lifecycle.SolveTrial();
		const auto first_trial = CopyVector(lifecycle.State());
		const auto first_trial_iterations = lifecycle.TrialLinearIterations();
		const auto trial_port_state = lifecycle.GetPortState(pressure_port);
		RequireNear(0.1, trial_port_state.time_s, 0.0, "trial port state time");
		assert(first_trial != committed_state);
		assert(lifecycle.Summary().linear_iterations == 0);
		RequireRejected([&lifecycle] { lifecycle.SolveTrial(); },
			"SolveTrial without rollback");
		lifecycle.RollbackTrial();
		assert(CopyVector(lifecycle.State()) == committed_state);
		assert(lifecycle.PressureTractionValue(1).has_value());
		RequireNear(0.0, *lifecycle.PressureTractionValue(1), 0.0,
			"committed pressure traction rollback");
		RequireRejected([&lifecycle] { lifecycle.CommitStep(); },
			"CommitStep after rollback");
		lifecycle.SetTrialBoundaryConfiguration(trial_flow_configuration);
		auto normal_traction_port = pressure_port;
		normal_traction_port.id = "normal_traction_input";
		normal_traction_port.requires = {iga::PortQuantity::MeanNormalTraction};
		iga::PortBoundaryData normal_traction_input;
		normal_traction_input.time_s = 0.1;
		normal_traction_input.mean_normal_traction_pa = -4.0;
		lifecycle.SetPortInput(normal_traction_port, normal_traction_input);
		RequireNear(4.0, *lifecycle.PressureTractionValue(1), 0.0,
			"outward normal traction to pressure conversion");
		lifecycle.SolveTrial();
		const auto second_trial = CopyVector(lifecycle.State());
		assert(second_trial == first_trial);
		assert(lifecycle.TrialLinearIterations() == first_trial_iterations);
		lifecycle.PrepareCommitStep();
		assert(lifecycle.Phase() == iga::FlowStepPhase::CommitPrepared);
		lifecycle.FinalizeCommitStep();
		assert(lifecycle.Phase() == iga::FlowStepPhase::Committed);
		assert(CopyVector(lifecycle.State()) == first_trial);
		assert(lifecycle.Summary().linear_iterations == first_trial_iterations);
		RequireRejected([&lifecycle] { lifecycle.CommitStep(); },
			"duplicate CommitStep");
		RequireRejected([&lifecycle, &pressure_port] {
			lifecycle.GetPortState(pressure_port);
		}, "GetPortState after commit");
		iga::TransientFlowRuntime advance_adapter(database, PETSC_COMM_WORLD, true, true,
			{1.0, 1.0, 0.1}, MakeLifecycleFlowBoundaries(lifecycle_labels),
			lifecycle_labels,
			std::vector<std::array<double, 3>>(64, {0.0, 0.0, 0.0}), {}, {});
		advance_adapter.InitializeState(initial_flow_configuration);
		advance_adapter.Advance(MakeFlowConfiguration(0.125, 4.0), 0, 0.1,
			12, 1e-8, 1e-8, 3.0);
		assert(CopyVector(advance_adapter.State()) == first_trial);
		assert(advance_adapter.Summary().linear_iterations == first_trial_iterations);
		iga::PortBoundaryData unsupported_flow_input;
		unsupported_flow_input.time_s = 0.2;
		unsupported_flow_input.outward_flow_m3_s = 1.0;
		lifecycle.BeginStep(1, 0.2, 12, 1e-8, 1e-8, 3.0);
		lifecycle.SetTrialBoundaryConfiguration(trial_flow_configuration);
		auto abort_pressure = pressure_input;
		abort_pressure.time_s = 0.2;
		abort_pressure.mean_pressure_pa = 9.0;
		lifecycle.SetPortInput(pressure_port, abort_pressure);
		RequireRejected([&lifecycle, &pressure_port, &unsupported_flow_input] {
			lifecycle.SetPortInput(pressure_port, unsupported_flow_input);
		}, "unsupported flow-profile trial input");
		lifecycle.AbortStep();
		assert(lifecycle.Phase() == iga::FlowStepPhase::Committed);
		assert(CopyVector(lifecycle.State()) == first_trial);
		RequireNear(4.0, *lifecycle.PressureTractionValue(1), 0.0,
			"committed pressure traction abort");
		lifecycle.BeginStep(1, 0.2, 12, 1e-8, 1e-8, 3.0);
		lifecycle.SetTrialBoundaryConfiguration(trial_flow_configuration);
		auto second_step_pressure = pressure_input;
		second_step_pressure.time_s = 0.2;
		lifecycle.SetPortInput(pressure_port, second_step_pressure);
		lifecycle.SolveTrial();
		lifecycle.AbortStep();
		assert(lifecycle.Phase() == iga::FlowStepPhase::Committed);
		assert(CopyVector(lifecycle.State()) == first_trial);
		assert(lifecycle.TrialLinearIterations() == 0);
		lifecycle.BeginStep(1, 0.2, 12, 1e-8, 1e-8, 3.0);
		lifecycle.SetTrialBoundaryConfiguration(trial_flow_configuration);
		lifecycle.SetPortInput(pressure_port, second_step_pressure);
		lifecycle.SolveTrial();
		lifecycle.PrepareCommitStep();
		lifecycle.AbortStep();
		assert(lifecycle.Phase() == iga::FlowStepPhase::Committed);
		assert(CopyVector(lifecycle.State()) == first_trial);
		assert(lifecycle.Summary().linear_iterations == first_trial_iterations);
		assert(lifecycle.TrialLinearIterations() == 0);
		lifecycle.AbortStep();
		{
			auto profile_configuration = MakeFlowConfiguration(0.0, 0.0);
			auto& profile = profile_configuration.boundaries.front().conditions.front();
			profile.value.clear();
			profile.profile = "initial_velocityfield.txt";
			profile.scale = 1.0;
			const std::vector<std::array<double, 3>> reference_velocity(
				64, {0.0, 0.0, 0.125});
			iga::TransientFlowRuntime adapter_native(database, PETSC_COMM_WORLD, true, true,
				{1.0, 1.0, 0.1},
				iga::ResolveFlowBoundaries(profile_configuration,
					profile_configuration.equation_systems.front(), lifecycle_labels,
					reference_velocity),
				lifecycle_labels, reference_velocity, {}, {});
			adapter_native.InitializeState(profile_configuration);
			const double reference_flow = adapter_native.ReferenceBoundaryFlow(1);
			assert(reference_flow != 0.0 && std::isfinite(reference_flow));
			iga::CouplingPort inlet;
			inlet.id = "inlet";
			inlet.subsystem_id = "adapter_three_d";
			inlet.locator_kind = "boundary_label";
			inlet.locator = "1";
			inlet.provides = {iga::PortQuantity::Area, iga::PortQuantity::FlowRate,
				iga::PortQuantity::MeanPressure};
			inlet.requires = {iga::PortQuantity::FlowRate};
			iga::ThreeDBodyFittedFlowDomainAdapter adapter("adapter_three_d",
				adapter_native, {inlet}, profile_configuration,
				std::filesystem::temp_directory_path(), {{"inlet", reference_flow}});
			adapter.BeginStep({0, 0.0, 0.1});
			iga::PortBoundaryData flow_input;
			flow_input.time_s = 0.1;
			flow_input.outward_flow_m3_s = reference_flow;
			adapter.SetPortInput("inlet", flow_input);
			adapter.SolveTrial();
			const auto adapter_state = adapter.GetPortState("inlet");
			assert(adapter_state.outward_flow_m3_s.has_value());
			adapter.PrepareCommitStep();
			adapter.FinalizeCommitStep();
			assert(adapter_native.Phase() == iga::FlowStepPhase::Committed);
		}

		auto outlet_model = iga::OutletModelState{};
		outlet_model.label = 1;
		outlet_model.kind = iga::FieldBoundaryKind::WindkesselRC;
		outlet_model.resistance = 1.0e-3;
		outlet_model.capacitance = 1.0;
		outlet_model.reference_pressure = 0.0;
		outlet_model.capacitor_pressure = 0.0;
		outlet_model.pressure = iga::EvaluateOutletModel(outlet_model, 0.0, 0.1).pressure;
		iga::TransientFlowRuntime outlet_lifecycle(database, PETSC_COMM_WORLD, true, true,
			{1.0, 1.0, 0.1}, MakeLifecycleFlowBoundaries(lifecycle_labels),
			lifecycle_labels,
			std::vector<std::array<double, 3>>(64, {0.0, 0.0, 0.0}), {}, {outlet_model});
		const auto outlet_configuration = MakeOutletFlowConfiguration();
		outlet_lifecycle.InitializeState(outlet_configuration);
		outlet_lifecycle.BeginStep(0, 0.1, 12, 1e-8, 1e-8, 3.0);
		outlet_lifecycle.SetTrialBoundaryConfiguration(outlet_configuration);
		RequireRejected([&outlet_lifecycle, &pressure_port, &pressure_input] {
			outlet_lifecycle.SetPortInput(pressure_port, pressure_input);
		}, "pressure override conflicts with outlet model");
		outlet_lifecycle.SolveTrial();
		const auto first_trial_capacitor =
			OutletCapacitorPressure(outlet_lifecycle);
		assert(first_trial_capacitor != outlet_model.capacitor_pressure);
		outlet_lifecycle.RollbackTrial();
		RequireNear(outlet_model.capacitor_pressure,
			OutletCapacitorPressure(outlet_lifecycle), 0.0,
			"outlet capacitor rollback");
		outlet_lifecycle.SetTrialBoundaryConfiguration(outlet_configuration);
		outlet_lifecycle.SolveTrial();
		RequireNear(first_trial_capacitor,
			OutletCapacitorPressure(outlet_lifecycle), 0.0,
			"outlet capacitor deterministic replay");
		outlet_lifecycle.CommitStep();
		RequireNear(first_trial_capacitor,
			OutletCapacitorPressure(outlet_lifecycle), 0.0,
			"outlet capacitor single commit");
		auto initial_configuration = MakeConfiguration(0.0);
		auto system = iga::CompileLinearSystem(initial_configuration, "oxygen_transport");
		system.velocity_source = "prescribed";
		iga::TransientTransportRuntime runtime(database, PETSC_COMM_WORLD,
			initial_configuration, system, std::vector<int>(64, 1));
		auto step_configuration = MakeConfiguration(2.0);
		assert(runtime.RequiredNodes().size() == 64);
		const auto transport_velocity
			= std::vector<std::array<double, 3>>(64, {0.0, 0.0, 0.0});
		const auto initial_transport_state = runtime.GatherState();
		runtime.BeginStep();
		RequireRejected([&runtime] { runtime.BeginStep(); },
			"nested transport BeginStep");
		runtime.SolveTrial(step_configuration, runtime.RequiredNodes(), transport_velocity);
		const auto first_transport_trial = runtime.GatherState();
		assert(first_transport_trial != initial_transport_state);
		assert(runtime.Steps() == 1);
		runtime.RollbackTrial();
		assert(runtime.GatherState() == initial_transport_state);
		assert(runtime.Steps() == 0);
		runtime.SolveTrial(step_configuration, runtime.RequiredNodes(), transport_velocity);
		assert(runtime.GatherState() == first_transport_trial);
		runtime.PrepareCommitStep();
		assert(runtime.Phase() == iga::TransportStepPhase::CommitPrepared);
		runtime.FinalizeCommitStep();
		assert(runtime.Phase() == iga::TransportStepPhase::Committed);
		assert(runtime.Steps() == 1);
		runtime.AbortStep();
		assert(runtime.GatherState() == first_transport_trial);
		runtime.BeginStep();
		RequireRejected([&runtime] { runtime.PrepareCommitStep(); },
			"transport prepare before solve");
		runtime.SolveTrial(step_configuration, runtime.RequiredNodes(), transport_velocity);
		runtime.PrepareCommitStep();
		runtime.AbortStep();
		assert(runtime.Phase() == iga::TransportStepPhase::Committed);
		assert(runtime.GatherState() == first_transport_trial);
		assert(runtime.Steps() == 1);
		iga::TransientTransportRuntime advance_transport(database, PETSC_COMM_WORLD,
			initial_configuration, system, std::vector<int>(64, 1));
		advance_transport.Advance(step_configuration, advance_transport.RequiredNodes(),
			transport_velocity);
		assert(advance_transport.GatherState() == first_transport_trial);
		runtime.BeginStep();
		runtime.AbortStep();
		assert(runtime.GatherState() == first_transport_trial);
		assert(runtime.Steps() == 1);
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
