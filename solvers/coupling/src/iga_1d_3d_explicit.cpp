#include "BoundarySupport.hpp"
#include "ExplicitOneDThreeDCoupling.hpp"
#include "IgaDatabase.hpp"
#include "OneDImplicit.hpp"
#include "OneDRuntime.hpp"
#include "ThreeDFlowCoupling.hpp"
#include "ThreeDVcaCoupling.hpp"
#include "TransientFlowRuntime.hpp"

#include <petscsys.h>

#include <algorithm>
#include <cmath>
#include <array>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr int kThreeDMaximumNewtonIterations = 30;
constexpr double kThreeDNonlinearRelativeTolerance = 1.0e-5;
constexpr double kThreeDNonlinearAbsoluteTolerance = 1.0e-10;
constexpr double kThreeDMassRelativeTolerance = 1.0e-3;

struct Options {
	fs::path database;
	fs::path three_d_case;
	fs::path upstream_case;
	fs::path downstream_case;
	int upstream_terminal_node = -1;
	fs::path output_directory;
	int stop_after_step = 0;
};

int PositiveInteger(const std::string& text, const std::string& option)
{
	std::size_t used = 0;
	int value = 0;
	try { value = std::stoi(text, &used); }
	catch (const std::exception&) { throw std::runtime_error(option+" requires a positive integer"); }
	if (used != text.size() || value < 1) throw std::runtime_error(option+" requires a positive integer");
	return value;
}

// Test-only failure injection.  It is deliberately read once before any
// runtime lifecycle transition and has no effect unless the environment is set.
int ExplicitCouplingFailureInjectionStep()
{
	const char* value = std::getenv("TUBULARFLOWIGA_INJECT_EXPLICIT_COUPLING_FAILURE_STEP");
	return value ? PositiveInteger(value, "TUBULARFLOWIGA_INJECT_EXPLICIT_COUPLING_FAILURE_STEP") : 0;
}

Options ParseOptions(int argc, char** argv)
{
	if (argc < 5) throw std::runtime_error(
		"usage: iga_1d_3d_explicit DB THREE_D_CASE UPSTREAM_1D_CASE DOWNSTREAM_1D_CASE "
		"--upstream-terminal-node ID --output-dir DIR [--stop-after-step N] [PETSc options]");
	Options options;
	options.database = argv[1];
	options.three_d_case = argv[2];
	options.upstream_case = argv[3];
	options.downstream_case = argv[4];
	for (int i = 5; i < argc; ++i) {
		const std::string argument(argv[i]);
		if (argument == "--upstream-terminal-node" || argument == "--output-dir"
			|| argument == "--stop-after-step") {
			if (++i >= argc) throw std::runtime_error(argument+" requires a value");
			const std::string value(argv[i]);
			if (argument == "--upstream-terminal-node") options.upstream_terminal_node = PositiveInteger(value, argument);
			else if (argument == "--output-dir") options.output_directory = value;
			else options.stop_after_step = PositiveInteger(value, argument);
			continue;
		}
		if (iga::IsExplicitCouplingPetscOption(argument)) {
			if (i+1 < argc && iga::ExplicitCouplingPetscOptionConsumesNextValue(argument, argv[i+1])) ++i;
			continue;
		}
		throw std::runtime_error("unexpected argument: "+argument);
	}
	if (options.upstream_terminal_node < 0 || options.output_directory.empty())
		throw std::runtime_error("--upstream-terminal-node and --output-dir are required");
	return options;
}

std::string ReadText(const fs::path& path)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot open file: "+path.string());
	std::ostringstream text;
	text << input.rdbuf();
	return text.str();
}

const iga::OneDFlowSystemDefinition& SelectFlow(const iga::OneDConfiguration& configuration)
{
	if (configuration.flow_systems.size() != 1)
		throw std::runtime_error("explicit 1D--3D coupling requires exactly one 1D flow system per case");
	return configuration.flow_systems.front();
}

void RequireOneDFlowOnly(const iga::OneDConfiguration& configuration, const char* name)
{
	if (configuration.coupling.mode != iga::SimulationScopeMode::FlowOnly
		|| !configuration.transport_systems.empty() || configuration.physiology.enabled)
		throw std::runtime_error(std::string("explicit 1D--3D coupling requires flow_only 1D ")+name
			+" case without transport or physiology");
}

iga::CouplingPort MakeThreeDPort(const std::string& id, int label,
	std::set<iga::PortQuantity> provides, std::set<iga::PortQuantity> requires)
{
	iga::CouplingPort port;
	port.id = id;
	port.subsystem_id = "three_d";
	port.locator_kind = "boundary_label";
	port.locator = std::to_string(label);
	port.provides = std::move(provides);
	port.requires = std::move(requires);
	iga::ValidateCouplingPort(port);
	return port;
}

double RequirePortValue(const std::optional<double>& value, const char* name)
{
	if (!value) throw std::runtime_error(std::string("explicit coupling port does not provide ")+name);
	iga::RequireFinitePortValue(name, *value);
	return *value;
}

void ApplyInitialPressure(iga::SimulationConfiguration& configuration, const std::string& pressure_field,
	int label, double pressure)
{
	iga::RequireFinitePortValue("initial 3D outlet pressure", pressure);
	iga::FieldBoundaryCondition* match = nullptr;
	for (auto& boundary : configuration.boundaries) {
		if (boundary.label != label) continue;
		for (auto& condition : boundary.conditions) {
			if (condition.field != pressure_field) continue;
			if (match || condition.kind != iga::FieldBoundaryKind::PressureTraction)
				throw std::runtime_error("explicit coupling 3D outlet requires exactly one pressure_traction declaration");
			match = &condition;
		}
	}
	if (!match) throw std::runtime_error("explicit coupling 3D outlet has no pressure_traction declaration");
	match->value = {pressure};
	match->waveform.clear();
}

std::string JsonEscape(const std::string& text)
{
	std::ostringstream escaped;
	for (const unsigned char character : text) {
		switch (character) {
		case '"': escaped << "\\\""; break;
		case '\\': escaped << "\\\\"; break;
		case '\b': escaped << "\\b"; break;
		case '\f': escaped << "\\f"; break;
		case '\n': escaped << "\\n"; break;
		case '\r': escaped << "\\r"; break;
		case '\t': escaped << "\\t"; break;
		default:
			if (character < 0x20) {
				escaped << "\\u" << std::hex << std::setw(4) << std::setfill('0')
					<< static_cast<unsigned int>(character) << std::dec << std::setfill(' ');
			} else escaped << character;
		}
	}
	return escaped.str();
}

void WriteExplicitCouplingManifest(const fs::path& path, int upstream_terminal_node,
	int three_d_inlet_label, int three_d_outlet_label, double dt_s, int configured_steps,
	int completed_steps,
	double density_kg_m3, double dynamic_viscosity_pa_s, double normalized_length_m,
	double reference_inlet_outward_flow_m3_s, double initial_three_d_inlet_pressure_pa,
	double initial_downstream_root_pressure_pa)
{
	std::ofstream output(path);
	if (!output) throw std::runtime_error("cannot create explicit coupling manifest");
	output << std::setprecision(17) << "{\n"
		<< "  \"scheme\": \"" << JsonEscape("explicit_staggered") << "\",\n"
		<< "  \"static_pressure_definition\": \"" << JsonEscape(
			"mean static pressure in Pa, applied through the existing 3D pressure_traction "
			"boundary condition and lagged from the downstream 1D root") << "\",\n"
		<< "  \"ports\": [\n"
		<< "    {\"id\": \"upstream_terminal\", \"locator\": \"outlet:"
		<< upstream_terminal_node << "\", \"native_to_outward_sign\": 1},\n"
		<< "    {\"id\": \"three_d_inlet\", \"locator_kind\": \"boundary_label\", \"locator\": \""
		<< three_d_inlet_label << "\", \"native_to_outward_sign\": 1},\n"
		<< "    {\"id\": \"three_d_outlet\", \"locator_kind\": \"boundary_label\", \"locator\": \""
		<< three_d_outlet_label << "\", \"native_to_outward_sign\": 1},\n"
		<< "    {\"id\": \"downstream_root\", \"locator\": \"root\", \"native_to_outward_sign\": -1}\n"
		<< "  ],\n"
		<< "  \"dt_s\": " << dt_s << ",\n"
		<< "  \"steps\": " << configured_steps << ",\n"
		<< "  \"completed_steps\": " << completed_steps << ",\n"
		<< "  \"density_kg_m3\": " << density_kg_m3 << ",\n"
		<< "  \"dynamic_viscosity_pa_s\": " << dynamic_viscosity_pa_s << ",\n"
		<< "  \"geometry_transform_product_m\": " << normalized_length_m << ",\n"
		<< "  \"reference_inlet_outward_flow_m3_s\": " << reference_inlet_outward_flow_m3_s << ",\n"
		<< "  \"newton_controls\": {\"maximum_iterations\": " << kThreeDMaximumNewtonIterations
		<< ", \"nonlinear_relative_tolerance\": " << kThreeDNonlinearRelativeTolerance
		<< ", \"nonlinear_absolute_tolerance\": " << kThreeDNonlinearAbsoluteTolerance
		<< ", \"mass_relative_tolerance\": " << kThreeDMassRelativeTolerance << "},\n"
		<< "  \"initial_lagged_pressures_pa\": {\"three_d_inlet\": "
		<< initial_three_d_inlet_pressure_pa << ", \"downstream_root\": "
		<< initial_downstream_root_pressure_pa << "}\n"
		<< "}\n";
	if (!output) throw std::runtime_error("cannot write explicit coupling manifest");
}

void RollbackSolved(iga::OneDFlowRuntime& runtime)
{
	if (runtime.CurrentPhase() == iga::OneDFlowRuntime::Phase::TrialSolved) runtime.RollbackTrial();
}

void RollbackSolved(iga::TransientFlowRuntime& runtime)
{
	if (runtime.Phase() == iga::FlowStepPhase::TrialSolved) runtime.RollbackTrial();
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, "TubularFlowIGA explicit 1D--3D flow-only coupling\n");
	int rank = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	int status = 0;
	try {
		const auto options = ParseOptions(argc, argv);
		const int injected_failure_step = ExplicitCouplingFailureInjectionStep();
		auto upstream_configuration = iga::ParseOneDConfiguration(
			ReadText(options.upstream_case/"simulation_config.json"));
		auto downstream_configuration = iga::ParseOneDConfiguration(
			ReadText(options.downstream_case/"simulation_config.json"));
		RequireOneDFlowOnly(upstream_configuration, "upstream");
		RequireOneDFlowOnly(downstream_configuration, "downstream");
		const auto upstream_flow = SelectFlow(upstream_configuration);
		const auto downstream_flow = SelectFlow(downstream_configuration);
		auto upstream_network = iga::ReadOneDNetwork(options.upstream_case/upstream_configuration.geometry.file,
			upstream_configuration.geometry.length_scale_to_m, upstream_flow.discretization.cells_per_segment,
			upstream_flow.dynamic_viscosity, upstream_configuration.geometry.root_node_id);
		auto downstream_network = iga::ReadOneDNetwork(options.downstream_case/downstream_configuration.geometry.file,
			downstream_configuration.geometry.length_scale_to_m, downstream_flow.discretization.cells_per_segment,
			downstream_flow.dynamic_viscosity, downstream_configuration.geometry.root_node_id);
		iga::ValidateOneDTopologyReferences(upstream_configuration, upstream_network);
		iga::ValidateOneDTopologyReferences(downstream_configuration, downstream_network);
		const auto upstream_inlet = iga::ResolveOneDInlet(upstream_configuration);
		const auto downstream_inlet = iga::ResolveOneDInlet(downstream_configuration);
		iga::OneDFlowRuntime upstream(upstream_configuration, upstream_flow, std::move(upstream_network),
			upstream_inlet, options.upstream_case, [](const iga::OneDNetwork& network,
				const iga::OneDFlowSystemDefinition& flow, iga::OneDFlowState& state, double inlet, double dt) {
				iga::AdvanceImplicitOneD(network, flow, state, inlet, dt);
			});
		iga::OneDFlowRuntime downstream(downstream_configuration, downstream_flow, std::move(downstream_network),
			downstream_inlet, options.downstream_case, [](const iga::OneDNetwork& network,
				const iga::OneDFlowSystemDefinition& flow, iga::OneDFlowState& state, double inlet, double dt) {
				iga::AdvanceImplicitOneD(network, flow, state, inlet, dt);
			});
		if (upstream.FlowState().outlets.size() != 1)
			throw std::runtime_error("explicit 1D--3D coupling requires exactly one upstream terminal closure");
		const auto upstream_terminal = "outlet:"+std::to_string(options.upstream_terminal_node);
		const auto upstream_node = upstream.Network().node_index.find(options.upstream_terminal_node);
		if (upstream_node == upstream.Network().node_index.end())
			throw std::runtime_error("unknown upstream terminal node");
		const auto& upstream_outlet = upstream.FlowState().outlets.front();
		if (upstream_outlet.node != upstream_node->second)
			throw std::runtime_error("--upstream-terminal-node must match the only upstream terminal closure");
		if (upstream_outlet.kind != iga::OneDOutletKind::Pressure)
			throw std::runtime_error("selected upstream terminal requires a pressure closure");
		if (downstream.FlowState().outlets.size() != 1)
			throw std::runtime_error("explicit 1D--3D coupling requires exactly one downstream terminal closure");
		const auto downstream_terminal = "outlet:"+std::to_string(
			downstream.Network().nodes[static_cast<std::size_t>(
				downstream.FlowState().outlets.front().node)].id);

		iga::Database database(options.database.string());
		auto three_d_configuration = iga::ReadSimulationConfiguration(
			(options.three_d_case/"simulation_config.json").string());
		if (three_d_configuration.equation_systems.size() != 1
			|| three_d_configuration.equation_systems.front().kind != iga::EquationKind::NavierStokes
			|| three_d_configuration.equation_systems.front().unknowns.size() != 2)
			throw std::runtime_error("explicit 1D--3D coupling requires exactly one two-unknown Navier-Stokes 3D system");
		for (const auto& field : three_d_configuration.fields)
			if (field.kind == iga::FieldKind::Scalar)
				throw std::runtime_error("explicit 1D--3D coupling rejects 3D transport fields");
		const auto& three_d_flow = three_d_configuration.equation_systems.front();
		const auto& three_d_pressure_unknown = three_d_flow.unknowns.at(1);
		if (three_d_configuration.coupling.mode != iga::SimulationScopeMode::FlowOnly
			|| three_d_configuration.physiology.enabled
			|| three_d_flow.time_integration != "backward_euler")
			throw std::runtime_error("explicit 1D--3D coupling requires flow_only backward_euler 3D flow without physiology");
		iga::RequireThreeDVascularPorts(three_d_configuration.coupling, "explicit 1D--3D coupling");
		const auto& ports = three_d_configuration.coupling.three_d_ports;
		if (ports.outlet_labels.size() != 1)
			throw std::runtime_error("explicit 1D--3D coupling requires exactly one declared 3D outlet");
		if (ports.inlet_label == 0 || ports.outlet_labels.front() == 0)
			throw std::runtime_error(
				"explicit 1D--3D coupling reserves boundary label 0 for the no-slip wall");
		const auto mesh = iga::ReadLabeledHexMesh((options.three_d_case/"controlmesh.vtk").string(),
			database.header().nodes, database.header().elements);
		const auto boundary_velocity = iga::ReadVelocity(
			(options.three_d_case/"initial_velocityfield.txt").string(), database.header().nodes);
		const auto wall_trace_basis = iga::WallTraceBasis(database, mesh, 0);
		auto initial_three_d = iga::MaterializeBoundaryWaveforms(three_d_configuration,
			options.three_d_case.string(), 0.0);
		auto outlet_models = iga::InitializeOutletModels(three_d_configuration, three_d_flow);
		if (!outlet_models.empty()) throw std::runtime_error("explicit 1D--3D coupling requires a static-pressure 3D outlet");
		const auto initial_boundaries = iga::ResolveFlowBoundaries(initial_three_d,
			initial_three_d.equation_systems.front(), mesh.labels, boundary_velocity);
		iga::TransientFlowRuntime three_d(database, PETSC_COMM_WORLD, true, true,
			{three_d_flow.density, three_d_flow.viscosity, three_d_configuration.time.dt},
			initial_boundaries, mesh.labels, boundary_velocity, wall_trace_basis, std::move(outlet_models));
		iga::RequireValidGeometry(three_d.Elements(), rank, PETSC_COMM_WORLD);
		const std::array<int, 3> required_three_d_labels{{0, ports.inlet_label,
			ports.outlet_labels.front()}};
		std::array<long long, 3> local_face_counts{{0, 0, 0}};
		int local_unknown_label = 0;
		for (const auto& element : three_d.OwnedElements())
			for (const int face_label : element.boundary_labels) {
				if (face_label < 0) continue;
				auto found = std::find(required_three_d_labels.begin(), required_three_d_labels.end(), face_label);
				if (found == required_three_d_labels.end()) {
					local_unknown_label = 1;
					continue;
				}
				++local_face_counts[static_cast<std::size_t>(found-required_three_d_labels.begin())];
			}
		int global_unknown_label = 0;
		MPI_Allreduce(&local_unknown_label, &global_unknown_label, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
		if (global_unknown_label)
			throw std::runtime_error("explicit coupling 3D database has a non-port non-wall boundary label");
		std::array<long long, 3> global_face_counts{{0, 0, 0}};
		MPI_Allreduce(local_face_counts.data(), global_face_counts.data(), 3, MPI_LONG_LONG, MPI_SUM,
			PETSC_COMM_WORLD);
		for (const auto count : global_face_counts)
			if (count == 0)
				throw std::runtime_error("explicit coupling 3D requires wall, inlet, and outlet .ntiga boundary faces");
		const double normalized_length = database.header().version == iga::kVersion
			? database.header().geometry_transform.source_units_per_normalized_unit
				*database.header().geometry_transform.source_length_scale_to_m : 0.0;
		const double reference_inlet_flow = three_d.ReferenceBoundaryFlow(ports.inlet_label);
		iga::ExplicitCouplingScalarPreflight scalar;
		scalar.dt_s = three_d_configuration.time.dt;
		scalar.steps = three_d_configuration.time.steps;
		scalar.density_kg_m3 = three_d_flow.density;
		scalar.dynamic_viscosity_pa_s = three_d_flow.viscosity;
		scalar.normalized_length_m = normalized_length;
		scalar.reference_inlet_outward_flow_m3_s = reference_inlet_flow;
		iga::ValidateExplicitCouplingScalarPreflight(scalar);
		if (upstream.Configuration().time.dt != scalar.dt_s || downstream.Configuration().time.dt != scalar.dt_s
			|| upstream.Configuration().time.steps != scalar.steps || downstream.Configuration().time.steps != scalar.steps
			|| upstream.FlowSystem().density != scalar.density_kg_m3
			|| downstream.FlowSystem().density != scalar.density_kg_m3
			|| upstream.FlowSystem().dynamic_viscosity != scalar.dynamic_viscosity_pa_s
			|| downstream.FlowSystem().dynamic_viscosity != scalar.dynamic_viscosity_pa_s)
			throw std::runtime_error("explicit 1D--3D coupling requires identical dt, steps, density, and viscosity");
		const auto three_d_inlet_profile = MakeThreeDPort("three_d_inlet", ports.inlet_label, {},
			{iga::PortQuantity::FlowRate});
		const auto three_d_inlet_measure = MakeThreeDPort("three_d_inlet", ports.inlet_label,
			{iga::PortQuantity::Area, iga::PortQuantity::FlowRate, iga::PortQuantity::MeanPressure}, {});
		const auto three_d_outlet_measure = MakeThreeDPort("three_d_outlet", ports.outlet_labels.front(),
			{iga::PortQuantity::Area, iga::PortQuantity::FlowRate, iga::PortQuantity::MeanPressure}, {});
		const auto three_d_wall_measure = MakeThreeDPort("three_d_wall", 0,
			{iga::PortQuantity::FlowRate}, {});
		const auto three_d_outlet_pressure = MakeThreeDPort("three_d_outlet_pressure", ports.outlet_labels.front(), {},
			{iga::PortQuantity::MeanPressure});

		const double initial_upstream_flow = iga::EvaluateOneDInlet(upstream.Configuration(), upstream.InletDefinition(),
			options.upstream_case, 0.0, upstream.Network().segments.front().area0);
		if (!(initial_upstream_flow > 0.0) || !std::isfinite(initial_upstream_flow))
			throw std::runtime_error("explicit 1D--3D coupling requires a positive initial upstream native flow");
		upstream.InitializeOpenLoop(initial_upstream_flow);
		downstream.InitializeOpenLoop(initial_upstream_flow);
		const auto initial_downstream_root = downstream.GetPortState("root");
		iga::PortBoundaryData initial_profile_input;
		initial_profile_input.time_s = 0.0;
		initial_profile_input.outward_flow_m3_s = -initial_upstream_flow;
		iga::ApplyThreeDReferenceProfileInput(initial_three_d,
			initial_three_d.equation_systems.front(), three_d_inlet_profile,
			initial_profile_input, reference_inlet_flow);
		ApplyInitialPressure(initial_three_d, three_d_pressure_unknown, ports.outlet_labels.front(),
			RequirePortValue(initial_downstream_root.mean_pressure_pa, "downstream root pressure"));
		three_d.InitializeState(initial_three_d);
		const auto initial_three_d_ports = three_d.MeasurePorts({three_d_inlet_measure, three_d_outlet_measure},
			0.0, {}, {});
		double lagged_three_d_inlet_pressure = RequirePortValue(
			initial_three_d_ports.at("three_d_inlet").mean_pressure_pa, "3D inlet pressure");
		double lagged_downstream_root_pressure = RequirePortValue(
			initial_downstream_root.mean_pressure_pa, "downstream root pressure");
		const double initial_lagged_three_d_inlet_pressure = lagged_three_d_inlet_pressure;
		const double initial_lagged_downstream_root_pressure = lagged_downstream_root_pressure;
		const int final_step = options.stop_after_step > 0 ? options.stop_after_step : scalar.steps;
		if (final_step > scalar.steps) throw std::runtime_error("--stop-after-step exceeds configured steps");
		std::vector<iga::ExplicitCouplingHistoryRow> history;
		history.reserve(static_cast<std::size_t>(final_step));
		for (int step = 1; step <= final_step; ++step) {
			const double time = step*scalar.dt_s;
			try {
				const double upstream_flow = iga::EvaluateOneDInlet(upstream.Configuration(), upstream.InletDefinition(),
					options.upstream_case, time, upstream.Network().segments.front().area0);
				upstream.BeginStep(upstream.FlowState().physical_time, scalar.dt_s);
				upstream.SetOpenLoopInlet(upstream.OpenLoopInlet(time, upstream_flow));
				iga::PortBoundaryData upstream_pressure;
				upstream_pressure.time_s = time;
				upstream_pressure.mean_pressure_pa = lagged_three_d_inlet_pressure;
				upstream.SetPortInput(upstream_terminal, upstream_pressure);
				upstream.SolveTrial();
				const auto upstream_port = upstream.GetPortState(upstream_terminal);
				const double upstream_q = RequirePortValue(upstream_port.outward_flow_m3_s, "upstream terminal flow");

				auto three_d_step = iga::MaterializeBoundaryWaveforms(three_d_configuration,
					options.three_d_case.string(), time);
				iga::PortBoundaryData three_d_profile_input;
				three_d_profile_input.time_s = time;
				three_d_profile_input.outward_flow_m3_s = -upstream_q;
				iga::ApplyThreeDReferenceProfileInput(three_d_step,
					three_d_step.equation_systems.front(), three_d_inlet_profile,
					three_d_profile_input, reference_inlet_flow);
				three_d.BeginStep(step-1, time, kThreeDMaximumNewtonIterations,
					kThreeDNonlinearRelativeTolerance, kThreeDNonlinearAbsoluteTolerance,
					kThreeDMassRelativeTolerance);
				three_d.SetTrialBoundaryConfiguration(three_d_step);
				iga::PortBoundaryData three_d_pressure_input;
				three_d_pressure_input.time_s = time;
				three_d_pressure_input.mean_pressure_pa = lagged_downstream_root_pressure;
				three_d.SetPortInput(three_d_outlet_pressure, three_d_pressure_input);
				three_d.SolveTrial();
				const auto three_d_inlet = three_d.GetPortState(three_d_inlet_measure);
				const auto three_d_outlet = three_d.GetPortState(three_d_outlet_measure);
				const auto three_d_wall = three_d.GetPortState(three_d_wall_measure);
				const double three_d_outlet_q = RequirePortValue(three_d_outlet.outward_flow_m3_s, "3D outlet flow");

				downstream.BeginStep(downstream.FlowState().physical_time, scalar.dt_s);
				iga::PortBoundaryData downstream_input;
				downstream_input.time_s = time;
				downstream_input.outward_flow_m3_s = -three_d_outlet_q;
				downstream.SetPortInput("root", downstream_input);
				downstream.SolveTrial();
				const auto downstream_root = downstream.GetPortState("root");
				const auto downstream_terminal_port = downstream.GetPortState(downstream_terminal);
				const auto upstream_root = upstream.GetPortState("root");

				iga::ExplicitCouplingHistoryRow row;
				row.time_s = time;
				row.upstream_root_pressure_pa = RequirePortValue(upstream_root.mean_pressure_pa, "upstream root pressure");
				row.upstream_root_outward_flow_m3_s = RequirePortValue(upstream_root.outward_flow_m3_s, "upstream root flow");
				row.upstream_root_area_m2 = RequirePortValue(upstream_root.area_m2, "upstream root area");
				row.upstream_terminal_pressure_pa = RequirePortValue(upstream_port.mean_pressure_pa, "upstream terminal pressure");
				row.upstream_terminal_outward_flow_m3_s = upstream_q;
				row.upstream_terminal_area_m2 = RequirePortValue(upstream_port.area_m2, "upstream terminal area");
				row.three_d_inlet_pressure_pa = RequirePortValue(three_d_inlet.mean_pressure_pa, "3D inlet pressure");
				row.three_d_inlet_outward_flow_m3_s = RequirePortValue(three_d_inlet.outward_flow_m3_s, "3D inlet flow");
				row.three_d_inlet_area_m2 = RequirePortValue(three_d_inlet.area_m2, "3D inlet area");
				row.three_d_outlet_pressure_pa = RequirePortValue(three_d_outlet.mean_pressure_pa, "3D outlet pressure");
				row.three_d_outlet_outward_flow_m3_s = three_d_outlet_q;
				row.three_d_outlet_area_m2 = RequirePortValue(three_d_outlet.area_m2, "3D outlet area");
				row.downstream_root_pressure_pa = RequirePortValue(downstream_root.mean_pressure_pa, "downstream root pressure");
				row.downstream_root_outward_flow_m3_s = RequirePortValue(downstream_root.outward_flow_m3_s, "downstream root flow");
				row.downstream_root_area_m2 = RequirePortValue(downstream_root.area_m2, "downstream root area");
				row.downstream_terminal_pressure_pa = RequirePortValue(downstream_terminal_port.mean_pressure_pa, "downstream terminal pressure");
				row.downstream_terminal_outward_flow_m3_s = RequirePortValue(downstream_terminal_port.outward_flow_m3_s, "downstream terminal flow");
				row.downstream_terminal_area_m2 = RequirePortValue(downstream_terminal_port.area_m2, "downstream terminal area");
				row.upstream_three_d_flow_residual_m3_s = row.upstream_terminal_outward_flow_m3_s+row.three_d_inlet_outward_flow_m3_s;
				row.three_d_downstream_flow_residual_m3_s = row.three_d_outlet_outward_flow_m3_s+row.downstream_root_outward_flow_m3_s;
				row.upstream_three_d_normalized_residual = iga::ExplicitCouplingNormalizedResidual(row.upstream_terminal_outward_flow_m3_s, row.three_d_inlet_outward_flow_m3_s);
				row.three_d_downstream_normalized_residual = iga::ExplicitCouplingNormalizedResidual(row.three_d_outlet_outward_flow_m3_s, row.downstream_root_outward_flow_m3_s);
				row.three_d_wall_outward_flow_m3_s = RequirePortValue(
					three_d_wall.outward_flow_m3_s, "3D wall flow");
				row.three_d_mass_imbalance_m3_s = row.three_d_inlet_outward_flow_m3_s
					+row.three_d_outlet_outward_flow_m3_s+row.three_d_wall_outward_flow_m3_s;
				row.net_external_outward_flow_m3_s = row.upstream_root_outward_flow_m3_s
					+row.downstream_terminal_outward_flow_m3_s;
				row.external_pressure_drop_pa = row.upstream_root_pressure_pa-row.downstream_terminal_pressure_pa;
				row.upstream_three_d_pressure_jump_pa = row.upstream_terminal_pressure_pa-row.three_d_inlet_pressure_pa;
				row.three_d_downstream_pressure_jump_pa = row.three_d_outlet_pressure_pa-row.downstream_root_pressure_pa;
				row.three_d_trial_linear_iterations = three_d.TrialLinearIterations();
				iga::ValidateExplicitCouplingHistoryRow(row);
				if (injected_failure_step == step)
					throw std::runtime_error("injected explicit coupling failure before commit");
				upstream.CommitStep();
				three_d.CommitStep();
				downstream.CommitStep();
				history.push_back(row);
				lagged_three_d_inlet_pressure = row.three_d_inlet_pressure_pa;
				lagged_downstream_root_pressure = row.downstream_root_pressure_pa;
			} catch (...) {
				RollbackSolved(upstream);
				RollbackSolved(three_d);
				RollbackSolved(downstream);
				throw;
			}
		}
		int output_failed = 0;
		std::string output_error;
		if (rank == 0) try {
			fs::create_directories(options.output_directory);
			std::ofstream output(options.output_directory/"explicit_coupling_history.csv");
			if (!output) throw std::runtime_error("cannot create explicit coupling history");
			iga::WriteExplicitCouplingHistoryHeader(output);
			for (const auto& row : history) iga::WriteExplicitCouplingHistoryRow(output, row);
			if (!output) throw std::runtime_error("cannot write explicit coupling history");
			WriteExplicitCouplingManifest(options.output_directory/"explicit_coupling_manifest.json",
				options.upstream_terminal_node, ports.inlet_label, ports.outlet_labels.front(), scalar.dt_s,
				scalar.steps, static_cast<int>(history.size()), scalar.density_kg_m3,
				scalar.dynamic_viscosity_pa_s, normalized_length,
				reference_inlet_flow, initial_lagged_three_d_inlet_pressure,
				initial_lagged_downstream_root_pressure);
		} catch (const std::exception& error) {
			output_failed = 1;
			output_error = error.what();
		}
		MPI_Bcast(&output_failed, 1, MPI_INT, 0, PETSC_COMM_WORLD);
		if (output_failed) {
			if (rank == 0) throw std::runtime_error(output_error);
			throw std::runtime_error("explicit coupling output failed on rank 0");
		}
		if (rank == 0) std::cout << "completed explicit 1D--3D coupling steps=" << history.size()
			<< " output=" << options.output_directory << '\n';
	} catch (const std::exception& error) {
		if (rank == 0) std::cerr << error.what() << '\n';
		status = 1;
	}
	PetscFinalize();
	return status;
}
