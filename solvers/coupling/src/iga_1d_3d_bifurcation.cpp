#include "CheckedText.hpp"
#include "GraphAcceptedHistory.hpp"
#include "NativeGraphCheckpoint.hpp"
#include "NativeGraphCheckpointSignal.hpp"
#include "ExecutionResources.hpp"
#include "MultidomainRunner.hpp"
#include "CollectiveFailure.hpp"
#include "CollectiveAssetInput.hpp"
#include "CollectivePetscOptions.hpp"
#include "BifurcationMultidomainCase.hpp"
#include "BoundarySupport.hpp"
#include "DomainRuntimeRegistry.hpp"
#include "CollectiveDomainRuntimeRegistry.hpp"
#include "ExplicitOneDThreeDCoupling.hpp"
#include "IgaDatabase.hpp"
#include "ImmersedFlowCase.hpp"
#include "OneDFlowDomainAdapter.hpp"
#include "OneDImplicit.hpp"
#include "OneDRuntime.hpp"
#include "PressureFlowComponentExecutor.hpp"
#include "CollectivePressureFlowExecution.hpp"
#include "CollectiveSpeciesPressureFlowExecution.hpp"
#include "SpeciesPressureFlowComponentExecutor.hpp"
#include "ThreeDBodyFittedFlowDomainAdapter.hpp"
#include "ThreeDBodyFittedFlowTransportDomainAdapter.hpp"
#include "ThreeDImmersedFlowDomain.hpp"
#include "ZeroDFlowDomain.hpp"
#include "ThreeDFlowCoupling.hpp"
#include "ThreeDVcaCoupling.hpp"
#include "TransientFlowRuntime.hpp"
#include "TransientTransportRuntime.hpp"

#include <petscsys.h>

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
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr int kMaximumNewtonIterations = 30;
constexpr double kNonlinearRelativeTolerance = 1.0e-5;
constexpr double kNonlinearAbsoluteTolerance = 1.0e-10;
constexpr double kMassRelativeTolerance = 1.0e-3;

struct Options {
	fs::path graph_case;
	fs::path output_directory;
	int stop_after_step = 0;
	int maximum_newton = kMaximumNewtonIterations;
	fs::path checkpoint_directory, restart_directory;
	int checkpoint_every = 1;
	bool checkpoint_interval_explicit = false;
};

int PositiveInteger(const std::string& text, const std::string& option)
{
	std::size_t used = 0;
	int value = 0;
	try { value = std::stoi(text, &used); }
	catch (const std::exception&) { throw std::runtime_error(option+" requires a positive integer"); }
	if (used != text.size() || value < 1)
		throw std::runtime_error(option+" requires a positive integer");
	return value;
}

Options ParseOptions(int argc, char** argv)
{
	Options options;
	for (int i = 1; i < argc; ++i) {
		const std::string argument(argv[i]);
		if (argument == "--graph-case" || argument == "--output-dir"
			|| argument == "--stop-after-step" || argument == "--three-d-max-newton"
			|| argument == "--checkpoint-dir" || argument == "--restart-dir" || argument == "--checkpoint-every") {
			if (++i >= argc) throw std::runtime_error(argument+" requires a value");
			if (argument == "--graph-case") options.graph_case = argv[i];
			else if (argument == "--output-dir") options.output_directory = argv[i];
			else if (argument == "--checkpoint-dir") {
				options.checkpoint_directory = argv[i];
				if (options.checkpoint_directory.empty()) throw std::runtime_error("--checkpoint-dir requires a nonempty path");
			}
			else if (argument == "--restart-dir") {
				options.restart_directory = argv[i];
				if (options.restart_directory.empty()) throw std::runtime_error("--restart-dir requires a nonempty path");
			}
			else if (argument == "--checkpoint-every") {
				options.checkpoint_every = PositiveInteger(argv[i], argument); options.checkpoint_interval_explicit = true;
			}
			else if (argument == "--stop-after-step")
				options.stop_after_step = PositiveInteger(argv[i], argument);
			else options.maximum_newton = PositiveInteger(argv[i], argument);
			continue;
		}
		if (iga::IsExplicitCouplingPetscOption(argument)) {
			if (i+1 < argc
				&& iga::ExplicitCouplingPetscOptionConsumesNextValue(argument, argv[i+1])) ++i;
			continue;
		}
		throw std::runtime_error("unexpected argument: "+argument);
	}
	if (options.graph_case.empty() || options.output_directory.empty())
		throw std::runtime_error("usage: iga_1d_3d_bifurcation --graph-case ROOT "
			"--output-dir DIR [--stop-after-step N] [--three-d-max-newton N] "
			"[--checkpoint-dir ROOT] [--checkpoint-every N] [--restart-dir ROOT] [PETSc options]");
	if (options.checkpoint_interval_explicit && options.checkpoint_directory.empty())
		throw std::runtime_error("--checkpoint-every requires --checkpoint-dir");
	return options;
}

std::string ReadText(const fs::path& path)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot open file: "+path.string());
	const auto text = iga::ReadCheckedText(input);
	return text;
}

void CanonicalizePeriodicTableAssets(std::vector<iga::TemporalFunctionDefinition>& functions,
	const fs::path& case_directory, const std::string& context)
{
	for (auto& function : functions) {
		if (function.kind != iga::TemporalFunctionKind::PeriodicTable) continue;
		const auto resolved = iga::ResolveContainedCaseFile(case_directory, function.file,
			context+" temporal function '"+function.name+"'");
		function.file = fs::relative(resolved, case_directory).generic_string();
	}
}

const iga::OneDFlowSystemDefinition& SelectFlow(const iga::OneDConfiguration& configuration)
{
	if (configuration.flow_systems.size() != 1)
		throw std::runtime_error("bifurcation coupling requires exactly one 1D flow system per case");
	return configuration.flow_systems.front();
}

void RequireOneDFlowOnly(const iga::OneDConfiguration& configuration)
{
	if (configuration.coupling.mode != iga::SimulationScopeMode::FlowOnly
		|| !configuration.transport_systems.empty() || configuration.physiology.enabled)
		throw std::runtime_error(
			"bifurcation coupling requires flow_only 1D cases without transport or physiology");
}

void ApplyInitialPressure(iga::SimulationConfiguration& configuration,
	const std::string& pressure_field, int label, double pressure)
{
	iga::RequireFinitePortValue("initial 3D outlet pressure", pressure);
	iga::FieldBoundaryCondition* match = nullptr;
	for (auto& boundary : configuration.boundaries) {
		if (boundary.label != label) continue;
		for (auto& condition : boundary.conditions) {
			if (condition.field != pressure_field) continue;
			if (match || condition.kind != iga::FieldBoundaryKind::PressureTraction)
				throw std::runtime_error(
					"bifurcation 3D outlet requires exactly one pressure_traction declaration");
			match = &condition;
		}
	}
	if (!match) throw std::runtime_error("bifurcation 3D outlet has no pressure_traction declaration");
	match->value = {pressure};
	match->waveform.clear();
}

double RequireValue(const std::optional<double>& value, const std::string& context)
{
	if (!value) throw std::runtime_error("bifurcation port omits "+context);
	iga::RequireFinitePortValue(context.c_str(), *value);
	return *value;
}

std::string JsonEscape(const std::string& value)
{
	std::ostringstream output;
	output.exceptions(std::ios::badbit | std::ios::failbit);
	for (const unsigned char character : value) {
		switch (character) {
		case '"': output << "\\\""; break;
		case '\\': output << "\\\\"; break;
		case '\b': output << "\\b"; break;
		case '\f': output << "\\f"; break;
		case '\n': output << "\\n"; break;
		case '\r': output << "\\r"; break;
		case '\t': output << "\\t"; break;
		default:
			if (character < 0x20)
				output << "\\u" << std::hex << std::setw(4) << std::setfill('0')
					<< static_cast<unsigned int>(character) << std::dec << std::setfill(' ');
			else output << character;
		}
	}
	return output.str();
}

std::string CsvEscape(const std::string& value)
{
	if (value.find_first_of(",\"\r\n") == std::string::npos) return value;
	std::string result = "\"";
	for (const char character : value) {
		if (character == '"') result.push_back('"');
		result.push_back(character);
	}
	return result+'"';
}

const char* ExecutionName(iga::GraphExecutionKind kind)
{
	if (kind == iga::GraphExecutionKind::Explicit) return "explicit";
	if (kind == iga::GraphExecutionKind::Fixed) return "fixed";
	if (kind == iga::GraphExecutionKind::Aitken) return "aitken";
	throw std::runtime_error("unknown graph execution kind");
}

struct NativeOneD {
	std::string domain_id;
	fs::path case_directory;
	iga::OneDInletPolicy inlet_policy = iga::OneDInletPolicy::ConfiguredOpenLoop;
	double initial_native_inlet_flow_m3_s = 0.0;
	std::unique_ptr<iga::OneDFlowRuntime> runtime;
};

struct NativeThreeD {
	std::string domain_id;
	fs::path case_directory;
	std::unique_ptr<iga::Database> database;
	iga::SimulationConfiguration configuration;
	iga::LabeledHexMesh mesh;
	std::vector<std::array<double, 3>> boundary_velocity;
	std::set<std::int32_t> wall_trace_basis;
	iga::SimulationConfiguration initial_configuration;
	iga::ResolvedBoundaryConditions initial_boundaries;
	std::vector<iga::OutletModelState> outlet_models;
	std::map<std::string, double> reference_outward_flow_m3_s;
	std::unique_ptr<iga::TransientFlowRuntime> runtime;
	std::unique_ptr<iga::TransientTransportRuntime> transport_runtime;
	std::optional<iga::CompiledLinearSystem> transport_system;
};

// The immersed production owner is itself a CoupledDomainRuntime and retains
// every catalog, runtime, and internal adapter needed by the registry.
using NativeImmersed = iga::ImmersedFlowCase;

std::size_t CountThreeDDomains(const iga::SimulationGraph& graph)
{
	std::size_t count = 0;
	for (const auto& domain : graph.Domains())
		if (iga::DomainDimensionOf(domain.second.kind) == 3) ++count;
	return count;
}

const char* ManifestDomainKind(iga::DomainKind kind)
{
	if (kind == iga::DomainKind::OneDFlow) return "network_flow";
	if (kind == iga::DomainKind::ThreeDBodyFittedFlow) return "body_fitted_iga_flow";
	if (kind == iga::DomainKind::ThreeDImmersedFlow) return "three_d_immersed_flow";
	if (kind == iga::DomainKind::ZeroDFlow) return "zero_d_flow";
	throw std::runtime_error("unknown domain kind in output manifest");
}

void RequireOneDStagedTransport(const iga::OneDConfiguration& configuration,
	const std::map<std::string, std::string>& bindings)
{
	if (configuration.flow_systems.size() != 1 || configuration.transport_systems.size() != 1
		|| configuration.coupling.mode != iga::SimulationScopeMode::FlowOnly
		|| configuration.physiology.enabled || configuration.physiology.vasodilation)
		throw std::runtime_error(
			"schema-v6 multidomain transport requires one 1D flow system, one transport system, flow_only coupling, and no physiology feedback");
	const auto& transport = configuration.transport_systems.front();
	if (transport.flow_system != configuration.flow_systems.front().name)
		throw std::runtime_error("schema-v6 1D transport must reference the configured flow system");
	std::set<std::string> native_fields;
	for (const auto& species : transport.species) native_fields.insert(species.field);
	for (const auto& binding : bindings)
		if (!native_fields.count(binding.second))
			throw std::runtime_error("schema-v6 1D species binding names no transported native field");
}

NativeOneD BuildOneD(const iga::MultidomainConfiguration& graph,
	const std::map<std::string, iga::ResolvedGraphDomainAssets>& assets,
	const std::string& domain_id, bool species_mode, MPI_Comm communicator,
	const std::string& configuration_text, const std::string& checkpoint_identity = {},
	iga::OneDPetscSolverContext* solver_context = nullptr)
{
	const auto& definition = iga::GraphDomainDefinitionFor(graph, domain_id);
	const auto case_directory = assets.at(domain_id).case_directory;
	auto configuration = iga::ParseOneDConfiguration(configuration_text);
	CanonicalizePeriodicTableAssets(configuration.temporal_functions, case_directory, domain_id);
	if (species_mode)
		RequireOneDStagedTransport(configuration, definition.species_bindings);
	else RequireOneDFlowOnly(configuration);
	const auto flow = SelectFlow(configuration);
	const auto geometry = iga::ResolveContainedCaseFile(case_directory,
		configuration.geometry.file, domain_id+" geometry");
	auto network = iga::ReadOneDNetwork(geometry,
		configuration.geometry.length_scale_to_m, flow.discretization.cells_per_segment,
		flow.dynamic_viscosity, configuration.geometry.root_node_id);
	iga::ValidateOneDTopologyReferences(configuration, network);
	auto inlet = iga::ResolveOneDInlet(configuration);
	auto runtime = std::make_unique<iga::OneDFlowRuntime>(configuration, flow,
		std::move(network), inlet, case_directory,
		[communicator, solver_context](const iga::OneDNetwork& network_definition,
			const iga::OneDFlowSystemDefinition& flow_definition,
			iga::OneDFlowState& state, double inlet_flow, double dt) {
			iga::AdvanceImplicitOneD(network_definition, flow_definition, state, inlet_flow, dt, communicator, solver_context);
		}, [communicator](const char* stage, std::exception_ptr error) {
			iga::CollectiveLocalStage(communicator, stage, [&] {
				if (error) std::rethrow_exception(error);
			});
		}, checkpoint_identity);
	(void)iga::MakeOneDSubcyclingPlan(configuration.time.dt, configuration.time.steps,
		graph.time.dt_s, graph.time.steps);
	const double seed = iga::EvaluateOneDInlet(runtime->Configuration(),
		runtime->InletDefinition(), case_directory, 0.0,
		runtime->Network().segments.front().area0);
	if (!(seed > 0.0) || !std::isfinite(seed))
		throw std::runtime_error("bifurcation 1D initial native inlet flow must be positive");
	runtime->InitializeOpenLoop(seed);
	return {domain_id, case_directory, definition.one_d_inlet_policy, seed,
		std::move(runtime)};
}

void ValidateOneDPort(const NativeOneD& native, const iga::CouplingPort& port)
{
	if (port.locator == "root") return;
	const int node_id = iga::ParseOneDOutletNodeLocator(port);
	const auto found = native.runtime->Network().node_index.find(node_id);
	if (found == native.runtime->Network().node_index.end())
		throw std::runtime_error("schema-v5 1D terminal locator does not match its native case");
	const auto outlet = std::find_if(native.runtime->FlowState().outlets.begin(),
		native.runtime->FlowState().outlets.end(), [&](const iga::OneDOutletState& candidate) {
			return candidate.node == found->second;
		});
	if (outlet == native.runtime->FlowState().outlets.end())
		throw std::runtime_error("schema-v5 1D terminal locator does not match its native case");
	if (port.requires.count(iga::PortQuantity::MeanPressure)
		&& outlet->kind != iga::OneDOutletKind::Pressure)
		throw std::runtime_error("graph-controlled 1D outlet requires a native pressure closure");
}

std::unique_ptr<NativeThreeD> BuildThreeDPreflight(
	const iga::MultidomainConfiguration& graph,
	const iga::PressureFlowComponentPlan& plan,
	const std::map<std::string, iga::ResolvedGraphDomainAssets>& assets,
	const std::string& domain_id, int rank, bool species_mode,
	const std::string& configuration_text)
{
	auto native = std::make_unique<NativeThreeD>();
	native->domain_id = domain_id;
	native->case_directory = assets.at(domain_id).case_directory;
	native->database = std::make_unique<iga::Database>(assets.at(domain_id).database.string());
	const auto required_elements = native->database->LoadRequired(rank);
	const auto owned_elements = native->database->LoadOwned(rank);
	native->configuration = iga::ParseSimulationConfiguration(configuration_text);
	CanonicalizePeriodicTableAssets(native->configuration.temporal_functions,
		native->case_directory, domain_id);
	int flow_count = 0;
	int transport_count = 0;
	for (const auto& system : native->configuration.equation_systems) {
		if (system.kind == iga::EquationKind::NavierStokes) ++flow_count;
		if (system.kind == iga::EquationKind::LinearTransport) ++transport_count;
	}
	const auto& flow_system = iga::FirstNavierStokesSystem(native->configuration);
	if (flow_count != 1 || flow_system.unknowns.size() != 2
		|| native->configuration.physiology.enabled
		|| flow_system.time_integration != "backward_euler"
		|| (!species_mode && (native->configuration.equation_systems.size() != 1
			|| native->configuration.coupling.mode != iga::SimulationScopeMode::FlowOnly))
		|| (species_mode && native->configuration.equation_systems.size() != 2))
		throw std::runtime_error(
			"multidomain flow requires exactly one backward_euler 3D Navier-Stokes system");
	if (species_mode) {
		if (transport_count != 1)
			throw std::runtime_error("schema-v6 multidomain transport requires exactly one 3D linear_transport system");
		if (native->configuration.coupling.mode != iga::SimulationScopeMode::FlowOnly)
			throw std::runtime_error(
				"schema-v6 multidomain transport rejects 3D coupling feedback");
		native->transport_system = iga::RequireThreeDVcaTransportSystem(native->configuration);
		if (native->transport_system->velocity_source != "prescribed")
			throw std::runtime_error(
				"schema-v6 multidomain transport requires prescribed 3D transport velocity");
		for (const auto& binding : graph.graph.Domain(domain_id).species_bindings)
			if (!native->transport_system->field_index.count(binding.second))
				throw std::runtime_error("schema-v6 3D species binding names no compiled scalar field");
	} else for (const auto& field : native->configuration.fields)
		if (field.kind == iga::FieldKind::Scalar)
			throw std::runtime_error("multidomain flow rejects 3D transport fields");
	if (graph.time.dt_s != native->configuration.time.dt
		|| graph.time.steps != native->configuration.time.steps)
		throw std::runtime_error("schema-v5 and every 3D time grid must match");
	native->mesh = iga::ReadLabeledHexMesh(
		iga::ResolveContainedCaseFile(native->case_directory, "controlmesh.vtk",
			domain_id+" 3D control mesh").string(), native->database->header().nodes,
		native->database->header().elements);
	native->boundary_velocity = iga::ReadVelocity(
		iga::ResolveContainedCaseFile(native->case_directory,
			"initial_velocityfield.txt", domain_id+" 3D initial velocity").string(),
		native->database->header().nodes);
	native->wall_trace_basis = iga::WallTraceBasis(*native->database, native->mesh, 0);
	native->initial_configuration = iga::MaterializeBoundaryWaveforms(
		native->configuration, native->case_directory.string(), 0.0);
	native->outlet_models = iga::InitializeOutletModels(native->configuration, flow_system);
	if (!native->outlet_models.empty())
		throw std::runtime_error(
			"multidomain graph-controlled 3D domains require static-pressure boundaries");
	std::set<int> declared_labels;
	for (const auto& port : graph.graph.Domain(domain_id).ports) {
		const int label = iga::ParseThreeDFlowBoundaryLabel(port);
		declared_labels.insert(label);
		if (!port.provides.count(iga::PortQuantity::FlowRate))
			throw std::runtime_error(
				"every multidomain 3D boundary port must provide flow for balance auditing");
		if (port.requires.count(iga::PortQuantity::MeanPressure))
		{
			const auto interface = std::find_if(plan.interfaces.begin(), plan.interfaces.end(),
				[&](const iga::PressureFlowInterfacePlan& candidate) {
					return candidate.pressure_receiver == iga::PortRef{domain_id, port.id};
				});
			if (interface == plan.interfaces.end())
				throw std::runtime_error("3D pressure receiver is not bound to a graph edge");
			ApplyInitialPressure(native->initial_configuration,
				flow_system.unknowns.at(1), label,
				graph.initial_pressure_pa.at(interface->edge_id));
		}
		if (port.requires.count(iga::PortQuantity::FlowRate)) {
			auto check = native->initial_configuration;
			iga::PortBoundaryData input;
			input.time_s = 0.0;
			input.outward_flow_m3_s = 1.0;
			iga::ApplyThreeDReferenceProfileInput(check,
				iga::FirstNavierStokesSystem(check), port, input, 1.0);
		}
	}
	for (const auto& element : required_elements)
		for (const int label : element.boundary_labels)
			if (label >= 0 && !declared_labels.count(label))
				throw std::runtime_error("3D database boundary label lacks a logical audit port");
	for (const auto& element : owned_elements)
		for (const int label : element.boundary_labels)
			if (label >= 0 && !declared_labels.count(label))
				throw std::runtime_error("owned 3D boundary label lacks a logical audit port");
	native->initial_boundaries = iga::ResolveFlowBoundaries(
		native->initial_configuration,
		iga::FirstNavierStokesSystem(native->initial_configuration), native->mesh.labels,
		native->boundary_velocity);
	return native;
}

using AcceptedStep = iga::GraphAcceptedFlowStep;
using AcceptedSpeciesStep = iga::GraphAcceptedSpeciesStep;

void WriteOutputs(const fs::path& directory,
	const iga::MultidomainConfiguration& configuration,
	const fs::path& graph_root,
	const std::map<std::string, iga::ResolvedGraphDomainAssets>& assets,
	const std::optional<iga::OneDThreeDBifurcationDefinition>& bifurcation,
	const std::map<std::string, NativeOneD>& one_d,
	const std::map<std::string, const NativeImmersed*>& immersed,
	const std::vector<AcceptedStep>& steps)
{
	if (!fs::create_directories(directory))
		throw std::runtime_error("bifurcation output directory must be newly created");
	std::ofstream summary(directory/"pressure_flow_steps.csv");
	std::ofstream edges(directory/"pressure_flow_edges.csv");
	std::ofstream iterations(directory/"pressure_flow_iterations.csv");
	std::ofstream ports(directory/"pressure_flow_ports.csv");
	std::ofstream initialization(directory/"one_d_initialization.csv");
	std::ofstream balances(directory/"pressure_flow_domain_balances.csv");
	std::ofstream zero_d_history(directory/"zero_d_flow_history.csv");
	if (!summary || !edges || !iterations || !ports || !initialization || !balances
		|| !zero_d_history)
		throw std::runtime_error("cannot create bifurcation long-form output");
	summary << "step,time_s,iterations,three_d_mass_imbalance_m3_s,external_outward_flow_m3_s\n";
	edges << "step,time_s,edge_id,applied_pressure_pa,measured_pressure_pa,pressure_residual_pa,"
		"normalized_pressure_residual,first_outward_flow_m3_s,second_outward_flow_m3_s,"
		"flow_residual_m3_s,normalized_flow_residual\n";
	iterations << "step,time_s,iteration,applied_relaxation,converged,edge_id,"
		"applied_pressure_pa,measured_pressure_pa,pressure_residual_pa,"
		"normalized_pressure_residual,first_outward_flow_m3_s,second_outward_flow_m3_s,"
		"flow_residual_m3_s,normalized_flow_residual\n";
	ports << "step,time_s,domain_id,port_id,area_m2,outward_flow_m3_s,mean_pressure_pa\n";
	initialization << "domain_id,inlet_policy,native_inlet_flow_m3_s\n";
	balances << "step,time_s,domain_id,boundary_flow_sum_m3_s,"
		"boundary_flow_absolute_sum_m3_s,normalized_mass_imbalance\n";
	zero_d_history << "step,time_s,domain_id,stored_pressure_pa,initial_stored_volume_m3,"
		"final_stored_volume_m3,prescribed_source_amount_m3,distal_sink_amount_m3,"
		"outward_graph_port_amount_m3,residual_m3\n";
	for (const auto& item : one_d)
		initialization << std::setprecision(17) << CsvEscape(item.first) << ','
			<< (item.second.inlet_policy == iga::OneDInletPolicy::ConfiguredOpenLoop
				? "configured_open_loop" : "coupled_root") << ','
			<< item.second.initial_native_inlet_flow_m3_s << '\n';
	for (const auto& step : steps) {
		summary << std::setprecision(17) << step.step << ',' << step.time_s << ','
			<< step.iterations << ',' << step.three_d_mass_m3_s << ','
			<< step.external_outward_flow_m3_s << '\n';
		for (const auto& edge : step.result.iterations.back().edges)
			edges << std::setprecision(17) << step.step << ',' << step.time_s << ','
				<< CsvEscape(edge.edge_id) << ',' << edge.applied_pressure_pa << ','
				<< edge.measured_pressure_pa << ',' << edge.pressure_residual_pa << ','
				<< edge.normalized_pressure_residual << ',' << edge.first_outward_flow_m3_s
				<< ',' << edge.second_outward_flow_m3_s << ',' << edge.flow_residual_m3_s
				<< ',' << edge.normalized_flow_residual << '\n';
		for (const auto& iteration : step.result.iterations)
			for (const auto& edge : iteration.edges)
				iterations << std::setprecision(17) << step.step << ',' << step.time_s << ','
					<< iteration.iteration << ',' << iteration.applied_relaxation << ','
					<< (iteration.converged ? 1 : 0) << ',' << CsvEscape(edge.edge_id) << ','
					<< edge.applied_pressure_pa << ',' << edge.measured_pressure_pa << ','
					<< edge.pressure_residual_pa << ',' << edge.normalized_pressure_residual
					<< ',' << edge.first_outward_flow_m3_s << ','
					<< edge.second_outward_flow_m3_s << ',' << edge.flow_residual_m3_s << ','
					<< edge.normalized_flow_residual << '\n';
		for (const auto& port : step.result.accepted_ports) {
			ports << std::setprecision(17) << step.step << ',' << step.time_s << ','
				<< CsvEscape(port.first.domain_id) << ',' << CsvEscape(port.first.port_id) << ',';
			if (port.second.area_m2) ports << *port.second.area_m2;
			ports << ',' << RequireValue(port.second.outward_flow_m3_s, "port flow") << ',';
			if (port.second.mean_pressure_pa) ports << *port.second.mean_pressure_pa;
			ports << '\n';
		}
		for (const auto& zero : step.zero_d_states) {
			const auto& accounting = step.zero_d_accounting.at(zero.first);
			zero_d_history << std::setprecision(17) << step.step << ',' << step.time_s << ','
				<< CsvEscape(zero.first) << ',' << zero.second.stored_pressure_pa << ','
				<< accounting.initial_stored_volume_m3 << ','
				<< accounting.final_stored_volume_m3 << ','
				<< accounting.prescribed_source_amount_m3 << ','
				<< accounting.distal_sink_amount_m3 << ','
				<< accounting.outward_graph_port_amount_m3 << ','
				<< accounting.residual_m3 << '\n';
		}
		for (const auto& balance : step.three_d_balance)
			balances << std::setprecision(17) << step.step << ',' << step.time_s << ','
				<< CsvEscape(balance.first) << ',' << balance.second.first << ','
				<< balance.second.second << ','
				<< (balance.second.second == 0.0 ? 0.0
					: 2.0*std::abs(balance.second.first)/balance.second.second) << '\n';
	}
	summary.close();
	edges.close();
	iterations.close();
	ports.close();
	initialization.close();
	balances.close();
	zero_d_history.close();
	if (!summary || !edges || !iterations || !ports || !initialization || !balances
		|| !zero_d_history)
		throw std::runtime_error("cannot finalize bifurcation long-form output");
	std::ofstream marker(directory/"graph_binding_manifest.json.tmp");
	if (!marker) throw std::runtime_error("cannot create bifurcation completion marker");
	marker << std::setprecision(17);
	if (bifurcation) {
		marker << "{\n  \"schema_version\": 5,\n  \"benchmark\": \"one_d_three_d_bifurcation\",\n"
			<< "  \"start_domain\": \"" << JsonEscape(configuration.start_domain_id) << "\",\n"
			<< "  \"three_d_domain\": \"" << JsonEscape(bifurcation->three_d_domain_id)
			<< "\",\n  \"branch_count\": " << bifurcation->branches.size() << ",\n"
			<< "  \"completed_steps\": " << steps.size() << "\n}\n";
	} else {
		marker << "{\n  \"schema_version\": 5,\n  \"benchmark\": \"acyclic_multidomain_flow\",\n"
			<< "  \"graph_root\": \"" << JsonEscape(graph_root.generic_string()) << "\",\n"
			<< "  \"start_domain\": \"" << JsonEscape(configuration.start_domain_id) << "\",\n"
			<< "  \"execution\": \"" << ExecutionName(configuration.execution.kind) << "\",\n"
			<< "  \"domain_count\": " << configuration.graph.Domains().size() << ",\n"
			<< "  \"three_d_domain_count\": " << CountThreeDDomains(configuration.graph) << ",\n"
			<< "  \"completed_steps\": " << steps.size() << ",\n  \"domains\": [\n";
		std::size_t domain_index = 0;
		for (const auto& domain : configuration.graph.Domains()) {
			const auto& resolved = assets.at(domain.first);
			marker << "    {\"id\":\"" << JsonEscape(domain.first) << "\",\"kind\":\""
				<< ManifestDomainKind(domain.second.kind) << "\",\"case\":\""
				<< JsonEscape(resolved.case_directory.generic_string()) << "\"";
			if (domain.second.kind == iga::DomainKind::ThreeDBodyFittedFlow)
				marker << ",\"database\":\""
					<< JsonEscape(resolved.database.generic_string()) << "\"";
			if (domain.second.kind == iga::DomainKind::ZeroDFlow) {
				const auto& definition = iga::GraphDomainDefinitionFor(configuration, domain.first);
				const auto& model = definition.zero_d_flow_model->model;
				marker << ",\"model_role\":\"" << iga::ZeroDFlowRoleName(model.role)
					<< "\",\"model_identity_sha256\":\""
					<< iga::BuildZeroDFlowModelIdentitySha256(model) << "\"";
				if (!steps.empty()) {
					const auto state = steps.back().zero_d_states.find(domain.first);
					const auto accounting = steps.back().zero_d_accounting.find(domain.first);
					if (state == steps.back().zero_d_states.end()
						|| accounting == steps.back().zero_d_accounting.end())
						throw std::runtime_error("accepted 0D output is missing final state or accounting");
					marker << ",\"final_state_identity_sha256\":\""
						<< iga::BuildZeroDFlowStateIdentitySha256(model, state->second)
						<< "\",\"final_step_accounting_identity_sha256\":\""
						<< iga::BuildZeroDFlowStepAccountingIdentitySha256(model, accounting->second)
						<< "\"";
				}
			}
			if (domain.second.kind == iga::DomainKind::ThreeDImmersedFlow) {
				const auto& native = *immersed.at(domain.first);
				const auto& grid = native.Grid();
				const auto& classification = native.ClassificationDiagnostics();
				const auto& volume = native.VolumeDiagnostics();
				const auto& surface = native.SurfaceDiagnostics();
				const auto& ghost = native.GhostDiagnostics();
				marker << ",\"surface_hash\":\"" << JsonEscape(native.SurfaceHash())
					<< "\",\"grid\":{\"lower_m\":[" << grid.lower_m[0] << ','
					<< grid.lower_m[1] << ',' << grid.lower_m[2] << "],\"upper_m\":["
					<< grid.upper_m[0] << ',' << grid.upper_m[1] << ',' << grid.upper_m[2]
					<< "],\"cells\":[" << grid.cells[0] << ',' << grid.cells[1] << ','
					<< grid.cells[2] << "]},\"catalog_audit\":{\"active_cells\":"
					<< classification.inside_count+classification.cut_count
					<< ",\"volume_points\":" << volume.output_points
					<< ",\"surface_points\":" << surface.output_points
					<< ",\"ghost_faces\":" << ghost.selected_faces << '}';
			}
			marker << '}'
				<< (++domain_index == configuration.graph.Domains().size() ? "\n" : ",\n");
		}
		marker << "  ],\n  \"couplings\": [\n";
		std::vector<const iga::CouplingEdge*> sorted_edges;
		for (const auto& edge : configuration.graph.Edges()) sorted_edges.push_back(&edge);
		std::sort(sorted_edges.begin(), sorted_edges.end(),
			[](const iga::CouplingEdge* first, const iga::CouplingEdge* second) {
				return first->id < second->id;
			});
		for (std::size_t index = 0; index < sorted_edges.size(); ++index) {
			const auto& edge = *sorted_edges[index];
			const iga::PortRef* first = &edge.first;
			const iga::PortRef* second = &edge.second;
			if (*second < *first) std::swap(first, second);
			marker << "    {\"id\":\"" << JsonEscape(edge.id)
				<< "\",\"first\":{\"domain\":\""
				<< JsonEscape(first->domain_id) << "\",\"port\":\""
				<< JsonEscape(first->port_id) << "\"},\"second\":{\"domain\":\""
				<< JsonEscape(second->domain_id) << "\",\"port\":\""
				<< JsonEscape(second->port_id) << "\"},\"initial_pressure_pa\":"
				<< std::setprecision(17) << configuration.initial_pressure_pa.at(edge.id) << '}'
				<< (index+1 == sorted_edges.size() ? "\n" : ",\n");
		}
		marker << "  ]\n}\n";
	}
	marker.close();
	if (!marker) throw std::runtime_error("cannot write bifurcation completion marker");
	fs::rename(directory/"graph_binding_manifest.json.tmp",
		directory/"graph_binding_manifest.json");
}

void WriteSpeciesOutputs(const fs::path& directory,
	const iga::MultidomainConfiguration& configuration,
	const fs::path& graph_root,
	const std::map<std::string, iga::ResolvedGraphDomainAssets>& assets,
	const std::map<std::string, NativeOneD>& one_d,
	const std::map<std::string, std::unique_ptr<NativeThreeD>>& three_d,
	const std::vector<AcceptedSpeciesStep>& steps)
{
	if (!fs::create_directories(directory))
		throw std::runtime_error("multidomain output directory must be newly created");
	std::ofstream hydraulic(directory/"species_hydraulic_steps.csv");
	std::ofstream edges(directory/"species_edge_amounts.csv");
	std::ofstream domains(directory/"species_domain_accounting.csv");
	std::ofstream global(directory/"species_global_balance.csv");
	std::ofstream ports(directory/"species_logical_ports.csv");
	if (!hydraulic || !edges || !domains || !global || !ports)
		throw std::runtime_error("cannot create schema-v6 multidomain long-form output");
	hydraulic << "step,time_s,hydraulic_iterations,three_d_mass_imbalance_m3_s,"
		"external_outward_flow_m3_s\n";
	edges << "step,time_s,edge_id,species_id,first_domain_id,first_port_id,"
		"second_domain_id,second_port_id,donor,first_outward_amount,"
		"second_outward_amount,amount_residual,normalized_amount_residual\n";
	domains << "step,time_s,domain_id,species_id,M0_amount,M1_amount,source_amount,"
		"total_outward_amount,recomputed_residual,normalized_residual\n";
	global << "step,time_s,species_id,M0_amount,M1_amount,source_amount,"
		"outward_amount,gross_activity,global_residual,normalized_residual\n";
	ports << "step,time_s,domain_id,port_id,species_id,concentration,"
		"end_step_outward_species_flux\n";
	for (const auto& step : steps) {
		hydraulic << std::setprecision(17) << step.step << ',' << step.time_s << ','
			<< step.hydraulic_iterations << ',' << step.three_d_mass_m3_s << ','
			<< step.external_outward_flow_m3_s << '\n';
		for (const auto& edge : step.result.edge_amounts)
			edges << std::setprecision(17) << step.step << ',' << step.time_s << ','
				<< CsvEscape(edge.edge_id) << ',' << CsvEscape(edge.species_id) << ','
				<< CsvEscape(edge.first.domain_id) << ',' << CsvEscape(edge.first.port_id)
				<< ',' << CsvEscape(edge.second.domain_id) << ','
				<< CsvEscape(edge.second.port_id) << ','
				<< (edge.donor == iga::SpeciesDonor::First ? "first" : "second") << ','
				<< edge.first_outward_amount << ',' << edge.second_outward_amount << ','
				<< edge.residual << ',' << edge.normalized_residual << '\n';
		for (const auto& domain : step.result.domain_balances)
		{
			double total_outward_amount = 0.0;
			for (const auto& amount : domain.accounting.outward_port_amount)
				total_outward_amount += amount.second;
			domains << std::setprecision(17) << step.step << ',' << step.time_s << ','
				<< CsvEscape(domain.domain_id) << ',' << CsvEscape(domain.species_id) << ','
				<< domain.accounting.initial_mass << ',' << domain.accounting.final_mass << ','
				<< domain.accounting.source_amount << ',' << total_outward_amount << ','
				<< domain.recomputed_residual << ','
				<< domain.normalized_residual << '\n';
		}
		for (const auto& item : step.result.global_balances) {
			const auto& balance = item.second;
			global << std::setprecision(17) << step.step << ',' << step.time_s << ','
				<< CsvEscape(item.first) << ',' << balance.initial_mass << ','
				<< balance.final_mass << ',' << balance.source_amount << ','
				<< balance.outward_amount << ',' << balance.gross_activity << ','
				<< balance.residual << ',' << balance.normalized_residual << '\n';
		}
		for (const auto& port : step.transport_ports) {
			std::set<std::string> species;
			for (const auto& concentration : port.second.concentration)
				species.insert(concentration.first);
			for (const auto& flux : port.second.outward_species_flux)
				species.insert(flux.first);
			for (const auto& species_id : species) {
				ports << std::setprecision(17) << step.step << ',' << step.time_s << ','
					<< CsvEscape(port.first.domain_id) << ',' << CsvEscape(port.first.port_id)
					<< ',' << CsvEscape(species_id) << ',';
				const auto concentration = port.second.concentration.find(species_id);
				if (concentration != port.second.concentration.end()) ports << concentration->second;
				ports << ',';
				const auto flux = port.second.outward_species_flux.find(species_id);
				if (flux != port.second.outward_species_flux.end()) ports << flux->second;
				ports << '\n';
			}
		}
	}
	hydraulic.close();
	edges.close();
	domains.close();
	global.close();
	ports.close();
	if (!hydraulic || !edges || !domains || !global || !ports)
		throw std::runtime_error("cannot finalize schema-v6 multidomain long-form output");
	std::ofstream marker(directory/"graph_binding_manifest.json.tmp");
	if (!marker) throw std::runtime_error("cannot create schema-v6 completion marker");
	marker << std::setprecision(17);
	marker << "{\n  \"schema_version\": " << configuration.schema_version
		<< ",\n  \"benchmark\": \"acyclic_multidomain_flow\",\n"
		<< "  \"species_mode\": true,\n  \"graph_root\": \""
		<< JsonEscape(graph_root.generic_string()) << "\",\n  \"completed_steps\": "
		<< steps.size() << ",\n  \"species\": [";
	std::size_t species_index = 0;
	for (const auto& species : configuration.graph.Species())
		marker << (species_index++ == 0 ? "" : ",") << "\""
			<< JsonEscape(species.first) << "\"";
	marker << "],\n  \"species_units\": {\n";
	std::size_t unit_index = 0;
	for (const auto& species : configuration.graph.Species()) {
		const auto& concentration = species.second.concentration_unit;
		marker << "    \"" << JsonEscape(species.first) << "\": {\"concentration\":\""
			<< JsonEscape(concentration) << "\",\"integrated_amount\":\"("
			<< JsonEscape(concentration) << ")*m^3\",\"outward_rate\":\""
			<< JsonEscape(iga::SpeciesFluxUnit(species.second)) << "\"}"
			<< (++unit_index == configuration.graph.Species().size() ? "\n" : ",\n");
	}
	marker << "  },\n  \"domains\": [\n";
	std::size_t domain_index = 0;
	for (const auto& domain : configuration.graph.Domains()) {
		const auto& resolved = assets.at(domain.first);
		marker << "    {\"id\":\"" << JsonEscape(domain.first) << "\",\"kind\":\""
			<< (domain.second.kind == iga::DomainKind::OneDFlow
				? "network_flow" : "body_fitted_iga_flow") << "\",\"case\":\""
			<< JsonEscape(resolved.case_directory.generic_string()) << "\"";
		if (domain.second.kind == iga::DomainKind::ThreeDBodyFittedFlow)
			marker << ",\"database\":\"" << JsonEscape(resolved.database.generic_string())
				<< "\"";
		marker << '}' << (++domain_index == configuration.graph.Domains().size()
			? "\n" : ",\n");
	}
	marker << "  ],\n  \"one_d_domain_count\": " << one_d.size()
		<< ",\n  \"three_d_domain_count\": " << three_d.size() << "\n}\n";
	marker.close();
	if (!marker) throw std::runtime_error("cannot write schema-v6 completion marker");
	fs::rename(directory/"graph_binding_manifest.json.tmp",
		directory/"graph_binding_manifest.json");
}

int FailureInjectionStep()
{
	// Share the explicit-coupling injection knob with the generic graph runner
	// so a Phase-6 closure test can exercise the same precommit-failure
	// contract through its production 1D--immersed--3D--1D entry point.
	const char* value = std::getenv("TUBULARFLOWIGA_INJECT_EXPLICIT_COUPLING_FAILURE_STEP");
	if (value) return PositiveInteger(value,
		"TUBULARFLOWIGA_INJECT_EXPLICIT_COUPLING_FAILURE_STEP");
	value = std::getenv("TUBULARFLOWIGA_INJECT_BIFURCATION_FAILURE_STEP");
	return value ? PositiveInteger(value,
		"TUBULARFLOWIGA_INJECT_BIFURCATION_FAILURE_STEP") : 0;
}


} // namespace

int iga::RunMultidomainFlow(int argc, char** argv, MPI_Comm communicator)
{
	int rank = 0;
	int mpi_size = 1;
	MPI_Comm_rank(communicator, &rank);
	MPI_Comm_size(communicator, &mpi_size);
	int status = 0;
	try {
		iga::RequireExecutionResources(communicator, &std::cout);
		Options options;
		int injected_failure_step = 0;
		std::string execution_controls;
		std::set<std::string> application_options;
		iga::CollectiveLocalStage(communicator, "graph arguments", [&] {
			options = ParseOptions(argc, argv);
			injected_failure_step = FailureInjectionStep();
			execution_controls = std::to_string(options.stop_after_step)+"\n"
				+std::to_string(options.maximum_newton)+"\n"+std::to_string(injected_failure_step)
				+"\n"+options.checkpoint_directory.generic_string()+"\n"+options.restart_directory.generic_string()
				+"\n"+std::to_string(options.checkpoint_every);
			// These application arguments are checked separately and are not
			// read by KSP/PC. options_file is an already-loaded input location;
			// its resulting entries, including prefixed/unused options, agree.
			application_options = {"--graph-case", "--output-dir", "--stop-after-step",
				"--three-d-max-newton", "--checkpoint-dir", "--checkpoint-every", "--restart-dir", "-options_file"};
		});
		iga::RequireCollectiveSameText(communicator, "graph execution controls", execution_controls);
		iga::RequireCollectivePetscOptions(communicator, nullptr, application_options);
		std::optional<iga::NativeGraphCheckpointSignal> checkpoint_signal;
		if (!options.checkpoint_directory.empty()) iga::CollectiveLocalStage(communicator, "graph checkpoint signal handler", [&] { checkpoint_signal.emplace(); });
		std::map<std::string, std::string> input_texts;
		std::optional<iga::MultidomainConfiguration> configuration_holder;
		std::optional<iga::PressureFlowComponentPlan> plan_holder;
		std::optional<iga::OneDThreeDBifurcationDefinition> bifurcation_holder;
		std::map<std::string, iga::ResolvedGraphDomainAssets> assets;
		std::map<std::string, std::unique_ptr<iga::OneDPetscSolverContext>> one_d_solvers;
		std::map<std::string, NativeOneD> one_d;
		std::map<std::string, std::unique_ptr<NativeThreeD>> three_d;
		std::map<std::string, std::unique_ptr<NativeImmersed>> immersed;
		std::map<std::string, const NativeImmersed*> immersed_audit;
		fs::path graph_root;
		iga::CollectiveLocalStage(communicator, "graph input", [&] {
			if (fs::exists(options.output_directory))
				throw std::runtime_error("bifurcation output directory must not already exist");
			graph_root = iga::CanonicalGraphCaseRoot(options.graph_case);
			configuration_holder.emplace(iga::ReadMultidomainConfiguration(
				iga::ResolveContainedCaseFile(graph_root, "simulation_config.json",
					"multidomain manifest").string(), input_texts));
#ifdef IGA_BIFURCATION_ENTRY
			if (configuration_holder->schema_version != 5)
				throw std::runtime_error("bifurcation runner supports schema version 5 only");
#else
			if (configuration_holder->schema_version != 5
				&& configuration_holder->schema_version != 6)
				throw std::runtime_error("multidomain flow runner supports schema version 5 or 6");
#endif
			plan_holder.emplace(iga::MakeAcyclicPressureFlowPlan(
				configuration_holder->graph, configuration_holder->start_domain_id));
#ifdef IGA_BIFURCATION_ENTRY
			bifurcation_holder.emplace(
				iga::ResolveOneDThreeDBifurcation(*configuration_holder));
#endif
			assets = iga::ResolveGraphDomainAssets(*configuration_holder, graph_root);

		});
		iga::RequireCollectiveSameText(communicator, "graph manifest", input_texts.at("graph manifest"));
		iga::CollectiveLocalStage(communicator, "graph configuration input", [&] {
			if (options.stop_after_step > configuration_holder->time.steps)
				throw std::runtime_error("--stop-after-step exceeds configured steps");
			for (const auto& domain : configuration_holder->domains) {
				if (domain.kind != iga::DomainKind::OneDFlow
					&& domain.kind != iga::DomainKind::ThreeDBodyFittedFlow) continue;
				input_texts["domain configuration "+domain.id] = ReadText(
					iga::ResolveContainedCaseFile(assets.at(domain.id).case_directory,
						"simulation_config.json", domain.id+" configuration"));
			}
		});
		// The manifest has already agreed, so all ranks visit the same logical
		// inputs in the same order. Parsed 1D/3D settings use these snapshots.
		for (const auto& input : input_texts)
			iga::RequireCollectiveSameText(communicator, input.first.c_str(), input.second);
		iga::AssetFileCatalog input_assets;
		iga::CollectiveLocalStage(communicator, "graph asset catalog", [&] {
			for (const auto& domain : configuration_holder->domains) {
				if (domain.kind != iga::DomainKind::OneDFlow
					&& domain.kind != iga::DomainKind::ThreeDBodyFittedFlow) continue;
				const auto& directory = assets.at(domain.id).case_directory;
				const auto prefix = "graph asset " + domain.id + " / ";
				const auto add_asset = [&](const std::string& role, const fs::path& path) {
					if (!input_assets.emplace(prefix + role, path).second)
						throw std::runtime_error("duplicate logical graph asset: " + prefix + role);
				};
				const auto add_tables = [&](const auto& functions) {
					for (const auto& function : functions)
						if (function.kind == iga::TemporalFunctionKind::PeriodicTable)
							add_asset("temporal " + function.name,
								iga::ResolveContainedCaseFile(directory, function.file,
									domain.id + " temporal " + function.name));
				};
				if (domain.kind == iga::DomainKind::OneDFlow) {
					const auto input = iga::ParseOneDConfiguration(input_texts.at("domain configuration " + domain.id));
					add_asset("network", iga::ResolveContainedCaseFile(
						directory, input.geometry.file, domain.id + " network"));
					add_tables(input.temporal_functions);
				} else {
					const auto input = iga::ParseSimulationConfiguration(input_texts.at("domain configuration " + domain.id));
					add_asset("database", assets.at(domain.id).database);
					add_asset("control mesh", iga::ResolveContainedCaseFile(
						directory, "controlmesh.vtk", domain.id + " control mesh"));
					add_asset("initial velocity", iga::ResolveContainedCaseFile(
						directory, "initial_velocityfield.txt", domain.id + " initial velocity"));
					add_tables(input.temporal_functions);
				}
			}
		});
		iga::RequireCollectiveAssetFiles(communicator, input_assets);
		const bool checkpoint_enabled = !options.checkpoint_directory.empty() || !options.restart_directory.empty();
		iga::CoupledCheckpointCompatibility checkpoint_identity;
		if (checkpoint_enabled) {
			iga::CollectiveLocalStage(communicator, "graph checkpoint capability", [&] {
				for (const auto& domain : configuration_holder->domains)
					if (domain.kind != iga::DomainKind::OneDFlow && domain.kind != iga::DomainKind::ThreeDBodyFittedFlow && domain.kind != iga::DomainKind::ZeroDFlow)
						throw std::runtime_error("native checkpoint currently supports 0D, 1D and body-fitted 3D graph domains");
			});
			checkpoint_identity = iga::BuildNativeGraphCheckpointIdentity(communicator, input_texts, input_assets, application_options, options.maximum_newton,
				{kNonlinearRelativeTolerance, kNonlinearAbsoluteTolerance, kMassRelativeTolerance});
		}
		for (const auto& domain_id : plan_holder->domain_order) {
			if (configuration_holder->graph.Domain(domain_id).kind != iga::DomainKind::OneDFlow) continue;
			std::string prefix;
			iga::CollectiveLocalStage(communicator, "graph 1D solver preparation", [&] {
				prefix = iga::PetscDomainOptionsPrefix(domain_id, "flow");
				one_d_solvers.emplace(domain_id, nullptr);
			});
			one_d_solvers.at(domain_id) = iga::AllocateCollectiveRuntime<iga::OneDPetscSolverContext>(communicator,
				communicator, prefix, application_options);
		}
		iga::CollectiveLocalStage(communicator, "graph 1D initialization", [&] {
			const bool species_mode = configuration_holder->schema_version == 6;
			for (const auto& domain_id : plan_holder->domain_order) {
				if (configuration_holder->graph.Domain(domain_id).kind
					!= iga::DomainKind::OneDFlow) continue;
				const auto& domain = iga::GraphDomainDefinitionFor(*configuration_holder,
					domain_id);
				const bool has_incoming = std::any_of(plan_holder->interfaces.begin(),
					plan_holder->interfaces.end(), [&](const iga::PressureFlowInterfacePlan& edge) {
						return edge.flow_receiver.domain_id == domain_id;
					});
				if ((!has_incoming
					&& domain.one_d_inlet_policy != iga::OneDInletPolicy::ConfiguredOpenLoop)
					|| (has_incoming
						&& domain.one_d_inlet_policy != iga::OneDInletPolicy::CoupledRoot))
					throw std::runtime_error(
						"1D inlet policy does not match its directed graph role");
				one_d.emplace(domain_id,
					BuildOneD(*configuration_holder, assets, domain_id, species_mode, communicator,
						input_texts.at("domain configuration "+domain_id),
						checkpoint_enabled ? iga::GraphCheckpointDomainIdentity(checkpoint_identity, domain_id, "one-d") : "",
						one_d_solvers.at(domain_id).get()));
			}
		});

		auto& configuration = *configuration_holder;
		const auto& plan = *plan_holder;
		iga::CollectiveLocalStage(communicator, "graph domain preflight", [&] {
			for (const auto& item : one_d)
				for (const auto& port : configuration.graph.Domain(item.first).ports)
					if (port.locator_kind == "runtime_port")
						ValidateOneDPort(item.second, port);
			for (const auto& domain_id : plan.domain_order)
				if (configuration.graph.Domain(domain_id).kind
					== iga::DomainKind::ThreeDBodyFittedFlow)
					three_d.emplace(domain_id, BuildThreeDPreflight(
						configuration, plan, assets, domain_id, rank,
						configuration.schema_version == 6,
						input_texts.at("domain configuration "+domain_id)));
			for (const auto& domain_id : plan.domain_order)
				if (configuration.graph.Domain(domain_id).kind
					== iga::DomainKind::ThreeDImmersedFlow) {
					if (configuration.schema_version == 6)
						throw std::runtime_error(
							"schema-v6 transport is unsupported by the immersed flow backend");
					immersed.emplace(domain_id, mpi_size == 1
						? iga::ImmersedFlowCase::Load(assets.at(domain_id).case_directory,domain_id,
							configuration.graph.Domain(domain_id).ports,mpi_size)
						: iga::ImmersedFlowCase::Preflight(assets.at(domain_id).case_directory,domain_id,
							configuration.graph.Domain(domain_id).ports));
				}
			for (const auto& volume : three_d)
				for (const auto& line : one_d)
					if (line.second.runtime->FlowSystem().density
							!= iga::FirstNavierStokesSystem(volume.second->configuration).density
						|| line.second.runtime->FlowSystem().dynamic_viscosity
							!= iga::FirstNavierStokesSystem(volume.second->configuration).viscosity)
						throw std::runtime_error(
							"multidomain flow requires identical density and viscosity");
			for (const auto& volume : immersed)
				if (volume.second->IsTransient()
					&& (volume.second->Configuration().time.dt != configuration.time.dt_s
						|| volume.second->Configuration().time.steps != configuration.time.steps))
					throw std::runtime_error("transient immersed domain and graph time grids must match");
			for (const auto& volume : immersed)
				for (const auto& line : one_d)
					if (line.second.runtime->FlowSystem().density
							!= volume.second->RuntimeParameters().density
						|| line.second.runtime->FlowSystem().dynamic_viscosity
							!= volume.second->RuntimeParameters().dynamic_viscosity)
						throw std::runtime_error(
							"multidomain flow requires identical density and viscosity");
		});
		for (const auto& domain_id : plan.domain_order) if (immersed.count(domain_id)) {
			auto& native = *immersed.at(domain_id);
			if (mpi_size > 1) native.InitializeDistributed(communicator);
			if (!native.IsDistributed()) continue;
			iga::CollectiveLocalStage(communicator,"immersed distribution diagnostics",[&] {
				const char* profile = std::getenv("IGA_PROFILE");
				if (!iga::CurrentPhaseProfile().Enabled() && (!profile || profile[0] != '1' || profile[1] != '\0')) return;
				const auto distribution = native.Distribution();
				std::cout << "hpc_immersed_distribution {\"domain\":\"" << JsonEscape(domain_id)
					<< "\",\"rank\":" << rank << ",\"ranks\":" << mpi_size
					<< ",\"time_integration\":\"" << (native.IsTransient() ? "backward_euler" : "steady") << "\""
					<< ",\"global_rows\":" << distribution.global_rows
					<< ",\"owned_rows\":" << distribution.owned_rows
					<< ",\"owned_stencils\":" << distribution.owned_stencils
					<< ",\"required_rows\":" << distribution.required_rows << "}\n";
				iga::FlushCheckedText(std::cout);
			});
		}
		for (const auto& domain_id : plan.domain_order) {
			if (!three_d.count(domain_id)) continue;
			auto& native = *three_d.at(domain_id);
			const auto& flow = iga::FirstNavierStokesSystem(native.configuration);
			std::string flow_identity, transport_identity, flow_prefix, transport_prefix;
			iga::CollectiveLocalStage(communicator, "graph checkpoint domain identities", [&] {
				flow_prefix = iga::PetscDomainOptionsPrefix(domain_id, "flow");
				transport_prefix = iga::PetscDomainOptionsPrefix(domain_id, "transport");
				if (checkpoint_enabled) {
					flow_identity = iga::GraphCheckpointDomainIdentity(checkpoint_identity, domain_id, "flow");
					transport_identity = iga::GraphCheckpointDomainIdentity(checkpoint_identity, domain_id, "transport");
				}
			});
			native.runtime = iga::AllocateCollectiveRuntime<iga::TransientFlowRuntime>(communicator, *native.database,
				communicator, true, true,
				iga::NavierStokesParameters{flow.density, flow.viscosity,
					native.configuration.time.dt}, native.initial_boundaries,
				native.mesh.labels, native.boundary_velocity, native.wall_trace_basis,
				std::move(native.outlet_models), application_options, flow_identity, configuration.time.dt_s, flow_prefix);
			iga::RequireValidGeometry(native.runtime->Elements(), rank, communicator);
			for (const auto& port : configuration.graph.Domain(domain_id).ports) {
				int label = 0;
				iga::RuntimeConstructionStage(communicator, "graph boundary label", [&] {
					label = iga::ParseThreeDFlowBoundaryLabel(port);
				});
				long long local_faces = 0;
				for (const auto& element : native.runtime->OwnedElements())
					local_faces += static_cast<long long>(std::count(
						element.boundary_labels.begin(), element.boundary_labels.end(), label));
				long long global_faces = 0;
				MPI_Allreduce(&local_faces, &global_faces, 1, MPI_LONG_LONG, MPI_SUM,
					communicator);
				iga::RuntimeConstructionStage(communicator, "graph boundary coverage", [&] {
					if (global_faces == 0)
						throw std::runtime_error("3D logical port has no database boundary face");
				});
				if (!port.requires.count(iga::PortQuantity::FlowRate)) continue;
				const double reference = port.orientation.ToOutward(
					native.runtime->ReferenceBoundaryFlow(label));
				iga::RuntimeConstructionStage(communicator, "graph initial port catalog", [&] {
					native.reference_outward_flow_m3_s.emplace(port.id, reference);
					const auto interface = std::find_if(plan.interfaces.begin(), plan.interfaces.end(),
						[&](const iga::PressureFlowInterfacePlan& edge) {
							return edge.flow_receiver.domain_id == domain_id
								&& edge.flow_receiver.port_id == port.id;
						});
					if (interface == plan.interfaces.end())
						throw std::runtime_error("3D flow receiver is not bound to a graph edge");
					iga::PortBoundaryData input;
					input.time_s = 0.0;
					if (one_d.count(interface->flow_provider.domain_id)) {
						const auto& provider = one_d.at(interface->flow_provider.domain_id);
						const auto provider_state = provider.runtime->GetPortState(
							configuration.graph.Port(interface->flow_provider).locator);
						input.outward_flow_m3_s = -RequireValue(
							provider_state.outward_flow_m3_s, "initial upstream provider flow");
					} else {
						const auto& provider = iga::GraphDomainDefinitionFor(configuration,
							interface->flow_provider.domain_id);
						if (provider.kind != iga::DomainKind::ZeroDFlow
							|| !provider.zero_d_flow_model
							|| provider.zero_d_flow_model->model.role
								!= iga::ZeroDFlowRole::SourceReservoir)
							throw std::runtime_error(
								"flow-controlled 3D inlet requires a 1D or source-0D flow provider");
						const auto& source = provider.zero_d_flow_model->model.source;
						const double edge_pressure = configuration.initial_pressure_pa.at(
							interface->edge_id);
						// This is an initialization-only boundary seed.  It uses the
						// committed compliant pressure and coupling-edge pressure directly;
						// do not publish or commit a speculative 0D trial here.
						input.outward_flow_m3_s = -(provider.zero_d_flow_model->initial_state
							.stored_pressure_pa-edge_pressure)/source.resistance_pa_s_m3;
					}
					iga::ApplyThreeDReferenceProfileInput(native.initial_configuration,
						iga::FirstNavierStokesSystem(native.initial_configuration), port, input,
						reference);
				});
			}
			native.runtime->InitializeState(native.initial_configuration);
			if (configuration.schema_version == 6) {
				iga::RuntimeConstructionStage(communicator, "graph transport preflight", [&] {
					if (!native.transport_system)
						throw std::runtime_error("schema-v6 3D transport system was not preflighted");
				});
				native.transport_runtime = iga::AllocateCollectiveRuntime<iga::TransientTransportRuntime>(communicator,
					*native.database, communicator, native.configuration,
					*native.transport_system, native.mesh.labels,
					std::map<std::uint64_t, iga::VolumeQuadratureRule>{}, application_options, transport_identity, transport_prefix);
				iga::RuntimeConstructionStage(communicator, "graph transport node agreement", [&] {
					if (native.runtime->RequiredNodes() != native.transport_runtime->RequiredNodes())
						throw std::runtime_error(
							"schema-v6 3D flow and transport runtimes require identical ordered nodes");
				});
			}
		}

		std::vector<std::unique_ptr<iga::CoupledDomainRuntime>> runtimes;
		std::map<std::string, iga::ZeroDFlowDomainRuntime*> zero_d;
		iga::ThreeDFlowDomainControls controls;
		controls.maximum_newton = options.maximum_newton;
		controls.nonlinear_relative_tolerance = kNonlinearRelativeTolerance;
		controls.nonlinear_absolute_tolerance = kNonlinearAbsoluteTolerance;
		controls.mass_relative_tolerance = kMassRelativeTolerance;
		iga::RuntimeConstructionStage(communicator, "graph adapter catalog", [&] {
			runtimes.reserve(plan.domain_order.size());
			for (const auto& domain_id : plan.domain_order) {
				if (one_d.count(domain_id)) {
					auto& native = one_d.at(domain_id);
					if (configuration.schema_version == 6)
					{
						iga::OneDFlowTransportDomainControls species_controls;
						species_controls.species_flow_epsilon_m3_s =
							configuration.execution.species_routing->flow_switch_m3_s;
						runtimes.push_back(std::make_unique<iga::OneDFlowTransportDomainAdapter>(
							domain_id, *native.runtime, configuration.graph.Domain(domain_id).ports,
							native.inlet_policy,
							configuration.graph.Domain(domain_id).species_bindings,
							species_controls));
					}
					else runtimes.push_back(std::make_unique<iga::OneDFlowDomainAdapter>(domain_id,
						*native.runtime, configuration.graph.Domain(domain_id).ports,
						native.inlet_policy));
				} else if (three_d.count(domain_id)) {
					auto& native = *three_d.at(domain_id);
					if (configuration.schema_version == 6) {
						iga::ThreeDFlowTransportDomainControls species_controls;
						species_controls.flow = controls;
						species_controls.species_flow_epsilon_m3_s =
							configuration.execution.species_routing->flow_switch_m3_s;
						runtimes.push_back(
							std::make_unique<iga::ThreeDBodyFittedFlowTransportDomainAdapter>(
								domain_id, *native.runtime, *native.transport_runtime,
								configuration.graph.Domain(domain_id).ports, native.configuration,
								native.case_directory,
								configuration.graph.Domain(domain_id).species_bindings,
								native.reference_outward_flow_m3_s, species_controls));
					} else runtimes.push_back(std::make_unique<iga::ThreeDBodyFittedFlowDomainAdapter>(
						domain_id, *native.runtime, configuration.graph.Domain(domain_id).ports,
						native.configuration, native.case_directory,
						native.reference_outward_flow_m3_s, controls));
				} else if (immersed.count(domain_id)) {
					// Allocate the audit entry before transferring the complete owner.
					immersed_audit.emplace(domain_id, immersed.at(domain_id).get());
					runtimes.push_back(std::move(immersed.at(domain_id)));
				} else {
					const auto& definition = iga::GraphDomainDefinitionFor(configuration, domain_id);
					if (definition.kind != iga::DomainKind::ZeroDFlow
						|| !definition.zero_d_flow_model)
						throw std::runtime_error("multidomain graph has no runtime owner for domain '"
							+domain_id+"'");
					auto runtime = std::make_unique<iga::ZeroDFlowDomainRuntime>(domain_id,
						definition.zero_d_flow_model->model,
						definition.zero_d_flow_model->initial_state,
						configuration.graph.Domain(domain_id).ports);
					zero_d.emplace(domain_id, runtime.get());
					runtimes.push_back(std::move(runtime));
				}
			}
		});
		auto registry_owner = iga::CreateCollectiveDomainRuntimeRegistry(communicator,
			configuration.graph, runtimes);
		auto& registry = *registry_owner;
		std::unique_ptr<iga::PressureFlowComponentExecutor> flow_executor;
		std::unique_ptr<iga::SpeciesPressureFlowComponentExecutor> species_executor;
		iga::CollectiveLocalStage(communicator, "graph executor construction", [&] {
			if (configuration.schema_version == 6)
				species_executor = std::make_unique<iga::SpeciesPressureFlowComponentExecutor>(registry,
					configuration.start_domain_id, iga::SpeciesPressureFlowControlsFor(configuration),
					iga::CollectiveSpeciesPressureFlowExecution(communicator));
			else flow_executor = std::make_unique<iga::PressureFlowComponentExecutor>(registry,
				configuration.start_domain_id, iga::PressureFlowControlsFor(configuration.execution),
					iga::CollectivePressureFlowExecution(communicator));
		});

		const int final_step = options.stop_after_step > 0
			? options.stop_after_step : configuration.time.steps;
		std::map<std::string, double> pressure;
		std::set<iga::PortRef> coupled_ports;
		iga::CollectiveLocalStage(communicator, "graph step catalog", [&] {
			if (final_step > configuration.time.steps)
				throw std::runtime_error("--stop-after-step exceeds configured steps");
			pressure = configuration.initial_pressure_pa;
			for (const auto& edge : configuration.graph.Edges()) {
				coupled_ports.insert(edge.first);
				coupled_ports.insert(edge.second);
			}
		});
		auto VerifyHydraulicBalance = [&](const auto& trial,
			std::map<std::string, std::pair<double, double>>& pending_balance,
			double& pending_mass, double& pending_external) {
			for (const auto& volume : configuration.graph.Domains()) {
				if (iga::DomainDimensionOf(volume.second.kind) != 3) continue;
				double sum = 0.0;
				double absolute_sum = 0.0;
				for (const auto& port : volume.second.ports) {
					const double flow = RequireValue(trial.accepted_ports.at(
						{volume.first, port.id}).outward_flow_m3_s, "3D boundary flow");
					sum += flow;
					absolute_sum += std::abs(flow);
				}
				if (absolute_sum > 0.0
					&& 2.0*std::abs(sum)/absolute_sum > kMassRelativeTolerance)
					throw std::runtime_error("3D domain balance exceeds coupling tolerance");
				pending_balance.emplace(volume.first, std::make_pair(sum, absolute_sum));
				pending_mass += sum;
			}
			for (const auto& domain : configuration.graph.Domains())
				for (const auto& port : domain.second.ports) {
					const iga::PortRef reference{domain.first, port.id};
					if (coupled_ports.count(reference)
						|| !port.provides.count(iga::PortQuantity::FlowRate)) continue;
					pending_external += RequireValue(trial.accepted_ports.at(reference)
						.outward_flow_m3_s, "external boundary flow");
				}
		};
		std::vector<AcceptedStep> accepted;
		std::vector<AcceptedSpeciesStep> accepted_species;
		// The start of every accepted step is the exact EndTime() that was
		// committed by the preceding context.  Do not reconstruct it as n*dt:
		// exact-clock runtimes deliberately reject that different rounding.
		double accepted_time_s = 0.0;
		long long first_step = 1;
		std::string previous_checkpoint;
		std::unique_ptr<iga::NativeGraphCheckpoint> checkpoint;
		if (checkpoint_enabled) {
			iga::NativeGraphCheckpointOwners owners;
			iga::CollectiveLocalStage(communicator, "graph checkpoint owner catalog", [&] {
				for (const auto& item : one_d) owners.one_d.emplace(item.first, item.second.runtime.get());
				for (const auto& item : three_d) {
					owners.flow.emplace(item.first, item.second->runtime.get());
					if (item.second->transport_runtime) owners.transport.emplace(item.first, item.second->transport_runtime.get());
				}
				owners.zero_d = zero_d;
			});
			checkpoint = iga::AllocateCollectiveRuntime<iga::NativeGraphCheckpoint>(communicator, communicator,
				configuration.graph, checkpoint_identity, std::move(owners), bool(species_executor));
			if (!options.restart_directory.empty()) {
				iga::PressureFlowCheckpointControls restored;
				const auto epoch = checkpoint->Restore(options.restart_directory, final_step, configuration.time.dt_s,
					restored, accepted, accepted_species, species_executor.get());
				iga::CollectiveLocalStage(communicator, "graph checkpoint activation", [&] {
					pressure = std::move(restored.next_pressure_pa); accepted_time_s = epoch.time_s;
					first_step = static_cast<long long>(epoch.accepted_steps)+1; previous_checkpoint = epoch.id;
					if (rank == 0) { std::cout << "restored graph checkpoint step=" << epoch.accepted_steps << " epoch=" << epoch.id << '\n'; iga::FlushCheckedText(std::cout); }
				});
			}
		}
		for (long long step_number = first_step; step_number <= final_step; ++step_number) {
			const int step = static_cast<int>(step_number);
			const iga::DomainStepContext step_context{step-1, accepted_time_s,
				configuration.time.dt_s};
			const double time = step_context.EndTime();
			std::map<std::string, std::pair<double, double>> pending_balance;
			double pending_mass = 0.0;
			double pending_external = 0.0;
			if (species_executor) {
				std::function<void(const iga::SpeciesPressureFlowStepResult&)> before_commit;
				iga::CollectiveLocalStage(communicator, "graph species precommit callback", [&] {
					before_commit = [&](const iga::SpeciesPressureFlowStepResult& trial) {
						VerifyHydraulicBalance(trial, pending_balance, pending_mass, pending_external);
						if (step == injected_failure_step)
							throw std::runtime_error("injected bifurcation failure before commit");
					};
				});
				auto result = species_executor->Advance(step_context, pressure, before_commit);
				iga::CollectiveLocalStage(communicator, "graph accepted species result", [&] {
					auto transport_ports = std::move(result.transport_ports);
					accepted_species.push_back({step, time,
						static_cast<int>(result.hydraulic_iterations.size()), pending_mass,
						pending_external, std::move(pending_balance),
						std::move(transport_ports), std::move(result)});
					for (const auto& edge : accepted_species.back().result.hydraulic_iterations.back().edges)
						pressure[edge.edge_id] = edge.measured_pressure_pa;
				});
			} else {
				std::function<void(const iga::PressureFlowStepResult&)> before_commit;
				iga::CollectiveLocalStage(communicator, "graph precommit callback", [&] {
					before_commit = [&](const iga::PressureFlowStepResult& trial) {
						VerifyHydraulicBalance(trial, pending_balance, pending_mass,
							pending_external);
						if (step == injected_failure_step)
							throw std::runtime_error("injected bifurcation failure before commit");
					};
				});
				auto result = flow_executor->Advance(step_context, pressure, before_commit);
				iga::CollectiveLocalStage(communicator, "graph accepted flow result", [&] {
					AcceptedStep accepted_step{step, time, static_cast<int>(result.iterations.size()),
						pending_mass, pending_external, std::move(pending_balance), {}, {},
						std::move(result)};
					for (const auto& zero : zero_d) {
						accepted_step.zero_d_states.emplace(zero.first, zero.second->CommittedState());
						accepted_step.zero_d_accounting.emplace(zero.first,
							*zero.second->CommittedStepAccounting());
					}
					accepted.push_back(std::move(accepted_step));
					for (const auto& edge : accepted.back().result.iterations.back().edges)
						pressure[edge.edge_id] = edge.measured_pressure_pa;
				});
			}
			for (const auto& entry : immersed_audit) if (entry.second->IsTransient()) {
				if (entry.second->IsMoving()) {
					const auto& runtime = entry.second->MovingRuntime();
					const auto conservation = runtime.ConservationDiagnostics();
					iga::CollectiveLocalStage(communicator,"immersed accepted moving diagnostics",[&] {
						const auto& clock=runtime.Clock();const auto& diagnostic=runtime.Diagnostics();
						if (clock.time_s!=time || clock.index!=static_cast<std::uint64_t>(step) || clock.trial_active
							|| diagnostic.commit_count!=1 || conservation.target_index!=clock.index)
							throw std::logic_error("immersed moving graph accepted backend clock differs");
						const char* profile=std::getenv("IGA_PROFILE");
						if (!iga::CurrentPhaseProfile().Enabled() && (!profile || profile[0]!='1' || profile[1]!='\0')) return;
						double minimum_damping=1.;
						for(const auto& solve_step:diagnostic.newton_steps)if(solve_step.damping>0)minimum_damping=std::min(minimum_damping,solve_step.damping);
						std::cout<<std::setprecision(17)<<"hpc_immersed_moving_step {\"domain\":\""<<JsonEscape(entry.first)
							<<"\",\"rank\":"<<rank<<",\"ranks\":"<<mpi_size<<",\"time_s\":"<<clock.time_s<<",\"index\":"<<clock.index
							<<",\"residual_norm\":"<<diagnostic.residual_norm
							<<",\"newton_iterations\":"<<diagnostic.nonlinear_iterations<<",\"minimum_accepted_damping\":"<<minimum_damping
							<<",\"assembly_s\":"<<diagnostic.aggregate_assembly_seconds<<",\"solve_s\":"<<diagnostic.aggregate_linear_solve_seconds
							<<",\"reynolds\":"<<conservation.normalized_reynolds_defect
							<<",\"moving_mass\":"<<conservation.normalized_moving_mass_defect
							<<",\"wall_relative_leakage\":"<<conservation.normalized_wall_relative_leakage<<"}\n";
						iga::FlushCheckedText(std::cout);
					});
					continue;
				}
				const auto& runtime = entry.second->TransientRuntime();
				const auto conservation = runtime.ConservationDiagnostics();
				iga::CollectiveLocalStage(communicator,"immersed accepted transient diagnostics",[&] {
					const auto& clock = runtime.Clock(); const auto& diagnostic = runtime.Diagnostics();
					if (clock.time_s != time || clock.index != static_cast<std::uint64_t>(step)
						|| clock.trial_active || runtime.History().Active() || diagnostic.commit_count != static_cast<std::size_t>(step))
						throw std::logic_error("immersed graph accepted backend clock differs");
					const char* profile = std::getenv("IGA_PROFILE");
					if (!iga::CurrentPhaseProfile().Enabled() && (!profile || profile[0] != '1' || profile[1] != '\0')) return;
					std::cout << std::setprecision(17) << "hpc_immersed_transient_step {\"domain\":\"" << JsonEscape(entry.first)
						<< "\",\"rank\":" << rank << ",\"ranks\":" << mpi_size << ",\"time_s\":" << clock.time_s
						<< ",\"index\":" << clock.index << ",\"commits\":" << diagnostic.commit_count
						<< ",\"residual_norm\":" << diagnostic.residual_norm
						<< ",\"assembly_s\":" << diagnostic.aggregate_assembly_seconds
						<< ",\"solve_s\":" << diagnostic.aggregate_linear_solve_seconds
						<< ",\"surface_flow\":" << conservation.total_surface_outward_flow_m3_s
						<< ",\"volume_divergence\":" << conservation.volume_divergence_integral_m3_s
						<< ",\"wall_flow\":" << conservation.wall_outward_flow_m3_s << "}\n";
					iga::FlushCheckedText(std::cout);
				});
			}
			// Advance only after the executor's transactional commit and all
			// accepted-result bookkeeping have completed.  A failed trial or
			// precommit callback therefore leaves this accepted clock unchanged.
			accepted_time_s = step_context.EndTime();
			iga::CollectiveLocalStage(communicator, "graph solver configuration output", [&] {
				if (rank != 0) return;
				for (const auto& domain : one_d_solvers)
					iga::WriteOneDPetscSolverConfiguration(std::cout, *domain.second, step);
				const auto write = [&](const std::string& id, const char* role, const iga::PetscKspConfiguration& solver) {
					const auto precision = std::cout.precision();
					std::cout << std::setprecision(17) << "solver_configuration {\"domain\":\"" << JsonEscape(id)
						<< "\",\"role\":\"" << role << "\",\"step\":" << step
						<< ",\"prefix\":\"" << JsonEscape(solver.prefix) << "\",\"ksp\":\"" << JsonEscape(solver.ksp)
						<< "\",\"pc\":\"" << JsonEscape(solver.pc) << "\",\"factor_backend\":\"" << JsonEscape(solver.factor_backend)
						<< "\",\"rtol\":" << solver.relative_tolerance << ",\"atol\":" << solver.absolute_tolerance
						<< ",\"max_iterations\":" << solver.maximum_iterations << ",\"last_iterations\":" << solver.last_iterations
						<< ",\"last_reason\":" << static_cast<int>(solver.last_reason) << "}\n";
					std::cout.precision(precision);
				};
				for (const auto& domain : three_d) {
					write(domain.first, "flow", domain.second->runtime->SolverConfiguration());
					if (domain.second->transport_runtime) write(domain.first, "transport", domain.second->transport_runtime->SolverConfiguration());
				}
				for (const auto& domain : immersed_audit) write(domain.first, "flow", domain.second->SolverConfiguration());
				iga::FlushCheckedText(std::cout);
			});

			int stop_requested = 0;
			if (checkpoint_signal) {
				const int local_request = checkpoint_signal->Requested();
				MPI_Allreduce(&local_request, &stop_requested, 1, MPI_INT, MPI_MAX, communicator);
			}
			if (checkpoint && !options.checkpoint_directory.empty() && (step % options.checkpoint_every == 0 || step == final_step || stop_requested)) {
				iga::PressureFlowCheckpointControls saved;
				iga::CollectiveLocalStage(communicator, "graph checkpoint accepted controls", [&] {
					saved = species_executor
						? iga::MakePressureFlowCheckpointControls(step_context, accepted_species.back().result.hydraulic_iterations.back(), species_executor->CaptureCheckpointDonors())
						: iga::MakePressureFlowCheckpointControls(step_context, accepted.back().result.iterations.back());
				});
				const auto epoch = checkpoint->Save(options.checkpoint_directory, step_context, saved, accepted, accepted_species, previous_checkpoint);
				iga::CollectiveLocalStage(communicator, "graph checkpoint completion", [&] {
					previous_checkpoint = epoch.id;
					if (rank == 0) {
						std::cout << "saved graph checkpoint step=" << epoch.accepted_steps << " epoch=" << epoch.id << '\n';
						if (stop_requested) std::cout << "stopped after checkpoint request step=" << step << '\n';
						iga::FlushCheckedText(std::cout);
					}
				});
			}
			if (stop_requested) break;
		}

		for (auto& entry : three_d) {
			auto& native = *entry.second;
			if (native.transport_runtime) native.transport_runtime->Close();
			native.runtime->Close();
		}
		for (const auto& entry : immersed_audit) entry.second->CloseDistributed();
		iga::CollectiveLocalStage(communicator, "graph output", [&] {
			if (rank != 0) return;
			if (species_executor)
				WriteSpeciesOutputs(options.output_directory, configuration, graph_root, assets,
					one_d, three_d, accepted_species);
			else WriteOutputs(options.output_directory, configuration, graph_root, assets,
				bifurcation_holder, one_d, immersed_audit, accepted);
		});
		iga::CollectiveLocalStage(communicator, "graph completion logging", [&] {
			if (rank != 0) return;
#ifdef IGA_BIFURCATION_ENTRY
			std::cout << "completed schema-v5 1D--3D bifurcation steps="
				<< accepted.size() << " branches=" << bifurcation_holder->branches.size()
				<< " output=" << options.output_directory << '\n';
#else
			std::cout << "completed schema-v" << configuration.schema_version
				<< " multidomain " << (species_executor ? "flow/transport" : "flow")
				<< " steps=" << (species_executor ? accepted_species.size() : accepted.size())
				<< " domains=" << configuration.graph.Domains().size()
				<< " output=" << options.output_directory << '\n';
#endif
			iga::FlushCheckedText(std::cout);
		});
	} catch (const std::exception& error) {
		if (rank == 0) std::cerr << error.what() << '\n';
		status = 1;
	}
	return status;
}

#ifndef IGA_MULTIDOMAIN_NO_MAIN
int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr,
		"TubularFlowIGA multidomain coupling\n");
	iga::CurrentPhaseProfile().EnableFromEnvironment();
	int rank = 0, ranks = 1; MPI_Comm_rank(PETSC_COMM_WORLD, &rank); MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
	int status = iga::RunMultidomainFlow(argc, argv, PETSC_COMM_WORLD);
	try {
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "graph profile output", [&] {
			iga::CurrentPhaseProfile().Write(std::cout, rank, ranks, status); iga::FlushCheckedText(std::cout);
		});
	} catch (const std::exception& error) { if (rank == 0) std::cerr << error.what() << '\n'; status = 1; }
	PetscFinalize();
	return status;
}
#endif
