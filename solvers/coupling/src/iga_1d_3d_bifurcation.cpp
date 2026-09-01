#include "BifurcationMultidomainCase.hpp"
#include "BoundarySupport.hpp"
#include "DomainRuntimeRegistry.hpp"
#include "ExplicitOneDThreeDCoupling.hpp"
#include "IgaDatabase.hpp"
#include "OneDFlowDomainAdapter.hpp"
#include "OneDImplicit.hpp"
#include "OneDRuntime.hpp"
#include "PressureFlowComponentExecutor.hpp"
#include "ThreeDBodyFittedFlowDomainAdapter.hpp"
#include "ThreeDFlowCoupling.hpp"
#include "ThreeDVcaCoupling.hpp"
#include "TransientFlowRuntime.hpp"

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
			|| argument == "--stop-after-step" || argument == "--three-d-max-newton") {
			if (++i >= argc) throw std::runtime_error(argument+" requires a value");
			if (argument == "--graph-case") options.graph_case = argv[i];
			else if (argument == "--output-dir") options.output_directory = argv[i];
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
			"--output-dir DIR [--stop-after-step N] [--three-d-max-newton N] [PETSc options]");
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
};

NativeOneD BuildOneD(const iga::MultidomainConfiguration& graph,
	const std::map<std::string, iga::ResolvedGraphDomainAssets>& assets,
	const std::string& domain_id)
{
	const auto& definition = iga::GraphDomainDefinitionFor(graph, domain_id);
	const auto case_directory = assets.at(domain_id).case_directory;
	auto configuration = iga::ParseOneDConfiguration(ReadText(
		iga::ResolveContainedCaseFile(case_directory, "simulation_config.json",
			"1D configuration")));
	CanonicalizePeriodicTableAssets(configuration.temporal_functions, case_directory, domain_id);
	RequireOneDFlowOnly(configuration);
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
		[](const iga::OneDNetwork& network_definition,
			const iga::OneDFlowSystemDefinition& flow_definition,
			iga::OneDFlowState& state, double inlet_flow, double dt) {
			iga::AdvanceImplicitOneD(network_definition, flow_definition, state, inlet_flow, dt);
		});
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
	const std::string& domain_id, int rank)
{
	auto native = std::make_unique<NativeThreeD>();
	native->domain_id = domain_id;
	native->case_directory = assets.at(domain_id).case_directory;
	native->database = std::make_unique<iga::Database>(assets.at(domain_id).database.string());
	const auto required_elements = native->database->LoadRequired(rank);
	const auto owned_elements = native->database->LoadOwned(rank);
	native->configuration = iga::ReadSimulationConfiguration(
		iga::ResolveContainedCaseFile(native->case_directory, "simulation_config.json",
			domain_id+" 3D configuration").string());
	CanonicalizePeriodicTableAssets(native->configuration.temporal_functions,
		native->case_directory, domain_id);
	if (native->configuration.equation_systems.size() != 1
		|| native->configuration.equation_systems.front().kind
			!= iga::EquationKind::NavierStokes
		|| native->configuration.equation_systems.front().unknowns.size() != 2
		|| native->configuration.coupling.mode != iga::SimulationScopeMode::FlowOnly
		|| native->configuration.physiology.enabled
		|| native->configuration.equation_systems.front().time_integration
			!= "backward_euler")
		throw std::runtime_error(
			"multidomain flow requires one flow_only backward_euler 3D Navier-Stokes system");
	for (const auto& field : native->configuration.fields)
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
	native->outlet_models = iga::InitializeOutletModels(native->configuration,
		native->configuration.equation_systems.front());
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
				native->configuration.equation_systems.front().unknowns.at(1), label,
				graph.initial_pressure_pa.at(interface->edge_id));
		}
		if (port.requires.count(iga::PortQuantity::FlowRate)) {
			auto check = native->initial_configuration;
			iga::PortBoundaryData input;
			input.time_s = 0.0;
			input.outward_flow_m3_s = 1.0;
			iga::ApplyThreeDReferenceProfileInput(check,
				check.equation_systems.front(), port, input, 1.0);
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
		native->initial_configuration.equation_systems.front(), native->mesh.labels,
		native->boundary_velocity);
	return native;
}

struct AcceptedStep {
	int step = 0;
	double time_s = 0.0;
	int iterations = 0;
	double three_d_mass_m3_s = 0.0;
	double external_outward_flow_m3_s = 0.0;
	std::map<std::string, std::pair<double, double>> three_d_balance;
	iga::PressureFlowStepResult result;
};

void WriteOutputs(const fs::path& directory,
	const iga::MultidomainConfiguration& configuration,
	const fs::path& graph_root,
	const std::map<std::string, iga::ResolvedGraphDomainAssets>& assets,
	const std::optional<iga::OneDThreeDBifurcationDefinition>& bifurcation,
	const std::map<std::string, NativeOneD>& one_d,
	const std::map<std::string, std::unique_ptr<NativeThreeD>>& three_d,
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
	if (!summary || !edges || !iterations || !ports || !initialization || !balances)
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
				<< CsvEscape(port.first.domain_id) << ',' << CsvEscape(port.first.port_id) << ','
				<< RequireValue(port.second.area_m2, "port area") << ','
				<< RequireValue(port.second.outward_flow_m3_s, "port flow") << ',';
			if (port.second.mean_pressure_pa) ports << *port.second.mean_pressure_pa;
			ports << '\n';
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
	if (!summary || !edges || !iterations || !ports || !initialization || !balances)
		throw std::runtime_error("cannot finalize bifurcation long-form output");
	std::ofstream marker(directory/"graph_binding_manifest.json.tmp");
	if (!marker) throw std::runtime_error("cannot create bifurcation completion marker");
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
			<< "  \"three_d_domain_count\": " << three_d.size() << ",\n"
			<< "  \"completed_steps\": " << steps.size() << ",\n  \"domains\": [\n";
		std::size_t domain_index = 0;
		for (const auto& domain : configuration.graph.Domains()) {
			const auto& resolved = assets.at(domain.first);
			marker << "    {\"id\":\"" << JsonEscape(domain.first) << "\",\"kind\":\""
				<< (domain.second.kind == iga::DomainKind::OneDFlow
					? "network_flow" : "body_fitted_iga_flow") << "\",\"case\":\""
				<< JsonEscape(resolved.case_directory.generic_string()) << "\"";
			if (domain.second.kind == iga::DomainKind::ThreeDBodyFittedFlow)
				marker << ",\"database\":\""
					<< JsonEscape(resolved.database.generic_string()) << "\"";
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

int FailureInjectionStep()
{
	const char* value = std::getenv("TUBULARFLOWIGA_INJECT_BIFURCATION_FAILURE_STEP");
	return value ? PositiveInteger(value,
		"TUBULARFLOWIGA_INJECT_BIFURCATION_FAILURE_STEP") : 0;
}

void RequireCollectivePreflight(const std::string& local_error)
{
	const int local_failed = local_error.empty() ? 0 : 1;
	int any_failed = 0;
	MPI_Allreduce(&local_failed, &any_failed, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
	if (!any_failed) return;
	if (!local_error.empty()) throw std::runtime_error(local_error);
	throw std::runtime_error("bifurcation graph preflight failed on another MPI rank");
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr,
		"TubularFlowIGA schema-v5 1D--3D bifurcation coupling\n");
	int rank = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	int status = 0;
	try {
		const auto options = ParseOptions(argc, argv);
		std::optional<iga::MultidomainConfiguration> configuration_holder;
		std::optional<iga::PressureFlowComponentPlan> plan_holder;
		std::optional<iga::OneDThreeDBifurcationDefinition> bifurcation_holder;
		std::map<std::string, iga::ResolvedGraphDomainAssets> assets;
		std::map<std::string, NativeOneD> one_d;
		std::map<std::string, std::unique_ptr<NativeThreeD>> three_d;
		fs::path graph_root;
		std::string local_error;
		try {
			if (fs::exists(options.output_directory))
				throw std::runtime_error("bifurcation output directory must not already exist");
			graph_root = iga::CanonicalGraphCaseRoot(options.graph_case);
			configuration_holder.emplace(iga::ReadMultidomainConfiguration(
				iga::ResolveContainedCaseFile(graph_root, "simulation_config.json",
					"schema-v5 manifest").string()));
			plan_holder.emplace(iga::MakeAcyclicPressureFlowPlan(
				configuration_holder->graph, configuration_holder->start_domain_id));
#ifdef IGA_BIFURCATION_ENTRY
			bifurcation_holder.emplace(
				iga::ResolveOneDThreeDBifurcation(*configuration_holder));
#endif
			assets = iga::ResolveGraphDomainAssets(*configuration_holder, graph_root);
			if (configuration_holder->graph.Domain(configuration_holder->start_domain_id).kind
				!= iga::DomainKind::OneDFlow)
				throw std::runtime_error("multidomain flow requires a 1D start domain");
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
					BuildOneD(*configuration_holder, assets, domain_id));
			}
		} catch (const std::exception& error) {
			local_error = error.what();
		}
		RequireCollectivePreflight(local_error);
		auto& configuration = *configuration_holder;
		const auto& plan = *plan_holder;
		local_error.clear();
		try {
			for (const auto& item : one_d)
				for (const auto& port : configuration.graph.Domain(item.first).ports)
					if (port.locator_kind == "runtime_port")
						ValidateOneDPort(item.second, port);
			for (const auto& domain_id : plan.domain_order)
				if (configuration.graph.Domain(domain_id).kind
					== iga::DomainKind::ThreeDBodyFittedFlow)
					three_d.emplace(domain_id, BuildThreeDPreflight(
						configuration, plan, assets, domain_id, rank));
			for (const auto& volume : three_d)
				for (const auto& line : one_d)
					if (line.second.runtime->FlowSystem().density
							!= volume.second->configuration.equation_systems.front().density
						|| line.second.runtime->FlowSystem().dynamic_viscosity
							!= volume.second->configuration.equation_systems.front().viscosity)
						throw std::runtime_error(
							"multidomain flow requires identical density and viscosity");
		} catch (const std::exception& error) {
			local_error = error.what();
		}
		RequireCollectivePreflight(local_error);
		for (const auto& domain_id : plan.domain_order) {
			if (!three_d.count(domain_id)) continue;
			auto& native = *three_d.at(domain_id);
			const auto& flow = native.configuration.equation_systems.front();
			native.runtime = std::make_unique<iga::TransientFlowRuntime>(*native.database,
				PETSC_COMM_WORLD, true, true,
				iga::NavierStokesParameters{flow.density, flow.viscosity,
					native.configuration.time.dt}, native.initial_boundaries,
				native.mesh.labels, native.boundary_velocity, native.wall_trace_basis,
				std::move(native.outlet_models));
			iga::RequireValidGeometry(native.runtime->Elements(), rank, PETSC_COMM_WORLD);
			for (const auto& port : configuration.graph.Domain(domain_id).ports) {
				const int label = iga::ParseThreeDFlowBoundaryLabel(port);
				long long local_faces = 0;
				for (const auto& element : native.runtime->OwnedElements())
					local_faces += static_cast<long long>(std::count(
						element.boundary_labels.begin(), element.boundary_labels.end(), label));
				long long global_faces = 0;
				MPI_Allreduce(&local_faces, &global_faces, 1, MPI_LONG_LONG, MPI_SUM,
					PETSC_COMM_WORLD);
				if (global_faces == 0)
					throw std::runtime_error("3D logical port has no database boundary face");
				if (!port.requires.count(iga::PortQuantity::FlowRate)) continue;
				const double reference = port.orientation.ToOutward(
					native.runtime->ReferenceBoundaryFlow(label));
				native.reference_outward_flow_m3_s.emplace(port.id, reference);
				const auto interface = std::find_if(plan.interfaces.begin(), plan.interfaces.end(),
					[&](const iga::PressureFlowInterfacePlan& edge) {
						return edge.flow_receiver == iga::PortRef{domain_id, port.id};
					});
				if (interface == plan.interfaces.end())
					throw std::runtime_error("3D flow receiver is not bound to a graph edge");
				const auto& provider = one_d.at(interface->flow_provider.domain_id);
				const auto provider_state = provider.runtime->GetPortState(
					configuration.graph.Port(interface->flow_provider).locator);
				iga::PortBoundaryData input;
				input.time_s = 0.0;
				input.outward_flow_m3_s = -RequireValue(
					provider_state.outward_flow_m3_s, "initial upstream provider flow");
				iga::ApplyThreeDReferenceProfileInput(native.initial_configuration,
					native.initial_configuration.equation_systems.front(), port, input,
					reference);
			}
			native.runtime->InitializeState(native.initial_configuration);
		}

		std::vector<std::unique_ptr<iga::CoupledDomainRuntime>> runtimes;
		iga::ThreeDFlowDomainControls controls;
		controls.maximum_newton = options.maximum_newton;
		controls.nonlinear_relative_tolerance = kNonlinearRelativeTolerance;
		controls.nonlinear_absolute_tolerance = kNonlinearAbsoluteTolerance;
		controls.mass_relative_tolerance = kMassRelativeTolerance;
		for (const auto& domain_id : plan.domain_order) {
			if (one_d.count(domain_id)) {
				auto& native = one_d.at(domain_id);
				runtimes.push_back(std::make_unique<iga::OneDFlowDomainAdapter>(domain_id,
					*native.runtime, configuration.graph.Domain(domain_id).ports,
					native.inlet_policy));
			} else {
				auto& native = *three_d.at(domain_id);
				runtimes.push_back(std::make_unique<iga::ThreeDBodyFittedFlowDomainAdapter>(
					domain_id, *native.runtime, configuration.graph.Domain(domain_id).ports,
					native.configuration, native.case_directory,
					native.reference_outward_flow_m3_s, controls));
			}
		}
		iga::DomainRuntimeRegistry registry(configuration.graph, std::move(runtimes));
		iga::PressureFlowComponentExecutor executor(registry,
			configuration.start_domain_id, iga::PressureFlowControlsFor(configuration.execution));

		const int final_step = options.stop_after_step > 0
			? options.stop_after_step : configuration.time.steps;
		if (final_step > configuration.time.steps)
			throw std::runtime_error("--stop-after-step exceeds configured steps");
		const int injected_failure_step = FailureInjectionStep();
		auto pressure = configuration.initial_pressure_pa;
		std::set<iga::PortRef> coupled_ports;
		for (const auto& edge : configuration.graph.Edges()) {
			coupled_ports.insert(edge.first);
			coupled_ports.insert(edge.second);
		}
		std::vector<AcceptedStep> accepted;
		for (int step = 1; step <= final_step; ++step) {
			const double time = step*configuration.time.dt_s;
			std::map<std::string, std::pair<double, double>> pending_balance;
			double pending_mass = 0.0;
			double pending_external = 0.0;
			auto result = executor.Advance({step-1, (step-1)*configuration.time.dt_s,
				configuration.time.dt_s}, pressure,
				[&](const iga::PressureFlowStepResult& trial) {
					for (const auto& volume : three_d) {
						double sum = 0.0;
						double absolute_sum = 0.0;
						for (const auto& port : configuration.graph.Domain(volume.first).ports) {
							const double flow = RequireValue(trial.accepted_ports.at(
								{volume.first, port.id}).outward_flow_m3_s, "3D boundary flow");
							sum += flow;
							absolute_sum += std::abs(flow);
						}
						if (absolute_sum > 0.0
							&& 2.0*std::abs(sum)/absolute_sum > kMassRelativeTolerance)
							throw std::runtime_error("3D domain balance exceeds coupling tolerance");
						pending_balance.emplace(volume.first,
							std::make_pair(sum, absolute_sum));
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
					if (step == injected_failure_step)
						throw std::runtime_error("injected bifurcation failure before commit");
				});
			accepted.push_back({step, time, static_cast<int>(result.iterations.size()),
				pending_mass, pending_external, std::move(pending_balance), std::move(result)});
			for (const auto& edge : accepted.back().result.iterations.back().edges)
				pressure[edge.edge_id] = edge.measured_pressure_pa;
		}

		int output_failed = 0;
		std::string output_error;
		if (rank == 0) try {
			WriteOutputs(options.output_directory, configuration, graph_root, assets,
				bifurcation_holder, one_d, three_d, accepted);
		} catch (const std::exception& error) {
			output_failed = 1;
			output_error = error.what();
		}
		MPI_Bcast(&output_failed, 1, MPI_INT, 0, PETSC_COMM_WORLD);
		if (output_failed) {
			if (rank == 0) throw std::runtime_error(output_error);
			throw std::runtime_error("bifurcation output failed on rank 0");
		}
		if (rank == 0) {
#ifdef IGA_BIFURCATION_ENTRY
			std::cout << "completed schema-v5 1D--3D bifurcation steps="
				<< accepted.size() << " branches=" << bifurcation_holder->branches.size()
				<< " output=" << options.output_directory << '\n';
#else
			std::cout << "completed schema-v5 multidomain flow steps="
				<< accepted.size() << " domains=" << configuration.graph.Domains().size()
				<< " output=" << options.output_directory << '\n';
#endif
		}
	} catch (const std::exception& error) {
		if (rank == 0) std::cerr << error.what() << '\n';
		status = 1;
	}
	PetscFinalize();
	return status;
}
