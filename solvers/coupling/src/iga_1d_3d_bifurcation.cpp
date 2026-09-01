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
	for (const char character : value) {
		if (character == '"' || character == '\\') output << '\\';
		output << character;
	}
	return output.str();
}

struct NativeOneD {
	std::string domain_id;
	fs::path case_directory;
	iga::OneDInletPolicy inlet_policy = iga::OneDInletPolicy::ConfiguredOpenLoop;
	double initial_native_inlet_flow_m3_s = 0.0;
	std::unique_ptr<iga::OneDFlowRuntime> runtime;
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
	if (runtime->FlowState().outlets.size() != 1
		|| runtime->FlowState().outlets.front().kind != iga::OneDOutletKind::Pressure)
		throw std::runtime_error(
			"bifurcation benchmark supports one pressure-closed terminal per 1D domain");
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
	if (found == native.runtime->Network().node_index.end()
		|| native.runtime->FlowState().outlets.front().node != found->second)
		throw std::runtime_error("schema-v5 1D terminal locator does not match its native case");
}

struct AcceptedStep {
	int step = 0;
	double time_s = 0.0;
	int iterations = 0;
	double three_d_mass_m3_s = 0.0;
	double external_outward_flow_m3_s = 0.0;
	iga::PressureFlowStepResult result;
};

void WriteOutputs(const fs::path& directory,
	const iga::MultidomainConfiguration& configuration,
	const iga::OneDThreeDBifurcationDefinition& definition,
	const std::map<std::string, NativeOneD>& one_d,
	const std::vector<AcceptedStep>& steps)
{
	if (!fs::create_directories(directory))
		throw std::runtime_error("bifurcation output directory must be newly created");
	std::ofstream summary(directory/"pressure_flow_steps.csv");
	std::ofstream edges(directory/"pressure_flow_edges.csv");
	std::ofstream iterations(directory/"pressure_flow_iterations.csv");
	std::ofstream ports(directory/"pressure_flow_ports.csv");
	std::ofstream initialization(directory/"one_d_initialization.csv");
	if (!summary || !edges || !iterations || !ports || !initialization)
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
	for (const auto& item : one_d)
		initialization << std::setprecision(17) << item.first << ','
			<< (item.second.inlet_policy == iga::OneDInletPolicy::ConfiguredOpenLoop
				? "configured_open_loop" : "coupled_root") << ','
			<< item.second.initial_native_inlet_flow_m3_s << '\n';
	for (const auto& step : steps) {
		summary << std::setprecision(17) << step.step << ',' << step.time_s << ','
			<< step.iterations << ',' << step.three_d_mass_m3_s << ','
			<< step.external_outward_flow_m3_s << '\n';
		for (const auto& edge : step.result.iterations.back().edges)
			edges << std::setprecision(17) << step.step << ',' << step.time_s << ','
				<< edge.edge_id << ',' << edge.applied_pressure_pa << ','
				<< edge.measured_pressure_pa << ',' << edge.pressure_residual_pa << ','
				<< edge.normalized_pressure_residual << ',' << edge.first_outward_flow_m3_s
				<< ',' << edge.second_outward_flow_m3_s << ',' << edge.flow_residual_m3_s
				<< ',' << edge.normalized_flow_residual << '\n';
		for (const auto& iteration : step.result.iterations)
			for (const auto& edge : iteration.edges)
				iterations << std::setprecision(17) << step.step << ',' << step.time_s << ','
					<< iteration.iteration << ',' << iteration.applied_relaxation << ','
					<< (iteration.converged ? 1 : 0) << ',' << edge.edge_id << ','
					<< edge.applied_pressure_pa << ',' << edge.measured_pressure_pa << ','
					<< edge.pressure_residual_pa << ',' << edge.normalized_pressure_residual
					<< ',' << edge.first_outward_flow_m3_s << ','
					<< edge.second_outward_flow_m3_s << ',' << edge.flow_residual_m3_s << ','
					<< edge.normalized_flow_residual << '\n';
		for (const auto& port : step.result.accepted_ports) {
			ports << std::setprecision(17) << step.step << ',' << step.time_s << ','
				<< port.first.domain_id << ',' << port.first.port_id << ','
				<< RequireValue(port.second.area_m2, "port area") << ','
				<< RequireValue(port.second.outward_flow_m3_s, "port flow") << ',';
			if (port.second.mean_pressure_pa) ports << *port.second.mean_pressure_pa;
			ports << '\n';
		}
	}
	summary.close();
	edges.close();
	iterations.close();
	ports.close();
	initialization.close();
	if (!summary || !edges || !iterations || !ports || !initialization)
		throw std::runtime_error("cannot finalize bifurcation long-form output");
	std::ofstream marker(directory/"graph_binding_manifest.json.tmp");
	if (!marker) throw std::runtime_error("cannot create bifurcation completion marker");
	marker << "{\n  \"schema_version\": 5,\n  \"benchmark\": \"one_d_three_d_bifurcation\",\n"
		<< "  \"start_domain\": \"" << JsonEscape(configuration.start_domain_id) << "\",\n"
		<< "  \"three_d_domain\": \"" << JsonEscape(definition.three_d_domain_id) << "\",\n"
		<< "  \"branch_count\": " << definition.branches.size() << ",\n"
		<< "  \"completed_steps\": " << steps.size() << "\n}\n";
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
		std::optional<iga::OneDThreeDBifurcationDefinition> definition_holder;
		std::map<std::string, iga::ResolvedGraphDomainAssets> assets;
		std::map<std::string, NativeOneD> one_d;
		std::string local_error;
		try {
			if (fs::exists(options.output_directory))
				throw std::runtime_error("bifurcation output directory must not already exist");
			const auto root = iga::CanonicalGraphCaseRoot(options.graph_case);
			configuration_holder.emplace(iga::ReadMultidomainConfiguration(
				iga::ResolveContainedCaseFile(root, "simulation_config.json",
					"schema-v5 manifest").string()));
			definition_holder.emplace(
				iga::ResolveOneDThreeDBifurcation(*configuration_holder));
			assets = iga::ResolveGraphDomainAssets(*configuration_holder, root);
			one_d.emplace(definition_holder->upstream_domain_id,
				BuildOneD(*configuration_holder, assets,
					definition_holder->upstream_domain_id));
			for (const auto& branch : definition_holder->branches)
				one_d.emplace(branch.downstream_domain_id,
					BuildOneD(*configuration_holder, assets,
						branch.downstream_domain_id));
		} catch (const std::exception& error) {
			local_error = error.what();
		}
		RequireCollectivePreflight(local_error);
		auto& configuration = *configuration_holder;
		const auto& definition = *definition_holder;
		const auto three_d_case = assets.at(definition.three_d_domain_id).case_directory;
		std::unique_ptr<iga::Database> database_holder;
		std::optional<iga::SimulationConfiguration> three_d_configuration_holder;
		std::optional<iga::LabeledHexMesh> mesh_holder;
		std::vector<std::array<double, 3>> boundary_velocity;
		std::set<std::int32_t> wall_trace_basis;
		std::optional<iga::SimulationConfiguration> initial_three_d_holder;
		std::vector<iga::OutletModelState> outlet_models;
		std::optional<iga::ResolvedBoundaryConditions> initial_boundaries_holder;
		std::set<int> graph_outlet_labels;
		local_error.clear();
		try {
			for (const auto& item : one_d)
				for (const auto& port : configuration.graph.Domain(item.first).ports)
					if (port.locator_kind == "runtime_port")
						ValidateOneDPort(item.second, port);
			database_holder = std::make_unique<iga::Database>(
				assets.at(definition.three_d_domain_id).database.string());
			(void)database_holder->LoadRequired(rank);
			(void)database_holder->LoadOwned(rank);
			three_d_configuration_holder.emplace(iga::ReadSimulationConfiguration(
				iga::ResolveContainedCaseFile(three_d_case, "simulation_config.json",
					"3D configuration").string()));
			auto& candidate = *three_d_configuration_holder;
			CanonicalizePeriodicTableAssets(candidate.temporal_functions,
				three_d_case, definition.three_d_domain_id);
			if (candidate.equation_systems.size() != 1
				|| candidate.equation_systems.front().kind != iga::EquationKind::NavierStokes
				|| candidate.equation_systems.front().unknowns.size() != 2
				|| candidate.coupling.mode != iga::SimulationScopeMode::FlowOnly
				|| candidate.physiology.enabled
				|| candidate.equation_systems.front().time_integration != "backward_euler")
				throw std::runtime_error(
					"bifurcation coupling requires one flow_only backward_euler 3D Navier-Stokes system");
			for (const auto& field : candidate.fields)
				if (field.kind == iga::FieldKind::Scalar)
					throw std::runtime_error("bifurcation coupling rejects 3D transport fields");
			iga::RequireThreeDVascularPorts(candidate.coupling, "bifurcation coupling");
			const auto& candidate_ports = candidate.coupling.three_d_ports;
			const auto& candidate_inlet = configuration.graph.Port(definition.three_d_inlet);
			if (iga::ParseThreeDFlowBoundaryLabel(candidate_inlet)
				!= candidate_ports.inlet_label)
				throw std::runtime_error("schema-v5 3D inlet label does not match the native case");
			for (const auto& branch : definition.branches)
				graph_outlet_labels.insert(iga::ParseThreeDFlowBoundaryLabel(
					configuration.graph.Port(branch.three_d_outlet)));
			const std::set<int> native_outlet_labels(candidate_ports.outlet_labels.begin(),
				candidate_ports.outlet_labels.end());
			if (graph_outlet_labels != native_outlet_labels
				|| iga::ParseThreeDFlowBoundaryLabel(configuration.graph.Port(
					definition.three_d_wall_observation)) != 0)
				throw std::runtime_error(
					"schema-v5 3D outlet/wall labels do not match the native case");
			const auto& candidate_flow = candidate.equation_systems.front();
			for (const auto& item : one_d)
				if (item.second.runtime->FlowSystem().density != candidate_flow.density
					|| item.second.runtime->FlowSystem().dynamic_viscosity
						!= candidate_flow.viscosity)
					throw std::runtime_error(
						"bifurcation coupling requires identical density and viscosity");
			if (configuration.time.dt_s != candidate.time.dt
				|| configuration.time.steps != candidate.time.steps)
				throw std::runtime_error("schema-v5 and 3D time grids must match");
			mesh_holder.emplace(iga::ReadLabeledHexMesh(
				iga::ResolveContainedCaseFile(three_d_case, "controlmesh.vtk",
					"3D control mesh").string(), database_holder->header().nodes,
				database_holder->header().elements));
			boundary_velocity = iga::ReadVelocity(
				iga::ResolveContainedCaseFile(three_d_case, "initial_velocityfield.txt",
					"3D initial velocity").string(), database_holder->header().nodes);
			wall_trace_basis = iga::WallTraceBasis(*database_holder, *mesh_holder, 0);
			initial_three_d_holder.emplace(iga::MaterializeBoundaryWaveforms(candidate,
				three_d_case.string(), 0.0));
			outlet_models = iga::InitializeOutletModels(candidate, candidate_flow);
			if (!outlet_models.empty())
				throw std::runtime_error(
					"bifurcation benchmark requires static-pressure 3D outlets");
			initial_boundaries_holder.emplace(iga::ResolveFlowBoundaries(
				*initial_three_d_holder,
				initial_three_d_holder->equation_systems.front(), mesh_holder->labels,
				boundary_velocity));
		} catch (const std::exception& error) {
			local_error = error.what();
		}
		RequireCollectivePreflight(local_error);
		auto& database = *database_holder;
		auto& three_d_configuration = *three_d_configuration_holder;
		auto& mesh = *mesh_holder;
		auto& initial_three_d = *initial_three_d_holder;
		const auto& flow = three_d_configuration.equation_systems.front();
		const auto& native_ports = three_d_configuration.coupling.three_d_ports;
		const auto& graph_inlet = configuration.graph.Port(definition.three_d_inlet);
		iga::TransientFlowRuntime three_d(database, PETSC_COMM_WORLD, true, true,
			{flow.density, flow.viscosity, three_d_configuration.time.dt},
			*initial_boundaries_holder, mesh.labels, boundary_velocity, wall_trace_basis,
			std::move(outlet_models));
		iga::RequireValidGeometry(three_d.Elements(), rank, PETSC_COMM_WORLD);

		std::set<int> required_labels = graph_outlet_labels;
		required_labels.insert(0);
		required_labels.insert(native_ports.inlet_label);
		std::map<int, long long> local_faces;
		int local_unknown_label = 0;
		for (const auto& element : three_d.OwnedElements())
			for (const int label : element.boundary_labels) {
				if (label < 0) continue;
				if (!required_labels.count(label)) local_unknown_label = 1;
				else ++local_faces[label];
			}
		int unknown_label = 0;
		MPI_Allreduce(&local_unknown_label, &unknown_label, 1, MPI_INT, MPI_MAX,
			PETSC_COMM_WORLD);
		if (unknown_label)
			throw std::runtime_error("bifurcation 3D database has an undeclared boundary label");
		for (const int label : required_labels) {
			long long local = local_faces[label];
			long long global = 0;
			MPI_Allreduce(&local, &global, 1, MPI_LONG_LONG, MPI_SUM, PETSC_COMM_WORLD);
			if (global == 0)
				throw std::runtime_error("bifurcation 3D database omits a declared boundary face");
		}

		const auto& source = one_d.at(definition.upstream_domain_id);
		const auto source_terminal = source.runtime->GetPortState(
			configuration.graph.Port(definition.upstream_terminal).locator);
		iga::PortBoundaryData initial_inlet;
		initial_inlet.time_s = 0.0;
		initial_inlet.outward_flow_m3_s = -RequireValue(
			source_terminal.outward_flow_m3_s, "source terminal flow");
		const double native_reference = three_d.ReferenceBoundaryFlow(native_ports.inlet_label);
		const double graph_reference = graph_inlet.orientation.ToOutward(native_reference);
		iga::ApplyThreeDReferenceProfileInput(initial_three_d,
			initial_three_d.equation_systems.front(), graph_inlet, initial_inlet,
			graph_reference);
		const auto& pressure_field = flow.unknowns.at(1);
		for (const auto& branch : definition.branches)
			ApplyInitialPressure(initial_three_d, pressure_field,
				iga::ParseThreeDFlowBoundaryLabel(configuration.graph.Port(branch.three_d_outlet)),
				configuration.initial_pressure_pa.at(branch.edge_id));
		three_d.InitializeState(initial_three_d);

		std::vector<std::unique_ptr<iga::CoupledDomainRuntime>> runtimes;
		for (auto& item : one_d)
			runtimes.push_back(std::make_unique<iga::OneDFlowDomainAdapter>(item.first,
				*item.second.runtime, configuration.graph.Domain(item.first).ports,
				item.second.inlet_policy));
		iga::ThreeDFlowDomainControls controls;
		controls.maximum_newton = options.maximum_newton;
		controls.nonlinear_relative_tolerance = kNonlinearRelativeTolerance;
		controls.nonlinear_absolute_tolerance = kNonlinearAbsoluteTolerance;
		controls.mass_relative_tolerance = kMassRelativeTolerance;
		runtimes.push_back(std::make_unique<iga::ThreeDBodyFittedFlowDomainAdapter>(
			definition.three_d_domain_id, three_d,
			configuration.graph.Domain(definition.three_d_domain_id).ports,
			three_d_configuration, three_d_case,
			std::map<std::string, double>{{graph_inlet.id, graph_reference}}, controls));
		iga::DomainRuntimeRegistry registry(configuration.graph, std::move(runtimes));
		iga::PressureFlowComponentExecutor executor(registry,
			configuration.start_domain_id, iga::PressureFlowControlsFor(configuration.execution));

		const int final_step = options.stop_after_step > 0
			? options.stop_after_step : configuration.time.steps;
		if (final_step > configuration.time.steps)
			throw std::runtime_error("--stop-after-step exceeds configured steps");
		const int injected_failure_step = FailureInjectionStep();
		auto pressure = configuration.initial_pressure_pa;
		std::vector<AcceptedStep> accepted;
		for (int step = 1; step <= final_step; ++step) {
			const double time = step*configuration.time.dt_s;
			auto result = executor.Advance({step-1, (step-1)*configuration.time.dt_s,
				configuration.time.dt_s}, pressure,
				[&](const iga::PressureFlowStepResult&) {
					if (step == injected_failure_step)
						throw std::runtime_error("injected bifurcation failure before commit");
				});
			double mass = RequireValue(result.accepted_ports.at(
				definition.three_d_inlet).outward_flow_m3_s, "3D inlet flow")
				+RequireValue(result.accepted_ports.at(
					definition.three_d_wall_observation).outward_flow_m3_s, "3D wall flow");
			for (const auto& branch : definition.branches)
				mass += RequireValue(result.accepted_ports.at(
					branch.three_d_outlet).outward_flow_m3_s, "3D outlet flow");
			double external = RequireValue(result.accepted_ports.at(
				definition.upstream_root_observation).outward_flow_m3_s, "source root flow");
			for (const auto& branch : definition.branches)
				external += RequireValue(result.accepted_ports.at(
					branch.downstream_terminal_observation).outward_flow_m3_s,
					"branch terminal flow");
			accepted.push_back({step, time, static_cast<int>(result.iterations.size()),
				mass, external, std::move(result)});
			for (const auto& edge : accepted.back().result.iterations.back().edges)
				pressure[edge.edge_id] = edge.measured_pressure_pa;
		}

		int output_failed = 0;
		std::string output_error;
		if (rank == 0) try {
			WriteOutputs(options.output_directory, configuration, definition, one_d, accepted);
		} catch (const std::exception& error) {
			output_failed = 1;
			output_error = error.what();
		}
		MPI_Bcast(&output_failed, 1, MPI_INT, 0, PETSC_COMM_WORLD);
		if (output_failed) {
			if (rank == 0) throw std::runtime_error(output_error);
			throw std::runtime_error("bifurcation output failed on rank 0");
		}
		if (rank == 0)
			std::cout << "completed schema-v5 1D--3D bifurcation steps="
				<< accepted.size() << " branches=" << definition.branches.size()
				<< " output=" << options.output_directory << '\n';
	} catch (const std::exception& error) {
		if (rank == 0) std::cerr << error.what() << '\n';
		status = 1;
	}
	PetscFinalize();
	return status;
}
