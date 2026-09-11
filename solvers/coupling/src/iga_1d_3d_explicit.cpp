#include "CheckedText.hpp"
#include "BoundarySupport.hpp"
#include "AitkenRelaxation.hpp"
#include "ExplicitOneDThreeDCoupling.hpp"
#include "IgaDatabase.hpp"
#include "DomainRuntimeRegistry.hpp"
#include "ExecutionResources.hpp"
#include "CollectiveDomainRuntimeRegistry.hpp"
#include "CollectiveAssetInput.hpp"
#include "CollectivePetscOptions.hpp"
#include "MultidomainRunner.hpp"
#include "OneDImplicit.hpp"
#include "OneDFlowDomainAdapter.hpp"
#include "OneDRuntime.hpp"
#include "PressureFlowComponentExecutor.hpp"
#include "CollectivePressureFlowExecution.hpp"
#include "SequentialMultidomainCase.hpp"
#include "SequentialExecution.hpp"
#include "SimulationGraph.hpp"
#include "StrongOneDThreeDCoupling.hpp"
#include "ThreeDBodyFittedFlowDomainAdapter.hpp"
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
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

namespace fs = std::filesystem;

namespace {

// Only local work belongs here. Values are retained until the common outcome
// is known, then moved without throwing; no distributed owner is constructed
// or destroyed inside this helper.
template<class Work>
auto SequentialLocalValue(MPI_Comm communicator, const char* stage, Work&& work)
{
	using Value = std::invoke_result_t<Work>;
	static_assert(std::is_nothrow_move_constructible<Value>::value,
		"prepared sequential values must move without throwing");
	std::optional<Value> result;
	iga::RuntimeConstructionStage(communicator, stage, [&] {
		result.emplace(std::forward<Work>(work)());
	});
	return std::move(*result);
}

constexpr int kThreeDMaximumNewtonIterations = 30;
constexpr double kThreeDNonlinearRelativeTolerance = 1.0e-5;
constexpr double kThreeDNonlinearAbsoluteTolerance = 1.0e-10;
constexpr double kThreeDMassRelativeTolerance = 1.0e-3;

struct Options {
	bool graph_case_mode = false;
	fs::path graph_case;
	fs::path database;
	fs::path three_d_case;
	fs::path upstream_case;
	fs::path downstream_case;
	int upstream_terminal_node = -1;
	fs::path output_directory;
	int stop_after_step = 0;
	int three_d_max_newton = kThreeDMaximumNewtonIterations;
	bool strong_fixed = false;
	bool strong_aitken = false;
	iga::StrongCouplingControls strong_controls;
	iga::AitkenRelaxationControls aitken_controls;
	bool strong_pressure_reference_set = false;
	bool coupling_mode_argument_seen = false;
	bool strong_argument_seen = false;
	bool aitken_argument_seen = false;
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

double FiniteDouble(const std::string& text, const std::string& option)
{
	std::size_t used = 0;
	double value = 0.0;
	try { value = std::stod(text, &used); }
	catch (const std::exception&) { throw std::runtime_error(option+" requires a finite number"); }
	if (used != text.size() || !std::isfinite(value)) throw std::runtime_error(option+" requires a finite number");
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
	const std::string usage =
		"usage: iga_1d_3d_explicit --graph-case ROOT --output-dir DIR "
		"[--stop-after-step N] [--three-d-max-newton N] [PETSc options]\n"
		"   or: iga_1d_3d_explicit DB THREE_D_CASE UPSTREAM_1D_CASE DOWNSTREAM_1D_CASE "
		"--upstream-terminal-node ID --output-dir DIR [--stop-after-step N] [--three-d-max-newton N] "
		"[--coupling-mode explicit|strong-fixed|strong-aitken --strong-max-iterations N "
		"--strong-pressure-relative-tol R --strong-pressure-reference-pa PA "
		"--strong-flow-relative-tol R --strong-relaxation W "
		"--strong-aitken-min-relaxation W --strong-aitken-max-relaxation W] [PETSc options]";
	Options options;
	int first_option = 0;
	if (argc >= 3 && std::string(argv[1]) == "--graph-case") {
		options.graph_case_mode = true;
		options.graph_case = argv[2];
		first_option = 3;
	} else {
		if (argc < 5) throw std::runtime_error(usage);
		options.database = argv[1];
		options.three_d_case = argv[2];
		options.upstream_case = argv[3];
		options.downstream_case = argv[4];
		first_option = 5;
	}
	for (int i = first_option; i < argc; ++i) {
		const std::string argument(argv[i]);
		if (argument == "--upstream-terminal-node" || argument == "--output-dir"
			|| argument == "--stop-after-step" || argument == "--three-d-max-newton" || argument == "--coupling-mode"
			|| argument == "--strong-max-iterations" || argument == "--strong-pressure-relative-tol"
			|| argument == "--strong-pressure-reference-pa" || argument == "--strong-flow-relative-tol"
			|| argument == "--strong-relaxation" || argument == "--strong-aitken-min-relaxation"
			|| argument == "--strong-aitken-max-relaxation") {
			if (++i >= argc) throw std::runtime_error(argument+" requires a value");
			const std::string value(argv[i]);
			if (argument == "--upstream-terminal-node") options.upstream_terminal_node = PositiveInteger(value, argument);
			else if (argument == "--output-dir") options.output_directory = value;
			else if (argument == "--stop-after-step") options.stop_after_step = PositiveInteger(value, argument);
			else if (argument == "--three-d-max-newton") options.three_d_max_newton = PositiveInteger(value, argument);
			if (argument == "--coupling-mode") {
				options.coupling_mode_argument_seen = true;
				if (value == "explicit") { options.strong_fixed = false; options.strong_aitken = false; }
				else if (value == "strong-fixed") { options.strong_fixed = true; options.strong_aitken = false; }
				else if (value == "strong-aitken") { options.strong_fixed = false; options.strong_aitken = true; }
				else throw std::runtime_error("--coupling-mode must be explicit, strong-fixed, or strong-aitken");
			} else if (argument == "--strong-max-iterations") {
				options.strong_argument_seen = true;
				options.strong_controls.maximum_iterations = PositiveInteger(value, argument);
			} else if (argument == "--strong-pressure-relative-tol") {
				options.strong_argument_seen = true;
				options.strong_controls.pressure_relative_tolerance = FiniteDouble(value, argument);
			} else if (argument == "--strong-pressure-reference-pa") {
				options.strong_argument_seen = true;
				options.strong_pressure_reference_set = true;
				options.strong_controls.pressure_reference_pa = FiniteDouble(value, argument);
			} else if (argument == "--strong-flow-relative-tol") {
				options.strong_argument_seen = true;
				options.strong_controls.flow_relative_tolerance = FiniteDouble(value, argument);
			} else if (argument == "--strong-relaxation") {
				options.strong_argument_seen = true;
				options.strong_controls.relaxation_factor = FiniteDouble(value, argument);
				options.aitken_controls.initial_relaxation = options.strong_controls.relaxation_factor;
			} else if (argument == "--strong-aitken-min-relaxation") {
				options.strong_argument_seen = true;
				options.aitken_argument_seen = true;
				options.aitken_controls.minimum_relaxation = FiniteDouble(value, argument);
			} else if (argument == "--strong-aitken-max-relaxation") {
				options.strong_argument_seen = true;
				options.aitken_argument_seen = true;
				options.aitken_controls.maximum_relaxation = FiniteDouble(value, argument);
			}
			continue;
		}
		if (iga::IsExplicitCouplingPetscOption(argument)) {
			if (i+1 < argc && iga::ExplicitCouplingPetscOptionConsumesNextValue(argument, argv[i+1])) ++i;
			continue;
		}
		throw std::runtime_error("unexpected argument: "+argument);
	}
	if (options.output_directory.empty())
		throw std::runtime_error("--output-dir is required");
	if (!options.graph_case_mode && options.upstream_terminal_node < 0)
		throw std::runtime_error("--upstream-terminal-node is required for legacy positional input");
	if (options.graph_case_mode && options.upstream_terminal_node >= 0)
		throw std::runtime_error("--graph-case derives the upstream terminal from schema v5");
	if (options.graph_case_mode && (options.coupling_mode_argument_seen
		|| options.strong_argument_seen))
		throw std::runtime_error("--graph-case takes coupling execution controls from schema v5");
	if (!options.strong_fixed && !options.strong_aitken && options.strong_argument_seen)
		throw std::runtime_error("strong coupling arguments require a strong coupling mode");
	if (!options.strong_aitken && options.aitken_argument_seen)
		throw std::runtime_error("Aitken bounds require --coupling-mode strong-aitken");
	if (options.strong_fixed || options.strong_aitken) {
		if (!options.strong_pressure_reference_set)
			throw std::runtime_error("strong coupling requires --strong-pressure-reference-pa");
		iga::ValidateStrongCouplingControls(options.strong_controls);
		if (options.strong_aitken) iga::ValidateAitkenRelaxationControls(options.aitken_controls);
	}
	return options;
}

void ApplyGraphExecution(Options& options, const iga::GraphExecutionDefinition& execution)
{
	options.strong_fixed = execution.kind == iga::GraphExecutionKind::Fixed;
	options.strong_aitken = execution.kind == iga::GraphExecutionKind::Aitken;
	if (!options.strong_fixed && !options.strong_aitken) return;
	options.strong_controls.maximum_iterations = execution.maximum_iterations;
	options.strong_controls.pressure_relative_tolerance
		= execution.pressure_relative_tolerance;
	options.strong_controls.pressure_reference_pa = execution.pressure_reference_pa;
	options.strong_controls.flow_relative_tolerance = execution.flow_relative_tolerance;
	options.strong_controls.relaxation_factor = execution.relaxation_factor;
	options.aitken_controls.initial_relaxation = execution.relaxation_factor;
	options.aitken_controls.minimum_relaxation = execution.minimum_relaxation;
	options.aitken_controls.maximum_relaxation = execution.maximum_relaxation;
	iga::ValidateStrongCouplingControls(options.strong_controls);
	if (options.strong_aitken)
		iga::ValidateAitkenRelaxationControls(options.aitken_controls);
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

iga::CouplingPort MakeThreeDPort(const std::string& subsystem_id, const std::string& id, int label,
	std::set<iga::PortQuantity> provides, std::set<iga::PortQuantity> requires)
{
	iga::CouplingPort port;
	port.id = id;
	port.subsystem_id = subsystem_id;
	port.locator_kind = "boundary_label";
	port.locator = std::to_string(label);
	port.provides = std::move(provides);
	port.requires = std::move(requires);
	iga::ValidateCouplingPort(port);
	return port;
}

iga::CouplingPort MakeOneDLogicalPort(const std::string& subsystem_id, const std::string& id,
	const std::string& runtime_port_id, iga::PortQuantity accepted)
{
	iga::CouplingPort port;
	port.id = id;
	port.subsystem_id = subsystem_id;
	port.locator_kind = "runtime_port";
	port.locator = runtime_port_id;
	port.provides = {iga::PortQuantity::Area, iga::PortQuantity::FlowRate,
		iga::PortQuantity::MeanPressure};
	port.requires = {accepted};
	iga::ValidateCouplingPort(port);
	return port;
}

iga::CouplingPort MakeOneDObservationPort(const std::string& subsystem_id,
	const std::string& id, const std::string& runtime_port_id)
{
	iga::CouplingPort port;
	port.id = id;
	port.subsystem_id = subsystem_id;
	port.locator_kind = "runtime_port";
	port.locator = runtime_port_id;
	port.provides = {iga::PortQuantity::Area, iga::PortQuantity::FlowRate,
		iga::PortQuantity::MeanPressure};
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
	escaped.exceptions(std::ios::badbit | std::ios::failbit);
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

const char* GraphExecutionKindName(iga::GraphExecutionKind kind)
{
	if (kind == iga::GraphExecutionKind::Explicit) return "explicit";
	if (kind == iga::GraphExecutionKind::Fixed) return "fixed";
	if (kind == iga::GraphExecutionKind::Aitken) return "aitken";
	return "unknown";
}

void WriteGraphBindingManifest(const fs::path& path,
	const iga::MultidomainConfiguration& configuration,
	const std::map<std::string, iga::ResolvedGraphDomainAssets>& assets,
	const fs::path& graph_case)
{
	auto temporary = path;
	temporary += ".tmp";
	std::ofstream output(temporary);
	if (!output) throw std::runtime_error("cannot create graph binding manifest");
	output << std::setprecision(17) << "{\n"
		<< "  \"schema_version\": 5,\n"
		<< "  \"graph_case\": \"" << JsonEscape(fs::canonical(graph_case).string())
		<< "\",\n"
		<< "  \"start_domain\": \"" << JsonEscape(configuration.start_domain_id)
		<< "\",\n"
		<< "  \"execution\": \"" << GraphExecutionKindName(configuration.execution.kind)
		<< "\",\n"
		<< "  \"domains\": [\n";
	for (std::size_t i = 0; i < configuration.domains.size(); ++i) {
		const auto& domain = configuration.domains[i];
		const auto& resolved = assets.at(domain.id);
		output << "    {\"id\":\"" << JsonEscape(domain.id) << "\",\"kind\":\""
			<< iga::DomainKindName(domain.kind) << "\",\"case\":\""
			<< JsonEscape(resolved.case_directory.string()) << "\"";
		if (!resolved.database.empty()) output << ",\"database\":\""
			<< JsonEscape(resolved.database.string()) << "\"";
		output << "}" << (i+1 == configuration.domains.size() ? "\n" : ",\n");
	}
	output << "  ],\n  \"couplings\": [\n";
	for (std::size_t i = 0; i < configuration.graph.Edges().size(); ++i) {
		const auto& edge = configuration.graph.Edges()[i];
		output << "    {\"id\":\"" << JsonEscape(edge.id)
			<< "\",\"a\":{\"domain\":\"" << JsonEscape(edge.first.domain_id)
			<< "\",\"port\":\"" << JsonEscape(edge.first.port_id)
			<< "\"},\"b\":{\"domain\":\"" << JsonEscape(edge.second.domain_id)
			<< "\",\"port\":\"" << JsonEscape(edge.second.port_id)
			<< "\"},\"initial_pressure_pa\":"
			<< configuration.initial_pressure_pa.at(edge.id) << "}"
			<< (i+1 == configuration.graph.Edges().size() ? "\n" : ",\n");
	}
	output << "  ]\n}\n";
	if (!output) throw std::runtime_error("cannot write graph binding manifest");
	output.close();
	if (!output) throw std::runtime_error("cannot close graph binding manifest");
	std::filesystem::rename(temporary, path);
}

struct OneDWorkTotals {
	long long accepted_configured = 0, all_configured = 0, rejected_configured = 0;
	long long accepted_cfl = 0, all_cfl = 0, rejected_cfl = 0;
};

OneDWorkTotals SumOneDWork(const std::vector<iga::ExplicitCouplingHistoryRow>& history, bool upstream)
{
	OneDWorkTotals totals;
	for (const auto& row : history) {
		const long long accepted_configured = upstream ? row.upstream_accepted_configured_substeps : row.downstream_accepted_configured_substeps;
		const long long all_configured = upstream ? row.upstream_all_configured_substeps : row.downstream_all_configured_substeps;
		const long long rejected_configured = upstream ? row.upstream_rejected_configured_substeps : row.downstream_rejected_configured_substeps;
		const long long accepted_cfl = upstream ? row.upstream_accepted_explicit_cfl_substeps : row.downstream_accepted_explicit_cfl_substeps;
		const long long all_cfl = upstream ? row.upstream_all_explicit_cfl_substeps : row.downstream_all_explicit_cfl_substeps;
		const long long rejected_cfl = upstream ? row.upstream_rejected_explicit_cfl_substeps : row.downstream_rejected_explicit_cfl_substeps;
		totals.accepted_configured += accepted_configured; totals.all_configured += all_configured; totals.rejected_configured += rejected_configured;
		totals.accepted_cfl += accepted_cfl; totals.all_cfl += all_cfl; totals.rejected_cfl += rejected_cfl;
	}
	return totals;
}

void WriteExplicitCouplingManifest(const fs::path& path, int upstream_terminal_node,
	int three_d_inlet_label, int three_d_outlet_label, double dt_s, int configured_steps,
	int completed_steps,
	double density_kg_m3, double dynamic_viscosity_pa_s, double normalized_length_m,
	double reference_inlet_outward_flow_m3_s, double initial_three_d_inlet_pressure_pa,
	double initial_downstream_root_pressure_pa, const iga::OneDSubcyclingPlan& upstream_subcycling,
	const iga::OneDSubcyclingPlan& downstream_subcycling, const std::vector<iga::ExplicitCouplingHistoryRow>& history,
	int three_d_max_newton)
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
		<< "  \"newton_controls\": {\"maximum_iterations\": " << three_d_max_newton
		<< ", \"nonlinear_relative_tolerance\": " << kThreeDNonlinearRelativeTolerance
		<< ", \"nonlinear_absolute_tolerance\": " << kThreeDNonlinearAbsoluteTolerance
		<< ", \"mass_relative_tolerance\": " << kThreeDMassRelativeTolerance << "},\n"
		<< "  \"initial_lagged_pressures_pa\": {\"three_d_inlet\": "
		<< initial_three_d_inlet_pressure_pa << ", \"downstream_root\": "
		<< initial_downstream_root_pressure_pa << "},\n"
		<< "  \"one_d_subcycling\": {\"macro_dt_s\": " << upstream_subcycling.macro_dt_s
		<< ", \"macro_steps\": " << upstream_subcycling.macro_steps
		<< ", \"ratio_relative_tolerance\": 1e-12, \"upstream\": {\"configured_dt_s\": " << upstream_subcycling.configured_dt_s
		<< ", \"configured_steps\": " << upstream_subcycling.configured_steps << ", \"N\": " << upstream_subcycling.configured_substeps_per_macro_step
		<< "}, \"downstream\": {\"configured_dt_s\": " << downstream_subcycling.configured_dt_s
		<< ", \"configured_steps\": " << downstream_subcycling.configured_steps << ", \"N\": " << downstream_subcycling.configured_substeps_per_macro_step
		<< "}, \"semantics\": \"zero-order held interface data; configured open-loop sampling at substep endpoints; endpoint residuals; internal_substeps is explicit CFL only\"},\n";
	const auto upstream_work = SumOneDWork(history, true);
	const auto downstream_work = SumOneDWork(history, false);
	output << "  \"one_d_work_totals\": {\"upstream\": {\"accepted_configured_substeps\": " << upstream_work.accepted_configured
		<< ", \"all_configured_substeps\": " << upstream_work.all_configured << ", \"rejected_configured_substeps\": " << upstream_work.rejected_configured
		<< ", \"accepted_explicit_cfl_substeps\": " << upstream_work.accepted_cfl << ", \"all_explicit_cfl_substeps\": " << upstream_work.all_cfl
		<< ", \"rejected_explicit_cfl_substeps\": " << upstream_work.rejected_cfl << "}, \"downstream\": {\"accepted_configured_substeps\": " << downstream_work.accepted_configured
		<< ", \"all_configured_substeps\": " << downstream_work.all_configured << ", \"rejected_configured_substeps\": " << downstream_work.rejected_configured
		<< ", \"accepted_explicit_cfl_substeps\": " << downstream_work.accepted_cfl << ", \"all_explicit_cfl_substeps\": " << downstream_work.all_cfl
		<< ", \"rejected_explicit_cfl_substeps\": " << downstream_work.rejected_cfl << "}}\n"
		<< "}\n";
	output.close();
	if (!output) throw std::runtime_error("cannot write explicit coupling manifest");
}

void WriteStrongCouplingManifest(const fs::path& path, const iga::StrongCouplingControls& controls,
	const iga::AitkenRelaxationControls& aitken_controls, bool strong_aitken,
	int configured_steps, int completed_steps, double dt_s, double density_kg_m3,
	double dynamic_viscosity_pa_s, double normalized_length_m, int upstream_terminal_node,
	int three_d_inlet_label, int three_d_outlet_label, double initial_upstream_pressure_pa,
	double initial_three_d_outlet_traction_pressure_pa, long long total_coupling_iterations,
	double reference_inlet_outward_flow_m3_s, const std::array<long long, 6>& aitken_status_counts,
	long long accepted_ksp, long long all_ksp, const iga::OneDSubcyclingPlan& upstream_subcycling,
	const iga::OneDSubcyclingPlan& downstream_subcycling, const std::vector<iga::ExplicitCouplingHistoryRow>& history,
	int three_d_max_newton)
{
	std::ofstream output(path);
	if (!output) throw std::runtime_error("cannot create strong coupling manifest");
	output << std::setprecision(17) << "{\n"
		<< "  \"scheme\": \"" << (strong_aitken ? "strong_aitken" : "strong_fixed") << "\",\n"
		<< "  \"pressure_fixed_point\": \"x=[upstream_terminal_pressure,three_d_outlet_pressure_traction_parameter]; G=[three_d_inlet_mean_static_pressure,downstream_root_mean_static_pressure]\",\n"
		<< "  \"pressure_traction_parameter\": \"The 3D outlet input is a pressure-traction parameter; measured 3D outlet static pressure is diagnostic only and is not the right fixed-point component.\",\n"
		<< "  \"restart\": \"unsupported\",\n"
		<< "  \"ports\": [{\"id\": \"upstream_terminal\", \"locator\": \"outlet:"
		<< upstream_terminal_node << "\", \"native_to_outward_sign\": 1}, "
		<< "{\"id\": \"three_d_inlet\", \"locator_kind\": \"boundary_label\", \"locator\": "
		<< three_d_inlet_label << ", \"native_to_outward_sign\": 1}, "
		<< "{\"id\": \"three_d_outlet\", \"locator_kind\": \"boundary_label\", \"locator\": "
		<< three_d_outlet_label << ", \"native_to_outward_sign\": 1}, "
		<< "{\"id\": \"downstream_root\", \"locator\": \"root\", \"native_to_outward_sign\": -1}],\n"
		<< "  \"dt_s\": " << dt_s << ",\n"
		<< "  \"configured_steps\": " << configured_steps << ",\n"
		<< "  \"completed_steps\": " << completed_steps << ",\n"
		<< "  \"density_kg_m3\": " << density_kg_m3 << ",\n"
		<< "  \"dynamic_viscosity_pa_s\": " << dynamic_viscosity_pa_s << ",\n"
		<< "  \"geometry_transform_product_m\": " << normalized_length_m << ",\n"
		<< "  \"reference_inlet_outward_flow_m3_s\": " << reference_inlet_outward_flow_m3_s << ",\n"
		<< "  \"initial_pressure_guesses_pa\": {\"upstream_terminal\": "
		<< initial_upstream_pressure_pa << ", \"three_d_outlet_traction\": "
		<< initial_three_d_outlet_traction_pressure_pa << "},\n"
		<< "  \"newton_controls\": {\"maximum_iterations\": " << three_d_max_newton
		<< ", \"nonlinear_relative_tolerance\": " << kThreeDNonlinearRelativeTolerance
		<< ", \"nonlinear_absolute_tolerance\": " << kThreeDNonlinearAbsoluteTolerance
		<< ", \"mass_relative_tolerance\": " << kThreeDMassRelativeTolerance << "},\n"
		<< "  \"controls\": {\"maximum_iterations\": " << controls.maximum_iterations
		<< ", \"pressure_relative_tolerance\": " << controls.pressure_relative_tolerance
		<< ", \"pressure_reference_pa\": " << controls.pressure_reference_pa
		<< ", \"flow_relative_tolerance\": " << controls.flow_relative_tolerance
		<< ", \"relaxation_factor\": " << controls.relaxation_factor << "},\n";
	if (strong_aitken) output
		<< "  \"aitken_controls\": {\"initial_relaxation\": " << aitken_controls.initial_relaxation
		<< ", \"minimum_relaxation\": " << aitken_controls.minimum_relaxation
		<< ", \"maximum_relaxation\": " << aitken_controls.maximum_relaxation
		<< ", \"scaled_difference_threshold\": " << aitken_controls.scaled_difference_threshold << "},\n"
		<< "  \"aitken_status_mapping\": {\"0\": \"Initial\", \"1\": \"Dynamic\", \"2\": \"ClampedMinimum\", \"3\": \"ClampedMaximum\", \"4\": \"TinyDifferenceFallback\", \"5\": \"NonfiniteCandidateFallback\"},\n"
		<< "  \"aitken_status_totals\": {\"initial\": " << aitken_status_counts[0]
		<< ", \"dynamic\": " << aitken_status_counts[1] << ", \"clamped_minimum\": " << aitken_status_counts[2]
		<< ", \"clamped_maximum\": " << aitken_status_counts[3] << ", \"tiny_difference_fallback\": " << aitken_status_counts[4]
		<< ", \"nonfinite_candidate_fallback\": " << aitken_status_counts[5] << "},\n";
	if (strong_aitken) output
		<< "  \"aitken_recurrence\": \"r=G-x ordered [upstream_terminal/three_d_inlet,three_d_outlet_pressure_traction_parameter/downstream_root]; scale=max(Pref,abs(all previous/current r components)); delta=(r_current/scale)-(r_previous/scale); omega_hat=-omega_previous*dot(r_previous/scale,delta)/dot(delta,delta); denominator <= threshold or nonfinite candidate retains prior clamped omega; otherwise exact min/max clamp; reset each macro-step; AcceptApplied occurs only after rollback; final converged proposal is unused\",\n";
	output
		<< "  \"formula\": \"pressure_raw_i=G_i-x_i; eta_i=abs(pressure_raw_i)/max(Pref,abs(G_i),abs(x_i)); convergence requires both nonnegative normalized pressure residuals eta_i <= pressure_relative_tolerance and both absolute signed outward-flow residuals abs((Q_a_out+Q_b_out)/max(abs(Q_a_out),abs(Q_b_out),1e-30)) <= flow_relative_tolerance; x_next=x+omega*(G-x)\",\n"
		<< "  \"output_semantics\": \"History, iteration CSV, and manifest are written only after the complete run succeeds; rejected trials produce no persistent coupling output.\",\n"
		<< "  \"work_semantics\": \"accepted 3D KSP work is the final committed attempt of each step; all_attempts includes rejected strong sweeps; rejected is all_attempts minus accepted. Per-domain configured and explicit-CFL 1D accepted/all/rejected work is serialized in every step row.\",\n"
		<< "  \"total_coupling_iterations\": " << total_coupling_iterations << ",\n"
		<< "  \"three_d_ksp_iterations\": {\"accepted\": " << accepted_ksp
		<< ", \"all_attempts\": " << all_ksp << ", \"rejected\": "
		<< all_ksp-accepted_ksp << "},\n";
	const auto upstream_work = SumOneDWork(history, true);
	const auto downstream_work = SumOneDWork(history, false);
	output << "  \"one_d_subcycling\": {\"macro_dt_s\": " << upstream_subcycling.macro_dt_s
		<< ", \"macro_steps\": " << upstream_subcycling.macro_steps << ", \"ratio_relative_tolerance\": 1e-12"
		<< ", \"upstream\": {\"configured_dt_s\": " << upstream_subcycling.configured_dt_s << ", \"configured_steps\": " << upstream_subcycling.configured_steps << ", \"N\": " << upstream_subcycling.configured_substeps_per_macro_step
		<< "}, \"downstream\": {\"configured_dt_s\": " << downstream_subcycling.configured_dt_s << ", \"configured_steps\": " << downstream_subcycling.configured_steps << ", \"N\": " << downstream_subcycling.configured_substeps_per_macro_step
		<< "}, \"semantics\": \"zero-order held interface data; configured open-loop sampling at substep endpoints; endpoint residuals; internal_substeps is explicit CFL only\"},\n"
		<< "  \"one_d_work_totals\": {\"upstream\": {\"accepted_configured_substeps\": " << upstream_work.accepted_configured << ", \"all_configured_substeps\": " << upstream_work.all_configured << ", \"rejected_configured_substeps\": " << upstream_work.rejected_configured << ", \"accepted_explicit_cfl_substeps\": " << upstream_work.accepted_cfl << ", \"all_explicit_cfl_substeps\": " << upstream_work.all_cfl << ", \"rejected_explicit_cfl_substeps\": " << upstream_work.rejected_cfl
		<< "}, \"downstream\": {\"accepted_configured_substeps\": " << downstream_work.accepted_configured << ", \"all_configured_substeps\": " << downstream_work.all_configured << ", \"rejected_configured_substeps\": " << downstream_work.rejected_configured << ", \"accepted_explicit_cfl_substeps\": " << downstream_work.accepted_cfl << ", \"all_explicit_cfl_substeps\": " << downstream_work.all_cfl << ", \"rejected_explicit_cfl_substeps\": " << downstream_work.rejected_cfl << "}}\n}\n";
	output.close();
	if (!output) throw std::runtime_error("cannot write strong coupling manifest");
}

std::string CsvLineWithoutNewline(const std::string& line)
{
	if (line.empty() || line.back() != '\n') throw std::runtime_error("invalid CSV serializer output");
	return line.substr(0, line.size()-1);
}

void WriteStrongFailureIteration(std::ostream& output, const iga::StrongCouplingIterationRow& row)
{
	output << std::setprecision(17) << "step=" << row.physical_step << " iteration=" << row.iteration
		<< " x=(" << row.applied_upstream_terminal_pressure_pa << ',' << row.applied_three_d_outlet_traction_pressure_pa << ')'
		<< " G=(" << row.measured_three_d_inlet_pressure_pa << ',' << row.measured_downstream_root_pressure_pa << ')'
		<< " raw_p=(" << row.signed_upstream_pressure_residual_pa << ',' << row.signed_downstream_pressure_residual_pa << ')'
		<< " normalized_p=(" << row.normalized_upstream_pressure_residual << ',' << row.normalized_downstream_pressure_residual << ')'
		<< " q=(" << row.upstream_terminal_outward_flow_m3_s << ',' << row.three_d_inlet_outward_flow_m3_s
		<< ',' << row.three_d_outlet_outward_flow_m3_s << ',' << row.downstream_root_outward_flow_m3_s << ')'
		<< " raw_q=(" << row.upstream_three_d_flow_residual_m3_s << ',' << row.three_d_downstream_flow_residual_m3_s << ')'
		<< " normalized_q=(" << row.normalized_upstream_three_d_flow_residual << ',' << row.normalized_three_d_downstream_flow_residual << ')'
		<< " wall=" << row.three_d_wall_outward_flow_m3_s << " mass=" << row.three_d_mass_imbalance_m3_s
		<< " omega=" << row.relaxation_factor_for_next_guess << " unclamped=" << row.unclamped_relaxation_factor
		<< " numerator=" << row.aitken_scaled_numerator << " denominator=" << row.aitken_scaled_denominator
		<< " status=" << row.aitken_status_code << " has_previous=" << row.aitken_has_previous_residual
		<< " applied=" << row.relaxation_update_applied << " converged=" << row.converged
		<< " ksp=(" << row.three_d_attempt_linear_iterations << ',' << row.three_d_cumulative_step_linear_iterations << ")\n";
	output << "upstream_1d_work=(" << row.upstream_1d_attempt_configured_substeps << ','
		<< row.upstream_1d_attempt_explicit_cfl_substeps << ',' << row.upstream_1d_cumulative_configured_substeps << ','
		<< row.upstream_1d_cumulative_explicit_cfl_substeps << ") downstream_1d_work=("
		<< row.downstream_1d_attempt_configured_substeps << ',' << row.downstream_1d_attempt_explicit_cfl_substeps << ','
		<< row.downstream_1d_cumulative_configured_substeps << ',' << row.downstream_1d_cumulative_explicit_cfl_substeps << ")\n";
}

void WriteStrongOneDWorkTuple(std::ostream& output, const char* name, const iga::OneDFlowRuntime::TrialDiagnostics& diagnostics,
	const iga::OneDTrialWorkAccumulator& cumulative)
{
	output << name << "_1d_work=(" << diagnostics.attempted_configured_substeps << ','
		<< diagnostics.explicit_cfl_substep_delta << ',' << cumulative.configured_substeps << ','
		<< cumulative.explicit_cfl_substeps << ')';
}

void RollbackSolved(iga::OneDFlowRuntime& runtime)
{
	if (runtime.CurrentPhase() == iga::OneDFlowRuntime::Phase::TrialSolved) runtime.RollbackTrial();
}

void RollbackSolved(iga::TransientFlowRuntime& runtime)
{
	if (runtime.Phase() == iga::FlowStepPhase::TrialSolved) runtime.RollbackTrial();
}

void SetOneDAttemptWork(iga::ExplicitCouplingHistoryRow& row, const iga::OneDFlowRuntime& upstream,
	const iga::OneDFlowRuntime& downstream)
{
	row.upstream_configured_substeps_attempted = upstream.Diagnostics().attempted_configured_substeps;
	row.upstream_explicit_cfl_substeps = upstream.Diagnostics().explicit_cfl_substep_delta;
	row.downstream_configured_substeps_attempted = downstream.Diagnostics().attempted_configured_substeps;
	row.downstream_explicit_cfl_substeps = downstream.Diagnostics().explicit_cfl_substep_delta;
	row.upstream_accepted_configured_substeps = row.upstream_all_configured_substeps = row.upstream_configured_substeps_attempted;
	row.upstream_accepted_explicit_cfl_substeps = row.upstream_all_explicit_cfl_substeps = row.upstream_explicit_cfl_substeps;
	row.downstream_accepted_configured_substeps = row.downstream_all_configured_substeps = row.downstream_configured_substeps_attempted;
	row.downstream_accepted_explicit_cfl_substeps = row.downstream_all_explicit_cfl_substeps = row.downstream_explicit_cfl_substeps;
}

} // namespace

int iga::RunSequentialFlow(int argc, char** argv, MPI_Comm communicator)
{
	int rank = 0;
	MPI_Comm_rank(communicator, &rank);
	int status = 0;
	try {
		iga::RequireExecutionResources(communicator, &std::cout);
		auto options = SequentialLocalValue(communicator, "sequential arguments", [&] {
			return ParseOptions(argc, argv);
		});
		iga::RequireCollectiveSameInt(communicator, "sequential input mode", options.graph_case_mode);
		std::optional<iga::MultidomainConfiguration> graph_configuration;
		std::optional<iga::SequentialOneDThreeDOneDDefinition> graph_definition;
		std::map<std::string, iga::ResolvedGraphDomainAssets> graph_assets;
		fs::path graph_case_root;
		if (options.graph_case_mode) {
			iga::AssetFileCatalog graph_files;
			iga::RuntimeConstructionStage(communicator, "sequential graph paths", [&] {
				graph_case_root = iga::CanonicalGraphCaseRoot(options.graph_case);
				graph_files.emplace("sequential graph manifest", iga::ResolveContainedCaseFile(
					graph_case_root, "simulation_config.json", "schema-v5 manifest"));
			});
			iga::RequireCollectiveAssetFiles(communicator, graph_files);
			iga::RuntimeConstructionStage(communicator, "sequential graph definition", [&] {
				if (fs::exists(options.output_directory))
					throw std::runtime_error(
						"schema-v5 graph output directory must not already exist");
				const auto& manifest = graph_files.at("sequential graph manifest");
				graph_configuration.emplace(
					iga::ReadMultidomainConfiguration(manifest.string()));
				if (graph_configuration->schema_version != 5)
					throw std::runtime_error(
						"sequential flow runner supports schema version 5; schema 6 species is not yet executable");
				graph_definition.emplace(
					iga::ResolveSequentialOneDThreeDOneD(*graph_configuration));
				graph_assets = iga::ResolveGraphDomainAssets(*graph_configuration,
					graph_case_root);
				options.upstream_case
					= graph_assets.at(graph_definition->upstream_domain_id).case_directory;
				options.three_d_case
					= graph_assets.at(graph_definition->three_d_domain_id).case_directory;
				options.database
					= graph_assets.at(graph_definition->three_d_domain_id).database;
				options.downstream_case
					= graph_assets.at(graph_definition->downstream_domain_id).case_directory;
				options.upstream_terminal_node = iga::ParseOneDOutletNodeLocator(
					graph_configuration->graph.Port(graph_definition->upstream_terminal));
				ApplyGraphExecution(options, graph_configuration->execution);
			});
		}
		const int injected_failure_step = SequentialLocalValue(communicator, "sequential failure control", [&] {
			return ExplicitCouplingFailureInjectionStep();
		});
		const auto execution_controls = SequentialLocalValue(communicator, "sequential controls", [&] {
			std::ostringstream text;
			text.exceptions(std::ios::badbit | std::ios::failbit);
			text << std::hexfloat << options.upstream_terminal_node << ' ' << options.stop_after_step
				<< ' ' << options.three_d_max_newton << ' ' << options.strong_fixed << ' ' << options.strong_aitken
				<< ' ' << options.strong_controls.maximum_iterations
				<< ' ' << options.strong_controls.pressure_relative_tolerance
				<< ' ' << options.strong_controls.pressure_reference_pa
				<< ' ' << options.strong_controls.flow_relative_tolerance
				<< ' ' << options.strong_controls.relaxation_factor
				<< ' ' << options.aitken_controls.initial_relaxation
				<< ' ' << options.aitken_controls.minimum_relaxation
				<< ' ' << options.aitken_controls.maximum_relaxation << ' ' << injected_failure_step;
			return text.str();
		});
		iga::RequireCollectiveSameText(communicator, "sequential control agreement", execution_controls);
		const auto application_options = SequentialLocalValue(communicator, "sequential option names", [&] {
			return std::set<std::string>{"--graph-case", "--output-dir", "--upstream-terminal-node",
				"--stop-after-step", "--three-d-max-newton", "--coupling-mode", "--strong-max-iterations",
				"--strong-pressure-relative-tol", "--strong-pressure-reference-pa", "--strong-flow-relative-tol",
				"--strong-relaxation", "--strong-aitken-min-relaxation", "--strong-aitken-max-relaxation", "-options_file"};
		});
		iga::RequireCollectivePetscOptions(communicator, nullptr, application_options);
		const auto case_file = [&](const fs::path& root, const fs::path& name, const char* context) {
			return options.graph_case_mode ? iga::ResolveContainedCaseFile(root, name, context) : root/name;
		};
		iga::AssetFileCatalog case_files;
		iga::RuntimeConstructionStage(communicator, "sequential case catalog", [&] {
			case_files = {{"sequential upstream configuration", case_file(options.upstream_case, "simulation_config.json", "upstream configuration")},
				{"sequential downstream configuration", case_file(options.downstream_case, "simulation_config.json", "downstream configuration")},
				{"sequential 3D configuration", case_file(options.three_d_case, "simulation_config.json", "3D configuration")},
				{"sequential database", options.database},
				{"sequential mesh", case_file(options.three_d_case, "controlmesh.vtk", "3D control mesh")},
				{"sequential velocity", case_file(options.three_d_case, "initial_velocityfield.txt", "3D initial velocity")}};
		});
		iga::RequireCollectiveAssetFiles(communicator, case_files);
		auto [upstream_configuration, downstream_configuration, three_d_configuration] =
			SequentialLocalValue(communicator, "sequential configurations", [&] {
				return std::make_tuple(
					iga::ParseOneDConfiguration(ReadText(case_files.at("sequential upstream configuration"))),
					iga::ParseOneDConfiguration(ReadText(case_files.at("sequential downstream configuration"))),
					iga::ReadSimulationConfiguration(case_files.at("sequential 3D configuration").string()));
			});
		iga::AssetFileCatalog selected_files;
		iga::RuntimeConstructionStage(communicator, "sequential selected assets", [&] {
			selected_files.emplace("sequential upstream geometry", case_file(options.upstream_case,
				upstream_configuration.geometry.file, "upstream geometry"));
			selected_files.emplace("sequential downstream geometry", case_file(options.downstream_case,
				downstream_configuration.geometry.file, "downstream geometry"));
			const auto inlet = iga::ResolveOneDInlet(upstream_configuration);
			if (!inlet.waveform.empty()) {
				const auto& function = iga::FindOneDTemporalFunction(upstream_configuration, inlet.waveform);
				if (function.kind == iga::TemporalFunctionKind::PeriodicTable)
					selected_files.emplace("sequential upstream temporal "+function.name,
						case_file(options.upstream_case, function.file, "upstream temporal table"));
			}
			for (const auto& boundary : three_d_configuration.boundaries)
				for (const auto& condition : boundary.conditions) {
					if (condition.waveform.empty()) continue;
					const auto& function = iga::FindTemporalFunction(three_d_configuration, condition.waveform);
					if (function.kind == iga::TemporalFunctionKind::PeriodicTable)
						selected_files.emplace("sequential 3D temporal "+function.name,
							case_file(options.three_d_case, function.file, "3D temporal table"));
				}
		});
		iga::RequireCollectiveAssetFiles(communicator, selected_files);
		auto [upstream_owner, downstream_owner, database_owner, mesh, boundary_velocity,
			wall_trace_basis, initial_three_d, outlet_models, initial_boundaries,
			upstream_terminal_runtime_id, downstream_terminal_runtime_id] =
			SequentialLocalValue(communicator, "sequential native inputs", [&] {
			if (options.graph_case_mode) {
				CanonicalizePeriodicTableAssets(upstream_configuration.temporal_functions,
					options.upstream_case, "upstream");
				CanonicalizePeriodicTableAssets(downstream_configuration.temporal_functions,
					options.downstream_case, "downstream");
			}
			RequireOneDFlowOnly(upstream_configuration, "upstream");
			RequireOneDFlowOnly(downstream_configuration, "downstream");
			const auto upstream_flow = SelectFlow(upstream_configuration);
			const auto downstream_flow = SelectFlow(downstream_configuration);
			const auto upstream_geometry_path = options.graph_case_mode
				? iga::ResolveContainedCaseFile(options.upstream_case,
					upstream_configuration.geometry.file, "upstream geometry")
				: options.upstream_case/upstream_configuration.geometry.file;
			const auto downstream_geometry_path = options.graph_case_mode
				? iga::ResolveContainedCaseFile(options.downstream_case,
					downstream_configuration.geometry.file, "downstream geometry")
				: options.downstream_case/downstream_configuration.geometry.file;
			auto upstream_network = iga::ReadOneDNetwork(upstream_geometry_path,
				upstream_configuration.geometry.length_scale_to_m, upstream_flow.discretization.cells_per_segment,
				upstream_flow.dynamic_viscosity, upstream_configuration.geometry.root_node_id);
			auto downstream_network = iga::ReadOneDNetwork(downstream_geometry_path,
				downstream_configuration.geometry.length_scale_to_m, downstream_flow.discretization.cells_per_segment,
				downstream_flow.dynamic_viscosity, downstream_configuration.geometry.root_node_id);
			iga::ValidateOneDTopologyReferences(upstream_configuration, upstream_network);
			iga::ValidateOneDTopologyReferences(downstream_configuration, downstream_network);
			const auto upstream_inlet = iga::ResolveOneDInlet(upstream_configuration);
			const auto downstream_inlet = iga::ResolveOneDInlet(downstream_configuration);
			auto upstream_owner = std::make_unique<iga::OneDFlowRuntime>(upstream_configuration, upstream_flow, std::move(upstream_network),
				upstream_inlet, options.upstream_case, [communicator](const iga::OneDNetwork& network,
					const iga::OneDFlowSystemDefinition& flow, iga::OneDFlowState& state, double inlet, double dt) {
					iga::AdvanceImplicitOneD(network, flow, state, inlet, dt, communicator);
				}, [communicator](const char* stage, std::exception_ptr error) {
					iga::RuntimeConstructionStage(communicator, stage, [&] {
						if (error) std::rethrow_exception(error);
					});
				});
			auto& upstream = *upstream_owner;
			auto downstream_owner = std::make_unique<iga::OneDFlowRuntime>(downstream_configuration, downstream_flow, std::move(downstream_network),
				downstream_inlet, options.downstream_case, [communicator](const iga::OneDNetwork& network,
					const iga::OneDFlowSystemDefinition& flow, iga::OneDFlowState& state, double inlet, double dt) {
					iga::AdvanceImplicitOneD(network, flow, state, inlet, dt, communicator);
				}, [communicator](const char* stage, std::exception_ptr error) {
					iga::RuntimeConstructionStage(communicator, stage, [&] {
						if (error) std::rethrow_exception(error);
					});
				});
			auto& downstream = *downstream_owner;
			if (upstream.FlowState().outlets.size() != 1)
				throw std::runtime_error("explicit 1D--3D coupling requires exactly one upstream terminal closure");
			auto upstream_terminal_runtime_id = "outlet:"+std::to_string(options.upstream_terminal_node);
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
			auto downstream_terminal_runtime_id = "outlet:"+std::to_string(
				downstream.Network().nodes[static_cast<std::size_t>(
					downstream.FlowState().outlets.front().node)].id);

			auto database_owner = std::make_unique<iga::Database>(options.database.string());
			auto& database = *database_owner;
			if (options.graph_case_mode)
				CanonicalizePeriodicTableAssets(three_d_configuration.temporal_functions,
					options.three_d_case, "3D");
			if (three_d_configuration.equation_systems.size() != 1
				|| three_d_configuration.equation_systems.front().kind != iga::EquationKind::NavierStokes
				|| three_d_configuration.equation_systems.front().unknowns.size() != 2)
				throw std::runtime_error("explicit 1D--3D coupling requires exactly one two-unknown Navier-Stokes 3D system");
			for (const auto& field : three_d_configuration.fields)
				if (field.kind == iga::FieldKind::Scalar)
					throw std::runtime_error("explicit 1D--3D coupling rejects 3D transport fields");
			const auto& three_d_flow = three_d_configuration.equation_systems.front();
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
			const auto mesh_path = options.graph_case_mode
				? iga::ResolveContainedCaseFile(options.three_d_case,
					"controlmesh.vtk", "3D control mesh")
				: options.three_d_case/"controlmesh.vtk";
			const auto velocity_path = options.graph_case_mode
				? iga::ResolveContainedCaseFile(options.three_d_case,
					"initial_velocityfield.txt", "3D initial velocity")
				: options.three_d_case/"initial_velocityfield.txt";
			auto mesh = iga::ReadLabeledHexMesh(mesh_path.string(),
				database.header().nodes, database.header().elements);
			auto boundary_velocity = iga::ReadVelocity(
				velocity_path.string(), database.header().nodes);
			auto wall_trace_basis = iga::WallTraceBasis(database, mesh, 0);
			auto initial_three_d = iga::MaterializeBoundaryWaveforms(three_d_configuration,
				options.three_d_case.string(), 0.0);
			auto outlet_models = iga::InitializeOutletModels(three_d_configuration, three_d_flow);
			if (!outlet_models.empty()) throw std::runtime_error("explicit 1D--3D coupling requires a static-pressure 3D outlet");
			auto initial_boundaries = iga::ResolveFlowBoundaries(initial_three_d,
				initial_three_d.equation_systems.front(), mesh.labels, boundary_velocity);
			return std::make_tuple(std::move(upstream_owner), std::move(downstream_owner),
				std::move(database_owner), std::move(mesh), std::move(boundary_velocity),
				std::move(wall_trace_basis), std::move(initial_three_d), std::move(outlet_models),
				std::move(initial_boundaries), std::move(upstream_terminal_runtime_id),
				std::move(downstream_terminal_runtime_id));
		});
		auto& upstream = *upstream_owner;
		auto& downstream = *downstream_owner;
		auto& database = *database_owner;
		const auto& three_d_flow = three_d_configuration.equation_systems.front();
		const auto& three_d_pressure_unknown = three_d_flow.unknowns.at(1);
		const auto& ports = three_d_configuration.coupling.three_d_ports;
		auto three_d_owner = iga::AllocateCollectiveRuntime<iga::TransientFlowRuntime>(communicator,
			database, communicator, true, true,
			iga::NavierStokesParameters{three_d_flow.density, three_d_flow.viscosity, three_d_configuration.time.dt},
			initial_boundaries, mesh.labels, boundary_velocity, wall_trace_basis, outlet_models, application_options);
		auto& three_d = *three_d_owner;
		iga::RequireValidGeometry(three_d.Elements(), rank, communicator);
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
		MPI_Allreduce(&local_unknown_label, &global_unknown_label, 1, MPI_INT, MPI_MAX, communicator);
		if (global_unknown_label)
			throw std::runtime_error("explicit coupling 3D database has a non-port non-wall boundary label");
		std::array<long long, 3> global_face_counts{{0, 0, 0}};
		MPI_Allreduce(local_face_counts.data(), global_face_counts.data(), 3, MPI_LONG_LONG, MPI_SUM,
			communicator);
		for (const auto count : global_face_counts)
			if (count == 0)
				throw std::runtime_error("explicit coupling 3D requires wall, inlet, and outlet .ntiga boundary faces");
		const double normalized_length = database.header().version == iga::kVersion
			? database.header().geometry_transform.source_units_per_normalized_unit
				*database.header().geometry_transform.source_length_scale_to_m : 0.0;
		const double native_reference_inlet_flow
			= three_d.ReferenceBoundaryFlow(ports.inlet_label);
		auto [reference_inlet_flow, scalar, upstream_subcycling, downstream_subcycling] =
			SequentialLocalValue(communicator, "sequential scalar preflight", [&] {
			const double reference_inlet_flow = graph_definition
				? graph_configuration->graph.Port(graph_definition->three_d_inlet)
					.orientation.ToOutward(native_reference_inlet_flow)
				: native_reference_inlet_flow;
			iga::ExplicitCouplingScalarPreflight scalar;
			scalar.dt_s = three_d_configuration.time.dt;
			scalar.steps = three_d_configuration.time.steps;
			scalar.density_kg_m3 = three_d_flow.density;
			scalar.dynamic_viscosity_pa_s = three_d_flow.viscosity;
			scalar.normalized_length_m = normalized_length;
			scalar.reference_inlet_outward_flow_m3_s = reference_inlet_flow;
			iga::ValidateExplicitCouplingScalarPreflight(scalar);
			if (graph_configuration) {
				const double time_scale = std::max({1.0, std::abs(scalar.dt_s),
					std::abs(graph_configuration->time.dt_s)});
				if (std::abs(scalar.dt_s-graph_configuration->time.dt_s)
					> 1.0e-12*time_scale || scalar.steps != graph_configuration->time.steps)
					throw std::runtime_error(
						"schema-v5 time must match the body-fitted 3D case time grid");
			}
			const auto upstream_subcycling = iga::MakeOneDSubcyclingPlan(upstream.Configuration().time.dt,
				upstream.Configuration().time.steps, scalar.dt_s, scalar.steps);
			const auto downstream_subcycling = iga::MakeOneDSubcyclingPlan(downstream.Configuration().time.dt,
				downstream.Configuration().time.steps, scalar.dt_s, scalar.steps);
			if (upstream.FlowSystem().density != scalar.density_kg_m3
				|| downstream.FlowSystem().density != scalar.density_kg_m3
				|| upstream.FlowSystem().dynamic_viscosity != scalar.dynamic_viscosity_pa_s
				|| downstream.FlowSystem().dynamic_viscosity != scalar.dynamic_viscosity_pa_s)
				throw std::runtime_error("explicit 1D--3D coupling requires identical density and viscosity");
			return std::make_tuple(reference_inlet_flow, scalar, upstream_subcycling, downstream_subcycling);
		});
		std::string upstream_domain_id;
		std::string three_d_domain_id;
		std::string downstream_domain_id;
		std::string upstream_edge_id;
		std::string downstream_edge_id;
		iga::PortRef upstream_terminal_ref;
		iga::PortRef upstream_root_ref;
		iga::PortRef three_d_inlet_ref;
		iga::PortRef three_d_outlet_ref;
		iga::PortRef three_d_wall_ref;
		iga::PortRef downstream_root_ref;
		iga::PortRef downstream_terminal_ref;
		iga::RuntimeConstructionStage(communicator, "sequential port identifiers", [&] {
			upstream_domain_id = "upstream";
			three_d_domain_id = "three_d";
			downstream_domain_id = "downstream";
			upstream_edge_id = "upstream_to_three_d";
			downstream_edge_id = "three_d_to_downstream";
			upstream_terminal_ref = {"upstream", "terminal"};
			upstream_root_ref = {"upstream", "root_state"};
			three_d_inlet_ref = {"three_d", "inlet"};
			three_d_outlet_ref = {"three_d", "outlet"};
			three_d_wall_ref = {"three_d", "wall"};
			downstream_root_ref = {"downstream", "root"};
			downstream_terminal_ref = {"downstream", "terminal_state"};
			if (graph_definition) {
				upstream_domain_id = graph_definition->upstream_domain_id;
				three_d_domain_id = graph_definition->three_d_domain_id;
				downstream_domain_id = graph_definition->downstream_domain_id;
				upstream_edge_id = graph_definition->upstream_edge_id;
				downstream_edge_id = graph_definition->downstream_edge_id;
				upstream_terminal_ref = graph_definition->upstream_terminal;
				upstream_root_ref = graph_definition->upstream_root_observation;
				three_d_inlet_ref = graph_definition->three_d_inlet;
				three_d_outlet_ref = graph_definition->three_d_outlet;
				three_d_wall_ref = graph_definition->three_d_wall_observation;
				downstream_root_ref = graph_definition->downstream_root;
				downstream_terminal_ref = graph_definition->downstream_terminal_observation;
			}
		});
		auto simulation_graph = SequentialLocalValue(communicator, "sequential graph binding", [&] {
			if (graph_configuration) return graph_configuration->graph;
			auto upstream_graph_port = MakeOneDLogicalPort("upstream", "terminal",
				upstream_terminal_runtime_id, iga::PortQuantity::MeanPressure);
			auto three_d_inlet_graph_port = MakeThreeDPort("three_d", "inlet",
				ports.inlet_label, {iga::PortQuantity::Area, iga::PortQuantity::FlowRate,
					iga::PortQuantity::MeanPressure}, {iga::PortQuantity::FlowRate});
			auto three_d_outlet_graph_port = MakeThreeDPort("three_d", "outlet",
				ports.outlet_labels.front(), {iga::PortQuantity::Area,
					iga::PortQuantity::FlowRate, iga::PortQuantity::MeanPressure},
				{iga::PortQuantity::MeanPressure});
			auto downstream_graph_port = MakeOneDLogicalPort("downstream", "root", "root",
				iga::PortQuantity::FlowRate);
			auto upstream_root_graph_port = MakeOneDObservationPort("upstream", "root_state",
				"root");
			auto downstream_terminal_graph_port = MakeOneDObservationPort("downstream",
				"terminal_state", downstream_terminal_runtime_id);
			auto three_d_wall_graph_port = MakeThreeDPort("three_d", "wall", 0,
				{iga::PortQuantity::FlowRate}, {});
			return iga::SimulationGraph({
				{"downstream", iga::DomainKind::OneDFlow, {std::move(downstream_graph_port),
					std::move(downstream_terminal_graph_port)}},
				{"three_d", iga::DomainKind::ThreeDBodyFittedFlow,
					{std::move(three_d_outlet_graph_port), std::move(three_d_inlet_graph_port),
					 std::move(three_d_wall_graph_port)}},
				{"upstream", iga::DomainKind::OneDFlow, {std::move(upstream_graph_port),
					std::move(upstream_root_graph_port)}}}, {
				{"upstream_to_three_d", {"upstream", "terminal"}, {"three_d", "inlet"},
					iga::CouplingLaw::PressureFlow},
				{"three_d_to_downstream", {"three_d", "outlet"}, {"downstream", "root"},
					iga::CouplingLaw::PressureFlow}});
		});
		iga::RuntimeConstructionStage(communicator, "sequential port validation", [&] {
			const auto sequential_plan = iga::MakeSequentialPressureFlowPlan(simulation_graph,
				upstream_domain_id);
			if (sequential_plan.domain_ids != std::vector<std::string>{upstream_domain_id,
				three_d_domain_id, downstream_domain_id})
				throw std::runtime_error("explicit coupling graph did not resolve the expected sequential chain");
			const auto& upstream_terminal_port = simulation_graph.Port(upstream_terminal_ref);
			const auto& three_d_inlet_port = simulation_graph.Port(three_d_inlet_ref);
			const auto& three_d_outlet_port = simulation_graph.Port(three_d_outlet_ref);
			const auto& downstream_root_port = simulation_graph.Port(downstream_root_ref);
			const auto& three_d_wall_measure = simulation_graph.Port(three_d_wall_ref);
			if (upstream_terminal_port.locator != upstream_terminal_runtime_id
				|| downstream_root_port.locator != "root"
				|| simulation_graph.Port(upstream_root_ref).locator != "root"
				|| simulation_graph.Port(downstream_terminal_ref).locator
					!= downstream_terminal_runtime_id)
				throw std::runtime_error("schema-v5 1D port locators do not match native cases");
			if (iga::ParseThreeDFlowBoundaryLabel(three_d_inlet_port) != ports.inlet_label
				|| iga::ParseThreeDFlowBoundaryLabel(three_d_outlet_port)
					!= ports.outlet_labels.front()
				|| iga::ParseThreeDFlowBoundaryLabel(three_d_wall_measure) != 0)
				throw std::runtime_error("schema-v5 3D port labels do not match native case coupling ports");
		});
		const auto& upstream_terminal_port = simulation_graph.Port(upstream_terminal_ref);
		const auto& three_d_inlet_port = simulation_graph.Port(three_d_inlet_ref);
		const auto& three_d_outlet_port = simulation_graph.Port(three_d_outlet_ref);
		const auto& downstream_root_port = simulation_graph.Port(downstream_root_ref);
		const auto& three_d_wall_measure = simulation_graph.Port(three_d_wall_ref);
		const auto& upstream_terminal = upstream_terminal_port.locator;
		const auto& downstream_root = downstream_root_port.locator;
		const auto& downstream_terminal = downstream_terminal_runtime_id;

		const double initial_outlet_pressure = SequentialLocalValue(communicator, "sequential initial state input", [&] {
			const double initial_upstream_flow = iga::EvaluateOneDInlet(upstream.Configuration(), upstream.InletDefinition(),
				options.upstream_case, 0.0, upstream.Network().segments.front().area0);
			if (!(initial_upstream_flow > 0.0) || !std::isfinite(initial_upstream_flow))
				throw std::runtime_error("explicit 1D--3D coupling requires a positive initial upstream native flow");
			upstream.InitializeOpenLoop(initial_upstream_flow);
			downstream.InitializeOpenLoop(initial_upstream_flow);
			const auto initial_downstream_root = downstream.GetPortState(downstream_root);
			iga::PortBoundaryData initial_profile_input;
			initial_profile_input.time_s = 0.0;
			initial_profile_input.outward_flow_m3_s = -initial_upstream_flow;
			iga::ApplyThreeDReferenceProfileInput(initial_three_d,
				initial_three_d.equation_systems.front(), three_d_inlet_port,
				initial_profile_input, reference_inlet_flow);
			const double initial_outlet_pressure = graph_configuration
				? graph_configuration->initial_pressure_pa.at(downstream_edge_id)
				: RequirePortValue(initial_downstream_root.mean_pressure_pa,
					"downstream root pressure");
			ApplyInitialPressure(initial_three_d, three_d_pressure_unknown, ports.outlet_labels.front(),
				initial_outlet_pressure);
			return initial_outlet_pressure;
		});
		three_d.InitializeState(initial_three_d);
		const auto initial_measurement_ports = SequentialLocalValue(communicator, "sequential initial measurements", [&] {
			return std::vector<iga::CouplingPort>{three_d_inlet_port, three_d_outlet_port};
		});
		const auto initial_three_d_ports = three_d.MeasurePorts(initial_measurement_ports, 0.0, {}, {});
		double lagged_three_d_inlet_pressure = SequentialLocalValue(communicator, "sequential initial lag", [&] {
			return graph_configuration
				? graph_configuration->initial_pressure_pa.at(upstream_edge_id)
				: RequirePortValue(initial_three_d_ports.at(three_d_inlet_port.id).mean_pressure_pa,
					"3D inlet pressure");
		});
		double lagged_downstream_root_pressure = initial_outlet_pressure;
		const double initial_lagged_three_d_inlet_pressure = lagged_three_d_inlet_pressure;
		const double initial_lagged_downstream_root_pressure = lagged_downstream_root_pressure;
		const int final_step = options.stop_after_step > 0 ? options.stop_after_step : scalar.steps;
		std::vector<iga::ExplicitCouplingHistoryRow> history;
		iga::RuntimeConstructionStage(communicator, "sequential history preparation", [&] {
			if (final_step > scalar.steps) throw std::runtime_error("--stop-after-step exceeds configured steps");
			history.reserve(static_cast<std::size_t>(final_step));
		});
		std::vector<iga::StrongCouplingIterationRow> strong_iterations;
		std::vector<long long> strong_accepted_ksp;
		std::vector<long long> strong_all_ksp;
		std::vector<int> strong_relaxation_updates;
		std::vector<double> strong_last_applied_relaxation;
		std::vector<double> strong_final_proposed_relaxation;
		std::string final_strong_diagnostics;
		if (options.strong_fixed || options.strong_aitken) {
			for (int step = 1; step <= final_step; ++step) {
				const double time = step*scalar.dt_s;
				double applied_upstream_pressure = lagged_three_d_inlet_pressure;
				double applied_three_d_pressure = lagged_downstream_root_pressure;
				long long all_ksp = 0;
				iga::OneDTrialWorkAccumulator upstream_work;
				iga::OneDTrialWorkAccumulator downstream_work;
				int relaxation_updates = 0;
				double last_applied_relaxation = 0.0;
				bool committed = false;
				auto aitken = SequentialLocalValue(communicator, "strong ready", [&] {
					iga::AitkenRelaxation<2> value(options.aitken_controls);
					value.Reset();
					iga::ObserveSequentialStrong("ready", step, upstream, three_d, downstream);
					return value;
				});
				std::array<std::exception_ptr, 3> abort_errors{};
				try {
					iga::SequentialRuntimeStage(communicator, "strong begin upstream", [&] {
						upstream.BeginStep(upstream.FlowState().physical_time, scalar.dt_s);
					});
					iga::SequentialRuntimeStage(communicator, "strong begin 3D", [&] {
						three_d.BeginStep(step-1, time, options.three_d_max_newton,
						kThreeDNonlinearRelativeTolerance, kThreeDNonlinearAbsoluteTolerance,
						kThreeDMassRelativeTolerance);
					});
					iga::SequentialRuntimeStage(communicator, "strong begin downstream", [&] {
						downstream.BeginStep(downstream.FlowState().physical_time, scalar.dt_s);
					});
					iga::SequentialRuntimeStage(communicator, "strong upstream schedule", [&] {
						upstream.SetConfiguredOpenLoopInlet();
					});
					for (int iteration = 1; iteration <= options.strong_controls.maximum_iterations; ++iteration) {
						iga::RuntimeConstructionStage(communicator, "strong upstream input", [&] {
							iga::PortBoundaryData upstream_pressure;
							upstream_pressure.time_s = time;
							upstream_pressure.mean_pressure_pa = applied_upstream_pressure;
							upstream.SetPortInput(upstream_terminal, upstream_pressure);
						});
						iga::SequentialRuntimeStage(communicator, "strong solve upstream", [&] {
							try { upstream.SolveTrial(); }
							catch (...) {
								upstream_work.Add(upstream.Diagnostics().attempted_configured_substeps,
									upstream.Diagnostics().explicit_cfl_substep_delta);
								std::ostringstream diagnostic;
							diagnostic.exceptions(std::ios::badbit | std::ios::failbit);
								diagnostic << std::setprecision(17) << "step=" << step << " iteration=" << iteration
									<< " x=(" << applied_upstream_pressure << ',' << applied_three_d_pressure << ')'
									<< " upstream SolveTrial failed ";
								WriteStrongOneDWorkTuple(diagnostic, "upstream", upstream.Diagnostics(), upstream_work);
								diagnostic << ' ';
								WriteStrongOneDWorkTuple(diagnostic, "downstream", downstream.Diagnostics(), downstream_work);
								final_strong_diagnostics = diagnostic.str();
								throw;
							}
						});
						iga::PortState upstream_port;
						double upstream_q = 0.0;
						iga::RuntimeConstructionStage(communicator, "strong 3D input", [&] {
							upstream_work.Add(upstream.Diagnostics().attempted_configured_substeps,
								upstream.Diagnostics().explicit_cfl_substep_delta);
							upstream_port = upstream.GetPortState(upstream_terminal);
							upstream_q = RequirePortValue(upstream_port.outward_flow_m3_s, "upstream terminal flow");

							auto three_d_step = iga::MaterializeBoundaryWaveforms(three_d_configuration,
								options.three_d_case.string(), time);
							iga::PortBoundaryData three_d_profile_input;
							three_d_profile_input.time_s = time;
							three_d_profile_input.outward_flow_m3_s = -upstream_q;
							iga::ApplyThreeDReferenceProfileInput(three_d_step,
								three_d_step.equation_systems.front(), three_d_inlet_port,
								three_d_profile_input, reference_inlet_flow);
							three_d.SetTrialBoundaryConfiguration(three_d_step);
							iga::PortBoundaryData three_d_pressure_input;
							three_d_pressure_input.time_s = time;
							three_d_pressure_input.mean_pressure_pa = applied_three_d_pressure;
							three_d.SetPortInput(three_d_outlet_port, three_d_pressure_input);
						});
						long long trial_ksp = 0;
						iga::SequentialRuntimeStage(communicator, "strong solve three_d", [&] {
							try {
								three_d.SolveTrial();
							} catch (...) {
								trial_ksp = static_cast<long long>(three_d.TrialLinearIterations());
								all_ksp += trial_ksp;
								std::ostringstream diagnostic;
							diagnostic.exceptions(std::ios::badbit | std::ios::failbit);
								diagnostic << std::setprecision(17) << "step=" << step << " iteration=" << iteration
									<< " x=(" << applied_upstream_pressure << ',' << applied_three_d_pressure << ')'
									<< " G=(unavailable,unavailable) signed_pressure_residual_pa=(unavailable,unavailable)"
									<< " normalized_pressure_residual=(unavailable,unavailable)"
									<< " flow_residuals_m3_s=(unavailable,unavailable) mass_imbalance_m3_s=unavailable"
									<< " three_d_ksp_attempt=" << trial_ksp << " three_d_ksp_cumulative=" << all_ksp << ' ';
								WriteStrongOneDWorkTuple(diagnostic, "upstream", upstream.Diagnostics(), upstream_work);
								diagnostic << ' ';
								WriteStrongOneDWorkTuple(diagnostic, "downstream", downstream.Diagnostics(), downstream_work);
								final_strong_diagnostics = diagnostic.str();
								throw;
							}
						});
						trial_ksp = static_cast<long long>(three_d.TrialLinearIterations());
						all_ksp += trial_ksp;
						iga::PortState three_d_inlet;
						iga::SequentialRuntimeStage(communicator, "strong measure inlet", [&] {
							three_d_inlet = three_d.GetPortState(three_d_inlet_port);
						});
						iga::PortState three_d_outlet;
						iga::SequentialRuntimeStage(communicator, "strong measure outlet", [&] {
							three_d_outlet = three_d.GetPortState(three_d_outlet_port);
						});
						iga::PortState three_d_wall;
						iga::SequentialRuntimeStage(communicator, "strong measure wall", [&] {
							three_d_wall = three_d.GetPortState(three_d_wall_measure);
						});
						double three_d_outlet_q = 0.0;
						iga::RuntimeConstructionStage(communicator, "strong downstream input", [&] {
							three_d_outlet_q = RequirePortValue(three_d_outlet.outward_flow_m3_s, "3D outlet flow");

							iga::PortBoundaryData downstream_input;
							downstream_input.time_s = time;
							downstream_input.outward_flow_m3_s = -three_d_outlet_q;
							downstream.SetPortInput(downstream_root, downstream_input);
						});
						iga::SequentialRuntimeStage(communicator, "strong solve downstream", [&] {
							try { downstream.SolveTrial(); }
							catch (...) {
								downstream_work.Add(downstream.Diagnostics().attempted_configured_substeps,
									downstream.Diagnostics().explicit_cfl_substep_delta);
								std::ostringstream diagnostic;
							diagnostic.exceptions(std::ios::badbit | std::ios::failbit);
								diagnostic << std::setprecision(17) << "step=" << step << " iteration=" << iteration
									<< " x=(" << applied_upstream_pressure << ',' << applied_three_d_pressure << ')'
									<< " downstream SolveTrial failed ";
								WriteStrongOneDWorkTuple(diagnostic, "upstream", upstream.Diagnostics(), upstream_work);
								diagnostic << ' ';
								WriteStrongOneDWorkTuple(diagnostic, "downstream", downstream.Diagnostics(), downstream_work);
								final_strong_diagnostics = diagnostic.str();
								throw;
							}
						});
						iga::ExplicitCouplingHistoryRow row;
						iga::StrongCouplingIterationRow iteration_row;
						std::array<double, 2> residual{};
						iga::RuntimeConstructionStage(communicator, "strong trial result", [&] {
							downstream_work.Add(downstream.Diagnostics().attempted_configured_substeps,
								downstream.Diagnostics().explicit_cfl_substep_delta);
							const auto downstream_root_state = downstream.GetPortState(downstream_root);
							const auto downstream_terminal_port = downstream.GetPortState(downstream_terminal);
							const auto upstream_root = upstream.GetPortState("root");

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
							row.downstream_root_pressure_pa = RequirePortValue(downstream_root_state.mean_pressure_pa, "downstream root pressure");
							row.downstream_root_outward_flow_m3_s = RequirePortValue(downstream_root_state.outward_flow_m3_s, "downstream root flow");
							row.downstream_root_area_m2 = RequirePortValue(downstream_root_state.area_m2, "downstream root area");
							row.downstream_terminal_pressure_pa = RequirePortValue(downstream_terminal_port.mean_pressure_pa, "downstream terminal pressure");
							row.downstream_terminal_outward_flow_m3_s = RequirePortValue(downstream_terminal_port.outward_flow_m3_s, "downstream terminal flow");
							row.downstream_terminal_area_m2 = RequirePortValue(downstream_terminal_port.area_m2, "downstream terminal area");
							row.upstream_three_d_flow_residual_m3_s = row.upstream_terminal_outward_flow_m3_s+row.three_d_inlet_outward_flow_m3_s;
							row.three_d_downstream_flow_residual_m3_s = row.three_d_outlet_outward_flow_m3_s+row.downstream_root_outward_flow_m3_s;
							row.upstream_three_d_normalized_residual = iga::ExplicitCouplingNormalizedResidual(row.upstream_terminal_outward_flow_m3_s, row.three_d_inlet_outward_flow_m3_s);
							row.three_d_downstream_normalized_residual = iga::ExplicitCouplingNormalizedResidual(row.three_d_outlet_outward_flow_m3_s, row.downstream_root_outward_flow_m3_s);
							row.three_d_wall_outward_flow_m3_s = RequirePortValue(three_d_wall.outward_flow_m3_s, "3D wall flow");
							row.three_d_mass_imbalance_m3_s = row.three_d_inlet_outward_flow_m3_s+row.three_d_outlet_outward_flow_m3_s+row.three_d_wall_outward_flow_m3_s;
							row.net_external_outward_flow_m3_s = row.upstream_root_outward_flow_m3_s+row.downstream_terminal_outward_flow_m3_s;
							row.external_pressure_drop_pa = row.upstream_root_pressure_pa-row.downstream_terminal_pressure_pa;
							row.upstream_three_d_pressure_jump_pa = row.upstream_terminal_pressure_pa-row.three_d_inlet_pressure_pa;
							row.three_d_downstream_pressure_jump_pa = row.three_d_outlet_pressure_pa-row.downstream_root_pressure_pa;
							row.iteration_count = iteration;
							row.relaxation_factor = options.strong_controls.relaxation_factor;
							row.three_d_trial_linear_iterations = trial_ksp;
							SetOneDAttemptWork(row, upstream, downstream);
							iga::ValidateExplicitCouplingHistoryRow(row);

							iteration_row.physical_step = step;
							iteration_row.time_s = time;
							iteration_row.iteration = iteration;
							iteration_row.applied_upstream_terminal_pressure_pa = applied_upstream_pressure;
							iteration_row.applied_three_d_outlet_traction_pressure_pa = applied_three_d_pressure;
							iteration_row.measured_three_d_inlet_pressure_pa = row.three_d_inlet_pressure_pa;
							iteration_row.measured_downstream_root_pressure_pa = row.downstream_root_pressure_pa;
							iteration_row.signed_upstream_pressure_residual_pa = row.three_d_inlet_pressure_pa-applied_upstream_pressure;
							iteration_row.signed_downstream_pressure_residual_pa = row.downstream_root_pressure_pa-applied_three_d_pressure;
							iteration_row.normalized_upstream_pressure_residual = iga::StrongCouplingPressureResidual(
								applied_upstream_pressure, row.three_d_inlet_pressure_pa,
								options.strong_controls.pressure_reference_pa);
							iteration_row.normalized_downstream_pressure_residual = iga::StrongCouplingPressureResidual(
								applied_three_d_pressure, row.downstream_root_pressure_pa,
								options.strong_controls.pressure_reference_pa);
							const std::array<double, 2> x{{applied_upstream_pressure, applied_three_d_pressure}};
							residual = {{iteration_row.signed_upstream_pressure_residual_pa,
								iteration_row.signed_downstream_pressure_residual_pa}};
							if (options.strong_aitken) {
								const auto proposal = aitken.Propose(x, residual, options.strong_controls.pressure_reference_pa);
								iteration_row.relaxation_factor_for_next_guess = proposal.relaxation_factor;
								iteration_row.unclamped_relaxation_factor = proposal.unclamped_relaxation_factor;
								iteration_row.aitken_scaled_numerator = proposal.scaled_numerator;
								iteration_row.aitken_scaled_denominator = proposal.scaled_denominator;
								iteration_row.aitken_status_code = iga::AitkenRelaxationStatusCode(proposal.status);
								iteration_row.aitken_has_previous_residual = proposal.has_previous_residual;
								iteration_row.next_upstream_terminal_pressure_pa = proposal.next[0];
								iteration_row.next_three_d_outlet_traction_pressure_pa = proposal.next[1];
							} else {
								iteration_row.relaxation_factor_for_next_guess = options.strong_controls.relaxation_factor;
								iteration_row.unclamped_relaxation_factor = options.strong_controls.relaxation_factor;
								iteration_row.next_upstream_terminal_pressure_pa = iga::StrongCouplingFixedUpdate(x[0], x[0]+residual[0], options.strong_controls.relaxation_factor);
								iteration_row.next_three_d_outlet_traction_pressure_pa = iga::StrongCouplingFixedUpdate(x[1], x[1]+residual[1], options.strong_controls.relaxation_factor);
							}
							iteration_row.upstream_terminal_outward_flow_m3_s = row.upstream_terminal_outward_flow_m3_s;
							iteration_row.three_d_inlet_outward_flow_m3_s = row.three_d_inlet_outward_flow_m3_s;
							iteration_row.three_d_outlet_outward_flow_m3_s = row.three_d_outlet_outward_flow_m3_s;
							iteration_row.downstream_root_outward_flow_m3_s = row.downstream_root_outward_flow_m3_s;
							iteration_row.upstream_three_d_flow_residual_m3_s = row.upstream_three_d_flow_residual_m3_s;
							iteration_row.three_d_downstream_flow_residual_m3_s = row.three_d_downstream_flow_residual_m3_s;
							iteration_row.normalized_upstream_three_d_flow_residual = row.upstream_three_d_normalized_residual;
							iteration_row.normalized_three_d_downstream_flow_residual = row.three_d_downstream_normalized_residual;
							iteration_row.three_d_wall_outward_flow_m3_s = row.three_d_wall_outward_flow_m3_s;
							iteration_row.three_d_mass_imbalance_m3_s = row.three_d_mass_imbalance_m3_s;
							iteration_row.three_d_attempt_linear_iterations = trial_ksp;
							iteration_row.three_d_cumulative_step_linear_iterations = all_ksp;
							iteration_row.upstream_1d_attempt_configured_substeps = upstream.Diagnostics().attempted_configured_substeps;
							iteration_row.upstream_1d_attempt_explicit_cfl_substeps = upstream.Diagnostics().explicit_cfl_substep_delta;
							iteration_row.upstream_1d_cumulative_configured_substeps = upstream_work.configured_substeps;
							iteration_row.upstream_1d_cumulative_explicit_cfl_substeps = upstream_work.explicit_cfl_substeps;
							iteration_row.downstream_1d_attempt_configured_substeps = downstream.Diagnostics().attempted_configured_substeps;
							iteration_row.downstream_1d_attempt_explicit_cfl_substeps = downstream.Diagnostics().explicit_cfl_substep_delta;
							iteration_row.downstream_1d_cumulative_configured_substeps = downstream_work.configured_substeps;
							iteration_row.downstream_1d_cumulative_explicit_cfl_substeps = downstream_work.explicit_cfl_substeps;
							iga::ValidateStrongCouplingIterationRow(iteration_row);
							std::ostringstream diagnostic;
							diagnostic.exceptions(std::ios::badbit | std::ios::failbit);
							diagnostic << std::setprecision(17) << "step=" << step << " iteration=" << iteration
								<< " x=(" << applied_upstream_pressure << ',' << applied_three_d_pressure << ')'
								<< " G=(" << row.three_d_inlet_pressure_pa << ',' << row.downstream_root_pressure_pa << ')'
								<< " signed_pressure_residual_pa=(" << iteration_row.signed_upstream_pressure_residual_pa
								<< ',' << iteration_row.signed_downstream_pressure_residual_pa << ')'
								<< " normalized_pressure_residual=(" << iteration_row.normalized_upstream_pressure_residual
								<< ',' << iteration_row.normalized_downstream_pressure_residual << ')'
								<< " flow_residuals_m3_s=(" << iteration_row.upstream_three_d_flow_residual_m3_s
								<< ',' << iteration_row.three_d_downstream_flow_residual_m3_s << ')'
								<< " mass_imbalance_m3_s=" << iteration_row.three_d_mass_imbalance_m3_s
								<< " three_d_ksp_attempt=" << trial_ksp << " three_d_ksp_cumulative=" << all_ksp;
							final_strong_diagnostics = diagnostic.str();
							const double maximum_flow_residual = std::max(
								std::abs(iteration_row.normalized_upstream_three_d_flow_residual),
								std::abs(iteration_row.normalized_three_d_downstream_flow_residual));
							if (maximum_flow_residual > options.strong_controls.flow_relative_tolerance)
							{
								iteration_row.relaxation_update_applied = false;
								strong_iterations.push_back(iteration_row);
								throw std::runtime_error("strong coupling transfer flow residual exceeds --strong-flow-relative-tol");
							}
							iteration_row.converged = iteration_row.normalized_upstream_pressure_residual <= options.strong_controls.pressure_relative_tolerance
								&& iteration_row.normalized_downstream_pressure_residual <= options.strong_controls.pressure_relative_tolerance
								&& std::abs(iteration_row.normalized_upstream_three_d_flow_residual) <= options.strong_controls.flow_relative_tolerance
								&& std::abs(iteration_row.normalized_three_d_downstream_flow_residual) <= options.strong_controls.flow_relative_tolerance;
						});
						iteration_row.converged = iga::SequentialStrongConverged(communicator, step,
							iteration, iteration_row.converged);
						if (iteration_row.converged) {
							iga::RuntimeConstructionStage(communicator, "strong precommit", [&] {
								iteration_row.relaxation_update_applied = false;
								row.relaxation_factor = relaxation_updates > 0 ? last_applied_relaxation : options.strong_controls.relaxation_factor;
								row.upstream_accepted_configured_substeps = iteration_row.upstream_1d_attempt_configured_substeps;
								row.upstream_all_configured_substeps = upstream_work.configured_substeps;
								row.upstream_rejected_configured_substeps = upstream_work.RejectedConfiguredSubsteps(row.upstream_accepted_configured_substeps);
								row.upstream_accepted_explicit_cfl_substeps = iteration_row.upstream_1d_attempt_explicit_cfl_substeps;
								row.upstream_all_explicit_cfl_substeps = upstream_work.explicit_cfl_substeps;
								row.upstream_rejected_explicit_cfl_substeps = upstream_work.RejectedExplicitCflSubsteps(row.upstream_accepted_explicit_cfl_substeps);
								row.downstream_accepted_configured_substeps = iteration_row.downstream_1d_attempt_configured_substeps;
								row.downstream_all_configured_substeps = downstream_work.configured_substeps;
								row.downstream_rejected_configured_substeps = downstream_work.RejectedConfiguredSubsteps(row.downstream_accepted_configured_substeps);
								row.downstream_accepted_explicit_cfl_substeps = iteration_row.downstream_1d_attempt_explicit_cfl_substeps;
								row.downstream_all_explicit_cfl_substeps = downstream_work.explicit_cfl_substeps;
								row.downstream_rejected_explicit_cfl_substeps = downstream_work.RejectedExplicitCflSubsteps(row.downstream_accepted_explicit_cfl_substeps);
								iga::ValidateExplicitCouplingHistoryRow(row);
								strong_iterations.push_back(iteration_row);
								if (injected_failure_step == step)
									throw std::runtime_error("injected strong coupling failure before commit");
							});
							iga::SequentialRuntimeStage(communicator, "strong prepare upstream", [&] {
								upstream.PrepareCommitStep();
							});
							iga::SequentialRuntimeStage(communicator, "strong prepare three_d", [&] {
								three_d.PrepareCommitStep();
							});
							iga::SequentialRuntimeStage(communicator, "strong prepare downstream", [&] {
								downstream.PrepareCommitStep();
							});
							upstream.FinalizeCommitStep();
							three_d.FinalizeCommitStep();
							downstream.FinalizeCommitStep();
							iga::RuntimeConstructionStage(communicator, "strong accepted result", [&] {
								history.push_back(row);
								strong_accepted_ksp.push_back(trial_ksp);
								strong_all_ksp.push_back(all_ksp);
								strong_relaxation_updates.push_back(relaxation_updates);
								strong_last_applied_relaxation.push_back(last_applied_relaxation);
								strong_final_proposed_relaxation.push_back(iteration_row.relaxation_factor_for_next_guess);
								lagged_three_d_inlet_pressure = row.three_d_inlet_pressure_pa;
								lagged_downstream_root_pressure = row.downstream_root_pressure_pa;
								committed = true;
								iga::ObserveSequentialStrong("accepted", step, upstream, three_d, downstream);
							});
							break;
						}
						iga::SequentialRuntimeStage(communicator, "strong rollback downstream", [&] {
							RollbackSolved(downstream);
						});
						iga::SequentialRuntimeStage(communicator, "strong rollback three_d", [&] {
							RollbackSolved(three_d);
						});
						iga::SequentialRuntimeStage(communicator, "strong rollback upstream", [&] {
							RollbackSolved(upstream);
						});
						iga::RuntimeConstructionStage(communicator, "strong relaxation update", [&] {
							if (options.strong_aitken) aitken.AcceptApplied(residual, iteration_row.relaxation_factor_for_next_guess);
							iteration_row.relaxation_update_applied = true;
							++relaxation_updates;
							last_applied_relaxation = iteration_row.relaxation_factor_for_next_guess;
							strong_iterations.push_back(iteration_row);
							applied_upstream_pressure = iteration_row.next_upstream_terminal_pressure_pa;
							applied_three_d_pressure = iteration_row.next_three_d_outlet_traction_pressure_pa;
						});
					}
					iga::RuntimeConstructionStage(communicator, "strong step completion", [&] {
						if (!committed)
							throw std::runtime_error("strong coupling did not converge within --strong-max-iterations: "+final_strong_diagnostics);
					});
				} catch (...) {
					const auto primary_error = std::current_exception();
					try {
						iga::SequentialRuntimeStage(communicator, "strong abort downstream", [&] { downstream.AbortStep(); });
					} catch (...) { abort_errors[0] = std::current_exception(); }
					try {
						iga::SequentialRuntimeStage(communicator, "strong abort three_d", [&] { three_d.AbortStep(); });
					} catch (...) { abort_errors[1] = std::current_exception(); }
					try {
						iga::SequentialRuntimeStage(communicator, "strong abort upstream", [&] { upstream.AbortStep(); });
					} catch (...) { abort_errors[2] = std::current_exception(); }
					iga::RuntimeConstructionStage(communicator, "strong failure diagnostics", [&] {
						iga::ObserveSequentialStrong("aborted", step, upstream, three_d, downstream);
						if (rank == 0 && !final_strong_diagnostics.empty()) {
							for (const auto& row : strong_iterations)
								if (row.physical_step == step)
									WriteStrongFailureIteration(std::cerr, row);
							std::cerr << "strong coupling final trial " << final_strong_diagnostics << '\n';
						}
						if (std::any_of(abort_errors.begin(), abort_errors.end(),
							[](const std::exception_ptr& error) { return bool(error); })) {
							std::ostringstream message;
							message.exceptions(std::ios::badbit | std::ios::failbit);
							const auto append = [&](const std::exception_ptr& error) {
								try { std::rethrow_exception(error); }
								catch (const std::exception& failure) { message << failure.what(); }
								catch (...) { message << "unknown exception"; }
							};
							message << "primary: "; append(primary_error);
							const std::array<const char*, 3> names{{"downstream", "three_d", "upstream"}};
							for (std::size_t i = 0; i < abort_errors.size(); ++i) if (abort_errors[i]) {
								message << "; abort " << names[i] << ": "; append(abort_errors[i]);
							}
							throw std::runtime_error(message.str());
						}
						std::rethrow_exception(primary_error);
					});
				}
			}
		} else {
			std::vector<std::unique_ptr<iga::CoupledDomainRuntime>> runtimes;
			iga::RuntimeConstructionStage(communicator, "sequential adapter catalog", [&] {
				runtimes.reserve(3);
				runtimes.push_back(std::make_unique<iga::OneDFlowDomainAdapter>(upstream_domain_id,
					upstream, simulation_graph.Domain(upstream_domain_id).ports,
					iga::OneDInletPolicy::ConfiguredOpenLoop));
				iga::ThreeDFlowDomainControls three_d_controls;
				three_d_controls.maximum_newton = options.three_d_max_newton;
				three_d_controls.nonlinear_relative_tolerance = kThreeDNonlinearRelativeTolerance;
				three_d_controls.nonlinear_absolute_tolerance = kThreeDNonlinearAbsoluteTolerance;
				three_d_controls.mass_relative_tolerance = kThreeDMassRelativeTolerance;
				runtimes.push_back(std::make_unique<iga::ThreeDBodyFittedFlowDomainAdapter>(
					three_d_domain_id, three_d, simulation_graph.Domain(three_d_domain_id).ports,
					three_d_configuration, options.three_d_case,
					std::map<std::string, double>{{three_d_inlet_port.id, reference_inlet_flow}},
					three_d_controls));
				runtimes.push_back(std::make_unique<iga::OneDFlowDomainAdapter>(downstream_domain_id,
					downstream, simulation_graph.Domain(downstream_domain_id).ports,
					iga::OneDInletPolicy::CoupledRoot));
			});
			auto registry_owner = iga::CreateCollectiveDomainRuntimeRegistry(communicator,
				simulation_graph, runtimes);
			auto& runtime_registry = *registry_owner;
			iga::PressureFlowExecutionControls execution_controls;
			std::unique_ptr<iga::PressureFlowComponentExecutor> executor;
			iga::CollectiveLocalStage(communicator, "sequential executor construction", [&] {
				executor = std::make_unique<iga::PressureFlowComponentExecutor>(runtime_registry,
					upstream_domain_id, execution_controls,
					iga::CollectivePressureFlowExecution(communicator));
			});
			for (int step = 1; step <= final_step; ++step) {
				const double time = step*scalar.dt_s;
				iga::ExplicitCouplingHistoryRow row;
				bool row_ready = false;
				std::map<std::string, double> pressure_input;
				std::function<void(const iga::PressureFlowStepResult&)> before_commit;
				iga::CollectiveLocalStage(communicator, "sequential executor step input", [&] {
					pressure_input = {{upstream_edge_id, lagged_three_d_inlet_pressure},
						{downstream_edge_id, lagged_downstream_root_pressure}};
					before_commit = [&](const iga::PressureFlowStepResult& result) {
						const auto& upstream_root = result.accepted_ports.at(upstream_root_ref);
						const auto& upstream_port = result.accepted_ports.at(upstream_terminal_ref);
						const auto& three_d_inlet = result.accepted_ports.at(three_d_inlet_ref);
						const auto& three_d_outlet = result.accepted_ports.at(three_d_outlet_ref);
						const auto& three_d_wall = result.accepted_ports.at(three_d_wall_ref);
						const auto& downstream_root_state = result.accepted_ports.at(downstream_root_ref);
						const auto& downstream_terminal_port
							= result.accepted_ports.at(downstream_terminal_ref);
						row.time_s = time;
						row.upstream_root_pressure_pa = RequirePortValue(upstream_root.mean_pressure_pa, "upstream root pressure");
						row.upstream_root_outward_flow_m3_s = RequirePortValue(upstream_root.outward_flow_m3_s, "upstream root flow");
						row.upstream_root_area_m2 = RequirePortValue(upstream_root.area_m2, "upstream root area");
						row.upstream_terminal_pressure_pa = RequirePortValue(upstream_port.mean_pressure_pa, "upstream terminal pressure");
						row.upstream_terminal_outward_flow_m3_s = RequirePortValue(upstream_port.outward_flow_m3_s, "upstream terminal flow");
						row.upstream_terminal_area_m2 = RequirePortValue(upstream_port.area_m2, "upstream terminal area");
						row.three_d_inlet_pressure_pa = RequirePortValue(three_d_inlet.mean_pressure_pa, "3D inlet pressure");
						row.three_d_inlet_outward_flow_m3_s = RequirePortValue(three_d_inlet.outward_flow_m3_s, "3D inlet flow");
						row.three_d_inlet_area_m2 = RequirePortValue(three_d_inlet.area_m2, "3D inlet area");
						row.three_d_outlet_pressure_pa = RequirePortValue(three_d_outlet.mean_pressure_pa, "3D outlet pressure");
						row.three_d_outlet_outward_flow_m3_s = RequirePortValue(three_d_outlet.outward_flow_m3_s, "3D outlet flow");
						row.three_d_outlet_area_m2 = RequirePortValue(three_d_outlet.area_m2, "3D outlet area");
						row.downstream_root_pressure_pa = RequirePortValue(downstream_root_state.mean_pressure_pa, "downstream root pressure");
						row.downstream_root_outward_flow_m3_s = RequirePortValue(downstream_root_state.outward_flow_m3_s, "downstream root flow");
						row.downstream_root_area_m2 = RequirePortValue(downstream_root_state.area_m2, "downstream root area");
						row.downstream_terminal_pressure_pa = RequirePortValue(downstream_terminal_port.mean_pressure_pa, "downstream terminal pressure");
						row.downstream_terminal_outward_flow_m3_s = RequirePortValue(downstream_terminal_port.outward_flow_m3_s, "downstream terminal flow");
						row.downstream_terminal_area_m2 = RequirePortValue(downstream_terminal_port.area_m2, "downstream terminal area");
						row.upstream_three_d_flow_residual_m3_s = row.upstream_terminal_outward_flow_m3_s+row.three_d_inlet_outward_flow_m3_s;
						row.three_d_downstream_flow_residual_m3_s = row.three_d_outlet_outward_flow_m3_s+row.downstream_root_outward_flow_m3_s;
						row.upstream_three_d_normalized_residual = iga::ExplicitCouplingNormalizedResidual(row.upstream_terminal_outward_flow_m3_s, row.three_d_inlet_outward_flow_m3_s);
						row.three_d_downstream_normalized_residual = iga::ExplicitCouplingNormalizedResidual(row.three_d_outlet_outward_flow_m3_s, row.downstream_root_outward_flow_m3_s);
						row.three_d_wall_outward_flow_m3_s = RequirePortValue(three_d_wall.outward_flow_m3_s, "3D wall flow");
						row.three_d_mass_imbalance_m3_s = row.three_d_inlet_outward_flow_m3_s+row.three_d_outlet_outward_flow_m3_s+row.three_d_wall_outward_flow_m3_s;
						row.net_external_outward_flow_m3_s = row.upstream_root_outward_flow_m3_s+row.downstream_terminal_outward_flow_m3_s;
						row.external_pressure_drop_pa = row.upstream_root_pressure_pa-row.downstream_terminal_pressure_pa;
						row.upstream_three_d_pressure_jump_pa = row.upstream_terminal_pressure_pa-row.three_d_inlet_pressure_pa;
						row.three_d_downstream_pressure_jump_pa = row.three_d_outlet_pressure_pa-row.downstream_root_pressure_pa;
						row.three_d_trial_linear_iterations = three_d.TrialLinearIterations();
						SetOneDAttemptWork(row, upstream, downstream);
						iga::ValidateExplicitCouplingHistoryRow(row);
						row_ready = true;
						if (injected_failure_step == step)
							throw std::runtime_error("injected explicit coupling failure before commit");
					};
				});
				executor->Advance({step-1, upstream.FlowState().physical_time, scalar.dt_s},
					pressure_input, before_commit);
				iga::CollectiveLocalStage(communicator, "sequential accepted flow result", [&] {
					if (!row_ready) throw std::runtime_error("explicit graph executor omitted its precommit observation");
					history.push_back(row);
					lagged_three_d_inlet_pressure = row.three_d_inlet_pressure_pa;
					lagged_downstream_root_pressure = row.downstream_root_pressure_pa;
				});
			}
		}
		three_d.Close();
		iga::CollectiveLocalStage(communicator, "sequential output", [&] {
			if (rank != 0) return;
			if (graph_configuration) {
				if (!fs::create_directories(options.output_directory))
					throw std::runtime_error(
						"schema-v5 graph output directory must be newly created");
			} else {
				fs::create_directories(options.output_directory);
			}
			if (!options.strong_fixed && !options.strong_aitken) {
				std::ofstream output(options.output_directory/"explicit_coupling_history.csv");
				if (!output) throw std::runtime_error("cannot create explicit coupling history");
				iga::WriteExplicitCouplingHistoryHeader(output);
				for (const auto& row : history) iga::WriteExplicitCouplingHistoryRow(output, row);
				output.close();
				if (!output) throw std::runtime_error("cannot write explicit coupling history");
				WriteExplicitCouplingManifest(options.output_directory/"explicit_coupling_manifest.json",
					options.upstream_terminal_node, ports.inlet_label, ports.outlet_labels.front(), scalar.dt_s,
					scalar.steps, static_cast<int>(history.size()), scalar.density_kg_m3,
					scalar.dynamic_viscosity_pa_s, normalized_length,
					reference_inlet_flow, initial_lagged_three_d_inlet_pressure,
					initial_lagged_downstream_root_pressure, upstream_subcycling, downstream_subcycling, history,
					options.three_d_max_newton);
			} else {
				std::ofstream history_output(options.output_directory/"strong_coupling_history.csv");
				std::ofstream iteration_output(options.output_directory/"strong_coupling_iterations.csv");
				if (!history_output || !iteration_output)
					throw std::runtime_error("cannot create strong coupling history output");
				history_output << std::setprecision(17);
				std::ostringstream serialized_header;
				serialized_header.exceptions(std::ios::badbit | std::ios::failbit);
				iga::WriteExplicitCouplingHistoryHeader(serialized_header);
				history_output << CsvLineWithoutNewline(serialized_header.str())
					<< ",final_signed_upstream_pressure_residual_pa,final_signed_downstream_pressure_residual_pa,"
					"final_normalized_upstream_pressure_residual,final_normalized_downstream_pressure_residual,"
					"final_max_normalized_pressure_residual,pressure_reference_pa,"
					"accepted_three_d_ksp_iterations,all_three_d_ksp_iterations,rejected_three_d_ksp_iterations,"
					"relaxation_updates_applied,last_applied_relaxation_factor,final_proposed_relaxation_factor,final_proposal_applied\n";
				for (std::size_t index = 0; index < history.size(); ++index) {
					const auto final_iteration = std::find_if(strong_iterations.rbegin(), strong_iterations.rend(),
						[index](const iga::StrongCouplingIterationRow& row) {
							return row.physical_step == static_cast<int>(index+1);
						});
					if (final_iteration == strong_iterations.rend())
						throw std::runtime_error("missing final strong coupling iteration row");
					std::ostringstream serialized_row;
					serialized_row.exceptions(std::ios::badbit | std::ios::failbit);
					iga::WriteExplicitCouplingHistoryRow(serialized_row, history[index]);
					history_output << CsvLineWithoutNewline(serialized_row.str()) << ','
						<< final_iteration->signed_upstream_pressure_residual_pa << ','
						<< final_iteration->signed_downstream_pressure_residual_pa << ','
						<< final_iteration->normalized_upstream_pressure_residual << ','
						<< final_iteration->normalized_downstream_pressure_residual << ','
						<< std::max(final_iteration->normalized_upstream_pressure_residual,
							final_iteration->normalized_downstream_pressure_residual) << ','
						<< options.strong_controls.pressure_reference_pa << ','
						<< strong_accepted_ksp.at(index) << ',' << strong_all_ksp.at(index) << ','
						<< strong_all_ksp.at(index)-strong_accepted_ksp.at(index) << ','
						<< strong_relaxation_updates.at(index) << ',' << strong_last_applied_relaxation.at(index) << ','
						<< strong_final_proposed_relaxation.at(index) << ",0\n";
				}
				iga::WriteStrongCouplingIterationHeader(iteration_output);
				for (const auto& row : strong_iterations) iga::WriteStrongCouplingIterationRow(iteration_output, row);
				history_output.close();
				iteration_output.close();
				if (!history_output || !iteration_output)
					throw std::runtime_error("cannot write strong coupling history output");
				const long long accepted_ksp = std::accumulate(strong_accepted_ksp.begin(), strong_accepted_ksp.end(), 0LL);
				const long long all_ksp = std::accumulate(strong_all_ksp.begin(), strong_all_ksp.end(), 0LL);
				std::array<long long, 6> aitken_status_counts{{0, 0, 0, 0, 0, 0}};
				for (const auto& row : strong_iterations)
					if (row.aitken_status_code >= 0 && row.aitken_status_code < static_cast<int>(aitken_status_counts.size()))
						++aitken_status_counts[static_cast<std::size_t>(row.aitken_status_code)];
				WriteStrongCouplingManifest(options.output_directory/"strong_coupling_manifest.json",
					options.strong_controls, options.aitken_controls, options.strong_aitken, scalar.steps, static_cast<int>(history.size()), scalar.dt_s,
					scalar.density_kg_m3, scalar.dynamic_viscosity_pa_s, normalized_length,
					options.upstream_terminal_node, ports.inlet_label, ports.outlet_labels.front(),
					initial_lagged_three_d_inlet_pressure, initial_lagged_downstream_root_pressure,
					static_cast<long long>(strong_iterations.size()), reference_inlet_flow, aitken_status_counts, accepted_ksp, all_ksp,
					upstream_subcycling, downstream_subcycling, history, options.three_d_max_newton);
			}
			if (graph_configuration)
				WriteGraphBindingManifest(options.output_directory/"graph_binding_manifest.json",
					*graph_configuration, graph_assets, graph_case_root);
		});
		iga::CollectiveLocalStage(communicator, "sequential completion logging", [&] {
			if (rank == 0) std::cout << "completed " << (options.strong_aitken ? "strong Aitken" : (options.strong_fixed ? "strong fixed" : "explicit"))
				<< " 1D--3D coupling steps=" << history.size()
				<< " output=" << options.output_directory << '\n';
			iga::FlushCheckedText(std::cout);
		});
	} catch (const std::exception& error) {
		if (rank == 0) std::cerr << error.what() << '\n';
		status = 1;
	}
	return status;
}

#ifndef IGA_SEQUENTIAL_NO_MAIN
int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, "TubularFlowIGA explicit 1D--3D flow-only coupling\n");
	const int status = iga::RunSequentialFlow(argc, argv, PETSC_COMM_WORLD);
	PetscFinalize();
	return status;
}
#endif
