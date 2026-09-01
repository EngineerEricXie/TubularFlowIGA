#include "BoundarySupport.hpp"
#include "AitkenRelaxation.hpp"
#include "ExplicitOneDThreeDCoupling.hpp"
#include "IgaDatabase.hpp"
#include "OneDImplicit.hpp"
#include "OneDRuntime.hpp"
#include "StrongOneDThreeDCoupling.hpp"
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
#include <numeric>
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
	bool strong_fixed = false;
	bool strong_aitken = false;
	iga::StrongCouplingControls strong_controls;
	iga::AitkenRelaxationControls aitken_controls;
	bool strong_pressure_reference_set = false;
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
	if (argc < 5) throw std::runtime_error(
		"usage: iga_1d_3d_explicit DB THREE_D_CASE UPSTREAM_1D_CASE DOWNSTREAM_1D_CASE "
		"--upstream-terminal-node ID --output-dir DIR [--stop-after-step N] "
		"[--coupling-mode explicit|strong-fixed|strong-aitken --strong-max-iterations N "
		"--strong-pressure-relative-tol R --strong-pressure-reference-pa PA "
		"--strong-flow-relative-tol R --strong-relaxation W "
		"--strong-aitken-min-relaxation W --strong-aitken-max-relaxation W] [PETSc options]");
	Options options;
	options.database = argv[1];
	options.three_d_case = argv[2];
	options.upstream_case = argv[3];
	options.downstream_case = argv[4];
	for (int i = 5; i < argc; ++i) {
		const std::string argument(argv[i]);
		if (argument == "--upstream-terminal-node" || argument == "--output-dir"
			|| argument == "--stop-after-step" || argument == "--coupling-mode"
			|| argument == "--strong-max-iterations" || argument == "--strong-pressure-relative-tol"
			|| argument == "--strong-pressure-reference-pa" || argument == "--strong-flow-relative-tol"
			|| argument == "--strong-relaxation" || argument == "--strong-aitken-min-relaxation"
			|| argument == "--strong-aitken-max-relaxation") {
			if (++i >= argc) throw std::runtime_error(argument+" requires a value");
			const std::string value(argv[i]);
			if (argument == "--upstream-terminal-node") options.upstream_terminal_node = PositiveInteger(value, argument);
			else if (argument == "--output-dir") options.output_directory = value;
			else if (argument == "--stop-after-step") options.stop_after_step = PositiveInteger(value, argument);
			if (argument == "--coupling-mode") {
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
	if (options.upstream_terminal_node < 0 || options.output_directory.empty())
		throw std::runtime_error("--upstream-terminal-node and --output-dir are required");
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

void WriteStrongCouplingManifest(const fs::path& path, const iga::StrongCouplingControls& controls,
	const iga::AitkenRelaxationControls& aitken_controls, bool strong_aitken,
	int configured_steps, int completed_steps, double dt_s, double density_kg_m3,
	double dynamic_viscosity_pa_s, double normalized_length_m, int upstream_terminal_node,
	int three_d_inlet_label, int three_d_outlet_label, double initial_upstream_pressure_pa,
	double initial_three_d_outlet_traction_pressure_pa, long long total_coupling_iterations,
	double reference_inlet_outward_flow_m3_s, const std::array<long long, 6>& aitken_status_counts,
	long long accepted_ksp, long long all_ksp)
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
		<< "  \"newton_controls\": {\"maximum_iterations\": " << kThreeDMaximumNewtonIterations
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
		<< "  \"work_semantics\": \"accepted 3D KSP work is the final committed attempt of each step; all_attempts includes rejected strong sweeps; rejected is all_attempts minus accepted. Rejected 1D computational work is not reported.\",\n"
		<< "  \"total_coupling_iterations\": " << total_coupling_iterations << ",\n"
		<< "  \"three_d_ksp_iterations\": {\"accepted\": " << accepted_ksp
		<< ", \"all_attempts\": " << all_ksp << ", \"rejected\": "
		<< all_ksp-accepted_ksp << "}\n}\n";
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
				int relaxation_updates = 0;
				double last_applied_relaxation = 0.0;
				bool committed = false;
				iga::AitkenRelaxation<2> aitken(options.aitken_controls);
				aitken.Reset();
				try {
					upstream.BeginStep(upstream.FlowState().physical_time, scalar.dt_s);
					three_d.BeginStep(step-1, time, kThreeDMaximumNewtonIterations,
						kThreeDNonlinearRelativeTolerance, kThreeDNonlinearAbsoluteTolerance,
						kThreeDMassRelativeTolerance);
					downstream.BeginStep(downstream.FlowState().physical_time, scalar.dt_s);
					for (int iteration = 1; iteration <= options.strong_controls.maximum_iterations; ++iteration) {
						const double upstream_flow = iga::EvaluateOneDInlet(upstream.Configuration(), upstream.InletDefinition(),
							options.upstream_case, time, upstream.Network().segments.front().area0);
						upstream.SetOpenLoopInlet(upstream.OpenLoopInlet(time, upstream_flow));
						iga::PortBoundaryData upstream_pressure;
						upstream_pressure.time_s = time;
						upstream_pressure.mean_pressure_pa = applied_upstream_pressure;
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
						three_d.SetTrialBoundaryConfiguration(three_d_step);
						iga::PortBoundaryData three_d_pressure_input;
						three_d_pressure_input.time_s = time;
						three_d_pressure_input.mean_pressure_pa = applied_three_d_pressure;
						three_d.SetPortInput(three_d_outlet_pressure, three_d_pressure_input);
						long long trial_ksp = 0;
						try {
							three_d.SolveTrial();
						} catch (...) {
							trial_ksp = static_cast<long long>(three_d.TrialLinearIterations());
							all_ksp += trial_ksp;
							std::ostringstream diagnostic;
							diagnostic << std::setprecision(17) << "step=" << step << " iteration=" << iteration
								<< " x=(" << applied_upstream_pressure << ',' << applied_three_d_pressure << ')'
								<< " G=(unavailable,unavailable) signed_pressure_residual_pa=(unavailable,unavailable)"
								<< " normalized_pressure_residual=(unavailable,unavailable)"
								<< " flow_residuals_m3_s=(unavailable,unavailable) mass_imbalance_m3_s=unavailable"
								<< " three_d_ksp_attempt=" << trial_ksp << " three_d_ksp_cumulative=" << all_ksp;
							final_strong_diagnostics = diagnostic.str();
							throw;
						}
						trial_ksp = static_cast<long long>(three_d.TrialLinearIterations());
						all_ksp += trial_ksp;
						const auto three_d_inlet = three_d.GetPortState(three_d_inlet_measure);
						const auto three_d_outlet = three_d.GetPortState(three_d_outlet_measure);
						const auto three_d_wall = three_d.GetPortState(three_d_wall_measure);
						const double three_d_outlet_q = RequirePortValue(three_d_outlet.outward_flow_m3_s, "3D outlet flow");

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
						row.three_d_wall_outward_flow_m3_s = RequirePortValue(three_d_wall.outward_flow_m3_s, "3D wall flow");
						row.three_d_mass_imbalance_m3_s = row.three_d_inlet_outward_flow_m3_s+row.three_d_outlet_outward_flow_m3_s+row.three_d_wall_outward_flow_m3_s;
						row.net_external_outward_flow_m3_s = row.upstream_root_outward_flow_m3_s+row.downstream_terminal_outward_flow_m3_s;
						row.external_pressure_drop_pa = row.upstream_root_pressure_pa-row.downstream_terminal_pressure_pa;
						row.upstream_three_d_pressure_jump_pa = row.upstream_terminal_pressure_pa-row.three_d_inlet_pressure_pa;
						row.three_d_downstream_pressure_jump_pa = row.three_d_outlet_pressure_pa-row.downstream_root_pressure_pa;
						row.iteration_count = iteration;
						row.relaxation_factor = options.strong_controls.relaxation_factor;
						row.three_d_trial_linear_iterations = trial_ksp;
						iga::ValidateExplicitCouplingHistoryRow(row);

						iga::StrongCouplingIterationRow iteration_row;
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
						const std::array<double, 2> residual{{iteration_row.signed_upstream_pressure_residual_pa,
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
						iga::ValidateStrongCouplingIterationRow(iteration_row);
						std::ostringstream diagnostic;
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
						if (iteration_row.converged) {
							iteration_row.relaxation_update_applied = false;
							row.relaxation_factor = relaxation_updates > 0 ? last_applied_relaxation : options.strong_controls.relaxation_factor;
							strong_iterations.push_back(iteration_row);
							if (injected_failure_step == step)
								throw std::runtime_error("injected strong coupling failure before commit");
							upstream.CommitStep();
							three_d.CommitStep();
							downstream.CommitStep();
							history.push_back(row);
							strong_accepted_ksp.push_back(trial_ksp);
							strong_all_ksp.push_back(all_ksp);
							strong_relaxation_updates.push_back(relaxation_updates);
							strong_last_applied_relaxation.push_back(last_applied_relaxation);
							strong_final_proposed_relaxation.push_back(iteration_row.relaxation_factor_for_next_guess);
							lagged_three_d_inlet_pressure = row.three_d_inlet_pressure_pa;
							lagged_downstream_root_pressure = row.downstream_root_pressure_pa;
							committed = true;
							break;
						}
						RollbackSolved(downstream);
						RollbackSolved(three_d);
						RollbackSolved(upstream);
						if (options.strong_aitken) aitken.AcceptApplied(residual, iteration_row.relaxation_factor_for_next_guess);
						iteration_row.relaxation_update_applied = true;
						++relaxation_updates;
						last_applied_relaxation = iteration_row.relaxation_factor_for_next_guess;
						strong_iterations.push_back(iteration_row);
						applied_upstream_pressure = iteration_row.next_upstream_terminal_pressure_pa;
						applied_three_d_pressure = iteration_row.next_three_d_outlet_traction_pressure_pa;
					}
					if (!committed)
						throw std::runtime_error("strong coupling did not converge within --strong-max-iterations: "+final_strong_diagnostics);
				} catch (...) {
					if (rank == 0 && !final_strong_diagnostics.empty()) {
						for (const auto& row : strong_iterations)
							if (row.physical_step == step)
								WriteStrongFailureIteration(std::cerr, row);
						std::cerr << "strong coupling final trial " << final_strong_diagnostics << '\n';
					}
					RollbackSolved(downstream);
					RollbackSolved(three_d);
					RollbackSolved(upstream);
					throw;
				}
			}
		} else for (int step = 1; step <= final_step; ++step) {
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
			if (!options.strong_fixed && !options.strong_aitken) {
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
			} else {
				std::ofstream history_output(options.output_directory/"strong_coupling_history.csv");
				std::ofstream iteration_output(options.output_directory/"strong_coupling_iterations.csv");
				if (!history_output || !iteration_output)
					throw std::runtime_error("cannot create strong coupling history output");
				history_output << std::setprecision(17);
				std::ostringstream serialized_header;
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
					static_cast<long long>(strong_iterations.size()), reference_inlet_flow, aitken_status_counts, accepted_ksp, all_ksp);
			}
		} catch (const std::exception& error) {
			output_failed = 1;
			output_error = error.what();
		}
		MPI_Bcast(&output_failed, 1, MPI_INT, 0, PETSC_COMM_WORLD);
		if (output_failed) {
			if (rank == 0) throw std::runtime_error(output_error);
			throw std::runtime_error("explicit coupling output failed on rank 0");
		}
		if (rank == 0) std::cout << "completed " << (options.strong_aitken ? "strong Aitken" : (options.strong_fixed ? "strong fixed" : "explicit"))
			<< " 1D--3D coupling steps=" << history.size()
			<< " output=" << options.output_directory << '\n';
	} catch (const std::exception& error) {
		if (rank == 0) std::cerr << error.what() << '\n';
		status = 1;
	}
	PetscFinalize();
	return status;
}
