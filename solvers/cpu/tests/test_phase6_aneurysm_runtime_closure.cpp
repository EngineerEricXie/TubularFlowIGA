#include "ImmersedFlowCase.hpp"
#include "SurfaceReaders.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr double kTargetFlowM3S = 1.0e-4;
constexpr double kControllerRelativeTolerance = 1.0e-10;
constexpr double kOpenBalanceTolerance = 1.0e-3;
constexpr double kWallLeakageTolerance = 1.0e-3;
constexpr double kTrueLinearRelativeTolerance = 1.0e-10;

void Require(bool condition, const std::string& message)
{
	if (!condition) throw std::runtime_error(message);
}

iga::CouplingPort Port(const char* id, int label, iga::PortQuantity requirement)
{
	iga::CouplingPort port;
	port.id = id;
	port.subsystem_id = "immersed";
	port.locator_kind = "boundary_label";
	port.locator = std::to_string(label);
	port.provides = {iga::PortQuantity::Area, iga::PortQuantity::FlowRate,
		iga::PortQuantity::MeanPressure};
	port.requires = {requirement};
	return port;
}

std::vector<iga::CouplingPort> Ports()
{
	return {Port("inlet", 1, iga::PortQuantity::FlowRate),
		Port("outlet", 2, iga::PortQuantity::MeanPressure)};
}

const iga::ImmersedStaticFlowDiagnostics::Port& DiagnosticPort(
	const iga::ImmersedStaticFlowDiagnostics& diagnostics, const char* id)
{
	const auto found = std::find_if(diagnostics.ports.begin(), diagnostics.ports.end(),
		[id](const auto& port) { return port.id == id; });
	if (found == diagnostics.ports.end()) throw std::runtime_error("missing immersed port diagnostic");
	return *found;
}

struct CapAudit {
	double maximum_area_error_m2 = 0.0;
	double maximum_normal_component_error = 0.0;
	std::size_t inlet_triangles = 0;
	std::size_t outlet_triangles = 0;
};

CapAudit AuditCaps(const fs::path& directory, const iga::ImmersedSurfaceQuadratureDiagnostics& diagnostics)
{
	constexpr double expected_area_m2 = 0.17677669529663687;
	CapAudit audit;
	for (const std::uint32_t label : {std::uint32_t{1}, std::uint32_t{2}}) {
		const auto area = diagnostics.area_by_boundary_id.find(label);
		Require(area != diagnostics.area_by_boundary_id.end(), "cap area is absent");
		audit.maximum_area_error_m2 = std::max(audit.maximum_area_error_m2,
			std::abs(area->second-expected_area_m2));
		Require(audit.maximum_area_error_m2 <= 2e-12,
			"cap area roundoff audit failed");
	}
	const auto surface = iga::SurfaceReaders::ReadVtpPath((directory/"surface.vtp").string());
	std::array<std::size_t, 3> triangles{{0, 0, 0}};
	for (const auto& triangle : surface.Triangles()) {
		if (triangle.boundary_id != 1 && triangle.boundary_id != 2) continue;
		const double expected_x = triangle.boundary_id == 1 ? -1.0 : 1.0;
		const double error = std::max({std::abs(triangle.outward_unit_normal[0]-expected_x),
			std::abs(triangle.outward_unit_normal[1]), std::abs(triangle.outward_unit_normal[2])});
		audit.maximum_normal_component_error = std::max(audit.maximum_normal_component_error, error);
		Require(audit.maximum_normal_component_error <= 2e-14,
			"cap outward-normal audit failed");
		++triangles[triangle.boundary_id];
	}
	Require(triangles[1] == 8 && triangles[2] == 8, "cap triangulation audit failed");
	audit.inlet_triangles = triangles[1];
	audit.outlet_triangles = triangles[2];
	return audit;
}

struct Sample {
	double flow_m3_s = 0.0;
	double resistance_pa_s_m3 = 0.0;
	double reynolds = 0.0;
	double controller_absolute_m3_s = 0.0;
	double controller_normalized = 0.0;
	double open_normalized = 0.0;
	double wall_normalized = 0.0;
	PetscInt nonlinear_iterations = 0;
	PetscInt ksp_iterations = 0;
	double nonlinear_residual_norm = 0.0;
	double maximum_true_linear_relative_residual = 0.0;
};

Sample SolveSample(iga::ImmersedStaticFlowRuntime& runtime, double flow_m3_s)
{
	const std::vector<PetscScalar> zero(runtime.Diagnostics().total_dofs, 0.0);
	runtime.SetCommittedState(zero);
	runtime.SetPortControlValue("inlet", -flow_m3_s);
	runtime.SetPortControlValue("outlet", 0.0);
	Require(runtime.SolveTrial(), "quasi-static immersed solve did not converge");
	runtime.Commit();
	const auto& diagnostics = runtime.Diagnostics();
	const auto& inlet = DiagnosticPort(diagnostics, "inlet");
	const auto& outlet = DiagnosticPort(diagnostics, "outlet");
	Require(diagnostics.converged && diagnostics.nonlinear_iterations > 0
		&& diagnostics.nonlinear_iterations <= 30, "immersed Newton gate failed");
	double maximum_true_linear_relative_residual = 0.0;
	for (const auto& step : diagnostics.newton_steps) {
		Require(step.ksp_reason > 0 && std::isfinite(step.linear_relative_residual)
			&& step.linear_relative_residual <= kTrueLinearRelativeTolerance,
			"immersed true linear residual gate failed");
		maximum_true_linear_relative_residual = std::max(maximum_true_linear_relative_residual,
			step.linear_relative_residual);
	}
	Require(std::isfinite(inlet.absolute_flow_residual_m3_s)
		&& inlet.absolute_flow_residual_m3_s <= inlet.flow_tolerance_m3_s
		&& inlet.absolute_flow_residual_m3_s/std::max(flow_m3_s, kTargetFlowM3S)
			<= kControllerRelativeTolerance,
		"immersed controller target gate failed");
	const auto conservation = runtime.ConservationDiagnostics();
	const double throughflow = std::abs(inlet.measurement.outward_flow_m3_s);
	Require(throughflow > 0.0 && std::isfinite(throughflow), "invalid immersed throughflow");
	const double open = std::abs(conservation.open_port_outward_flow_m3_s)/throughflow;
	const double wall = std::abs(conservation.wall_outward_flow_m3_s)/throughflow;
	Require(std::isfinite(open) && std::isfinite(wall) && open <= kOpenBalanceTolerance
		&& wall <= kWallLeakageTolerance, "immersed conservation gate failed");
	const double delta_pressure = inlet.measurement.mean_pressure_pa-outlet.measurement.mean_pressure_pa;
	const double resistance = delta_pressure/throughflow;
	Require(std::isfinite(resistance) && resistance > 0.0, "immersed hydraulic resistance is not positive");
	const double area = inlet.area_m2;
	const double hydraulic_diameter = 2.0*std::sqrt(area/M_PI);
	(void)hydraulic_diameter;
	return {throughflow, resistance, 0.0, inlet.absolute_flow_residual_m3_s,
		inlet.absolute_flow_residual_m3_s/std::max(flow_m3_s, kTargetFlowM3S), open, wall,
		diagnostics.nonlinear_iterations, diagnostics.ksp_iterations, diagnostics.residual_norm,
		maximum_true_linear_relative_residual};
}

double Reynolds(const iga::ImmersedFlowCase& owner, double flow_m3_s)
{
	const auto& inlet = DiagnosticPort(owner.Runtime().Diagnostics(), "inlet");
	const double diameter = 2.0*std::sqrt(inlet.area_m2/M_PI);
	return owner.RuntimeParameters().density*4.0*flow_m3_s
		/(M_PI*owner.RuntimeParameters().dynamic_viscosity*diameter);
}

fs::path DepthTwoCopy(const fs::path& source)
{
	const auto nonce = std::chrono::steady_clock::now().time_since_epoch().count();
	const fs::path copy = fs::temp_directory_path()/
		("tubularflowiga-phase6-depth2-"+std::to_string(nonce));
	fs::create_directory(copy);
	fs::copy_file(source/"simulation_config.json", copy/"simulation_config.json");
	fs::copy_file(source/"surface.vtp", copy/"surface.vtp");
	std::ifstream input(source/"immersed_geometry.json");
	std::ofstream output(copy/"immersed_geometry.json");
	std::string text((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
	const auto pos = text.find("\"max_depth\":3");
	Require(pos != std::string::npos, "depth-3 geometry fixture is malformed");
	text.replace(pos, std::string("\"max_depth\":3").size(), "\"max_depth\":2");
	output << text;
	Require(input.good() || input.eof(), "cannot read depth-3 geometry fixture");
	Require(static_cast<bool>(output), "cannot write depth-2 geometry fixture");
	return copy;
}

void SetInputs(iga::ThreeDImmersedFlowDomain& domain, double time_s, double flow_m3_s)
{
	iga::PortBoundaryData inlet;
	inlet.time_s = time_s;
	inlet.outward_flow_m3_s = -flow_m3_s;
	domain.SetPortInput("inlet", inlet);
	iga::PortBoundaryData outlet;
	outlet.time_s = time_s;
	outlet.mean_pressure_pa = 0.0;
	domain.SetPortInput("outlet", outlet);
}

void RequireEqualState(const std::vector<PetscScalar>& first, const std::vector<PetscScalar>& second)
{
	Require(first.size() == second.size(), "transaction replay state size differs");
	for (std::size_t i = 0; i < first.size(); ++i)
		Require(PetscRealPart(first[i]) == PetscRealPart(second[i]), "transaction retry did not exactly replay committed state");
}

void RequireEqualPortStates(const std::map<std::string, iga::PortState>& first,
	const std::map<std::string, iga::PortState>& second)
{
	Require(first.size() == second.size(), "transaction replay port count differs");
	for (const auto& entry : first) {
		const auto found = second.find(entry.first);
		Require(found != second.end(), "transaction replay port is absent");
		const auto& a = entry.second;
		const auto& b = found->second;
		Require(a.time_s == b.time_s && a.area_m2 == b.area_m2
			&& a.outward_flow_m3_s == b.outward_flow_m3_s
			&& a.mean_pressure_pa == b.mean_pressure_pa
			&& a.mean_normal_traction_pa == b.mean_normal_traction_pa
			&& a.total_pressure_pa == b.total_pressure_pa
			&& a.concentration == b.concentration && a.outward_species_flux == b.outward_species_flux,
			"transaction retry did not exactly replay committed port state");
	}
}

void RequireEqualWork(const iga::ImmersedStaticFlowDiagnostics& first,
	const iga::ImmersedStaticFlowDiagnostics& second)
{
	Require(first.nonlinear_iterations == second.nonlinear_iterations
		&& first.ksp_iterations == second.ksp_iterations && first.ksp_reason == second.ksp_reason
		&& first.newton_steps.size() == second.newton_steps.size(),
		"transaction retry did not exactly replay nonlinear/KSP work");
	for (std::size_t i = 0; i < first.newton_steps.size(); ++i) {
		const auto& a = first.newton_steps[i];
		const auto& b = second.newton_steps[i];
		Require(a.iteration == b.iteration && a.ksp_iterations == b.ksp_iterations
			&& a.ksp_reason == b.ksp_reason,
			"transaction retry did not exactly replay Newton/KSP step work");
	}
}

struct TransactionEvidence {
	std::size_t baseline_domain_committed_steps = 0;
	std::size_t baseline_domain_prepared_count = 0;
	std::size_t baseline_runtime_commit_count = 0;
	std::size_t baseline_runtime_prepare_count = 0;
	std::size_t baseline_runtime_finalize_count = 0;
	std::size_t retry_domain_committed_steps = 0;
	std::size_t retry_domain_prepared_count = 0;
	std::size_t retry_domain_rollback_count = 0;
	std::size_t retry_domain_abort_count = 0;
	std::size_t retry_runtime_commit_count = 0;
	std::size_t retry_runtime_prepare_count = 0;
	std::size_t retry_runtime_finalize_count = 0;
	std::size_t retry_runtime_rollback_count = 0;
	PetscInt baseline_nonlinear_iterations = 0;
	PetscInt baseline_ksp_iterations = 0;
	PetscInt retry_nonlinear_iterations = 0;
	PetscInt retry_ksp_iterations = 0;
};

TransactionEvidence AuditPrecommitFailure(const fs::path& root)
{
	auto baseline_owner = iga::ImmersedFlowCase::Load(root, "immersed", Ports(), 1);
	iga::ThreeDImmersedFlowDomain baseline("immersed", baseline_owner->Runtime(), Ports());
	baseline.BeginStep({0, 0.0, 0.01});
	SetInputs(baseline, 0.01, kTargetFlowM3S);
	baseline.SolveTrial();
	baseline.PrepareCommitStep();
	baseline.FinalizeCommitStep();
	const auto expected = baseline.CommittedBackendState();
	const auto expected_ports = baseline.CommittedPortStates();
	const auto expected_work = baseline_owner->Runtime().Diagnostics();
	const auto baseline_domain = baseline.Diagnostics();
	Require(baseline_domain.committed_steps == 1 && baseline_domain.prepared_count == 1
		&& baseline_domain.rollback_count == 0 && baseline_domain.abort_count == 0,
		"baseline transaction counters are invalid");
	Require(expected_work.commit_count == 1 && expected_work.prepare_count == 1
		&& expected_work.finalize_count == 1 && expected_work.rollback_count == 0,
		"baseline runtime transaction counters are invalid");

	auto retry_owner = iga::ImmersedFlowCase::Load(root, "immersed", Ports(), 1);
	iga::ThreeDImmersedFlowDomain retry("immersed", retry_owner->Runtime(), Ports());
	retry.BeginStep({0, 0.0, 0.01});
	SetInputs(retry, 0.01, kTargetFlowM3S);
	retry.SolveTrial();
	retry.PrepareCommitStep();
	retry.FinalizeCommitStep();
	RequireEqualState(expected, retry.CommittedBackendState());
	RequireEqualPortStates(expected_ports, retry.CommittedPortStates());
	RequireEqualWork(expected_work, retry_owner->Runtime().Diagnostics());
	// Reject a later transaction so these snapshots contain the publicly
	// committed inlet and outlet values, rather than an initially empty map.
	const auto before_state = retry.CommittedBackendState();
	const auto before_ports = retry.CommittedPortStates();
	Require(before_ports.size() == Ports().size(),
		"retry did not publish every configured committed port state");
	const auto before = retry.Diagnostics();
	retry.BeginStep({1, 0.01, 0.01});
	SetInputs(retry, 0.02, kTargetFlowM3S);
	retry.SolveTrial();
	retry_owner->Runtime().FailNextPrepareForTesting();
	bool failed = false;
	try { retry.PrepareCommitStep(); }
	catch (const std::runtime_error&) { failed = true; }
	Require(failed, "injected precommit failure was not observed");
	const auto after_failure = retry.Diagnostics();
	RequireEqualState(before_state, retry.CommittedBackendState());
	RequireEqualPortStates(before_ports, retry.CommittedPortStates());
	Require(after_failure.committed_time_s == before.committed_time_s
		&& after_failure.committed_step_index == before.committed_step_index
		&& after_failure.committed_steps == before.committed_steps
		&& after_failure.prepared_count == before.prepared_count,
		"precommit failure mutated immersed committed time or counters");
	const auto runtime_after_failure = retry_owner->Runtime().Diagnostics();
	Require(runtime_after_failure.commit_count == expected_work.commit_count
		&& runtime_after_failure.prepare_count == expected_work.prepare_count
		&& runtime_after_failure.finalize_count == expected_work.finalize_count
		&& runtime_after_failure.rollback_count == expected_work.rollback_count,
		"precommit failure mutated runtime commit counters");
	retry.RollbackTrial();
	retry.AbortStep();
	const auto retry_work = retry_owner->Runtime().Diagnostics();
	const auto retry_domain = retry.Diagnostics();
	Require(retry_domain.committed_steps == 1 && retry_domain.prepared_count == 1
		&& retry_domain.rollback_count == 1 && retry_domain.abort_count == 1,
		"retry transaction counters are invalid");
	Require(retry_work.commit_count == 1 && retry_work.prepare_count == 1
		&& retry_work.finalize_count == 1 && retry_work.rollback_count == 1,
		"retry runtime transaction counters are invalid");
	return {baseline_domain.committed_steps, baseline_domain.prepared_count,
		expected_work.commit_count, expected_work.prepare_count, expected_work.finalize_count,
		retry_domain.committed_steps, retry_domain.prepared_count, retry_domain.rollback_count,
		retry_domain.abort_count, retry_work.commit_count, retry_work.prepare_count,
		retry_work.finalize_count, retry_work.rollback_count, expected_work.nonlinear_iterations,
		expected_work.ksp_iterations, retry_work.nonlinear_iterations, retry_work.ksp_iterations};
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int status = 0;
	fs::path depth_two;
	try {
		const fs::path depth_three = "../../examples/vascular_flow/immersed_aneurysm_chain/immersed";
		auto three = iga::ImmersedFlowCase::Load(depth_three, "immersed", Ports(), 1);
		const CapAudit cap_audit = AuditCaps(depth_three, three->SurfaceDiagnostics());
		std::vector<Sample> samples;
		for (const double flow : {0.5*kTargetFlowM3S, 0.75*kTargetFlowM3S, kTargetFlowM3S}) {
			auto sample = SolveSample(three->Runtime(), flow);
			sample.reynolds = Reynolds(*three, flow);
			Require(sample.reynolds < 1.0, "quasi-static aneurysm Reynolds gate failed");
			samples.push_back(sample);
		}
		depth_two = DepthTwoCopy(depth_three);
		auto two = iga::ImmersedFlowCase::Load(depth_two, "immersed", Ports(), 1);
		const Sample depth_two_sample = SolveSample(two->Runtime(), kTargetFlowM3S);
		const double resistance = samples.back().resistance_pa_s_m3;
		const double refined_difference = std::abs(resistance-depth_two_sample.resistance_pa_s_m3)
			/std::max(resistance, depth_two_sample.resistance_pa_s_m3);
		// Cell-centre lubrication audit of the explicit [0,.25,.5,.75,1] m
		// radius schedule.  The 0.5 m bulge cell has r=.3 m; the two 0.25 m
		// end cells have r=.25 m, yielding 8 mu/pi sum(dx/r^4).
		const double lubrication = 8.0*three->RuntimeParameters().dynamic_viscosity/M_PI
			*(0.25/std::pow(0.25, 4)+0.5/std::pow(0.3, 4)+0.25/std::pow(0.25, 4));
		const double lubrication_difference = std::abs(resistance-lubrication)/lubrication;
		Require(lubrication_difference <= 0.25, "hydraulic resistance differs from lubrication estimate by more than 25%");
		Require(refined_difference <= 0.10, "depth-2/depth-3 hydraulic resistance differs by more than 10%");
		for (std::size_t i = 0; i < samples.size(); ++i) {
			const auto& sample = samples[i];
			std::cout << "{\"phase6_immersed_sample\":{\"flow_m3_s\":" << sample.flow_m3_s
				<< ",\"depth\":3,\"sample_index\":" << i
				<< ",\"resistance_pa_s_m3\":" << sample.resistance_pa_s_m3
				<< ",\"reynolds\":" << sample.reynolds
				<< ",\"controller_absolute_m3_s\":" << sample.controller_absolute_m3_s
				<< ",\"controller_normalized\":" << sample.controller_normalized
				<< ",\"open_port_normalized_balance\":" << sample.open_normalized
				<< ",\"wall_leakage_normalized\":" << sample.wall_normalized
				<< ",\"nonlinear_iterations\":" << sample.nonlinear_iterations
				<< ",\"ksp_iterations\":" << sample.ksp_iterations
				<< ",\"nonlinear_residual_norm\":" << sample.nonlinear_residual_norm
				<< ",\"maximum_true_linear_relative_residual\":" << sample.maximum_true_linear_relative_residual
				<< ",\"nonlinear_iteration_gate\":true,\"controller_gate\":true"
				<< ",\"open_balance_gate\":true,\"wall_leakage_gate\":true,\"true_linear_gate\":true}}\n";
		}
		std::cout << "{\"phase6_immersed_sample\":{\"flow_m3_s\":" << depth_two_sample.flow_m3_s
			<< ",\"depth\":2,\"sample_index\":0"
			<< ",\"resistance_pa_s_m3\":" << depth_two_sample.resistance_pa_s_m3
			<< ",\"reynolds\":" << Reynolds(*two, depth_two_sample.flow_m3_s)
			<< ",\"controller_absolute_m3_s\":" << depth_two_sample.controller_absolute_m3_s
			<< ",\"controller_normalized\":" << depth_two_sample.controller_normalized
			<< ",\"open_port_normalized_balance\":" << depth_two_sample.open_normalized
			<< ",\"wall_leakage_normalized\":" << depth_two_sample.wall_normalized
			<< ",\"nonlinear_iterations\":" << depth_two_sample.nonlinear_iterations
			<< ",\"ksp_iterations\":" << depth_two_sample.ksp_iterations
			<< ",\"nonlinear_residual_norm\":" << depth_two_sample.nonlinear_residual_norm
			<< ",\"maximum_true_linear_relative_residual\":" << depth_two_sample.maximum_true_linear_relative_residual
			<< ",\"nonlinear_iteration_gate\":true,\"controller_gate\":true"
			<< ",\"open_balance_gate\":true,\"wall_leakage_gate\":true,\"true_linear_gate\":true}}\n";
		std::cout << "{\"phase6_caps\":{\"expected_area_m2\":0.17677669529663687"
			<< ",\"maximum_area_error_m2\":" << cap_audit.maximum_area_error_m2
			<< ",\"area_error_gate_m2\":2e-12"
			<< ",\"maximum_normal_component_error\":" << cap_audit.maximum_normal_component_error
			<< ",\"normal_component_error_gate\":2e-14"
			<< ",\"inlet_triangles\":" << cap_audit.inlet_triangles
			<< ",\"outlet_triangles\":" << cap_audit.outlet_triangles
			<< ",\"cap_area_gate\":true,\"cap_normal_gate\":true}}\n";
		std::cout << "{\"phase6_resistance\":{\"depth2_pa_s_m3\":" << depth_two_sample.resistance_pa_s_m3
			<< ",\"depth3_pa_s_m3\":" << resistance << ",\"lubrication_pa_s_m3\":" << lubrication
			<< ",\"lubrication_relative_difference\":" << lubrication_difference
			<< ",\"lubrication_difference_gate\":0.25,\"depth_relative_difference\":" << refined_difference
			<< ",\"depth_difference_gate\":0.1}}\n";
		const TransactionEvidence transaction = AuditPrecommitFailure(depth_two);
		std::cout << "{\"phase6_precommit_failure\":{\"completion_published\":false"
			<< ",\"baseline_domain_committed_steps\":" << transaction.baseline_domain_committed_steps
			<< ",\"baseline_domain_prepared_count\":" << transaction.baseline_domain_prepared_count
			<< ",\"baseline_runtime_commit_count\":" << transaction.baseline_runtime_commit_count
			<< ",\"baseline_runtime_prepare_count\":" << transaction.baseline_runtime_prepare_count
			<< ",\"baseline_runtime_finalize_count\":" << transaction.baseline_runtime_finalize_count
			<< ",\"retry_domain_committed_steps\":" << transaction.retry_domain_committed_steps
			<< ",\"retry_domain_prepared_count\":" << transaction.retry_domain_prepared_count
			<< ",\"retry_domain_rollback_count\":" << transaction.retry_domain_rollback_count
			<< ",\"retry_domain_abort_count\":" << transaction.retry_domain_abort_count
			<< ",\"retry_runtime_commit_count\":" << transaction.retry_runtime_commit_count
			<< ",\"retry_runtime_prepare_count\":" << transaction.retry_runtime_prepare_count
			<< ",\"retry_runtime_finalize_count\":" << transaction.retry_runtime_finalize_count
			<< ",\"retry_runtime_rollback_count\":" << transaction.retry_runtime_rollback_count
			<< ",\"baseline_nonlinear_iterations\":" << transaction.baseline_nonlinear_iterations
			<< ",\"baseline_ksp_iterations\":" << transaction.baseline_ksp_iterations
			<< ",\"retry_nonlinear_iterations\":" << transaction.retry_nonlinear_iterations
			<< ",\"retry_ksp_iterations\":" << transaction.retry_ksp_iterations
			<< ",\"exact_committed_state_replay\":true,\"exact_committed_port_replay\":true"
			<< ",\"exact_nonlinear_ksp_work_replay\":true}}\n";
	} catch (const std::exception& error) {
		std::cerr << "phase6 aneurysm runtime closure failure: " << error.what() << '\n';
		status = 1;
	}
	std::error_code ignored;
	if (!depth_two.empty()) fs::remove_all(depth_two, ignored);
	PetscFinalize();
	return status;
}
