#include "CouplingFixture.hpp"
#include "OneDRuntime.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr double kPeriodS = 0.16;
// Advance four periods, compare the last two, and measure the final one.  The
// first periods damp the start from the steady mean-flow profile without
// making this modest-mesh opt-in test prohibitively expensive.
constexpr double kDurationS = 4.0*kPeriodS;
constexpr double kMeanFlowM3S = 1.0e-3;
// Retain low Reynolds number while making the unsteady term visible over one
// 0.16 s pulse; a nearly zero density would reduce this to a sampled steady
// resistance check rather than a temporal verification.
constexpr double kDensityKgM3 = 1.0;
constexpr double kDynamicViscosityPaS = 1.0;
constexpr int kTransverseElements = 2;
constexpr int kAxialElements = 2;
constexpr int kOneDSubstepsPerMacro = 2;
constexpr double kTightPressureTolerance = 1.0e-10;

struct TemporaryDirectory {
	fs::path path;
	~TemporaryDirectory()
	{
		if (std::getenv("TUBULARFLOWIGA_KEEP_TEST_OUTPUT")) return;
		std::error_code error;
		fs::remove_all(path, error);
	}
};

struct CsvTable {
	std::vector<std::map<std::string, double>> rows;
};

struct ReferenceSample {
	double time_s = 0.0;
	double root_pressure_pa = 0.0;
	double terminal_pressure_pa = 0.0;
	double root_outward_flow_m3_s = 0.0;
	double terminal_outward_flow_m3_s = 0.0;
};

struct RunResult {
	double macro_dt_s = 0.0;
	double pressure_tolerance = 0.0;
	CsvTable history;
	CsvTable iterations;
	std::vector<ReferenceSample> reference;
};

struct WaveformMetrics {
	double linf_relative = 0.0;
	double l2_relative = 0.0;
	double amplitude = 0.0;
	double phase_rad = 0.0;
};

std::string ReadText(const fs::path& path)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot read "+path.string());
	std::ostringstream text;
	text << input.rdbuf();
	return text.str();
}

std::string Number(double value)
{
	std::ostringstream output;
	output << std::setprecision(17) << value;
	return output.str();
}

std::string Quote(const fs::path& path)
{
	return "'"+path.string()+"'";
}

CsvTable ReadCsv(const fs::path& path)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("missing CSV "+path.string());
	std::string header;
	if (!std::getline(input, header) || header.empty()) throw std::runtime_error("empty CSV "+path.string());
	std::vector<std::string> names;
	std::istringstream headings(header);
	for (std::string name; std::getline(headings, name, ',');) names.push_back(name);
	CsvTable table;
	for (std::string line; std::getline(input, line);) {
		if (line.empty()) continue;
		std::istringstream values(line);
		std::map<std::string, double> row;
		for (const auto& name : names) {
			std::string value;
			if (!std::getline(values, value, ',')) throw std::runtime_error("truncated CSV row in "+path.string());
			const double parsed = std::stod(value);
			if (!std::isfinite(parsed)) throw std::runtime_error("nonfinite CSV value in "+path.string());
			row.emplace(name, parsed);
		}
		std::string extra;
		if (std::getline(values, extra, ',')) throw std::runtime_error("overlong CSV row in "+path.string());
		table.rows.push_back(std::move(row));
	}
	return table;
}

double Value(const std::map<std::string, double>& row, const std::string& name)
{
	const auto found = row.find(name);
	if (found == row.end()) throw std::runtime_error("CSV is missing column "+name);
	return found->second;
}

long long IntegerValue(const std::map<std::string, double>& row, const std::string& name)
{
	const double value = Value(row, name);
	const auto integer = static_cast<long long>(std::llround(value));
	if (std::abs(value-integer) > 1.0e-12) throw std::runtime_error("expected integer CSV column "+name);
	return integer;
}

void WriteHydraulicEquivalentNetwork(const fs::path& directory, double radius_m)
{
	std::ofstream output(directory/"tree.swc", std::ios::trunc);
	if (!output) throw std::runtime_error("cannot create all-1D hydraulic-equivalent network");
	const double middle_end = 1.0+iga::test::kEquivalentSquareDuctHydraulicLengthM;
	output << std::setprecision(17)
		<< "1 2 0 0 0 " << radius_m << " -1\n"
		<< "2 2 1 0 0 " << radius_m << " 1\n"
		<< "3 2 " << middle_end << " 0 0 " << radius_m << " 2\n"
		<< "4 2 " << middle_end+1.0 << " 0 0 " << radius_m << " 3\n";
}

std::vector<ReferenceSample> AdvanceAllOneDReference(const fs::path& directory, double dt_s, int steps)
{
	const auto configuration = iga::ParseOneDConfiguration(ReadText(directory/"simulation_config.json"));
	if (configuration.flow_systems.size() != 1) throw std::runtime_error("all-1D reference requires one flow system");
	const auto& flow = configuration.flow_systems.front();
	const auto network = iga::ReadOneDNetwork(directory/"tree.swc", 1.0, 1, flow.dynamic_viscosity);
	iga::OneDFlowRuntime runtime(configuration, flow, network, iga::ResolveOneDInlet(configuration), directory);
	runtime.InitializeOpenLoop(kMeanFlowM3S);
	std::vector<ReferenceSample> history;
	for (int step = 1; step <= steps; ++step) {
		runtime.BeginStep(runtime.FlowState().physical_time, dt_s);
		runtime.SetConfiguredOpenLoopInlet();
		runtime.SolveTrial();
		const auto root = runtime.GetPortState("root");
		const auto terminal = runtime.GetPortState("outlet:4");
		if (!root.mean_pressure_pa || !terminal.mean_pressure_pa || !root.outward_flow_m3_s || !terminal.outward_flow_m3_s)
			throw std::runtime_error("all-1D reference does not provide required port data");
		history.push_back({runtime.FlowState().physical_time, *root.mean_pressure_pa, *terminal.mean_pressure_pa,
			*root.outward_flow_m3_s, *terminal.outward_flow_m3_s});
		runtime.CommitStep();
	}
	return history;
}

std::vector<double> Extract(const CsvTable& table, const std::string& name)
{
	std::vector<double> values;
	values.reserve(table.rows.size());
	for (const auto& row : table.rows) values.push_back(Value(row, name));
	return values;
}

std::vector<double> TailPeriod(const std::vector<double>& values, int samples_per_period)
{
	if (samples_per_period < 4 || values.size() < static_cast<std::size_t>(samples_per_period))
		throw std::runtime_error("temporal metric needs at least four samples in the final period");
	return std::vector<double>(values.end()-samples_per_period, values.end());
}

std::vector<double> PreviousPeriod(const std::vector<double>& values, int samples_per_period)
{
	if (samples_per_period < 4 || values.size() < static_cast<std::size_t>(2*samples_per_period))
		throw std::runtime_error("cycle audit requires at least two complete periods");
	return std::vector<double>(values.end()-2*samples_per_period, values.end()-samples_per_period);
}

WaveformMetrics MeasureWaveform(const std::vector<double>& values, double dt_s, double scale)
{
	if (!(dt_s > 0.0) || !(scale > 0.0) || values.size() < 4) throw std::runtime_error("invalid waveform metric input");
	double mean = 0.0;
	for (double value : values) mean += value;
	mean /= static_cast<double>(values.size());
	std::complex<double> harmonic = 0.0;
	for (std::size_t index = 0; index < values.size(); ++index) {
		// History rows are macro-step endpoints: the final-period tail begins at
		// period_start+dt, not at period_start.  Use that physical endpoint phase
		// so grids with different dt share the same period origin.
		const double angle = -2.0*iga::test::kPi*static_cast<double>(index+1)
			/static_cast<double>(values.size());
		harmonic += (values[index]-mean)*std::complex<double>(std::cos(angle), std::sin(angle));
	}
	harmonic *= 2.0/static_cast<double>(values.size());
	WaveformMetrics metric;
	metric.amplitude = std::abs(harmonic);
	metric.phase_rad = std::atan2(harmonic.imag(), harmonic.real());
	metric.linf_relative = 0.0;
	for (double value : values) metric.linf_relative = std::max(metric.linf_relative, std::abs(value-mean)/scale);
	double squared = 0.0;
	for (double value : values) squared += std::pow((value-mean)/scale, 2);
	metric.l2_relative = std::sqrt(squared/static_cast<double>(values.size()));
	return metric;
}

WaveformMetrics CompareWaveforms(const std::vector<double>& first, const std::vector<double>& second,
	double dt_s, double scale)
{
	if (first.size() != second.size()) throw std::runtime_error("unaligned waveform comparison");
	std::vector<double> difference(first.size());
	for (std::size_t index = 0; index < first.size(); ++index) difference[index] = first[index]-second[index];
	WaveformMetrics metric = MeasureWaveform(difference, dt_s, scale);
	double squared = 0.0;
	metric.linf_relative = 0.0;
	for (double value : difference) {
		metric.linf_relative = std::max(metric.linf_relative, std::abs(value)/scale);
		squared += std::pow(value/scale, 2);
	}
	metric.l2_relative = std::sqrt(squared/static_cast<double>(difference.size()));
	return metric;
}

void CheckMetricUtilities()
{
	constexpr double amplitude = 0.37;
	constexpr double input_phase = 0.41;
	const double expected_phase = input_phase-0.5*iga::test::kPi;
	for (int samples : {16, 32, 64}) {
		std::vector<double> endpoint_values;
		for (int index = 0; index < samples; ++index) {
			const double endpoint_angle = 2.0*iga::test::kPi*static_cast<double>(index+1)
				/static_cast<double>(samples);
			endpoint_values.push_back(1.25+amplitude*std::sin(endpoint_angle+input_phase));
		}
		const auto metric = MeasureWaveform(endpoint_values, kPeriodS/samples, 1.0);
		if (std::abs(metric.amplitude-amplitude) > 64.0*std::numeric_limits<double>::epsilon()
			|| std::abs(iga::test::HarmonicPhaseDifferenceRadians(metric.phase_rad, expected_phase))
				> 64.0*std::numeric_limits<double>::epsilon())
			throw std::runtime_error("endpoint harmonic metric has a grid-dependent phase origin");
	}
	std::vector<double> cycles;
	for (int cycle = 0; cycle < 4; ++cycle)
		for (int sample = 0; sample < 16; ++sample) cycles.push_back(100.0*cycle+sample);
	const auto previous = PreviousPeriod(cycles, 16);
	const auto tail = TailPeriod(cycles, 16);
	if (previous.front() != 200.0 || previous.back() != 215.0
		|| tail.front() != 300.0 || tail.back() != 315.0)
		throw std::runtime_error("cycle audit did not select the final two periods");
}

std::vector<double> RestrictToCoarseTimes(const std::vector<double>& fine, int ratio)
{
	if (ratio < 1 || fine.size()%static_cast<std::size_t>(ratio) != 0)
		throw std::runtime_error("invalid nested temporal grid ratio");
	std::vector<double> restricted;
	for (std::size_t index = static_cast<std::size_t>(ratio-1); index < fine.size(); index += static_cast<std::size_t>(ratio))
		restricted.push_back(fine[index]);
	return restricted;
}

void CheckStrongArtifacts(const RunResult& result, int macro_steps, const fs::path& output_directory)
{
	if (result.history.rows.size() != static_cast<std::size_t>(macro_steps))
		throw std::runtime_error("strong history does not contain every macro endpoint");
	if (result.iterations.rows.empty()) throw std::runtime_error("strong iteration history is empty");
	const std::string manifest = ReadText(output_directory/"strong_coupling_manifest.json");
	if (manifest.find("\"scheme\": \"strong_aitken\"") == std::string::npos
		|| manifest.find("\"N\": 2") == std::string::npos
		|| manifest.find("endpoint residuals") == std::string::npos)
		throw std::runtime_error("strong manifest does not record Aitken/subcycling semantics");
	std::map<int, long long> iteration_rows;
	std::map<int, long long> converged_rows;
	for (const auto& row : result.iterations.rows) {
		const int step = static_cast<int>(IntegerValue(row, "physical_step"));
		if (step < 1 || step > macro_steps
			|| std::abs(Value(row, "time_s")-step*result.macro_dt_s)
				> 64.0*std::numeric_limits<double>::epsilon()*std::max(1.0, step*result.macro_dt_s))
			throw std::runtime_error("strong iteration history is not aligned to its macro step");
		++iteration_rows[step];
		converged_rows[step] += IntegerValue(row, "converged");
		if (IntegerValue(row, "upstream_1d_attempt_configured_substeps") != kOneDSubstepsPerMacro
			|| IntegerValue(row, "downstream_1d_attempt_configured_substeps") != kOneDSubstepsPerMacro)
			throw std::runtime_error("one-dimensional subcycling attempt accounting is incorrect");
		if (std::abs(Value(row, "normalized_upstream_three_d_flow_residual")) > 1.0e-10
			|| std::abs(Value(row, "normalized_three_d_downstream_flow_residual")) > 1.0e-10)
			throw std::runtime_error("strong transfer-flow residual exceeds configured gate");
		const double proposed_relaxation = Value(row, "relaxation_factor_for_next_guess");
		if (!(proposed_relaxation >= 0.0 && proposed_relaxation <= 1.0)
			|| IntegerValue(row, "aitken_status_code") < 0
			|| IntegerValue(row, "aitken_status_code") > 5)
			throw std::runtime_error("strong iteration Aitken diagnostics are invalid");
	}
	for (int step = 1; step <= macro_steps; ++step) {
		const auto& row = result.history.rows.at(static_cast<std::size_t>(step-1));
		const double expected_time = step*result.macro_dt_s;
		if (std::abs(Value(row, "time_s")-expected_time) > 64.0*std::numeric_limits<double>::epsilon()*std::max(1.0, expected_time))
			throw std::runtime_error("coupling history time is not aligned to the 3D macro grid");
		const long long sweeps = IntegerValue(row, "iteration_count");
		if (sweeps != iteration_rows[step] || sweeps < 1 || converged_rows[step] != 1)
			throw std::runtime_error("iteration history disagrees with step summary");
		for (const char* prefix : {"upstream", "downstream"}) {
			const std::string base(prefix);
			const long long accepted = IntegerValue(row, base+"_accepted_configured_substeps");
			const long long all = IntegerValue(row, base+"_all_configured_substeps");
			const long long rejected = IntegerValue(row, base+"_rejected_configured_substeps");
			if (accepted != kOneDSubstepsPerMacro || all != sweeps*kOneDSubstepsPerMacro || rejected != all-accepted)
				throw std::runtime_error("one-dimensional configured work totals are inconsistent");
		}
		if (Value(row, "final_max_normalized_pressure_residual") > result.pressure_tolerance*(1.0+64.0*std::numeric_limits<double>::epsilon()))
			throw std::runtime_error("committed strong step exceeds its requested pressure tolerance");
		if (std::abs(Value(row, "three_d_mass_imbalance_m3_s"))/kMeanFlowM3S > 1.0e-6
			|| std::abs(Value(row, "three_d_wall_outward_flow_m3_s"))/kMeanFlowM3S > 1.0e-10)
			throw std::runtime_error("3D mass or wall-flow diagnostic is outside its established gate");
		const double relaxation = Value(row, "relaxation_factor");
		if (!(relaxation >= 0.0 && relaxation <= 1.0) || !std::isfinite(relaxation))
			throw std::runtime_error("Aitken relaxation factor is invalid");
	}
}

RunResult RunPulsatileCase(const fs::path& root, double macro_dt_s, double pressure_tolerance)
{
	const int macro_steps = static_cast<int>(std::llround(kDurationS/macro_dt_s));
	if (std::abs(macro_steps*macro_dt_s-kDurationS) > 1.0e-14) throw std::runtime_error("temporal duration is not on macro grid");
	const double one_d_dt_s = macro_dt_s/kOneDSubstepsPerMacro;
	const int one_d_steps = macro_steps*kOneDSubstepsPerMacro;
	const auto case_directory = root/("dt-"+Number(macro_dt_s)+"-tol-"+Number(pressure_tolerance));
	const auto output = case_directory/"output";
	fs::create_directories(case_directory/"three_d");
	fs::create_directories(case_directory/"upstream");
	fs::create_directories(case_directory/"downstream");
	fs::create_directories(case_directory/"all_one_d");
	iga::test::WriteC2SquareDuctDatabase(case_directory/"case.ntiga", kTransverseElements, 1,
		kAxialElements, 1.0);
	iga::test::WriteC2SquareDuctThreeDCase(case_directory/"three_d", kTransverseElements, macro_dt_s,
		macro_steps, kDensityKgM3, kDynamicViscosityPaS, kMeanFlowM3S, kAxialElements, 1.0);
	const double radius_m = std::sqrt(1.0/iga::test::kPi);
	iga::test::WriteRigidOneDStraightCase(case_directory/"upstream", 1.0, radius_m, one_d_dt_s,
		one_d_steps, kDensityKgM3, kDynamicViscosityPaS, kMeanFlowM3S, "sinusoid", kPeriodS);
	iga::test::WriteRigidOneDStraightCase(case_directory/"downstream", 1.0, radius_m, one_d_dt_s,
		one_d_steps, kDensityKgM3, kDynamicViscosityPaS, kMeanFlowM3S);
	iga::test::WriteRigidOneDStraightCase(case_directory/"all_one_d", 1.0, radius_m, one_d_dt_s,
		one_d_steps, kDensityKgM3, kDynamicViscosityPaS, kMeanFlowM3S, "sinusoid", kPeriodS, 4);
	WriteHydraulicEquivalentNetwork(case_directory/"all_one_d", radius_m);
	const double pressure_reference = (16.0*iga::test::kPi+iga::test::kSquareDuctResistanceCoefficient)*kMeanFlowM3S;
	const auto log = case_directory/"driver.log";
	const std::string command = "mpiexec -np 1 ./iga_1d_3d_explicit "+Quote(case_directory/"case.ntiga")+" "
		+Quote(case_directory/"three_d")+" "+Quote(case_directory/"upstream")+" "+Quote(case_directory/"downstream")
		+" --upstream-terminal-node 2 --output-dir "+Quote(output)
		+" --coupling-mode strong-aitken --strong-pressure-reference-pa "+Number(pressure_reference)
		+" --strong-pressure-relative-tol "+Number(pressure_tolerance)
		+" --strong-flow-relative-tol 1e-10 --strong-max-iterations 30 --three-d-max-newton 60"
		+" -ksp_type preonly -pc_type lu > "+Quote(log)+" 2>&1";
	if (std::system(command.c_str()) != 0)
		throw std::runtime_error("pulsatile strong-Aitken driver failed for dt="+Number(macro_dt_s)+"\n"+ReadText(log));
	RunResult result;
	result.macro_dt_s = macro_dt_s;
	result.pressure_tolerance = pressure_tolerance;
	result.history = ReadCsv(output/"strong_coupling_history.csv");
	result.iterations = ReadCsv(output/"strong_coupling_iterations.csv");
	result.reference = AdvanceAllOneDReference(case_directory/"all_one_d", one_d_dt_s, one_d_steps);
	CheckStrongArtifacts(result, macro_steps, output);
	return result;
}

std::vector<double> ReferenceAtMacroEndpoints(const std::vector<ReferenceSample>& samples,
	int macro_steps, double ReferenceSample::*member)
{
	if (samples.size() != static_cast<std::size_t>(macro_steps*kOneDSubstepsPerMacro))
		throw std::runtime_error("all-1D reference has an unexpected substep count");
	std::vector<double> values;
	for (int step = 1; step <= macro_steps; ++step)
		values.push_back(samples.at(static_cast<std::size_t>(step*kOneDSubstepsPerMacro-1)).*member);
	return values;
}

void CheckReferenceAgreement(const RunResult& result)
{
	const int macro_steps = static_cast<int>(result.history.rows.size());
	const auto coupled_root_q = Extract(result.history, "upstream_root_outward_flow_m3_s");
	const auto coupled_terminal_q = Extract(result.history, "downstream_terminal_outward_flow_m3_s");
	const auto reference_root_q = ReferenceAtMacroEndpoints(result.reference, macro_steps, &ReferenceSample::root_outward_flow_m3_s);
	const auto reference_terminal_q = ReferenceAtMacroEndpoints(result.reference, macro_steps, &ReferenceSample::terminal_outward_flow_m3_s);
	const auto root_q_difference = CompareWaveforms(coupled_root_q, reference_root_q, result.macro_dt_s, kMeanFlowM3S);
	const auto terminal_q_difference = CompareWaveforms(coupled_terminal_q, reference_terminal_q, result.macro_dt_s, kMeanFlowM3S);
	if (root_q_difference.linf_relative > 1.0e-10 || terminal_q_difference.linf_relative > 1.0e-10)
		throw std::runtime_error("pulsatile external flow differs from independently advanced all-1D reference");
	const auto coupled_drop = Extract(result.history, "external_pressure_drop_pa");
	const auto reference_root = ReferenceAtMacroEndpoints(result.reference, macro_steps, &ReferenceSample::root_pressure_pa);
	const auto reference_terminal = ReferenceAtMacroEndpoints(result.reference, macro_steps, &ReferenceSample::terminal_pressure_pa);
	std::vector<double> reference_drop(reference_root.size());
	for (std::size_t index = 0; index < reference_drop.size(); ++index) reference_drop[index] = reference_root[index]-reference_terminal[index];
	const double pressure_scale = (16.0*iga::test::kPi+iga::test::kSquareDuctResistanceCoefficient)*kMeanFlowM3S;
	const int samples_per_period = static_cast<int>(std::llround(kPeriodS/result.macro_dt_s));
	const auto coupled_period = TailPeriod(coupled_drop, samples_per_period);
	const auto reference_period = TailPeriod(reference_drop, samples_per_period);
	const auto pressure_difference = CompareWaveforms(coupled_period, reference_period, result.macro_dt_s, pressure_scale);
	const auto coupled_metric = MeasureWaveform(coupled_period, result.macro_dt_s, pressure_scale);
	const auto reference_metric = MeasureWaveform(reference_period, result.macro_dt_s, pressure_scale);
	std::cout << std::setprecision(17) << "temporal_reference dt_s=" << result.macro_dt_s
		<< " q_linf=" << root_q_difference.linf_relative
		<< " pressure_drop_l2=" << pressure_difference.l2_relative
		<< " coupled_drop_amplitude_pa=" << coupled_metric.amplitude
		<< " reference_drop_amplitude_pa=" << reference_metric.amplitude
		<< " phase_difference_rad=" << iga::test::HarmonicPhaseDifferenceRadians(coupled_metric.phase_rad, reference_metric.phase_rad) << '\n';
}

void CheckPeriodicSettling(const RunResult& result)
{
	const int samples_per_period = static_cast<int>(std::llround(kPeriodS/result.macro_dt_s));
	const double pressure_scale = (16.0*iga::test::kPi+iga::test::kSquareDuctResistanceCoefficient)*kMeanFlowM3S;
	double maximum_pressure_l2 = 0.0, maximum_pressure_linf = 0.0;
	for (const char* field : {"external_pressure_drop_pa", "upstream_root_pressure_pa",
		"three_d_inlet_pressure_pa", "three_d_outlet_pressure_pa", "downstream_terminal_pressure_pa"}) {
		const auto values = Extract(result.history, field);
		const auto difference = CompareWaveforms(PreviousPeriod(values, samples_per_period),
			TailPeriod(values, samples_per_period), result.macro_dt_s, pressure_scale);
		maximum_pressure_l2 = std::max(maximum_pressure_l2, difference.l2_relative);
		maximum_pressure_linf = std::max(maximum_pressure_linf, difference.linf_relative);
	}
	double maximum_flow_l2 = 0.0, maximum_flow_linf = 0.0;
	for (const char* field : {"upstream_terminal_outward_flow_m3_s", "three_d_inlet_outward_flow_m3_s",
		"three_d_outlet_outward_flow_m3_s", "downstream_root_outward_flow_m3_s"}) {
		const auto values = Extract(result.history, field);
		const auto difference = CompareWaveforms(PreviousPeriod(values, samples_per_period),
			TailPeriod(values, samples_per_period), result.macro_dt_s, kMeanFlowM3S);
		maximum_flow_l2 = std::max(maximum_flow_l2, difference.l2_relative);
		maximum_flow_linf = std::max(maximum_flow_linf, difference.linf_relative);
	}
	std::cout << std::setprecision(17) << "periodic_settling dt_s=" << result.macro_dt_s
		<< " pressure_l2=" << maximum_pressure_l2 << " pressure_linf=" << maximum_pressure_linf
		<< " flow_l2=" << maximum_flow_l2 << " flow_linf=" << maximum_flow_linf << '\n';
	if (maximum_pressure_linf > 1.0e-10 || maximum_flow_linf > 1.0e-10)
		throw std::runtime_error("pulsatile case did not reach a repeatable final cycle");
}

void CheckTemporalSelfConvergence(const RunResult& coarse, const RunResult& medium, const RunResult& fine)
{
	const int coarse_to_medium = static_cast<int>(std::llround(coarse.macro_dt_s/medium.macro_dt_s));
	const int medium_to_fine = static_cast<int>(std::llround(medium.macro_dt_s/fine.macro_dt_s));
	const int coarse_samples = static_cast<int>(std::llround(kPeriodS/coarse.macro_dt_s));
	const int medium_samples = static_cast<int>(std::llround(kPeriodS/medium.macro_dt_s));
	const int fine_samples = static_cast<int>(std::llround(kPeriodS/fine.macro_dt_s));
	const double pressure_scale = (16.0*iga::test::kPi+iga::test::kSquareDuctResistanceCoefficient)*kMeanFlowM3S;
	double maximum_pressure_l2_coarse_medium = 0.0, maximum_pressure_l2_medium_fine = 0.0;
	double maximum_pressure_linf_coarse_medium = 0.0, maximum_pressure_linf_medium_fine = 0.0;
	for (const char* field : {"external_pressure_drop_pa", "upstream_root_pressure_pa",
		"three_d_inlet_pressure_pa", "three_d_outlet_pressure_pa", "downstream_terminal_pressure_pa"}) {
		const auto coarse_values = TailPeriod(Extract(coarse.history, field), coarse_samples);
		const auto medium_values = TailPeriod(Extract(medium.history, field), medium_samples);
		const auto fine_values = TailPeriod(Extract(fine.history, field), fine_samples);
		const auto coarse_medium = CompareWaveforms(coarse_values,
			RestrictToCoarseTimes(medium_values, coarse_to_medium), coarse.macro_dt_s, pressure_scale);
		const auto medium_fine = CompareWaveforms(medium_values,
			RestrictToCoarseTimes(fine_values, medium_to_fine), medium.macro_dt_s, pressure_scale);
		maximum_pressure_l2_coarse_medium = std::max(maximum_pressure_l2_coarse_medium, coarse_medium.l2_relative);
		maximum_pressure_l2_medium_fine = std::max(maximum_pressure_l2_medium_fine, medium_fine.l2_relative);
		maximum_pressure_linf_coarse_medium = std::max(maximum_pressure_linf_coarse_medium, coarse_medium.linf_relative);
		maximum_pressure_linf_medium_fine = std::max(maximum_pressure_linf_medium_fine, medium_fine.linf_relative);
		const double roundoff = 128.0*std::numeric_limits<double>::epsilon();
		if (medium_fine.l2_relative > coarse_medium.l2_relative+roundoff
			|| medium_fine.linf_relative > coarse_medium.linf_relative+roundoff)
			throw std::runtime_error("temporal pressure field did not converge: "+std::string(field));
		std::cout << std::setprecision(17) << "temporal_pressure_field name=" << field
			<< " l2_coarse_medium=" << coarse_medium.l2_relative
			<< " l2_medium_fine=" << medium_fine.l2_relative
			<< " linf_coarse_medium=" << coarse_medium.linf_relative
			<< " linf_medium_fine=" << medium_fine.linf_relative << '\n';
	}
	double maximum_flow_l2_coarse_medium = 0.0, maximum_flow_l2_medium_fine = 0.0;
	double maximum_flow_linf_coarse_medium = 0.0, maximum_flow_linf_medium_fine = 0.0;
	for (const char* field : {"upstream_terminal_outward_flow_m3_s", "three_d_inlet_outward_flow_m3_s",
		"three_d_outlet_outward_flow_m3_s", "downstream_root_outward_flow_m3_s"}) {
		const auto coarse_values = TailPeriod(Extract(coarse.history, field), coarse_samples);
		const auto medium_values = TailPeriod(Extract(medium.history, field), medium_samples);
		const auto fine_values = TailPeriod(Extract(fine.history, field), fine_samples);
		const auto coarse_medium = CompareWaveforms(coarse_values,
			RestrictToCoarseTimes(medium_values, coarse_to_medium), coarse.macro_dt_s, kMeanFlowM3S);
		const auto medium_fine = CompareWaveforms(medium_values,
			RestrictToCoarseTimes(fine_values, medium_to_fine), medium.macro_dt_s, kMeanFlowM3S);
		maximum_flow_l2_coarse_medium = std::max(maximum_flow_l2_coarse_medium, coarse_medium.l2_relative);
		maximum_flow_l2_medium_fine = std::max(maximum_flow_l2_medium_fine, medium_fine.l2_relative);
		maximum_flow_linf_coarse_medium = std::max(maximum_flow_linf_coarse_medium, coarse_medium.linf_relative);
		maximum_flow_linf_medium_fine = std::max(maximum_flow_linf_medium_fine, medium_fine.linf_relative);
		const double roundoff = 128.0*std::numeric_limits<double>::epsilon();
		if (medium_fine.l2_relative > coarse_medium.l2_relative+roundoff
			|| medium_fine.linf_relative > coarse_medium.linf_relative+roundoff)
			throw std::runtime_error("temporal flow field did not converge: "+std::string(field));
		std::cout << std::setprecision(17) << "temporal_flow_field name=" << field
			<< " l2_coarse_medium=" << coarse_medium.l2_relative
			<< " l2_medium_fine=" << medium_fine.l2_relative
			<< " linf_coarse_medium=" << coarse_medium.linf_relative
			<< " linf_medium_fine=" << medium_fine.linf_relative << '\n';
	}
	const auto coarse_drop = TailPeriod(Extract(coarse.history, "external_pressure_drop_pa"), coarse_samples);
	const auto medium_drop = TailPeriod(Extract(medium.history, "external_pressure_drop_pa"), medium_samples);
	const auto fine_drop = TailPeriod(Extract(fine.history, "external_pressure_drop_pa"), fine_samples);
	const auto coarse_metric = MeasureWaveform(coarse_drop, coarse.macro_dt_s, pressure_scale);
	const auto medium_metric = MeasureWaveform(medium_drop, medium.macro_dt_s, pressure_scale);
	const auto fine_metric = MeasureWaveform(fine_drop, fine.macro_dt_s, pressure_scale);
	const double roundoff = 128.0*std::numeric_limits<double>::epsilon();
	if (maximum_pressure_l2_medium_fine > maximum_pressure_l2_coarse_medium+roundoff
		|| maximum_pressure_linf_medium_fine > maximum_pressure_linf_coarse_medium+roundoff
		|| maximum_flow_l2_medium_fine > maximum_flow_l2_coarse_medium+roundoff
		|| maximum_flow_linf_medium_fine > maximum_flow_linf_coarse_medium+roundoff)
		throw std::runtime_error("temporal refinement did not reduce nested waveform differences");
	const double amplitude_coarse_medium = std::abs(coarse_metric.amplitude-medium_metric.amplitude);
	const double amplitude_medium_fine = std::abs(medium_metric.amplitude-fine_metric.amplitude);
	const double phase_coarse_medium = std::abs(iga::test::HarmonicPhaseDifferenceRadians(
		coarse_metric.phase_rad, medium_metric.phase_rad));
	const double phase_medium_fine = std::abs(iga::test::HarmonicPhaseDifferenceRadians(
		medium_metric.phase_rad, fine_metric.phase_rad));
	if (!(fine_metric.amplitude > 1.0e-6*pressure_scale)
		|| amplitude_medium_fine > amplitude_coarse_medium+roundoff*pressure_scale
		|| phase_medium_fine > phase_coarse_medium+roundoff)
		throw std::runtime_error("temporal refinement did not converge pressure amplitude and phase");
	std::cout << std::setprecision(17) << "temporal_self_convergence pressure_l2_coarse_medium="
		<< maximum_pressure_l2_coarse_medium << " pressure_l2_medium_fine=" << maximum_pressure_l2_medium_fine
		<< " pressure_linf_coarse_medium=" << maximum_pressure_linf_coarse_medium
		<< " pressure_linf_medium_fine=" << maximum_pressure_linf_medium_fine
		<< " internal_q_l2_coarse_medium=" << maximum_flow_l2_coarse_medium
		<< " internal_q_l2_medium_fine=" << maximum_flow_l2_medium_fine
		<< " internal_q_linf_coarse_medium=" << maximum_flow_linf_coarse_medium
		<< " internal_q_linf_medium_fine=" << maximum_flow_linf_medium_fine
		<< " amplitude_pa=" << coarse_metric.amplitude << ',' << medium_metric.amplitude << ',' << fine_metric.amplitude
		<< " phase_rad=" << coarse_metric.phase_rad << ',' << medium_metric.phase_rad << ',' << fine_metric.phase_rad
		<< " amplitude_differences_pa=" << amplitude_coarse_medium << ',' << amplitude_medium_fine
		<< " phase_differences_rad=" << phase_coarse_medium << ',' << phase_medium_fine << '\n';
}

double ToleranceDifference(const RunResult& candidate, const RunResult& tight)
{
	const int samples_per_period = static_cast<int>(std::llround(kPeriodS/candidate.macro_dt_s));
	const auto candidate_drop = TailPeriod(Extract(candidate.history, "external_pressure_drop_pa"), samples_per_period);
	const auto tight_drop = TailPeriod(Extract(tight.history, "external_pressure_drop_pa"), samples_per_period);
	const auto candidate_root = TailPeriod(Extract(candidate.history, "upstream_root_pressure_pa"), samples_per_period);
	const auto tight_root = TailPeriod(Extract(tight.history, "upstream_root_pressure_pa"), samples_per_period);
	const double pressure_scale = (16.0*iga::test::kPi+iga::test::kSquareDuctResistanceCoefficient)*kMeanFlowM3S;
	return std::max(CompareWaveforms(candidate_drop, tight_drop, candidate.macro_dt_s, pressure_scale).linf_relative,
		CompareWaveforms(candidate_root, tight_root, candidate.macro_dt_s, pressure_scale).linf_relative);
}

int Run()
{
	CheckMetricUtilities();
	const double radius_m = std::sqrt(1.0/iga::test::kPi);
	assert(radius_m > 0.0);
	assert(std::abs(kDurationS/kPeriodS-4.0) < 1.0e-14);
	assert(std::abs(kPeriodS/0.01-16.0) < 1.0e-14);
	assert(std::abs(kPeriodS/0.005-32.0) < 1.0e-14);
	assert(std::abs(kPeriodS/0.0025-64.0) < 1.0e-14);
	const auto fixture_root = fs::temp_directory_path()/("tubularflowiga-pulsatile-temporal-fixture-"
		+std::to_string(static_cast<long long>(getpid())));
	std::error_code fixture_error;
	fs::remove_all(fixture_root, fixture_error);
	TemporaryDirectory fixture_cleanup{fixture_root};
	iga::test::WriteC2SquareDuctThreeDCase(fixture_root/"three_d", kTransverseElements, 0.0025, 256,
		kDensityKgM3, kDynamicViscosityPaS, kMeanFlowM3S, kAxialElements, 1.0);
	iga::test::WriteRigidOneDStraightCase(fixture_root/"upstream", 1.0, radius_m, 0.00125, 512,
		kDensityKgM3, kDynamicViscosityPaS, kMeanFlowM3S, "sinusoid", kPeriodS);
	const std::string upstream_configuration = ReadText(fixture_root/"upstream/simulation_config.json");
	assert(upstream_configuration.find("\"kind\":\"sinusoid\"") != std::string::npos);
	assert(upstream_configuration.find("\"period\":") != std::string::npos);
	if (std::getenv("TUBULARFLOWIGA_TEMPORAL_FIXTURE_ONLY")) return 0;
	const auto root = fs::temp_directory_path()/(
		"tubularflowiga-pulsatile-temporal-"+std::to_string(static_cast<long long>(getpid())));
	TemporaryDirectory cleanup{root};
	if (std::getenv("TUBULARFLOWIGA_KEEP_TEST_OUTPUT")) std::cout << "retained_test_output=" << root << '\n';
	const auto coarse = RunPulsatileCase(root, 0.01, kTightPressureTolerance);
	if (std::getenv("TUBULARFLOWIGA_TEMPORAL_SINGLE_CASE")) {
		CheckReferenceAgreement(coarse);
		CheckPeriodicSettling(coarse);
		return 0;
	}
	const auto medium = RunPulsatileCase(root, 0.005, kTightPressureTolerance);
	const auto tight = RunPulsatileCase(root, 0.0025, kTightPressureTolerance);
	CheckReferenceAgreement(coarse);
	CheckReferenceAgreement(medium);
	CheckReferenceAgreement(tight);
	CheckPeriodicSettling(coarse);
	CheckPeriodicSettling(medium);
	CheckPeriodicSettling(tight);
	CheckTemporalSelfConvergence(coarse, medium, tight);
	const auto loose = RunPulsatileCase(root, 0.0025, 1.0e-4);
	const auto intermediate = RunPulsatileCase(root, 0.0025, 1.0e-6);
	const double loose_error = ToleranceDifference(loose, tight);
	const double intermediate_error = ToleranceDifference(intermediate, tight);
	const double roundoff = 128.0*std::numeric_limits<double>::epsilon();
	if (intermediate_error > loose_error+roundoff)
		throw std::runtime_error("tighter coupling tolerance did not move the waveform toward the tight reference");
	std::cout << std::setprecision(17) << "coupling_tolerance_sensitivity loose_to_tight=" << loose_error
		<< " intermediate_to_tight=" << intermediate_error << " tight_to_tight=0"
		<< " tolerances=0.0001,9.9999999999999995e-07,1e-10\n";
	return 0;
}

} // namespace

int main()
{
	try {
		return Run();
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n';
		return 1;
	}
}
