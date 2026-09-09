#include "OneDAcceptedCheckpoint.hpp"
#include "CoupledCheckpointBundle.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <sys/wait.h>

namespace fs = std::filesystem;
using namespace iga;
using checkpoint_metadata::Require;

namespace {
int rejections = 0;
using Clock = std::chrono::steady_clock;
std::array<double, 5> solve_seconds{};
std::array<int, 5> solve_calls{};
double Seconds(Clock::time_point start) { return std::chrono::duration<double>(Clock::now()-start).count(); }
void PrintSolveTimes(const char* process)
{
	for (std::size_t mode = 0; mode < solve_seconds.size(); ++mode)
		std::cout << "one_d_checkpoint process=" << process << " mode=" << mode << " solve_calls=" << solve_calls[mode]
			<< " solve_seconds=" << solve_seconds[mode] << '\n';
}
template<class Work> void Reject(Work&& work)
{
	try { work(); } catch (const std::exception&) { ++rejections; return; }
	throw std::runtime_error("expected rejected 1d checkpoint");
}
std::string Digest(const std::string& value) { return coupled_checkpoint_detail::Hash(value); }
std::string FileDigest(const fs::path& path)
{
	std::ifstream input(path); Require(bool(input), "fixture input missing"); return Digest(ReadCheckedText(input));
}
std::string Configuration()
{
	return R"json({
 "schema_version":3,"dimension":"1d",
 "geometry":{"kind":"swc_network","file":"tree.swc","length_scale_to_m":1.0},
 "fields":[{"name":"area","kind":"scalar"},{"name":"flow_rate","kind":"scalar"},{"name":"pressure","kind":"pressure"},{"name":"signal","kind":"scalar","initial_value":1.0}],
 "time":{"dt":0.001,"steps":20,"output_every":1},
 "temporal_functions":[{"name":"signal_wave","kind":"periodic_table","units":"1","file":"signal.csv","period":2.0,"interpolation":"linear"}],
 "equation_systems":[
  {"name":"flow","kind":"network_flow_1d","unknowns":["area","flow_rate","pressure"],"model":"rigid","scheme":"steady_poiseuille","dynamic_viscosity":0.004,"density":1060.0,"discretization":{"cells_per_segment":1}},
  {"name":"transport","kind":"network_transport_1d","unknowns":["signal"],"flow_system":"flow","species":[{"field":"signal","diffusivity":0.0}]}
 ],
 "boundaries":[
  {"name":"inlet","role":"inlet","conditions":[{"field":"flow_rate","type":"dirichlet","value":1e-9},{"field":"signal","type":"dirichlet","value":1.0,"waveform":"signal_wave"}]},
  {"name":"outlet","role":"outlet","conditions":[{"field":"pressure","type":"pressure","value":0.0}]}
 ],
 "physiology":{"enabled":true,"oxygen_capacity":{"enabled":true,"hematocrit_percent":40.0,"hemoglobin_g_dl":13.6},"vasodilation":{"enabled":true,"field":"signal","emax_radius_fraction":0.1,"ec50":1.0,"relaxation_tau":0.01}}
})json";
}


OneDFlowRuntime Runtime(const fs::path& root, int mode, bool bound = true)
{
	auto config = ParseOneDConfiguration(Configuration());
	config.physiology.vasodilation = mode == 0;
	config.coupling.perfusate.oxygen.hematocrit_percent = 40;
	config.coupling.perfusate.oxygen.hemoglobin_g_dl = 14;
	if (mode == 2) {
		auto& rcr = config.boundaries.back().conditions.front(); rcr.type = "windkessel_rcr";
		rcr.proximal_resistance = 1e8; rcr.distal_resistance = 1e9; rcr.capacitance = 1e-10;
		rcr.reference_pressure = 2; rcr.initial_pressure = 3;
	}
	if (mode == 3) {
		config.flow_systems.front().model = OneDFlowModel::Compliant;
		config.flow_systems.front().scheme = OneDFlowScheme::ExplicitRusanov;
	}
	const auto flow = config.flow_systems.front();
	auto network = ReadOneDNetwork(root/"tree.swc", 1.0, 2, flow.dynamic_viscosity);
	OneDFlowRuntime runtime(config, flow, network, ResolveOneDInlet(config), root, {}, {},
		bound ? Digest(Configuration()+std::to_string(mode)+"mode-recipe-v1"+FileDigest(root/"tree.swc")+FileDigest(root/"signal.csv")) : "");
	runtime.InitializeOpenLoop(1e-9); return runtime;
}

void Trial(OneDFlowRuntime& runtime, int index, int mode)
{
	const double time = runtime.FlowState().physical_time;
	runtime.BeginStep(time, 0.002);
	VascularInletState inlet; inlet.time_s = time+0.002; inlet.has_flow = true;
	inlet.flow_m3_s = index%4 == 3 && mode != 2 ? -1e-9 : 1e-9*(1+0.2*index);
	inlet.has_pressure = true; inlet.pressure_pa = 1.25;
	inlet.has_temperature = index%2 == 0; inlet.temperature_c = inlet.has_temperature ? 37 : NAN;
	inlet.has_hematocrit = index == 0 || index == 1 || index == 5;
	inlet.hematocrit_percent = index == 5 ? 28 : 35;
	if (index%3 == 0) inlet.species = {{"signal", 2.0+0.1*index}};
	inlet.metadata = {{"", "unnamed"}, {"source", "checkpoint fixture"}, {"step", std::to_string(index)}};
	if (mode == 4) runtime.SetConfiguredOpenLoopInlet();
	else runtime.SetCoupledInlet(inlet);
	if (mode != 2) {
		PortBoundaryData pressure; pressure.time_s = inlet.time_s; pressure.mean_pressure_pa = mode == 3 ? 0.0 : 0.5*index; pressure.concentration = {{"signal", 3.0}};
		runtime.SetPortInput("outlet:2", pressure); runtime.SetPortInput("outlet:3", pressure);
	}
	const auto solve_start = Clock::now();
	if (mode == 1) {
		runtime.SolveHydraulicTrial();
		std::map<int, std::map<std::string, double>> outlets;
		for (int node : runtime.Network().outlet_nodes) outlets[node] = {{"signal", 3.0}};
		runtime.SolveStagedTransportTrial({{"signal", 2.0+0.1*index}}, outlets);
	} else runtime.SolveTrial();
	solve_seconds[mode] += Seconds(solve_start); ++solve_calls[mode];
}
void Advance(OneDFlowRuntime& runtime, int index, int mode)
{
	Trial(runtime, index, mode); runtime.PrepareCommitStep(); runtime.FinalizeCommitStep();
}
std::string Fields(const OneDFlowCheckpointState& state)
{
	std::string bytes;
	WriteOneDAcceptedCheckpointFields(state, [&](const void* data, std::size_t size) { bytes.append(static_cast<const char*>(data), size); });
	Require(bytes.size() == OneDAcceptedCheckpointFieldBytes(state), "wrong field byte count"); return bytes;
}
std::string Fingerprint(const OneDFlowRuntime& runtime)
{
	const auto state = runtime.CaptureCheckpointState();
	checkpoint_metadata::Writer extras;
	for (const auto& segment : runtime.Network().segments) { extras.Real(segment.area0); extras.Real(segment.resistance); }
	checkpoint_metadata::WritePort(extras, runtime.GetPortState("root"));
	checkpoint_metadata::WritePort(extras, runtime.GetPortState("outlet:2"));
	checkpoint_metadata::WritePort(extras, runtime.GetPortState("outlet:3"));
	for (const auto& entry : runtime.GetSpeciesStepAccounting()) {
		extras.Text(entry.first); extras.Real(entry.second.initial_mass); extras.Real(entry.second.final_mass);
		extras.Real(entry.second.root_outward_amount); extras.Real(entry.second.source_amount); extras.Real(entry.second.balance_residual);
	}
	return Digest(SerializeOneDAcceptedCheckpointMetadata(state)+Fields(state)+extras.Bytes());
}
CoupledCheckpointEpoch Epoch(int mode)
{
	return {Digest("one-d-step-3-mode-"+std::to_string(mode)), {}, {Digest("one-d-case"), Digest("one-d-config"), Digest("one-d-execution"), 1}, 3, 0.006, 0.002};
}
std::vector<CoupledCheckpointShardSpec> Catalog()
{
	return {{"fields", "one-d-fields-v1", 0}, {"metadata", "one-d-accepted-v1", 0}};
}

void StateTests(const fs::path& root, int mode)
{
	auto original = Runtime(root, mode); Reject([&] { original.CaptureCheckpointState(); });
	for (int i = 0; i < 3; ++i) Advance(original, i, mode);
	const auto state = original.CaptureCheckpointState(); const auto fingerprint = Fingerprint(original);
	Require(state.accepted_macro_steps == 3 && state.flow.completed_step == 6, "integer subcycling counter differs");
	if (mode == 0) {
		Require(state.segment_radii_m.front() != original.Network().segments.front().baseline_radius0, "fixture did not dilate");
		Require(!state.last_inlet.has_hematocrit && state.blood_state[0] == 35 && state.blood_state[3] == 12.25,
			"fixture did not retain previous hematocrit/hemoglobin");
		Require(state.transports.front().species.front().inlet_waveform.empty() && state.last_inlet.species.empty(), "fixture did not clear inlet waveform");
	}
	auto fresh = Runtime(root, mode);
	const auto metadata = SerializeOneDAcceptedCheckpointMetadata(state), fields = Fields(state);
	for (std::size_t n = 0; n < metadata.size(); ++n)
		Reject([&] { ParseOneDAcceptedCheckpointMetadata(std::string_view(metadata).substr(0, n), fresh, Epoch(mode)); });
	Reject([&] { ParseOneDAcceptedCheckpointMetadata(metadata+"x", fresh, Epoch(mode)); });
	for (int change = 0; change < 3; ++change) {
		auto epoch = Epoch(mode); if (change == 0) ++epoch.accepted_steps; if (change == 1) epoch.time_s += 1; if (change == 2) epoch.dt_s += 1;
		Reject([&] { ParseOneDAcceptedCheckpointMetadata(metadata, fresh, epoch); });
	}
	for (int change = 0; change < 25; ++change) {
		auto bad = state;
		if (change == 0) bad.configuration_identity_sha256 = Digest("wrong");
		if (change == 1) ++bad.accepted_macro_steps;
		if (change == 2) ++bad.flow.completed_step;
		if (change == 3) bad.last_macro_dt_s = 0;
		if (change == 4) bad.last_macro_start_s = -1;
		if (change == 5) bad.flow.physical_time += 1;
		if (change == 6) bad.flow.internal_substeps = -1;
		if (change == 7) bad.flow.area.pop_back();
		if (change == 8) bad.flow.area.front() = 0;
		if (change == 9) bad.flow.flow.front() = NAN;
		if (change == 10) bad.segment_radii_m.front() = 0;
		if (change == 11) bad.flow.outlets.front().resistance += 1;
		if (change == 12) bad.flow.outlets.front().node = -1;
		if (change == 13) bad.last_inlet.time_s += 1;
		if (change == 14) bad.blood_state[0] = 101;
		if (change == 15) bad.transports.front().species.front().inlet_waveform = "wrong.csv";
		if (change == 16) bad.transports.front().species.front().step_accounting_valid = false;
		if (change == 17) bad.transports.front().species.front().outlet_native_flux.clear();
		if (change == 18) bad.transports.front().species.front().definition.volume_source += 1;
		if (change == 19) bad.last_inlet.species["unknown"] = 1;
		if (change == 20) bad.last_macro_dt_s = std::numeric_limits<double>::max();
		if (change == 21) bad.transports.front().species.front().step_outlet_native_amount.begin()->second = INFINITY;
		if (change == 22) bad.flow.outlets.front().capacitor_pressure = INFINITY;
		if (change == 23) bad.segment_radii_m.front() = std::numeric_limits<double>::min();
		if (change == 24) bad.transports.front().species.front().concentration.front() = NAN;
		Reject([&] { fresh.RestoreCheckpointState(bad); });
		Require(fresh.FlowState().completed_step == 0 && fresh.FlowState().physical_time == 0
			&& fresh.Network().segments.front().radius0 == fresh.Network().segments.front().baseline_radius0
			&& fresh.Transports().front().species.front().inlet_waveform == "signal_wave", "failed restore changed fresh runtime");
	}
	for (std::size_t n = 0; n < fields.size(); ++n) {
		auto candidate = ParseOneDAcceptedCheckpointMetadata(metadata, fresh, Epoch(mode));
		OneDAcceptedCheckpointFieldsReader reader(candidate); reader.Consume(fields.data(), n); Reject([&] { reader.Finish(); });
	}
	{
		auto candidate = ParseOneDAcceptedCheckpointMetadata(metadata, fresh, Epoch(mode));
		OneDAcceptedCheckpointFieldsReader reader(candidate); const auto extra = fields+"x";
		Reject([&] { reader.Consume(extra.data(), extra.size()); }); Reject([&] { reader.Finish(); });
	}
	{
		auto candidate = ParseOneDAcceptedCheckpointMetadata(metadata, fresh, Epoch(mode));
		OneDAcceptedCheckpointFieldsReader reader(candidate);
		const unsigned char nonfinite[8] = {0,0,0,0,0,0,0xf0,0x7f};
		Reject([&] { reader.Consume(nonfinite, 8); }); Reject([&] { reader.Consume(fields.data(), fields.size()); });
	}
	auto candidate = ParseOneDAcceptedCheckpointMetadata(metadata, fresh, Epoch(mode));
	OneDAcceptedCheckpointFieldsReader reader(candidate);
	for (std::size_t offset = 0; offset < fields.size();) {
		const auto n = std::min<std::size_t>(1+offset%19, fields.size()-offset); reader.Consume(fields.data()+offset, n); offset += n;
	}
	reader.Finish(); fresh.RestoreCheckpointState(std::move(candidate));
	Require(Fingerprint(fresh) == fingerprint, "accepted snapshot differs after roundtrip");
	Reject([&] { fresh.RestoreCheckpointState(state); });
	original.BeginStep(original.FlowState().physical_time, 0.001);
	Reject([&] { original.CaptureCheckpointState(); }); Reject([&] { original.RestoreCheckpointState(state); }); original.AbortStep();
	Require(Fingerprint(original) == fingerprint, "aborted different dt altered accepted clock");
	Trial(original, 3, mode); Reject([&] { original.CaptureCheckpointState(); });
	original.PrepareCommitStep(); Reject([&] { original.CaptureCheckpointState(); }); original.AbortStep();
	Require(Fingerprint(original) == fingerprint, "aborted prepared state changed accepted snapshot");
	for (int i = 3; i < 10; ++i) {
		Advance(original, i, mode); Advance(fresh, i, mode);
		Require(Fingerprint(fresh) == Fingerprint(original), "subsequent state/ports/accounting differ");
	}
	// Existing transport has no positivity clamp. Interior finite values and
	// accounting must roundtrip without introducing a new numerical gate.
	if (mode == 0) {
		auto signed_state = state; signed_state.transports.front().species.front().concentration.front() = -0.25;
		signed_state.transports.front().species.front().step_initial_mass = -1e-8;
		auto signed_target = Runtime(root, mode); signed_target.RestoreCheckpointState(signed_state);
		Require(Fields(signed_target.CaptureCheckpointState()) == Fields(signed_state), "signed interior field was changed");
	}
	// Preserve the 64-bit explicit counter beyond the old legacy int encoding.
	auto large_count = state; large_count.flow.internal_substeps = static_cast<long long>(INT_MAX)+17;
	auto counter = Runtime(root, mode); counter.RestoreCheckpointState(large_count);
	auto counter_target = Runtime(root, mode);
	Require(ParseOneDAcceptedCheckpointMetadata(SerializeOneDAcceptedCheckpointMetadata(counter.CaptureCheckpointState()),
		counter_target, Epoch(mode)).flow.internal_substeps == large_count.flow.internal_substeps, "64-bit internal counter lost");
}

void Save(const fs::path& root, int mode)
{
	const auto directory = root/("bundle-"+std::to_string(mode)); Require(fs::create_directory(directory), "bundle directory exists");
	auto runtime = Runtime(root, mode); for (int i = 0; i < 3; ++i) Advance(runtime, i, mode);
	const auto snapshot = runtime.CaptureCheckpointState(); const auto metadata = SerializeOneDAcceptedCheckpointMetadata(snapshot);
	const auto write_start = Clock::now();
	const auto epoch = Epoch(mode); CreateCoupledCheckpointEpoch(directory, epoch); const auto catalog = Catalog();
	std::vector<CoupledCheckpointShard> receipts;
	receipts.push_back(WriteCoupledCheckpointShard(directory, epoch, catalog[0], OneDAcceptedCheckpointFieldBytes(snapshot), [&](auto& output) {
		WriteOneDAcceptedCheckpointFields(snapshot, [&](const void* bytes, std::size_t size) { output.Write(bytes, size); });
	}));
	receipts.push_back(WriteCoupledCheckpointShard(directory, epoch, catalog[1], metadata.size(), [&](auto& output) { output.Write(metadata.data(), metadata.size()); }));
	PublishCoupledCheckpoint(directory, epoch, catalog, receipts);
	std::cout << "one_d_checkpoint mode=" << mode << " bundle_write_publish_seconds=" << Seconds(write_start) << '\n';
	std::ofstream reference(root/("reference-"+std::to_string(mode)+".txt")); reference << Fingerprint(runtime) << '\n';
	for (int i = 3; i < 10; ++i) { Advance(runtime, i, mode); reference << Fingerprint(runtime) << '\n'; }
	Require(bool(reference), "reference write failed");
}
void Resume(const fs::path& root)
{
	for (int mode = 0; mode < 5; ++mode) {
		const auto directory = root/("bundle-"+std::to_string(mode)); const auto epoch = Epoch(mode);
		const auto load_start = Clock::now();
		const auto manifest = LoadCoupledCheckpoint(directory, epoch.id, epoch.compatibility, Catalog());
		std::string metadata;
		ReadCoupledCheckpointShard(directory, manifest, 1, [&](const void* bytes, std::size_t size) {
			Require(size <= checkpoint_metadata::maximum_bytes-metadata.size(), "metadata shard exceeds limit"); metadata.append(static_cast<const char*>(bytes), size);
		});
		auto runtime = Runtime(root, mode); auto candidate = ParseOneDAcceptedCheckpointMetadata(metadata, runtime, manifest.epoch);
		Require(manifest.shards[0].payload_bytes == OneDAcceptedCheckpointFieldBytes(candidate), "field shard size differs");
		OneDAcceptedCheckpointFieldsReader reader(candidate);
		ReadCoupledCheckpointShard(directory, manifest, 0, [&](const void* bytes, std::size_t size) { reader.Consume(bytes, size); });
		reader.Finish(); runtime.RestoreCheckpointState(std::move(candidate));
		std::cout << "one_d_checkpoint mode=" << mode << " load_verify_restore_seconds=" << Seconds(load_start) << '\n';
		std::ifstream reference(root/("reference-"+std::to_string(mode)+".txt")); std::string expected;
		for (int i = 2; i < 10; ++i) {
			if (i > 2) Advance(runtime, i, mode);
			Require(bool(reference >> expected) && Fingerprint(runtime) == expected, "fresh-process continuation differs");
		}
		Require(!(reference >> expected), "unexpected reference steps");
	}
	PrintSolveTimes("resume");
	std::cout << "one_d_checkpoint resume=passed modes=5 subsequent_steps=7 exact_fields_ports_accounting=1\n";
}
void LargeStreamTest()
{
	OneDFlowCheckpointState source; source.flow.area.assign(3*1024*1024, 1.25);
	auto target = source; std::fill(target.flow.area.begin(), target.flow.area.end(), 0);
	OneDAcceptedCheckpointFieldsReader reader(target); std::size_t maximum_chunk = 0;
	const auto start = Clock::now();
	WriteOneDAcceptedCheckpointFields(source, [&](const void* bytes, std::size_t size) {
		maximum_chunk = std::max(maximum_chunk, size); reader.Consume(bytes, size);
	}); reader.Finish();
	Require(source.flow.area == target.flow.area && maximum_chunk <= 65536, "large streaming roundtrip differs");
	std::cout << "one_d_checkpoint large_field_bytes=" << OneDAcceptedCheckpointFieldBytes(source) << " maximum_chunk=" << maximum_chunk << " encode_decode_seconds=" << Seconds(start) << '\n';
}
} // namespace

int main(int argc, char** argv)
{
	try {
		if (argc == 3 && std::string(argv[1]) == "--resume") { Resume(argv[2]); return 0; }
		Require(argc == 2, "test requires unused output directory"); const fs::path root(argv[1]);
		Require(fs::create_directory(root), "test directory exists");
		{ std::ofstream tree(root/"tree.swc"); tree << "1 2 0 0 0 0.001 -1\n2 2 0.01 0.01 0 0.001 1\n3 2 0.01 -0.01 0 0.001 1\n"; }
		{ std::ofstream waveform(root/"signal.csv"); waveform << "time,value\n0,1\n1,2\n"; }
		for (int mode = 0; mode < 5; ++mode) { StateTests(root, mode); Save(root, mode); }
		auto unbound = Runtime(root, 0, false); Advance(unbound, 0, 0); Reject([&] { unbound.CaptureCheckpointState(); });
		LargeStreamTest();
		const auto pid = ::fork(); Require(pid >= 0, "fork failed");
		if (pid == 0) { ::execl(argv[0], argv[0], "--resume", argv[1], static_cast<char*>(nullptr)); ::_exit(99); }
		int status = 0; Require(::waitpid(pid, &status, 0) == pid && WIFEXITED(status) && WEXITSTATUS(status) == 0, "fresh executable failed");
		PrintSolveTimes("parent");
		std::cout << "one_d_checkpoint status=passed rejections=" << rejections << '\n'; return 0;
	} catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
