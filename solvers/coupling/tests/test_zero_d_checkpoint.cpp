#include "ZeroDFlowCheckpoint.hpp"
#include "CoupledCheckpointBundle.hpp"
#include <iostream>
#include <sys/wait.h>

namespace fs = std::filesystem;
using namespace iga;
using checkpoint_metadata::Require;

namespace {
int rejections = 0;
template<class Work> void Reject(Work&& work)
{
	try { work(); } catch (const std::exception&) { ++rejections; return; }
	throw std::runtime_error("expected rejected 0D checkpoint");
}

ZeroDFlowModel Model(bool source)
{
	ZeroDFlowModel result; result.role = source ? ZeroDFlowRole::SourceReservoir : ZeroDFlowRole::TerminalRcr;
	if (source) result.source = {2.0, 4.0, 3.0};
	else result.terminal = {1.5, 5.0, 2.0, 7.0};
	return result;
}

ZeroDFlowDomainRuntime Runtime(bool source)
{
	const std::string id = source ? "source" : "terminal"; const auto model = Model(source);
	return ZeroDFlowDomainRuntime(id, model, {source ? 11.0 : 9.0}, {MakeZeroDFlowPort(id, model.role)});
}

void Trial(ZeroDFlowDomainRuntime& runtime, int index, double dt = 0.1)
{
	const DomainStepContext step{index, runtime.CommittedTime(), dt}; runtime.BeginStep(step);
	PortBoundaryData input; input.time_s = step.EndTime();
	if (runtime.Model().role == ZeroDFlowRole::SourceReservoir) input.mean_pressure_pa = index%3 == 0 ? 15.0 : 4.0+index*0.2;
	else input.outward_flow_m3_s = index%3 == 0 ? 0.75 : -2.0-index*0.1;
	runtime.SetPortInput("port", input); runtime.SolveTrial();
}

void Advance(ZeroDFlowDomainRuntime& runtime, int index)
{
	Trial(runtime, index); runtime.PrepareCommitStep(); runtime.FinalizeCommitStep();
}

std::string Digest(const std::string& value) { return coupled_checkpoint_detail::Hash(value); }
CoupledCheckpointEpoch Epoch()
{
	return {Digest("zero-d-step-4"), {}, {Digest("zero-d-case"), Digest("zero-d-configuration"), Digest("zero-d-execution"), 1}, 4, 0.4, 0.1};
}
std::vector<CoupledCheckpointShardSpec> Catalog()
{
	return {{"source", "zero-d-accepted-v1", 0}, {"terminal", "zero-d-accepted-v1", 0}};
}

void CodecTests()
{
	using namespace checkpoint_metadata;
	PortState port; port.time_s = 0.3; port.area_m2 = 1.5; port.outward_flow_m3_s = -0.0;
	port.mean_normal_traction_pa = -5.0; port.concentration = {{"dye", 2.0}, {"tracer", 3.0}};
	port.outward_species_flux = {{"tracer", -4.0}};
	Writer output; WritePort(output, port); Reader input(output.Bytes()); const auto restored = ReadPort(input); input.Finish();
	Writer copy; WritePort(copy, restored); Require(copy.Bytes() == output.Bytes(), "port codec did not preserve exact bytes");
	Require(!restored.mean_pressure_pa && !restored.total_pressure_pa && std::signbit(*restored.outward_flow_m3_s), "optional or signed zero was lost");
	Writer invalid; invalid.Unsigned(2); Reader presence(invalid.Bytes()); Reject([&] { presence.Optional(); });
	Writer too_long; too_long.Unsigned(UINT64_MAX); Reader text(too_long.Bytes()); Reject([&] { text.Text(); });
	Reader count(too_long.Bytes()); Reject([&] { count.Reals(); });
	Writer duplicate; duplicate.Unsigned(2); duplicate.Text("dye"); duplicate.Real(1); duplicate.Text("dye"); duplicate.Real(2);
	Reader duplicate_reader(duplicate.Bytes()); Reject([&] { duplicate_reader.Reals(); });
	Reject([&] { Writer large; large.Text(std::string(maximum_string_bytes+1, 'x')); });
	Reject([&] { Reader large(std::string(maximum_bytes+1, 'x')); });
}

void StateTests(bool source)
{
	auto original = Runtime(source); Reject([&] { original.CaptureCheckpointState(); });
	for (int i = 0; i < 4; ++i) Advance(original, i);
	const auto state = original.CaptureCheckpointState(); const auto encoded = SerializeZeroDFlowCheckpoint(state);
	Require(SerializeZeroDFlowCheckpoint(ParseZeroDFlowCheckpoint(encoded)) == encoded, "0D codec roundtrip differs");
	for (std::size_t n = 0; n < encoded.size(); ++n) Reject([&] { ParseZeroDFlowCheckpoint(std::string_view(encoded).substr(0, n)); });
	Reject([&] { ParseZeroDFlowCheckpoint(encoded+"x"); });
	for (int changed = 0; changed < 3; ++changed) {
		auto epoch = Epoch();
		if (changed == 0) ++epoch.accepted_steps;
		if (changed == 1) epoch.time_s += 0.1;
		if (changed == 2) epoch.dt_s += 0.1;
		Reject([&] { ParseZeroDFlowCheckpoint(encoded, epoch); });
	}
	auto fresh = Runtime(source); const auto initial_hash = fresh.CommittedStateIdentitySha256();
	for (int mutation = 0; mutation < 12; ++mutation) {
		auto bad = state;
		if (mutation == 0) bad.domain_id += "wrong";
		if (mutation == 1) bad.model_identity_sha256 = Digest("wrong");
		if (mutation == 2) bad.accepted_step.step_index = -1;
		if (mutation == 3) bad.port.time_s += 1;
		if (mutation == 4) bad.port.mean_pressure_pa.reset();
		if (mutation == 5) bad.port.area_m2 = 1;
		if (mutation == 6) bad.state.stored_pressure_pa = std::numeric_limits<double>::quiet_NaN();
		if (mutation == 7) bad.accounting.final_stored_volume_m3 += 1;
		if (mutation == 8) bad.accounting.outward_graph_port_amount_m3 += 1;
		if (mutation == 9) bad.accounting.initial_stored_volume_m3 += 1;
		if (mutation == 10) bad.accounting.residual_m3 += 1;
		if (mutation == 11) bad.port.concentration["dye"] = 1;
		Reject([&] { fresh.RestoreCheckpointState(bad); });
		Require(fresh.CommittedStateIdentitySha256() == initial_hash && fresh.CommittedTime() == 0
			&& fresh.CommittedStepIndex() == -1 && fresh.CommittedStepCount() == 0
			&& !fresh.CommittedPortState() && !fresh.CommittedStepAccounting(), "failed restore changed accepted state");
	}
	fresh.RestoreCheckpointState(ParseZeroDFlowCheckpoint(encoded));
	Require(SerializeZeroDFlowCheckpoint(fresh.CaptureCheckpointState()) == encoded
		&& fresh.CommittedTime() == original.CommittedTime() && fresh.CommittedStepIndex() == 3
		&& fresh.CommittedStepCount() == 4, "complete accepted restore differs");
	Reject([&] { fresh.RestoreCheckpointState(state); });
	const DomainStepContext pending{4, original.CommittedTime(), 0.2}; original.BeginStep(pending);
	Reject([&] { original.CaptureCheckpointState(); }); original.AbortStep();
	Trial(original, 4, 0.2); Reject([&] { original.CaptureCheckpointState(); });
	original.PrepareCommitStep(); Reject([&] { original.CaptureCheckpointState(); }); original.AbortStep();
	Require(SerializeZeroDFlowCheckpoint(original.CaptureCheckpointState()) == encoded, "aborted dt replaced accepted context");
	for (int i = 4; i < 10; ++i) {
		Advance(original, i); Advance(fresh, i);
		Require(SerializeZeroDFlowCheckpoint(original.CaptureCheckpointState())
			== SerializeZeroDFlowCheckpoint(fresh.CaptureCheckpointState()), "restored continuation differs");
	}
	auto last = state; last.accepted_step.step_index = std::numeric_limits<int>::max();
	auto exhausted = Runtime(source); exhausted.RestoreCheckpointState(last);
	Reject([&] { exhausted.BeginStep({std::numeric_limits<int>::max(), exhausted.CommittedTime(), 0.1}); });
}

void Save(const fs::path& root)
{
	Require(fs::create_directory(root), "test needs an unused directory"); const auto epoch = Epoch();
	CreateCoupledCheckpointEpoch(root, epoch); std::vector<CoupledCheckpointShard> receipts;
	for (const bool source : {true, false}) {
		auto runtime = Runtime(source); for (int i = 0; i < 4; ++i) Advance(runtime, i);
		const auto snapshot = runtime.CaptureCheckpointState();
		Require(snapshot.accepted_step.EndTime() == epoch.time_s, "0D and bundle clock differ");
		const auto bytes = SerializeZeroDFlowCheckpoint(snapshot); const auto spec = Catalog().at(source ? 0 : 1);
		receipts.push_back(WriteCoupledCheckpointShard(root, epoch, spec, bytes.size(), [&](auto& output) { output.Write(bytes.data(), bytes.size()); }));
	}
	PublishCoupledCheckpoint(root, epoch, Catalog(), receipts);
}

void Resume(const fs::path& root)
{
	const auto epoch = Epoch(); const auto manifest = LoadCoupledCheckpoint(root, epoch.id, epoch.compatibility, Catalog());
	for (const bool source : {true, false}) {
		std::string bytes; ReadCoupledCheckpointShard(root, manifest, source ? 0 : 1,
			[&](const void* data, std::size_t count) { bytes.append(static_cast<const char*>(data), count); });
		auto resumed = Runtime(source); resumed.RestoreCheckpointState(ParseZeroDFlowCheckpoint(bytes, manifest.epoch));
		auto reference = Runtime(source);
		for (int i = 0; i < 10; ++i) {
			Advance(reference, i); if (i >= 4) Advance(resumed, i);
			if (i >= 3) Require(SerializeZeroDFlowCheckpoint(reference.CaptureCheckpointState())
				== SerializeZeroDFlowCheckpoint(resumed.CaptureCheckpointState()), "new-process continuation differs from uninterrupted state");
		}
	}
	std::cout << "zero_d_checkpoint resume=passed roles=2 subsequent_steps=6 exact_state_and_accounting=1\n";
}
} // namespace

int main(int argc, char** argv)
{
	try {
		if (argc == 3 && std::string(argv[1]) == "--resume") { Resume(argv[2]); return 0; }
		Require(argc == 2, "test requires an unused output directory");
		CodecTests(); StateTests(true); StateTests(false); Save(argv[1]);
		const auto pid = ::fork(); Require(pid >= 0, "fork failed");
		if (pid == 0) { ::execl(argv[0], argv[0], "--resume", argv[1], static_cast<char*>(nullptr)); ::_exit(99); }
		int status = 0; Require(::waitpid(pid, &status, 0) == pid && WIFEXITED(status) && WEXITSTATUS(status) == 0, "fresh executable restart failed");
		std::cout << "zero_d_checkpoint status=passed rejections=" << rejections << '\n'; return 0;
	} catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
