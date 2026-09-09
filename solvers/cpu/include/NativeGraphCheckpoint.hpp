#ifndef IGA_NATIVE_GRAPH_CHECKPOINT_HPP
#define IGA_NATIVE_GRAPH_CHECKPOINT_HPP

#include "BodyFittedCheckpointBundle.hpp"
#include "CollectiveCheckpointReceipts.hpp"
#include "CollectiveAssetInput.hpp"
#include "CollectivePetscOptions.hpp"
#include "GraphAcceptedHistory.hpp"
#include "CheckpointRecordStream.hpp"
#include "OneDAcceptedCheckpoint.hpp"
#include "PressureFlowCheckpointControls.hpp"
#include "ZeroDFlowCheckpoint.hpp"
#include <numeric>

namespace iga {
inline std::string GraphCheckpointDomainPrefix(const std::string& id)
{
	return coupled_checkpoint_detail::Hash("native-graph-domain-v1:"+id);
}
inline std::string GraphCheckpointDomainIdentity(const CoupledCheckpointCompatibility& identity, const std::string& id, const char* role)
{
	checkpoint_metadata::Writer out; out.Text(identity.configuration_sha256); out.Text(identity.execution_sha256); out.Text(id); out.Text(role);
	return coupled_checkpoint_detail::Hash(out.Bytes());
}
inline CoupledCheckpointCompatibility BuildNativeGraphCheckpointIdentity(MPI_Comm comm,
	const std::map<std::string, std::string>& texts, const AssetFileCatalog& assets,
	const std::set<std::string>& application_options, int maximum_newton, const std::array<double, 3>& gates)
{
	CoupledCheckpointCompatibility result; int size = 0, world_rank = 0; MPI_Comm_size(comm, &size); MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	std::vector<int> mapping;
	CollectiveLocalStage(comm, "graph checkpoint identity allocation", [&] { mapping.resize(size); });
	MPI_Allgather(&world_rank, 1, MPI_INT, mapping.data(), 1, MPI_INT, comm);
	CollectiveLocalStage(comm, "graph checkpoint identity", [&] {
		Sha256 input; const auto append = [&](const std::string& value) { input.AppendLittleEndian64(value.size()); input.Append(value.data(), value.size()); };
		append("native-graph-inputs-v1"); append(std::to_string(texts.size()));
		for (const auto& item : texts) { append(item.first); append(item.second); }
		append(std::to_string(assets.size()));
		for (const auto& item : assets) { append(item.first); append(ReadAssetFingerprint(item.second)); }
		result.case_sha256 = input.Hex();
		checkpoint_metadata::Writer config; config.Text("native-graph-numerics-v1"); config.Text(result.case_sha256); config.Unsigned(maximum_newton);
		// Includes the native runner's fixed nonlinear and conservation gates.
		for (double gate : gates) config.Real(gate);
		result.configuration_sha256 = coupled_checkpoint_detail::Hash(config.Bytes());
		checkpoint_metadata::Writer execution; execution.Text("native-graph-execution-v1"); execution.Unsigned(size);
		for (int rank : mapping) execution.Unsigned(rank);
		// Content identity also covers dirty builds and source archives without Git.
		std::string source_identity;
#ifdef IGA_NATIVE_CHECKPOINT_SOURCE_SHA256
		source_identity = IGA_NATIVE_CHECKPOINT_SOURCE_SHA256;
#endif
#ifdef IGA_COUPLED_CHECKPOINT_TESTING
		if (const char* override_identity = std::getenv("IGA_NATIVE_CHECKPOINT_TEST_SOURCE_SHA256")) source_identity = override_identity;
#endif
		checkpoint_metadata::Require(source_identity.size() == 64 && source_identity.find_first_not_of("0123456789abcdef") == std::string::npos,
			"checkpoint requires a build-time IGA_NATIVE_CHECKPOINT_SOURCE_SHA256; rebuild using the native Makefile");
		execution.Text(source_identity);
#ifdef _OPENMP
		execution.Unsigned(_OPENMP);
#else
		execution.Unsigned(0);
#endif
#ifdef __FAST_MATH__
		execution.Unsigned(1);
#else
		execution.Unsigned(0);
#endif
		execution.Unsigned(sizeof(PetscInt)); execution.Unsigned(sizeof(PetscScalar)); execution.Text(__VERSION__);
		char petsc[256] = {}; if (PetscGetVersion(petsc, sizeof(petsc))) throw std::runtime_error("cannot inspect PETSc version"); execution.Text(petsc);
		char mpi[MPI_MAX_LIBRARY_VERSION_STRING] = {}; int length = 0; MPI_Get_library_version(mpi, &length); execution.Text(std::string(mpi, length));
		for (const char* name : {"OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS"}) {
			execution.Text(name); const char* value = std::getenv(name); execution.Text(value ? value : "");
		}
		const auto petsc_options = CapturePetscOptions(nullptr, application_options);
		execution.Text(coupled_checkpoint_detail::Hash(petsc_options)); result.execution_sha256 = coupled_checkpoint_detail::Hash(execution.Bytes()); result.ranks = size;
	});
	RequireCollectiveSameText(comm, "graph checkpoint input identity", result.case_sha256);
	RequireCollectiveSameText(comm, "graph checkpoint configuration identity", result.configuration_sha256);
	RequireCollectiveSameText(comm, "graph checkpoint execution identity", result.execution_sha256);
	return result;
}

struct NativeGraphCheckpointOwners {
	std::map<std::string, OneDFlowRuntime*> one_d;
	std::map<std::string, TransientFlowRuntime*> flow;
	std::map<std::string, TransientTransportRuntime*> transport;
	std::map<std::string, ZeroDFlowDomainRuntime*> zero_d;
};

// The runner creates these owners outside any running step. Restore may dirty
// them after a validated PETSc staging failure; the whole unpublished graph is
// then discarded. No step/output is exposed before Restore returns on all ranks.
class NativeGraphCheckpoint {
public:
	NativeGraphCheckpoint(MPI_Comm comm, const SimulationGraph& graph, const CoupledCheckpointCompatibility& identity,
		NativeGraphCheckpointOwners&& owners, bool species)
		: comm_(comm), graph_(graph), owners_(std::move(owners)), species_(species)
	{
		MPI_Comm_rank(comm_, &rank_);
		CollectiveLocalStage(comm_, "native graph checkpoint catalog", [&] {
			identity_ = identity;
			for (const auto& domain : graph.Domains()) {
				const auto prefix = GraphCheckpointDomainPrefix(domain.first);
				ids_.emplace(domain.first, std::array<std::string, 3>{prefix+".zero-d", prefix+".one-d.metadata", prefix+".one-d.fields"});
			}
			catalog_ = {{"graph.controls", "pressure-flow-controls-v1", 0}, {"graph.history", species_ ? "graph-species-history-v1" : "graph-flow-history-v1", 0}};
			for (const auto& item : owners_.one_d) {
				const auto prefix = GraphCheckpointDomainPrefix(item.first);
				catalog_.push_back({prefix+".one-d.metadata", "one-d-accepted-v1", 0}); catalog_.push_back({prefix+".one-d.fields", "one-d-fields-v1", 0});
			}
			for (const auto& item : owners_.zero_d) catalog_.push_back({GraphCheckpointDomainPrefix(item.first)+".zero-d", "zero-d-accepted-v1", 0});
			for (const auto& item : owners_.flow) {
				BodyFittedCheckpointLayout layout{GraphCheckpointDomainPrefix(item.first), {}, owners_.transport.count(item.first) != 0};
				layout.world_ranks.resize(identity_.ranks); std::iota(layout.world_ranks.begin(), layout.world_ranks.end(), 0);
				const auto domain = BodyFittedCheckpointCatalog(layout, identity_.ranks); catalog_.insert(catalog_.end(), domain.begin(), domain.end()); layouts_.emplace(item.first, std::move(layout));
			}
			checkpoint_metadata::Require(owners_.one_d.size()+owners_.zero_d.size()+owners_.flow.size() == graph.Domains().size(), "checkpoint owner catalog incomplete");
			std::sort(catalog_.begin(), catalog_.end(), [](const auto& a, const auto& b) { return a.id < b.id; });
			coupled_checkpoint_detail::ValidateCatalog(catalog_, identity_.ranks);
		});
	}
	CoupledCheckpointEpoch Save(const std::filesystem::path& root, const DomainStepContext& step,
		const PressureFlowCheckpointControls& controls, const std::vector<GraphAcceptedFlowStep>& flow_history,
		const std::vector<GraphAcceptedSpeciesStep>& species_history, const std::string& previous)
	{
		PhaseScope checkpoint_phase(ProfilePhase::Output);
		CoupledCheckpointEpoch epoch;
		CollectiveLocalStage(comm_, "native checkpoint epoch", [&] {
			epoch = {std::string(64, '0'), previous, identity_, static_cast<std::uint64_t>(step.step_index)+1, step.EndTime(), step.dt_s};
			epoch.id = coupled_checkpoint_detail::Hash(CheckpointEpochAgreementBytes(epoch));
			ValidatePressureFlowCheckpointControls(controls, graph_, epoch);
			checkpoint_metadata::Require((species_ ? species_history.size() : flow_history.size()) == epoch.accepted_steps, "checkpoint history prefix is incomplete");
			if (rank_ == 0) {
				std::vector<std::filesystem::path> created;
				auto parent = std::filesystem::absolute(root).lexically_normal();
				while (!std::filesystem::exists(parent)) { created.push_back(parent); parent = parent.parent_path(); }
				std::filesystem::create_directories(root);
				for (const auto& path : created) coupled_checkpoint_detail::Root(path).Sync();
				if (!created.empty()) coupled_checkpoint_detail::Root(parent).Sync();
				const auto initial = epoch.id;
				std::uint64_t attempt = 0;
				while (std::filesystem::exists(root/epoch.id)) {
					checkpoint_metadata::Require(++attempt <= 65536, "too many checkpoint epoch collisions");
					epoch.id = coupled_checkpoint_detail::Hash(initial+":"+std::to_string(attempt));
				}
				CreateCoupledCheckpointEpoch(root, epoch);
			}
		});
		MPI_Bcast(epoch.id.data(), 64, MPI_CHAR, 0, comm_);
		std::vector<CoupledCheckpointShard> receipts;
		auto metadata = [&](const std::string& id, const std::string& bytes) {
			std::string digest;
			CollectiveLocalStage(comm_, "native checkpoint metadata digest", [&] { digest = coupled_checkpoint_detail::Hash(bytes); });
			RequireCollectiveSameText(comm_, "native checkpoint replicated metadata", digest);
			CollectiveLocalStage(comm_, "native checkpoint metadata write", [&] {
				if (rank_ == 0) receipts.push_back(WriteCoupledCheckpointShard(root, epoch, Spec(id), bytes.size(), [&](auto& out) { out.Write(bytes.data(), bytes.size()); }));
			});
		};
		std::string bytes;
		CollectiveLocalStage(comm_, "native checkpoint controls", [&] { bytes = SerializePressureFlowCheckpointControls(controls); }); metadata("graph.controls", bytes);
		for (const auto& item : owners_.zero_d) {
			CollectiveLocalStage(comm_, "native checkpoint 0D capture", [&] { bytes = SerializeZeroDFlowCheckpoint(item.second->CaptureCheckpointState()); (void)ParseZeroDFlowCheckpoint(bytes, epoch); });
			metadata(ids_.at(item.first)[0], bytes);
		}
		for (const auto& item : owners_.one_d) {
			OneDFlowCheckpointState state; std::string digest;
			CollectiveLocalStage(comm_, "native checkpoint 1D capture", [&] {
				state = item.second->CaptureCheckpointState();
				checkpoint_metadata::RequireAcceptedStep({state.accepted_macro_steps-1, state.last_macro_start_s, state.last_macro_dt_s}, epoch.accepted_steps, epoch.time_s, epoch.dt_s);
				bytes = SerializeOneDAcceptedCheckpointMetadata(state);
				Sha256 hash; WriteOneDAcceptedCheckpointFields(state, [&](const void* p, std::size_t n) { hash.Append(p, n); }); digest = hash.Hex();
			});
			metadata(ids_.at(item.first)[1], bytes);
			RequireCollectiveSameText(comm_, "native checkpoint replicated 1D field", digest);
			CollectiveLocalStage(comm_, "native checkpoint 1D field write", [&] {
				if (rank_ == 0) receipts.push_back(WriteCoupledCheckpointShard(root, epoch, Spec(ids_.at(item.first)[2]), OneDAcceptedCheckpointFieldBytes(state), [&](auto& out) {
					WriteOneDAcceptedCheckpointFields(state, [&](const void* p, std::size_t n) { out.Write(p, n); });
				}));
			});
		}
		for (const auto& item : owners_.flow) {
			auto state = item.second->CaptureCheckpointState(); TransportAcceptedCheckpointState scalar;
			if (owners_.transport.count(item.first)) scalar = owners_.transport.at(item.first)->CaptureCheckpointState();
			CollectiveLocalStage(comm_, "native checkpoint 3D shard write", [&] {
				auto local = WriteBodyFittedCheckpointShards(root, epoch, layouts_.at(item.first), rank_, state, owners_.transport.count(item.first) ? &scalar : nullptr);
				receipts.insert(receipts.end(), local.begin(), local.end());
			});
		}
		auto history = [&](const auto& records) {
			std::uint64_t count = 0; std::string digest;
			CollectiveLocalStage(comm_, "native checkpoint history capture", [&] {
				count = CheckpointRecordBytes(records); Sha256 hash;
				WriteCheckpointRecords(records, [&](const void* p, std::size_t n) { hash.Append(p, n); }); digest = hash.Hex();
			});
			RequireCollectiveSameText(comm_, "native checkpoint history agreement", digest);
			CollectiveLocalStage(comm_, "native checkpoint history write", [&] {
				if (rank_ == 0) receipts.push_back(WriteCoupledCheckpointShard(root, epoch, Spec("graph.history"), count, [&](auto& out) {
					WriteCheckpointRecords(records, [&](const void* p, std::size_t n) { out.Write(p, n); });
				}));
			});
		};
		if (species_) history(species_history); else history(flow_history);
		CollectiveLocalStage(comm_, "native checkpoint receipt order", [&] { std::sort(receipts.begin(), receipts.end(), [](const auto& a, const auto& b) { return a.spec.id < b.spec.id; }); });
		const auto gathered = GatherCheckpointReceipts(comm_, epoch, receipts);
		CollectiveLocalStage(comm_, "native checkpoint publish", [&] { if (rank_ == 0) PublishCoupledCheckpoint(root, epoch, catalog_, gathered); });
		return epoch;
	}

	CoupledCheckpointEpoch Restore(const std::filesystem::path& root, int final_step, double dt,
		PressureFlowCheckpointControls& controls, std::vector<GraphAcceptedFlowStep>& flow_history,
		std::vector<GraphAcceptedSpeciesStep>& species_history, SpeciesPressureFlowComponentExecutor* executor)
	{
		PhaseScope checkpoint_phase(ProfilePhase::Output);
		std::string bytes;
		CollectiveLocalStage(comm_, "native checkpoint discovery", [&] {
			if (rank_ != 0) return;
			const auto discovered = FindLatestCoupledCheckpoint(root, identity_, catalog_);
			checkpoint_metadata::Require(bool(discovered.latest), "no compatible complete graph checkpoint"); bytes = SerializeCoupledCheckpointManifest(*discovered.latest);
		});
		std::uint64_t length = bytes.size(); MPI_Bcast(&length, 1, MPI_UINT64_T, 0, comm_);
		CollectiveLocalStage(comm_, "native checkpoint manifest allocation", [&] { checkpoint_metadata::Require(length <= coupled_checkpoint_detail::maximum_manifest_bytes, "graph manifest too large"); bytes.resize(length); });
		MPI_Bcast(bytes.data(), static_cast<int>(length), MPI_CHAR, 0, comm_);
		CoupledCheckpointManifest manifest;
		CollectiveLocalStage(comm_, "native checkpoint manifest", [&] {
			manifest = ParseCoupledCheckpointManifest(bytes); coupled_checkpoint_detail::RequireCompatible(manifest.epoch.compatibility, identity_);
			checkpoint_metadata::Require(manifest.epoch.accepted_steps <= static_cast<std::uint64_t>(final_step) && manifest.epoch.dt_s == dt, "checkpoint exceeds requested horizon or dt differs");
		});
		auto index = [&](const std::string& id) {
			const auto found = std::find_if(manifest.shards.begin(), manifest.shards.end(), [&](const auto& s) { return s.spec.id == id; });
			checkpoint_metadata::Require(found != manifest.shards.end(), "missing graph checkpoint shard"); return static_cast<std::size_t>(found-manifest.shards.begin());
		};
		auto metadata = [&](const std::string& id) {
			const auto i = index(id); checkpoint_metadata::Require(manifest.shards[i].payload_bytes <= checkpoint_metadata::maximum_bytes, "graph metadata too large");
			std::string text; ReadCoupledCheckpointShard(root, manifest, i, [&](const void* p, std::size_t n) { text.append(static_cast<const char*>(p), n); }); return text;
		};
		std::map<std::string, ZeroDFlowCheckpointState> zero;
		std::map<std::string, OneDFlowCheckpointState> one;
		std::map<std::string, FlowAcceptedCheckpointState> flow;
		std::map<std::string, TransportAcceptedCheckpointState> transport;
		// Allocate all PETSc decoder candidates collectively before any local I/O.
		for (const auto& item : owners_.flow) {
			auto value = item.second->CreateCheckpointRestoreCandidate();
			CollectiveLocalStage(comm_, "native checkpoint flow candidate", [&] { flow.emplace(item.first, std::move(value)); });
		}
		for (const auto& item : owners_.transport) {
			auto value = item.second->CreateCheckpointRestoreCandidate();
			CollectiveLocalStage(comm_, "native checkpoint transport candidate", [&] { transport.emplace(item.first, std::move(value)); });
		}
		CollectiveLocalStage(comm_, "native checkpoint candidate read", [&] {
			controls = ParsePressureFlowCheckpointControls(metadata("graph.controls"), graph_, manifest.epoch);
			for (const auto& item : owners_.zero_d) zero.emplace(item.first, ParseZeroDFlowCheckpoint(metadata(GraphCheckpointDomainPrefix(item.first)+".zero-d"), manifest.epoch));
			for (const auto& item : owners_.one_d) {
				const auto prefix = GraphCheckpointDomainPrefix(item.first);
				auto value = ParseOneDAcceptedCheckpointMetadata(metadata(prefix+".one-d.metadata"), *item.second, manifest.epoch);
				const auto i = index(prefix+".one-d.fields"); checkpoint_metadata::Require(manifest.shards[i].payload_bytes == OneDAcceptedCheckpointFieldBytes(value), "1D graph field size differs");
				OneDAcceptedCheckpointFieldsReader reader(value); ReadCoupledCheckpointShard(root, manifest, i, [&](const void* p, std::size_t n) { reader.Consume(p, n); }); reader.Finish(); one.emplace(item.first, std::move(value));
			}
			for (const auto& item : owners_.flow) LoadBodyFittedCheckpointShards(root, manifest, layouts_.at(item.first), rank_, flow.at(item.first),
				transport.count(item.first) ? &transport.at(item.first) : nullptr, dt);
			double time = 0.0;
			CheckpointRecordReader reader(manifest.epoch.accepted_steps, [&](std::uint64_t record, std::string_view image) {
				time += dt;
				if (species_) {
					auto value = ParseGraphSpeciesHistoryRecord(image); ValidateHistory(value, record, time, value.result.hydraulic_iterations, value.result.accepted_ports); ValidateSpeciesHistory(value); species_history.push_back(std::move(value));
				} else {
					auto value = ParseGraphFlowHistoryRecord(image); ValidateHistory(value, record, time, value.result.iterations, value.result.accepted_ports); ValidateFlowHistory(value); flow_history.push_back(std::move(value));
				}
			});
			ReadCoupledCheckpointShard(root, manifest, index("graph.history"), [&](const void* p, std::size_t n) { reader.Consume(p, n); }); reader.Finish();
			checkpoint_metadata::Require(time == manifest.epoch.time_s, "graph history end time differs");
			const auto& last = species_ ? species_history.back().result.hydraulic_iterations.back() : flow_history.back().result.iterations.back();
			const auto reconstructed = MakePressureFlowCheckpointControls(controls.accepted_step, last, species_ ? species_history.back().result.donor_ownership : decltype(controls.donors){});
			checkpoint_metadata::Require(SerializePressureFlowCheckpointControls(controls) == SerializePressureFlowCheckpointControls(reconstructed), "graph history and controls differ");
		});
		// This candidate graph has not executed a step or published output. A
		// failure below rejects the entire graph, even if an earlier owner restored.
		CollectiveLocalStage(comm_, "native checkpoint local owner restore", [&] {
			for (const auto& item : owners_.zero_d) item.second->RestoreCheckpointState(std::move(zero.at(item.first)));
			for (const auto& item : owners_.one_d) item.second->RestoreCheckpointState(std::move(one.at(item.first)));
			if (executor) executor->RestoreCheckpointDonors(controls.donors);
		});
		for (const auto& item : owners_.flow) item.second->RestoreCheckpointState(flow.at(item.first));
		for (const auto& item : owners_.transport) item.second->RestoreCheckpointState(transport.at(item.first));
		CollectiveLocalStage(comm_, "native checkpoint candidate agreement", [] {}); return manifest.epoch;
	}
private:
	const CoupledCheckpointShardSpec& Spec(const std::string& id) const
	{
		const auto found = std::find_if(catalog_.begin(), catalog_.end(), [&](const auto& s) { return s.id == id; });
		checkpoint_metadata::Require(found != catalog_.end(), "unknown graph shard"); return *found;
	}
	template<class Step> void ValidateHistory(const Step& value, std::uint64_t index, double time,
		const std::vector<PressureFlowIterationState>& iterations, const std::map<PortRef, PortState>& ports) const
	{
		checkpoint_metadata::Require(static_cast<std::uint64_t>(value.step) == index+1 && value.time_s == time, "history prefix clock differs");
		std::size_t expected_ports = 0;
		for (const auto& domain : graph_.Domains()) expected_ports += domain.second.ports.size();
		checkpoint_metadata::Require(ports.size() == expected_ports, "incomplete history port catalog");
		for (const auto& port : ports) { (void)graph_.Port(port.first); checkpoint_metadata::Require(bool(port.second.outward_flow_m3_s), "history port has no flow"); }
		checkpoint_metadata::Require(value.three_d_balance.size() == owners_.flow.size(), "history 3D balance catalog differs");
		for (const auto& item : owners_.flow) checkpoint_metadata::Require(value.three_d_balance.count(item.first), "missing 3D history balance");
		int number = 0;
		for (const auto& iteration : iterations) {
			checkpoint_metadata::Require(iteration.iteration == ++number, "history iteration order differs");
			std::set<std::string> ids;
			for (const auto& edge : iteration.edges) checkpoint_metadata::Require(ids.insert(edge.edge_id).second, "duplicate history edge");
			checkpoint_metadata::Require(ids.size() == graph_.Edges().size(), "incomplete history edge catalog");
			for (const auto& edge : graph_.Edges()) checkpoint_metadata::Require(ids.count(edge.id), "unknown history edge");
		}
	}
	void ValidateFlowHistory(const GraphAcceptedFlowStep& value) const
	{
		checkpoint_metadata::Require(value.zero_d_states.size() == owners_.zero_d.size() && value.zero_d_accounting.size() == owners_.zero_d.size(), "history 0D catalog differs");
		for (const auto& item : owners_.zero_d) {
			checkpoint_metadata::Require(value.zero_d_states.count(item.first) && value.zero_d_accounting.count(item.first), "missing 0D history");
			ValidateZeroDFlowState(value.zero_d_states.at(item.first)); ValidateZeroDFlowStepAccounting(value.zero_d_accounting.at(item.first));
		}
	}
	void ValidateSpeciesHistory(const GraphAcceptedSpeciesStep& value) const
	{
		const auto& result = value.result;
		checkpoint_metadata::Require(value.transport_ports.size() == result.accepted_ports.size(), "history transport port catalog differs");
		for (const auto& port : value.transport_ports) (void)graph_.Port(port.first);
		std::set<std::string> domains(result.transport_domain_order.begin(), result.transport_domain_order.end());
		checkpoint_metadata::Require(domains.size() == graph_.Domains().size() && domains.size() == result.transport_domain_order.size(), "history transport order differs");
		for (const auto& item : graph_.Domains()) checkpoint_metadata::Require(domains.count(item.first), "missing transport domain");
		std::set<std::pair<std::string, std::string>> expected, actual;
		for (const auto& edge : graph_.Edges()) for (const auto& species : edge.species) expected.emplace(edge.id, species);
		checkpoint_metadata::Require(result.donor_ownership.size() == expected.size(), "history donor catalog differs");
		for (const auto& key : expected) checkpoint_metadata::Require(result.donor_ownership.count(key), "missing history donor");
		for (const auto& amount : result.edge_amounts) {
			const auto key = std::make_pair(amount.edge_id, amount.species_id); const auto& edge = graph_.Edge(amount.edge_id);
			checkpoint_metadata::Require(actual.insert(key).second && expected.count(key), "history edge amount catalog differs");
			const bool swapped = edge.second < edge.first;
			auto donor = result.donor_ownership.at(key);
			if (swapped) donor = donor == SpeciesDonor::First ? SpeciesDonor::Second : SpeciesDonor::First;
			checkpoint_metadata::Require(amount.first == (swapped ? edge.second : edge.first)
				&& amount.second == (swapped ? edge.first : edge.second) && amount.donor == donor, "history canonical edge orientation differs");
		}
		checkpoint_metadata::Require(actual == expected, "incomplete history edge amounts");
		expected.clear(); actual.clear();
		for (const auto& domain : graph_.Domains()) for (const auto& binding : domain.second.species_bindings) expected.emplace(domain.first, binding.first);
		for (const auto& balance : result.domain_balances) {
			const auto key = std::make_pair(balance.domain_id, balance.species_id);
			checkpoint_metadata::Require(actual.insert(key).second && expected.count(key), "history domain balance catalog differs");
			const auto& domain = graph_.Domain(balance.domain_id);
			checkpoint_metadata::Require(balance.accounting.outward_port_amount.size() == domain.ports.size(), "history amount port catalog differs");
			for (const auto& port : domain.ports) checkpoint_metadata::Require(balance.accounting.outward_port_amount.count(port.id), "missing history amount port");
		}
		checkpoint_metadata::Require(actual == expected && result.global_balances.size() == graph_.Species().size(), "incomplete species balance catalog");
		for (const auto& species : graph_.Species()) {
			const auto found = result.global_balances.find(species.first);
			checkpoint_metadata::Require(found != result.global_balances.end() && found->second.species_id == species.first, "history global species differs");
		}
	}
	MPI_Comm comm_;
	const SimulationGraph& graph_;
	CoupledCheckpointCompatibility identity_;
	NativeGraphCheckpointOwners owners_;
	bool species_ = false;
	int rank_ = 0;
	std::vector<CoupledCheckpointShardSpec> catalog_;
	std::map<std::string, BodyFittedCheckpointLayout> layouts_;
	std::map<std::string, std::array<std::string, 3>> ids_;
};
} // namespace iga
#endif
