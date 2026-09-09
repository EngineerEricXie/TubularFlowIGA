#ifndef IGA_BODY_FITTED_CHECKPOINT_BUNDLE_HPP
#define IGA_BODY_FITTED_CHECKPOINT_BUNDLE_HPP

#include "BodyFittedAcceptedCheckpoint.hpp"
#include "CoupledCheckpointBundle.hpp"

namespace iga {

struct BodyFittedCheckpointLayout {
	std::string domain_id;
	// Entry i is the bundle/world rank corresponding to solver local rank i.
	std::vector<std::uint32_t> world_ranks;
	bool transport = false;
};

inline std::vector<CoupledCheckpointShardSpec> BodyFittedCheckpointCatalog(
	const BodyFittedCheckpointLayout& layout, std::uint32_t world_size)
{
	using namespace coupled_checkpoint_detail;
	Require(Identifier(layout.domain_id) && !layout.world_ranks.empty(), "invalid body-fitted checkpoint layout");
	std::set<std::uint32_t> unique;
	for (const auto rank : layout.world_ranks) Require(rank < world_size && unique.insert(rank).second, "invalid body-fitted rank mapping");
	const auto shared_owner = layout.world_ranks.front(); const auto prefix = layout.domain_id+".";
	std::vector<CoupledCheckpointShardSpec> catalog{{prefix+"flow.metadata", "body-flow-v1", shared_owner},
		{prefix+"flow.boundaries", "body-boundaries-v1", shared_owner}};
	if (layout.transport) catalog.push_back({prefix+"transport.metadata", "body-transport-v1", shared_owner});
	for (const auto rank : layout.world_ranks) {
		catalog.push_back({prefix+"flow.rank-"+std::to_string(rank), "owned-real-field-v1", rank});
		if (layout.transport) catalog.push_back({prefix+"transport.rank-"+std::to_string(rank), "owned-real-field-v1", rank});
	}
	std::sort(catalog.begin(), catalog.end(), [](const auto& a, const auto& b) { return a.id < b.id; });
	ValidateCatalog(catalog, world_size); return catalog;
}

// LOCAL producer/consumer operations. The caller surrounds the entire local
// operation with group failure agreement; never call these from a worker that
// could enter MPI while another group member is still performing local I/O.
inline std::vector<CoupledCheckpointShard> WriteBodyFittedCheckpointShards(const std::filesystem::path& root,
	const CoupledCheckpointEpoch& epoch, const BodyFittedCheckpointLayout& layout, std::size_t local_rank,
	const FlowAcceptedCheckpointState& flow, const TransportAcceptedCheckpointState* transport = nullptr)
{
	using namespace coupled_checkpoint_detail;
	Validate(epoch); const auto catalog = BodyFittedCheckpointCatalog(layout, epoch.compatibility.ranks);
	Require(local_rank < layout.world_ranks.size() && layout.transport == (transport != nullptr), "body-fitted producer layout differs");
	body_fitted_checkpoint_detail::Epoch(flow.accepted_steps, flow.accepted_time_s, flow.macro_dt_s, epoch);
	if (transport) Require(transport->accepted_steps == flow.accepted_steps, "flow and transport accepted counts differ");
	const auto owner = layout.world_ranks[local_rank]; std::vector<CoupledCheckpointShard> receipts;
	const auto prefix = layout.domain_id+".";
	for (const auto& spec : catalog) {
		if (spec.owner_rank != owner) continue;
		if (spec.id == prefix+"flow.metadata" || spec.id == prefix+"transport.metadata") {
			const auto bytes = spec.id == prefix+"flow.metadata" ? SerializeBodyFittedFlowMetadata(flow)
				: SerializeBodyFittedTransportMetadata(*transport, epoch.dt_s);
			receipts.push_back(WriteCoupledCheckpointShard(root, epoch, spec, bytes.size(), [&](auto& output) { output.Write(bytes.data(), bytes.size()); }));
		} else if (spec.id == prefix+"flow.boundaries") {
			receipts.push_back(WriteCoupledCheckpointShard(root, epoch, spec, BodyFittedBoundaryBytes(flow), [&](auto& output) {
				WriteBodyFittedBoundaries(flow, [&](const void* bytes, std::size_t size) { output.Write(bytes, size); });
			}));
		} else {
			const auto& field = spec.id == prefix+"flow.rank-"+std::to_string(owner) ? flow.field : transport->field;
			receipts.push_back(WriteCoupledCheckpointShard(root, epoch, spec, OwnedCheckpointFieldBytes(field), [&](auto& output) {
				WriteOwnedCheckpointField(field, [&](const void* bytes, std::size_t size) { output.Write(bytes, size); });
			}));
		}
	}
	return receipts;
}

inline void LoadBodyFittedCheckpointShards(const std::filesystem::path& root, const CoupledCheckpointManifest& manifest,
	const BodyFittedCheckpointLayout& layout, std::size_t local_rank, FlowAcceptedCheckpointState& flow,
	TransportAcceptedCheckpointState* transport = nullptr, double transport_dt = 0.0)
{
	using namespace coupled_checkpoint_detail;
	Validate(manifest.epoch); const auto catalog = BodyFittedCheckpointCatalog(layout, manifest.epoch.compatibility.ranks);
	Require(local_rank < layout.world_ranks.size() && layout.transport == (transport != nullptr), "body-fitted consumer layout differs");
	std::map<std::string, std::size_t> indices;
	for (const auto& spec : catalog) {
		const auto found = std::find_if(manifest.shards.begin(), manifest.shards.end(), [&](const auto& shard) { return shard.spec.id == spec.id; });
		Require(found != manifest.shards.end() && Same(found->spec, spec), "body-fitted manifest shard catalog differs");
		indices.emplace(spec.id, static_cast<std::size_t>(found-manifest.shards.begin()));
	}
	const auto prefix = layout.domain_id+".";
	auto metadata = [&](const std::string& id) {
		const auto index = indices.at(prefix+id);
		Require(manifest.shards[index].payload_bytes <= checkpoint_metadata::maximum_bytes, "body-fitted metadata exceeds limit");
		std::string bytes;
		ReadCoupledCheckpointShard(root, manifest, index, [&](const void* data, std::size_t size) {
			Require(size <= checkpoint_metadata::maximum_bytes-bytes.size(), "body-fitted metadata exceeds limit");
			bytes.append(static_cast<const char*>(data), size);
		}); return bytes;
	};
	ParseBodyFittedFlowMetadata(metadata("flow.metadata"), flow, manifest.epoch);
	if (transport) ParseBodyFittedTransportMetadata(metadata("transport.metadata"), *transport, transport_dt, manifest.epoch);
	const auto boundary = indices.at(prefix+"flow.boundaries");
	Require(manifest.shards[boundary].payload_bytes == BodyFittedBoundaryBytes(flow), "boundary payload size differs");
	auto boundaries = ReadBodyFittedBoundaries(flow);
	ReadCoupledCheckpointShard(root, manifest, boundary, [&](const void* data, std::size_t size) { boundaries.Consume(data, size); });
	boundaries.Finish();
	auto field = [&](const std::string& name, OwnedCheckpointVector& value) {
		const auto index = indices.at(prefix+name+".rank-"+std::to_string(layout.world_ranks[local_rank]));
		Require(manifest.shards[index].payload_bytes == OwnedCheckpointFieldBytes(value), "owned field payload size differs");
		auto reader = ReadOwnedCheckpointField(value);
		ReadCoupledCheckpointShard(root, manifest, index, [&](const void* data, std::size_t size) { reader.Consume(data, size); });
		reader.Finish();
	};
	field("flow", flow.field); if (transport) field("transport", transport->field);
}

} // namespace iga
#endif
