#ifndef IGA_COLLECTIVE_CHECKPOINT_RECEIPTS_HPP
#define IGA_COLLECTIVE_CHECKPOINT_RECEIPTS_HPP

#include "CollectiveFailure.hpp"
#include "CheckpointMetadataCodec.hpp"
#include "CoupledCheckpointManifest.hpp"

namespace iga {
inline std::string CheckpointEpochAgreementBytes(const CoupledCheckpointEpoch& epoch)
{
	coupled_checkpoint_detail::Validate(epoch); checkpoint_metadata::Writer output;
	output.Text(epoch.id); output.Text(epoch.previous_id); output.Text(epoch.compatibility.case_sha256);
	output.Text(epoch.compatibility.configuration_sha256); output.Text(epoch.compatibility.execution_sha256);
	output.Unsigned(epoch.compatibility.ranks); output.Unsigned(epoch.accepted_steps); output.Real(epoch.time_s); output.Real(epoch.dt_s);
	return output.Bytes();
}

// Bounded receipt exchange only; field bytes never pass through MPI here.
// Returns the sorted union on communicator rank 0, empty elsewhere. Empty
// producers are allowed so a later graph coordinator can include idle ranks.
inline std::vector<CoupledCheckpointShard> GatherCheckpointReceipts(MPI_Comm communicator,
	const CoupledCheckpointEpoch& epoch, const std::vector<CoupledCheckpointShard>& local_receipts)
{
	int rank = 0, ranks = 1; MPI_Comm_rank(communicator, &rank); MPI_Comm_size(communicator, &ranks);
	std::string epoch_bytes, local_bytes;
	CollectiveLocalStage(communicator, "checkpoint receipt preparation", [&] {
		epoch_bytes = CheckpointEpochAgreementBytes(epoch);
		if (!local_receipts.empty()) local_bytes = SerializeCoupledCheckpointManifest({epoch, local_receipts});
	});
	RequireCollectiveSameText(communicator, "checkpoint receipt epoch agreement", epoch_bytes);
	std::vector<CoupledCheckpointShard> receipts;
	for (int source = 0; source < ranks; ++source) {
		std::uint64_t bytes = rank == source ? local_bytes.size() : 0;
		MPI_Bcast(&bytes, 1, MPI_UINT64_T, source, communicator);
		std::string buffer;
		CollectiveLocalStage(communicator, "checkpoint receipt receive preparation", [&] {
			coupled_checkpoint_detail::Require(bytes <= coupled_checkpoint_detail::maximum_manifest_bytes, "receipt payload exceeds limit");
			if (rank == source) buffer = local_bytes; else buffer.resize(static_cast<std::size_t>(bytes));
		});
		MPI_Bcast(buffer.data(), static_cast<int>(bytes), MPI_CHAR, source, communicator);
		CollectiveLocalStage(communicator, "checkpoint receipt merge", [&] {
			if (rank != 0 || !bytes) return;
			auto partial = ParseCoupledCheckpointManifest(buffer);
			coupled_checkpoint_detail::Require(CheckpointEpochAgreementBytes(partial.epoch) == epoch_bytes, "receipt epoch differs");
			coupled_checkpoint_detail::Require(partial.shards.size() <= coupled_checkpoint_detail::maximum_shards-receipts.size(), "too many combined receipts");
			receipts.insert(receipts.end(), std::make_move_iterator(partial.shards.begin()), std::make_move_iterator(partial.shards.end()));
		});
	}
	CollectiveLocalStage(communicator, "checkpoint receipt catalog validation", [&] {
		if (rank != 0) return;
		std::sort(receipts.begin(), receipts.end(), [](const auto& a, const auto& b) { return a.spec.id < b.spec.id; });
		if (!receipts.empty()) SerializeCoupledCheckpointManifest({epoch, receipts});
	});
	return receipts;
}
} // namespace iga
#endif
