#ifndef IGA_PARALLEL_OWNERSHIP_VALIDATION_HPP
#define IGA_PARALLEL_OWNERSHIP_VALIDATION_HPP

#include "ParallelOwnership.hpp"
#include <petscsys.h>

#include <algorithm>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace iga {

// All members must call each validation stage. This synchronizes controlled
// validation failures, not process loss or arbitrary exceptions in a solver.
inline void RequireOwnershipCondition(bool local_valid, MPI_Comm communicator,
	const std::string& message)
{
	int local = local_valid ? 1 : 0, global = 0;
	MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, communicator);
	if (!global) throw std::runtime_error("parallel ownership: "+message);
}

inline void ValidateDistributedNodeRanges(const OwnedNodeRange& local,
	MPI_Comm communicator)
{
	int ranks = 1;
	MPI_Comm_size(communicator, &ranks);
	const std::uint64_t values[3] = {local.begin, local.end, local.global_nodes};
	std::vector<std::uint64_t> ranges(static_cast<std::size_t>(ranks)*3);
	MPI_Allgather(values, 3, MPI_UINT64_T, ranges.data(), 3, MPI_UINT64_T, communicator);
	bool valid = true;
	std::uint64_t end = 0;
	for (int rank = 0; rank < ranks; ++rank) {
		const auto offset = static_cast<std::size_t>(rank)*3;
		valid = valid && ranges[offset] == end && ranges[offset] <= ranges[offset+1]
			&& ranges[offset+1] <= local.global_nodes && ranges[offset+2] == local.global_nodes;
		end = ranges[offset+1];
	}
	RequireOwnershipCondition(valid && end == local.global_nodes, communicator,
		"node ranges must cover the global range once in communicator rank order");
}

namespace ownership_detail {

// Fixed-width records are routed by their first ID. Width is a call-site
// constant shared by every group member (one for IDs, three for incidences).
inline std::vector<std::uint64_t> RedistributeRecords(const std::vector<std::uint64_t>& records,
	std::size_t width, MPI_Comm communicator, const std::string& context)
{
	int ranks = 1;
	MPI_Comm_size(communicator, &ranks);
	const auto limit = static_cast<std::size_t>(std::numeric_limits<int>::max());
	RequireOwnershipCondition(width > 0 && records.size()%width == 0 && records.size() <= limit,
		communicator, context+": invalid record or MPI count range");
	std::vector<int> sent(ranks, 0), received(ranks, 0);
	for (std::size_t i = 0; i < records.size(); i += width)
		sent[static_cast<std::size_t>(records[i]%ranks)] += static_cast<int>(width);
	MPI_Alltoall(sent.data(), 1, MPI_INT, received.data(), 1, MPI_INT, communicator);
	std::vector<int> send_offsets(ranks, 0), receive_offsets(ranks, 0);
	std::size_t send_total = 0, receive_total = 0;
	bool valid = true;
	for (int peer = 0; peer < ranks; ++peer) {
		if (send_total > limit || receive_total > limit || received[peer] < 0) {
			valid = false;
			break;
		}
		send_offsets[peer] = static_cast<int>(send_total);
		receive_offsets[peer] = static_cast<int>(receive_total);
		send_total += static_cast<std::size_t>(sent[peer]);
		receive_total += static_cast<std::size_t>(received[peer]);
	}
	RequireOwnershipCondition(valid && send_total <= limit && receive_total <= limit && receive_total%width == 0,
		communicator, context+": audit shard exceeds MPI count range");
	std::vector<std::uint64_t> outgoing(send_total), incoming(receive_total);
	auto positions = send_offsets;
	for (std::size_t i = 0; i < records.size(); i += width) {
		auto& position = positions[records[i]%ranks];
		for (std::size_t column = 0; column < width; ++column)
			outgoing[static_cast<std::size_t>(position++)] = records[i+column];
	}
	MPI_Alltoallv(outgoing.data(), sent.data(), send_offsets.data(), MPI_UINT64_T,
		incoming.data(), received.data(), receive_offsets.data(), MPI_UINT64_T, communicator);
	return incoming;
}

} // namespace ownership_detail

// Check exact coverage of [0, global_count), not merely a sum or hash of IDs.
// Each rank stores only its incoming shard plus local selected IDs; there is
// no full-catalog gather.
inline void ValidateUniqueEntityCoverage(const std::vector<std::uint64_t>& local_ids,
	std::uint64_t global_count, MPI_Comm communicator, const std::string& context)
{
	int rank = 0, ranks = 1;
	MPI_Comm_rank(communicator, &rank);
	MPI_Comm_size(communicator, &ranks);
	std::uint64_t minimum = 0, maximum = 0;
	MPI_Allreduce(&global_count, &minimum, 1, MPI_UINT64_T, MPI_MIN, communicator);
	MPI_Allreduce(&global_count, &maximum, 1, MPI_UINT64_T, MPI_MAX, communicator);
	bool valid = minimum == maximum;
	for (auto id : local_ids) valid = valid && id < global_count;
	RequireOwnershipCondition(valid, communicator, context+": invalid ID or global count");
	auto incoming = ownership_detail::RedistributeRecords(local_ids, 1, communicator, context);
	std::sort(incoming.begin(), incoming.end());
	const auto expected = global_count/static_cast<std::uint64_t>(ranks)
		+(static_cast<std::uint64_t>(rank) < global_count%ranks ? 1 : 0);
	valid = incoming.size() == expected;
	for (std::size_t i = 0; valid && i < incoming.size(); ++i)
		valid = incoming[i] == static_cast<std::uint64_t>(rank)
			+static_cast<std::uint64_t>(i)*static_cast<std::uint64_t>(ranks);
	RequireOwnershipCondition(valid, communicator, context+": missing or duplicate physical entity");
}

using OwnershipIncidence = std::pair<std::uint64_t, std::uint64_t>;

// The authoritative catalog and the consumed catalog may use different
// owners. Compare each (entity, node) pair exactly, including multiplicity.
inline void ValidateOwnershipIncidences(const std::vector<OwnershipIncidence>& authoritative,
	const std::vector<OwnershipIncidence>& consumed, MPI_Comm communicator, const std::string& context)
{
	const auto limit = static_cast<std::size_t>(std::numeric_limits<int>::max())/3;
	RequireOwnershipCondition(authoritative.size() <= limit && consumed.size() <= limit
		&& authoritative.size() <= limit-consumed.size(), communicator, context+": too many incidences");
	std::vector<std::uint64_t> records;
	records.reserve(3*(authoritative.size()+consumed.size()));
	for (const auto& pair : authoritative) records.insert(records.end(), {pair.first, pair.second, 0});
	for (const auto& pair : consumed) records.insert(records.end(), {pair.first, pair.second, 1});
	const auto incoming = ownership_detail::RedistributeRecords(records, 3, communicator, context);
	std::vector<OwnershipIncidence> expected, actual;
	for (std::size_t i = 0; i < incoming.size(); i += 3)
		(incoming[i+2] == 0 ? expected : actual).emplace_back(incoming[i], incoming[i+1]);
	std::sort(expected.begin(), expected.end());
	std::sort(actual.begin(), actual.end());
	RequireOwnershipCondition(expected == actual
		&& std::adjacent_find(expected.begin(), expected.end()) == expected.end(),
		communicator, context+": missing, duplicate or incorrect entity/node incidence");
}

} // namespace iga

#endif
