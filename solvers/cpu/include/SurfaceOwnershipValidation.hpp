#ifndef IGA_SURFACE_OWNERSHIP_VALIDATION_HPP
#define IGA_SURFACE_OWNERSHIP_VALIDATION_HPP

#include "DistributedSurfaceInterface.hpp"
#include "ParallelOwnershipValidation.hpp"

namespace iga {

// Validate publication ownership on an existing surface communicator. This
// adapts the replicated-reference layout contract; it does not distribute a
// membrane solve, surface quadrature, force transfer or reference geometry.
inline void ValidateSurfacePublicationOwnership(const SurfaceInterfaceRef& reference,
	const DistributedSurfaceLayout& layout, MPI_Comm communicator)
{
	int rank = 0, ranks = 1;
	MPI_Comm_rank(communicator, &rank);
	MPI_Comm_size(communicator, &ranks);
	bool valid = true;
	std::string identity;
	try {
		ValidateSurfaceInterfaceRef(reference);
		ValidateDistributedSurfaceLayout(layout);
		valid = layout.partition_count == static_cast<std::uint64_t>(ranks)
			&& layout.partition_rank == static_cast<std::uint64_t>(rank);
		Sha256 hash;
		for (const auto* value : {&reference.domain_id, &reference.subsystem_id,
			&reference.interface_id, &layout.layout_identity_sha256}) {
			hash.AppendLittleEndian64(static_cast<std::uint64_t>(value->size()));
			hash.Append(value->data(), value->size());
		}
		identity = hash.Hex();
	} catch (const std::exception&) {
		valid = false;
	}
	RequireOwnershipCondition(valid, communicator, "invalid surface layout or communicator partition");
	std::array<char, 64> root_identity{};
	if (rank == 0) std::copy(identity.begin(), identity.end(), root_identity.begin());
	MPI_Bcast(root_identity.data(), static_cast<int>(root_identity.size()), MPI_CHAR, 0, communicator);
	RequireOwnershipCondition(std::equal(identity.begin(), identity.end(), root_identity.begin()),
		communicator, "surface interface or reference geometry differs across ranks");
	std::vector<OwnershipIncidence> expected, published;
	// The current layout already replicates reference positions. Rank zero
	// contributes that authoritative ID catalog once to the sharded audit.
	if (rank == 0)
		for (const auto& position : layout.reference_positions)
			expected.emplace_back(position.global_node_id, 0);
	for (auto id : layout.owned_global_node_ids) published.emplace_back(id, 0);
	ValidateOwnershipIncidences(expected, published, communicator, "surface publication");
}

} // namespace iga

#endif
