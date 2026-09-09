#ifndef IGA_COLLECTIVE_FAILURE_HPP
#define IGA_COLLECTIVE_FAILURE_HPP

#include <petscsys.h>
#include "PhaseProfile.hpp"

#include <array>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <exception>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace iga {

// Every member of the borrowed communicator must call the same stage in the
// same order. The callback must contain LOCAL operations only: wrapping a
// collective cannot rescue peers already waiting inside that collective.
// This protocol coordinates caught errors while every MPI process is alive;
// it does not recover from process loss or errors inside MPI itself.
template <class Function>
void CollectiveLocalStage(MPI_Comm communicator, const char* stage, Function&& function)
{
	int rank = 0, ranks = 1;
	MPI_Comm_rank(communicator, &rank);
	MPI_Comm_size(communicator, &ranks);
	std::array<char, 1024> message{};
	int failed_rank = ranks;
	try {
		std::forward<Function>(function)();
	} catch (const std::exception& error) {
		failed_rank = rank;
		// Bounded stack storage also handles a caught allocation failure without
		// allocating another string before agreement. Long diagnostics truncate.
		std::snprintf(message.data(), message.size(), "%s: rank %d: %s",
			stage, rank, error.what());
	} catch (...) {
		failed_rank = rank;
		std::snprintf(message.data(), message.size(), "%s: rank %d: unknown exception",
			stage, rank);
	}
	int first_failed_rank = ranks;
	PhaseScope communication_phase(ProfilePhase::Communication);
	MPI_Allreduce(&failed_rank, &first_failed_rank, 1, MPI_INT, MPI_MIN, communicator);
	if (first_failed_rank == ranks) return;
	MPI_Bcast(message.data(), static_cast<int>(message.size()), MPI_CHAR,
		first_failed_rank, communicator);
	throw std::runtime_error(message.data());
}

// The PETSc operation must finish before entering this helper. Its status is
// coordinated here; this cannot recover a peer stuck inside the operation.
inline void RequireCollectivePetscSuccess(MPI_Comm communicator, const char* stage,
	PetscErrorCode status)
{
	CollectiveLocalStage(communicator, stage, [&] {
		if (status) throw std::runtime_error("PETSc returned error "+std::to_string(status));
	});
}

inline void RequireCollectiveSameInt(MPI_Comm communicator, const char* stage, int value)
{
	int minimum = 0, maximum = 0;
	{
		PhaseScope communication_phase(ProfilePhase::Communication);
		MPI_Allreduce(&value, &minimum, 1, MPI_INT, MPI_MIN, communicator);
		MPI_Allreduce(&value, &maximum, 1, MPI_INT, MPI_MAX, communicator);
	}
	CollectiveLocalStage(communicator, stage, [&] {
		if (minimum != maximum) throw std::runtime_error("value differs within communicator");
	});
}

// Exact agreement on a previously captured input, without allocating a copy
// of rank 0's text on every process. All ranks complete the same broadcasts,
// including when lengths differ. Local paths need not match between replicas.
inline void RequireCollectiveSameText(MPI_Comm communicator, const char* stage,
	std::string_view local)
{
	int rank = 0;
	MPI_Comm_rank(communicator, &rank);
	std::uint64_t root_size = local.size();
	bool same = true;
	{
		PhaseScope communication_phase(ProfilePhase::Communication);
		MPI_Bcast(&root_size, 1, MPI_UINT64_T, 0, communicator);
		same = root_size == local.size();
		std::array<char, 4096> buffer{};
		for (std::uint64_t offset = 0; offset < root_size;) {
			const auto count = static_cast<std::size_t>(std::min<std::uint64_t>(
				buffer.size(), root_size-offset));
			if (rank == 0) std::memcpy(buffer.data(), local.data()+offset, count);
			MPI_Bcast(buffer.data(), static_cast<int>(count), MPI_CHAR, 0, communicator);
			if (same && std::memcmp(buffer.data(), local.data()+offset, count) != 0)
				same = false;
			offset += count;
		}
	}
	CollectiveLocalStage(communicator, stage, [&] {
		if (!same) throw std::runtime_error("input differs from communicator rank 0");
	});
}

} // namespace iga

#endif
