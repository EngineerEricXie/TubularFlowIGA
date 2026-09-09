#ifndef IGA_OWNED_CHECKPOINT_VECTOR_HPP
#define IGA_OWNED_CHECKPOINT_VECTOR_HPP

#include "PetscReadArray.hpp"
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

struct OwnedCheckpointVector {
	std::uint64_t global_rows = 0, row_begin = 0, row_end = 0;
	std::vector<double> values;
};

// LOCAL operations only. Call them in a CollectiveLocalStage before any
// subsequent collective. The runtime owns the communicator agreement protocol.
inline void ValidateOwnedCheckpointVector(Vec target, const OwnedCheckpointVector& value)
{
	if (!target) throw std::runtime_error("checkpoint target vector is closed");
	PetscInt global = 0, begin = 0, end = 0, local = 0;
	if (VecGetSize(target, &global) || VecGetOwnershipRange(target, &begin, &end) || VecGetLocalSize(target, &local))
		throw std::runtime_error("cannot query checkpoint vector ownership");
	if (global < 0 || begin < 0 || end < begin || local != end-begin
		|| value.global_rows != static_cast<std::uint64_t>(global)
		|| value.row_begin != static_cast<std::uint64_t>(begin)
		|| value.row_end != static_cast<std::uint64_t>(end)
		|| value.values.size() != static_cast<std::uint64_t>(local))
		throw std::runtime_error("checkpoint vector ownership or field shape differs");
	for (double scalar : value.values)
		if (!std::isfinite(scalar)) throw std::runtime_error("nonfinite checkpoint vector value");
}

inline OwnedCheckpointVector CaptureOwnedCheckpointVector(Vec source)
{
	if (!source) throw std::runtime_error("checkpoint source vector is closed");
	PetscInt global = 0, begin = 0, end = 0;
	if (VecGetSize(source, &global) || VecGetOwnershipRange(source, &begin, &end)
		|| global < 0 || begin < 0 || end < begin)
		throw std::runtime_error("cannot query checkpoint source ownership");
	OwnedCheckpointVector value; value.global_rows = static_cast<std::uint64_t>(global);
	value.row_begin = static_cast<std::uint64_t>(begin); value.row_end = static_cast<std::uint64_t>(end);
	value.values.resize(static_cast<std::size_t>(end-begin));
	PetscReadArray view; view.Acquire(source);
	for (std::size_t i = 0; i < value.values.size(); ++i) {
		if (PetscImaginaryPart(view.Data()[i]) != 0.0)
			throw std::runtime_error("checkpoint requires a real-valued field");
		value.values[i] = PetscRealPart(view.Data()[i]);
	}
	view.Restore(); ValidateOwnedCheckpointVector(source, value); return value;
}

inline void StageOwnedCheckpointVector(Vec scratch, const OwnedCheckpointVector& value)
{
	ValidateOwnedCheckpointVector(scratch, value);
	PetscScalar* values = nullptr;
	if (VecGetArray(scratch, &values)) throw std::runtime_error("cannot acquire checkpoint staging array");
	for (std::size_t i = 0; i < value.values.size(); ++i) values[i] = value.values[i];
	if (VecRestoreArray(scratch, &values)) throw std::runtime_error("cannot restore checkpoint staging array");
}

inline void ValidateCheckpointConfigurationIdentity(const std::string& identity)
{
	if (identity.size() != 64 || identity.find_first_not_of("0123456789abcdef") != std::string::npos)
		throw std::runtime_error("checkpoint configuration identity must be lowercase SHA-256");
}

} // namespace iga
#endif
