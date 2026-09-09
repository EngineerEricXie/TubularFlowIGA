#ifndef IGA_PETSC_CHECKPOINT_READ_HPP
#define IGA_PETSC_CHECKPOINT_READ_HPP

#include "CollectiveFailure.hpp"
#include "CollectiveAssetInput.hpp"

#include <petscvec.h>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <limits>
#include <string>
#include <vector>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

namespace iga {
namespace petsc_checkpoint_detail {

inline void Check(PetscErrorCode error, const char* operation)
{
	if (error) throw std::runtime_error(std::string(operation)+" returned PETSc error "+std::to_string(error));
}

// All file operations below are local. PetscBinaryRead performs PETSc's
// endian conversion without entering a viewer's communicator collectives.
class LocalBinaryInput {
public:
	explicit LocalBinaryInput(const std::string& path)
	{
		fd_ = ::open(path.c_str(), O_RDONLY | O_NONBLOCK);
		if (fd_ < 0) throw std::runtime_error("cannot open checkpoint state: "+std::string(std::strerror(errno)));
	}
	LocalBinaryInput(const LocalBinaryInput&) = delete;
	LocalBinaryInput& operator=(const LocalBinaryInput&) = delete;
	~LocalBinaryInput() { if (fd_ >= 0) ::close(fd_); }
	std::uint64_t Size() const
	{
		struct stat status{};
		if (::fstat(fd_, &status) || !S_ISREG(status.st_mode) || status.st_size < 0)
			throw std::runtime_error("checkpoint state must be a readable regular file");
		return static_cast<std::uint64_t>(status.st_size);
	}
	void Read(void* values, PetscInt count, PetscDataType type)
	{
		if (count) Check(PetscBinaryRead(fd_, values, count, nullptr, type), "checkpoint binary read");
	}
	void Seek(std::uint64_t offset)
	{
		if (offset > static_cast<std::uint64_t>(std::numeric_limits<off_t>::max()))
			throw std::runtime_error("checkpoint byte offset exceeds off_t");
		off_t actual = 0;
		Check(PetscBinarySeek(fd_, static_cast<off_t>(offset), PETSC_BINARY_SEEK_SET, &actual), "checkpoint binary seek");
		if (actual != static_cast<off_t>(offset)) throw std::runtime_error("checkpoint seek returned wrong offset");
	}
	void Close()
	{
		const auto fd = fd_; fd_ = -1;
		if (::close(fd)) throw std::runtime_error("cannot close checkpoint state");
	}
private:
	int fd_ = -1;
};

struct ReturnErrors {
	ReturnErrors() { Check(PetscPushErrorHandler(PetscReturnErrorHandler, nullptr), "checkpoint error handler"); }
	~ReturnErrors() { PetscPopErrorHandler(); }
	ReturnErrors(const ReturnErrors&) = delete;
	ReturnErrors& operator=(const ReturnErrors&) = delete;
};

struct CandidateVector {
	Vec value = nullptr;
	CandidateVector() = default;
	CandidateVector(const CandidateVector&) = delete;
	CandidateVector& operator=(const CandidateVector&) = delete;
	~CandidateVector() { if (value) VecDestroy(&value); }
};

} // namespace petsc_checkpoint_detail

// Read one ordinary, uncompressed PETSc binary Vec using the current build's
// PetscInt/PetscScalar ABI. Every rank reads only its owned contiguous rows.
namespace petsc_checkpoint_detail {

// Callers either require a shared path or have already verified equal contents
// for immutable local replicas. Keep candidate validation/publication identical.
inline void ReadVector(Vec target, MPI_Comm communicator,
	const std::filesystem::path& path, bool compare_path)
{
	PetscInt rows = 0, begin = 0, end = 0;
	std::size_t integer_bytes = 0, scalar_bytes = 0;
	std::uint64_t expected_bytes = 0, data_offset = 0;
	std::string filename, description;
	std::vector<PetscScalar> values;
	CollectiveLocalStage(communicator, "checkpoint read preparation", [&] {
		using petsc_checkpoint_detail::Check;
		filename = path.string();
		if (filename.empty() || filename.find('\0') != std::string::npos)
			throw std::runtime_error("checkpoint path must be nonempty and contain no NUL");
		filename = std::filesystem::absolute(path).lexically_normal().string();
		Check(VecGetSize(target, &rows), "checkpoint vector size");
		Check(VecGetOwnershipRange(target, &begin, &end), "checkpoint vector ownership");
		if (rows < 0 || begin < 0 || end < begin || end > rows)
			throw std::runtime_error("invalid checkpoint vector ownership");
		Check(PetscDataTypeGetSize(PETSC_INT, &integer_bytes), "checkpoint integer size");
		Check(PetscDataTypeGetSize(PETSC_SCALAR, &scalar_bytes), "checkpoint scalar size");
		if (!scalar_bytes || static_cast<std::uint64_t>(rows)
			> (std::numeric_limits<std::uint64_t>::max()-2*integer_bytes)/scalar_bytes)
			throw std::runtime_error("checkpoint byte count overflows");
		expected_bytes = 2*integer_bytes+static_cast<std::uint64_t>(rows)*scalar_bytes;
		data_offset = 2*integer_bytes+static_cast<std::uint64_t>(begin)*scalar_bytes;
		values.resize(static_cast<std::size_t>(end-begin));
		description = (compare_path ? filename+"\n" : "")
			+std::to_string(rows)+" "+std::to_string(integer_bytes)+" "+std::to_string(scalar_bytes);
	});
	RequireCollectiveSameText(communicator, "checkpoint read agreement", description);
	CollectiveLocalStage(communicator, "checkpoint local read", [&] {
		petsc_checkpoint_detail::ReturnErrors errors;
		petsc_checkpoint_detail::LocalBinaryInput input(filename);
		if (input.Size() != expected_bytes) throw std::runtime_error("checkpoint state is truncated or has trailing data");
		PetscInt header[2]{};
		input.Read(header, 2, PETSC_INT);
		if (header[0] != VEC_FILE_CLASSID || header[1] != rows)
			throw std::runtime_error("checkpoint Vec header does not match the target layout");
		input.Seek(data_offset);
		input.Read(values.data(), end-begin, PETSC_SCALAR);
		input.Close();
		for (const auto value : values)
			if (PetscIsInfOrNanScalar(value)) throw std::runtime_error("checkpoint state contains a nonfinite value");
	});
	std::vector<PetscInt> indices;
	CollectiveLocalStage(communicator, "checkpoint candidate preparation", [&] {
		indices.resize(values.size());
		for (PetscInt i = 0; i < end-begin; ++i) indices[static_cast<std::size_t>(i)] = begin+i;
	});
	petsc_checkpoint_detail::CandidateVector candidate;
	RequireCollectivePetscSuccess(communicator, "checkpoint candidate create", VecDuplicate(target, &candidate.value));
	CollectiveLocalStage(communicator, "checkpoint candidate insertion", [&] {
		petsc_checkpoint_detail::Check(VecSetValues(candidate.value, end-begin, indices.data(), values.data(), INSERT_VALUES),
			"checkpoint candidate VecSetValues");
	});
	RequireCollectivePetscSuccess(communicator, "checkpoint candidate assembly begin", VecAssemblyBegin(candidate.value));
	RequireCollectivePetscSuccess(communicator, "checkpoint candidate assembly end", VecAssemblyEnd(candidate.value));
	RequireCollectivePetscSuccess(communicator, "checkpoint publish vector", VecCopy(candidate.value, target));
	RequireCollectivePetscSuccess(communicator, "checkpoint candidate destroy", VecDestroy(&candidate.value));
}

} // namespace petsc_checkpoint_detail

// The path must name the same immutable regular file on the shared filesystem.
// Adjacent .info files do not configure the solver or change this fixed layout.
// Publication happens only after all local file/format/value checks succeed.
inline void ReadPetscCheckpointVector(Vec target, MPI_Comm communicator,
	const std::filesystem::path& path)
{
	petsc_checkpoint_detail::ReadVector(target, communicator, path, true);
}

// Rank-local replicas need identical bytes, not identical absolute paths.
// Revalidate content here rather than trusting an unchecked caller assertion.
inline void ReadReplicatedPetscCheckpointVector(Vec target, MPI_Comm communicator,
	const std::filesystem::path& path)
{
	AssetFileCatalog files;
	CollectiveLocalStage(communicator, "checkpoint replica catalog", [&] {
		files.emplace("checkpoint asset state", path);
	});
	RequireCollectiveAssetFiles(communicator, files);
	petsc_checkpoint_detail::ReadVector(target, communicator, path, false);
}

} // namespace iga

#endif
