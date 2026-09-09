#ifndef IGA_PETSC_CHECKPOINT_WRITE_HPP
#define IGA_PETSC_CHECKPOINT_WRITE_HPP

#include "PetscCheckpointRead.hpp"
#include "PetscReadArray.hpp"
#include <climits>
#include <cstdlib>

namespace iga {
namespace petsc_checkpoint_detail {

class LocalBinaryOutput {
public:
	explicit LocalBinaryOutput(int fd) : fd_(fd) {}
	LocalBinaryOutput(const LocalBinaryOutput&) = delete;
	LocalBinaryOutput& operator=(const LocalBinaryOutput&) = delete;
	~LocalBinaryOutput() { if (fd_ >= 0) ::close(fd_); }
	void Resize(std::uint64_t bytes)
	{
		if (::ftruncate(fd_, static_cast<off_t>(bytes))) throw std::runtime_error("cannot size checkpoint temporary file");
	}
	void Seek(std::uint64_t offset)
	{
		off_t actual = 0;
		Check(PetscBinarySeek(fd_, static_cast<off_t>(offset), PETSC_BINARY_SEEK_SET, &actual), "checkpoint write seek");
		if (actual != static_cast<off_t>(offset)) throw std::runtime_error("checkpoint write seek returned wrong offset");
	}
	void Write(void* values, PetscInt count, PetscDataType type)
	{
		// PETSc may byte-swap this buffer in place, including on an error path.
		if (count) Check(PetscBinaryWrite(fd_, values, count, type), "checkpoint binary write");
	}
	void Close()
	{
		if (::fsync(fd_)) throw std::runtime_error("cannot sync checkpoint temporary file");
		const auto fd = fd_; fd_ = -1;
		if (::close(fd)) throw std::runtime_error("cannot close checkpoint temporary file");
	}
private:
	int fd_ = -1;
};

class TemporaryCheckpoint {
public:
	TemporaryCheckpoint() = default;
	TemporaryCheckpoint(const TemporaryCheckpoint&) = delete;
	TemporaryCheckpoint& operator=(const TemporaryCheckpoint&) = delete;
	~TemporaryCheckpoint() { if (created_) ::unlink(path_.c_str()); }
	void Create(const std::string& target, std::uint64_t bytes, PetscInt rows)
	{
		struct stat previous{};
		const auto status = ::lstat(target.c_str(), &previous);
		if (status == 0 && !S_ISREG(previous.st_mode))
			throw std::runtime_error("checkpoint destination must be a regular file or absent");
		if (status != 0 && errno != ENOENT) throw std::runtime_error("cannot inspect checkpoint destination");
		path_ = target+".tmp.XXXXXX";
		const auto fd = ::mkstemp(path_.data());
		if (fd < 0) throw std::runtime_error("cannot create checkpoint temporary file: "+std::string(std::strerror(errno)));
		created_ = true;
		LocalBinaryOutput output(fd);
		if (status == 0 && ::fchmod(fd, previous.st_mode & 0777))
			throw std::runtime_error("cannot preserve checkpoint permissions");
		output.Resize(bytes);
		PetscInt header[2] = {VEC_FILE_CLASSID, rows};
		output.Write(header, 2, PETSC_INT);
		output.Close();
	}
	const std::string& Path() const { return path_; }
	void Publish(const std::string& target)
	{
		if (::rename(path_.c_str(), target.c_str())) throw std::runtime_error("cannot publish checkpoint state");
		created_ = false;
	}
private:
	std::string path_;
	bool created_ = false;
};

} // namespace petsc_checkpoint_detail

// Collective over the source vector's group. Requires shared storage and an
// exclusive writer for target. Local owned-row I/O completes before root renames
// the sibling temporary file. This publishes one Vec, not a coupled checkpoint.
inline void WritePetscCheckpointVector(Vec source, MPI_Comm communicator,
	const std::filesystem::path& path)
{
	PetscInt rows = 0, begin = 0, end = 0;
	std::size_t integer_bytes = 0, scalar_bytes = 0;
	std::uint64_t bytes = 0, offset = 0;
	std::string filename, description;
	std::vector<PetscScalar> values;
	CollectiveLocalStage(communicator, "checkpoint write preparation", [&] {
		using petsc_checkpoint_detail::Check;
		filename = path.string();
		if (filename.empty() || filename.find('\0') != std::string::npos)
			throw std::runtime_error("checkpoint path must be nonempty and contain no NUL");
		filename = std::filesystem::absolute(path).lexically_normal().string();
		Check(VecGetSize(source, &rows), "checkpoint source size");
		Check(VecGetOwnershipRange(source, &begin, &end), "checkpoint source ownership");
		Check(PetscDataTypeGetSize(PETSC_INT, &integer_bytes), "checkpoint integer size");
		Check(PetscDataTypeGetSize(PETSC_SCALAR, &scalar_bytes), "checkpoint scalar size");
		if (rows < 0 || begin < 0 || end < begin || end > rows || !scalar_bytes
			|| static_cast<std::uint64_t>(rows) > (std::numeric_limits<std::uint64_t>::max()-2*integer_bytes)/scalar_bytes)
			throw std::runtime_error("invalid checkpoint source layout or byte count");
		bytes = 2*integer_bytes+static_cast<std::uint64_t>(rows)*scalar_bytes;
		offset = 2*integer_bytes+static_cast<std::uint64_t>(begin)*scalar_bytes;
		if (bytes > static_cast<std::uint64_t>(std::numeric_limits<off_t>::max()))
			throw std::runtime_error("checkpoint byte count exceeds off_t");
		values.resize(static_cast<std::size_t>(end-begin));
		PetscReadArray view; view.Acquire(source);
		for (PetscInt i = 0; i < end-begin; ++i) {
			values[static_cast<std::size_t>(i)] = view.Data()[i];
			if (PetscIsInfOrNanScalar(values[static_cast<std::size_t>(i)]))
				throw std::runtime_error("checkpoint source contains a nonfinite value");
		}
		view.Restore();
		description = filename+"\n"+std::to_string(rows)+" "+std::to_string(integer_bytes)+" "+std::to_string(scalar_bytes);
	});
	RequireCollectiveSameText(communicator, "checkpoint write agreement", description);
	int rank = 0;
	MPI_Comm_rank(communicator, &rank);
	petsc_checkpoint_detail::TemporaryCheckpoint temporary;
	std::string temporary_path;
	CollectiveLocalStage(communicator, "checkpoint write file preparation", [&] {
		if (rank == 0) {
			petsc_checkpoint_detail::ReturnErrors errors;
			temporary.Create(filename, bytes, rows);
			temporary_path = temporary.Path();
			if (temporary_path.size() > INT_MAX) throw std::runtime_error("checkpoint temporary path too long");
		}
	});
	int length = static_cast<int>(temporary_path.size());
	MPI_Bcast(&length, 1, MPI_INT, 0, communicator);
	CollectiveLocalStage(communicator, "checkpoint temporary path preparation", [&] {
		temporary_path.resize(static_cast<std::size_t>(length));
	});
	MPI_Bcast(temporary_path.data(), length, MPI_CHAR, 0, communicator);
	CollectiveLocalStage(communicator, "checkpoint local write", [&] {
		petsc_checkpoint_detail::ReturnErrors errors;
		const auto fd = ::open(temporary_path.c_str(), O_WRONLY | O_NONBLOCK);
		if (fd < 0) throw std::runtime_error("cannot open checkpoint temporary file on this rank");
		petsc_checkpoint_detail::LocalBinaryOutput output(fd);
		output.Seek(offset);
		output.Write(values.data(), end-begin, PETSC_SCALAR);
		output.Close();
	});
	CollectiveLocalStage(communicator, "checkpoint state publication", [&] {
		if (rank == 0) temporary.Publish(filename);
	});
}

} // namespace iga

#endif
