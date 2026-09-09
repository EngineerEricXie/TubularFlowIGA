#ifndef IGA_EXECUTION_ENVIRONMENT_HPP
#define IGA_EXECUTION_ENVIRONMENT_HPP

#include <cerrno>
#include <climits>
#include <cstdlib>
#include <cctype>
#include <cstring>
#include <stdexcept>
#include <string>

namespace iga {

// Parse documented numeric resource settings without changing environment or
// runtime policy. -1 means unset; a BLAS zero is retained as a library default.
inline int ResourceThreadSetting(const char* name, bool list = false, bool allow_zero = false)
{
	const auto* value = std::getenv(name);
	if (!value) return -1;
	int first = -1;
	for (;;) {
		while (std::isspace(static_cast<unsigned char>(*value))) ++value;
		char* end = nullptr;
		errno = 0;
		const auto number = std::strtol(value, &end, 10);
		if (end == value || errno == ERANGE || number < (allow_zero ? 0 : 1) || number > INT_MAX)
			throw std::runtime_error(std::string(name)+" requires "+(allow_zero ? "nonnegative" : "positive")+" integer thread counts");
		if (first < 0) first = static_cast<int>(number);
		while (std::isspace(static_cast<unsigned char>(*end))) ++end;
		if (!*end) return first;
		if (!list || *end != ',') throw std::runtime_error(std::string(name)+" has an invalid thread count list");
		value = end+1;
	}
}

inline int ExecutionInteger(const char* name, bool allow_zero = false)
{
	const auto* value = std::getenv(name);
	if (!value || !*value) throw std::runtime_error(std::string(name)+" requires an integer");
	int result = 0;
	for (const auto* digit = value; *digit; ++digit) {
		if (*digit < '0' || *digit > '9' || result > (INT_MAX-(*digit-'0'))/10)
			throw std::runtime_error(std::string(name)+" requires an integer within INT_MAX");
		result = result*10+(*digit-'0');
	}
	if (!allow_zero && !result) throw std::runtime_error(std::string(name)+" must be positive");
	return result;
}

inline void RequireSingleProcessPair(const char* rank_name, const char* size_name)
{
	const auto rank = ExecutionInteger(rank_name, true);
	const auto size = ExecutionInteger(size_name);
	if (rank >= size) throw std::runtime_error(std::string(rank_name)+" is outside the launcher size");
	if (size != 1) throw std::runtime_error(std::string("single-process executable requires one task; ")
		+size_name+"="+std::to_string(size));
}

// MPI launcher metadata takes precedence over the enclosing Slurm allocation.
// A batch shell is one process even when its allocation requests many tasks.
// This detects the supported launchers; it is not an output-directory lock.
inline const char* RequireSingleProcessEnvironment()
{
	if (std::getenv("OMPI_COMM_WORLD_RANK") || std::getenv("OMPI_COMM_WORLD_SIZE")) {
		RequireSingleProcessPair("OMPI_COMM_WORLD_RANK", "OMPI_COMM_WORLD_SIZE");
		return "openmpi";
	}
	if (std::getenv("PMI_RANK") || std::getenv("PMI_SIZE")) {
		RequireSingleProcessPair("PMI_RANK", "PMI_SIZE");
		return "pmi";
	}
	const char* step = std::getenv("SLURM_STEP_ID");
	if (!step) step = std::getenv("SLURM_STEPID");
	if (step && (!std::strcmp(step, "batch") || !std::strcmp(step, "extern"))) return "slurm_allocation";
	if (std::getenv("SLURM_PROCID")) {
		const char* size = std::getenv("SLURM_STEP_NUM_TASKS") ? "SLURM_STEP_NUM_TASKS" : "SLURM_NTASKS";
		RequireSingleProcessPair("SLURM_PROCID", size);
		return "slurm_step";
	}
	return "direct";
}

} // namespace iga
#endif
