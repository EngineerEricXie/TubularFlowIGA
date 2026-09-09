#ifndef IGA_ELEMENT_ASSEMBLY_EXECUTION_HPP
#define IGA_ELEMENT_ASSEMBLY_EXECUTION_HPP

#include "ExecutionResources.hpp"
#include "ParallelElementBatch.hpp"
#include <thread>

namespace iga {

// Execution settings are frozen at construction and are not physical inputs.
// This object must be constructed and used by the MPI initialization thread.
class ElementAssemblyExecution {
public:
	ElementAssemblyExecution()
		: caller_(std::this_thread::get_id())
	{
#ifdef _OPENMP
		if (omp_in_parallel()) throw std::invalid_argument("element assembly cannot start inside an OpenMP region");
		options_.threads = omp_get_max_threads();
#endif
		const int requested = ResourceThreadSetting("IGA_ASSEMBLY_THREADS");
		if (requested != -1) options_.threads = requested;
#ifndef _OPENMP
		if (options_.threads != 1) throw std::invalid_argument("parallel element assembly requires an OpenMP build");
#endif
		options_.capacity = static_cast<std::size_t>(std::min(options_.threads, 8));
		const int capacity = ResourceThreadSetting("IGA_ASSEMBLY_BATCH_SIZE");
		if (capacity != -1) options_.capacity = static_cast<std::size_t>(capacity);
		int main_thread = 0, provided = MPI_THREAD_SINGLE;
		if (MPI_Is_thread_main(&main_thread) != MPI_SUCCESS || !main_thread)
			throw std::invalid_argument("element assembly requires the MPI initialization thread");
		if (MPI_Query_thread(&provided) != MPI_SUCCESS)
			throw std::runtime_error("cannot query element assembly MPI thread support");
		if (options_.threads > 1 && provided < MPI_THREAD_FUNNELED)
			throw std::invalid_argument("parallel element assembly requires MPI_THREAD_FUNNELED or stronger");
	}

	void RequireCaller() const
	{
		// Reject a foreign caller before any PETSc or MPI operation.
		if (std::this_thread::get_id() != caller_)
			throw std::logic_error("element assembly called from a different thread");
#ifdef _OPENMP
		if (omp_in_parallel()) throw std::logic_error("element assembly called inside an OpenMP region");
#endif
	}

	ElementBatchOptions Options() const noexcept { return options_; }

private:
	std::thread::id caller_;
	ElementBatchOptions options_;
};

} // namespace iga

#endif
