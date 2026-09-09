#ifndef IGA_CUDA_EXECUTION_HPP
#define IGA_CUDA_EXECUTION_HPP

#include "ExecutionEnvironment.hpp"
#include "CheckedText.hpp"
#include <algorithm>
#include <cstdint>
#include <limits>
#include <ostream>

namespace iga::cuda {

inline void RequireExecutionEnvironment(std::ostream& output)
{
	const char* launcher = RequireSingleProcessEnvironment();
	const int omp = ResourceThreadSetting("OMP_NUM_THREADS", true);
	const int limit = ResourceThreadSetting("OMP_THREAD_LIMIT");
	const int openblas = ResourceThreadSetting("OPENBLAS_NUM_THREADS", false, true);
	const int mkl = ResourceThreadSetting("MKL_NUM_THREADS", false, true);
	const int blis = ResourceThreadSetting("BLIS_NUM_THREADS", false, true);
#ifdef _OPENMP
	constexpr int openmp = 1;
#else
	constexpr int openmp = 0;
#endif
	output << "execution_resources backend=cuda processes=1 gpu_limit=1 launcher=" << launcher
		<< " scalar=real64 index_bits=" << sizeof(int)*CHAR_BIT
		<< " host_openmp_compiled=" << openmp << " omp_requested=" << omp << " omp_thread_limit=" << limit
		<< " openblas_requested=" << openblas << " mkl_requested=" << mkl
		<< " blis_requested=" << blis << '\n';
	FlushCheckedText(output);
}

inline void RequireMeshIndexCapacity(std::uint64_t nodes, std::uint64_t elements)
{
	if (!nodes || nodes > static_cast<std::uint64_t>(INT_MAX))
		throw std::runtime_error("CUDA node count exceeds int32 capacity");
	if (!elements || elements > static_cast<std::uint64_t>(INT_MAX))
		throw std::runtime_error("CUDA element count exceeds int32 capacity");
}

inline void RequireFieldIndexCapacity(std::uint64_t nodes, std::size_t fields)
{
	// Both nodal fields and the three-component prescribed velocity use int
	// subscripts; GMRES/cuBLAS also receive the field vector length as int.
	if (!fields || fields > 8 || !nodes
		|| nodes > static_cast<std::uint64_t>(INT_MAX)/std::max<std::size_t>(fields, 3))
		throw std::runtime_error("CUDA node/field rows exceed int32 capacity or supported field count");
}

inline std::size_t CheckedBlockValueCount(std::size_t blocks, unsigned fields)
{
	// Kernels already widen block-value offsets to size_t. Preserve that range;
	// only the block IDs themselves are int32. Check bytes before multiplying.
	if (!fields || fields > 8 || blocks > static_cast<std::size_t>(INT_MAX)
		|| blocks > std::numeric_limits<std::size_t>::max()/sizeof(double)/(fields*fields))
		throw std::runtime_error("CUDA block matrix values exceed allocation/index capacity");
	return blocks*fields*fields;
}

} // namespace iga::cuda
#endif
