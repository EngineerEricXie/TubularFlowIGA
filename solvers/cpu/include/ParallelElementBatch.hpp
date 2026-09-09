#ifndef IGA_PARALLEL_ELEMENT_BATCH_HPP
#define IGA_PARALLEL_ELEMENT_BATCH_HPP

#include <algorithm>
#include <cstddef>
#include <exception>
#include <limits>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace iga {

struct ElementBatchOptions {
	int threads = 1;
	std::size_t capacity = 8;
};

struct ElementBatchStatistics {
	std::size_t items = 0, batches = 0, maximum_resident_items = 0;
	int maximum_team_size = 1;
};

// Prepare and consume execute on the calling thread. Compute may run in
// OpenMP workers and must use immutable inputs and private scratch, with no
// PETSc/MPI calls. Every batch joins and checks all errors before consumption;
// insertion order is independent of worker scheduling. The caller remains
// responsible for rollback of any previously consumed batches after failure.
template <class Prepare, class Compute, class Consume>
ElementBatchStatistics ForEachElementBatch(std::size_t count, ElementBatchOptions options,
	Prepare prepare, Compute compute, Consume consume)
{
	if (options.threads < 1 || options.capacity == 0
		|| options.capacity > static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()))
		throw std::invalid_argument("invalid element batch threads or capacity");
#ifdef _OPENMP
	if (omp_in_parallel()) throw std::invalid_argument("element batches cannot start inside an OpenMP region");
#else
	if (options.threads != 1) throw std::invalid_argument("parallel element assembly requires an OpenMP build");
#endif
	using Prepared = std::decay_t<decltype(prepare(std::size_t{}))>;
	using Result = std::decay_t<decltype(compute(std::declval<const Prepared&>(), std::size_t{}))>;
	const auto capacity = std::min(count, options.capacity);
	std::vector<std::optional<Prepared>> prepared(capacity);
	std::vector<std::optional<Result>> results(capacity);
	std::vector<std::exception_ptr> errors(capacity);
	ElementBatchStatistics statistics;
	for (std::size_t begin = 0; begin < count;) {
		const auto length = std::min(capacity, count-begin);
		for (std::size_t i = 0; i < length; ++i) prepared[i].emplace(prepare(begin+i));
		int team_size = 1;
#ifdef _OPENMP
#pragma omp parallel num_threads(options.threads) if(options.threads > 1)
		{
#pragma omp single
			team_size = omp_get_num_threads();
#pragma omp for schedule(dynamic, 1)
#endif
			for (std::ptrdiff_t index = 0; index < static_cast<std::ptrdiff_t>(length); ++index) {
				const auto i = static_cast<std::size_t>(index);
				try { results[i].emplace(compute(std::as_const(*prepared[i]), begin+i)); }
				catch (...) { errors[i] = std::current_exception(); }
			}
#ifdef _OPENMP
		}
#endif
		for (std::size_t i = 0; i < length; ++i)
			if (errors[i]) std::rethrow_exception(errors[i]);
		for (std::size_t i = 0; i < length; ++i) {
			consume(*prepared[i], *results[i], begin+i);
			results[i].reset(); prepared[i].reset();
		}
		statistics.items += length;
		++statistics.batches;
		statistics.maximum_resident_items = std::max(statistics.maximum_resident_items, length);
		statistics.maximum_team_size = std::max(statistics.maximum_team_size, team_size);
		begin += length;
	}
	return statistics;
}

} // namespace iga

#endif
