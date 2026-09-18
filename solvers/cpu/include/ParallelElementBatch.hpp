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
#include "AssemblyDetailProfile.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace iga {

struct ElementBatchOptions {
	int threads = 1;
	std::size_t capacity = 8;
	bool profile_detail = false;
};

struct ElementBatchStatistics {
	std::size_t items = 0, batches = 0, maximum_resident_items = 0;
	std::size_t maximum_resident_result_payload_bytes = 0;
	int maximum_team_size = 1;
	double prepare_wall_seconds=0.,compute_wall_seconds=0.,consume_wall_seconds=0.;
	double worker_work_seconds=0.,maximum_worker_seconds=0.,maximum_worker_completion_spread_seconds=0.;
};

namespace parallel_element_batch_detail {
template<class T> auto ResultPayloadBytes(const T& value,int) noexcept -> decltype(value.ResultPayloadBytes()) { return value.ResultPayloadBytes(); }
template<class T> std::size_t ResultPayloadBytes(const T&,long) noexcept { return sizeof(T); }
}

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
	std::vector<double> item_seconds(options.profile_detail?capacity:0);
	std::vector<int> item_worker(options.profile_detail?capacity:0);
	std::vector<std::size_t> item_result_bytes(options.profile_detail?capacity:0);
	ElementBatchStatistics statistics;
	for (std::size_t begin = 0; begin < count;) {
		const auto length = std::min(capacity, count-begin);
		{
			DetailTimer timer(statistics.prepare_wall_seconds,options.profile_detail);
			for (std::size_t i = 0; i < length; ++i) prepared[i].emplace(prepare(begin+i));
		}
		int team_size = 1;
		const double compute_start=options.profile_detail?DetailNow():0.;
#ifdef _OPENMP
#pragma omp parallel num_threads(options.threads) if(options.threads > 1)
		{
#pragma omp single
			team_size = omp_get_num_threads();
#pragma omp for schedule(dynamic, 1)
#endif
			for (std::ptrdiff_t index = 0; index < static_cast<std::ptrdiff_t>(length); ++index) {
				const auto i = static_cast<std::size_t>(index);
				const double item_start=options.profile_detail?DetailNow():0.;
				try { results[i].emplace(compute(std::as_const(*prepared[i]), begin+i)); }
				catch (...) { errors[i] = std::current_exception(); }
				if(options.profile_detail) {
					if(results[i]) item_result_bytes[i]=parallel_element_batch_detail::ResultPayloadBytes(*results[i],0);
					item_seconds[i]=DetailNow()-item_start;
#ifdef _OPENMP
					item_worker[i]=omp_get_thread_num();
#else
					item_worker[i]=0;
#endif
				}
			}
#ifdef _OPENMP
		}
#endif
		for (std::size_t i = 0; i < length; ++i)
			if (errors[i]) std::rethrow_exception(errors[i]);
		if(options.profile_detail) {
			statistics.compute_wall_seconds+=DetailNow()-compute_start;
			std::vector<double> work(static_cast<std::size_t>(team_size),0.);
			for(std::size_t i=0;i<length;++i) { work[item_worker[i]]+=item_seconds[i]; statistics.worker_work_seconds+=item_seconds[i]; }
			const auto limits=std::minmax_element(work.begin(),work.end());
			statistics.maximum_worker_seconds+=*limits.second;
			statistics.maximum_worker_completion_spread_seconds=std::max(statistics.maximum_worker_completion_spread_seconds,*limits.second-*limits.first);
			std::size_t resident=0;
			for(std::size_t i=0;i<length;++i) {
				if(item_result_bytes[i]>std::numeric_limits<std::size_t>::max()-resident) throw std::overflow_error("element batch result payload bytes overflow");
				resident+=item_result_bytes[i];
			}
			statistics.maximum_resident_result_payload_bytes=std::max(statistics.maximum_resident_result_payload_bytes,resident);
		}
		{
		DetailTimer timer(statistics.consume_wall_seconds,options.profile_detail);
		for (std::size_t i = 0; i < length; ++i) {
			consume(*prepared[i], *results[i], begin+i);
			results[i].reset(); prepared[i].reset();
		}
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
