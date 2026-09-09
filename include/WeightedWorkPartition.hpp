#ifndef IGA_WEIGHTED_WORK_PARTITION_HPP
#define IGA_WEIGHTED_WORK_PARTITION_HPP

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace iga {

inline std::uint64_t CheckedWorkSum(std::uint64_t a,std::uint64_t b)
{
	if (b > std::numeric_limits<std::uint64_t>::max()-a) throw std::overflow_error("partition work sum overflows");
	return a+b;
}
inline std::uint64_t CheckedWorkProduct(std::uint64_t a,std::uint64_t b)
{
	if (a && b > std::numeric_limits<std::uint64_t>::max()/a) throw std::overflow_error("partition work product overflows");
	return a*b;
}

// Preserve the input spatial order. Binary search the smallest maximum work
// that permits at most ranks contiguous groups; greedy placement then realizes
// that bound. Empty ranks are allowed, but every supplied item must cost > 0.
inline std::vector<int> PartitionContiguousWork(const std::vector<std::uint64_t>& weights,int ranks)
{
	if (ranks <= 0) throw std::invalid_argument("partition ranks must be positive");
	std::uint64_t lower = 0,upper = 0;
	for (auto weight : weights) {
		if (!weight) throw std::invalid_argument("partition work must be positive");
		lower = std::max(lower,weight); upper = CheckedWorkSum(upper,weight);
	}
	while (lower < upper) {
		const auto limit = lower+(upper-lower)/2;
		std::size_t groups = 1; std::uint64_t used = 0;
		for (auto weight : weights) {
			if (weight > limit-used) { ++groups; used = 0; }
			used += weight;
		}
		if (groups <= static_cast<std::size_t>(ranks)) upper = limit;
		else lower = limit+1;
	}
	std::vector<int> owners; owners.reserve(weights.size());
	int rank = 0; std::uint64_t used = 0;
	for (auto weight : weights) {
		if (weight > lower-used) { ++rank; used = 0; }
		owners.push_back(rank); used += weight;
	}
	return owners;
}

} // namespace iga
#endif
