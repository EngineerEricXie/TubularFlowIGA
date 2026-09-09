#include "WeightedWorkPartition.hpp"
#include <functional>
#include <iostream>

namespace {
std::uint64_t Brute(const std::vector<std::uint64_t>& weights,std::size_t begin,int ranks)
{
	if (begin == weights.size()) return 0;
	std::uint64_t result = std::numeric_limits<std::uint64_t>::max(),sum = 0;
	for (std::size_t end = begin+1; end <= weights.size(); ++end) {
		sum += weights[end-1];
		if (ranks > 1 || end == weights.size())
			result = std::min(result,std::max(sum,end == weights.size() ? 0 : Brute(weights,end,ranks-1)));
	}
	return result;
}
void Check(const std::vector<std::uint64_t>& weights,int ranks)
{
	const auto owners = iga::PartitionContiguousWork(weights,ranks);
	std::vector<std::uint64_t> totals(ranks,0);
	for (std::size_t i = 0; i < weights.size(); ++i) {
		if (owners[i] < 0 || owners[i] >= ranks || (i && owners[i] < owners[i-1])) throw std::runtime_error("invalid work ownership");
		totals[owners[i]] += weights[i];
	}
	if (*std::max_element(totals.begin(),totals.end()) != Brute(weights,0,ranks)) throw std::runtime_error("partition does not attain the optimal contiguous bound");
}
}
int main()
{
	try {
		for (std::size_t n = 0; n <= 7; ++n) {
			std::vector<std::uint64_t> weights(n,1);
			std::function<void(std::size_t)> enumerate = [&](std::size_t index) {
				if (index == n) { for (int ranks = 1; ranks <= 5; ++ranks) Check(weights,ranks); return; }
				for (std::uint64_t weight : {1u,3u,11u}) { weights[index] = weight; enumerate(index+1); }
			};
			enumerate(0);
		}
		const auto maximum = std::numeric_limits<std::uint64_t>::max();
		Check({maximum},4);
		for (int scenario = 0; scenario < 4; ++scenario) {
			bool rejected = false;
			try {
				if (scenario == 0) iga::PartitionContiguousWork({1},0);
				if (scenario == 1) iga::PartitionContiguousWork({0},2);
				if (scenario == 2) iga::PartitionContiguousWork({maximum,1},2);
				if (scenario == 3) iga::CheckedWorkProduct(maximum,2);
			} catch (const std::exception&) { rejected = true; }
			if (!rejected) throw std::runtime_error("invalid work was accepted");
		}
		std::cout << "weighted contiguous work partition passed\n";
	} catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
