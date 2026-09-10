#ifndef IGA_OWNED_VECTOR_CONTRIBUTIONS_HPP
#define IGA_OWNED_VECTOR_CONTRIBUTIONS_HPP

#include "OwnedScalarContributions.hpp"
#include <array>

namespace iga {

// Vector force contributions retain material IDs without encoding a component
// into the ID (which would overflow the UInt64 catalog). The three components
// use consecutive scalar exchanges. This is a correctness bridge, with three
// routing passes; limits apply to each pass, not aggregate traffic or RSS.
inline std::vector<std::array<double,3>> SumOwnedVectorContributions(MPI_Comm comm,
	const std::vector<std::uint64_t>& owned,
	const std::vector<std::pair<std::uint64_t,std::array<double,3>>>& contributions,
	PointIdentityLimits limits={})
{
	std::vector<std::array<double,3>> result;
	std::vector<std::pair<std::uint64_t,double>> scalar;
	CollectiveLocalStage(comm,"owned vector input",[&] {
		if(!limits.max_local_occurrences||owned.size()>limits.max_local_occurrences
			||contributions.size()>limits.max_local_occurrences-owned.size())
			throw std::invalid_argument("owned vector contributions exceed record limit");
		for(const auto& contribution:contributions)for(double value:contribution.second)
			if(!std::isfinite(value))throw std::invalid_argument("nonfinite vector contribution");
		result.resize(owned.size());scalar.resize(contributions.size());
	});
	for(int component=0;component<3;++component) {
		for(std::size_t row=0;row<contributions.size();++row)
			scalar[row]={contributions[row].first,contributions[row].second[component]};
		const auto reduced=SumOwnedScalarContributions(comm,owned,scalar,limits);
		for(std::size_t row=0;row<owned.size();++row)result[row][component]=reduced[row];
	}
	return result;
}
} // namespace iga
#endif
