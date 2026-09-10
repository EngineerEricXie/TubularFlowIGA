#ifndef IGA_DISTRIBUTED_POINT_VALUES_HPP
#define IGA_DISTRIBUTED_POINT_VALUES_HPP

#include "DistributedPointIdentity.hpp"
#include <cmath>
#include <cstring>

namespace iga {

// Values are interleaved by local occurrence, with a common component count.
// Every result row is copied from its elected representative, including local
// representatives. Coordinates and fields can be exchanged in the same tuple.
// Input arrays remain unchanged on every caught collective failure.
inline std::vector<double> ExchangePointRepresentativeValues(MPI_Comm comm,
	const std::vector<std::int64_t>& local_ids,const std::vector<PointIdentityRepresentative>& representatives,
	const std::vector<double>& values,int components,PointIdentityLimits limits={})
{
	using namespace point_identity_detail;
	static_assert(sizeof(double)==sizeof(std::uint64_t),"point value exchange requires 64-bit double");
	int ranks=0;MPI_Comm_size(comm,&ranks);
	std::map<std::int64_t,std::size_t> rows;
	std::vector<std::string> buckets;
	CollectiveLocalStage(comm,"point values input validation",[&] {
		if(components<1||!limits.max_wire_bytes||limits.max_wire_bytes>static_cast<std::size_t>(std::numeric_limits<int>::max())
			||!limits.max_local_occurrences||local_ids.size()>limits.max_local_occurrences
			||representatives.size()!=local_ids.size()
			||local_ids.size()>values.max_size()/static_cast<std::size_t>(components)
			||values.size()!=local_ids.size()*static_cast<std::size_t>(components))
			throw std::invalid_argument("invalid point value dimensions or limits");
		for(double value:values)if(!std::isfinite(value))throw std::invalid_argument("nonfinite point representative value");
		buckets.resize(ranks);std::size_t used=0;
		for(std::size_t i=0;i<local_ids.size();++i) {
			const auto& representative=representatives[i];
			if(local_ids[i]<0||!rows.emplace(local_ids[i],i).second||representative.occurrence<0
				||representative.rank<0||representative.rank>=ranks)
				throw std::invalid_argument("invalid point value occurrence or owner");
			Charge(used,16,limits.max_wire_bytes);
			Put(buckets[representative.rank],i);Put(buckets[representative.rank],static_cast<std::uint64_t>(representative.occurrence));
		}
	});
	RequireCollectiveSameInt(comm,"point values component agreement",components);
	{
		const auto requests=Exchange(comm,buckets,limits.max_wire_bytes);
		CollectiveLocalStage(comm,"point values response construction",[&] {
			buckets.assign(ranks,{});std::size_t used=0;
			for(int peer=0;peer<ranks;++peer) {
				auto packet=requests.Peer(peer);
				while(!packet.empty()) {
					const auto query=Get(packet),id=Get(packet);
					if(id>static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()))throw std::runtime_error("invalid representative occurrence id");
					const auto found=rows.find(static_cast<std::int64_t>(id));
					if(found==rows.end())throw std::runtime_error("point representative is absent from its owner");
					Charge(used,8,limits.max_wire_bytes);Put(buckets[peer],query);
					const auto offset=found->second*static_cast<std::size_t>(components);
					for(int component=0;component<components;++component) {
						Charge(used,8,limits.max_wire_bytes);
						std::uint64_t bits=0;std::memcpy(&bits,&values[offset+component],sizeof(bits));Put(buckets[peer],bits);
					}
				}
			}
		});
	}
	const auto replies=Exchange(comm,buckets,limits.max_wire_bytes);
	std::vector<double> result;
	CollectiveLocalStage(comm,"point values response validation",[&] {
		result.resize(values.size());std::vector<bool> seen(local_ids.size(),false);
		for(int peer=0;peer<ranks;++peer) {
			auto packet=replies.Peer(peer);
			while(!packet.empty()) {
				const auto query=Get(packet);
				if(query>=local_ids.size()||seen[query]||representatives[query].rank!=peer)
					throw std::runtime_error("invalid point value response source or row");
				seen[query]=true;const auto offset=static_cast<std::size_t>(query)*static_cast<std::size_t>(components);
				for(int component=0;component<components;++component) {
					const auto bits=Get(packet);std::memcpy(&result[offset+component],&bits,sizeof(bits));
					if(!std::isfinite(result[offset+component]))throw std::runtime_error("nonfinite received point value");
				}
			}
		}
		for(bool value:seen)if(!value)throw std::runtime_error("missing point value response");
	});
	return result;
}
} // namespace iga
#endif
