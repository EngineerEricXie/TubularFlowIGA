#ifndef IGA_OWNED_POINT_VALUES_HPP
#define IGA_OWNED_POINT_VALUES_HPP

#include "DistributedPointIdentity.hpp"
#include <cmath>
#include <cstring>

namespace iga {

// Fetch interleaved tuples by full UInt64 node ID without a global owner map.
// Owners publish to ID shards; queries use the same shards. Duplicate queries
// are allowed and retain their original order. Publication authority and field
// epoch binding belong to the caller; duplicate owners are rejected here.
// All live ranks enter, including empty owners/requesters. Wire caps do not
// bound allocator RSS and skew can concentrate a shard. Inputs are immutable.
inline std::vector<double> FetchOwnedPointValues(MPI_Comm comm,
	const std::vector<std::uint64_t>& owned,const std::vector<double>& values,
	const std::vector<std::uint64_t>& queries,int components,PointIdentityLimits limits={})
{
	using namespace point_identity_detail;
	static_assert(sizeof(double)==sizeof(std::uint64_t),"point exchange requires 64-bit double");
	int ranks=0;MPI_Comm_size(comm,&ranks);
	std::vector<std::string> buckets;
	CollectiveLocalStage(comm,"owned point input",[&] {
		if(components<1||!limits.max_wire_bytes||limits.max_wire_bytes>static_cast<std::size_t>(std::numeric_limits<int>::max())
			||!limits.max_local_occurrences||owned.size()>limits.max_local_occurrences
			||queries.size()>limits.max_local_occurrences-owned.size()
			||owned.size()>values.max_size()/static_cast<std::size_t>(components)
			||queries.size()>values.max_size()/static_cast<std::size_t>(components)
			||values.size()!=owned.size()*static_cast<std::size_t>(components))
			throw std::invalid_argument("invalid owned point dimensions or limits");
		std::set<std::uint64_t> unique;
		buckets.resize(ranks);std::size_t used=0;
		for(std::size_t row=0;row<owned.size();++row) {
			if(!unique.insert(owned[row]).second)throw std::invalid_argument("duplicate local point owner");
			Charge(used,16,limits.max_wire_bytes);auto& target=buckets[owned[row]%ranks];
			Put(target,0);Put(target,owned[row]);
			for(int c=0;c<components;++c) {
				const double value=values[row*static_cast<std::size_t>(components)+c];
				if(!std::isfinite(value))throw std::invalid_argument("nonfinite owned point value");
				Charge(used,8,limits.max_wire_bytes);std::uint64_t bits=0;
				std::memcpy(&bits,&value,sizeof(bits));Put(target,bits);
			}
		}
		for(std::size_t row=0;row<queries.size();++row) {
			Charge(used,24,limits.max_wire_bytes);auto& target=buckets[queries[row]%ranks];
			Put(target,1);Put(target,queries[row]);Put(target,row);
		}
	});
	RequireCollectiveSameInt(comm,"owned point components",components);
	{
		const auto incoming=Exchange(comm,buckets,limits.max_wire_bytes);
		CollectiveLocalStage(comm,"owned point shard replies",[&] {
			struct Query { std::uint64_t id,row;int peer; };
			std::map<std::uint64_t,std::vector<std::uint64_t>> tuples;
			std::vector<Query> requests;
			for(int peer=0;peer<ranks;++peer) {
				auto packet=incoming.Peer(peer);
				while(!packet.empty()) {
					const auto kind=Get(packet),id=Get(packet);
					if(kind==0) {
						auto added=tuples.emplace(id,std::vector<std::uint64_t>{});
						if(!added.second)throw std::runtime_error("duplicate global point owner");
						auto& tuple=added.first->second;tuple.reserve(components);
						for(int c=0;c<components;++c)tuple.push_back(Get(packet));
					} else if(kind==1)requests.push_back({id,Get(packet),peer});
					else throw std::runtime_error("invalid owned point packet kind");
				}
			}
			buckets.assign(ranks,{});std::size_t used=0;
			for(const auto& query:requests) {
				const auto found=tuples.find(query.id);
				if(found==tuples.end())throw std::runtime_error("queried point has no owner");
				Charge(used,8,limits.max_wire_bytes);auto& target=buckets[query.peer];Put(target,query.row);
				for(auto bits:found->second){Charge(used,8,limits.max_wire_bytes);Put(target,bits);}
			}
		});
	}
	const auto incoming=Exchange(comm,buckets,limits.max_wire_bytes);
	std::vector<double> result;
	CollectiveLocalStage(comm,"owned point response",[&] {
		result.resize(queries.size()*static_cast<std::size_t>(components));std::vector<bool> seen(queries.size(),false);
		for(int peer=0;peer<ranks;++peer) {
			auto packet=incoming.Peer(peer);
			while(!packet.empty()) {
				const auto row=Get(packet);
				if(row>=queries.size()||seen[row]||queries[row]%ranks!=static_cast<std::uint64_t>(peer))
					throw std::runtime_error("invalid owned point response row or source");
				seen[row]=true;
				for(int c=0;c<components;++c) {
					const auto bits=Get(packet);auto& value=result[row*static_cast<std::size_t>(components)+c];
					std::memcpy(&value,&bits,sizeof(bits));
					if(!std::isfinite(value))throw std::runtime_error("nonfinite fetched point value");
				}
			}
		}
		for(bool value:seen)if(!value)throw std::runtime_error("missing owned point response");
	});
	return result;
}
} // namespace iga
#endif
