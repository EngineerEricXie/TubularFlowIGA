#ifndef IGA_OWNED_SCALAR_CONTRIBUTIONS_HPP
#define IGA_OWNED_SCALAR_CONTRIBUTIONS_HPP

#include "DistributedPointIdentity.hpp"
#include <cmath>
#include <cstring>

namespace iga {

// Sparse scalar reduction to publication owners. IDs use the full UInt64
// range. A routing shard collects each ID's owner and contributions, then
// returns one sum per owned ID. No full global field/catalog is gathered.
// Duplicate owners and contributions to an unowned ID are collective errors;
// repeated contributions are intentional. The caller certifies their physical
// integration authority. Wire limits are not RSS limits; skew may concentrate
// a shard. Live-rank errors are coordinated; process loss is not supported.
inline std::vector<double> SumOwnedScalarContributions(MPI_Comm comm,
	const std::vector<std::uint64_t>& owned,
	const std::vector<std::pair<std::uint64_t,double>>& contributions,
	PointIdentityLimits limits={})
{
	using namespace point_identity_detail;
	static_assert(sizeof(double)==sizeof(std::uint64_t),"scalar reduction requires 64-bit double");
	int ranks=0;MPI_Comm_size(comm,&ranks);
	std::map<std::uint64_t,std::size_t> rows;
	std::vector<std::string> buckets;
	CollectiveLocalStage(comm,"owned scalar input",[&] {
		if(!limits.max_wire_bytes||limits.max_wire_bytes>static_cast<std::size_t>(std::numeric_limits<int>::max())
			||!limits.max_local_occurrences||owned.size()>limits.max_local_occurrences
			||contributions.size()>limits.max_local_occurrences-owned.size())
			throw std::invalid_argument("invalid owned scalar limits");
		buckets.resize(ranks);std::size_t used=0;
		for(std::size_t i=0;i<owned.size();++i) {
			if(!rows.emplace(owned[i],i).second)throw std::invalid_argument("duplicate local scalar owner");
			Charge(used,24,limits.max_wire_bytes);auto& target=buckets[owned[i]%ranks];
			Put(target,owned[i]);Put(target,0);Put(target,0);
		}
		for(const auto& item:contributions) {
			if(!std::isfinite(item.second))throw std::invalid_argument("nonfinite scalar contribution");
			Charge(used,24,limits.max_wire_bytes);auto& target=buckets[item.first%ranks];
			std::uint64_t bits=0;std::memcpy(&bits,&item.second,sizeof(bits));
			Put(target,item.first);Put(target,1);Put(target,bits);
		}
	});
	{
		const auto incoming=Exchange(comm,buckets,limits.max_wire_bytes);
		CollectiveLocalStage(comm,"owned scalar shard reduction",[&] {
			struct Entry { int owner=-1;long double sum=0.; };
			std::map<std::uint64_t,Entry> entries;
			for(int peer=0;peer<ranks;++peer) {
				auto packet=incoming.Peer(peer);
				while(!packet.empty()) {
					const auto id=Get(packet),kind=Get(packet),bits=Get(packet);auto& entry=entries[id];
					if(kind==0) {
						if(entry.owner!=-1)throw std::runtime_error("duplicate global scalar owner");
						entry.owner=peer;
					} else if(kind==1) {
						double value=0.;std::memcpy(&value,&bits,sizeof(bits));entry.sum+=static_cast<long double>(value);
						if(!std::isfinite(entry.sum))throw std::overflow_error("scalar sum overflow");
					} else throw std::runtime_error("invalid scalar packet kind");
				}
			}
			buckets.assign(ranks,{});std::size_t used=0;
			for(const auto& item:entries) {
				if(item.second.owner<0)throw std::runtime_error("scalar contribution has no owner");
				if(std::abs(item.second.sum)>static_cast<long double>(std::numeric_limits<double>::max()))
					throw std::overflow_error("scalar sum exceeds double range");
				const double value=static_cast<double>(item.second.sum);
				if(!std::isfinite(value))throw std::overflow_error("scalar sum exceeds double range");
				Charge(used,16,limits.max_wire_bytes);auto& target=buckets[item.second.owner];
				std::uint64_t bits=0;std::memcpy(&bits,&value,sizeof(bits));Put(target,item.first);Put(target,bits);
			}
		});
	}
	const auto incoming=Exchange(comm,buckets,limits.max_wire_bytes);
	std::vector<double> result;
	CollectiveLocalStage(comm,"owned scalar result",[&] {
		result.resize(owned.size());std::vector<bool> seen(owned.size(),false);
		for(int peer=0;peer<ranks;++peer) {
			auto packet=incoming.Peer(peer);
			while(!packet.empty()) {
				const auto id=Get(packet),bits=Get(packet);const auto found=rows.find(id);
				if(found==rows.end()||id%ranks!=static_cast<std::uint64_t>(peer)||seen[found->second])
					throw std::runtime_error("invalid owned scalar response");
				seen[found->second]=true;std::memcpy(&result[found->second],&bits,sizeof(bits));
				if(!std::isfinite(result[found->second]))throw std::runtime_error("nonfinite owned scalar response");
			}
		}
		for(bool value:seen)if(!value)throw std::runtime_error("missing owned scalar response");
	});
	return result;
}
} // namespace iga
#endif
