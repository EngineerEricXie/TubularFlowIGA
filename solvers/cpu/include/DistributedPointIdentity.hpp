#ifndef IGA_DISTRIBUTED_POINT_IDENTITY_HPP
#define IGA_DISTRIBUTED_POINT_IDENTITY_HPP

#include "CollectiveFailure.hpp"
#include <limits>
#include <map>
#include <set>
#include <vector>

namespace iga {
struct PointIdentityOccurrence {
	std::int64_t occurrence=-1;
	std::string key;
};
struct PointIdentityRepresentative {
	std::int64_t occurrence=-1;
	int rank=-1;
};
struct PointIdentityLimits {
	// Wire payload per rank, in each direction of each exchange. These are
	// not allocator/RSS limits; buckets and flattened buffers may coexist.
	std::size_t max_wire_bytes=64u*1024u*1024u;
	std::size_t max_local_occurrences=1000000;
};
namespace point_identity_detail {
inline std::uint64_t BucketHash(std::string_view key)
{
	std::uint64_t hash=1469598103934665603ULL;
	for(unsigned char byte:key) { hash^=byte;hash*=1099511628211ULL; }
	return hash;
}
inline void Put(std::string& buffer,std::uint64_t value)
{
	for(unsigned q=0;q<8;++q)buffer.push_back(static_cast<char>((value>>(8*q))&255u));
}
inline std::uint64_t Get(std::string_view& buffer)
{
	if(buffer.size()<8)throw std::runtime_error("truncated point identity packet");
	std::uint64_t value=0;
	for(unsigned q=0;q<8;++q)value|=static_cast<std::uint64_t>(static_cast<unsigned char>(buffer[q]))<<(8*q);
	buffer.remove_prefix(8);return value;
}
inline std::string_view Key(std::string_view& buffer)
{
	const auto size=Get(buffer);
	if(!size||size>buffer.size())throw std::runtime_error("invalid point identity key length");
	const auto key=buffer.substr(0,static_cast<std::size_t>(size));buffer.remove_prefix(static_cast<std::size_t>(size));return key;
}
inline void Charge(std::size_t& used,std::size_t count,std::size_t cap)
{
	if(used>cap||count>cap-used)throw std::runtime_error("point identity wire cap reached");
	used+=count;
}
struct Packets {
	std::vector<char> bytes;
	std::vector<int> counts,offsets;
	std::string_view Peer(int rank) const
	{
		return counts[rank]?std::string_view(bytes.data()+offsets[rank],static_cast<std::size_t>(counts[rank])):std::string_view();
	}
};
inline Packets Exchange(MPI_Comm comm,const std::vector<std::string>& buckets,std::size_t cap)
{
	int ranks=0;MPI_Comm_size(comm,&ranks);
	std::vector<int> counts,offsets;std::vector<char> outgoing;Packets received;
	CollectiveLocalStage(comm,"point identity exchange preparation",[&] {
		counts.resize(ranks);offsets.resize(ranks);received.counts.resize(ranks);received.offsets.resize(ranks);
		std::size_t total=0;
		for(int peer=0;peer<ranks;++peer) {
			offsets[peer]=static_cast<int>(total);Charge(total,buckets.at(peer).size(),cap);
			counts[peer]=static_cast<int>(buckets[peer].size());
		}
		outgoing.reserve(total);
		for(const auto& bucket:buckets)outgoing.insert(outgoing.end(),bucket.begin(),bucket.end());
	});
	MPI_Alltoall(counts.data(),1,MPI_INT,received.counts.data(),1,MPI_INT,comm);
	CollectiveLocalStage(comm,"point identity receive allocation",[&] {
		std::size_t total=0;
		for(int peer=0;peer<ranks;++peer) {
			if(received.counts[peer]<0)throw std::runtime_error("negative point identity MPI count");
			received.offsets[peer]=static_cast<int>(total);Charge(total,static_cast<std::size_t>(received.counts[peer]),cap);
		}
		received.bytes.resize(total);
	});
	MPI_Alltoallv(outgoing.data(),counts.data(),offsets.data(),MPI_BYTE,
		received.bytes.data(),received.counts.data(),received.offsets.data(),MPI_BYTE,comm);
	return received;
}
}

// Exact full-key equality determines equivalence; a hash only routes packets.
// Unique occurrence ids must be independent of rank layout. The minimum id
// wins, and the returned rank identifies its source in the given communicator.
// No full-key catalog is gathered to root or replicated on every rank. Skewed
// buckets can hit the explicit wire cap and then all live ranks reject.
inline std::vector<PointIdentityRepresentative> ResolveDistributedPointIdentities(MPI_Comm comm,
	const std::vector<PointIdentityOccurrence>& local,PointIdentityLimits limits={})
{
	using namespace point_identity_detail;
	int ranks=0;MPI_Comm_size(comm,&ranks);
	std::vector<std::string> buckets;
	CollectiveLocalStage(comm,"point identity input validation",[&] {
		if(!limits.max_wire_bytes||limits.max_wire_bytes>static_cast<std::size_t>(std::numeric_limits<int>::max())
			||!limits.max_local_occurrences||local.size()>limits.max_local_occurrences)
			throw std::invalid_argument("invalid point identity limits");
		buckets.resize(ranks);std::size_t used=0;
		for(const auto& value:local) {
			if(value.occurrence<0||value.key.empty())throw std::invalid_argument("invalid point identity occurrence");
			Charge(used,8,limits.max_wire_bytes);
			Put(buckets[static_cast<std::uint64_t>(value.occurrence)%ranks],static_cast<std::uint64_t>(value.occurrence));
		}
	});
	{
		const auto ids=Exchange(comm,buckets,limits.max_wire_bytes);
		CollectiveLocalStage(comm,"point identity occurrence uniqueness",[&] {
			std::set<std::uint64_t> unique;
			for(int peer=0;peer<ranks;++peer) {
				auto packet=ids.Peer(peer);
				while(!packet.empty())if(!unique.insert(Get(packet)).second)
					throw std::runtime_error("duplicate global point occurrence");
			}
		});
	}
	CollectiveLocalStage(comm,"point identity key routing",[&] {
		buckets.assign(ranks,{});std::size_t used=0;
		for(const auto& value:local) {
			Charge(used,16,limits.max_wire_bytes);Charge(used,value.key.size(),limits.max_wire_bytes);
			auto& bucket=buckets[BucketHash(value.key)%ranks];Put(bucket,static_cast<std::uint64_t>(value.occurrence));
			Put(bucket,value.key.size());bucket.append(value.key);
		}
	});
	{
		const auto keys=Exchange(comm,buckets,limits.max_wire_bytes);
		CollectiveLocalStage(comm,"point identity representative election",[&] {
			std::map<std::string,PointIdentityRepresentative> representatives;
			for(int peer=0;peer<ranks;++peer) {
				auto packet=keys.Peer(peer);
				while(!packet.empty()) {
					const auto id=static_cast<std::int64_t>(Get(packet));const std::string key(Key(packet));
					const auto found=representatives.find(key);
					if(found==representatives.end())representatives.emplace(key,PointIdentityRepresentative{id,peer});
					else if(id<found->second.occurrence)found->second={id,peer};
				}
			}
			buckets.assign(ranks,{});std::size_t used=0;
			for(int peer=0;peer<ranks;++peer) {
				auto packet=keys.Peer(peer);
				while(!packet.empty()) {
					const auto id=Get(packet);const auto& chosen=representatives.at(std::string(Key(packet)));
					Charge(used,24,limits.max_wire_bytes);Put(buckets[peer],id);
					Put(buckets[peer],static_cast<std::uint64_t>(chosen.occurrence));Put(buckets[peer],static_cast<std::uint64_t>(chosen.rank));
				}
			}
		});
	}
	const auto replies=Exchange(comm,buckets,limits.max_wire_bytes);
	std::vector<PointIdentityRepresentative> result;
	CollectiveLocalStage(comm,"point identity response validation",[&] {
		std::map<std::int64_t,std::size_t> indices;
		for(std::size_t i=0;i<local.size();++i)indices.emplace(local[i].occurrence,i);
		result.resize(local.size());std::vector<bool> seen(local.size(),false);
		for(int peer=0;peer<ranks;++peer) {
			auto packet=replies.Peer(peer);
			while(!packet.empty()) {
				const auto id=Get(packet),chosen=Get(packet),owner=Get(packet);
				const auto found=indices.find(static_cast<std::int64_t>(id));
				if(found==indices.end()||seen[found->second]||chosen>id||owner>=static_cast<std::uint64_t>(ranks))
					throw std::runtime_error("invalid point identity response");
				seen[found->second]=true;result[found->second]={static_cast<std::int64_t>(chosen),static_cast<int>(owner)};
			}
		}
		for(bool value:seen)if(!value)throw std::runtime_error("missing point identity response");
	});
	return result;
}
} // namespace iga
#endif
