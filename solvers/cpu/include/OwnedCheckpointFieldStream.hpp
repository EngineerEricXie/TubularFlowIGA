#ifndef IGA_OWNED_CHECKPOINT_FIELD_STREAM_HPP
#define IGA_OWNED_CHECKPOINT_FIELD_STREAM_HPP

#include "OwnedCheckpointVector.hpp"
#include "CheckpointWordStream.hpp"
#include <map>

namespace iga {
inline std::uint64_t OwnedCheckpointFieldBytes(const OwnedCheckpointVector& field)
{
	return checkpoint_stream::Bytes(checkpoint_stream::Add(3, field.values.size()));
}
template<class Sink> void WriteOwnedCheckpointField(const OwnedCheckpointVector& field, Sink sink)
{
	checkpoint_stream::Writer<Sink> output(std::move(sink));
	output.Word(field.global_rows); output.Word(field.row_begin); output.Word(field.row_end);
	for (double value : field.values) output.Real(value);
	output.Finish();
}
inline auto ReadOwnedCheckpointField(OwnedCheckpointVector& field)
{
	return checkpoint_stream::Reader(checkpoint_stream::Add(3, field.values.size()), [&field](std::uint64_t index, std::uint64_t word) {
		if (index == 0) checkpoint_stream::Require(word == field.global_rows, "global rows differ");
		else if (index == 1) checkpoint_stream::Require(word == field.row_begin, "owned row begin differs");
		else if (index == 2) checkpoint_stream::Require(word == field.row_end, "owned row end differs");
		else field.values.at(static_cast<std::size_t>(index-3)) = checkpoint_stream::Real(word);
	});
}
// Read source shards in any order into an unpublished target partition.
// Each source callback streams one authenticated file; only target-owned rows
// are retained. Memory is O(target rows + source shard count), not global rows.
class OwnedCheckpointRepartitionReader {
public:
	explicit OwnedCheckpointRepartitionReader(OwnedCheckpointVector& target) : target_(target)
	{
		checkpoint_stream::Require(target.row_begin<=target.row_end&&target.row_end<=target.global_rows
			&&target.values.size()==target.row_end-target.row_begin,"invalid repartition target shape");
	}
	template<class Source> void ReadShard(std::uint64_t payload_bytes,Source source)
	{
		checkpoint_stream::Require(!failed_&&!finished_,"repartition reader is closed");
		try {
			checkpoint_stream::Require(payload_bytes>=24&&payload_bytes%8==0,"invalid owned shard byte count");
			const auto count=payload_bytes/8-3;
			checkpoint_stream::Require(count<=target_.global_rows,"owned shard exceeds global shape");
			std::uint64_t begin=0,end=0;
			checkpoint_stream::Reader reader(count+3,[&](std::uint64_t index,std::uint64_t word) {
				if(index==0)checkpoint_stream::Require(word==target_.global_rows,"repartition global rows differ");
				else if(index==1)begin=word;
				else if(index==2) {
					end=word;checkpoint_stream::Require(begin<=end&&end<=target_.global_rows&&end-begin==count,
						"repartition source range or size differs");
				} else {
					const double value=checkpoint_stream::Real(word);const auto row=begin+index-3;
					if(row>=target_.row_begin&&row<target_.row_end)target_.values[static_cast<std::size_t>(row-target_.row_begin)]=value;
				}
			});
			source([&](const void* bytes,std::size_t size) { reader.Consume(bytes,size); });reader.Finish();
			if(begin!=end) {
				const auto next=ranges_.lower_bound(begin);
				checkpoint_stream::Require(next==ranges_.end()||end<=next->first,"overlapping checkpoint source shards");
				checkpoint_stream::Require(next==ranges_.begin()||std::prev(next)->second<=begin,"overlapping checkpoint source shards");
				ranges_.emplace(begin,end);
			}
		} catch(...) { failed_=true;throw; }
	}
	void Finish()
	{
		checkpoint_stream::Require(!failed_&&!finished_,"repartition reader is closed");
		try {
			std::uint64_t next=0;
			for(const auto& range:ranges_) { checkpoint_stream::Require(range.first==next,"missing checkpoint source rows");next=range.second; }
			checkpoint_stream::Require(next==target_.global_rows,"incomplete checkpoint global coverage");finished_=true;
		} catch(...) { failed_=true;throw; }
	}
private:
	OwnedCheckpointVector& target_;std::map<std::uint64_t,std::uint64_t> ranges_;
	bool failed_=false,finished_=false;
};
} // namespace iga
#endif
