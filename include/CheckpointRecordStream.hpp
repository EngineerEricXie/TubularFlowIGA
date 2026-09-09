#ifndef IGA_CHECKPOINT_RECORD_STREAM_HPP
#define IGA_CHECKPOINT_RECORD_STREAM_HPP
#include "CheckpointMetadataCodec.hpp"
#include "CheckpointWordStream.hpp"

namespace iga {
// A history can exceed the metadata limit; each individual record is bounded.
// The bundle supplies the total record count. Records carry uint64 byte lengths.
template<class Records> std::uint64_t CheckpointRecordBytes(const Records& records)
{
	std::uint64_t bytes = 0;
	for (const auto& record : records) bytes = checkpoint_stream::Add(bytes, checkpoint_stream::Add(8, SerializeGraphHistoryRecord(record).size()));
	return bytes;
}
template<class Records, class Sink> void WriteCheckpointRecords(const Records& records, Sink sink)
{
	for (const auto& record : records) {
		const auto bytes = SerializeGraphHistoryRecord(record);
		checkpoint_stream::Writer length(sink); length.Word(bytes.size()); length.Finish();
		for (std::size_t offset = 0; offset < bytes.size(); offset += 65536)
			sink(bytes.data()+offset, std::min(std::size_t{65536}, bytes.size()-offset));
	}
}
template<class Consumer> class CheckpointRecordReader {
public:
	CheckpointRecordReader(std::uint64_t count, Consumer consumer) : expected_(count), consumer_(std::move(consumer)) {}
	void Consume(const void* data, std::size_t size)
	{
		checkpoint_metadata::Require(!failed_, "history reader already failed");
		try {
			checkpoint_metadata::Require(!size || data, "null history bytes");
			const auto* bytes = static_cast<const unsigned char*>(data);
			while (size) {
				checkpoint_metadata::Require(records_ < expected_, "trailing history records");
				if (header_ < 8) {
					length_ |= std::uint64_t(*bytes++) << (8*header_++); --size;
					if (header_ < 8) continue;
					checkpoint_metadata::Require(length_ > 0 && length_ <= checkpoint_metadata::maximum_bytes, "invalid history record length");
				}
				const auto count = std::min(size, static_cast<std::size_t>(length_)-record_.size());
				record_.append(reinterpret_cast<const char*>(bytes), count); bytes += count; size -= count;
				if (record_.size() == length_) {
					consumer_(records_, std::string_view(record_)); ++records_; record_.clear(); header_ = 0; length_ = 0;
				}
			}
		} catch (...) { failed_ = true; throw; }
	}
	void Finish() const
	{
		checkpoint_metadata::Require(!failed_ && records_ == expected_ && header_ == 0 && record_.empty(), "incomplete history records");
	}
private:
	std::uint64_t expected_ = 0, records_ = 0, length_ = 0;
	unsigned header_ = 0;
	bool failed_ = false;
	Consumer consumer_;
	std::string record_;
};
} // namespace iga
#endif
