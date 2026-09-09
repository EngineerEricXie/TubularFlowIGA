#ifndef IGA_CHECKPOINT_WORD_STREAM_HPP
#define IGA_CHECKPOINT_WORD_STREAM_HPP

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace iga {
namespace checkpoint_stream {
inline void Require(bool condition, const char* message)
{
	if (!condition) throw std::runtime_error(std::string("checkpoint stream: ")+message);
}
inline std::uint64_t RealBits(double value)
{
	static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559, "checkpoint requires IEEE binary64");
	Require(std::isfinite(value), "nonfinite field"); std::uint64_t bits; std::memcpy(&bits, &value, 8); return bits;
}
inline double Real(std::uint64_t bits)
{
	double value; std::memcpy(&value, &bits, 8); Require(std::isfinite(value), "nonfinite field"); return value;
}
inline std::uint64_t Bytes(std::uint64_t words)
{
	Require(words <= UINT64_MAX/8, "field byte count overflows"); return words*8;
}
inline std::uint64_t Add(std::uint64_t first, std::uint64_t second)
{
	Require(second <= UINT64_MAX-first, "field word count overflows"); return first+second;
}

template<class Sink> class Writer {
public:
	explicit Writer(Sink sink) : sink_(std::move(sink)) {}
	void Word(std::uint64_t value)
	{
		Require(!finished_, "writer is finished");
		for (unsigned i = 0; i < 8; ++i) buffer_[used_++] = static_cast<char>((value>>(8*i)) & 255);
		if (used_ == buffer_.size()) Flush();
	}
	void Real(double value) { Word(RealBits(value)); }
	void Finish() { Require(!finished_, "writer is finished"); Flush(); finished_ = true; }
private:
	void Flush()
	{
		if (!used_) return;
		try { sink_(buffer_.data(), used_); used_ = 0; }
		catch (...) { finished_ = true; throw; }
	}
	Sink sink_;
	std::array<char, 65536> buffer_;
	std::size_t used_ = 0;
	bool finished_ = false;
};

// Consume arbitrary byte fragments into an UNPUBLISHED candidate. The caller
// must also verify the bundle SHA before publishing it. No payload buffer is
// accumulated here; a consumer failure makes the reader permanently unusable.
template<class Consumer> class Reader {
public:
	Reader(std::uint64_t words, Consumer consumer) : expected_(Bytes(words)), consumer_(std::move(consumer)) {}
	Reader(const Reader&) = delete;
	Reader& operator=(const Reader&) = delete;
	void Consume(const void* data, std::size_t count)
	{
		Require(!failed_, "reader already failed");
		try {
			Require(count <= expected_-consumed_, "trailing field bytes");
			Require(!count || data, "null field input");
			const auto* bytes = static_cast<const unsigned char*>(data);
			for (std::size_t i = 0; i < count; ++i) {
				bits_ |= std::uint64_t(bytes[i])<<(8*partial_); ++partial_; ++consumed_;
				if (partial_ == 8) { consumer_(words_++, bits_); bits_ = 0; partial_ = 0; }
			}
		} catch (...) { failed_ = true; throw; }
	}
	void Finish() const { Require(!failed_ && consumed_ == expected_ && partial_ == 0, "incomplete field stream"); }
private:
	std::uint64_t expected_ = 0, consumed_ = 0, words_ = 0, bits_ = 0;
	Consumer consumer_;
	unsigned partial_ = 0;
	bool failed_ = false;
};

} // namespace checkpoint_stream
} // namespace iga
#endif
