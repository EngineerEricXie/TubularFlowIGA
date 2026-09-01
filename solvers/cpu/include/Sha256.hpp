#ifndef IGA_SHA256_HPP
#define IGA_SHA256_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>

namespace iga {

class Sha256 {
public:
	void Append(const void* data, std::size_t count)
	{
		if (count && !data) throw std::invalid_argument("SHA-256 input is null");
		constexpr auto maximum_bytes = std::numeric_limits<std::uint64_t>::max()/8;
		if (bytes_ > maximum_bytes || count > maximum_bytes-bytes_)
			throw std::overflow_error("SHA-256 length overflows");
		const auto* input = static_cast<const std::uint8_t*>(data);
		bytes_ += count;
		while (count) {
			const auto copied = std::min<std::size_t>(count, block_.size()-used_);
			std::memcpy(block_.data()+used_, input, copied);
			used_ += copied;
			input += copied;
			count -= copied;
			if (used_ == block_.size()) {
				Transform();
				used_ = 0;
			}
		}
	}

	void AppendLittleEndian32(std::uint32_t value)
	{
		std::array<std::uint8_t, 4> bytes{};
		for (int i = 0; i < 4; ++i)
			bytes[static_cast<std::size_t>(i)] = static_cast<std::uint8_t>(value>>(8*i));
		Append(bytes.data(), bytes.size());
	}

	void AppendLittleEndian64(std::uint64_t value)
	{
		std::array<std::uint8_t, 8> bytes{};
		for (int i = 0; i < 8; ++i)
			bytes[static_cast<std::size_t>(i)] = static_cast<std::uint8_t>(value>>(8*i));
		Append(bytes.data(), bytes.size());
	}

	void AppendNormalizedDouble(double value)
	{
		if (!std::isfinite(value))
			throw std::invalid_argument("SHA-256 canonical double must be finite");
		if (value == 0.0) value = 0.0;
		std::uint64_t bits = 0;
		static_assert(sizeof(bits) == sizeof(value), "unsupported double representation");
		std::memcpy(&bits, &value, sizeof(bits));
		AppendLittleEndian64(bits);
	}

	std::array<std::uint8_t, 32> Digest() const
	{
		auto copy = *this;
		copy.block_[copy.used_++] = 0x80;
		if (copy.used_ > 56) {
			while (copy.used_ < copy.block_.size()) copy.block_[copy.used_++] = 0;
			copy.Transform();
			copy.used_ = 0;
		}
		while (copy.used_ < 56) copy.block_[copy.used_++] = 0;
		const auto bits = copy.bytes_*8;
		for (int i = 7; i >= 0; --i)
			copy.block_[copy.used_++] = static_cast<std::uint8_t>(bits>>(8*i));
		copy.Transform();
		std::array<std::uint8_t, 32> result{};
		for (int i = 0; i < 8; ++i)
			for (int j = 0; j < 4; ++j)
				result[static_cast<std::size_t>(4*i+j)] =
					static_cast<std::uint8_t>(copy.state_[static_cast<std::size_t>(i)]>>(24-8*j));
		return result;
	}

	std::string Hex() const
	{
		static constexpr char hexadecimal[] = "0123456789abcdef";
		std::string result;
		result.reserve(64);
		for (const auto byte : Digest()) {
			result.push_back(hexadecimal[byte>>4]);
			result.push_back(hexadecimal[byte&0x0f]);
		}
		return result;
	}

private:
	static std::uint32_t RotateRight(std::uint32_t value, int shift)
	{
		return (value>>shift) | (value<<(32-shift));
	}

	void Transform()
	{
		static constexpr std::array<std::uint32_t, 64> constants{{
			0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1,
			0x923f82a4, 0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
			0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786,
			0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
			0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147,
			0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
			0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
			0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
			0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a,
			0x5b9cca4f, 0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
			0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2}};
		std::array<std::uint32_t, 64> words{};
		for (int i = 0; i < 16; ++i) {
			const auto offset = static_cast<std::size_t>(4*i);
			words[static_cast<std::size_t>(i)] =
				(static_cast<std::uint32_t>(block_[offset])<<24)
				| (static_cast<std::uint32_t>(block_[offset+1])<<16)
				| (static_cast<std::uint32_t>(block_[offset+2])<<8)
				| static_cast<std::uint32_t>(block_[offset+3]);
		}
		for (int i = 16; i < 64; ++i) {
			const auto x = words[static_cast<std::size_t>(i-15)];
			const auto y = words[static_cast<std::size_t>(i-2)];
			const auto sigma0 = RotateRight(x, 7)^RotateRight(x, 18)^(x>>3);
			const auto sigma1 = RotateRight(y, 17)^RotateRight(y, 19)^(y>>10);
			words[static_cast<std::size_t>(i)] = words[static_cast<std::size_t>(i-16)]
				+sigma0+words[static_cast<std::size_t>(i-7)]+sigma1;
		}
		auto a = state_[0], b = state_[1], c = state_[2], d = state_[3];
		auto e = state_[4], f = state_[5], g = state_[6], h = state_[7];
		for (int i = 0; i < 64; ++i) {
			const auto sum1 = RotateRight(e, 6)^RotateRight(e, 11)^RotateRight(e, 25);
			const auto choose = (e&f)^((~e)&g);
			const auto temporary1 = h+sum1+choose+constants[static_cast<std::size_t>(i)]
				+words[static_cast<std::size_t>(i)];
			const auto sum0 = RotateRight(a, 2)^RotateRight(a, 13)^RotateRight(a, 22);
			const auto majority = (a&b)^(a&c)^(b&c);
			const auto temporary2 = sum0+majority;
			h = g;
			g = f;
			f = e;
			e = d+temporary1;
			d = c;
			c = b;
			b = a;
			a = temporary1+temporary2;
		}
		state_[0] += a; state_[1] += b; state_[2] += c; state_[3] += d;
		state_[4] += e; state_[5] += f; state_[6] += g; state_[7] += h;
	}

	std::array<std::uint8_t, 64> block_{};
	std::size_t used_ = 0;
	std::uint64_t bytes_ = 0;
	std::array<std::uint32_t, 8> state_{{0x6a09e667, 0xbb67ae85, 0x3c6ef372,
		0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19}};
};

} // namespace iga

#endif
