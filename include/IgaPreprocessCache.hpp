#ifndef IGA_PREPROCESS_CACHE_HPP
#define IGA_PREPROCESS_CACHE_HPP

#include <array>
#include <cstdint>
#include <fstream>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <string>

namespace igacache {

constexpr std::array<char, 8> kMagic{{'I', 'G', 'A', 'C', 'A', 'C', 'H', '1'}};
constexpr std::uint32_t kVersion = 1;
constexpr std::uint32_t kBezierPointCount = 64;

struct Header
{
	std::array<char, 8> magic;
	std::uint32_t version;
	std::uint64_t elements;
	std::uint64_t nodes;
	std::uint64_t mesh_hash;
	std::uint32_t bezier_points;
	std::uint32_t reserved;
};

template <class T>
inline void Write(std::ostream& out, const T& value)
{
	out.write(reinterpret_cast<const char*>(&value), sizeof(T));
	if (!out) throw std::runtime_error("preprocessing cache write failed");
}

template <class T>
inline void Read(std::istream& in, T& value)
{
	in.read(reinterpret_cast<char*>(&value), sizeof(T));
	if (!in) throw std::runtime_error("preprocessing cache is truncated");
}

inline std::uint64_t HashFile(const std::string& path)
{
	std::ifstream input(path, std::ios::binary);
	if (!input) throw std::runtime_error("cannot hash control mesh: " + path);
	std::uint64_t hash = 14695981039346656037ull;
	std::array<char, 65536> buffer{};
	while (input) {
		input.read(buffer.data(), buffer.size());
		const std::streamsize count = input.gcount();
		for (std::streamsize i = 0; i < count; ++i) {
			hash ^= static_cast<unsigned char>(buffer[static_cast<std::size_t>(i)]);
			hash *= 1099511628211ull;
		}
	}
	if (!input.eof()) throw std::runtime_error("failed while hashing control mesh: " + path);
	return hash;
}

inline void WriteHeader(
	std::ostream& out,
	std::uint64_t elements,
	std::uint64_t nodes,
	std::uint64_t mesh_hash)
{
	out.write(kMagic.data(), kMagic.size());
	Write(out, kVersion);
	Write(out, elements);
	Write(out, nodes);
	Write(out, mesh_hash);
	Write(out, kBezierPointCount);
	Write(out, std::uint32_t{0});
}

inline Header ReadHeader(std::istream& in)
{
	Header header{};
	in.read(header.magic.data(), header.magic.size());
	Read(in, header.version);
	Read(in, header.elements);
	Read(in, header.nodes);
	Read(in, header.mesh_hash);
	Read(in, header.bezier_points);
	Read(in, header.reserved);
	if (header.magic != kMagic) throw std::runtime_error("invalid preprocessing cache magic");
	if (header.version != kVersion) throw std::runtime_error("unsupported preprocessing cache version");
	if (header.elements == 0 || header.nodes == 0) throw std::runtime_error("empty preprocessing cache");
	if (header.bezier_points != kBezierPointCount)
		throw std::runtime_error("unsupported preprocessing cache Bezier point count");
	return header;
}

} // namespace igacache

#endif
