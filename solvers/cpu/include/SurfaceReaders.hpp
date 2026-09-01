#ifndef IGA_SURFACE_READERS_HPP
#define IGA_SURFACE_READERS_HPP

#include "SurfaceGeometry.hpp"

#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

class SurfaceReaders {
public:
	static ClosedTriangulatedSurface ReadStlBuffer(const std::vector<std::uint8_t>& bytes,
		const SurfaceValidationOptions& options = SurfaceValidationOptions(), const std::string& context = "buffer")
	{
		try {
			return ClosedTriangulatedSurface::Build(IsExactBinaryStl(bytes)
				? ParseBinary(bytes, options) : ParseAscii(bytes, options), options);
		} catch (const std::exception& error) {
			throw std::invalid_argument("STL " + context + ": " + error.what());
		}
	}

	static ClosedTriangulatedSurface ReadStlText(const std::string& text,
		const SurfaceValidationOptions& options = SurfaceValidationOptions(), const std::string& context = "text")
	{
		return ReadStlBuffer(std::vector<std::uint8_t>(text.begin(), text.end()), options, context);
	}

	static ClosedTriangulatedSurface ReadStlPath(const std::string& path,
		const SurfaceValidationOptions& options = SurfaceValidationOptions())
	{
		std::ifstream input(path, std::ios::binary);
		if (!input) throw std::invalid_argument("STL " + path + ": cannot open file");
		input.seekg(0, std::ios::end);
		const auto size = input.tellg();
		if (size < 0) throw std::invalid_argument("STL " + path + ": cannot determine file size");
		if (static_cast<std::uintmax_t>(size) > std::numeric_limits<std::size_t>::max()
			|| static_cast<std::uintmax_t>(size) > std::numeric_limits<std::streamsize>::max())
			throw std::invalid_argument("STL " + path + ": file is too large");
		input.seekg(0, std::ios::beg);
		std::vector<std::uint8_t> bytes(static_cast<std::size_t>(size));
		if (!bytes.empty()) input.read(reinterpret_cast<char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
		if (!input && !bytes.empty()) throw std::invalid_argument("STL " + path + ": cannot read file");
		return ReadStlBuffer(bytes, options, path);
	}

private:
	static std::uint32_t ReadLe32(const std::vector<std::uint8_t>& bytes, std::size_t offset)
	{
		return static_cast<std::uint32_t>(bytes[offset]) | (static_cast<std::uint32_t>(bytes[offset+1])<<8)
			| (static_cast<std::uint32_t>(bytes[offset+2])<<16) | (static_cast<std::uint32_t>(bytes[offset+3])<<24);
	}
	static std::uint16_t ReadLe16(const std::vector<std::uint8_t>& bytes, std::size_t offset)
	{
		return static_cast<std::uint16_t>(bytes[offset]) | static_cast<std::uint16_t>(bytes[offset+1]<<8);
	}
	static float ReadLeFloat(const std::vector<std::uint8_t>& bytes, std::size_t offset)
	{
		const auto bits = ReadLe32(bytes, offset); float value = 0.0F;
		static_assert(sizeof(bits) == sizeof(value), "unsupported float representation");
		static_assert(std::numeric_limits<float>::is_iec559, "binary STL requires IEEE-754 float32");
		std::memcpy(&value, &bits, sizeof(value)); return value;
	}
	static bool IsExactBinaryStl(const std::vector<std::uint8_t>& bytes)
	{
		constexpr std::size_t prefix = 84U, record = 50U;
		if (bytes.size() < prefix) return false;
		const auto count = static_cast<std::size_t>(ReadLe32(bytes, 80));
		return count <= (std::numeric_limits<std::size_t>::max()-prefix)/record
			&& bytes.size() == prefix+record*count;
	}
	static void ValidateCount(std::size_t count, const SurfaceValidationOptions& options)
	{
		if (options.max_triangles < 0) throw std::invalid_argument("configured maximum triangle count is negative");
		if (static_cast<std::uintmax_t>(count) > static_cast<std::uintmax_t>(options.max_triangles))
			throw std::invalid_argument("triangle count exceeds configured maximum");
	}
	static RawSurfaceSoup ParseBinary(const std::vector<std::uint8_t>& bytes, const SurfaceValidationOptions& options)
	{
		const auto count = static_cast<std::size_t>(ReadLe32(bytes, 80)); ValidateCount(count, options);
		if (count > std::numeric_limits<std::size_t>::max()/3U) throw std::overflow_error("vertex count overflows");
		RawSurfaceSoup soup; soup.vertices.reserve(3*count); soup.triangles.reserve(count);
		for (std::size_t face = 0, offset = 84; face < count; ++face, offset += 50) {
			for (std::size_t component = 0; component < 3; ++component)
				if (!std::isfinite(ReadLeFloat(bytes, offset+4*component))) throw std::invalid_argument("binary facet normal is nonfinite");
			RawSurfaceTriangle triangle;
			for (std::size_t corner = 0; corner < 3; ++corner) {
				std::array<double, 3> vertex{};
				for (std::size_t component = 0; component < 3; ++component) {
					const float value = ReadLeFloat(bytes, offset+12+12*corner+4*component);
					if (!std::isfinite(value)) throw std::invalid_argument("binary vertex is nonfinite");
					vertex[component] = value;
				}
				triangle.indices[corner] = static_cast<std::int64_t>(soup.vertices.size()); soup.vertices.push_back(vertex);
			}
			if (ReadLe16(bytes, offset+48) != 0) throw std::invalid_argument("binary attribute word must be zero");
			soup.triangles.push_back(triangle);
		}
		return soup;
	}
	static double ParseFinite(const std::string& token)
	{
		errno = 0; char* end = nullptr; const double value = std::strtod(token.c_str(), &end);
		if (end != token.c_str()+token.size() || errno == ERANGE || !std::isfinite(value))
			throw std::invalid_argument("ASCII numeric token is invalid");
		return value;
	}
	static RawSurfaceSoup ParseAscii(const std::vector<std::uint8_t>& bytes, const SurfaceValidationOptions& options)
	{
		if (bytes.size() > 64U*1024U*1024U) throw std::invalid_argument("ASCII STL exceeds parser limit");
		const std::string text(bytes.begin(), bytes.end()); std::istringstream input(text); std::vector<std::string> tokens;
		for (std::string token; input >> token;) { if (tokens.size() == 10000000U) throw std::invalid_argument("ASCII STL has too many tokens"); tokens.push_back(token); }
		std::size_t at = 0;
		auto need = [&]() -> const std::string& { if (at == tokens.size()) throw std::invalid_argument("ASCII STL is truncated"); return tokens[at++]; };
		auto expect = [&](const char* word) { if (need() != word) throw std::invalid_argument(std::string("expected '")+word+"'"); };
		expect("solid");
		if (at < tokens.size() && tokens[at] != "facet" && tokens[at] != "endsolid") ++at;
		RawSurfaceSoup soup;
		while (at < tokens.size() && tokens[at] == "facet") {
			expect("facet"); expect("normal"); for (int i = 0; i < 3; ++i) (void)ParseFinite(need());
			expect("outer"); expect("loop"); RawSurfaceTriangle triangle;
			for (std::size_t corner = 0; corner < 3; ++corner) {
				expect("vertex"); std::array<double, 3> vertex{{ParseFinite(need()), ParseFinite(need()), ParseFinite(need())}};
				triangle.indices[corner] = static_cast<std::int64_t>(soup.vertices.size()); soup.vertices.push_back(vertex);
			}
			expect("endloop"); expect("endfacet"); soup.triangles.push_back(triangle);
			ValidateCount(soup.triangles.size(), options);
		}
		expect("endsolid");
		if (at < tokens.size()) ++at;
		if (at != tokens.size()) throw std::invalid_argument("ASCII STL has trailing tokens");
		return soup;
	}
};

} // namespace iga

#endif
