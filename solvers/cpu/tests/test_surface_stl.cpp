#include "SurfaceReaders.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

template <class Function>
void RequireRejected(Function&& function)
{
	bool rejected = false;
	try { function(); } catch (const std::exception&) { rejected = true; }
	assert(rejected);
}

void Append32(std::vector<std::uint8_t>& bytes, std::uint32_t value)
{
	for (int shift = 0; shift < 32; shift += 8) bytes.push_back(static_cast<std::uint8_t>(value>>shift));
}

void AppendFloat(std::vector<std::uint8_t>& bytes, float value)
{
	std::uint32_t bits = 0; std::memcpy(&bits, &value, sizeof(bits)); Append32(bytes, bits);
}

std::string TetraAscii()
{
	return "solid tetra\n"
		"facet normal 0 0 0 outer loop vertex 0 0 0 vertex 0 1 0 vertex 1 0 0 endloop endfacet\n"
		"facet normal 0 0 0 outer loop vertex 0 0 0 vertex 1 0 0 vertex 0 0 1 endloop endfacet\n"
		"facet normal 0 0 0 outer loop vertex 0 0 0 vertex 0 0 1 vertex 0 1 0 endloop endfacet\n"
		"facet normal 0 0 0 outer loop vertex 1 0 0 vertex 0 1 0 vertex 0 0 1 endloop endfacet\n"
		"endsolid tetra\n";
}

std::vector<std::uint8_t> TetraBinary()
{
	std::vector<std::uint8_t> bytes(80, 0); std::memcpy(bytes.data(), "solid binary tetra", 18); Append32(bytes, 4);
	const float faces[4][3][3] = {{{0,0,0},{0,1,0},{1,0,0}}, {{0,0,0},{1,0,0},{0,0,1}}, {{0,0,0},{0,0,1},{0,1,0}}, {{1,0,0},{0,1,0},{0,0,1}}};
	for (const auto& face : faces) { for (int i = 0; i < 3; ++i) AppendFloat(bytes, 0.0F); for (const auto& vertex : face) for (const auto coordinate : vertex) AppendFloat(bytes, coordinate); bytes.push_back(0); bytes.push_back(0); }
	return bytes;
}

using Point = std::array<float, 3>;
using Face = std::array<int, 3>;

const std::array<Point, 8> kCubePoints{{{{0,0,0}}, {{1,0,0}}, {{1,1,0}}, {{0,1,0}}, {{0,0,1}}, {{1,0,1}}, {{1,1,1}}, {{0,1,1}}}};
const std::array<Face, 12> kCubeFaces{{{{0,2,1}}, {{0,3,2}}, {{4,5,6}}, {{4,6,7}}, {{0,1,5}}, {{0,5,4}}, {{1,2,6}}, {{1,6,5}}, {{2,3,7}}, {{2,7,6}}, {{3,0,4}}, {{3,4,7}}}};

std::string AsciiMesh(const std::array<Point, 8>& points, const std::array<Face, 12>& faces)
{
	std::ostringstream output; output << "solid cube\n";
	for (const auto& face : faces) {
		output << "facet normal 0 0 0 outer loop ";
		for (const auto index : face) output << "vertex " << points[static_cast<std::size_t>(index)][0] << ' ' << points[static_cast<std::size_t>(index)][1] << ' ' << points[static_cast<std::size_t>(index)][2] << ' ';
		output << "endloop endfacet\n";
	}
	output << "endsolid cube\n"; return output.str();
}

std::vector<std::uint8_t> BinaryMesh(const std::array<Point, 8>& points, const std::array<Face, 12>& faces,
	bool inward = false, bool reverse_faces = false)
{
	std::vector<std::uint8_t> bytes(80, 0); std::memcpy(bytes.data(), "solid binary cube", 17); Append32(bytes, 12);
	for (std::size_t ordinal = 0; ordinal < faces.size(); ++ordinal) {
		const auto& face = faces[reverse_faces ? faces.size()-1-ordinal : ordinal];
		for (int i = 0; i < 3; ++i) AppendFloat(bytes, 0.0F);
		for (std::size_t corner = 0; corner < 3; ++corner) {
			const auto index = face[inward ? 2-corner : corner];
			for (const auto coordinate : points[static_cast<std::size_t>(index)]) AppendFloat(bytes, coordinate);
		}
		bytes.push_back(0); bytes.push_back(0);
	}
	return bytes;
}

} // namespace

int main()
{
	const auto ascii = iga::SurfaceReaders::ReadStlText(TetraAscii());
	const auto binary_bytes = TetraBinary();
	const auto binary = iga::SurfaceReaders::ReadStlBuffer(binary_bytes, {}, "solid-header-binary");
	assert(ascii.CanonicalSha256() == binary.CanonicalSha256());
	assert(std::fabs(binary.Diagnostics().volume_m3-1.0/6.0) < 1.0e-12);
	assert(binary.Diagnostics().input_triangle_count == 4);
	const auto cube_ascii = iga::SurfaceReaders::ReadStlText(AsciiMesh(kCubePoints, kCubeFaces));
	const auto cube_binary = iga::SurfaceReaders::ReadStlBuffer(BinaryMesh(kCubePoints, kCubeFaces));
	assert(cube_ascii.CanonicalSha256() == cube_binary.CanonicalSha256());
	assert(std::fabs(cube_ascii.Diagnostics().area_m2-6.0) < 1.0e-12);
	assert(std::fabs(cube_ascii.Diagnostics().volume_m3-1.0) < 1.0e-12);
	assert(iga::SurfaceReaders::ReadStlBuffer(BinaryMesh(kCubePoints, kCubeFaces, true, true)).CanonicalSha256()
		== cube_binary.CanonicalSha256());

	iga::SurfaceValidationOptions millimetres; millimetres.length_scale_to_m = 0.001;
	auto mm = TetraAscii();
	const std::array<std::pair<const char*, const char*>, 3> replacements{{{"1 0 0", "1000 0 0"}, {"0 1 0", "0 1000 0"}, {"0 0 1", "0 0 1000"}}};
	for (const auto& replacement : replacements) {
		std::size_t position = 0; while ((position = mm.find(replacement.first, position)) != std::string::npos) { mm.replace(position, std::strlen(replacement.first), replacement.second); position += std::strlen(replacement.second); }
	}
	assert(iga::SurfaceReaders::ReadStlText(mm, millimetres).CanonicalSha256() == ascii.CanonicalSha256());

	RequireRejected([] { iga::SurfaceReaders::ReadStlText("solid x facet normal 0 0 0 endsolid x"); });
	RequireRejected([] { iga::SurfaceReaders::ReadStlText("solid x facet normal nan 0 0 outer loop vertex 0 0 0 vertex 0 1 0 vertex 1 0 0 endloop endfacet endsolid x"); });
	auto truncated = binary_bytes; truncated.pop_back(); RequireRejected([&] { iga::SurfaceReaders::ReadStlBuffer(truncated); });
	auto trailing = binary_bytes; trailing.push_back(0); RequireRejected([&] { iga::SurfaceReaders::ReadStlBuffer(trailing); });
	auto bad_count = binary_bytes; bad_count[80] = 5; RequireRejected([&] { iga::SurfaceReaders::ReadStlBuffer(bad_count); });
	auto huge_count = binary_bytes; huge_count[80] = 0xff; huge_count[81] = 0xff; huge_count[82] = 0xff; huge_count[83] = 0xff; RequireRejected([&] { iga::SurfaceReaders::ReadStlBuffer(huge_count); });
	auto attribute = binary_bytes; attribute[84+48] = 1; RequireRejected([&] { iga::SurfaceReaders::ReadStlBuffer(attribute); });
	auto nonfinite_vertex = binary_bytes; const std::uint32_t nan = 0x7fc00000U; for (int i = 0; i < 4; ++i) nonfinite_vertex[84+12+i] = static_cast<std::uint8_t>(nan>>(8*i)); RequireRejected([&] { iga::SurfaceReaders::ReadStlBuffer(nonfinite_vertex); });
	auto nonfinite_normal = binary_bytes; for (int i = 0; i < 4; ++i) nonfinite_normal[84+i] = static_cast<std::uint8_t>(nan>>(8*i)); RequireRejected([&] { iga::SurfaceReaders::ReadStlBuffer(nonfinite_normal); });
	iga::SurfaceValidationOptions max; max.max_triangles = 3; RequireRejected([&] { iga::SurfaceReaders::ReadStlBuffer(binary_bytes, max); });
	RequireRejected([] { iga::SurfaceReaders::ReadStlText("solid x facet normal inf 0 0 outer loop vertex 0 0 0 vertex 0 1 0 vertex 1 0 0 endloop endfacet endsolid x"); });
	RequireRejected([] { iga::SurfaceReaders::ReadStlText("solid x endsolid x trailing"); });
	RequireRejected([] { iga::SurfaceReaders::ReadStlText("solid x facet normal 0 0 0 outer loop vertex 0 0 0 vertex 0 1 0 vertex 1 0 0 endloop endfacet endsolid x"); });
	const std::string path = "/tmp/tubularflowiga_surface_stl_test.stl";
	{ std::ofstream output(path, std::ios::binary); output.write(reinterpret_cast<const char*>(binary_bytes.data()), static_cast<std::streamsize>(binary_bytes.size())); }
	assert(iga::SurfaceReaders::ReadStlPath(path).CanonicalSha256() == binary.CanonicalSha256());
	assert(std::remove(path.c_str()) == 0);
	RequireRejected([] { iga::SurfaceReaders::ReadStlPath("/tmp/tubularflowiga_missing_surface.stl"); });
	std::cout << "surface STL tests passed\n";
}
