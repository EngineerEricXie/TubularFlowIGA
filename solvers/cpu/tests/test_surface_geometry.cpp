#include "SurfaceGeometry.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace {

using iga::RawSurfaceSoup;
using iga::RawSurfaceTriangle;

static_assert(!std::is_default_constructible<iga::ClosedTriangulatedSurface>::value,
	"closed surfaces must be created through validated canonical construction");

bool Near(double left, double right, double tolerance = 1.0e-12)
{
	return std::fabs(left-right) <= tolerance;
}

template <class Function>
void RequireRejected(Function&& function)
{
	bool rejected = false;
	try {
		function();
	} catch (const std::exception&) {
		rejected = true;
	}
	assert(rejected);
}

RawSurfaceTriangle Face(std::int64_t a, std::int64_t b, std::int64_t c, std::int64_t id = 7)
{
	RawSurfaceTriangle result;
	result.indices = {{a, b, c}};
	result.boundary_id = id;
	return result;
}

RawSurfaceSoup Tetrahedron()
{
	RawSurfaceSoup result;
	result.vertices = {{{{0.0, 0.0, 0.0}}, {{1.0, 0.0, 0.0}}, {{0.0, 1.0, 0.0}}, {{0.0, 0.0, 1.0}}}};
	result.triangles = {Face(0, 2, 1), Face(0, 1, 3), Face(0, 3, 2), Face(1, 2, 3)};
	return result;
}

RawSurfaceSoup Cube()
{
	RawSurfaceSoup result;
	result.vertices = {{{{0, 0, 0}}, {{1, 0, 0}}, {{1, 1, 0}}, {{0, 1, 0}}, {{0, 0, 1}}, {{1, 0, 1}}, {{1, 1, 1}}, {{0, 1, 1}}}};
	result.triangles = {Face(0, 2, 1), Face(0, 3, 2), Face(4, 5, 6), Face(4, 6, 7),
		Face(0, 1, 5), Face(0, 5, 4), Face(1, 2, 6), Face(1, 6, 5),
		Face(2, 3, 7), Face(2, 7, 6), Face(3, 0, 4), Face(3, 4, 7)};
	return result;
}

void CheckOutwardNormals(const iga::ClosedTriangulatedSurface& surface)
{
	std::array<double, 3> interior{{0.0, 0.0, 0.0}};
	for (const auto& vertex : surface.Vertices()) for (std::size_t axis = 0; axis < 3; ++axis)
		interior[axis] += vertex[axis]/static_cast<double>(surface.Vertices().size());
	for (const auto& triangle : surface.Triangles()) {
		std::array<double, 3> center{{0.0, 0.0, 0.0}};
		for (const auto index : triangle.indices) for (std::size_t axis = 0; axis < 3; ++axis)
			center[axis] += surface.Vertices()[index][axis]/3.0;
		const double dot = (center[0]-interior[0])*triangle.outward_unit_normal[0]
			+(center[1]-interior[1])*triangle.outward_unit_normal[1]
			+(center[2]-interior[2])*triangle.outward_unit_normal[2];
		assert(dot > 0.0);
	}
}

} // namespace

int main()
{
	const auto tetra = iga::ClosedTriangulatedSurface::Build(Tetrahedron());
	assert(Near(tetra.Diagnostics().volume_m3, 1.0/6.0));
	assert(tetra.Diagnostics().unique_edge_count == 6);
	assert(tetra.Diagnostics().unique_edge_count == 3*tetra.Triangles().size()/2);
	CheckOutwardNormals(tetra);
	constexpr double tiny_length_m = 0x1p-300;
	auto tiny_metres = Tetrahedron();
	for (auto& vertex : tiny_metres.vertices) for (auto& coordinate : vertex) coordinate *= tiny_length_m;
	iga::SurfaceValidationOptions tiny_unit_options;
	tiny_unit_options.length_scale_to_m = tiny_length_m;
	const auto tiny_direct = iga::ClosedTriangulatedSurface::Build(tiny_metres);
	const auto tiny_from_units = iga::ClosedTriangulatedSurface::Build(Tetrahedron(), tiny_unit_options);
	assert(tiny_direct.Vertices() == tiny_from_units.Vertices());
	assert(tiny_direct.CanonicalSha256() == tiny_from_units.CanonicalSha256());
	assert(std::isfinite(tiny_direct.Diagnostics().area_m2) && tiny_direct.Diagnostics().area_m2 > 0.0);
	assert(std::isfinite(tiny_direct.Diagnostics().signed_volume_m3)
		&& tiny_direct.Diagnostics().signed_volume_m3 > 0.0);
	assert(std::isfinite(tiny_direct.Diagnostics().volume_m3) && tiny_direct.Diagnostics().volume_m3 > 0.0);
	auto subdouble_volume = Tetrahedron();
	for (auto& vertex : subdouble_volume.vertices) for (auto& coordinate : vertex) coordinate *= 0x1p-500;
	RequireRejected([&] { iga::ClosedTriangulatedSurface::Build(subdouble_volume); });

	const auto cube = iga::ClosedTriangulatedSurface::Build(Cube());
	assert(Near(cube.Diagnostics().area_m2, 6.0));
	assert(Near(cube.Diagnostics().volume_m3, 1.0));
	assert(cube.Diagnostics().unique_edge_count == 3*cube.Triangles().size()/2);
	assert(cube.Diagnostics().component_count == 1);
	assert(cube.Diagnostics().vertex_fan_count == cube.Vertices().size());
	assert(cube.Diagnostics().boundary_id_histogram.at(7) == 12);
	CheckOutwardNormals(cube);

	auto inward = Cube();
	for (auto& face : inward.triangles) {
		std::swap(face.indices[1], face.indices[2]);
		std::rotate(face.indices.begin(), face.indices.begin()+1, face.indices.end());
	}
	std::reverse(inward.vertices.begin(), inward.vertices.end());
	for (auto& face : inward.triangles) for (auto& index : face.indices) index = 7-index;
	std::reverse(inward.triangles.begin(), inward.triangles.end());
	const auto canonical_inward = iga::ClosedTriangulatedSurface::Build(inward);
	assert(canonical_inward.Diagnostics().flipped_inward_shell);
	assert(canonical_inward.CanonicalSha256() == cube.CanonicalSha256());
	assert(canonical_inward.Vertices() == cube.Vertices());

	auto millimetres = Cube();
	for (auto& vertex : millimetres.vertices) for (auto& coordinate : vertex) coordinate *= 1000.0;
	iga::SurfaceValidationOptions mm_options;
	mm_options.length_scale_to_m = 0.001;
	const auto millimetre_cube = iga::ClosedTriangulatedSurface::Build(millimetres, mm_options);
	assert(millimetre_cube.CanonicalSha256() == cube.CanonicalSha256());
	assert(millimetre_cube.Vertices() == cube.Vertices());
	auto exact_duplicate_vertex = Tetrahedron();
	exact_duplicate_vertex.vertices.push_back(exact_duplicate_vertex.vertices[0]);
	exact_duplicate_vertex.triangles[0].indices[0] = 4;
	const auto collapsed_duplicate = iga::ClosedTriangulatedSurface::Build(exact_duplicate_vertex);
	assert(collapsed_duplicate.Vertices().size() == tetra.Vertices().size());
	assert(collapsed_duplicate.CanonicalSha256() == tetra.CanonicalSha256());

	auto changed_label = Cube();
	changed_label.triangles.front().boundary_id = 8;
	const auto changed = iga::ClosedTriangulatedSurface::Build(changed_label);
	assert(changed.CanonicalSha256() != cube.CanonicalSha256());
	assert(changed.Diagnostics().boundary_id_histogram.at(8) == 1);
	auto default_label = Tetrahedron();
	default_label.triangles.front().boundary_id = -1;
	iga::SurfaceValidationOptions default_options;
	default_options.default_boundary_id = 13;
	assert(iga::ClosedTriangulatedSurface::Build(default_label, default_options)
		.Diagnostics().boundary_id_histogram.at(13) == 1);

	auto nonfinite = Tetrahedron(); nonfinite.vertices[0][0] = std::numeric_limits<double>::quiet_NaN();
	RequireRejected([&] { iga::ClosedTriangulatedSurface::Build(nonfinite); });
	iga::SurfaceValidationOptions zero_scale; zero_scale.length_scale_to_m = 0.0;
	RequireRejected([&] { iga::ClosedTriangulatedSurface::Build(Tetrahedron(), zero_scale); });
	iga::SurfaceValidationOptions nan_scale; nan_scale.length_scale_to_m = std::numeric_limits<double>::quiet_NaN();
	RequireRejected([&] { iga::ClosedTriangulatedSurface::Build(Tetrahedron(), nan_scale); });
	auto scaled_overflow = Tetrahedron(); scaled_overflow.vertices[1][0] = 2.0;
	iga::SurfaceValidationOptions overflow_scale;
	overflow_scale.length_scale_to_m = std::numeric_limits<double>::max();
	RequireRejected([&] { iga::ClosedTriangulatedSurface::Build(scaled_overflow, overflow_scale); });
	auto bad_index = Tetrahedron(); bad_index.triangles[0].indices[0] = 99;
	RequireRejected([&] { iga::ClosedTriangulatedSurface::Build(bad_index); });
	auto negative_index = Tetrahedron(); negative_index.triangles[0].indices[0] = -1;
	RequireRejected([&] { iga::ClosedTriangulatedSurface::Build(negative_index); });
	auto bad_label = Tetrahedron(); bad_label.triangles[0].boundary_id = -2;
	RequireRejected([&] { iga::ClosedTriangulatedSurface::Build(bad_label); });
	auto repeated = Tetrahedron(); repeated.triangles[0].indices[1] = 0;
	RequireRejected([&] { iga::ClosedTriangulatedSurface::Build(repeated); });
	auto degenerate = Tetrahedron(); degenerate.vertices[2] = {{2.0, 0.0, 0.0}};
	RequireRejected([&] { iga::ClosedTriangulatedSurface::Build(degenerate); });
	auto duplicate = Tetrahedron(); duplicate.triangles.push_back(Face(1, 0, 2, 99));
	RequireRejected([&] { iga::ClosedTriangulatedSurface::Build(duplicate); });
	auto open = Tetrahedron(); open.triangles.pop_back();
	RequireRejected([&] { iga::ClosedTriangulatedSurface::Build(open); });
	auto nonmanifold = Tetrahedron(); nonmanifold.vertices.push_back({{2.0, 0.0, 0.0}}); nonmanifold.triangles.push_back(Face(0, 1, 4));
	RequireRejected([&] { iga::ClosedTriangulatedSurface::Build(nonmanifold); });
	auto disconnected = Cube(); const auto second = Tetrahedron();
	for (const auto vertex : second.vertices) disconnected.vertices.push_back({{vertex[0]+3.0, vertex[1], vertex[2]}});
	for (auto face : second.triangles) { for (auto& index : face.indices) index += 8; disconnected.triangles.push_back(face); }
	RequireRejected([&] { iga::ClosedTriangulatedSurface::Build(disconnected); });
	auto local_flip = Cube(); std::swap(local_flip.triangles[0].indices[1], local_flip.triangles[0].indices[2]);
	RequireRejected([&] { iga::ClosedTriangulatedSurface::Build(local_flip); });

	// Two closed tetrahedra touching at one vertex are a pinched, non-manifold vertex fan.
	auto pinched = Tetrahedron();
	pinched.vertices.push_back({{-1.0, 0.0, 0.0}});
	pinched.vertices.push_back({{0.0, -1.0, 0.0}});
	pinched.vertices.push_back({{0.0, 0.0, -1.0}});
	pinched.triangles.push_back(Face(0, 4, 5));
	pinched.triangles.push_back(Face(0, 6, 4));
	pinched.triangles.push_back(Face(0, 5, 6));
	pinched.triangles.push_back(Face(4, 6, 5));
	RequireRejected([&] { iga::ClosedTriangulatedSurface::Build(pinched); });

	RawSurfaceSoup zero_volume;
	zero_volume.vertices = {{{{0, 0, 0}}, {{0.2, 0.3, 0}}, {{1, 0, 0}}, {{0, 1, 0}}, {{-1, 0, 0}}, {{0, -1, 0}}}};
	zero_volume.triangles = {Face(0, 2, 3), Face(0, 3, 4), Face(0, 4, 5), Face(0, 5, 2),
		Face(1, 3, 2), Face(1, 4, 3), Face(1, 5, 4), Face(1, 2, 5)};
	RequireRejected([&] { iga::ClosedTriangulatedSurface::Build(zero_volume); });
	iga::SurfaceValidationOptions overflow; overflow.max_triangles = 3;
	RequireRejected([&] { iga::ClosedTriangulatedSurface::Build(Tetrahedron(), overflow); });
	iga::SurfaceValidationOptions invalid_options; invalid_options.default_boundary_id = -1;
	RequireRejected([&] { iga::ClosedTriangulatedSurface::Build(Tetrahedron(), invalid_options); });
	const auto above_boundary_id = static_cast<std::int64_t>(std::numeric_limits<std::uint32_t>::max())+1;
	auto oversized_label = Tetrahedron(); oversized_label.triangles[0].boundary_id = above_boundary_id;
	RequireRejected([&] { iga::ClosedTriangulatedSurface::Build(oversized_label); });
	iga::SurfaceValidationOptions oversized_default; oversized_default.default_boundary_id = above_boundary_id;
	RequireRejected([&] { iga::ClosedTriangulatedSurface::Build(Tetrahedron(), oversized_default); });

	std::cout << "surface geometry tests passed\n";
}
