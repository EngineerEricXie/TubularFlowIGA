#include "PrescribedSurfaceMotion.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using iga::PrescribedSurfaceFrame;
using iga::RawSurfaceSoup;
using iga::RawSurfaceTriangle;

using Evaluation = iga::PrescribedSurfaceMotion::Evaluation;
static_assert(std::is_same<decltype(std::declval<const Evaluation&>().Surface()),
	const iga::ClosedTriangulatedSurface&>::value, "evaluation surface must be read-only");
static_assert(std::is_same<decltype(std::declval<const Evaluation&>().SourceVerticesM()),
	const std::vector<std::array<double, 3>>&>::value, "evaluation positions must be read-only");
static_assert(std::is_same<decltype(std::declval<const Evaluation&>().SourceVertexVelocitiesMPerS()),
	const std::vector<std::array<double, 3>>&>::value, "evaluation velocities must be read-only");
static_assert(std::is_same<decltype(std::declval<const Evaluation&>().CanonicalTriangleProvenance()),
	const std::vector<iga::SourceTriangleProvenance>&>::value, "evaluation provenance must be read-only");
static_assert(!std::is_assignable<Evaluation&, Evaluation>::value,
	"validated evaluation state must not be replaceable");

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

RawSurfaceTriangle Face(std::int64_t a, std::int64_t b, std::int64_t c, std::int64_t label)
{
	RawSurfaceTriangle result;
	result.indices = {{a, b, c}};
	result.boundary_id = label;
	return result;
}

RawSurfaceSoup Tetrahedron(double shift_x)
{
	RawSurfaceSoup result;
	result.vertices = {{{{shift_x, 0.0, 0.0}}, {{shift_x+1.0, 0.0, 0.0}},
		{{shift_x, 1.0, 0.0}}, {{shift_x, 0.0, 1.0}}}};
	result.triangles = {Face(0, 2, 1, 17), Face(0, 1, 3, 23), Face(0, 3, 2, 17), Face(1, 2, 3, 23)};
	return result;
}

RawSurfaceSoup Cube()
{
	RawSurfaceSoup result;
	result.vertices = {{{{0, 0, 0}}, {{1, 0, 0}}, {{1, 1, 0}}, {{0, 1, 0}},
		{{0, 0, 1}}, {{1, 0, 1}}, {{1, 1, 1}}, {{0, 1, 1}}}};
	result.triangles = {Face(0, 2, 1, 1), Face(0, 3, 2, 1), Face(4, 5, 6, 1), Face(4, 6, 7, 1),
		Face(0, 1, 5, 1), Face(0, 5, 4, 1), Face(1, 2, 6, 1), Face(1, 6, 5, 1),
		Face(2, 3, 7, 1), Face(2, 7, 6, 1), Face(3, 0, 4, 1), Face(3, 4, 7, 1)};
	return result;
}

std::vector<PrescribedSurfaceFrame> ThreeFrames()
{
	return {{0.0, Tetrahedron(0.0)}, {1.0, Tetrahedron(1.0)}, {2.0, Tetrahedron(3.0)}};
}

std::vector<PrescribedSurfaceFrame> DistinctVelocityFrames()
{
	auto end = Tetrahedron(0.0);
	end.vertices[0] = {{0.10, 0.00, 0.00}};
	end.vertices[1] = {{1.00, 0.20, 0.00}};
	end.vertices[2] = {{0.00, 1.00, 0.30}};
	end.vertices[3] = {{0.20, 0.00, 1.00}};
	return {{0.0, Tetrahedron(0.0)}, {1.0, end}};
}

RawSurfaceSoup RotatedCube()
{
	auto result = Cube();
	for (auto& vertex : result.vertices) {
		vertex[0] = 1.0-vertex[0];
		vertex[1] = 1.0-vertex[1];
	}
	return result;
}

} // namespace

int main()
{
	const iga::PrescribedSurfaceMotion motion(ThreeFrames());
	const auto midpoint = motion.Evaluate(0.5, 0.0, 0.5);
	assert(Near(midpoint.SourceVerticesM()[0][0], 0.5));
	assert(Near(midpoint.SourceVerticesM()[1][0], 1.5));
	assert(Near(midpoint.SourceVertexVelocitiesMPerS()[0][0], 1.0));
	const auto endpoint = motion.Evaluate(1.0, 0.5, 1.0);
	assert(Near(endpoint.SourceVerticesM()[0][0], 1.0));
	assert(Near(endpoint.SourceVertexVelocitiesMPerS()[0][0], 1.0));
	const auto next_interval = motion.Evaluate(1.0, 1.0, 1.5);
	assert(Near(next_interval.SourceVertexVelocitiesMPerS()[0][0], 2.0));
	RequireRejected([&] { motion.Evaluate(0.75, 0.25, 1.25); });

	const auto repeat = motion.Evaluate(0.5, 0.0, 0.5);
	assert(midpoint.IdentitySha256() == repeat.IdentitySha256());
	assert(midpoint.IdentitySha256() == "52ba4de4c7f05fb58c8663dbe6025ea41eec5adfec05c4d586f1ce2682539d42");
	assert(midpoint.SourceVerticesM() == repeat.SourceVerticesM());
	assert(midpoint.SourceVertexVelocitiesMPerS() == repeat.SourceVertexVelocitiesMPerS());
	assert(midpoint.Surface().CanonicalSha256() == repeat.Surface().CanonicalSha256());
	assert(midpoint.CanonicalTriangleProvenance() == repeat.CanonicalTriangleProvenance());
	assert(midpoint.CanonicalTriangleProvenance().size() == midpoint.Surface().Triangles().size());
	for (const auto& provenance : midpoint.CanonicalTriangleProvenance()) {
		assert(provenance.source_triangle < motion.SourceTriangles().size());
		assert(provenance.source_vertex_indices == motion.SourceTriangles()[provenance.source_triangle].source_vertex_indices);
		assert(provenance.boundary_id == motion.SourceTriangles()[provenance.source_triangle].boundary_id);
	}
	const auto wall_velocity = midpoint.WallVelocity(0, {{0.2, 0.3, 0.5}});
	assert(Near(wall_velocity[0], 1.0) && Near(wall_velocity[1], 0.0) && Near(wall_velocity[2], 0.0));
	RequireRejected([&] { midpoint.WallVelocity(static_cast<std::uint32_t>(midpoint.Surface().Triangles().size()), {{0.2, 0.3, 0.5}}); });
	RequireRejected([&] { midpoint.WallVelocity(0, {{-1.0e-6, 0.3, 0.700001}}); });
	RequireRejected([&] { midpoint.WallVelocity(0, {{0.2, 0.3, 0.6}}); });
	RequireRejected([&] { midpoint.WallVelocity(0, {{std::numeric_limits<double>::quiet_NaN(), 0.3, 0.5}}); });

	const iga::PrescribedSurfaceMotion distinct_motion(DistinctVelocityFrames());
	const auto distinct = distinct_motion.Evaluate(0.5, 0.0, 0.5);
	std::size_t cyclic_canonical_triangle = distinct.Surface().Triangles().size();
	for (std::size_t triangle = 0; triangle < distinct.CanonicalTriangleProvenance().size(); ++triangle)
		if (distinct.CanonicalTriangleProvenance()[triangle].canonical_corner_to_source_corner
			!= std::array<std::uint32_t, 3>{{0, 1, 2}}) { cyclic_canonical_triangle = triangle; break; }
	assert(cyclic_canonical_triangle < distinct.Surface().Triangles().size());
	const std::array<double, 3> canonical_barycentric{{0.2, 0.3, 0.5}};
	for (double weight : canonical_barycentric) assert(std::isfinite(weight) && weight >= 0.0 && weight <= 1.0);
	assert(Near(canonical_barycentric[0]+canonical_barycentric[1]+canonical_barycentric[2], 1.0));
	std::array<double, 3> expected_velocity{{0.0, 0.0, 0.0}};
	const auto& canonical_triangle = distinct.Surface().Triangles()[cyclic_canonical_triangle];
	for (std::size_t canonical_corner = 0; canonical_corner < 3; ++canonical_corner) {
		std::size_t source_vertex = distinct.SourceVerticesM().size();
		for (std::size_t vertex = 0; vertex < distinct.SourceVerticesM().size(); ++vertex)
			if (distinct.SourceVerticesM()[vertex] == distinct.Surface().Vertices()[canonical_triangle.indices[canonical_corner]]) {
				source_vertex = vertex; break;
			}
		assert(source_vertex < distinct.SourceVerticesM().size());
		for (std::size_t axis = 0; axis < 3; ++axis)
			expected_velocity[axis] += canonical_barycentric[canonical_corner]
				*distinct.SourceVertexVelocitiesMPerS()[source_vertex][axis];
	}
	const auto mapped_velocity = distinct.WallVelocity(static_cast<std::uint32_t>(cyclic_canonical_triangle), canonical_barycentric);
	for (std::size_t axis = 0; axis < 3; ++axis) {
		assert(std::isfinite(mapped_velocity[axis]));
		assert(Near(mapped_velocity[axis], expected_velocity[axis]));
	}

	const auto same_state_different_envelope = motion.Evaluate(0.5, 0.25, 0.75);
	assert(midpoint.IdentitySha256() == same_state_different_envelope.IdentitySha256());

	auto nonfinite_coordinate = ThreeFrames();
	nonfinite_coordinate[1].surface.vertices[0][0] = std::numeric_limits<double>::quiet_NaN();
	RequireRejected([&] { iga::PrescribedSurfaceMotion rejected(nonfinite_coordinate); });
	auto nonfinite_time = ThreeFrames(); nonfinite_time[1].time_s = std::numeric_limits<double>::infinity();
	RequireRejected([&] { iga::PrescribedSurfaceMotion rejected(nonfinite_time); });
	auto time_order = ThreeFrames(); time_order[1].time_s = 0.0;
	RequireRejected([&] { iga::PrescribedSurfaceMotion rejected(time_order); });
	auto overflowing_frame_duration = ThreeFrames();
	overflowing_frame_duration[0].time_s = -std::numeric_limits<double>::max();
	overflowing_frame_duration[1].time_s = std::numeric_limits<double>::max();
	RequireRejected([&] { iga::PrescribedSurfaceMotion rejected(overflowing_frame_duration); });
	RequireRejected([&] { motion.Evaluate(0.5, -std::numeric_limits<double>::max(), std::numeric_limits<double>::max()); });
	auto bad_connectivity = ThreeFrames(); std::swap(bad_connectivity[1].surface.triangles[0].indices[0], bad_connectivity[1].surface.triangles[0].indices[1]);
	RequireRejected([&] { iga::PrescribedSurfaceMotion rejected(bad_connectivity); });
	auto bad_label = ThreeFrames(); bad_label[1].surface.triangles[0].boundary_id = 99;
	RequireRejected([&] { iga::PrescribedSurfaceMotion rejected(bad_label); });
	auto bad_winding = ThreeFrames(); std::swap(bad_winding[1].surface.triangles[0].indices[1], bad_winding[1].surface.triangles[0].indices[2]);
	RequireRejected([&] { iga::PrescribedSurfaceMotion rejected(bad_winding); });
	auto degenerate = ThreeFrames(); degenerate[1].surface.vertices[2] = {{2.0, 0.0, 0.0}};
	RequireRejected([&] { iga::PrescribedSurfaceMotion rejected(degenerate); });

	std::vector<PrescribedSurfaceFrame> inward{{0.0, Cube()}, {1.0, Cube()}};
	for (auto& frame : inward) for (auto& face : frame.surface.triangles) std::swap(face.indices[1], face.indices[2]);
	RequireRejected([&] { iga::PrescribedSurfaceMotion rejected(inward); });
	std::vector<PrescribedSurfaceFrame> self_intersecting{{0.0, Cube()}, {1.0, Cube()}};
	self_intersecting[1].surface.vertices[6] = {{1.0, 1.0, -1.0}};
	RequireRejected([&] { iga::PrescribedSurfaceMotion rejected(self_intersecting); });

	iga::PrescribedSurfaceMotionOptions containment;
	containment.require_containment_in_fixed_bounds = true;
	containment.fixed_bounds_m.minimum = {{0.0, 0.0, 0.0}};
	containment.fixed_bounds_m.maximum = {{1.25, 1.0, 1.0}};
	RequireRejected([&] { iga::PrescribedSurfaceMotion rejected(ThreeFrames(), containment); });
	iga::PrescribedSurfaceMotionOptions displacement;
	displacement.maximum_displacement_m = 0.5;
	RequireRejected([&] { iga::PrescribedSurfaceMotion rejected(ThreeFrames(), displacement); });
	iga::PrescribedSurfaceMotionOptions velocity;
	velocity.maximum_velocity_m_per_s = 0.5;
	RequireRejected([&] { iga::PrescribedSurfaceMotion rejected(ThreeFrames(), velocity); });
	iga::PrescribedSurfaceMotionOptions extension_band;
	extension_band.extension_band_m = std::numeric_limits<double>::min();
	extension_band.maximum_extension_band_cfl = 1.0;
	RequireRejected([&] { iga::PrescribedSurfaceMotion rejected(ThreeFrames(), extension_band); });
	std::vector<PrescribedSurfaceFrame> denormal_displacement{{0.0, Tetrahedron(0.0)}, {1.0, Tetrahedron(0.0)}};
	denormal_displacement[1].surface.vertices[0][0] = std::numeric_limits<double>::denorm_min();
	extension_band.extension_band_m = 2.0;
	extension_band.maximum_extension_band_cfl = 0.0;
	RequireRejected([&] { iga::PrescribedSurfaceMotion rejected(denormal_displacement, extension_band); });

	std::vector<PrescribedSurfaceFrame> invalid_interpolation{{0.0, Cube()}, {1.0, RotatedCube()}};
	const iga::PrescribedSurfaceMotion interpolation_motion(invalid_interpolation);
	RequireRejected([&] { interpolation_motion.Evaluate(0.5, 0.0, 1.0); });

	std::cout << "prescribed surface motion tests passed\n";
	return 0;
}
