#include "SurfaceSpatialIndex.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {

iga::RawSurfaceTriangle Face(std::int64_t a, std::int64_t b, std::int64_t c)
{
	iga::RawSurfaceTriangle result; result.indices = {{a, b, c}}; result.boundary_id = 7; return result;
}
iga::RawSurfaceSoup Cube()
{
	iga::RawSurfaceSoup result;
	result.vertices = {{{{0,0,0}}, {{1,0,0}}, {{1,1,0}}, {{0,1,0}}, {{0,0,1}}, {{1,0,1}}, {{1,1,1}}, {{0,1,1}}}};
	result.triangles = {Face(0,2,1), Face(0,3,2), Face(4,5,6), Face(4,6,7), Face(0,1,5), Face(0,5,4), Face(1,2,6), Face(1,6,5), Face(2,3,7), Face(2,7,6), Face(3,0,4), Face(3,4,7)};
	return result;
}
iga::RawSurfaceSoup Tetrahedron()
{
	iga::RawSurfaceSoup result;
	result.vertices = {{{{0,0,0}}, {{1,0,0}}, {{0,1,0}}, {{0,0,1}}}};
	result.triangles = {Face(0,2,1), Face(0,1,3), Face(0,3,2), Face(1,2,3)};
	return result;
}
iga::SurfaceAabb Box(double x0, double y0, double z0, double x1, double y1, double z1)
{ iga::SurfaceAabb result; result.minimum = {{x0,y0,z0}}; result.maximum = {{x1,y1,z1}}; return result; }
template <class Function> void RequireRejected(Function&& function) { bool rejected = false; try { function(); } catch (const std::exception&) { rejected = true; } assert(rejected); }

} // namespace

int main()
{
	const auto cube = iga::ClosedTriangulatedSurface::Build(Cube());
	iga::SurfaceSpatialIndex index(cube);
	assert(index.Surface().CanonicalSha256() == cube.CanonicalSha256());
	assert(index.LocatePoint({{0.5, 0.5, 0.5}}) == iga::PointLocation::Inside);
	assert(index.LocatePoint({{2.0, 0.5, 0.5}}) == iga::PointLocation::Outside);
	assert(index.LocatePoint({{0.0, 0.5, 0.5}}) == iga::PointLocation::Boundary);
	assert(index.LocatePoint({{0.0, 0.0, 0.5}}) == iga::PointLocation::Boundary);
	assert(index.LocatePoint({{0.0, 0.0, 0.0}}) == iga::PointLocation::Boundary);
	const iga::SurfaceSpatialIndex tetra(iga::ClosedTriangulatedSurface::Build(Tetrahedron()));
	assert(tetra.LocatePoint({{0.1,0.1,0.1}}) == iga::PointLocation::Inside);
	assert(tetra.LocatePoint({{0.7,0.7,0.1}}) == iga::PointLocation::Outside);
	assert(tetra.LocatePoint({{0.2,0.2,0.6}}) == iga::PointLocation::Boundary);
	// The oblique face x+y+z=1 slices this box without containing any tetra
	// vertex, exercising the edge-cross-box SAT axes rather than vertex tests.
	assert(tetra.IntersectBox(Box(0.30,0.30,0.30,0.40,0.40,0.40)) == iga::BoxContact::Interior);
	// The lower corner lies exactly on the same oblique face; shifting it by one
	// representable step separates all points of the closed box from the tetra.
	assert(tetra.IntersectBox(Box(0.40,0.30,0.30,0.50,0.40,0.40)) == iga::BoxContact::BoundaryOnly);
	assert(tetra.IntersectBox(Box(std::nextafter(0.40,1.0),0.30,0.30,0.50,0.40,0.40)) == iga::BoxContact::None);
	assert(tetra.IntersectBox(Box(0.40,0.30,0.30,0.50,0.40,0.30)) == iga::BoxContact::BoundaryOnly);
	// The box overlaps the tetra AABB and the oblique face's supporting plane,
	// but an edge×box SAT axis separates it from the triangular face.
	assert(tetra.IntersectBox(Box(0.80,0.30,-0.20,0.90,0.40,0.10)) == iga::BoxContact::None);
	assert(tetra.IntersectBox(Box(0.70,0.20,-0.05,0.85,0.35,0.10)) == iga::BoxContact::Interior);
	assert(index.IntersectBox(Box(0.2,0.2,0.2,0.8,0.8,0.8)) == iga::BoxContact::None);
	assert(index.IntersectBox(Box(-0.1,0.2,0.2,0.1,0.8,0.8)) == iga::BoxContact::Interior);
	assert(index.IntersectBox(Box(1.0,0.2,0.2,1.0,0.8,0.8)) == iga::BoxContact::BoundaryOnly);
	assert(index.IntersectBox(Box(1.0,1.0,1.0,1.0,1.0,1.0)) == iga::BoxContact::BoundaryOnly);
	assert(index.IntersectBox(Box(std::nextafter(1.0, 2.0),0.2,0.2,2.0,0.8,0.8)) == iga::BoxContact::None);
	iga::exact_dyadic::SetTestMagnitudeCap(1);
	std::vector<std::size_t> ambiguous_ids;
	assert(index.IntersectBox(Box(-0.1,0.2,0.2,0.1,0.8,0.8), &ambiguous_ids) == iga::BoxContact::Ambiguous);
	assert(!ambiguous_ids.empty()); for (std::size_t i = 1; i < ambiguous_ids.size(); ++i) assert(ambiguous_ids[i-1] < ambiguous_ids[i]);
	iga::exact_dyadic::SetTestMagnitudeCap(128);
	std::vector<std::size_t> ids; assert(index.IntersectBox(Box(-0.1,0.2,0.2,0.1,0.8,0.8), &ids) == iga::BoxContact::Interior); for (std::size_t i = 1; i < ids.size(); ++i) assert(ids[i-1] < ids[i]);
	RequireRejected([&] { index.LocatePoint({{std::numeric_limits<double>::quiet_NaN(), 0, 0}}); });
	RequireRejected([&] { index.IntersectBox(Box(1,0,0,0,1,1)); });
	auto millimetres = Cube(); for (auto& v : millimetres.vertices) for (auto& x : v) x *= 1000.0;
	iga::SurfaceValidationOptions options; options.length_scale_to_m = 0.001;
	const iga::SurfaceSpatialIndex scaled(iga::ClosedTriangulatedSurface::Build(millimetres, options));
	assert(scaled.Surface().CanonicalSha256() == index.Surface().CanonicalSha256());
	assert(scaled.LocatePoint({{0.5,0.5,0.5}}) == iga::PointLocation::Inside);
	auto translated_soup = Cube(); for (auto& vertex : translated_soup.vertices) { vertex[0] += 1.0e6; vertex[1] -= 2.0e6; vertex[2] += 3.0e6; }
	const iga::SurfaceSpatialIndex translated(iga::ClosedTriangulatedSurface::Build(translated_soup));
	assert(translated.LocatePoint({{1000000.5,-1999999.5,3000000.5}}) == iga::PointLocation::Inside);
	assert(translated.IntersectBox(Box(999999.9,-1999999.8,3000000.2,1000000.1,-1999999.2,3000000.8)) == iga::BoxContact::Interior);
	assert(index.Diagnostics().node_count > 0 && index.Diagnostics().narrow_phase_count > 0);
	const auto before_box = index.Diagnostics().box_candidate_count;
	index.LocatePoint({{0.25,0.25,0.25}});
	assert(index.Diagnostics().ray_candidate_count > 0 && index.Diagnostics().box_candidate_count == before_box);
	iga::SurfaceSpatialIndex moved(std::move(index));
	assert(moved.LocatePoint({{0.5,0.5,0.5}}) == iga::PointLocation::Inside);
	iga::SurfaceSpatialIndex assigned(iga::ClosedTriangulatedSurface::Build(Cube())); assigned = std::move(moved);
	assert(assigned.LocatePoint({{2.0,0.5,0.5}}) == iga::PointLocation::Outside);
	std::cout << "surface spatial index tests passed\n";
}
