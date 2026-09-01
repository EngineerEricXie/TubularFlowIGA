#include "CartesianDomainClassification.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <algorithm>

namespace {
iga::RawSurfaceTriangle Face(std::int64_t a, std::int64_t b, std::int64_t c) { iga::RawSurfaceTriangle r; r.indices={{a,b,c}}; r.boundary_id=7; return r; }
iga::RawSurfaceSoup Cube(double lo, double hi)
{
	iga::RawSurfaceSoup s; s.vertices={{{{lo,lo,lo}},{{hi,lo,lo}},{{hi,hi,lo}},{{lo,hi,lo}},{{lo,lo,hi}},{{hi,lo,hi}},{{hi,hi,hi}},{{lo,hi,hi}}}};
	s.triangles={Face(0,2,1),Face(0,3,2),Face(4,5,6),Face(4,6,7),Face(0,1,5),Face(0,5,4),Face(1,2,6),Face(1,6,5),Face(2,3,7),Face(2,7,6),Face(3,0,4),Face(3,4,7)}; return s;
}
template <class F> void Reject(F&& f) { bool yes=false; try { f(); } catch (const std::exception&) { yes=true; } assert(yes); }
}
int main()
{
	iga::CubicCartesianGridSpec spec{{{0,0,0}},{{1,1,1}},{{6,6,6}}};
	iga::CubicCartesianBackground background(spec);
	iga::CartesianDomainClassification catalog(background, iga::SurfaceSpatialIndex(iga::ClosedTriangulatedSurface::Build(Cube(.25,.75))));
	assert(catalog.Cells().size() == 216 && catalog.Diagnostics().inside_count == 8 && catalog.Diagnostics().cut_count == 56 && catalog.Diagnostics().outside_count == 152);
	assert(catalog.ActiveCellIds().size() == 64 && catalog.ActiveCellIds().front() == 43);
	for (std::size_t i=1;i<catalog.ActiveCellIds().size();++i) assert(catalog.ActiveCellIds()[i-1] < catalog.ActiveCellIds()[i]);
	for (std::size_t position=0; position<catalog.Cells().size(); ++position) { const auto& record=catalog.Cells()[position]; assert(record.id==position && record.index[0]+6*(record.index[1]+6*record.index[2])==position); for (std::size_t i=1;i<record.triangle_ids.size();++i) assert(record.triangle_ids[i-1]<record.triangle_ids[i]); }
	for (const auto id : catalog.ActiveCellIds()) assert(catalog.Cells()[id].classification != iga::CellClassification::Outside);
	iga::CartesianDomainClassification eight(iga::CubicCartesianBackground({{{0,0,0}},{{1,1,1}},{{8,8,8}}}), iga::SurfaceSpatialIndex(iga::ClosedTriangulatedSurface::Build(Cube(.25,.75))));
	assert(eight.Diagnostics().inside_count==8 && eight.Diagnostics().cut_count==208 && eight.Diagnostics().outside_count==296 && eight.ActiveCellIds().size()==216);
	iga::CartesianDomainClassification translated(iga::CubicCartesianBackground({{{10,10,10}},{{11,11,11}},{{6,6,6}}}), iga::SurfaceSpatialIndex(iga::ClosedTriangulatedSurface::Build(Cube(10.25,10.75))));
	assert(translated.Diagnostics().inside_count==8 && translated.Diagnostics().cut_count==56 && translated.Diagnostics().outside_count==152);
	auto permuted=Cube(.25,.75); std::reverse(permuted.vertices.begin(),permuted.vertices.end()); for(auto& f:permuted.triangles) for(auto& index:f.indices) index=7-index; std::reverse(permuted.triangles.begin(),permuted.triangles.end());
	iga::CartesianDomainClassification permuted_catalog(background,iga::SurfaceSpatialIndex(iga::ClosedTriangulatedSurface::Build(permuted))); assert(permuted_catalog.SurfaceCanonicalHash()==catalog.SurfaceCanonicalHash());
	for (std::uint64_t id=0; id<background.ElementCount(); ++id) { const auto cell=background.Cell(id); assert(cell.lower_m[0] == background.Plane(0,cell.index[0]) && cell.upper_m[0] == background.Plane(0,cell.index[0]+1)); }
	iga::CubicCartesianGridSpec coarse{{{0,0,0}},{{1,1,1}},{{2,2,2}}};
	iga::CartesianDomainClassification closed_plane(iga::CubicCartesianBackground(coarse), iga::SurfaceSpatialIndex(iga::ClosedTriangulatedSurface::Build(Cube(.5,1.0))));
	assert(closed_plane.Diagnostics().cut_count == 8 && closed_plane.Diagnostics().inside_count == 0 && closed_plane.Diagnostics().outside_count == 0);
	assert(closed_plane.ActiveCellIds().size() == 8);
	for (std::size_t i=0;i<closed_plane.Cells().size();++i) { const auto& cell=closed_plane.Cells()[i]; assert(cell.id==i && cell.classification==iga::CellClassification::Cut && cell.boundary_only_contact && !cell.triangle_ids.empty()); if (i) assert(closed_plane.ActiveCellIds()[i-1] == i-1); }
	assert(background.Plane(0,0) == 0.0 && background.Plane(0,6) == 1.0 && background.Cell(0).upper_m[0] == background.Cell(1).lower_m[0]);
	const auto final_element = background.MaterializeElement(background.ElementCount()-1);
	assert((final_element.bezier_points[63] == std::array<double,3>{{1.0,1.0,1.0}}));
	Reject([&] { iga::CubicCartesianBackground({{{1.0,0,0}},{{std::nextafter(1.0,2.0),1,1}},{{2,1,1}}}); });
	{
		struct ResetCap { ~ResetCap() { iga::exact_dyadic::SetTestMagnitudeCap(128); } } reset;
		iga::CubicCartesianBackground ambiguous_background(spec);
		iga::SurfaceSpatialIndex ambiguous_surface(iga::ClosedTriangulatedSurface::Build(Cube(.25,.75)));
		iga::exact_dyadic::SetTestMagnitudeCap(1);
		iga::CartesianDomainClassification ambiguous(std::move(ambiguous_background), std::move(ambiguous_surface));
		assert(ambiguous.Diagnostics().ambiguous_count > 0);
		bool retained = false; for (const auto& record : ambiguous.Cells()) if (record.ambiguous) { assert(record.classification == iga::CellClassification::Cut); if (!record.triangle_ids.empty()) { retained = true; for (std::size_t i=1;i<record.triangle_ids.size();++i) assert(record.triangle_ids[i-1] < record.triangle_ids[i]); } }
		assert(retained);
	}
	Reject([&] { iga::CartesianDomainClassification bad(iga::CubicCartesianBackground({{{0,0,0}},{{.5,.5,.5}},{{2,2,2}}}), iga::SurfaceSpatialIndex(iga::ClosedTriangulatedSurface::Build(Cube(.25,.75)))); });
	std::cout << "cartesian domain classification tests passed\n";
}
