#include "CutCellVolumeQuadrature.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {

iga::RawSurfaceTriangle Face(std::int64_t a, std::int64_t b, std::int64_t c)
{
	iga::RawSurfaceTriangle result; result.indices = {{a, b, c}}; result.boundary_id = 7; return result;
}
iga::RawSurfaceSoup Cube(double lower, double upper)
{
	iga::RawSurfaceSoup result;
	result.vertices = {{{{lower,lower,lower}}, {{upper,lower,lower}}, {{upper,upper,lower}}, {{lower,upper,lower}},
		{{lower,lower,upper}}, {{upper,lower,upper}}, {{upper,upper,upper}}, {{lower,upper,upper}}}};
	result.triangles = {Face(0,2,1), Face(0,3,2), Face(4,5,6), Face(4,6,7), Face(0,1,5), Face(0,5,4),
		Face(1,2,6), Face(1,6,5), Face(2,3,7), Face(2,7,6), Face(3,0,4), Face(3,4,7)};
	return result;
}
iga::RawSurfaceSoup Tetrahedron()
{
	iga::RawSurfaceSoup result;
	result.vertices = {{{{0,0,0}}, {{1,0,0}}, {{0,1,0}}, {{0,0,1}}}};
	result.triangles = {Face(0,2,1), Face(0,1,3), Face(0,3,2), Face(1,2,3)};
	return result;
}
iga::RawSurfaceSoup Translate(iga::RawSurfaceSoup soup, double shift)
{
	for (auto& vertex : soup.vertices) for (double& coordinate : vertex) coordinate += shift;
	return soup;
}
iga::RawSurfaceSoup PermuteCube()
{
	auto result = Cube(0.0, 1.0);
	std::reverse(result.vertices.begin(), result.vertices.end());
	for (auto& triangle : result.triangles) for (auto& index : triangle.indices) index = 7-index;
	std::reverse(result.triangles.begin(), result.triangles.end());
	return result;
}
iga::RawSurfaceSoup PermuteTetrahedron()
{
	auto result = Tetrahedron();
	std::reverse(result.vertices.begin(), result.vertices.end());
	for (auto& triangle : result.triangles) for (auto& index : triangle.indices) index = 3-index;
	std::reverse(result.triangles.begin(), result.triangles.end());
	return result;
}
iga::RawSurfaceSoup Affine(iga::RawSurfaceSoup soup, const std::array<double, 3>& origin,
	const std::array<double, 3>& scale)
{
	for (auto& vertex : soup.vertices) for (std::size_t axis = 0; axis < 3; ++axis)
		vertex[axis] = origin[axis]+scale[axis]*vertex[axis];
	return soup;
}
iga::CartesianDomainClassification Domain(iga::RawSurfaceSoup soup, iga::CubicCartesianGridSpec spec)
{
	return iga::CartesianDomainClassification(iga::CubicCartesianBackground(spec),
		iga::SurfaceSpatialIndex(iga::ClosedTriangulatedSurface::Build(std::move(soup))));
}
bool Near(double left, double right, double tolerance = 2.0e-12)
{
	return std::abs(left-right) <= tolerance*std::max({1.0, std::abs(left), std::abs(right)});
}
template <class Function> void Reject(Function&& function)
{
	bool rejected = false; try { function(); } catch (const std::exception&) { rejected = true; } assert(rejected);
}
double Moment(const iga::VolumeQuadratureRule& rule, int x, int y, int z)
{
	double result = 0.0;
	for (const auto& point : rule.Points()) result += point.weight*std::pow(point.parametric[0], x)*std::pow(point.parametric[1], y)*std::pow(point.parametric[2], z);
	return result;
}
double WeightSum(const iga::VolumeQuadratureRule& rule)
{
	double result = 0.0;
	for (const auto& point : rule.Points()) result += point.weight;
	return result;
}
double PhysicalWeightSum(const iga::Element& element, const iga::VolumeQuadratureRule& rule)
{
	double result = 0.0;
	for (const auto& point : rule.Points()) result += point.weight*iga::EvaluateElementGeometry(element, point.parametric).raw_determinant;
	return result;
}
void RequireSameReferenceRule(const iga::CutCellVolumeQuadratureCatalog& first,
	const iga::CutCellVolumeQuadratureCatalog& second)
{
	const auto& left = first.Cell(0); const auto& right = second.Cell(0);
	assert(left.rule.Points().size() == right.rule.Points().size());
	assert(left.diagnostics.nodes == right.diagnostics.nodes && left.diagnostics.leaves == right.diagnostics.leaves
		&& left.diagnostics.reached_depth == right.diagnostics.reached_depth);
	for (std::size_t i = 0; i < left.rule.Points().size(); ++i)
		assert(left.rule.Points()[i].parametric == right.rule.Points()[i].parametric
			&& left.rule.Points()[i].weight == right.rule.Points()[i].weight);
}

} // namespace

int main()
{
	const iga::CubicCartesianGridSpec root{{{0,0,0}}, {{1,1,1}}, {{1,1,1}}};
	const auto cube_domain = Domain(Cube(0.0, 1.0), root);
	const iga::CutCellVolumeQuadratureCatalog cube(cube_domain, iga::OctreeCutQuadratureOptions{4,10000,10000,100000});
	const iga::CutCellVolumeQuadratureCatalog cube_repeat(cube_domain, iga::OctreeCutQuadratureOptions{4,10000,10000,100000});
	const auto& cube_cell = cube.Cell(0);
	assert(cube_cell.usable && cube_cell.rule.Points().size() == 64);
	assert(Near(cube_cell.diagnostics.estimated_reference_volume, 1.0));
	assert(Near(Moment(cube_cell.rule, 1, 0, 0), .5));
	assert(Near(Moment(cube_cell.rule, 2, 0, 0), 1.0/3.0));
	assert(Near(Moment(cube_cell.rule, 1, 1, 1), .125));
	for (const auto& point : cube_cell.rule.Points()) {
		assert(point.weight > 0.0 && point.parametric[0] >= 0.0 && point.parametric[0] <= 1.0);
	}
	cube.ValidateUsableRule(cube_domain, 0);
	assert(&cube.UsableRule(cube_domain, 0) == &cube_cell.rule);
	assert(Near(WeightSum(cube_cell.rule), cube_cell.diagnostics.estimated_reference_volume));
	assert(cube_repeat.Cell(0).diagnostics.nodes == cube_cell.diagnostics.nodes);
	assert(cube_repeat.Cell(0).rule.Points().size() == cube_cell.rule.Points().size());

	const auto inner_domain = Domain(Cube(.25, .75), root);
	const iga::CutCellVolumeQuadratureCatalog inner(inner_domain, iga::OctreeCutQuadratureOptions{5,500000,500000,3000000});
	assert(inner.Cell(0).diagnostics.unresolved_reference_volume == 0.0);
	assert(Near(inner.Cell(0).diagnostics.estimated_reference_volume, .125));
	assert(Near(Moment(inner.Cell(0).rule, 1, 0, 0), 1.0/16.0));
	assert(Near(Moment(inner.Cell(0).rule, 2, 0, 0), 13.0/384.0));
	assert(Near(Moment(inner.Cell(0).rule, 1, 1, 1), 1.0/64.0));
	assert(Near(WeightSum(inner.Cell(0).rule), inner.Cell(0).diagnostics.estimated_reference_volume));

	const iga::CubicCartesianGridSpec six{{{0,0,0}}, {{1,1,1}}, {{6,6,6}}};
	const auto mixed_domain = Domain(Cube(.25, .75), six);
	const iga::CutCellVolumeQuadratureCatalog mixed(mixed_domain, iga::OctreeCutQuadratureOptions{4,100000,100000,1000000});
	bool checked_inside = false, checked_outside = false;
	for (const auto& record : mixed.Cells()) {
		if (record.classification == iga::CellClassification::Inside) {
			const iga::FullCell4x4x4VolumeQuadratureProvider full(mixed_domain.Background().MaterializeElement(record.id));
			assert(record.rule.Points().size() == full.Rule().Points().size());
			for (std::size_t i = 0; i < record.rule.Points().size(); ++i) {
				assert(record.rule.Points()[i].parametric == full.Rule().Points()[i].parametric);
				assert(record.rule.Points()[i].weight == full.Rule().Points()[i].weight);
			}
			checked_inside = true;
		}
		if (record.classification == iga::CellClassification::Outside) {
			assert(record.rule.Points().empty() && record.usable);
			mixed.ValidateUsableRule(mixed_domain, record.id);
			checked_outside = true;
		}
	}
	assert(checked_inside && checked_outside);

	double previous_lower = -1.0, previous_upper = 2.0, previous_unresolved = 2.0, first_error = 0.0, final_error = 0.0;
	double first_moment_error = 0.0, final_moment_error = 0.0;
	for (std::uint32_t depth = 0; depth <= 5; ++depth) {
		const auto tetra_domain = Domain(Tetrahedron(), root);
		const iga::CutCellVolumeQuadratureCatalog tetra(tetra_domain, iga::OctreeCutQuadratureOptions{depth,200000,200000,2000000});
		const auto& record = tetra.Cell(0);
		const auto& d = record.diagnostics;
		assert(record.usable && d.lower_reference_volume >= previous_lower && d.upper_reference_volume <= previous_upper);
		assert(d.unresolved_reference_volume <= previous_unresolved && d.lower_reference_volume <= 1.0/6.0 && d.upper_reference_volume >= 1.0/6.0);
		const double error = std::abs(d.estimated_reference_volume-1.0/6.0);
		if (depth == 0) first_error = error;
		const double moment_error = std::abs(Moment(record.rule, 1, 0, 0)-1.0/24.0);
		assert(Near(Moment(record.rule, 1, 0, 0), Moment(record.rule, 0, 1, 0), 1.0e-10));
		assert(Near(Moment(record.rule, 1, 0, 0), Moment(record.rule, 0, 0, 1), 1.0e-10));
		if (depth == 0) first_moment_error = moment_error;
		if (depth == 5) { final_error = error; final_moment_error = moment_error; }
		assert(Near(WeightSum(record.rule), d.estimated_reference_volume));
		previous_lower = d.lower_reference_volume; previous_upper = d.upper_reference_volume; previous_unresolved = d.unresolved_reference_volume;
	}
	assert(final_error < first_error && final_moment_error < first_moment_error);
	const auto tetra_domain = Domain(Tetrahedron(), root);
	const iga::CutCellVolumeQuadratureCatalog tetra_reference(tetra_domain, iga::OctreeCutQuadratureOptions{5,200000,200000,2000000});
	const iga::CutCellVolumeQuadratureCatalog tetra_repeat(tetra_domain, iga::OctreeCutQuadratureOptions{5,200000,200000,2000000});
	const auto tetra_permuted_domain = Domain(PermuteTetrahedron(), root);
	const iga::CutCellVolumeQuadratureCatalog tetra_permuted(tetra_permuted_domain, iga::OctreeCutQuadratureOptions{5,200000,200000,2000000});
	RequireSameReferenceRule(tetra_reference, tetra_repeat);
	RequireSameReferenceRule(tetra_reference, tetra_permuted);
	assert(tetra_reference.SurfaceCanonicalHash() == tetra_permuted.SurfaceCanonicalHash());
	const iga::CubicCartesianGridSpec tetra_translated_grid{{{10,-3,7}}, {{11,-2,8}}, {{1,1,1}}};
	const auto tetra_translated_domain = Domain(Affine(Tetrahedron(), std::array<double,3>{{10,-3,7}}, std::array<double,3>{{1,1,1}}), tetra_translated_grid);
	const iga::CutCellVolumeQuadratureCatalog tetra_translated(tetra_translated_domain, iga::OctreeCutQuadratureOptions{5,200000,200000,2000000});
	const iga::CubicCartesianGridSpec tetra_scaled_grid{{{0,0,0}}, {{2,2,2}}, {{1,1,1}}};
	const auto tetra_scaled_domain = Domain(Affine(Tetrahedron(), std::array<double,3>{{0,0,0}}, std::array<double,3>{{2,2,2}}), tetra_scaled_grid);
	const iga::CutCellVolumeQuadratureCatalog tetra_scaled(tetra_scaled_domain, iga::OctreeCutQuadratureOptions{5,200000,200000,2000000});
	RequireSameReferenceRule(tetra_reference, tetra_translated);
	RequireSameReferenceRule(tetra_reference, tetra_scaled);
	const auto& tetra_diagnostics = tetra_reference.Cell(0).diagnostics;
	const auto& tetra_translated_diagnostics = tetra_translated.Cell(0).diagnostics;
	const auto& tetra_scaled_diagnostics = tetra_scaled.Cell(0).diagnostics;
	assert(Near(tetra_diagnostics.estimated_physical_volume, tetra_translated_diagnostics.estimated_physical_volume));
	assert(Near(tetra_diagnostics.lower_physical_volume, tetra_translated_diagnostics.lower_physical_volume));
	assert(Near(tetra_diagnostics.upper_physical_volume, tetra_translated_diagnostics.upper_physical_volume));
	assert(Near(tetra_scaled_diagnostics.estimated_physical_volume, 8.0*tetra_diagnostics.estimated_physical_volume));
	assert(Near(tetra_scaled_diagnostics.lower_physical_volume, 8.0*tetra_diagnostics.lower_physical_volume));
	assert(Near(tetra_scaled_diagnostics.upper_physical_volume, 8.0*tetra_diagnostics.upper_physical_volume));
	const iga::CubicCartesianGridSpec tetra_anisotropic_grid{{{4,-5,8}}, {{6,-2,12}}, {{1,1,1}}};
	const auto tetra_anisotropic_domain = Domain(Affine(Tetrahedron(), std::array<double,3>{{4,-5,8}}, std::array<double,3>{{2,3,4}}), tetra_anisotropic_grid);
	const iga::CutCellVolumeQuadratureCatalog tetra_anisotropic(tetra_anisotropic_domain, iga::OctreeCutQuadratureOptions{5,200000,200000,2000000});
	const auto anisotropic_element = tetra_anisotropic_domain.Background().MaterializeElement(0);
	for (const auto& point : tetra_anisotropic.Cell(0).rule.Points())
		assert(tetra_anisotropic_domain.SurfaceIndex().LocatePoint(iga::EvaluateElementGeometry(anisotropic_element, point.parametric).physical) == iga::PointLocation::Inside);
	assert(Near(PhysicalWeightSum(anisotropic_element, tetra_anisotropic.Cell(0).rule), tetra_anisotropic.Cell(0).diagnostics.estimated_physical_volume));

	const iga::CubicCartesianGridSpec twice{{{0,0,0}}, {{2,2,2}}, {{1,1,1}}};
	const auto scaled_domain = Domain(Cube(0.0, 2.0), twice);
	const iga::CutCellVolumeQuadratureCatalog scaled(scaled_domain, iga::OctreeCutQuadratureOptions{4,10000,10000,100000});
	assert(scaled.Cell(0).rule.Points().size() == cube.Cell(0).rule.Points().size());
	for (std::size_t i = 0; i < cube.Cell(0).rule.Points().size(); ++i) {
		assert(scaled.Cell(0).rule.Points()[i].parametric == cube.Cell(0).rule.Points()[i].parametric);
		assert(scaled.Cell(0).rule.Points()[i].weight == cube.Cell(0).rule.Points()[i].weight);
	}
	assert(Near(scaled.Cell(0).diagnostics.estimated_physical_volume, 8.0));
	const iga::CubicCartesianGridSpec translated_grid{{{10,10,10}}, {{11,11,11}}, {{1,1,1}}};
	const auto translated_domain = Domain(Translate(Cube(0.0, 1.0), 10.0), translated_grid);
	const iga::CutCellVolumeQuadratureCatalog translated(translated_domain, iga::OctreeCutQuadratureOptions{4,10000,10000,100000});
	assert(translated.Cell(0).rule.Points().size() == cube.Cell(0).rule.Points().size());
	for (std::size_t i = 0; i < cube.Cell(0).rule.Points().size(); ++i)
		assert(translated.Cell(0).rule.Points()[i].parametric == cube.Cell(0).rule.Points()[i].parametric
			&& translated.Cell(0).rule.Points()[i].weight == cube.Cell(0).rule.Points()[i].weight);
	assert(translated.Cell(0).diagnostics.nodes == cube.Cell(0).diagnostics.nodes
		&& translated.Cell(0).diagnostics.leaves == cube.Cell(0).diagnostics.leaves);
	const auto permuted_domain = Domain(PermuteCube(), root);
	const iga::CutCellVolumeQuadratureCatalog permuted(permuted_domain, iga::OctreeCutQuadratureOptions{4,10000,10000,100000});
	assert(permuted.SurfaceCanonicalHash() == cube.SurfaceCanonicalHash());
	assert(permuted.Cell(0).rule.Points().size() == cube.Cell(0).rule.Points().size());
	for (std::size_t i = 0; i < cube.Cell(0).rule.Points().size(); ++i)
		assert(permuted.Cell(0).rule.Points()[i].parametric == cube.Cell(0).rule.Points()[i].parametric
			&& permuted.Cell(0).rule.Points()[i].weight == cube.Cell(0).rule.Points()[i].weight);
	const auto wrong_surface_domain = Domain(Tetrahedron(), root);
	const iga::CubicCartesianGridSpec wrong_grid{{{1,1,1}}, {{2,2,2}}, {{1,1,1}}};
	const auto wrong_grid_domain = Domain(Cube(1.0, 2.0), wrong_grid);
	Reject([&] { cube.UsableRule(wrong_surface_domain, 0); });
	Reject([&] { cube.ValidateUsableRule(wrong_grid_domain, 0); });
	const iga::CubicCartesianGridSpec tangent_grid{{{0,0,0}}, {{1,1,1}}, {{2,2,2}}};
	const auto tangent_domain = Domain(Cube(.5, 1.0), tangent_grid);
	const iga::CutCellVolumeQuadratureCatalog tangent(tangent_domain, iga::OctreeCutQuadratureOptions{4,10000,10000,100000});
	assert(tangent.Cell(0).classification == iga::CellClassification::Cut);
	assert(tangent.Cell(0).rule.Points().empty() && tangent.Cell(0).usable);
	assert(tangent.Cell(0).diagnostics.certified_reference_volume == 0.0);
	assert(tangent.Cell(4).rule.Points().empty() && tangent.Cell(4).usable); // edge contact
	assert(tangent.Cell(6).rule.Points().empty() && tangent.Cell(6).usable); // face contact
	const auto gauss_boundary_domain = Domain(Cube(iga::kGaussFourPoints[1], 1.0), root);
	const iga::CutCellVolumeQuadratureCatalog gauss_boundary(gauss_boundary_domain, iga::OctreeCutQuadratureOptions{0,100,100,10000});
	assert(gauss_boundary.Cell(0).diagnostics.boundary_samples > 0);
	assert(gauss_boundary.Cell(0).rule.Points().size() < 64);

	Reject([&] { iga::CutCellVolumeQuadratureCatalog bad(cube_domain, iga::OctreeCutQuadratureOptions{21,1,1,1}); });
	Reject([&] { iga::CutCellVolumeQuadratureCatalog bad(cube_domain, iga::OctreeCutQuadratureOptions{1,0,1,1}); });
	Reject([&] { iga::CutCellVolumeQuadratureCatalog bad(cube_domain, iga::OctreeCutQuadratureOptions{1,1,0,1}); });
	Reject([&] { iga::CutCellVolumeQuadratureCatalog bad(cube_domain, iga::OctreeCutQuadratureOptions{1,1,1,0}); });
	const auto cap_domain = Domain(Tetrahedron(), root);
	Reject([&] { iga::CutCellVolumeQuadratureCatalog bad(cap_domain, iga::OctreeCutQuadratureOptions{1,1,100,100}); });
	Reject([&] { iga::CutCellVolumeQuadratureCatalog bad(cap_domain, iga::OctreeCutQuadratureOptions{1,100,1,1000}); });
	Reject([&] { iga::CutCellVolumeQuadratureCatalog bad(cap_domain, iga::OctreeCutQuadratureOptions{0,100,100,1}); });

#ifdef IGA_EXACT_DYADIC_TESTING
	{
		struct Reset { ~Reset() { iga::exact_dyadic::SetTestMagnitudeCap(128); } } reset;
		const auto ambiguous_domain = Domain(Cube(0.0, 1.0), root);
		iga::exact_dyadic::SetTestMagnitudeCap(1);
		const iga::CutCellVolumeQuadratureCatalog ambiguous(ambiguous_domain, iga::OctreeCutQuadratureOptions{0,100,100,10000});
		assert(!ambiguous.Cell(0).usable && ambiguous.Cell(0).diagnostics.predicate_ambiguities > 0);
		Reject([&] { ambiguous.UsableRule(ambiguous_domain, 0); });
	}
#endif

	const double narrow_lower = 1.0e16;
	const double narrow_upper = std::nextafter(narrow_lower, std::numeric_limits<double>::infinity());
	const iga::CubicCartesianGridSpec narrow_grid{{{narrow_lower,narrow_lower,narrow_lower}}, {{narrow_upper,narrow_upper,narrow_upper}}, {{1,1,1}}};
	const auto narrow_domain = Domain(Cube(narrow_lower, narrow_upper), narrow_grid);
	const iga::CutCellVolumeQuadratureCatalog narrow(narrow_domain, iga::OctreeCutQuadratureOptions{4,10000,10000,100000});
	assert(narrow.Cell(0).diagnostics.precision_limited_leaves > 0);

	std::cout << "cut-cell volume quadrature tests passed\n";
}
