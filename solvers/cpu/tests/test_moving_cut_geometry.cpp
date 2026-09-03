#include "MovingCutGeometry.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
using namespace iga;

RawSurfaceTriangle Face(std::int64_t a, std::int64_t b, std::int64_t c, std::int64_t label)
{ RawSurfaceTriangle result; result.indices = {{a,b,c}}; result.boundary_id = label; return result; }

RawSurfaceSoup Shell(double x)
{
	RawSurfaceSoup result;
	result.vertices = {{{{x, .45, .45}}, {{x+.72, .45, .45}}, {{x, 1.17, .45}}, {{x, .45, 1.17}}}};
	result.triangles = {Face(0,2,1,5), Face(0,1,3,9), Face(0,3,2,5), Face(1,2,3,9)};
	return result;
}

PrescribedSurfaceMotion Motion()
{
	auto middle = Shell(.85);
	middle.vertices[1][1] += .06; middle.vertices[2][2] += .04; middle.vertices[3][1] -= .03;
	return PrescribedSurfaceMotion({{0.0, Shell(.15)}, {1.0, middle}, {2.0, Shell(.15)}});
}

void Rejected(const std::function<void()>& action)
{ bool rejected = false; try { action(); } catch (const std::exception&) { rejected = true; } assert(rejected); }

bool Has(const MovingCutGeometryDiagnostics& diagnostics, CellClassification old_value, CellClassification new_value)
{
	auto found = diagnostics.transition_counts.find({old_value,new_value});
	return found != diagnostics.transition_counts.end() && found->second != 0;
}

bool Active(CellClassification value)
{ return value == CellClassification::Inside || value == CellClassification::Cut; }

void AssertTransitions(const MovingCutGeometry& geometry, const MovingCutGeometry& previous)
{
	const auto& records = geometry.Diagnostics().transitions;
	const auto& old_cells = previous.Domain().Cells();
	const auto& new_cells = geometry.Domain().Cells();
	assert(old_cells.size() == new_cells.size());
	std::map<std::pair<CellClassification, CellClassification>, std::size_t> counted;
	std::size_t unchanged = 0, active = 0;
	std::size_t record = 0;
	for (std::size_t id = 0; id < new_cells.size(); ++id) {
		const auto old_value = old_cells[id].classification;
		const auto new_value = new_cells[id].classification;
		if (old_value == new_value) {
			++unchanged;
		} else {
			assert(record < records.size());
			const auto& transition = records[record++];
			assert(transition.id == id);
			assert(transition.old_classification == old_value);
			assert(transition.new_classification == new_value);
			++counted[{old_value, new_value}];
		}
		if (Active(new_value)) ++active;
	}
	assert(record == records.size());
	assert(counted == geometry.Diagnostics().transition_counts);
	assert(unchanged == geometry.Diagnostics().unchanged_cells);
	assert(unchanged + records.size() == new_cells.size());
	assert(active == geometry.Diagnostics().active_cells);
	assert(active == geometry.Domain().ActiveCellIds().size());
}

void AssertInverseTransitions(const MovingCutGeometry& forward, const MovingCutGeometry& reverse)
{
	const auto& first = forward.Diagnostics().transitions;
	const auto& second = reverse.Diagnostics().transitions;
	assert(first.size() == second.size());
	for (std::size_t id = 0; id < first.size(); ++id) {
		assert(first[id].id == second[id].id);
		assert(first[id].old_classification == second[id].new_classification);
		assert(first[id].new_classification == second[id].old_classification);
	}
}

RawSurfaceSoup LargeTetra(double origin)
{
	RawSurfaceSoup result;
	const double inset = 16.0, edge = 128.0;
	result.vertices = {{{{origin+inset, origin+inset, origin+inset}},
		{{origin+inset+edge, origin+inset, origin+inset}},
		{{origin+inset, origin+inset+edge, origin+inset}},
		{{origin+inset, origin+inset, origin+inset+edge}}}};
	result.triangles = {Face(0,2,1,11), Face(0,1,3,12), Face(0,3,2,13), Face(1,2,3,14)};
	return result;
}

PrescribedSurfaceMotion LargeMotion(double origin)
{
	auto first = LargeTetra(origin), second = first;
	const std::array<std::array<double, 3>, 4> displacement{{{{4.,8.,12.}}, {{12.,4.,8.}}, {{8.,12.,4.}}, {{4.,12.,8.}}}};
	for (std::size_t vertex = 0; vertex < second.vertices.size(); ++vertex)
		for (std::size_t axis = 0; axis < 3; ++axis) second.vertices[vertex][axis] += displacement[vertex][axis];
	PrescribedSurfaceMotionOptions motion_options;
	motion_options.require_containment_in_fixed_bounds = true;
	motion_options.fixed_bounds_m.minimum = {{origin, origin, origin}};
	motion_options.fixed_bounds_m.maximum = {{origin+1024., origin+1024., origin+1024.}};
	return PrescribedSurfaceMotion({{0.0, first}, {2.0, second}}, motion_options);
}

std::array<double, 3> LocalBarycentric(const ClosedSurfaceTriangle& triangle,
	const ClosedTriangulatedSurface& surface, const std::array<double, 3>& local_point,
	const std::array<double, 3>& cell_lower)
{
	std::array<double, 3> u{}, v{}, w{};
	const auto& a = surface.Vertices()[triangle.indices[0]];
	const auto& b = surface.Vertices()[triangle.indices[1]];
	const auto& c = surface.Vertices()[triangle.indices[2]];
	for (std::size_t axis = 0; axis < 3; ++axis) {
		u[axis] = (b[axis]-cell_lower[axis])-(a[axis]-cell_lower[axis]);
		v[axis] = (c[axis]-cell_lower[axis])-(a[axis]-cell_lower[axis]);
		w[axis] = local_point[axis]-(a[axis]-cell_lower[axis]);
	}
	const auto dot = [](const std::array<double, 3>& left, const std::array<double, 3>& right) {
		return left[0]*right[0]+left[1]*right[1]+left[2]*right[2];
	};
	const double d00 = dot(u,u), d01 = dot(u,v), d11 = dot(v,v), d20 = dot(w,u), d21 = dot(w,v);
	const double determinant = d00*d11-d01*d01;
	assert(determinant > 0.0);
	const double beta = (d11*d20-d01*d21)/determinant;
	const double gamma = (d00*d21-d01*d20)/determinant;
	return {{1.0-beta-gamma, beta, gamma}};
}

void AssertLargeOriginProvenance()
{
	const double origin = std::ldexp(1.0, 53);
	const auto motion = LargeMotion(origin);
	const CubicCartesianGridSpec grid{{{origin,origin,origin}}, {{origin+1024.,origin+1024.,origin+1024.}}, {{1,1,1}}};
	MovingCutGeometryOptions options;
	options.volume.max_depth = 3;
	options.volume.max_nodes = options.volume.max_leaves = options.volume.max_points = 1000000;
	options.volume.max_records = options.volume.max_logical_points = 1000000;
	options.volume.max_retained_bytes = 100000000;
	const auto geometry = MovingCutGeometry::Build(grid, motion.Evaluate(1.0, 0.0, 2.0), options);
	assert(geometry->Surface().Usable() && geometry->Surface().Diagnostics().output_points > 1);
	std::size_t checked = 0;
	for (std::uint64_t id = 0; id < geometry->Domain().Cells().size(); ++id) {
		if (geometry->Domain().Cells()[id].classification != CellClassification::Cut) continue;
		const auto& rule = geometry->Surface().UsableRule(geometry->Domain(), id);
		const auto& provenance = geometry->Surface().UsableProvenance(geometry->Domain(), id);
		for (std::size_t point = 0; point < rule.Points().size(); ++point) {
			const auto& material = provenance[point];
			const auto& triangle = geometry->Evaluation().Surface().Triangles()[material.canonical_triangle];
			std::array<double, 3> local_point{};
			for (std::size_t axis = 0; axis < 3; ++axis)
				local_point[axis] = rule.Points()[point].parametric[axis]
					*(geometry->Domain().Cells()[id].bounds.maximum[axis]-geometry->Domain().Cells()[id].bounds.minimum[axis]);
			const auto expected_barycentric = LocalBarycentric(triangle, geometry->Evaluation().Surface(),
				local_point, geometry->Domain().Cells()[id].bounds.minimum);
			for (std::size_t corner = 0; corner < 3; ++corner)
				assert(std::abs(expected_barycentric[corner]-material.canonical_barycentric[corner]) < 2.0e-12);
			const auto& triangle_provenance = geometry->Evaluation().CanonicalTriangleProvenance()[material.canonical_triangle];
			std::array<double, 3> expected_velocity{{0.,0.,0.}};
			for (std::size_t corner = 0; corner < 3; ++corner) {
				const auto source_corner = triangle_provenance.canonical_corner_to_source_corner[corner];
				const auto source_vertex = triangle_provenance.source_vertex_indices[source_corner];
				for (std::size_t axis = 0; axis < 3; ++axis)
					expected_velocity[axis] += expected_barycentric[corner]
						*geometry->Evaluation().SourceVertexVelocitiesMPerS()[source_vertex][axis];
			}
			const auto actual_velocity = geometry->Evaluation().WallVelocity(material.canonical_triangle, material.canonical_barycentric);
			for (std::size_t axis = 0; axis < 3; ++axis)
				assert(std::abs(actual_velocity[axis]-expected_velocity[axis]) < 2.0e-12);
			++checked;
		}
	}
	assert(checked > 1);
}
}

int main()
{
	const auto motion = Motion();
	const CubicCartesianGridSpec grid{{{0.,0.,0.}}, {{2.,2.,2.}}, {{6,6,6}}};
	MovingCutGeometryOptions options; options.volume.max_depth = 3; options.volume.max_nodes = 1000000; options.volume.max_leaves = 1000000;
	options.volume.max_points = 1000000; options.volume.max_records = 1000000; options.volume.max_retained_bytes = 100000000; options.volume.max_logical_points = 1000000;
	auto initial = MovingCutGeometry::Build(grid, motion.Evaluate(0.0, 0.0, 1.0), options);
	auto forward = MovingCutGeometry::Build(grid, motion.Evaluate(1.0, 0.0, 1.0), options, initial.get());
	auto backward = MovingCutGeometry::Build(grid, motion.Evaluate(2.0, 1.0, 2.0), options, forward.get());
	auto repeat = MovingCutGeometry::Build(grid, motion.Evaluate(1.0, 0.0, 1.0), options, initial.get());
	auto no_previous = MovingCutGeometry::Build(grid, motion.Evaluate(1.0, 0.0, 1.0), options);
	auto same_current_different_previous = MovingCutGeometry::Build(grid, motion.Evaluate(1.0, 0.0, 1.0), options, backward.get());
	MovingCutGeometryOptions changed_surface_options = options;
	++changed_surface_options.surface.max_candidates;
	auto changed_options = MovingCutGeometry::Build(grid, motion.Evaluate(1.0, 0.0, 1.0), changed_surface_options, initial.get());
	assert(forward->IdentitySha256() == repeat->IdentitySha256());
	assert(forward->GeometryIdentitySha256() == repeat->GeometryIdentitySha256());
	assert(forward->PublicationIdentitySha256() == repeat->PublicationIdentitySha256());
	assert(forward->GeometryIdentitySha256() == no_previous->GeometryIdentitySha256());
	assert(forward->PublicationIdentitySha256() != no_previous->PublicationIdentitySha256());
	assert(forward->GeometryIdentitySha256() == same_current_different_previous->GeometryIdentitySha256());
	assert(forward->PublicationIdentitySha256() != same_current_different_previous->PublicationIdentitySha256());
	assert(forward->GeometryIdentitySha256() != changed_options->GeometryIdentitySha256());
	assert(forward->PublicationIdentitySha256() != changed_options->PublicationIdentitySha256());
	assert(forward->Diagnostics().time_s == forward->Evaluation().EvaluatedTimeS());
	assert(initial->Diagnostics().time_s == 0.0 && backward->Diagnostics().time_s == 2.0);
	assert(initial->GeometryIdentitySha256() != backward->GeometryIdentitySha256());
	assert(forward->Diagnostics().transition_counts == repeat->Diagnostics().transition_counts);
	assert(forward->Diagnostics().unchanged_cells == repeat->Diagnostics().unchanged_cells);
	assert(backward->GeometryIdentitySha256() != forward->GeometryIdentitySha256());
	AssertTransitions(*forward, *initial);
	AssertTransitions(*backward, *forward);
	AssertInverseTransitions(*forward, *backward);
	assert(Has(forward->Diagnostics(), CellClassification::Outside, CellClassification::Cut)
		|| Has(forward->Diagnostics(), CellClassification::Cut, CellClassification::Inside));
	assert(Has(backward->Diagnostics(), CellClassification::Cut, CellClassification::Outside)
		|| Has(backward->Diagnostics(), CellClassification::Inside, CellClassification::Cut));
	assert(forward->Diagnostics().active_cells > 0 && forward->Surface().Usable());
	assert(std::abs(forward->Diagnostics().surface_area_residual_m2) < 1.0e-8);
	assert(!forward->Surface().Diagnostics().area_by_boundary_id.empty());
	for (std::uint64_t id = 0; id < forward->Domain().Cells().size(); ++id) {
		if (forward->Domain().Cells()[id].classification != CellClassification::Cut) continue;
		const auto& rule = forward->Surface().UsableRule(forward->Domain(), id);
		const auto& provenance = forward->Surface().UsableProvenance(forward->Domain(), id);
		assert(rule.Points().size() == provenance.size());
		for (std::size_t point = 0; point < rule.Points().size(); ++point) {
			const auto& p = provenance[point]; const auto& triangle = forward->Evaluation().Surface().Triangles()[p.canonical_triangle];
			std::array<double,3> reconstructed{{0.,0.,0.}};
			for (std::size_t c = 0; c < 3; ++c) for (std::size_t axis = 0; axis < 3; ++axis)
				reconstructed[axis] += p.canonical_barycentric[c]*forward->Evaluation().Surface().Vertices()[triangle.indices[c]][axis];
			for (std::size_t axis = 0; axis < 3; ++axis) assert(std::abs(reconstructed[axis]-rule.Points()[point].physical[axis]) < 1.0e-9);
			assert(static_cast<std::uint32_t>(rule.Points()[point].boundary_id) == triangle.boundary_id);
			const auto velocity = forward->Evaluation().WallVelocity(p.canonical_triangle, p.canonical_barycentric);
			for (double value : velocity) assert(std::isfinite(value));
		}
	}
	CubicCartesianGridSpec wrong_grid = grid; wrong_grid.cells[0] = 5;
	Rejected([&] { (void)MovingCutGeometry::Build(wrong_grid, motion.Evaluate(1.0, 0.0, 1.0), options, initial.get()); });
	// Topology mismatch is represented by a different source label sequence.
	auto different_first = Shell(.15); auto different = Shell(.85);
	different_first.triangles[0].boundary_id = 99; different.triangles[0].boundary_id = 99;
	const auto label_motion = PrescribedSurfaceMotion({{0.0, different_first}, {1.0, different}});
	Rejected([&] { (void)MovingCutGeometry::Build(grid, label_motion.Evaluate(1.0, 0.0, 1.0), options, initial.get()); });
	std::unique_ptr<MovingCutGeometry> committed = std::move(initial);
	const std::string committed_identity = committed->IdentitySha256();
	MovingCutGeometryOptions invalid = options; invalid.surface.max_points = 1;
	Rejected([&] { auto temporary = MovingCutGeometry::Build(grid, motion.Evaluate(1.0, 0.0, 1.0), invalid, committed.get()); (void)temporary; });
	assert(committed->IdentitySha256() == committed_identity);
	MovingCutGeometryOptions compact = options;
	compact.volume_storage = CutCellVolumeQuadratureStorageMode::Compact;
	const auto compact_geometry = MovingCutGeometry::Build(grid, motion.Evaluate(1.0, 0.0, 1.0), compact, initial.get());
	assert(compact_geometry->GeometryIdentitySha256() != forward->GeometryIdentitySha256());
	AssertLargeOriginProvenance();
	std::cout << "moving cut geometry tests passed\n";
}
