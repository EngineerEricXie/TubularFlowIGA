#ifndef IGA_SURFACE_SPATIAL_INDEX_HPP
#define IGA_SURFACE_SPATIAL_INDEX_HPP

#include "ExactDyadicPredicates.hpp"
#include "SurfaceGeometry.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace iga {

enum class PointLocation { Outside, Inside, Boundary, Ambiguous };
enum class BoxContact { None, BoundaryOnly, Interior, Ambiguous };

struct SurfaceSpatialIndexDiagnostics {
	std::size_t node_count = 0;
	std::size_t point_candidate_count = 0;
	std::size_t ray_candidate_count = 0;
	std::size_t box_candidate_count = 0;
	std::size_t narrow_phase_count = 0;
	std::size_t ray_crossing_count = 0;
};

class SurfaceSpatialIndex {
public:
	SurfaceSpatialIndex(const SurfaceSpatialIndex&) = delete;
	SurfaceSpatialIndex& operator=(const SurfaceSpatialIndex&) = delete;
	explicit SurfaceSpatialIndex(ClosedTriangulatedSurface surface) : surface_(std::move(surface))
	{
		const auto& triangles = surface_.Triangles(); ids_.resize(triangles.size()); bounds_.reserve(triangles.size());
		for (std::size_t i = 0; i < triangles.size(); ++i) { ids_[i] = i; bounds_.push_back(triangles[i].bounds); }
		if (triangles.empty()) throw std::invalid_argument("surface spatial index needs triangles");
		if (triangles.size() > std::numeric_limits<std::size_t>::max()/2) throw std::overflow_error("surface spatial index node count overflows");
		nodes_.reserve(2*triangles.size()); BuildNode(0, ids_.size());
	}
	SurfaceSpatialIndex(SurfaceSpatialIndex&& other) noexcept : surface_(std::move(other.surface_)), bounds_(std::move(other.bounds_)), ids_(std::move(other.ids_)), nodes_(std::move(other.nodes_))
	{ point_candidates_.store(other.point_candidates_.load()); ray_candidates_.store(other.ray_candidates_.load()); box_candidates_.store(other.box_candidates_.load()); narrow_.store(other.narrow_.load()); crossings_.store(other.crossings_.load()); }
	SurfaceSpatialIndex& operator=(SurfaceSpatialIndex&& other) noexcept
	{ if (this != &other) { surface_ = std::move(other.surface_); bounds_ = std::move(other.bounds_); ids_ = std::move(other.ids_); nodes_ = std::move(other.nodes_); point_candidates_.store(other.point_candidates_.load()); ray_candidates_.store(other.ray_candidates_.load()); box_candidates_.store(other.box_candidates_.load()); narrow_.store(other.narrow_.load()); crossings_.store(other.crossings_.load()); } return *this; }

	const ClosedTriangulatedSurface& Surface() const { return surface_; }
	SurfaceSpatialIndexDiagnostics Diagnostics() const
	{
		SurfaceSpatialIndexDiagnostics result; result.node_count = nodes_.size();
		result.point_candidate_count = point_candidates_.load(); result.ray_candidate_count = ray_candidates_.load(); result.box_candidate_count = box_candidates_.load();
		result.narrow_phase_count = narrow_.load(); result.ray_crossing_count = crossings_.load(); return result;
	}

	PointLocation LocatePoint(const std::array<double, 3>& point) const
	{
		for (const auto coordinate : point) if (!std::isfinite(coordinate)) throw std::invalid_argument("surface point query is nonfinite");
		try {
			const SurfaceAabb point_box{{point[0], point[1], point[2]}, {point[0], point[1], point[2]}};
			const auto candidates = Candidates(point_box, true);
			for (const auto id : candidates) { Increment(narrow_); if (PointOnTriangle(point, id)) return PointLocation::Boundary; }
			const auto& bounds = surface_.Bounds(); for (std::size_t axis = 0; axis < 3; ++axis) if (exact_dyadic::CoordinateDifferenceSign(point[axis], bounds.minimum[axis]) < 0 || exact_dyadic::CoordinateDifferenceSign(point[axis], bounds.maximum[axis]) > 0) return PointLocation::Outside;
			const SurfaceAabb ray_box{{point[0], point[1], point[2]}, {bounds.maximum[0], point[1], point[2]}};
			const auto candidates_x = Candidates(ray_box, false, true); std::size_t crossings = 0;
			for (const auto id : candidates_x) { Increment(narrow_); if (RayCrossesPositiveX(point, id)) ++crossings; }
			Increment(crossings_, crossings); return crossings%2 == 0 ? PointLocation::Outside : PointLocation::Inside;
		} catch (const std::overflow_error&) { return PointLocation::Ambiguous; }
	}

	BoxContact IntersectBox(const SurfaceAabb& box, std::vector<std::size_t>* canonical_triangle_ids = nullptr) const
	{
		for (std::size_t axis = 0; axis < 3; ++axis) if (!std::isfinite(box.minimum[axis]) || !std::isfinite(box.maximum[axis]) || box.minimum[axis] > box.maximum[axis]) throw std::invalid_argument("surface box query is invalid");
		const auto candidates = Candidates(box, false); BoxContact result = BoxContact::None;
		if (canonical_triangle_ids) canonical_triangle_ids->clear();
		for (const auto id : candidates) {
			Increment(narrow_); BoxContact contact = BoxContact::Ambiguous;
			try { contact = TriangleBoxContact(id, box); }
			catch (const std::overflow_error&) { if (canonical_triangle_ids) canonical_triangle_ids->push_back(id); result = BoxContact::Ambiguous; continue; }
			if (contact == BoxContact::None) continue;
			if (canonical_triangle_ids) canonical_triangle_ids->push_back(id);
			if (contact == BoxContact::Ambiguous) result = BoxContact::Ambiguous;
			else if (result != BoxContact::Ambiguous && contact == BoxContact::Interior) result = BoxContact::Interior;
			else if (result == BoxContact::None) result = BoxContact::BoundaryOnly;
		}
		return result;
	}

private:
	struct Node { SurfaceAabb bounds; std::size_t begin = 0, end = 0, left = 0, right = 0; bool leaf = false; };
	static std::array<double, 3> SubtractD(const std::array<double, 3>& left, const std::array<double, 3>& right)
	{ return {{left[0]-right[0], left[1]-right[1], left[2]-right[2]}}; }
	static std::array<double, 3> CrossD(const std::array<double, 3>& left, const std::array<double, 3>& right)
	{ return {{left[1]*right[2]-left[2]*right[1], left[2]*right[0]-left[0]*right[2], left[0]*right[1]-left[1]*right[0]}}; }
	static std::array<long double, 3> SubtractL(const std::array<long double, 3>& left, const std::array<long double, 3>& right)
	{ return {{left[0]-right[0], left[1]-right[1], left[2]-right[2]}}; }
	static std::array<long double, 3> CrossL(const std::array<long double, 3>& left, const std::array<long double, 3>& right)
	{ return {{left[1]*right[2]-left[2]*right[1], left[2]*right[0]-left[0]*right[2], left[0]*right[1]-left[1]*right[0]}}; }
	static void Increment(std::atomic<std::size_t>& count, std::size_t amount = 1)
	{
		std::size_t old = count.load();
		while (old != std::numeric_limits<std::size_t>::max()) {
			const std::size_t next = amount > std::numeric_limits<std::size_t>::max()-old ? std::numeric_limits<std::size_t>::max() : old+amount;
			if (count.compare_exchange_weak(old, next)) return;
		}
	}
	static bool Overlaps(const SurfaceAabb& left, const SurfaceAabb& right)
	{
		for (std::size_t axis = 0; axis < 3; ++axis) if (left.maximum[axis] < right.minimum[axis] || right.maximum[axis] < left.minimum[axis]) return false;
		return true;
	}
	std::size_t BuildNode(std::size_t begin, std::size_t end)
	{
		Node node; node.begin = begin; node.end = end; node.bounds = bounds_[ids_[begin]];
		for (std::size_t i = begin+1; i < end; ++i) for (std::size_t axis = 0; axis < 3; ++axis) { node.bounds.minimum[axis] = std::min(node.bounds.minimum[axis], bounds_[ids_[i]].minimum[axis]); node.bounds.maximum[axis] = std::max(node.bounds.maximum[axis], bounds_[ids_[i]].maximum[axis]); }
		const auto result = nodes_.size(); nodes_.push_back(node); if (end-begin <= 4) { nodes_[result].leaf = true; return result; }
		std::size_t axis = 0; if (node.bounds.maximum[1]-node.bounds.minimum[1] > node.bounds.maximum[axis]-node.bounds.minimum[axis]) axis = 1; if (node.bounds.maximum[2]-node.bounds.minimum[2] > node.bounds.maximum[axis]-node.bounds.minimum[axis]) axis = 2;
		std::sort(ids_.begin()+static_cast<std::ptrdiff_t>(begin), ids_.begin()+static_cast<std::ptrdiff_t>(end), [&](std::size_t left, std::size_t right) { const double lc = (bounds_[left].minimum[axis]+bounds_[left].maximum[axis])/2.0, rc = (bounds_[right].minimum[axis]+bounds_[right].maximum[axis])/2.0; return lc == rc ? left < right : lc < rc; });
		const auto middle = begin+(end-begin)/2; nodes_[result].left = BuildNode(begin, middle); nodes_[result].right = BuildNode(middle, end); return result;
	}
	std::vector<std::size_t> Candidates(const SurfaceAabb& box, bool point_query, bool ray_query = false) const
	{
		std::vector<std::size_t> result, stack(1, 0);
		while (!stack.empty()) { const auto index = stack.back(); stack.pop_back(); const auto& node = nodes_[index]; if (!Overlaps(node.bounds, box)) continue; if (node.leaf) for (std::size_t i = node.begin; i < node.end; ++i) if (Overlaps(bounds_[ids_[i]], box)) result.push_back(ids_[i]); else {} else { stack.push_back(node.right); stack.push_back(node.left); } }
		std::sort(result.begin(), result.end()); if (point_query) Increment(point_candidates_, result.size()); else if (ray_query) Increment(ray_candidates_, result.size()); else Increment(box_candidates_, result.size()); return result;
	}
	bool PointOnTriangle(const std::array<double, 3>& point, std::size_t id) const
	{
		const auto& triangle = surface_.Triangles()[id]; const auto& vertices = surface_.Vertices(); const auto& a = vertices[triangle.indices[0]], b = vertices[triangle.indices[1]], c = vertices[triangle.indices[2]];
		if (exact_dyadic::Orient3D(a, b, c, point) != 0) return false;
		const auto normal = CrossD(SubtractD(b, a), SubtractD(c, a)); std::size_t drop = 0; if (std::fabs(normal[1]) > std::fabs(normal[drop])) drop = 1; if (std::fabs(normal[2]) > std::fabs(normal[drop])) drop = 2;
		const std::size_t x = (drop+1)%3, y = (drop+2)%3;
		const int ab = exact_dyadic::Orient2D(a, b, point, x, y), bc = exact_dyadic::Orient2D(b, c, point, x, y), ca = exact_dyadic::Orient2D(c, a, point, x, y);
		return (ab >= 0 && bc >= 0 && ca >= 0) || (ab <= 0 && bc <= 0 && ca <= 0);
	}
	bool RayCrossesPositiveX(const std::array<double, 3>& point, std::size_t id) const
	{
		const auto& tri = surface_.Triangles()[id]; const auto& v = surface_.Vertices(); const auto& a = v[tri.indices[0]], b = v[tri.indices[1]], c = v[tri.indices[2]];
		const int area = exact_dyadic::Orient2D(a, b, c, 1, 2); if (area == 0) return false;
		const int ab = exact_dyadic::Orient2D(a, b, point, 1, 2), bc = exact_dyadic::Orient2D(b, c, point, 1, 2), ca = exact_dyadic::Orient2D(c, a, point, 1, 2);
		auto top_left = [](const std::array<double, 3>& first, const std::array<double, 3>& second) {
			const int dy = exact_dyadic::CoordinateDifferenceSign(second[2], first[2]);
			const int dz = exact_dyadic::CoordinateDifferenceSign(second[1], first[1]);
			return dy > 0 || (dy == 0 && dz < 0);
		};
		auto covered = [&](int side, const std::array<double, 3>& first, const std::array<double, 3>& second) {
			const int normalized = area > 0 ? side : -side; if (normalized > 0) return true; if (normalized < 0) return false;
			return area > 0 ? top_left(first, second) : top_left(second, first);
		};
		if (!covered(ab, a, b) || !covered(bc, b, c) || !covered(ca, c, a)) return false;
		const int plane = exact_dyadic::Orient3D(a, b, c, point); return plane != 0 && (plane < 0) != (area < 0);
	}
	BoxContact TriangleBoxContact(std::size_t id, const SurfaceAabb& box) const
	{
		using exact_dyadic::Add; using exact_dyadic::CompareSigned; using exact_dyadic::Cross; using exact_dyadic::Dot; using exact_dyadic::FromDouble; using exact_dyadic::Multiply; using exact_dyadic::Negate; using exact_dyadic::Number; using exact_dyadic::Point; using exact_dyadic::PointFromDouble; using exact_dyadic::Sign; using exact_dyadic::Subtract; using exact_dyadic::SubtractPoint;
		const auto& tri = surface_.Triangles()[id]; const auto& vertices = surface_.Vertices(); const std::array<Point, 3> p{{PointFromDouble(vertices[tri.indices[0]]), PointFromDouble(vertices[tri.indices[1]]), PointFromDouble(vertices[tri.indices[2]])}};
		const Point e0 = SubtractPoint(p[1], p[0]), e1 = SubtractPoint(p[2], p[1]), e2 = SubtractPoint(p[0], p[2]);
		std::array<Point, 13> axes{}; axes[0] = {{FromDouble(1),FromDouble(0),FromDouble(0)}}; axes[1] = {{FromDouble(0),FromDouble(1),FromDouble(0)}}; axes[2] = {{FromDouble(0),FromDouble(0),FromDouble(1)}}; axes[3] = Cross(e0,e1); std::size_t next = 4;
		for (const auto& edge : {e0,e1,e2}) for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) { Point unit{{FromDouble(0),FromDouble(0),FromDouble(0)}}; unit[coordinate]=FromDouble(1); axes[next++] = Cross(edge,unit); }
		const Point lower = PointFromDouble(box.minimum), upper = PointFromDouble(box.maximum);
		const Number two = FromDouble(2); bool strict_all = true, positive_volume = true;
		for (std::size_t coordinate = 0; coordinate < 3; ++coordinate)
			if (CompareSigned(upper[coordinate], lower[coordinate]) <= 0) positive_volume = false;
		for (const auto& axis : axes) {
			if (Sign(axis[0]) == 0 && Sign(axis[1]) == 0 && Sign(axis[2]) == 0) continue;
			Number radius, minimum, maximum;
			for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
				const Number extent = Subtract(upper[coordinate], lower[coordinate]);
				const Number magnitude = Sign(axis[coordinate]) < 0 ? Negate(axis[coordinate]) : axis[coordinate];
				radius = Add(radius, Multiply(magnitude, extent));
			}
			for (std::size_t vertex = 0; vertex < 3; ++vertex) {
				Number value;
				for (std::size_t coordinate = 0; coordinate < 3; ++coordinate)
					value = Add(value, Multiply(axis[coordinate], Subtract(Multiply(two, p[vertex][coordinate]), Add(lower[coordinate], upper[coordinate]))));
				if (vertex == 0 || CompareSigned(value, minimum) < 0) minimum = value;
				if (vertex == 0 || CompareSigned(value, maximum) > 0) maximum = value;
			}
			const Number negative_radius = Negate(radius);
			if (CompareSigned(maximum, negative_radius) < 0 || CompareSigned(minimum, radius) > 0) return BoxContact::None;
			if (CompareSigned(maximum, negative_radius) <= 0 || CompareSigned(minimum, radius) >= 0) strict_all = false;
		}
		return positive_volume && strict_all ? BoxContact::Interior : BoxContact::BoundaryOnly;
	}

	ClosedTriangulatedSurface surface_;
	std::vector<SurfaceAabb> bounds_; std::vector<std::size_t> ids_; std::vector<Node> nodes_;
	mutable std::atomic<std::size_t> point_candidates_{0}, ray_candidates_{0}, box_candidates_{0}, narrow_{0}, crossings_{0};
};

} // namespace iga

#endif
