#ifndef IGA_SURFACE_GEOMETRY_HPP
#define IGA_SURFACE_GEOMETRY_HPP

#include "Sha256.hpp"
#include "ExactDyadicPredicates.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace iga {

struct RawSurfaceTriangle {
	std::array<std::int64_t, 3> indices{{0, 0, 0}};
	// -1 selects SurfaceValidationOptions::default_boundary_id.
	std::int64_t boundary_id = -1;
};

struct RawSurfaceSoup {
	std::vector<std::array<double, 3>> vertices;
	std::vector<RawSurfaceTriangle> triangles;
};

struct SurfaceValidationOptions {
	double length_scale_to_m = 1.0;
	double weld_tolerance_m = 0.0;
	std::int64_t default_boundary_id = 0;
	std::int64_t max_triangles = std::numeric_limits<std::int64_t>::max();
};

struct SurfaceAabb {
	std::array<double, 3> minimum{{0.0, 0.0, 0.0}};
	std::array<double, 3> maximum{{0.0, 0.0, 0.0}};
};

struct ClosedSurfaceTriangle {
	std::array<std::uint32_t, 3> indices{{0, 0, 0}};
	std::uint32_t boundary_id = 0;
	std::array<double, 3> outward_unit_normal{{0.0, 0.0, 0.0}};
	double area_m2 = 0.0;
	SurfaceAabb bounds;
};

struct SurfaceDiagnostics {
	std::size_t input_vertex_count = 0;
	std::size_t input_triangle_count = 0;
	std::size_t canonical_vertex_count = 0;
	std::size_t welded_vertex_count = 0;
	double applied_length_scale_to_m = 1.0;
	double applied_weld_tolerance_m = 0.0;
	std::size_t vertex_count = 0;
	std::size_t triangle_count = 0;
	std::size_t component_count = 0;
	std::size_t unique_edge_count = 0;
	std::size_t vertex_fan_count = 0;
	std::size_t bvh_broad_phase_candidate_count = 0;
	std::size_t bvh_narrow_phase_candidate_count = 0;
	std::size_t accepted_self_intersection_count = 0;
	bool flipped_inward_shell = false;
	SurfaceAabb bounds;
	double area_m2 = 0.0;
	double signed_volume_m3 = 0.0;
	double volume_m3 = 0.0;
	double minimum_edge_length_m = 0.0;
	double minimum_triangle_area_m2 = 0.0;
	std::map<std::uint32_t, std::size_t> boundary_id_histogram;
	std::string canonical_sha256;
};

class ClosedTriangulatedSurface {
public:
	static ClosedTriangulatedSurface Build(const RawSurfaceSoup& soup,
		const SurfaceValidationOptions& options = SurfaceValidationOptions())
	{
		ValidateOptions(options);
		if (soup.vertices.empty()) throw std::invalid_argument("surface has no vertices");
		if (soup.triangles.empty()) throw std::invalid_argument("surface has no triangles");
		if (static_cast<std::uintmax_t>(soup.triangles.size())
			> static_cast<std::uintmax_t>(options.max_triangles))
			throw std::invalid_argument("surface triangle count exceeds configured maximum");
		if (soup.vertices.size() > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()))
			throw std::overflow_error("surface vertex count exceeds canonical index range");

		std::vector<std::array<double, 3>> scaled;
		scaled.reserve(soup.vertices.size());
		for (const auto& vertex : soup.vertices) {
			std::array<double, 3> result{};
			for (std::size_t axis = 0; axis < 3; ++axis) {
				if (!std::isfinite(vertex[axis])) throw std::invalid_argument("surface vertex is nonfinite");
				result[axis] = NormalizeZero(vertex[axis]*options.length_scale_to_m);
				if (!std::isfinite(result[axis])) throw std::invalid_argument("scaled surface vertex is nonfinite");
			}
			scaled.push_back(result);
		}

		auto welded = WeldVertices(scaled, options.weld_tolerance_m);
		std::vector<std::array<double, 3>> vertices = std::move(welded.vertices);
		const auto& remap = welded.remap;

		std::vector<WorkingTriangle> triangles;
		triangles.reserve(soup.triangles.size());
		std::set<std::array<std::uint32_t, 3>> facets;
		for (const auto& raw : soup.triangles) {
			WorkingTriangle triangle;
			for (std::size_t corner = 0; corner < 3; ++corner) {
				if (raw.indices[corner] < 0
					|| static_cast<std::uintmax_t>(raw.indices[corner]) >= scaled.size())
					throw std::invalid_argument("surface triangle index is out of range");
				triangle.indices[corner] = remap[static_cast<std::size_t>(raw.indices[corner])];
			}
			const auto label = raw.boundary_id == -1 ? options.default_boundary_id : raw.boundary_id;
			if (label < 0 || static_cast<std::uintmax_t>(label)
				> static_cast<std::uintmax_t>(std::numeric_limits<std::uint32_t>::max()))
				throw std::invalid_argument("surface boundary id is invalid");
			triangle.boundary_id = static_cast<std::uint32_t>(label);
			if (triangle.indices[0] == triangle.indices[1] || triangle.indices[1] == triangle.indices[2]
				|| triangle.indices[2] == triangle.indices[0])
				throw std::invalid_argument("surface triangle has repeated vertices");
			auto facet = triangle.indices;
			std::sort(facet.begin(), facet.end());
			if (!facets.insert(facet).second)
				throw std::invalid_argument("surface has duplicate facet (winding and label ignored)");
			triangles.push_back(triangle);
		}

		const auto bounds = ComputeBounds(vertices);
		const double characteristic_length = CharacteristicLength(bounds);
		const double edge_tolerance = 128.0*std::numeric_limits<double>::epsilon()*characteristic_length;
		const double area_tolerance = edge_tolerance*characteristic_length;
		for (auto& triangle : triangles) ValidateTriangleGeometry(triangle, vertices, edge_tolerance, area_tolerance);

		const auto edge_uses = BuildEdges(triangles);
		ValidateEdges(edge_uses);
		ValidateVertexFans(vertices.size(), triangles, edge_uses);
		const auto components = CountComponents(triangles, edge_uses);
		if (components != 1) throw std::invalid_argument("surface face adjacency is disconnected");

		CanonicalizeTriangleOrder(triangles);
		const auto self_intersection = ValidateSelfIntersection(vertices, triangles, bounds,
			characteristic_length);
		bool flipped = false;
		const auto volume = SignedVolume(vertices, triangles, bounds);
		const auto length = static_cast<long double>(characteristic_length);
		const auto length_cubed = length*length*length;
		const long double volume_tolerance = 1024.0L*std::numeric_limits<long double>::epsilon()
			*(length_cubed+volume.absolute_terms);
		if (!std::isfinite(volume.signed_sum) || !std::isfinite(volume.absolute_terms)
			|| std::fabs(volume.signed_sum) <= volume_tolerance)
			throw std::invalid_argument("surface enclosed volume is numerically zero or nonfinite");
		if (volume.signed_sum < 0.0L) {
			for (auto& triangle : triangles) std::swap(triangle.indices[1], triangle.indices[2]);
			flipped = true;
		}
		CanonicalizeTriangleOrder(triangles);

		ClosedTriangulatedSurface result;
		result.vertices_ = std::move(vertices);
		result.bounds_ = ComputeBounds(result.vertices_);
		result.triangles_.reserve(triangles.size());
		result.diagnostics_.minimum_edge_length_m = std::numeric_limits<double>::infinity();
		result.diagnostics_.minimum_triangle_area_m2 = std::numeric_limits<double>::infinity();
		for (const auto& triangle : triangles) {
			ClosedSurfaceTriangle canonical;
			canonical.indices = triangle.indices;
			canonical.boundary_id = triangle.boundary_id;
			PopulateGeometry(canonical, result.vertices_);
			result.triangles_.push_back(canonical);
			result.diagnostics_.area_m2 += canonical.area_m2;
			result.diagnostics_.minimum_triangle_area_m2 = std::min(result.diagnostics_.minimum_triangle_area_m2, canonical.area_m2);
			result.diagnostics_.minimum_edge_length_m = std::min(result.diagnostics_.minimum_edge_length_m,
				MinimumEdgeLength(canonical.indices, result.vertices_));
			++result.diagnostics_.boundary_id_histogram[canonical.boundary_id];
		}
		const auto final_volume = SignedVolume(result.vertices_, triangles, result.bounds_);
		const double final_signed_volume_m3 = static_cast<double>(final_volume.signed_sum);
		const double final_volume_m3 = static_cast<double>(std::fabs(final_volume.signed_sum));
		if (!std::isfinite(final_signed_volume_m3) || !std::isfinite(final_volume_m3)
			|| final_signed_volume_m3 <= 0.0 || final_volume_m3 <= 0.0)
			throw std::invalid_argument("surface volume is outside public double geometry range");
		result.diagnostics_.vertex_count = result.vertices_.size();
		result.diagnostics_.input_vertex_count = soup.vertices.size();
		result.diagnostics_.input_triangle_count = soup.triangles.size();
		result.diagnostics_.canonical_vertex_count = result.vertices_.size();
		result.diagnostics_.welded_vertex_count = welded.welded_vertex_count;
		result.diagnostics_.applied_length_scale_to_m = options.length_scale_to_m;
		result.diagnostics_.applied_weld_tolerance_m = options.weld_tolerance_m;
		result.diagnostics_.triangle_count = result.triangles_.size();
		result.diagnostics_.component_count = components;
		result.diagnostics_.unique_edge_count = edge_uses.size();
		result.diagnostics_.vertex_fan_count = result.vertices_.size();
		result.diagnostics_.bvh_broad_phase_candidate_count = self_intersection.broad_candidates;
		result.diagnostics_.bvh_narrow_phase_candidate_count = self_intersection.narrow_candidates;
		result.diagnostics_.accepted_self_intersection_count = 0;
		result.diagnostics_.flipped_inward_shell = flipped;
		result.diagnostics_.bounds = result.bounds_;
		result.diagnostics_.signed_volume_m3 = final_signed_volume_m3;
		result.diagnostics_.volume_m3 = final_volume_m3;
		result.diagnostics_.canonical_sha256 = Hash(result.vertices_, result.triangles_);
		return result;
	}

	const std::vector<std::array<double, 3>>& Vertices() const { return vertices_; }
	const std::vector<ClosedSurfaceTriangle>& Triangles() const { return triangles_; }
	const SurfaceAabb& Bounds() const { return bounds_; }
	const SurfaceDiagnostics& Diagnostics() const { return diagnostics_; }
	const std::string& CanonicalSha256() const { return diagnostics_.canonical_sha256; }

private:
	ClosedTriangulatedSurface() = default;

	struct WorkingTriangle { std::array<std::uint32_t, 3> indices{}; std::uint32_t boundary_id = 0; };
	struct EdgeUse { std::size_t face = 0; std::uint32_t from = 0; std::uint32_t to = 0; };
	struct VolumeSum { long double signed_sum = 0.0L; long double absolute_terms = 0.0L; };
	struct WeldingResult {
		std::vector<std::array<double, 3>> vertices;
		std::vector<std::uint32_t> remap;
		std::size_t welded_vertex_count = 0;
	};
	struct SelfIntersectionCounts { std::size_t broad_candidates = 0; std::size_t narrow_candidates = 0; };
	using EdgeKey = std::pair<std::uint32_t, std::uint32_t>;

	static double NormalizeZero(double value) { return value == 0.0 ? 0.0 : value; }
	static bool LexicographicLess(const std::array<double, 3>& left, const std::array<double, 3>& right)
	{
		return left < right;
	}
	static void ValidateOptions(const SurfaceValidationOptions& options)
	{
		if (!std::isfinite(options.length_scale_to_m) || options.length_scale_to_m <= 0.0)
			throw std::invalid_argument("length_scale_to_m must be finite and positive");
		if (!std::isfinite(options.weld_tolerance_m) || options.weld_tolerance_m < 0.0)
			throw std::invalid_argument("weld_tolerance_m must be finite and nonnegative");
		if (options.default_boundary_id < 0 || options.max_triangles < 0)
			throw std::invalid_argument("surface validation ids and counts must be nonnegative");
		if (static_cast<std::uintmax_t>(options.default_boundary_id)
			> static_cast<std::uintmax_t>(std::numeric_limits<std::uint32_t>::max()))
			throw std::invalid_argument("surface default boundary id exceeds canonical range");
	}
	static WeldingResult ExactCanonicalVertices(const std::vector<std::array<double, 3>>& points)
	{
		WeldingResult result;
		result.vertices = points;
		std::sort(result.vertices.begin(), result.vertices.end(), LexicographicLess);
		result.vertices.erase(std::unique(result.vertices.begin(), result.vertices.end()), result.vertices.end());
		result.remap.resize(points.size());
		for (std::size_t source = 0; source < points.size(); ++source) {
			const auto found = std::lower_bound(result.vertices.begin(), result.vertices.end(), points[source], LexicographicLess);
			result.remap[source] = static_cast<std::uint32_t>(found-result.vertices.begin());
		}
		result.welded_vertex_count = points.size()-result.vertices.size();
		return result;
	}
	static double Distance(const std::array<double, 3>& left, const std::array<double, 3>& right)
	{
		return std::hypot(std::hypot(left[0]-right[0], left[1]-right[1]), left[2]-right[2]);
	}
	static WeldingResult WeldVertices(const std::vector<std::array<double, 3>>& points, double tolerance)
	{
		if (tolerance == 0.0) return ExactCanonicalVertices(points);
		SurfaceAabb bounds;
		bounds.minimum = points.front();
		for (const auto& point : points) for (std::size_t axis = 0; axis < 3; ++axis)
			bounds.minimum[axis] = std::min(bounds.minimum[axis], point[axis]);
		struct DisjointSet {
			std::vector<std::size_t> parent;
			explicit DisjointSet(std::size_t count) : parent(count) { for (std::size_t i = 0; i < count; ++i) parent[i] = i; }
			std::size_t Find(std::size_t value) { while (parent[value] != value) { parent[value] = parent[parent[value]]; value = parent[value]; } return value; }
			void Join(std::size_t left, std::size_t right) { left = Find(left); right = Find(right); if (left != right) parent[std::max(left, right)] = std::min(left, right); }
		};
		using Cell = std::tuple<std::int64_t, std::int64_t, std::int64_t>;
		auto cell_for = [tolerance, &bounds](const std::array<double, 3>& point) {
			std::array<std::int64_t, 3> result{};
			for (std::size_t axis = 0; axis < 3; ++axis) {
				const long double relative = static_cast<long double>(point[axis])-static_cast<long double>(bounds.minimum[axis]);
				const long double quotient = std::floor(relative/static_cast<long double>(tolerance));
				if (!std::isfinite(relative) || !std::isfinite(quotient) || quotient < static_cast<long double>(std::numeric_limits<std::int64_t>::min())+1.0L
					|| quotient > static_cast<long double>(std::numeric_limits<std::int64_t>::max())-1.0L)
					throw std::overflow_error("surface weld spatial-hash cell coordinate overflows");
				result[axis] = static_cast<std::int64_t>(quotient);
			}
			return Cell(result[0], result[1], result[2]);
		};
		std::map<Cell, std::vector<std::size_t>> cells;
		std::vector<Cell> point_cells;
		point_cells.reserve(points.size());
		for (std::size_t index = 0; index < points.size(); ++index) {
			const auto cell = cell_for(points[index]);
			point_cells.push_back(cell);
			cells[cell].push_back(index);
		}
		DisjointSet sets(points.size());
		for (std::size_t index = 0; index < points.size(); ++index) {
			const auto& cell = point_cells[index];
			for (int dx = -1; dx <= 1; ++dx) for (int dy = -1; dy <= 1; ++dy) for (int dz = -1; dz <= 1; ++dz) {
				const auto x = std::get<0>(cell), y = std::get<1>(cell), z = std::get<2>(cell);
				if ((dx < 0 && x == std::numeric_limits<std::int64_t>::min()) || (dx > 0 && x == std::numeric_limits<std::int64_t>::max())
					|| (dy < 0 && y == std::numeric_limits<std::int64_t>::min()) || (dy > 0 && y == std::numeric_limits<std::int64_t>::max())
					|| (dz < 0 && z == std::numeric_limits<std::int64_t>::min()) || (dz > 0 && z == std::numeric_limits<std::int64_t>::max()))
					throw std::overflow_error("surface weld spatial-hash neighbor coordinate overflows");
				const Cell neighbor(x+dx, y+dy, z+dz);
				const auto found = cells.find(neighbor);
				if (found == cells.end()) continue;
				for (const auto other : found->second) if (other < index && Distance(points[index], points[other]) <= tolerance)
					sets.Join(index, other);
			}
		}
		std::map<std::size_t, std::vector<std::size_t>> clusters;
		for (std::size_t index = 0; index < points.size(); ++index) clusters[sets.Find(index)].push_back(index);
		std::vector<std::array<double, 3>> representatives(points.size());
		for (const auto& cluster : clusters) {
			const auto& members = cluster.second;
			for (std::size_t i = 0; i < members.size(); ++i) for (std::size_t j = i+1; j < members.size(); ++j)
				if (Distance(points[members[i]], points[members[j]]) > tolerance)
					throw std::invalid_argument("transitive surface weld cluster exceeds tolerance diameter");
			auto representative = points[members.front()];
			for (const auto member : members) if (LexicographicLess(points[member], representative)) representative = points[member];
			for (const auto member : members) representatives[member] = representative;
		}
		return ExactCanonicalVertices(representatives);
	}
	static SurfaceAabb ComputeBounds(const std::vector<std::array<double, 3>>& vertices)
	{
		SurfaceAabb result;
		result.minimum = vertices.front(); result.maximum = vertices.front();
		for (const auto& vertex : vertices) for (std::size_t axis = 0; axis < 3; ++axis) {
			result.minimum[axis] = std::min(result.minimum[axis], vertex[axis]);
			result.maximum[axis] = std::max(result.maximum[axis], vertex[axis]);
		}
		return result;
	}
	static double CharacteristicLength(const SurfaceAabb& bounds)
	{
		const double dx = bounds.maximum[0]-bounds.minimum[0];
		const double dy = bounds.maximum[1]-bounds.minimum[1];
		const double dz = bounds.maximum[2]-bounds.minimum[2];
		const double result = std::hypot(std::hypot(dx, dy), dz);
		if (!std::isfinite(result) || result <= 0.0)
			throw std::invalid_argument("surface bounding-box diagonal is numerically invalid");
		return result;
	}
	static std::array<double, 3> Subtract(const std::array<double, 3>& a, const std::array<double, 3>& b)
	{
		return {{a[0]-b[0], a[1]-b[1], a[2]-b[2]}};
	}
	static std::array<double, 3> Cross(const std::array<double, 3>& a, const std::array<double, 3>& b)
	{
		return {{a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]}};
	}
	static double Dot(const std::array<double, 3>& a, const std::array<double, 3>& b)
	{
		return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
	}
	static double Norm(const std::array<double, 3>& vector)
	{
		return std::hypot(std::hypot(vector[0], vector[1]), vector[2]);
	}
	static void ValidateTriangleGeometry(const WorkingTriangle& triangle,
		const std::vector<std::array<double, 3>>& vertices, double edge_tolerance, double area_tolerance)
	{
		const auto& a = vertices[triangle.indices[0]]; const auto& b = vertices[triangle.indices[1]];
		const auto& c = vertices[triangle.indices[2]];
		if (Norm(Subtract(b, a)) <= edge_tolerance || Norm(Subtract(c, b)) <= edge_tolerance
			|| Norm(Subtract(a, c)) <= edge_tolerance)
			throw std::invalid_argument("surface triangle edge is numerically zero");
		if (0.5*Norm(Cross(Subtract(b, a), Subtract(c, a))) <= area_tolerance)
			throw std::invalid_argument("surface triangle area is numerically zero");
	}
	static std::map<EdgeKey, std::vector<EdgeUse>> BuildEdges(const std::vector<WorkingTriangle>& triangles)
	{
		std::map<EdgeKey, std::vector<EdgeUse>> result;
		for (std::size_t face = 0; face < triangles.size(); ++face) for (std::size_t edge = 0; edge < 3; ++edge) {
			const auto from = triangles[face].indices[edge];
			const auto to = triangles[face].indices[(edge+1)%3];
			result[std::minmax(from, to)].push_back({face, from, to});
		}
		return result;
	}
	static void ValidateEdges(const std::map<EdgeKey, std::vector<EdgeUse>>& edges)
	{
		for (const auto& entry : edges) {
			const auto& uses = entry.second;
			if (uses.size() != 2) throw std::invalid_argument("surface edge is not used by exactly two faces");
			if (uses[0].from != uses[1].to || uses[0].to != uses[1].from)
				throw std::invalid_argument("surface adjacent faces do not have opposite edge directions");
		}
	}
	static std::size_t CountComponents(const std::vector<WorkingTriangle>& triangles,
		const std::map<EdgeKey, std::vector<EdgeUse>>& edges)
	{
		std::vector<std::vector<std::size_t>> adjacency(triangles.size());
		for (const auto& edge : edges) { const auto& use = edge.second; adjacency[use[0].face].push_back(use[1].face); adjacency[use[1].face].push_back(use[0].face); }
		std::vector<bool> seen(triangles.size(), false); std::size_t components = 0;
		for (std::size_t start = 0; start < triangles.size(); ++start) if (!seen[start]) {
			++components; std::queue<std::size_t> pending; pending.push(start); seen[start] = true;
			while (!pending.empty()) { const auto face = pending.front(); pending.pop(); for (const auto next : adjacency[face]) if (!seen[next]) { seen[next] = true; pending.push(next); } }
		}
		return components;
	}
	static void ValidateVertexFans(std::size_t vertex_count, const std::vector<WorkingTriangle>& triangles,
		const std::map<EdgeKey, std::vector<EdgeUse>>& edges)
	{
		std::vector<std::vector<std::size_t>> incident(vertex_count);
		for (std::size_t face = 0; face < triangles.size(); ++face) for (const auto vertex : triangles[face].indices) incident[vertex].push_back(face);
		std::vector<std::vector<std::size_t>> neighbors(triangles.size());
		for (const auto& entry : edges) { const auto& use = entry.second; neighbors[use[0].face].push_back(use[1].face); neighbors[use[1].face].push_back(use[0].face); }
		for (std::size_t vertex = 0; vertex < vertex_count; ++vertex) {
			if (incident[vertex].empty()) throw std::invalid_argument("surface has unused vertex");
			std::set<std::size_t> local(incident[vertex].begin(), incident[vertex].end());
			std::set<std::size_t> reached; std::queue<std::size_t> pending; pending.push(incident[vertex].front()); reached.insert(incident[vertex].front());
			while (!pending.empty()) { const auto face = pending.front(); pending.pop(); std::size_t degree = 0; for (const auto next : neighbors[face]) if (local.count(next)) { ++degree; if (reached.insert(next).second) pending.push(next); } if (degree != 2) throw std::invalid_argument("surface vertex fan is not cyclic"); }
			if (reached.size() != local.size()) throw std::invalid_argument("surface vertex fan is pinched (bow-tie)");
		}
	}
	static VolumeSum SignedVolume(const std::vector<std::array<double, 3>>& vertices,
		const std::vector<WorkingTriangle>& triangles, const SurfaceAabb& bounds)
	{
		std::array<long double, 3> center{{(static_cast<long double>(bounds.minimum[0])+bounds.maximum[0])/2.0L,
			(static_cast<long double>(bounds.minimum[1])+bounds.maximum[1])/2.0L,
			(static_cast<long double>(bounds.minimum[2])+bounds.maximum[2])/2.0L}};
		long double sum = 0.0L, compensation = 0.0L;
		long double absolute_sum = 0.0L, absolute_compensation = 0.0L;
		for (const auto& triangle : triangles) {
			auto point = [&](std::uint32_t index) { return std::array<long double, 3>{{static_cast<long double>(vertices[index][0])-center[0], static_cast<long double>(vertices[index][1])-center[1], static_cast<long double>(vertices[index][2])-center[2]}}; };
			const auto a = point(triangle.indices[0]), b = point(triangle.indices[1]), c = point(triangle.indices[2]);
			const long double term = (a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0]))/6.0L;
			const long double adjusted = term-compensation; const long double next = sum+adjusted; compensation = (next-sum)-adjusted; sum = next;
			const long double magnitude = std::fabs(term);
			const long double absolute_adjusted = magnitude-absolute_compensation;
			const long double absolute_next = absolute_sum+absolute_adjusted;
			absolute_compensation = (absolute_next-absolute_sum)-absolute_adjusted;
			absolute_sum = absolute_next;
		}
		return {sum, absolute_sum};
	}
	static void RotateToSmallestIndex(std::array<std::uint32_t, 3>& indices)
	{
		const auto smallest = std::min_element(indices.begin(), indices.end());
		std::rotate(indices.begin(), smallest, indices.end());
	}
	static void CanonicalizeTriangleOrder(std::vector<WorkingTriangle>& triangles)
	{
		for (auto& triangle : triangles) RotateToSmallestIndex(triangle.indices);
		std::sort(triangles.begin(), triangles.end(), [](const WorkingTriangle& left,
			const WorkingTriangle& right) {
			if (left.indices != right.indices) return left.indices < right.indices;
			return left.boundary_id < right.boundary_id;
		});
	}
	static SurfaceAabb TriangleBounds(const WorkingTriangle& triangle,
		const std::vector<std::array<double, 3>>& vertices)
	{
		SurfaceAabb result;
		result.minimum = vertices[triangle.indices[0]]; result.maximum = result.minimum;
		for (std::size_t corner = 1; corner < 3; ++corner) for (std::size_t axis = 0; axis < 3; ++axis) {
			const auto coordinate = vertices[triangle.indices[corner]][axis];
			result.minimum[axis] = std::min(result.minimum[axis], coordinate);
			result.maximum[axis] = std::max(result.maximum[axis], coordinate);
		}
		return result;
	}
	static bool BoundsOverlap(const SurfaceAabb& left, const SurfaceAabb& right, double tolerance)
	{
		for (std::size_t axis = 0; axis < 3; ++axis)
			if (left.maximum[axis]+tolerance < right.minimum[axis]
				|| right.maximum[axis]+tolerance < left.minimum[axis]) return false;
		return true;
	}
	struct BvhNode { SurfaceAabb bounds; std::size_t begin = 0, end = 0, left = 0, right = 0; bool leaf = false; };
	static std::size_t BuildBvh(std::vector<BvhNode>& nodes, std::vector<std::size_t>& ids,
		const std::vector<SurfaceAabb>& triangle_bounds, std::size_t begin, std::size_t end)
	{
		BvhNode node; node.begin = begin; node.end = end; node.bounds = triangle_bounds[ids[begin]];
		for (std::size_t index = begin+1; index < end; ++index) for (std::size_t axis = 0; axis < 3; ++axis) {
			node.bounds.minimum[axis] = std::min(node.bounds.minimum[axis], triangle_bounds[ids[index]].minimum[axis]);
			node.bounds.maximum[axis] = std::max(node.bounds.maximum[axis], triangle_bounds[ids[index]].maximum[axis]);
		}
		const auto node_index = nodes.size(); nodes.push_back(node);
		if (end-begin <= 4) { nodes[node_index].leaf = true; return node_index; }
		std::size_t axis = 0;
		if (node.bounds.maximum[1]-node.bounds.minimum[1] > node.bounds.maximum[axis]-node.bounds.minimum[axis]) axis = 1;
		if (node.bounds.maximum[2]-node.bounds.minimum[2] > node.bounds.maximum[axis]-node.bounds.minimum[axis]) axis = 2;
		std::sort(ids.begin()+static_cast<std::ptrdiff_t>(begin), ids.begin()+static_cast<std::ptrdiff_t>(end),
			[&triangle_bounds, axis](std::size_t left, std::size_t right) {
				const double left_center = (triangle_bounds[left].minimum[axis]+triangle_bounds[left].maximum[axis])/2.0;
				const double right_center = (triangle_bounds[right].minimum[axis]+triangle_bounds[right].maximum[axis])/2.0;
				return left_center == right_center ? left < right : left_center < right_center;
			});
		const auto middle = begin+(end-begin)/2;
		nodes[node_index].left = BuildBvh(nodes, ids, triangle_bounds, begin, middle);
		nodes[node_index].right = BuildBvh(nodes, ids, triangle_bounds, middle, end);
		return node_index;
	}
	static std::array<long double, 3> ScaledPoint(const std::array<double, 3>& point,
		const SurfaceAabb& bounds, double length)
	{
		return {{(static_cast<long double>(point[0])-(static_cast<long double>(bounds.minimum[0])+bounds.maximum[0])/2.0L)/length,
			(static_cast<long double>(point[1])-(static_cast<long double>(bounds.minimum[1])+bounds.maximum[1])/2.0L)/length,
			(static_cast<long double>(point[2])-(static_cast<long double>(bounds.minimum[2])+bounds.maximum[2])/2.0L)/length}};
	}
	static std::array<long double, 3> SubtractL(const std::array<long double, 3>& a, const std::array<long double, 3>& b)
	{ return {{a[0]-b[0], a[1]-b[1], a[2]-b[2]}}; }
	static std::array<long double, 3> CrossL(const std::array<long double, 3>& a, const std::array<long double, 3>& b)
	{ return {{a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]}}; }
	static long double DotL(const std::array<long double, 3>& a, const std::array<long double, 3>& b)
	{ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
	static long double NormL(const std::array<long double, 3>& vector)
	{ return std::hypot(std::hypot(vector[0], vector[1]), vector[2]); }
	// Coordinates enter as finite binary64 values.  Ambiguous long-double signs
	// are recomputed exactly as signed dyadic integers, not approximated again.
	template <class ExactCalculation>
	static int FilteredSign(long double value, long double operand_magnitude, ExactCalculation exact_calculation)
	{
		const long double bound = 128.0L*std::numeric_limits<long double>::epsilon()*operand_magnitude;
		if (value != 0.0L && std::fabs(value) > bound) return value < 0.0L ? -1 : 1;
		return exact_calculation();
	}
	static int Orient2D(const std::array<long double, 3>& a, const std::array<long double, 3>& b,
		const std::array<long double, 3>& c, std::uint32_t ai, std::uint32_t bi, std::uint32_t ci,
		std::size_t drop, const std::vector<std::array<double, 3>>& vertices)
	{
		const std::size_t x = (drop+1)%3, y = (drop+2)%3;
		const long double first = (b[x]-a[x])*(c[y]-a[y]), second = (b[y]-a[y])*(c[x]-a[x]);
		return FilteredSign(first-second, std::fabs(first)+std::fabs(second), [&] { return exact_dyadic::Orient2D(vertices[ai], vertices[bi], vertices[ci], x, y); });
	}
	static int Orient3D(const std::array<long double, 3>& a, const std::array<long double, 3>& b,
		const std::array<long double, 3>& c, const std::array<long double, 3>& p,
		std::uint32_t ai, std::uint32_t bi, std::uint32_t ci, std::uint32_t pi,
		const std::vector<std::array<double, 3>>& vertices)
	{
		const auto u = SubtractL(b, a), v = SubtractL(c, a), w = SubtractL(p, a);
		const long double q0 = v[1]*w[2], q1 = v[2]*w[1], q2 = v[0]*w[2], q3 = v[2]*w[0], q4 = v[0]*w[1], q5 = v[1]*w[0];
		const long double t0 = u[0]*(q0-q1), t1 = u[1]*(q2-q3), t2 = u[2]*(q4-q5);
		const long double magnitude = std::fabs(u[0])*(std::fabs(q0)+std::fabs(q1))+
			std::fabs(u[1])*(std::fabs(q2)+std::fabs(q3))+std::fabs(u[2])*(std::fabs(q4)+std::fabs(q5));
		return FilteredSign(t0-t1+t2, magnitude, [&] { return exact_dyadic::Orient3D(vertices[ai], vertices[bi], vertices[ci], vertices[pi]); });
	}
	static int DirectionCrossComponent(const std::array<long double, 3>& a, const std::array<long double, 3>& b,
		const std::array<long double, 3>& c, const std::array<long double, 3>& d,
		std::uint32_t ai, std::uint32_t bi, std::uint32_t ci, std::uint32_t di, std::size_t axis,
		const std::vector<std::array<double, 3>>& vertices)
	{
		const auto u = SubtractL(b, a), v = SubtractL(d, c);
		const std::size_t x = (axis+1)%3, y = (axis+2)%3;
		const long double first = u[x]*v[y], second = u[y]*v[x];
		return FilteredSign(first-second, std::fabs(first)+std::fabs(second), [&] { return exact_dyadic::DirectionCrossComponent(vertices[ai], vertices[bi], vertices[ci], vertices[di], x, y); });
	}
	struct SegmentTriangleHit { bool hit = false; bool at_start = false; bool at_end = false; std::uint32_t triangle_vertex = std::numeric_limits<std::uint32_t>::max(); };
	static SegmentTriangleHit SegmentTriangleContact(const std::array<long double, 3>& start,
		const std::array<long double, 3>& end, std::uint32_t start_id, std::uint32_t end_id,
		const std::array<long double, 3>& a, const std::array<long double, 3>& b, const std::array<long double, 3>& c,
		const std::array<std::uint32_t, 3>& ids, const std::vector<std::array<double, 3>>& vertices)
	{
		const auto direction = SubtractL(end, start), e1 = SubtractL(b, a), e2 = SubtractL(c, a), tvec = SubtractL(start, a);
		const auto p = CrossL(direction, e2), q = CrossL(tvec, e1);
		const long double determinant = DotL(e1, p), u = DotL(tvec, p), v = DotL(direction, q), t = DotL(e2, q);
		auto exact_values = [&] { return exact_dyadic::SegmentTriangleValues(vertices[start_id], vertices[end_id], vertices[ids[0]], vertices[ids[1]], vertices[ids[2]]); };
		auto cross_magnitude = [](const std::array<long double, 3>& left, const std::array<long double, 3>& right, const std::array<long double, 3>& multiplier) {
			return std::fabs(multiplier[0])*(std::fabs(left[1]*right[2])+std::fabs(left[2]*right[1]))+
				std::fabs(multiplier[1])*(std::fabs(left[2]*right[0])+std::fabs(left[0]*right[2]))+
				std::fabs(multiplier[2])*(std::fabs(left[0]*right[1])+std::fabs(left[1]*right[0]));
		};
		const long double determinant_magnitude = cross_magnitude(direction, e2, e1), u_magnitude = cross_magnitude(direction, e2, tvec);
		const long double v_magnitude = cross_magnitude(tvec, e1, direction), t_magnitude = cross_magnitude(tvec, e1, e2);
		auto sign = [&](long double value, long double magnitude, std::size_t item) { return FilteredSign(value, magnitude, [&] { return exact_dyadic::Sign(exact_values()[item]); }); };
		int determinant_sign = sign(determinant, determinant_magnitude, 0); if (determinant_sign == 0) return {};
		int u_sign = sign(u, u_magnitude, 1), v_sign = sign(v, v_magnitude, 2), t_sign = sign(t, t_magnitude, 3);
		int duv_sign = FilteredSign(determinant-u-v, determinant_magnitude+u_magnitude+v_magnitude, [&] { const auto values = exact_values(); return exact_dyadic::Sign(exact_dyadic::Subtract(exact_dyadic::Subtract(values[0], values[1]), values[2])); });
		int dt_sign = FilteredSign(determinant-t, determinant_magnitude+t_magnitude, [&] { const auto values = exact_values(); return exact_dyadic::Sign(exact_dyadic::Subtract(values[0], values[3])); });
		if (determinant_sign < 0) { u_sign = -u_sign; v_sign = -v_sign; t_sign = -t_sign; duv_sign = -duv_sign; dt_sign = -dt_sign; }
		if (u_sign < 0 || v_sign < 0 || t_sign < 0 || duv_sign < 0 || dt_sign < 0) return {};
		SegmentTriangleHit result; result.hit = true; result.at_start = t_sign == 0; result.at_end = dt_sign == 0;
		if (u_sign == 0 && v_sign == 0) result.triangle_vertex = ids[0];
		else if (duv_sign == 0 && v_sign == 0) result.triangle_vertex = ids[1];
		else if (u_sign == 0 && duv_sign == 0) result.triangle_vertex = ids[2];
		return result;
	}
	static bool EdgeContains(std::uint32_t first, std::uint32_t second, std::uint32_t vertex)
	{ return first == vertex || second == vertex; }
	static bool IsAllowedPointContact(const SegmentTriangleHit& hit, std::uint32_t start, std::uint32_t end,
		const std::array<std::uint32_t, 3>& shared, std::size_t shared_count)
	{
		for (std::size_t i = 0; i < shared_count; ++i)
			if (hit.triangle_vertex == shared[i]
				&& ((hit.at_start && start == shared[i]) || (hit.at_end && end == shared[i]))) return true;
		return false;
	}
	enum class Segment2DContact { none, point, overlap };
	static bool IsAllowedEdgeContact(Segment2DContact contact, std::uint32_t ai, std::uint32_t aj,
		std::uint32_t bi, std::uint32_t bj, const std::array<std::uint32_t, 3>& shared, std::size_t shared_count)
	{
		if (contact == Segment2DContact::point) for (std::size_t i = 0; i < shared_count; ++i)
			if (EdgeContains(ai, aj, shared[i]) && EdgeContains(bi, bj, shared[i])) return true;
		return contact == Segment2DContact::overlap && shared_count == 2
			&& EdgeContains(ai, aj, shared[0]) && EdgeContains(ai, aj, shared[1])
			&& EdgeContains(bi, bj, shared[0]) && EdgeContains(bi, bj, shared[1]);
	}
	static Segment2DContact SegmentContact2D(const std::array<long double, 3>& a, const std::array<long double, 3>& b,
		const std::array<long double, 3>& c, const std::array<long double, 3>& d,
		std::uint32_t ai, std::uint32_t bi, std::uint32_t ci, std::uint32_t di, std::size_t drop,
		const std::vector<std::array<double, 3>>& vertices)
	{
		const int ab_c = Orient2D(a, b, c, ai, bi, ci, drop, vertices), ab_d = Orient2D(a, b, d, ai, bi, di, drop, vertices);
		const int cd_a = Orient2D(c, d, a, ci, di, ai, drop, vertices), cd_b = Orient2D(c, d, b, ci, di, bi, drop, vertices);
		if (ab_c == 0 && ab_d == 0 && cd_a == 0 && cd_b == 0) {
			const auto contact = exact_dyadic::ClassifyCollinearIntervals(vertices[ai], vertices[bi], vertices[ci], vertices[di]);
			return contact == exact_dyadic::CollinearIntervalContact::none ? Segment2DContact::none
				: (contact == exact_dyadic::CollinearIntervalContact::point ? Segment2DContact::point : Segment2DContact::overlap);
		}
		if ((ab_c == 0 || ab_d == 0 || ab_c != ab_d) && (cd_a == 0 || cd_b == 0 || cd_a != cd_b)) return Segment2DContact::point;
		return Segment2DContact::none;
	}
	static bool PointInTriangle2D(std::uint32_t point, const std::array<long double, 3>& point_l,
		const std::array<std::uint32_t, 3>& triangle, const std::array<std::array<long double, 3>, 3>& points, std::size_t drop,
		const std::vector<std::array<double, 3>>& vertices)
	{
		const int a = Orient2D(points[0], points[1], point_l, triangle[0], triangle[1], point, drop, vertices);
		const int b = Orient2D(points[1], points[2], point_l, triangle[1], triangle[2], point, drop, vertices);
		const int c = Orient2D(points[2], points[0], point_l, triangle[2], triangle[0], point, drop, vertices);
		return (a >= 0 && b >= 0 && c >= 0) || (a <= 0 && b <= 0 && c <= 0);
	}
	static bool PairIntersectsForbidden(const WorkingTriangle& left, const WorkingTriangle& right,
		const std::vector<std::array<double, 3>>& vertices, const SurfaceAabb& bounds, double length)
	{
		std::array<std::uint32_t, 3> shared{}; std::size_t shared_count = 0;
		for (const auto a : left.indices) for (const auto b : right.indices) if (a == b) shared[shared_count++] = a;
		std::array<std::array<long double, 3>, 3> a{}, b{};
		for (std::size_t i = 0; i < 3; ++i) { a[i] = ScaledPoint(vertices[left.indices[i]], bounds, length); b[i] = ScaledPoint(vertices[right.indices[i]], bounds, length); }
		const auto na = CrossL(SubtractL(a[1], a[0]), SubtractL(a[2], a[0]));
		const auto nb = CrossL(SubtractL(b[1], b[0]), SubtractL(b[2], b[0]));
		const bool coplanar = Orient3D(a[0], a[1], a[2], b[0], left.indices[0], left.indices[1], left.indices[2], right.indices[0], vertices) == 0
			&& Orient3D(a[0], a[1], a[2], b[1], left.indices[0], left.indices[1], left.indices[2], right.indices[1], vertices) == 0
			&& Orient3D(a[0], a[1], a[2], b[2], left.indices[0], left.indices[1], left.indices[2], right.indices[2], vertices) == 0;
		if (coplanar) {
			std::size_t drop = 0; if (std::fabs(na[1]) > std::fabs(na[drop])) drop = 1; if (std::fabs(na[2]) > std::fabs(na[drop])) drop = 2;
			for (std::size_t i = 0; i < 3; ++i) for (std::size_t j = 0; j < 3; ++j) {
				const auto ai = left.indices[i], aj = left.indices[(i+1)%3], bi = right.indices[j], bj = right.indices[(j+1)%3];
				const auto contact = SegmentContact2D(a[i], a[(i+1)%3], b[j], b[(j+1)%3], ai, aj, bi, bj, drop, vertices);
				if (contact == Segment2DContact::none) continue;
				bool allowed = false;
				if (contact == Segment2DContact::point) for (std::size_t shared_index = 0; shared_index < shared_count; ++shared_index)
					if (EdgeContains(ai, aj, shared[shared_index]) && EdgeContains(bi, bj, shared[shared_index])) allowed = true;
				if (contact == Segment2DContact::overlap && shared_count == 2
					&& EdgeContains(ai, aj, shared[0]) && EdgeContains(ai, aj, shared[1])
					&& EdgeContains(bi, bj, shared[0]) && EdgeContains(bi, bj, shared[1])) allowed = true;
				if (!allowed) return true;
			}
			for (std::size_t i = 0; i < 3; ++i) if (PointInTriangle2D(left.indices[i], a[i], right.indices, b, drop, vertices)) {
				bool allowed = false; for (std::size_t k = 0; k < shared_count; ++k) if (left.indices[i] == shared[k]) allowed = true;
				if (!allowed) return true;
			}
			for (std::size_t i = 0; i < 3; ++i) if (PointInTriangle2D(right.indices[i], b[i], left.indices, a, drop, vertices)) {
				bool allowed = false; for (std::size_t k = 0; k < shared_count; ++k) if (right.indices[i] == shared[k]) allowed = true;
				if (!allowed) return true;
			}
			return false;
		}
		bool forbidden = false;
		for (std::size_t edge = 0; edge < 3; ++edge) {
			const auto first = SegmentTriangleContact(a[edge], a[(edge+1)%3], left.indices[edge], left.indices[(edge+1)%3],
				b[0], b[1], b[2], right.indices, vertices);
			const auto second = SegmentTriangleContact(b[edge], b[(edge+1)%3], right.indices[edge], right.indices[(edge+1)%3],
				a[0], a[1], a[2], left.indices, vertices);
			if (first.hit && !IsAllowedPointContact(first, left.indices[edge], left.indices[(edge+1)%3], shared, shared_count)) forbidden = true;
			if (second.hit && !IsAllowedPointContact(second, right.indices[edge], right.indices[(edge+1)%3], shared, shared_count)) forbidden = true;
		}
		auto edge_in_other_plane = [&](const auto& edge_points, const auto& edge_ids, const auto& triangle_points,
			const auto& triangle_ids, const auto& normal) {
			if (Orient3D(triangle_points[0], triangle_points[1], triangle_points[2], edge_points[0], triangle_ids[0], triangle_ids[1], triangle_ids[2], edge_ids[0], vertices) != 0
				|| Orient3D(triangle_points[0], triangle_points[1], triangle_points[2], edge_points[1], triangle_ids[0], triangle_ids[1], triangle_ids[2], edge_ids[1], vertices) != 0) return;
			std::size_t drop = 0; if (std::fabs(normal[1]) > std::fabs(normal[drop])) drop = 1; if (std::fabs(normal[2]) > std::fabs(normal[drop])) drop = 2;
			for (std::size_t endpoint = 0; endpoint < 2; ++endpoint) if (PointInTriangle2D(edge_ids[endpoint], edge_points[endpoint], triangle_ids, triangle_points, drop, vertices)) {
				bool allowed = false; for (std::size_t k = 0; k < shared_count; ++k) if (edge_ids[endpoint] == shared[k]) allowed = true;
				if (!allowed) forbidden = true;
			}
		};
		for (std::size_t edge = 0; edge < 3; ++edge) {
			edge_in_other_plane(std::array<std::array<long double, 3>, 2>{{a[edge], a[(edge+1)%3]}},
				std::array<std::uint32_t, 2>{{left.indices[edge], left.indices[(edge+1)%3]}}, b, right.indices, nb);
			edge_in_other_plane(std::array<std::array<long double, 3>, 2>{{b[edge], b[(edge+1)%3]}},
				std::array<std::uint32_t, 2>{{right.indices[edge], right.indices[(edge+1)%3]}}, a, left.indices, na);
		}
		// A line segment contained in the other triangle's plane gives a zero
		// Moller determinant.  Test every edge pair separately so collinear
		// overlap on the line where the two triangle planes meet is still seen.
		for (std::size_t i = 0; i < 3; ++i) for (std::size_t j = 0; j < 3; ++j) {
			const auto ai = left.indices[i], aj = left.indices[(i+1)%3], bi = right.indices[j], bj = right.indices[(j+1)%3];
			if (Orient3D(a[i], a[(i+1)%3], b[j], b[(j+1)%3], ai, aj, bi, bj, vertices) != 0) continue;
			std::array<int, 3> cross_sign{};
			for (std::size_t axis = 0; axis < 3; ++axis) cross_sign[axis] = DirectionCrossComponent(a[i], a[(i+1)%3], b[j], b[(j+1)%3], ai, aj, bi, bj, axis, vertices);
			std::size_t drop = 0;
			const auto direction_cross = CrossL(SubtractL(a[(i+1)%3], a[i]), SubtractL(b[(j+1)%3], b[j]));
			if (cross_sign[0] != 0 || cross_sign[1] != 0 || cross_sign[2] != 0) {
				for (std::size_t axis = 1; axis < 3; ++axis) if (cross_sign[axis] != 0
					&& (cross_sign[drop] == 0 || std::fabs(direction_cross[axis]) > std::fabs(direction_cross[drop]))) drop = axis;
			} else {
				std::array<int, 3> plane_sign{};
				for (std::size_t axis = 0; axis < 3; ++axis) plane_sign[axis] = DirectionCrossComponent(a[i], a[(i+1)%3], a[i], b[j], ai, aj, ai, bi, axis, vertices);
				const auto plane_normal = CrossL(SubtractL(a[(i+1)%3], a[i]), SubtractL(b[j], a[i]));
				if (plane_sign[0] != 0 || plane_sign[1] != 0 || plane_sign[2] != 0) {
					for (std::size_t axis = 1; axis < 3; ++axis) if (plane_sign[axis] != 0
						&& (plane_sign[drop] == 0 || std::fabs(plane_normal[axis]) > std::fabs(plane_normal[drop]))) drop = axis;
				} else {
					const auto direction = SubtractL(a[(i+1)%3], a[i]); std::size_t longest = 0;
					if (std::fabs(direction[1]) > std::fabs(direction[longest])) longest = 1;
					if (std::fabs(direction[2]) > std::fabs(direction[longest])) longest = 2;
					drop = longest == 0 ? 1 : 0;
				}
			}
			const auto contact = SegmentContact2D(a[i], a[(i+1)%3], b[j], b[(j+1)%3], ai, aj, bi, bj, drop, vertices);
			if (contact != Segment2DContact::none && !IsAllowedEdgeContact(contact, ai, aj, bi, bj, shared, shared_count)) forbidden = true;
		}
		return forbidden;
	}
	static SelfIntersectionCounts ValidateSelfIntersection(const std::vector<std::array<double, 3>>& vertices,
		const std::vector<WorkingTriangle>& triangles, const SurfaceAabb& bounds, double length)
	{
		std::vector<SurfaceAabb> triangle_bounds; triangle_bounds.reserve(triangles.size());
		for (const auto& triangle : triangles) triangle_bounds.push_back(TriangleBounds(triangle, vertices));
		std::vector<std::size_t> ids(triangles.size()); for (std::size_t i = 0; i < ids.size(); ++i) ids[i] = i;
		if (triangles.size() > std::numeric_limits<std::size_t>::max()/2)
			throw std::overflow_error("surface BVH node reserve overflows");
		std::vector<BvhNode> nodes; nodes.reserve(2*triangles.size()); BuildBvh(nodes, ids, triangle_bounds, 0, ids.size());
		SelfIntersectionCounts result; const double tolerance = 256.0*std::numeric_limits<double>::epsilon()*length;
		auto increment = [](std::size_t& count) { if (count == std::numeric_limits<std::size_t>::max()) throw std::overflow_error("surface BVH candidate count overflows"); ++count; };
		auto test_pair = [&](std::size_t first, std::size_t second) {
			increment(result.broad_candidates);
			if (!BoundsOverlap(triangle_bounds[first], triangle_bounds[second], tolerance)) return;
			increment(result.narrow_candidates);
			if (PairIntersectsForbidden(triangles[first], triangles[second], vertices, bounds, length)) {
				const auto& a = vertices[triangles[first].indices[0]];
				const auto& b = vertices[triangles[second].indices[0]];
				throw std::invalid_argument("surface self-intersection detected between canonical triangles "
					+std::to_string(first)+" and "+std::to_string(second)
					+" near ("+std::to_string(a[0])+","+std::to_string(a[1])+","+std::to_string(a[2])+")"
					+" and ("+std::to_string(b[0])+","+std::to_string(b[1])+","+std::to_string(b[2])+") m");
			}
		};
		std::function<void(std::size_t, std::size_t)> visit = [&](std::size_t left, std::size_t right) {
			if (!BoundsOverlap(nodes[left].bounds, nodes[right].bounds, tolerance)) return;
			const auto& a = nodes[left]; const auto& b = nodes[right];
			if (left == right) { if (a.leaf) { for (std::size_t i = a.begin; i < a.end; ++i) for (std::size_t j = i+1; j < a.end; ++j) test_pair(ids[i], ids[j]); } else { visit(a.left, a.left); visit(a.left, a.right); visit(a.right, a.right); } return; }
			if (a.leaf && b.leaf) { for (std::size_t i = a.begin; i < a.end; ++i) for (std::size_t j = b.begin; j < b.end; ++j) test_pair(ids[i], ids[j]); return; }
			if (b.leaf || (!a.leaf && a.end-a.begin >= b.end-b.begin)) { visit(a.left, right); visit(a.right, right); } else { visit(left, b.left); visit(left, b.right); }
		};
		visit(0, 0); return result;
	}
	static double MinimumEdgeLength(const std::array<std::uint32_t, 3>& indices, const std::vector<std::array<double, 3>>& vertices)
	{
		return std::min({Norm(Subtract(vertices[indices[1]], vertices[indices[0]])), Norm(Subtract(vertices[indices[2]], vertices[indices[1]])), Norm(Subtract(vertices[indices[0]], vertices[indices[2]]))});
	}
	static void PopulateGeometry(ClosedSurfaceTriangle& triangle, const std::vector<std::array<double, 3>>& vertices)
	{
		const auto& a = vertices[triangle.indices[0]]; const auto& b = vertices[triangle.indices[1]]; const auto& c = vertices[triangle.indices[2]];
		const auto cross = Cross(Subtract(b, a), Subtract(c, a)); const double magnitude = Norm(cross);
		triangle.area_m2 = 0.5*magnitude;
		triangle.outward_unit_normal = {{cross[0]/magnitude, cross[1]/magnitude, cross[2]/magnitude}};
		triangle.bounds.minimum = a; triangle.bounds.maximum = a;
		for (const auto& point : {b, c}) for (std::size_t axis = 0; axis < 3; ++axis) { triangle.bounds.minimum[axis] = std::min(triangle.bounds.minimum[axis], point[axis]); triangle.bounds.maximum[axis] = std::max(triangle.bounds.maximum[axis], point[axis]); }
	}
	static std::string Hash(const std::vector<std::array<double, 3>>& vertices, const std::vector<ClosedSurfaceTriangle>& triangles)
	{
		Sha256 hash; static constexpr char version[] = "TubularFlowIGA.ClosedTriangulatedSurface.v1";
		hash.Append(version, sizeof(version)-1); hash.AppendLittleEndian64(vertices.size()); hash.AppendLittleEndian64(triangles.size());
		for (const auto& vertex : vertices) for (const auto coordinate : vertex) hash.AppendNormalizedDouble(coordinate);
		for (const auto& triangle : triangles) { for (const auto index : triangle.indices) hash.AppendLittleEndian32(index); hash.AppendLittleEndian32(triangle.boundary_id); }
		return hash.Hex();
	}

	std::vector<std::array<double, 3>> vertices_;
	std::vector<ClosedSurfaceTriangle> triangles_;
	SurfaceAabb bounds_;
	SurfaceDiagnostics diagnostics_;
};

} // namespace iga

#endif
