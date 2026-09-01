#ifndef IGA_SURFACE_GEOMETRY_HPP
#define IGA_SURFACE_GEOMETRY_HPP

#include "Sha256.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
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
	std::size_t vertex_count = 0;
	std::size_t triangle_count = 0;
	std::size_t component_count = 0;
	std::size_t unique_edge_count = 0;
	std::size_t vertex_fan_count = 0;
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

		std::vector<std::array<double, 3>> vertices = scaled;
		std::sort(vertices.begin(), vertices.end(), LexicographicLess);
		vertices.erase(std::unique(vertices.begin(), vertices.end()), vertices.end());
		std::vector<std::uint32_t> remap(scaled.size());
		for (std::size_t source = 0; source < scaled.size(); ++source) {
			const auto found = std::lower_bound(vertices.begin(), vertices.end(), scaled[source], LexicographicLess);
			remap[source] = static_cast<std::uint32_t>(found-vertices.begin());
		}

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
		result.diagnostics_.triangle_count = result.triangles_.size();
		result.diagnostics_.component_count = components;
		result.diagnostics_.unique_edge_count = edge_uses.size();
		result.diagnostics_.vertex_fan_count = result.vertices_.size();
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
		if (options.default_boundary_id < 0 || options.max_triangles < 0)
			throw std::invalid_argument("surface validation ids and counts must be nonnegative");
		if (static_cast<std::uintmax_t>(options.default_boundary_id)
			> static_cast<std::uintmax_t>(std::numeric_limits<std::uint32_t>::max()))
			throw std::invalid_argument("surface default boundary id exceeds canonical range");
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
