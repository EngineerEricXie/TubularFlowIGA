#ifndef IGA_PRESCRIBED_SURFACE_MOTION_HPP
#define IGA_PRESCRIBED_SURFACE_MOTION_HPP

#include "SurfaceGeometry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

// Source ordering is deliberately separate from ClosedTriangulatedSurface's
// canonical ordering.  A moving surface must retain its material vertices and
// directed source facets; canonicalizing independently at each frame would
// destroy that correspondence.
struct PrescribedSurfaceFrame {
	double time_s = 0.0;
	RawSurfaceSoup surface;
};

struct PrescribedSurfaceMotionOptions {
	SurfaceValidationOptions surface_validation;
	bool require_containment_in_fixed_bounds = false;
	SurfaceAabb fixed_bounds_m;
	double maximum_displacement_m = std::numeric_limits<double>::infinity();
	double maximum_velocity_m_per_s = std::numeric_limits<double>::infinity();
	// An interval displacement may not exceed this fraction of the configured
	// fixed-background extension band.  Infinity disables the check.
	double extension_band_m = std::numeric_limits<double>::infinity();
	double maximum_extension_band_cfl = 1.0;
};

struct SourceTriangleProvenance {
	std::uint32_t source_triangle = 0;
	std::array<std::uint32_t, 3> source_vertex_indices{{0, 0, 0}};
	std::uint32_t boundary_id = 0;
	// For canonical-triangle provenance only, maps a canonical triangle corner
	// to the corresponding source/material triangle corner.
	std::array<std::uint32_t, 3> canonical_corner_to_source_corner{{0, 1, 2}};
	bool operator==(const SourceTriangleProvenance& other) const
	{
		return source_triangle == other.source_triangle && source_vertex_indices == other.source_vertex_indices
			&& boundary_id == other.boundary_id
			&& canonical_corner_to_source_corner == other.canonical_corner_to_source_corner;
	}
};

class PrescribedSurfaceMotion {
public:
	struct Evaluation {
		Evaluation(const Evaluation&) = default;
		Evaluation(Evaluation&&) noexcept = default;
		Evaluation& operator=(const Evaluation&) = delete;
		Evaluation& operator=(Evaluation&&) = delete;

		const ClosedTriangulatedSurface& Surface() const noexcept { return surface_; }
		const std::vector<std::array<double, 3>>& SourceVerticesM() const noexcept { return source_vertices_m_; }
		const std::vector<std::array<double, 3>>& SourceVertexVelocitiesMPerS() const noexcept
		{ return source_vertex_velocities_m_per_s_; }
		// Indexed in canonical surface triangle order.  It provides the unique
		// material/source facet and its canonical-to-material corner permutation.
		const std::vector<SourceTriangleProvenance>& CanonicalTriangleProvenance() const noexcept
		{ return canonical_triangle_provenance_; }
		const std::string& IdentitySha256() const noexcept { return identity_sha256_; }
		double EvaluatedTimeS() const noexcept { return evaluated_time_s_; }

		// The input triangle and barycentric weights use the canonical surface
		// triangle order returned by Surface().
		std::array<double, 3> WallVelocity(std::uint32_t canonical_triangle,
			const std::array<double, 3>& barycentric) const
		{
			if (canonical_triangle >= surface_.Triangles().size()
				|| canonical_triangle >= canonical_triangle_provenance_.size())
				throw std::out_of_range("canonical triangle is out of range");
			const auto& provenance = canonical_triangle_provenance_[canonical_triangle];
			if (provenance.source_triangle >= source_triangles_.size())
				throw std::out_of_range("canonical triangle source provenance is out of range");
			const auto& source_triangle = source_triangles_[provenance.source_triangle];
			if (source_triangle.source_triangle != provenance.source_triangle
				|| source_triangle.source_vertex_indices != provenance.source_vertex_indices
				|| source_triangle.boundary_id != provenance.boundary_id)
				throw std::invalid_argument("canonical triangle source provenance is inconsistent");
			constexpr double barycentric_tolerance = 64.0*std::numeric_limits<double>::epsilon();
			double weight_sum = 0.0;
			for (double weight : barycentric) {
				if (!std::isfinite(weight)) throw std::invalid_argument("barycentric weight is nonfinite");
				if (weight < -barycentric_tolerance || weight > 1.0+barycentric_tolerance)
					throw std::invalid_argument("barycentric weight is outside the simplex");
				weight_sum += weight;
				if (!std::isfinite(weight_sum)) throw std::invalid_argument("barycentric weight sum is nonfinite");
			}
			if (std::fabs(weight_sum-1.0) > barycentric_tolerance)
				throw std::invalid_argument("barycentric weights do not sum to one");
			std::array<double, 3> result{{0.0, 0.0, 0.0}};
			for (std::size_t canonical_corner = 0; canonical_corner < 3; ++canonical_corner) {
				const std::uint32_t source_corner = provenance.canonical_corner_to_source_corner[canonical_corner];
				if (source_corner >= 3) throw std::out_of_range("canonical triangle corner provenance is out of range");
				const std::uint32_t source_vertex = provenance.source_vertex_indices[source_corner];
				if (source_vertex >= source_vertex_velocities_m_per_s_.size())
					throw std::out_of_range("source vertex velocity is out of range");
				for (std::size_t axis = 0; axis < 3; ++axis)
					result[axis] += barycentric[canonical_corner]
						*source_vertex_velocities_m_per_s_[source_vertex][axis];
			}
			for (double value : result)
				if (!std::isfinite(value)) throw std::invalid_argument("wall velocity is nonfinite");
			return result;
		}

	private:
		friend class PrescribedSurfaceMotion;
		explicit Evaluation(ClosedTriangulatedSurface checked) : surface_(std::move(checked)) {}
		double evaluated_time_s_ = 0.0;
		ClosedTriangulatedSurface surface_;
		std::vector<std::array<double, 3>> source_vertices_m_;
		std::vector<std::array<double, 3>> source_vertex_velocities_m_per_s_;
		std::vector<SourceTriangleProvenance> canonical_triangle_provenance_;
		std::string identity_sha256_;
		std::vector<SourceTriangleProvenance> source_triangles_;
	};

	explicit PrescribedSurfaceMotion(std::vector<PrescribedSurfaceFrame> frames,
		const PrescribedSurfaceMotionOptions& options = PrescribedSurfaceMotionOptions())
		: frames_(std::move(frames)), options_(options)
	{
		ValidateOptions();
		if (frames_.size() < 2) throw std::invalid_argument("prescribed surface motion requires at least two frames");
		if (frames_.front().surface.vertices.size() > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())
			|| frames_.front().surface.triangles.size() > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()))
			throw std::overflow_error("prescribed surface source count exceeds uint32 provenance range");
		for (std::size_t frame = 0; frame < frames_.size(); ++frame) {
			if (!std::isfinite(frames_[frame].time_s)) throw std::invalid_argument("prescribed surface frame time is nonfinite");
			if (frame && !(frames_[frame].time_s > frames_[frame-1].time_s))
				throw std::invalid_argument("prescribed surface frame times are not strictly increasing");
			ValidateSourceTopology(frames_[frame].surface);
			const auto checked = BuildValidated(frames_[frame].surface);
			ValidateContainment(frames_[frame].surface.vertices);
			ValidateSourceCanMapToCanonical(frames_[frame].surface, checked);
		}
		BuildSourceProvenance();
		ValidateFrameDurations();
		ValidateIntervalLimits();
	}

	const std::vector<PrescribedSurfaceFrame>& Frames() const noexcept { return frames_; }
	const std::vector<SourceTriangleProvenance>& SourceTriangles() const noexcept { return source_triangles_; }

	Evaluation Evaluate(double time_s, double step_start_s, double step_end_s) const
	{
		if (!std::isfinite(time_s) || !std::isfinite(step_start_s) || !std::isfinite(step_end_s))
			throw std::invalid_argument("prescribed surface trial time is nonfinite");
		if (!(step_end_s > step_start_s) || time_s < step_start_s || time_s > step_end_s)
			throw std::invalid_argument("prescribed surface trial time is outside its positive step");
		(void)PositiveFiniteDifference(step_start_s, step_end_s, "prescribed surface step duration is invalid");
		const std::size_t interval = StepInterval(step_start_s, step_end_s);
		const double left_time = frames_[interval].time_s;
		const double duration = FrameDuration(interval);
		const double theta = time_s == left_time ? 0.0 : time_s == frames_[interval+1].time_s ? 1.0
			: PositiveFiniteDifference(left_time, time_s, "prescribed surface interpolation time is invalid")/duration;
		if (!std::isfinite(theta) || theta < 0.0 || theta > 1.0)
			throw std::invalid_argument("prescribed surface interpolation parameter is invalid");

		std::vector<std::array<double, 3>> positions;
		std::vector<std::array<double, 3>> velocities;
		positions.reserve(frames_.front().surface.vertices.size());
		velocities.reserve(frames_.front().surface.vertices.size());
		for (std::size_t vertex = 0; vertex < frames_.front().surface.vertices.size(); ++vertex) {
			std::array<double, 3> position{};
			std::array<double, 3> velocity{};
			for (std::size_t axis = 0; axis < 3; ++axis) {
				const double first = frames_[interval].surface.vertices[vertex][axis];
				const double second = frames_[interval+1].surface.vertices[vertex][axis];
				const long double delta = FiniteDifference(first, second,
					"prescribed surface coordinate difference is invalid");
				const long double interpolated = static_cast<long double>(first)+delta*static_cast<long double>(theta);
				const long double velocity_value = delta/static_cast<long double>(duration);
				if (!std::isfinite(interpolated) || !std::isfinite(velocity_value)
					|| std::fabs(interpolated) > std::numeric_limits<double>::max()
					|| std::fabs(velocity_value) > std::numeric_limits<double>::max())
					throw std::invalid_argument("prescribed surface interpolation is unrepresentable");
				position[axis] = theta == 0.0 ? first : theta == 1.0 ? second : static_cast<double>(interpolated);
				velocity[axis] = static_cast<double>(velocity_value);
				if (!std::isfinite(position[axis]) || !std::isfinite(velocity[axis]))
					throw std::invalid_argument("prescribed surface interpolation is nonfinite");
			}
			positions.push_back(position);
			velocities.push_back(velocity);
		}
		RawSurfaceSoup soup;
		soup.vertices = positions;
		soup.triangles = frames_.front().surface.triangles;
		Evaluation result(BuildValidated(soup));
		result.evaluated_time_s_ = time_s;
		result.source_vertices_m_ = std::move(positions);
		result.source_vertex_velocities_m_per_s_ = std::move(velocities);
		ValidateContainment(result.source_vertices_m_);
		result.canonical_triangle_provenance_ = CanonicalProvenance(soup, result.surface_);
		result.source_triangles_ = source_triangles_;
		result.identity_sha256_ = HashEvaluation(result);
		return result;
	}

private:
	void ValidateOptions() const
	{
		if (options_.surface_validation.length_scale_to_m != 1.0 || options_.surface_validation.weld_tolerance_m != 0.0)
			throw std::invalid_argument("prescribed surface motion requires metre coordinates and zero welding");
		for (double value : {options_.maximum_displacement_m, options_.maximum_velocity_m_per_s,
			options_.extension_band_m})
			if (std::isnan(value) || value < 0.0) throw std::invalid_argument("prescribed surface motion limit is invalid");
		if (!std::isfinite(options_.maximum_extension_band_cfl) || options_.maximum_extension_band_cfl < 0.0)
			throw std::invalid_argument("prescribed surface extension-band CFL is invalid");
		if (options_.require_containment_in_fixed_bounds)
			for (std::size_t axis = 0; axis < 3; ++axis)
				if (!std::isfinite(options_.fixed_bounds_m.minimum[axis]) || !std::isfinite(options_.fixed_bounds_m.maximum[axis])
					|| options_.fixed_bounds_m.minimum[axis] > options_.fixed_bounds_m.maximum[axis])
					throw std::invalid_argument("prescribed surface fixed bounds are invalid");
	}

	void ValidateSourceTopology(const RawSurfaceSoup& surface) const
	{
		const auto& source = frames_.front().surface;
		if (surface.vertices.size() != source.vertices.size() || surface.triangles.size() != source.triangles.size())
			throw std::invalid_argument("prescribed surface source counts differ between frames");
		for (const auto& triangle : surface.triangles)
			if (triangle.boundary_id < 0) throw std::invalid_argument("prescribed surface requires explicit nonnegative labels");
		for (std::size_t triangle = 0; triangle < source.triangles.size(); ++triangle)
			if (surface.triangles[triangle].indices != source.triangles[triangle].indices
				|| surface.triangles[triangle].boundary_id != source.triangles[triangle].boundary_id)
				throw std::invalid_argument("prescribed surface source connectivity, labels, or winding differs between frames");
	}

	ClosedTriangulatedSurface BuildValidated(const RawSurfaceSoup& soup) const
	{
		auto result = ClosedTriangulatedSurface::Build(soup, options_.surface_validation);
		if (result.Diagnostics().flipped_inward_shell)
			throw std::invalid_argument("prescribed surface source winding is inward");
		return result;
	}

	void ValidateContainment(const std::vector<std::array<double, 3>>& vertices) const
	{
		if (!options_.require_containment_in_fixed_bounds) return;
		for (const auto& vertex : vertices)
			for (std::size_t axis = 0; axis < 3; ++axis)
				if (vertex[axis] < options_.fixed_bounds_m.minimum[axis]
					|| vertex[axis] > options_.fixed_bounds_m.maximum[axis])
					throw std::invalid_argument("prescribed surface leaves fixed Cartesian bounds");
	}

	void ValidateSourceCanMapToCanonical(const RawSurfaceSoup& source,
		const ClosedTriangulatedSurface& canonical) const
	{
		if (source.vertices.size() != canonical.Vertices().size())
			throw std::invalid_argument("prescribed surface source vertices are not one-to-one");
		(void)CanonicalProvenance(source, canonical);
	}

	void BuildSourceProvenance()
	{
		source_triangles_.reserve(frames_.front().surface.triangles.size());
		for (std::size_t triangle = 0; triangle < frames_.front().surface.triangles.size(); ++triangle) {
			const auto& raw = frames_.front().surface.triangles[triangle];
			SourceTriangleProvenance result;
			result.source_triangle = static_cast<std::uint32_t>(triangle);
			result.boundary_id = static_cast<std::uint32_t>(raw.boundary_id);
			for (std::size_t corner = 0; corner < 3; ++corner)
				result.source_vertex_indices[corner] = static_cast<std::uint32_t>(raw.indices[corner]);
			source_triangles_.push_back(result);
		}
	}

	void ValidateIntervalLimits() const
	{
		for (std::size_t frame = 0; frame+1 < frames_.size(); ++frame) {
			const double duration = FrameDuration(frame);
			for (std::size_t vertex = 0; vertex < frames_[frame].surface.vertices.size(); ++vertex) {
				std::array<long double, 3> delta{};
				long double scale = 0.0L;
				for (std::size_t axis = 0; axis < 3; ++axis) {
					delta[axis] = FiniteDifference(frames_[frame].surface.vertices[vertex][axis],
						frames_[frame+1].surface.vertices[vertex][axis],
						"prescribed surface coordinate difference is invalid");
					scale = std::max(scale, std::fabs(delta[axis]));
				}
				long double displacement_value = 0.0L;
				if (scale != 0.0L) {
					long double normalized_squared = 0.0L;
					for (long double component : delta) normalized_squared += (component/scale)*(component/scale);
					displacement_value = scale*std::sqrt(normalized_squared);
				}
				const long double velocity_value = displacement_value/static_cast<long double>(duration);
				if (!std::isfinite(displacement_value) || !std::isfinite(velocity_value)
					|| displacement_value > std::numeric_limits<double>::max()
					|| velocity_value > std::numeric_limits<double>::max())
					throw std::invalid_argument("prescribed surface displacement or velocity is unrepresentable");
				const double displacement = static_cast<double>(displacement_value);
				const double velocity = static_cast<double>(velocity_value);
				if (!std::isfinite(displacement) || !std::isfinite(velocity)
					|| displacement > options_.maximum_displacement_m || velocity > options_.maximum_velocity_m_per_s)
					throw std::invalid_argument("prescribed surface displacement or velocity exceeds configured limit");
				if (std::isfinite(options_.extension_band_m)) {
					if (options_.extension_band_m == 0.0 || options_.maximum_extension_band_cfl == 0.0) {
						if (displacement != 0.0)
							throw std::invalid_argument("prescribed surface motion exceeds extension-band CFL limit");
					} else {
						// With a positive finite CFL, a quotient that underflows is below
						// every representable positive CFL; an overflow exceeds every one.
						const double extension_ratio = displacement/options_.extension_band_m;
						if (!std::isfinite(extension_ratio) || extension_ratio > options_.maximum_extension_band_cfl)
							throw std::invalid_argument("prescribed surface motion exceeds extension-band CFL limit");
					}
				}
			}
		}
	}

	static long double FiniteDifference(double left, double right, const char* message)
	{
		const long double difference = static_cast<long double>(right)-static_cast<long double>(left);
		if (!std::isfinite(difference)) throw std::invalid_argument(message);
		return difference;
	}

	static double PositiveFiniteDifference(double left, double right, const char* message)
	{
		const long double difference = FiniteDifference(left, right, message);
		if (!(difference > 0.0L) || difference > std::numeric_limits<double>::max())
			throw std::invalid_argument(message);
		return static_cast<double>(difference);
	}

	void ValidateFrameDurations() const
	{
		for (std::size_t frame = 0; frame+1 < frames_.size(); ++frame)
			(void)FrameDuration(frame);
	}

	double FrameDuration(std::size_t frame) const
	{
		if (frame+1 >= frames_.size()) throw std::out_of_range("prescribed surface frame interval is out of range");
		return PositiveFiniteDifference(frames_[frame].time_s, frames_[frame+1].time_s,
			"prescribed surface frame duration is invalid");
	}

	std::size_t StepInterval(double start, double end) const
	{
		if (start < frames_.front().time_s || end > frames_.back().time_s)
			throw std::invalid_argument("prescribed surface step is outside frame time range");
		auto upper = std::upper_bound(frames_.begin(), frames_.end(), end,
			[](double value, const PrescribedSurfaceFrame& frame) { return value < frame.time_s; });
		std::size_t right = static_cast<std::size_t>(upper-frames_.begin());
		if (right && right <= frames_.size() && end == frames_[right-1].time_s) --right;
		if (!right || right >= frames_.size()) throw std::invalid_argument("prescribed surface step has no interpolation interval");
		const std::size_t interval = right-1;
		if (start < frames_[interval].time_s || end > frames_[interval+1].time_s)
			throw std::invalid_argument("prescribed surface step spans an interpolation knot");
		return interval;
	}

	std::vector<SourceTriangleProvenance> CanonicalProvenance(const RawSurfaceSoup& source,
		const ClosedTriangulatedSurface& canonical) const
	{
		std::map<std::array<double, 3>, std::uint32_t> vertex_ids;
		for (std::size_t vertex = 0; vertex < canonical.Vertices().size(); ++vertex)
			vertex_ids.emplace(canonical.Vertices()[vertex], static_cast<std::uint32_t>(vertex));
		struct FacetMapping {
			SourceTriangleProvenance provenance;
			std::array<std::uint32_t, 3> canonical_vertex_by_source_corner{};
		};
		std::map<std::array<std::uint32_t, 3>, FacetMapping> facets;
		for (std::size_t triangle = 0; triangle < source.triangles.size(); ++triangle) {
			FacetMapping mapping;
			mapping.provenance.source_triangle = static_cast<std::uint32_t>(triangle);
			mapping.provenance.boundary_id = static_cast<std::uint32_t>(source.triangles[triangle].boundary_id);
			std::array<std::uint32_t, 3> canonical_indices{};
			for (std::size_t corner = 0; corner < 3; ++corner) {
				mapping.provenance.source_vertex_indices[corner] = static_cast<std::uint32_t>(source.triangles[triangle].indices[corner]);
				const auto found = vertex_ids.find(source.vertices[mapping.provenance.source_vertex_indices[corner]]);
				if (found == vertex_ids.end()) throw std::invalid_argument("prescribed source vertex cannot map to canonical surface");
				mapping.canonical_vertex_by_source_corner[corner] = found->second;
				canonical_indices[corner] = found->second;
			}
			const auto smallest = std::min_element(canonical_indices.begin(), canonical_indices.end());
			std::rotate(canonical_indices.begin(), smallest, canonical_indices.end());
			if (!facets.emplace(canonical_indices, mapping).second)
				throw std::invalid_argument("prescribed source triangle is not uniquely mappable");
		}
		std::vector<SourceTriangleProvenance> result;
		result.reserve(canonical.Triangles().size());
		for (const auto& triangle : canonical.Triangles()) {
			const auto found = facets.find(triangle.indices);
			if (found == facets.end() || found->second.provenance.boundary_id != triangle.boundary_id)
				throw std::invalid_argument("prescribed source triangle winding differs from validated surface");
			auto provenance = found->second.provenance;
			provenance.canonical_corner_to_source_corner = {{0, 1, 2}};
			for (std::size_t canonical_corner = 0; canonical_corner < 3; ++canonical_corner) {
				const auto source_corner = std::find(found->second.canonical_vertex_by_source_corner.begin(),
					found->second.canonical_vertex_by_source_corner.end(), triangle.indices[canonical_corner]);
				if (source_corner == found->second.canonical_vertex_by_source_corner.end())
					throw std::invalid_argument("prescribed source triangle cannot map canonical corner");
				provenance.canonical_corner_to_source_corner[canonical_corner] =
					static_cast<std::uint32_t>(source_corner-found->second.canonical_vertex_by_source_corner.begin());
			}
			result.push_back(provenance);
		}
		return result;
	}

	static void AppendCount(Sha256& hash, std::size_t count)
	{
		if (count > static_cast<std::size_t>(std::numeric_limits<std::uint64_t>::max()))
			throw std::overflow_error("prescribed surface digest count is unrepresentable");
		hash.AppendLittleEndian64(static_cast<std::uint64_t>(count));
	}

	static void AppendString(Sha256& hash, const std::string& value)
	{
		AppendCount(hash, value.size());
		hash.Append(value.data(), value.size());
	}

	std::string HashEvaluation(const Evaluation& evaluation) const
	{
		Sha256 hash;
		static constexpr char domain[] = "iga-prescribed-surface-motion-v2";
		hash.Append(domain, sizeof(domain)-1);
		hash.AppendNormalizedDouble(evaluation.evaluated_time_s_);
		AppendCount(hash, evaluation.source_vertices_m_.size());
		for (const auto& vertex : evaluation.source_vertices_m_)
			for (double value : vertex) hash.AppendNormalizedDouble(value);
		AppendCount(hash, evaluation.source_vertex_velocities_m_per_s_.size());
		for (const auto& velocity : evaluation.source_vertex_velocities_m_per_s_)
			for (double value : velocity) hash.AppendNormalizedDouble(value);
		AppendCount(hash, source_triangles_.size());
		for (const auto& triangle : source_triangles_) {
			hash.AppendLittleEndian32(triangle.source_triangle);
			hash.AppendLittleEndian32(triangle.boundary_id);
			for (auto index : triangle.source_vertex_indices) hash.AppendLittleEndian32(index);
		}
		AppendCount(hash, evaluation.Surface().Vertices().size());
		for (const auto& vertex : evaluation.Surface().Vertices())
			for (double value : vertex) hash.AppendNormalizedDouble(value);
		AppendCount(hash, evaluation.Surface().Triangles().size());
		for (const auto& triangle : evaluation.Surface().Triangles()) {
			for (auto index : triangle.indices) hash.AppendLittleEndian32(index);
			hash.AppendLittleEndian32(triangle.boundary_id);
		}
		AppendString(hash, evaluation.Surface().CanonicalSha256());
		AppendCount(hash, evaluation.canonical_triangle_provenance_.size());
		for (const auto& triangle : evaluation.canonical_triangle_provenance_) {
			hash.AppendLittleEndian32(triangle.source_triangle);
			hash.AppendLittleEndian32(triangle.boundary_id);
			for (auto index : triangle.source_vertex_indices) hash.AppendLittleEndian32(index);
			for (auto corner : triangle.canonical_corner_to_source_corner) hash.AppendLittleEndian32(corner);
		}
		return hash.Hex();
	}

	std::vector<PrescribedSurfaceFrame> frames_;
	PrescribedSurfaceMotionOptions options_;
	std::vector<SourceTriangleProvenance> source_triangles_;
};

} // namespace iga

#endif
