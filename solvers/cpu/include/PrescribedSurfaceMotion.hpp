#ifndef IGA_PRESCRIBED_SURFACE_MOTION_HPP
#define IGA_PRESCRIBED_SURFACE_MOTION_HPP

#include "MaterialSurfaceKinematics.hpp"

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

inline double PrescribedClockRoundoffAllowance(std::uint64_t additions)
{
	// gamma_n bounds repeated nonnegative clock additions. Epsilon is twice
	// unit roundoff; the extra margin covers the endpoint comparisons.
	const long double error=static_cast<long double>(additions)*std::numeric_limits<double>::epsilon();
	if (!(error<1)) throw std::invalid_argument("prescribed clock accumulation bound is invalid");
	const double result=static_cast<double>(error/(1-error)+8*std::numeric_limits<double>::epsilon());
	if (!(result<.25)) throw std::invalid_argument("prescribed clock accumulation allowance is too large");
	return result;
}

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
	// Relative allowance for the caller's floating-point clock arithmetic.
	// Interpolation snapping is also capped by adjacent frame durations.
	double clock_roundoff_relative_tolerance = 8*std::numeric_limits<double>::epsilon();
};

class PrescribedSurfaceMotion {
public:
	// Compatibility spelling for existing prescribed-motion callers.  New
	// consumers take MaterialSurfaceKinematics directly.
	using Evaluation = MaterialSurfaceKinematics;

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
		double interpolation_start = SnapFrameTime(step_start_s), interpolation_end = SnapFrameTime(step_end_s);
		// Do not collapse a genuinely short, representable step near a knot.
		const bool snap = interpolation_start < interpolation_end;
		if (!snap) { interpolation_start=step_start_s; interpolation_end=step_end_s; }
		const std::size_t interval = StepInterval(interpolation_start, interpolation_end);
		const double interpolation_time = snap ? SnapFrameTime(time_s) : time_s;
		const double left_time = frames_[interval].time_s;
		const double duration = FrameDuration(interval);
		const double theta = interpolation_time == left_time ? 0.0 : interpolation_time == frames_[interval+1].time_s ? 1.0
			: PositiveFiniteDifference(left_time, interpolation_time, "prescribed surface interpolation time is invalid")/duration;
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
		result.step_start_s_ = step_start_s;
		result.step_end_s_ = step_end_s;
		result.reference_material_vertices_m_ = frames_.front().surface.vertices;
		result.source_vertices_m_ = std::move(positions);
		result.source_vertex_velocities_m_per_s_ = std::move(velocities);
		ValidateContainment(result.source_vertices_m_);
		result.canonical_triangle_provenance_ = CanonicalProvenance(soup, result.surface_);
		result.source_triangles_ = source_triangles_;
		result.identity_sha256_ = HashEvaluation(result);
		// Preserve the prescribed producer's reference-coordinate material
		// identity while the neutral owner independently recomputes it in Validate.
		result.material_identity_sha256_ = HashMaterialIdentity();
		result.topology_identity_sha256_ = HashTopologyIdentity();
		result.content_identity_sha256_ = result.HashContentIdentity();
		result.Validate();
		return result;
	}

private:
	void ValidateOptions() const
	{
		if (!std::isfinite(options_.clock_roundoff_relative_tolerance)
			|| options_.clock_roundoff_relative_tolerance<0 || options_.clock_roundoff_relative_tolerance>=.25)
			throw std::invalid_argument("prescribed clock roundoff allowance is invalid");
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
		(void)CanonicalProvenance(source, canonical);
		if (source.vertices.size() != canonical.Vertices().size())
			throw std::invalid_argument("prescribed surface source vertices are not one-to-one");
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

	// Only interpolation queries are snapped; published material clocks and
	// input frame times retain their original bits. The neighbour-duration cap
	// prevents the tolerance from merging closely spaced physical knots.
	double SnapFrameTime(double time) const
	{
		const auto next=std::lower_bound(frames_.begin(),frames_.end(),time,
			[](const PrescribedSurfaceFrame& frame,double value) { return frame.time_s<value; });
		const std::size_t right=static_cast<std::size_t>(next-frames_.begin());
		double result=time;long double best=std::numeric_limits<long double>::infinity();
		for (const auto index : {right,right ? right-1 : frames_.size()}) {
			if (index>=frames_.size()) continue;
			const double knot=frames_[index].time_s;
			long double tolerance=options_.clock_roundoff_relative_tolerance
				*std::max(std::abs(static_cast<long double>(time)),std::abs(static_cast<long double>(knot)));
			if (index) tolerance=std::min(tolerance,static_cast<long double>(FrameDuration(index-1))/4);
			if (index+1<frames_.size()) tolerance=std::min(tolerance,static_cast<long double>(FrameDuration(index))/4);
			const long double difference=std::abs(static_cast<long double>(time)-knot);
			if (difference<=tolerance && difference<best) { result=knot;best=difference; }
		}
		return result;
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
		const auto unmapped_vertex = std::numeric_limits<std::uint32_t>::max();
		std::vector<std::uint32_t> canonical_vertex_for_source(source.vertices.size(), unmapped_vertex);
		std::vector<std::uint32_t> source_vertex_for_canonical(canonical.Vertices().size(), unmapped_vertex);
		for (std::size_t source_vertex = 0; source_vertex < source.vertices.size(); ++source_vertex) {
			const auto found = vertex_ids.find(source.vertices[source_vertex]);
			if (found == vertex_ids.end())
				throw std::invalid_argument("prescribed source vertex cannot map to canonical surface");
			const auto canonical_vertex = found->second;
			if (source_vertex_for_canonical[canonical_vertex] != unmapped_vertex)
				throw std::invalid_argument("prescribed source/canonical vertex mapping is not one-to-one");
			canonical_vertex_for_source[source_vertex] = canonical_vertex;
			source_vertex_for_canonical[canonical_vertex] = static_cast<std::uint32_t>(source_vertex);
		}
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
				const auto source_vertex = mapping.provenance.source_vertex_indices[corner];
				mapping.canonical_vertex_by_source_corner[corner] = canonical_vertex_for_source[source_vertex];
				canonical_indices[corner] = canonical_vertex_for_source[source_vertex];
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
			const auto& permutation = provenance.canonical_corner_to_source_corner;
			const unsigned inversions = (permutation[0] > permutation[1]) + (permutation[0] > permutation[2])
				+ (permutation[1] > permutation[2]);
			if (inversions%2 != 0)
				throw std::invalid_argument("prescribed source triangle orientation is reversed");
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

	std::string HashTopologyIdentity() const
	{
		Sha256 hash;
		static constexpr char domain[] = "MaterialSurfaceKinematics/topology/v1";
		hash.Append(domain, sizeof(domain)-1);
		AppendCount(hash, source_triangles_.size());
		for (const auto& triangle : source_triangles_) {
			hash.AppendLittleEndian32(triangle.source_triangle);
			hash.AppendLittleEndian32(triangle.boundary_id);
			for (const auto index : triangle.source_vertex_indices) hash.AppendLittleEndian32(index);
		}
		return hash.Hex();
	}

	std::string HashMaterialIdentity() const
	{
		Sha256 hash;
		static constexpr char domain[] = "MaterialSurfaceKinematics/material/v1";
		hash.Append(domain, sizeof(domain)-1);
		AppendCount(hash, frames_.front().surface.vertices.size());
		for (const auto& vertex : frames_.front().surface.vertices)
			for (const auto value : vertex) hash.AppendNormalizedDouble(value);
		AppendString(hash, HashTopologyIdentity());
		return hash.Hex();
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
