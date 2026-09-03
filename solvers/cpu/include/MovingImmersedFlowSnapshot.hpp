#ifndef IGA_MOVING_IMMERSED_FLOW_SNAPSHOT_HPP
#define IGA_MOVING_IMMERSED_FLOW_SNAPSHOT_HPP

// Immutable, local post-processing for one already-committed moving immersed
// flow state.  This deliberately has no runtime/PETSc objects and performs no
// publication or file I/O.
#include "MovingCutGeometry.hpp"
#include "ImmersedTransientState.hpp"
#include "TransportElement.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

struct MovingImmersedFlowSnapshotOptions {
	double stagnant_speed_threshold_m_per_s = 0.0;
	std::size_t maximum_points = 1024u*1024u;
	std::size_t maximum_output_bytes = 256u*1024u*1024u;
};

struct MovingImmersedPortFlow {
	std::uint64_t label = 0;
	double outward_flow_m3_s = 0.0;
};

struct MovingImmersedFlowSnapshotRequest {
	double time_s = 0.0;
	std::uint64_t index = 0;
	// A transition permits a step-level, well-mixed replacement proxy.  With
	// no transition dt must be canonical zero and that proxy is zero.
	bool transition_available = false;
	double dt_s = 0.0;
	// These are the complete configured open-boundary labels, independent of
	// flow-controller rows (which exist only for flow-rate-controlled ports).
	std::vector<std::uint32_t> port_labels;
	// False means the endpoint port measurements are unavailable, not zero.
	// In that case port_flows must be empty.
	bool port_flows_available = false;
	std::vector<MovingImmersedPortFlow> port_flows;
	std::vector<std::uint32_t> wall_labels;
};

enum class MovingImmersedFieldCellKind : std::uint8_t { Inside = 0, Cut = 1 };

struct MovingImmersedFieldPoint {
	std::array<double, 3> physical_m{};
	std::array<double, 3> parametric{};
	std::array<double, 3> velocity_m_per_s{};
	std::array<double, 3> vorticity_per_s{};
	double pressure = 0.0;
	double speed_m_per_s = 0.0;
	double q_criterion_per_s2 = 0.0;
	double enstrophy_density_per_s2 = 0.0;
	double physical_integration_weight_m3 = 0.0;
	std::uint64_t background_cell_id = 0;
	std::uint64_t cell_quadrature_ordinal = 0;
	MovingImmersedFieldCellKind cell_kind = MovingImmersedFieldCellKind::Inside;
	bool stagnant = false;
};

struct MovingImmersedFlowSnapshotMetrics {
	double quadrature_volume_m3 = 0.0;
	double audited_volume_m3 = 0.0;
	double enstrophy_integral_m3_per_s2 = 0.0;
	double mean_enstrophy_per_s2 = 0.0;
	double mean_q_criterion_per_s2 = 0.0;
	double q_positive_volume_m3 = 0.0;
	double q_positive_volume_fraction = 0.0;
	double stagnant_volume_m3 = 0.0;
	double stagnant_volume_fraction = 0.0;
	double inlet_flow_m3_s = 0.0;
	double outlet_flow_m3_s = 0.0;
	bool endpoint_turnover_available = false;
	double endpoint_turnover_rate_per_s = 0.0;
	std::optional<double> endpoint_turnover_time_s;
	double well_mixed_replacement_fraction_over_step = 0.0;
	// Integral of |u-w|^2 over wall area: (m/s)^2 m^2 = m^4/s^2.
	double wall_relative_velocity_squared_area_integral_m4_per_s2 = 0.0;
	double wall_relative_velocity_rms_m_per_s = 0.0;
	double wall_relative_velocity_max_m_per_s = 0.0;
	double wall_area_m2 = 0.0;
};

class MovingImmersedFlowSnapshot {
public:
	static MovingImmersedFlowSnapshot Build(const MovingCutGeometry& geometry,
		const ImmersedActiveLayout& layout, const ImmersedGlobalFlowState& state,
		MovingImmersedFlowSnapshotRequest request = {}, MovingImmersedFlowSnapshotOptions options = {})
	{
		ValidateIdentity(geometry, layout, state, request, options);
		const auto count = CountPoints(geometry, layout, options);
		MovingImmersedFlowSnapshot result;
		result.options_ = options; result.request_ = std::move(request); result.storage_mode_ = geometry.Volume().StorageMode();
		result.points_.reserve(count);
		long double volume = 0.0L, enstrophy = 0.0L, q_sum = 0.0L;
		long double q_positive = 0.0L, stagnant = 0.0L;
		for (std::uint64_t id = 0; id < geometry.Domain().Cells().size(); ++id) {
			const auto& classified = geometry.Domain().Cells()[static_cast<std::size_t>(id)];
			if (!Usable(classified, geometry.Volume().Cell(id))) continue;
			const auto element = geometry.Domain().Background().MaterializeElement(id);
			std::uint64_t ordinal = 0;
			const auto add = [&](const VolumeQuadraturePoint& point) {
				const auto basis = EvaluateBasis(element, point.parametric[0], point.parametric[1], point.parametric[2], false);
				const auto mapped = EvaluateElementGeometry(element, point.parametric);
				MovingImmersedFieldPoint item;
				item.physical_m = mapped.physical; item.parametric = point.parametric;
				item.background_cell_id = id; item.cell_quadrature_ordinal = ordinal++;
				item.cell_kind = classified.classification == CellClassification::Inside
					? MovingImmersedFieldCellKind::Inside : MovingImmersedFieldCellKind::Cut;
				const double weight = point.weight*basis.raw_determinant;
				RequireFinitePositive(weight, "moving immersed snapshot volume weight is invalid");
				item.physical_integration_weight_m3 = weight;
				std::array<std::array<double, 3>, 3> gradient{};
				for (std::size_t a = 0; a < element.connectivity.size(); ++a) {
					const auto local = layout.LocalNode(element.connectivity[a]);
					const auto& coefficient = state.Coefficients()[local];
					for (std::size_t i = 0; i < 3; ++i) {
						item.velocity_m_per_s[i] += basis.value[a]*coefficient[i];
						for (std::size_t j = 0; j < 3; ++j) gradient[i][j] += coefficient[i]*basis.gradient[a][j];
					}
					item.pressure += basis.value[a]*coefficient[3];
				}
				item.vorticity_per_s = {{gradient[2][1]-gradient[1][2], gradient[0][2]-gradient[2][0], gradient[1][0]-gradient[0][1]}};
				item.speed_m_per_s = Norm(item.velocity_m_per_s);
				double s2 = 0.0, w2 = 0.0;
				for (std::size_t i = 0; i < 3; ++i) for (std::size_t j = 0; j < 3; ++j) {
					const double symmetric = 0.5*(gradient[i][j]+gradient[j][i]);
					const double skew = 0.5*(gradient[i][j]-gradient[j][i]);
					s2 += symmetric*symmetric; w2 += skew*skew;
				}
				item.q_criterion_per_s2 = 0.5*(w2-s2);
				item.enstrophy_density_per_s2 = 0.5*Dot(item.vorticity_per_s, item.vorticity_per_s);
				item.stagnant = item.speed_m_per_s < options.stagnant_speed_threshold_m_per_s;
				ValidatePoint(item);
				AddFinite(volume, weight, "moving immersed snapshot volume sum is nonfinite");
				AddFinite(enstrophy, weight*item.enstrophy_density_per_s2, "moving immersed snapshot enstrophy sum is nonfinite");
				AddFinite(q_sum, weight*item.q_criterion_per_s2, "moving immersed snapshot Q sum is nonfinite");
				if (item.q_criterion_per_s2 > 0.0) AddFinite(q_positive, weight, "moving immersed snapshot Q-positive volume is nonfinite");
				if (item.stagnant) AddFinite(stagnant, weight, "moving immersed snapshot stagnant volume is nonfinite");
				result.points_.push_back(item);
			};
			ForEachUsableVolumePoint(geometry, id, add);
		}
		if (result.points_.size() != count || !(volume > 0.0L)) throw std::runtime_error("moving immersed snapshot volume point coverage is invalid");
		result.metrics_.quadrature_volume_m3 = FiniteCast(volume, "moving immersed snapshot quadrature volume is nonfinite");
		result.metrics_.audited_volume_m3 = geometry.Diagnostics().catalog_estimated_physical_volume_m3;
		RequireFinitePositive(result.metrics_.audited_volume_m3, "moving immersed snapshot audited volume is invalid");
		const double tolerance = 1.0e-12*std::max({1.0, std::abs(result.metrics_.quadrature_volume_m3), std::abs(result.metrics_.audited_volume_m3)});
		if (std::abs(result.metrics_.quadrature_volume_m3-result.metrics_.audited_volume_m3) > tolerance)
			throw std::runtime_error("moving immersed snapshot quadrature volume does not match audited volume");
		result.metrics_.enstrophy_integral_m3_per_s2 = FiniteCast(enstrophy, "moving immersed snapshot enstrophy integral is nonfinite");
		result.metrics_.mean_enstrophy_per_s2 = Divide(result.metrics_.enstrophy_integral_m3_per_s2, result.metrics_.quadrature_volume_m3, "moving immersed snapshot mean enstrophy is nonfinite");
		result.metrics_.mean_q_criterion_per_s2 = Divide(FiniteCast(q_sum, "moving immersed snapshot Q integral is nonfinite"), result.metrics_.quadrature_volume_m3, "moving immersed snapshot mean Q is nonfinite");
		result.metrics_.q_positive_volume_m3 = FiniteCast(q_positive, "moving immersed snapshot Q-positive volume is nonfinite");
		result.metrics_.q_positive_volume_fraction = Divide(result.metrics_.q_positive_volume_m3, result.metrics_.quadrature_volume_m3, "moving immersed snapshot Q-positive fraction is nonfinite");
		result.metrics_.stagnant_volume_m3 = FiniteCast(stagnant, "moving immersed snapshot stagnant volume is nonfinite");
		result.metrics_.stagnant_volume_fraction = Divide(result.metrics_.stagnant_volume_m3, result.metrics_.quadrature_volume_m3, "moving immersed snapshot stagnant fraction is nonfinite");
		SetFlows(result.metrics_, result.request_, result.metrics_.quadrature_volume_m3);
		SetWallMetrics(result.metrics_, geometry, layout, state, result.request_.wall_labels);
		result.content_hash_sha256_ = result.ContentHash();
		result.snapshot_identity_sha256_ = result.SnapshotIdentity(geometry, layout, state);
		return result;
	}

	const std::vector<MovingImmersedFieldPoint>& Points() const noexcept { return points_; }
	const MovingImmersedFlowSnapshotMetrics& Metrics() const noexcept { return metrics_; }
	const MovingImmersedFlowSnapshotOptions& Options() const noexcept { return options_; }
	const MovingImmersedFlowSnapshotRequest& Request() const noexcept { return request_; }
	const std::string& ContentHashSha256() const noexcept { return content_hash_sha256_; }
	const std::string& SnapshotIdentitySha256() const noexcept { return snapshot_identity_sha256_; }

private:
	static bool Usable(const CartesianDomainCell& domain, const CutCellVolumeQuadratureCell& volume)
	{
		return domain.classification == CellClassification::Inside || (domain.classification == CellClassification::Cut && volume.diagnostics.estimated_reference_volume > 0.0);
	}
	template <class Callback> static void ForEachUsableVolumePoint(const MovingCutGeometry& geometry, std::uint64_t id, Callback&& callback)
	{
		if (geometry.Volume().StorageMode() == CutCellVolumeQuadratureStorageMode::Expanded) {
			for (const auto& point : geometry.Volume().UsableRule(geometry.Domain(), id).Points()) callback(point);
		} else ForEachVolumePoint(geometry.Volume().UsableCompactRule(geometry.Domain(), id), std::forward<Callback>(callback));
	}
	static void ValidateIdentity(const MovingCutGeometry& geometry, const ImmersedActiveLayout& layout,
		const ImmersedGlobalFlowState& state, const MovingImmersedFlowSnapshotRequest& request,
		const MovingImmersedFlowSnapshotOptions& options)
	{
		if (!layout.Valid() || !state.Valid() || layout.GeometryIdentity() != geometry.GeometryIdentitySha256()
			|| state.GeometryIdentity() != geometry.GeometryIdentitySha256() || state.NodeIds() != layout.NodeIds()
			|| state.Coefficients().size() != layout.NodeIds().size() || state.PortIds() != layout.PortIds()
			|| state.HasGaugeMultiplier() != layout.HasGaugeRow())
			throw std::invalid_argument("moving immersed snapshot geometry, layout, and state identities do not match");
		if (!std::isfinite(request.time_s) || request.time_s != geometry.Evaluation().EvaluatedTimeS()
			|| request.time_s != state.TimeS() || request.index != state.Index())
			throw std::invalid_argument("moving immersed snapshot time or index does not match state");
		if (!std::isfinite(options.stagnant_speed_threshold_m_per_s) || options.stagnant_speed_threshold_m_per_s < 0.0
			|| !options.maximum_points || !options.maximum_output_bytes)
			throw std::invalid_argument("moving immersed snapshot options are invalid");
		if (!std::isfinite(request.dt_s) || request.dt_s < 0.0 || (request.transition_available && !(request.dt_s > 0.0)) || (!request.transition_available && request.dt_s != 0.0))
			throw std::invalid_argument("moving immersed snapshot transition is invalid");
		if (!std::is_sorted(request.port_labels.begin(), request.port_labels.end())
			|| std::adjacent_find(request.port_labels.begin(), request.port_labels.end()) != request.port_labels.end())
			throw std::invalid_argument("moving immersed snapshot port labels must be sorted and unique");
		if ((!request.port_flows_available && !request.port_flows.empty())
			|| (request.port_flows_available && request.port_flows.size() != request.port_labels.size()))
			throw std::invalid_argument("moving immersed snapshot port flow availability is invalid");
		for (std::size_t i = 0; i < request.port_flows.size(); ++i) {
			if (request.port_flows[i].label != request.port_labels[i] || !std::isfinite(request.port_flows[i].outward_flow_m3_s)
				|| (i && request.port_flows[i-1].label >= request.port_flows[i].label))
				throw std::invalid_argument("moving immersed snapshot port labels are invalid");
		}
		if (!std::is_sorted(request.wall_labels.begin(), request.wall_labels.end())
			|| std::adjacent_find(request.wall_labels.begin(), request.wall_labels.end()) != request.wall_labels.end())
			throw std::invalid_argument("moving immersed snapshot wall labels must be sorted and unique");
		ValidateSurfacePartition(geometry, request.wall_labels, request.port_labels);
	}
	static void ValidateSurfacePartition(const MovingCutGeometry& geometry, const std::vector<std::uint32_t>& walls,
		const std::vector<std::uint32_t>& ports)
	{
		std::vector<std::uint32_t> all;
		std::vector<long double> areas;
		for (std::uint64_t id = 0; id < geometry.Surface().Cells().size(); ++id) {
			const auto& rule = geometry.Surface().UsableRule(geometry.Domain(), id);
			for (const auto& point : rule.Points()) {
				if (point.boundary_id <= 0 || static_cast<std::uint64_t>(point.boundary_id) > std::numeric_limits<std::uint32_t>::max())
					throw std::invalid_argument("moving immersed snapshot surface boundary label is invalid");
				const auto boundary_id = static_cast<std::uint32_t>(point.boundary_id);
				RequireFinitePositive(point.weight, "moving immersed snapshot wall weight is invalid");
				const auto found = std::lower_bound(all.begin(), all.end(), boundary_id);
				const auto offset = static_cast<std::size_t>(found-all.begin());
				if (found == all.end() || *found != boundary_id) { all.insert(found, boundary_id); areas.insert(areas.begin()+static_cast<std::ptrdiff_t>(offset), 0.0L); }
				AddFinite(areas[offset], point.weight, "moving immersed snapshot surface area is nonfinite");
			}
		}
		std::vector<std::uint32_t> requested=walls; requested.insert(requested.end(),ports.begin(),ports.end()); std::sort(requested.begin(),requested.end());
		if (std::adjacent_find(requested.begin(),requested.end()) != requested.end() || requested != all)
			throw std::invalid_argument("moving immersed snapshot wall and port labels must partition surface boundary labels");
		for (const auto area : areas) if (!(FiniteCast(area, "moving immersed snapshot surface area is nonfinite") > 0.0))
			throw std::invalid_argument("moving immersed snapshot requested surface label has no finite positive surface area");
	}
	static std::size_t CountPoints(const MovingCutGeometry& geometry, const ImmersedActiveLayout& layout, const MovingImmersedFlowSnapshotOptions& options)
	{
		std::size_t count = 0;
		for (std::uint64_t id = 0; id < geometry.Domain().Cells().size(); ++id) {
			const auto& cell = geometry.Domain().Cells()[static_cast<std::size_t>(id)]; const auto& rule = geometry.Volume().Cell(id);
			if (!Usable(cell, rule)) continue;
			const std::size_t add = geometry.Volume().StorageMode() == CutCellVolumeQuadratureStorageMode::Expanded
				? geometry.Volume().UsableRule(geometry.Domain(), id).Points().size() : CompactCutCellVolumeLogicalPointCount(geometry.Volume().UsableCompactRule(geometry.Domain(), id));
			if (add > options.maximum_points-count) throw std::overflow_error("moving immersed snapshot point cap exceeded");
			count += add;
		}
		if (!count) throw std::runtime_error("moving immersed snapshot volume is empty");
		if (count > std::numeric_limits<std::size_t>::max()/sizeof(MovingImmersedFieldPoint) || count*sizeof(MovingImmersedFieldPoint) > options.maximum_output_bytes)
			throw std::overflow_error("moving immersed snapshot output byte cap exceeded");
		(void)layout; return count;
	}
	static double Dot(const std::array<double, 3>& a, const std::array<double, 3>& b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
	static double Norm(const std::array<double, 3>& value) { const double result=std::sqrt(Dot(value,value)); if (!std::isfinite(result)) throw std::runtime_error("moving immersed snapshot vector norm is nonfinite"); return result; }
	static void RequireFinitePositive(double value, const char* message) { if (!std::isfinite(value) || !(value > 0.0)) throw std::runtime_error(message); }
	static void ValidatePoint(const MovingImmersedFieldPoint& point)
	{
		for (double value : point.physical_m) if (!std::isfinite(value)) throw std::runtime_error("moving immersed snapshot physical point is nonfinite");
		for (double value : point.parametric) if (!std::isfinite(value)) throw std::runtime_error("moving immersed snapshot parametric point is nonfinite");
		for (double value : point.velocity_m_per_s) if (!std::isfinite(value)) throw std::runtime_error("moving immersed snapshot velocity is nonfinite");
		for (double value : point.vorticity_per_s) if (!std::isfinite(value)) throw std::runtime_error("moving immersed snapshot vorticity is nonfinite");
		for (double value : {point.pressure, point.speed_m_per_s, point.q_criterion_per_s2, point.enstrophy_density_per_s2, point.physical_integration_weight_m3}) if (!std::isfinite(value)) throw std::runtime_error("moving immersed snapshot derived field is nonfinite");
	}
	static void AddFinite(long double& sum, double value, const char* message) { if (!std::isfinite(value) || !std::isfinite(sum += static_cast<long double>(value))) throw std::runtime_error(message); }
	static double FiniteCast(long double value, const char* message) { const double result=static_cast<double>(value); if (!std::isfinite(result)) throw std::runtime_error(message); return result; }
	static double Divide(double top, double bottom, const char* message) { const double value=top/bottom; if (!std::isfinite(value)) throw std::runtime_error(message); return value; }
	static void SetFlows(MovingImmersedFlowSnapshotMetrics& metrics, const MovingImmersedFlowSnapshotRequest& request, double volume)
	{
		if (!request.port_flows_available) return;
		metrics.endpoint_turnover_available = true;
		for (const auto& port : request.port_flows) { metrics.inlet_flow_m3_s += std::max(0.0, -port.outward_flow_m3_s); metrics.outlet_flow_m3_s += std::max(0.0, port.outward_flow_m3_s); }
		if (!std::isfinite(metrics.inlet_flow_m3_s) || !std::isfinite(metrics.outlet_flow_m3_s)) throw std::runtime_error("moving immersed snapshot endpoint flow is nonfinite");
		metrics.endpoint_turnover_rate_per_s = Divide(metrics.inlet_flow_m3_s, volume, "moving immersed snapshot turnover rate is nonfinite");
		if (metrics.inlet_flow_m3_s > 0.0) metrics.endpoint_turnover_time_s = Divide(volume, metrics.inlet_flow_m3_s, "moving immersed snapshot turnover time is nonfinite");
		if (request.transition_available) { metrics.well_mixed_replacement_fraction_over_step = -std::expm1(-metrics.inlet_flow_m3_s*request.dt_s/volume); if (!std::isfinite(metrics.well_mixed_replacement_fraction_over_step)) throw std::runtime_error("moving immersed snapshot replacement fraction is nonfinite"); }
	}
	static void SetWallMetrics(MovingImmersedFlowSnapshotMetrics& metrics, const MovingCutGeometry& geometry,
		const ImmersedActiveLayout& layout, const ImmersedGlobalFlowState& state, const std::vector<std::uint32_t>& labels)
	{
		if (labels.empty()) return;
		long double l2 = 0.0L, area = 0.0L; double maximum = 0.0;
		for (std::uint64_t id = 0; id < geometry.Surface().Cells().size(); ++id) {
			const auto element = geometry.Domain().Background().MaterializeElement(id);
			const auto& rule = geometry.Surface().UsableRule(geometry.Domain(), id); const auto& provenance = geometry.Surface().UsableProvenance(geometry.Domain(), id);
			for (std::size_t q = 0; q < rule.Points().size(); ++q) {
				const auto& point = rule.Points()[q]; if (!std::binary_search(labels.begin(), labels.end(), static_cast<std::uint32_t>(point.boundary_id))) continue;
				const auto basis = EvaluateBasis(element, point.parametric[0], point.parametric[1], point.parametric[2], false); std::array<double,3> u{};
				for (std::size_t a=0;a<element.connectivity.size();++a) { const auto& c=state.Coefficients()[layout.LocalNode(element.connectivity[a])]; for(std::size_t d=0;d<3;++d)u[d]+=basis.value[a]*c[d]; }
				const auto w=geometry.Evaluation().WallVelocity(provenance[q].canonical_triangle, provenance[q].canonical_barycentric); std::array<double,3> difference{{u[0]-w[0],u[1]-w[1],u[2]-w[2]}}; const double norm=Norm(difference);
				RequireFinitePositive(point.weight,"moving immersed snapshot wall weight is invalid"); AddFinite(area,point.weight,"moving immersed snapshot wall area is nonfinite"); AddFinite(l2,point.weight*norm*norm,"moving immersed snapshot wall L2 is nonfinite"); maximum=std::max(maximum,norm);
			}
		}
		metrics.wall_area_m2=FiniteCast(area,"moving immersed snapshot wall area is nonfinite"); RequireFinitePositive(metrics.wall_area_m2,"moving immersed snapshot configured wall area is invalid"); metrics.wall_relative_velocity_squared_area_integral_m4_per_s2=FiniteCast(l2,"moving immersed snapshot wall squared-area integral is nonfinite"); metrics.wall_relative_velocity_rms_m_per_s=std::sqrt(Divide(metrics.wall_relative_velocity_squared_area_integral_m4_per_s2,metrics.wall_area_m2,"moving immersed snapshot wall RMS is nonfinite")); metrics.wall_relative_velocity_max_m_per_s=maximum;
	}
	std::string ContentHash() const
	{
		Sha256 hash; immersed_transient_detail::AppendString(hash,"MovingImmersedFlowSnapshot/content/v3"); hash.AppendNormalizedDouble(options_.stagnant_speed_threshold_m_per_s); hash.AppendLittleEndian64(options_.maximum_points); hash.AppendLittleEndian64(options_.maximum_output_bytes); hash.AppendLittleEndian32(static_cast<std::uint32_t>(storage_mode_)); hash.AppendLittleEndian32(request_.transition_available?1u:0u); hash.AppendNormalizedDouble(request_.dt_s); hash.AppendLittleEndian32(request_.port_flows_available?1u:0u);
		hash.AppendLittleEndian64(request_.port_labels.size()); for(auto label:request_.port_labels)hash.AppendLittleEndian32(label); hash.AppendLittleEndian64(request_.port_flows.size()); for(const auto& port:request_.port_flows){hash.AppendLittleEndian64(port.label);hash.AppendNormalizedDouble(port.outward_flow_m3_s);} hash.AppendLittleEndian64(request_.wall_labels.size());for(auto label:request_.wall_labels)hash.AppendLittleEndian32(label);
		for(const auto& point:points_){for(double x:point.physical_m)hash.AppendNormalizedDouble(x);for(double x:point.parametric)hash.AppendNormalizedDouble(x);for(double x:point.velocity_m_per_s)hash.AppendNormalizedDouble(x);for(double x:point.vorticity_per_s)hash.AppendNormalizedDouble(x);for(double x:{point.pressure,point.speed_m_per_s,point.q_criterion_per_s2,point.enstrophy_density_per_s2,point.physical_integration_weight_m3})hash.AppendNormalizedDouble(x);hash.AppendLittleEndian64(point.background_cell_id);hash.AppendLittleEndian64(point.cell_quadrature_ordinal);hash.AppendLittleEndian32(static_cast<std::uint32_t>(point.cell_kind));hash.AppendLittleEndian32(point.stagnant?1u:0u);}
		const auto& m=metrics_;for(double x:{m.quadrature_volume_m3,m.audited_volume_m3,m.enstrophy_integral_m3_per_s2,m.mean_enstrophy_per_s2,m.mean_q_criterion_per_s2,m.q_positive_volume_m3,m.q_positive_volume_fraction,m.stagnant_volume_m3,m.stagnant_volume_fraction,m.inlet_flow_m3_s,m.outlet_flow_m3_s,m.endpoint_turnover_rate_per_s,m.well_mixed_replacement_fraction_over_step,m.wall_relative_velocity_squared_area_integral_m4_per_s2,m.wall_relative_velocity_rms_m_per_s,m.wall_relative_velocity_max_m_per_s,m.wall_area_m2})hash.AppendNormalizedDouble(x);hash.AppendLittleEndian32(m.endpoint_turnover_available?1u:0u);hash.AppendLittleEndian32(m.endpoint_turnover_time_s?1u:0u);if(m.endpoint_turnover_time_s)hash.AppendNormalizedDouble(*m.endpoint_turnover_time_s);return hash.Hex();
	}
	std::string SnapshotIdentity(const MovingCutGeometry& geometry, const ImmersedActiveLayout& layout, const ImmersedGlobalFlowState& state) const
	{ Sha256 hash; immersed_transient_detail::AppendString(hash,"MovingImmersedFlowSnapshot/identity/v1"); immersed_transient_detail::AppendString(hash,content_hash_sha256_); immersed_transient_detail::AppendString(hash,geometry.GeometryIdentitySha256()); immersed_transient_detail::AppendString(hash,geometry.PublicationIdentitySha256()); immersed_transient_detail::AppendString(hash,layout.HashSha256()); immersed_transient_detail::AppendString(hash,state.HashSha256()); hash.AppendNormalizedDouble(request_.time_s); hash.AppendLittleEndian64(request_.index); return hash.Hex(); }
	MovingImmersedFlowSnapshotOptions options_{}; MovingImmersedFlowSnapshotRequest request_{}; CutCellVolumeQuadratureStorageMode storage_mode_ = CutCellVolumeQuadratureStorageMode::Expanded; std::vector<MovingImmersedFieldPoint> points_; MovingImmersedFlowSnapshotMetrics metrics_{}; std::string content_hash_sha256_, snapshot_identity_sha256_;
};

} // namespace iga

#endif
