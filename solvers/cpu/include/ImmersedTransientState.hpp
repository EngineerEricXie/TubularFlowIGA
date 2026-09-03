#ifndef IGA_IMMERSED_TRANSIENT_STATE_HPP
#define IGA_IMMERSED_TRANSIENT_STATE_HPP

// Value-only, global-ID keyed transient data.  This header deliberately owns
// no mutable geometry or solver state; PR7.3a uses it only on fixed geometry.
#include "CutCellVolumeQuadrature.hpp"
#include "Sha256.hpp"

#include <petscsys.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace iga {

namespace immersed_transient_detail {
inline void AppendString(Sha256& hash, const std::string& value)
{
	hash.AppendLittleEndian64(static_cast<std::uint64_t>(value.size()));
	hash.Append(value.data(), value.size());
}
inline void AppendIds(Sha256& hash, const std::vector<std::int32_t>& values)
{
	hash.AppendLittleEndian64(static_cast<std::uint64_t>(values.size()));
	for (const auto value : values) hash.AppendLittleEndian32(static_cast<std::uint32_t>(value));
}
inline void RequireFinite(double value, const char* message)
{
	if (!std::isfinite(value)) throw std::invalid_argument(message);
}
inline void RequireSortedUnique(const std::vector<std::int32_t>& values, const char* message)
{
	for (std::size_t i = 0; i < values.size(); ++i) {
		if (values[i] < 0 || (i && values[i-1] >= values[i])) throw std::invalid_argument(message);
	}
}
inline bool SameGridSpec(const CubicCartesianGridSpec& left, const CubicCartesianGridSpec& right)
{
	return left.lower_m == right.lower_m && left.upper_m == right.upper_m && left.cells == right.cells;
}
inline std::size_t CheckedMultiply(std::size_t left, std::size_t right, const char* message)
{
	if (left != 0 && right > std::numeric_limits<std::size_t>::max()/left)
		throw std::overflow_error(message);
	return left*right;
}
inline std::size_t CheckedAdd(std::size_t left, std::size_t right, const char* message)
{
	if (right > std::numeric_limits<std::size_t>::max()-left) throw std::overflow_error(message);
	return left+right;
}
inline void RequirePetscIntRange(std::size_t value, const char* message)
{
	if (value > static_cast<std::size_t>(std::numeric_limits<PetscInt>::max()))
		throw std::overflow_error(message);
}
} // namespace immersed_transient_detail

class ImmersedActiveLayout {
public:
	ImmersedActiveLayout() = default;
	ImmersedActiveLayout(const ImmersedActiveLayout&) = default;
	ImmersedActiveLayout& operator=(const ImmersedActiveLayout&) = default;
	ImmersedActiveLayout(ImmersedActiveLayout&&) noexcept = default;
	ImmersedActiveLayout& operator=(ImmersedActiveLayout&&) noexcept = default;

	static ImmersedActiveLayout Build(const CartesianDomainClassification& domain,
		const CutCellVolumeQuadratureCatalog& volume, std::string geometry_identity,
		std::vector<std::uint64_t> port_ids = {}, bool has_gauge_row = false)
	{
		if (geometry_identity.empty()) throw std::invalid_argument("immersed active layout geometry identity is empty");
		const auto& cells = domain.Cells();
		if (cells.size() != volume.Cells().size() || cells.size() != domain.Background().ElementCount())
			throw std::invalid_argument("immersed active layout needs complete matching domain and volume catalogs");
		if (!immersed_transient_detail::SameGridSpec(domain.Background().Spec(), volume.GridSpec())
			|| domain.SurfaceCanonicalHash() != volume.SurfaceCanonicalHash())
			throw std::invalid_argument("immersed active layout domain and volume bindings differ");
		ImmersedActiveLayout result;
		result.geometry_identity_ = std::move(geometry_identity);
		result.port_ids_ = std::move(port_ids);
		if (!std::is_sorted(result.port_ids_.begin(), result.port_ids_.end())
			|| std::adjacent_find(result.port_ids_.begin(), result.port_ids_.end()) != result.port_ids_.end())
			throw std::invalid_argument("immersed active layout port IDs must be sorted and unique");
		result.has_gauge_row_ = has_gauge_row;
		for (std::size_t id = 0; id < cells.size(); ++id) {
			const auto& cell = cells[id]; const auto& quadrature = volume.Cells()[id];
			if (cell.id != id || quadrature.id != id || cell.classification != quadrature.classification
				|| cell.ambiguous || !quadrature.usable)
				throw std::invalid_argument("immersed active layout catalog ordering differs");
			const double fraction = quadrature.diagnostics.estimated_reference_volume;
			if (!std::isfinite(fraction) || fraction < 0.0 || fraction > 1.0)
				throw std::invalid_argument("immersed active layout cut fraction is invalid");
			const bool usable = cell.classification == CellClassification::Inside
				|| (cell.classification == CellClassification::Cut && quadrature.usable && fraction > 0.0);
			if (!usable) continue;
			const auto element = domain.Background().MaterializeElement(static_cast<std::uint64_t>(id));
			if (element.connectivity.size() != 64) throw std::runtime_error("immersed active layout requires cubic 64-node elements");
			result.node_ids_.insert(result.node_ids_.end(), element.connectivity.begin(), element.connectivity.end());
		}
		std::sort(result.node_ids_.begin(), result.node_ids_.end());
		result.node_ids_.erase(std::unique(result.node_ids_.begin(), result.node_ids_.end()), result.node_ids_.end());
		immersed_transient_detail::RequireSortedUnique(result.node_ids_, "immersed active layout node IDs are invalid");
		const auto node_rows = immersed_transient_detail::CheckedMultiply(result.node_ids_.size(), 4u,
			"immersed active layout node row count overflows");
		immersed_transient_detail::RequirePetscIntRange(node_rows,
			"immersed active layout node row count exceeds PetscInt");
		for (const auto node_id : result.node_ids_)
			if (static_cast<std::uint64_t>(node_id) > static_cast<std::uint64_t>(std::numeric_limits<PetscInt>::max()))
				throw std::overflow_error("immersed active layout global node ID exceeds PetscInt");
		const auto rows = immersed_transient_detail::CheckedAdd(
			immersed_transient_detail::CheckedAdd(node_rows, result.port_ids_.size(),
				"immersed active layout controller row count overflows"), has_gauge_row ? 1u : 0u,
			"immersed active layout gauge row count overflows");
		immersed_transient_detail::RequirePetscIntRange(rows, "immersed active layout row count exceeds PetscInt");
		result.global_to_local_.reserve(result.node_ids_.size());
		for (std::size_t i = 0; i < result.node_ids_.size(); ++i)
			result.global_to_local_.emplace(result.node_ids_[i], i);
		result.hash_ = result.ComputeHash();
		return result;
	}

	const std::string& GeometryIdentity() const noexcept { return geometry_identity_; }
	bool Valid() const noexcept { return !geometry_identity_.empty() && global_to_local_.size() == node_ids_.size(); }
	const std::vector<std::int32_t>& NodeIds() const noexcept { return node_ids_; }
	const std::vector<std::uint64_t>& PortIds() const noexcept { return port_ids_; }
	bool HasGaugeRow() const noexcept { return has_gauge_row_; }
	std::size_t NodeFieldRows() const { return immersed_transient_detail::CheckedMultiply(node_ids_.size(), 4u, "immersed active layout node row count overflows"); }
	std::size_t ControllerRow(std::uint64_t port_id) const
	{
		const auto pos = std::lower_bound(port_ids_.begin(), port_ids_.end(), port_id);
		if (pos == port_ids_.end() || *pos != port_id) throw std::out_of_range("immersed active layout port ID is absent");
		return NodeFieldRows()+static_cast<std::size_t>(pos-port_ids_.begin());
	}
	std::size_t GaugeRow() const
	{
		if (!has_gauge_row_) throw std::logic_error("immersed active layout has no gauge row");
		return NodeFieldRows()+port_ids_.size();
	}
	std::size_t Rows() const { return immersed_transient_detail::CheckedAdd(immersed_transient_detail::CheckedAdd(NodeFieldRows(), port_ids_.size(), "immersed active layout controller row count overflows"), has_gauge_row_ ? 1u : 0u, "immersed active layout gauge row count overflows"); }
	std::size_t LocalNode(std::int32_t global_id) const
	{
		const auto found = global_to_local_.find(global_id);
		if (found == global_to_local_.end()) throw std::out_of_range("immersed active layout global node ID is absent");
		return found->second;
	}
	const std::string& HashSha256() const noexcept { return hash_; }

private:
	std::string ComputeHash() const
	{
		Sha256 hash; immersed_transient_detail::AppendString(hash, "ImmersedActiveLayout/v1");
		immersed_transient_detail::AppendString(hash, geometry_identity_); immersed_transient_detail::AppendIds(hash, node_ids_);
		hash.AppendLittleEndian64(static_cast<std::uint64_t>(port_ids_.size()));
		for (const auto id : port_ids_) hash.AppendLittleEndian64(id);
		hash.AppendLittleEndian32(has_gauge_row_ ? 1u : 0u); return hash.Hex();
	}
	std::string geometry_identity_, hash_; std::vector<std::int32_t> node_ids_; std::vector<std::uint64_t> port_ids_;
	bool has_gauge_row_ = false; std::unordered_map<std::int32_t, std::size_t> global_to_local_;
};

class ImmersedGlobalFlowState {
public:
	ImmersedGlobalFlowState(double time_s, std::uint64_t index, const ImmersedActiveLayout& layout,
		std::vector<std::array<double, 4>> coefficients, std::vector<double> port_multipliers = {},
		bool has_gauge_multiplier = false, double gauge_multiplier = 0.0)
		: time_s_(time_s), index_(index), geometry_identity_(layout.GeometryIdentity()), node_ids_(layout.NodeIds()),
		coefficients_(std::move(coefficients)), port_ids_(layout.PortIds()), port_multipliers_(std::move(port_multipliers)),
		has_gauge_multiplier_(has_gauge_multiplier), gauge_multiplier_(gauge_multiplier)
	{
		if (!std::isfinite(time_s_) || time_s_ < 0.0) throw std::invalid_argument("immersed global flow state time is invalid");
		if (coefficients_.size() != node_ids_.size() || port_multipliers_.size() != port_ids_.size()
			|| has_gauge_multiplier_ != layout.HasGaugeRow()) throw std::invalid_argument("immersed global flow state does not match layout");
		for (const auto& q : coefficients_) for (double value : q) immersed_transient_detail::RequireFinite(value, "immersed global flow state coefficient is not finite");
		for (double value : port_multipliers_) immersed_transient_detail::RequireFinite(value, "immersed global flow state port multiplier is not finite");
		immersed_transient_detail::RequireFinite(gauge_multiplier_, "immersed global flow state gauge multiplier is not finite");
		if (!has_gauge_multiplier_ && gauge_multiplier_ != 0.0)
			throw std::invalid_argument("immersed global flow state without a gauge must use canonical zero");
		if (!has_gauge_multiplier_) gauge_multiplier_ = 0.0;
		hash_ = ComputeHash();
	}
	const std::string& GeometryIdentity() const noexcept { return geometry_identity_; }
	bool Valid() const noexcept { return std::isfinite(time_s_) && !geometry_identity_.empty() && coefficients_.size() == node_ids_.size(); }
	double TimeS() const noexcept { return time_s_; } std::uint64_t Index() const noexcept { return index_; }
	const std::vector<std::int32_t>& NodeIds() const noexcept { return node_ids_; }
	const std::vector<std::array<double, 4>>& Coefficients() const noexcept { return coefficients_; }
	const std::vector<std::uint64_t>& PortIds() const noexcept { return port_ids_; }
	const std::vector<double>& PortMultipliers() const noexcept { return port_multipliers_; }
	bool HasGaugeMultiplier() const noexcept { return has_gauge_multiplier_; } double GaugeMultiplier() const noexcept { return gauge_multiplier_; }
	const std::string& HashSha256() const noexcept { return hash_; }
private:
	std::string ComputeHash() const { Sha256 h; immersed_transient_detail::AppendString(h, "ImmersedGlobalFlowState/v1"); h.AppendNormalizedDouble(time_s_); h.AppendLittleEndian64(index_); immersed_transient_detail::AppendString(h, geometry_identity_); immersed_transient_detail::AppendIds(h,node_ids_); for(const auto&q:coefficients_)for(double x:q)h.AppendNormalizedDouble(x); h.AppendLittleEndian64(port_ids_.size()); for(std::size_t i=0;i<port_ids_.size();++i){h.AppendLittleEndian64(port_ids_[i]);h.AppendNormalizedDouble(port_multipliers_[i]);} h.AppendLittleEndian32(has_gauge_multiplier_?1u:0u); if(has_gauge_multiplier_)h.AppendNormalizedDouble(gauge_multiplier_); return h.Hex(); }
	double time_s_; std::uint64_t index_; std::string geometry_identity_,hash_; std::vector<std::int32_t> node_ids_; std::vector<std::array<double,4>> coefficients_; std::vector<std::uint64_t> port_ids_; std::vector<double> port_multipliers_; bool has_gauge_multiplier_; double gauge_multiplier_;
};

enum class ImmersedVelocityHistoryProvenance : std::uint8_t { Committed = 0, Extended = 1 };
class ImmersedVelocityHistory {
public:
	ImmersedVelocityHistory(double source_time_s, double target_time_s, std::string source_geometry_identity,
		std::string target_geometry_identity, std::vector<std::int32_t> node_ids,
		std::vector<std::array<double, 3>> velocities, std::vector<ImmersedVelocityHistoryProvenance> provenance)
		: source_time_s_(source_time_s), target_time_s_(target_time_s), source_geometry_identity_(std::move(source_geometry_identity)), target_geometry_identity_(std::move(target_geometry_identity)), node_ids_(std::move(node_ids)), velocities_(std::move(velocities)), provenance_(std::move(provenance))
	{
		if (!std::isfinite(source_time_s_) || source_time_s_ < 0.0 || !std::isfinite(target_time_s_) || target_time_s_ < 0.0 || source_geometry_identity_.empty() || target_geometry_identity_.empty()) throw std::invalid_argument("immersed velocity history identity or time is invalid");
		immersed_transient_detail::RequireSortedUnique(node_ids_, "immersed velocity history node IDs are not sorted and unique");
		if (velocities_.size()!=node_ids_.size() || provenance_.size()!=node_ids_.size()) throw std::invalid_argument("immersed velocity history coverage is inconsistent");
		for(const auto&q:velocities_)for(double x:q) immersed_transient_detail::RequireFinite(x,"immersed velocity history velocity is not finite");
		for(const auto p:provenance_) if(p!=ImmersedVelocityHistoryProvenance::Committed&&p!=ImmersedVelocityHistoryProvenance::Extended) throw std::invalid_argument("immersed velocity history provenance is invalid");
		hash_=ComputeHash();
	}
	void ValidateCoverage(const ImmersedActiveLayout& layout) const { if(target_geometry_identity_!=layout.GeometryIdentity() || node_ids_!=layout.NodeIds()) throw std::invalid_argument("immersed velocity history does not exactly cover target layout"); }
	bool Valid() const noexcept{return std::isfinite(source_time_s_)&&std::isfinite(target_time_s_)&&!source_geometry_identity_.empty()&&!target_geometry_identity_.empty()&&node_ids_.size()==velocities_.size()&&node_ids_.size()==provenance_.size();} double SourceTimeS()const noexcept{return source_time_s_;} double TargetTimeS()const noexcept{return target_time_s_;} const std::vector<std::int32_t>& NodeIds() const noexcept{return node_ids_;} const std::vector<std::array<double,3>>& Velocities()const noexcept{return velocities_;} const std::vector<ImmersedVelocityHistoryProvenance>& Provenance()const noexcept{return provenance_;} const std::string& SourceGeometryIdentity()const noexcept{return source_geometry_identity_;} const std::string& TargetGeometryIdentity()const noexcept{return target_geometry_identity_;} const std::string& HashSha256()const noexcept{return hash_;}
private:
	std::string ComputeHash()const{Sha256 h;immersed_transient_detail::AppendString(h,"ImmersedVelocityHistory/v1");h.AppendNormalizedDouble(source_time_s_);h.AppendNormalizedDouble(target_time_s_);immersed_transient_detail::AppendString(h,source_geometry_identity_);immersed_transient_detail::AppendString(h,target_geometry_identity_);immersed_transient_detail::AppendIds(h,node_ids_);for(std::size_t i=0;i<node_ids_.size();++i){for(double x:velocities_[i])h.AppendNormalizedDouble(x);h.AppendLittleEndian32(static_cast<std::uint32_t>(provenance_[i]));}return h.Hex();}
	double source_time_s_,target_time_s_;std::string source_geometry_identity_,target_geometry_identity_,hash_;std::vector<std::int32_t>node_ids_;std::vector<std::array<double,3>>velocities_;std::vector<ImmersedVelocityHistoryProvenance>provenance_;
};

// The element's connectivity is global.  Localizing by it here prevents a
// same-sized positional vector from being silently applied to another cell.
inline std::vector<std::array<double, 3>> LocalizeImmersedVelocityHistory(
	const Element& element, const ImmersedActiveLayout& layout,
	const ImmersedVelocityHistory& history, double target_time_s)
{
	if (!layout.Valid() || !history.Valid() || !std::isfinite(target_time_s) || target_time_s < 0.0)
		throw std::invalid_argument("immersed velocity history localization identity or time is invalid");
	history.ValidateCoverage(layout);
	if (history.TargetTimeS() != target_time_s)
		throw std::invalid_argument("immersed velocity history target time does not match assembly time");
	if (element.connectivity.empty()) throw std::invalid_argument("immersed velocity history cannot localize an empty element");
	std::vector<std::array<double, 3>> local;
	local.reserve(element.connectivity.size());
	for (const auto global_id : element.connectivity) {
		if (global_id < 0) throw std::out_of_range("immersed velocity history element global node ID is invalid");
		const auto index = layout.LocalNode(global_id);
		if (index >= history.Velocities().size() || index >= history.Provenance().size())
			throw std::out_of_range("immersed velocity history local node index is out of range");
		if (history.Provenance()[index] != ImmersedVelocityHistoryProvenance::Committed
			&& history.Provenance()[index] != ImmersedVelocityHistoryProvenance::Extended)
			throw std::invalid_argument("immersed velocity history provenance is invalid");
		for (const auto value : history.Velocities()[index])
			immersed_transient_detail::RequireFinite(value, "immersed velocity history velocity is not finite");
		local.push_back(history.Velocities()[index]);
	}
	return local;
}

inline ImmersedVelocityHistory BuildIdentityImmersedVelocityHistory(const ImmersedGlobalFlowState& source,
	const ImmersedActiveLayout& target, double target_time_s)
{
	if (source.GeometryIdentity()!=target.GeometryIdentity() || source.NodeIds()!=target.NodeIds()) throw std::invalid_argument("identity immersed velocity history requires identical fixed geometry and layout coverage");
	if (!std::isfinite(target_time_s) || target_time_s<0.0) throw std::invalid_argument("identity immersed velocity history target time is invalid");
	std::vector<std::array<double,3>> velocities(source.Coefficients().size());
	for(std::size_t i=0;i<velocities.size();++i) for(int c=0;c<3;++c) velocities[i][c]=source.Coefficients()[i][c];
	return ImmersedVelocityHistory(source.TimeS(),target_time_s,source.GeometryIdentity(),target.GeometryIdentity(),source.NodeIds(),std::move(velocities),std::vector<ImmersedVelocityHistoryProvenance>(source.NodeIds().size(),ImmersedVelocityHistoryProvenance::Committed));
}

} // namespace iga

#endif
