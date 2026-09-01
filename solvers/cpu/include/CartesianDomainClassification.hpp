#ifndef IGA_CARTESIAN_DOMAIN_CLASSIFICATION_HPP
#define IGA_CARTESIAN_DOMAIN_CLASSIFICATION_HPP

#include "CartesianBackground.hpp"
#include "SurfaceSpatialIndex.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

enum class CellClassification { Outside, Inside, Cut };
struct CartesianDomainCell {
	std::uint64_t id = 0;
	std::array<std::uint32_t, 3> index{};
	SurfaceAabb bounds;
	CellClassification classification = CellClassification::Outside;
	std::vector<std::size_t> triangle_ids;
	bool boundary_only_contact = false;
	bool ambiguous = false;
};
struct CartesianDomainDiagnostics {
	std::size_t outside_count = 0, inside_count = 0, cut_count = 0, ambiguous_count = 0;
	std::string surface_canonical_hash;
};

class CartesianDomainClassification {
public:
	CartesianDomainClassification(const CartesianDomainClassification&) = delete;
	CartesianDomainClassification& operator=(const CartesianDomainClassification&) = delete;
	CartesianDomainClassification(CubicCartesianBackground background, SurfaceSpatialIndex surface)
		: background_(std::move(background)), surface_(std::move(surface))
	{
		const auto& surface_bounds = surface_.Surface().Bounds(); const auto& spec = background_.Spec();
		try { for (std::size_t axis = 0; axis < 3; ++axis) if (exact_dyadic::CoordinateDifferenceSign(surface_bounds.minimum[axis], spec.lower_m[axis]) < 0 || exact_dyadic::CoordinateDifferenceSign(surface_bounds.maximum[axis], spec.upper_m[axis]) > 0) throw std::invalid_argument("Cartesian background does not contain surface bounds"); }
		catch (const std::overflow_error&) {
#ifndef IGA_EXACT_DYADIC_TESTING
			throw;
#endif
		}
		if (background_.ElementCount() > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) throw std::overflow_error("Cartesian domain cell count exceeds size range");
		cells_.reserve(static_cast<std::size_t>(background_.ElementCount())); diagnostics_.surface_canonical_hash = surface_.Surface().CanonicalSha256();
		for (std::uint64_t id = 0; id < background_.ElementCount(); ++id) Classify(id);
	}
	CartesianDomainClassification(CartesianDomainClassification&&) noexcept = default;
	CartesianDomainClassification& operator=(CartesianDomainClassification&&) noexcept = default;
	const CubicCartesianBackground& Background() const { return background_; }
	const SurfaceSpatialIndex& SurfaceIndex() const { return surface_; }
	const std::vector<CartesianDomainCell>& Cells() const { return cells_; }
	const std::vector<std::uint64_t>& ActiveCellIds() const { return active_ids_; }
	const CartesianDomainDiagnostics& Diagnostics() const { return diagnostics_; }
	const std::string& SurfaceCanonicalHash() const { return diagnostics_.surface_canonical_hash; }

private:
	void Classify(std::uint64_t id)
	{
		const auto cell = background_.Cell(id); CartesianDomainCell record; record.id = id; record.index = cell.index; record.bounds.minimum = cell.lower_m; record.bounds.maximum = cell.upper_m;
		const auto contact = surface_.IntersectBox(record.bounds, &record.triangle_ids);
		if (contact == BoxContact::Interior || contact == BoxContact::BoundaryOnly || contact == BoxContact::Ambiguous) {
			record.classification = CellClassification::Cut; record.boundary_only_contact = contact == BoxContact::BoundaryOnly; record.ambiguous = contact == BoxContact::Ambiguous;
			++diagnostics_.cut_count; if (record.ambiguous) ++diagnostics_.ambiguous_count; active_ids_.push_back(id);
		} else {
			std::array<double, 3> center{}; for (std::size_t axis = 0; axis < 3; ++axis) center[axis] = static_cast<double>((static_cast<long double>(cell.lower_m[axis])+cell.upper_m[axis])/2.0L);
			const auto location = surface_.LocatePoint(center);
			if (location == PointLocation::Inside) { record.classification = CellClassification::Inside; ++diagnostics_.inside_count; active_ids_.push_back(id); }
			else if (location == PointLocation::Outside) { record.classification = CellClassification::Outside; ++diagnostics_.outside_count; }
			else { record.classification = CellClassification::Cut; record.ambiguous = location == PointLocation::Ambiguous; ++diagnostics_.cut_count; if (record.ambiguous) ++diagnostics_.ambiguous_count; active_ids_.push_back(id); }
		}
		cells_.push_back(std::move(record));
	}

	CubicCartesianBackground background_;
	SurfaceSpatialIndex surface_;
	std::vector<CartesianDomainCell> cells_;
	std::vector<std::uint64_t> active_ids_;
	CartesianDomainDiagnostics diagnostics_;
};

} // namespace iga

#endif
