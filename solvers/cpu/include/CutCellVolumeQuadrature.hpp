#ifndef IGA_CUT_CELL_VOLUME_QUADRATURE_HPP
#define IGA_CUT_CELL_VOLUME_QUADRATURE_HPP

#include "CartesianDomainClassification.hpp"
#include "Quadrature.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

// All limits are explicit: an adaptive integration must fail predictably rather
// than silently replacing a cut cell with a different rule.
struct OctreeCutQuadratureOptions {
	std::uint32_t max_depth = 6;
	std::size_t max_nodes = 1000000;
	std::size_t max_leaves = 875000;
	std::size_t max_points = 4000000;
};

struct CutCellVolumeQuadratureDiagnostics {
	double certified_reference_volume = 0.0;
	double estimated_reference_volume = 0.0;
	double unresolved_reference_volume = 0.0;
	double lower_reference_volume = 0.0;
	double upper_reference_volume = 0.0;
	double certified_physical_volume = 0.0;
	double estimated_physical_volume = 0.0;
	double unresolved_physical_volume = 0.0;
	double lower_physical_volume = 0.0;
	double upper_physical_volume = 0.0;
	std::size_t nodes = 0;
	std::size_t leaves = 0;
	std::size_t output_points = 0;
	std::size_t samples = 0;
	std::size_t boundary_samples = 0;
	std::size_t ambiguous_samples = 0;
	std::size_t predicate_ambiguities = 0;
	std::size_t precision_limited_leaves = 0;
	std::uint32_t reached_depth = 0;
};

struct CutCellVolumeQuadratureCell {
	std::uint64_t id = 0;
	CellClassification classification = CellClassification::Outside;
	VolumeQuadratureRule rule;
	CutCellVolumeQuadratureDiagnostics diagnostics;
	bool usable = true;
};

class CutCellVolumeQuadratureCatalog {
public:
	CutCellVolumeQuadratureCatalog(const CutCellVolumeQuadratureCatalog&) = delete;
	CutCellVolumeQuadratureCatalog& operator=(const CutCellVolumeQuadratureCatalog&) = delete;

	explicit CutCellVolumeQuadratureCatalog(const CartesianDomainClassification& domain,
		OctreeCutQuadratureOptions options = {})
		: grid_spec_(domain.Background().Spec()), surface_canonical_hash_(domain.SurfaceCanonicalHash()),
		options_(ValidateOptions(options))
	{
		const auto& cells = domain.Cells();
		if (cells.size() != domain.Background().ElementCount())
			throw std::runtime_error("Cartesian domain catalog is incomplete");
		cells_.reserve(cells.size());
		for (const auto& cell : cells) {
			if (cell.id != cells_.size()) throw std::runtime_error("Cartesian domain cell ids are not x-fast");
			CutCellVolumeQuadratureCell result;
			result.id = cell.id;
			result.classification = cell.classification;
			BuildCell(domain, cell, result);
			Accumulate(diagnostics_, result.diagnostics);
			cells_.push_back(std::move(result));
		}
		ValidateDiagnostics(diagnostics_);
	}

	CutCellVolumeQuadratureCatalog(CutCellVolumeQuadratureCatalog&&) noexcept = default;
	CutCellVolumeQuadratureCatalog& operator=(CutCellVolumeQuadratureCatalog&&) noexcept = default;

	const CubicCartesianGridSpec& GridSpec() const noexcept { return grid_spec_; }
	const std::string& SurfaceCanonicalHash() const noexcept { return surface_canonical_hash_; }
	const OctreeCutQuadratureOptions& Options() const noexcept { return options_; }
	const std::vector<CutCellVolumeQuadratureCell>& Cells() const noexcept { return cells_; }
	const CutCellVolumeQuadratureDiagnostics& Diagnostics() const noexcept { return diagnostics_; }

	const CutCellVolumeQuadratureCell& Cell(std::uint64_t id) const
	{
		if (id >= cells_.size()) throw std::out_of_range("cut-cell quadrature id is out of range");
		return cells_[static_cast<std::size_t>(id)];
	}

	// Rules are bound to the exact classified surface and Cartesian grid that
	// produced them; a matching id by itself is not a safe assembly binding.
	const VolumeQuadratureRule& UsableRule(const CartesianDomainClassification& domain,
		std::uint64_t id) const
	{
		ValidateBinding(domain, id);
		const auto& cell = Cell(id);
		if (!cell.usable) throw std::runtime_error("cut-cell quadrature is unusable after predicate ambiguity");
		if (cell.rule.Points().empty()) throw std::runtime_error("cut-cell quadrature is certified empty");
		return cell.rule;
	}

	void ValidateUsableRule(const CartesianDomainClassification& domain, std::uint64_t id) const
	{
		ValidateBinding(domain, id);
		const auto& cell = Cell(id);
		if (!cell.usable) throw std::runtime_error("cut-cell quadrature is unusable after predicate ambiguity");
		if (cell.rule.Points().empty()) {
			if (cell.diagnostics.certified_reference_volume != 0.0
				|| cell.diagnostics.unresolved_reference_volume != 0.0)
				throw std::runtime_error("empty cut-cell rule is not certified empty");
			return;
		}
		ValidateVolumeQuadratureRule(domain.Background().MaterializeElement(id), cell.rule);
	}

private:
	struct Node {
		std::array<double, 3> reference_lower{};
		std::array<double, 3> reference_upper{};
		SurfaceAabb physical_bounds;
		std::uint32_t depth = 0;
	};

	static OctreeCutQuadratureOptions ValidateOptions(OctreeCutQuadratureOptions options)
	{
		if (options.max_depth > 20) throw std::invalid_argument("octree cut quadrature depth exceeds supported cap");
		if (!options.max_nodes || !options.max_leaves || !options.max_points)
			throw std::invalid_argument("octree cut quadrature caps must be positive");
		return options;
	}
	static void CheckAdd(std::size_t& value, std::size_t add, std::size_t cap, const char* message)
	{
		if (add > cap-value) throw std::runtime_error(message);
		value += add;
	}
	static void AddFinite(double& total, double value, const char* message)
	{
		if (!std::isfinite(value) || value < 0.0 || !std::isfinite(total)
			|| total > std::numeric_limits<double>::max()-value)
			throw std::runtime_error(message);
		total += value;
	}
	static double Volume(const std::array<double, 3>& lower, const std::array<double, 3>& upper)
	{
		double volume = 1.0;
		for (std::size_t axis = 0; axis < 3; ++axis) {
			const double width = upper[axis]-lower[axis];
			if (!std::isfinite(width) || !(width > 0.0)) throw std::runtime_error("octree node has non-positive volume");
			if (volume > std::numeric_limits<double>::max()/width) throw std::overflow_error("octree volume overflows");
			volume *= width;
		}
		return volume;
	}
	static void Accumulate(CutCellVolumeQuadratureDiagnostics& target,
		const CutCellVolumeQuadratureDiagnostics& source)
	{
		AddFinite(target.certified_reference_volume, source.certified_reference_volume, "cut quadrature diagnostic overflows");
		AddFinite(target.estimated_reference_volume, source.estimated_reference_volume, "cut quadrature diagnostic overflows");
		AddFinite(target.unresolved_reference_volume, source.unresolved_reference_volume, "cut quadrature diagnostic overflows");
		AddFinite(target.lower_reference_volume, source.lower_reference_volume, "cut quadrature diagnostic overflows");
		AddFinite(target.upper_reference_volume, source.upper_reference_volume, "cut quadrature diagnostic overflows");
		AddFinite(target.certified_physical_volume, source.certified_physical_volume, "cut quadrature diagnostic overflows");
		AddFinite(target.estimated_physical_volume, source.estimated_physical_volume, "cut quadrature diagnostic overflows");
		AddFinite(target.unresolved_physical_volume, source.unresolved_physical_volume, "cut quadrature diagnostic overflows");
		AddFinite(target.lower_physical_volume, source.lower_physical_volume, "cut quadrature diagnostic overflows");
		AddFinite(target.upper_physical_volume, source.upper_physical_volume, "cut quadrature diagnostic overflows");
		const auto add_count = [](std::size_t& left, std::size_t right) {
			if (right > std::numeric_limits<std::size_t>::max()-left)
				throw std::overflow_error("cut quadrature count overflows");
			left += right;
		};
		add_count(target.nodes, source.nodes); add_count(target.leaves, source.leaves);
		add_count(target.output_points, source.output_points); add_count(target.samples, source.samples);
		add_count(target.boundary_samples, source.boundary_samples); add_count(target.ambiguous_samples, source.ambiguous_samples);
		add_count(target.predicate_ambiguities, source.predicate_ambiguities); add_count(target.precision_limited_leaves, source.precision_limited_leaves);
		target.reached_depth = std::max(target.reached_depth, source.reached_depth);
	}
	static void ValidateDiagnostics(const CutCellVolumeQuadratureDiagnostics& diagnostics)
	{
		const auto close = [](double a, double b) { return std::abs(a-b) <= 2.0e-12*std::max({1.0, std::abs(a), std::abs(b)}); };
		if (!close(diagnostics.lower_reference_volume, diagnostics.certified_reference_volume)
			|| !close(diagnostics.lower_physical_volume, diagnostics.certified_physical_volume)
			|| diagnostics.upper_reference_volume+2.0e-12 < diagnostics.lower_reference_volume
			|| diagnostics.upper_physical_volume+2.0e-12 < diagnostics.lower_physical_volume
			|| !close(diagnostics.upper_reference_volume, diagnostics.certified_reference_volume+diagnostics.unresolved_reference_volume)
			|| !close(diagnostics.upper_physical_volume, diagnostics.certified_physical_volume+diagnostics.unresolved_physical_volume)
			|| diagnostics.estimated_reference_volume+2.0e-12 < diagnostics.lower_reference_volume
			|| diagnostics.estimated_reference_volume > diagnostics.upper_reference_volume+2.0e-12
			|| diagnostics.estimated_physical_volume+2.0e-12 < diagnostics.lower_physical_volume
			|| diagnostics.estimated_physical_volume > diagnostics.upper_physical_volume+2.0e-12)
			throw std::runtime_error("cut quadrature diagnostic bounds are inconsistent");
	}
	static void ValidateStoredRule(const CutCellVolumeQuadratureCell& cell)
	{
		double weight_sum = 0.0;
		for (const auto& point : cell.rule.Points()) {
			if (!QuadratureFinite(point.parametric) || point.parametric[0] < 0.0 || point.parametric[0] > 1.0
				|| point.parametric[1] < 0.0 || point.parametric[1] > 1.0
				|| point.parametric[2] < 0.0 || point.parametric[2] > 1.0
				|| !std::isfinite(point.weight) || !(point.weight > 0.0))
				throw std::runtime_error("stored cut-cell quadrature point is invalid");
			AddFinite(weight_sum, point.weight, "stored cut-cell quadrature weight sum overflows");
		}
		if (cell.diagnostics.output_points != cell.rule.Points().size()
			|| !QuadratureClose(weight_sum, cell.diagnostics.estimated_reference_volume))
			throw std::runtime_error("stored cut-cell quadrature weights are inconsistent");
	}
	static bool SameGridSpec(const CubicCartesianGridSpec& left, const CubicCartesianGridSpec& right)
	{
		return left.lower_m == right.lower_m && left.upper_m == right.upper_m && left.cells == right.cells;
	}
	void ValidateBinding(const CartesianDomainClassification& domain, std::uint64_t id) const
	{
		if (!SameGridSpec(grid_spec_, domain.Background().Spec())
			|| domain.Background().ElementCount() != cells_.size()
			|| domain.Cells().size() != cells_.size()
			|| domain.SurfaceCanonicalHash() != surface_canonical_hash_)
			throw std::runtime_error("cut-cell quadrature domain binding does not match catalog");
		if (id >= cells_.size() || domain.Cells()[static_cast<std::size_t>(id)].id != id
			|| cells_[static_cast<std::size_t>(id)].id != id)
			throw std::runtime_error("cut-cell quadrature id is not x-fast in bound domain");
	}
	static double PhysicalRuleVolume(const Element& element, const VolumeQuadratureRule& rule)
	{
		double result = 0.0;
		for (const auto& point : rule.Points()) {
			const auto geometry = EvaluateElementGeometry(element, point.parametric);
			if (!std::isfinite(geometry.raw_determinant) || !(geometry.raw_determinant > 0.0)
				|| !std::isfinite(point.weight) || point.weight > std::numeric_limits<double>::max()/geometry.raw_determinant)
				throw std::runtime_error("cut-cell quadrature has invalid physical Jacobian");
			AddFinite(result, point.weight*geometry.raw_determinant, "cut quadrature physical estimate overflows");
		}
		return result;
	}
	static VolumeQuadratureRule FullRule(const Element& element)
	{
		return FullCell4x4x4VolumeQuadratureProvider(element).Rule();
	}
	static bool Midpoint(const double lower, const double upper, double& midpoint)
	{
		midpoint = lower+(upper-lower)/2.0;
		return std::isfinite(midpoint) && midpoint > lower && midpoint < upper;
	}
	static void AppendScaledGauss(const Node& node, std::vector<VolumeQuadraturePoint>& points,
		CutCellVolumeQuadratureDiagnostics& diagnostics, const OctreeCutQuadratureOptions& options)
	{
		CheckAdd(diagnostics.samples, 64, std::numeric_limits<std::size_t>::max(), "cut quadrature sample count overflows");
		CheckAdd(diagnostics.output_points, 64, options.max_points, "octree cut quadrature point cap reached");
		const double volume = Volume(node.reference_lower, node.reference_upper);
		for (std::size_t qz = 0; qz < 4; ++qz)
			for (std::size_t qy = 0; qy < 4; ++qy)
				for (std::size_t qx = 0; qx < 4; ++qx) {
					std::array<double, 3> point{};
					for (std::size_t axis = 0; axis < 3; ++axis) {
						const std::size_t q = axis == 0 ? qx : (axis == 1 ? qy : qz);
						point[axis] = node.reference_lower[axis]
							+(node.reference_upper[axis]-node.reference_lower[axis])*kGaussFourPoints[q];
					}
					points.push_back({point, volume*kGaussFourWeights[qx]*kGaussFourWeights[qy]*kGaussFourWeights[qz]/8.0});
				}
	}
	static void AddCertified(Node const& node, CutCellVolumeQuadratureDiagnostics& diagnostics,
		double physical_cell_volume)
	{
		const double reference_volume = Volume(node.reference_lower, node.reference_upper);
		AddFinite(diagnostics.certified_reference_volume, reference_volume, "cut quadrature certified volume overflows");
		AddFinite(diagnostics.certified_physical_volume, reference_volume*physical_cell_volume, "cut quadrature physical volume overflows");
	}
	static void AddUnresolved(Node const& node, CutCellVolumeQuadratureDiagnostics& diagnostics,
		double physical_cell_volume)
	{
		const double reference_volume = Volume(node.reference_lower, node.reference_upper);
		AddFinite(diagnostics.unresolved_reference_volume, reference_volume, "cut quadrature unresolved volume overflows");
		AddFinite(diagnostics.unresolved_physical_volume, reference_volume*physical_cell_volume, "cut quadrature physical volume overflows");
	}

	void BuildCell(const CartesianDomainClassification& domain, const CartesianDomainCell& source,
		CutCellVolumeQuadratureCell& result) const
	{
		const auto cell = domain.Background().Cell(source.id);
		const Element element = domain.Background().MaterializeElement(source.id);
		const double physical_volume = Volume(cell.lower_m, cell.upper_m);
		if (source.classification == CellClassification::Inside) {
			result.rule = FullRule(element);
			result.diagnostics.certified_reference_volume = 1.0;
			result.diagnostics.estimated_reference_volume = 1.0;
			result.diagnostics.lower_reference_volume = 1.0;
			result.diagnostics.upper_reference_volume = 1.0;
			result.diagnostics.certified_physical_volume = physical_volume;
			result.diagnostics.estimated_physical_volume = PhysicalRuleVolume(element, result.rule);
			result.diagnostics.lower_physical_volume = physical_volume;
			result.diagnostics.upper_physical_volume = physical_volume;
			CheckAdd(result.diagnostics.nodes, 1, options_.max_nodes, "octree cut quadrature node cap reached");
			CheckAdd(result.diagnostics.leaves, 1, options_.max_leaves, "octree cut quadrature leaf cap reached");
			CheckAdd(result.diagnostics.output_points, result.rule.Points().size(), options_.max_points, "octree cut quadrature point cap reached");
			CheckAdd(result.diagnostics.samples, result.rule.Points().size(), std::numeric_limits<std::size_t>::max(), "cut quadrature sample count overflows");
			ValidateVolumeQuadratureRule(element, result.rule);
			ValidateStoredRule(result);
			ValidateDiagnostics(result.diagnostics);
			return;
		}
		if (source.classification == CellClassification::Outside) {
			CheckAdd(result.diagnostics.nodes, 1, options_.max_nodes, "octree cut quadrature node cap reached");
			CheckAdd(result.diagnostics.leaves, 1, options_.max_leaves, "octree cut quadrature leaf cap reached");
			ValidateStoredRule(result);
			ValidateDiagnostics(result.diagnostics);
			return;
		}

		std::vector<VolumeQuadraturePoint> points;
		std::vector<Node> stack;
		stack.push_back({{{0.0, 0.0, 0.0}}, {{1.0, 1.0, 1.0}}, source.bounds, 0});
		if (source.ambiguous) { result.usable = false; ++result.diagnostics.predicate_ambiguities; }
		while (!stack.empty()) {
			Node node = stack.back(); stack.pop_back();
			CheckAdd(result.diagnostics.nodes, 1, options_.max_nodes, "octree cut quadrature node cap reached");
			result.diagnostics.reached_depth = std::max(result.diagnostics.reached_depth, node.depth);
			bool ambiguous = false;
			const auto contact = domain.SurfaceIndex().IntersectBox(node.physical_bounds);
			if (contact == BoxContact::Ambiguous) ambiguous = true;
			else if (contact == BoxContact::None || contact == BoxContact::BoundaryOnly) {
				std::array<double, 3> center{};
				for (std::size_t axis = 0; axis < 3; ++axis) center[axis] = node.reference_lower[axis]+(node.reference_upper[axis]-node.reference_lower[axis])/2.0;
				const auto location = domain.SurfaceIndex().LocatePoint(EvaluateElementGeometry(element, center).physical);
				if (location == PointLocation::Inside) {
					CheckAdd(result.diagnostics.leaves, 1, options_.max_leaves, "octree cut quadrature leaf cap reached");
					AddCertified(node, result.diagnostics, physical_volume);
					AppendScaledGauss(node, points, result.diagnostics, options_);
					continue;
				}
				if (location == PointLocation::Outside) {
					CheckAdd(result.diagnostics.leaves, 1, options_.max_leaves, "octree cut quadrature leaf cap reached");
					continue;
				}
				ambiguous = location == PointLocation::Ambiguous;
			}
			if (ambiguous) { result.usable = false; ++result.diagnostics.predicate_ambiguities; }
			std::array<double, 3> ref_mid{}, physical_mid{};
			bool midpoint_ok = node.depth < options_.max_depth;
			for (std::size_t axis = 0; axis < 3 && midpoint_ok; ++axis)
				midpoint_ok = Midpoint(node.reference_lower[axis], node.reference_upper[axis], ref_mid[axis])
					&& Midpoint(node.physical_bounds.minimum[axis], node.physical_bounds.maximum[axis], physical_mid[axis]);
			if (!midpoint_ok) {
				CheckAdd(result.diagnostics.leaves, 1, options_.max_leaves, "octree cut quadrature leaf cap reached");
				if (node.depth < options_.max_depth) ++result.diagnostics.precision_limited_leaves;
				AddUnresolved(node, result.diagnostics, physical_volume);
				for (std::size_t qz = 0; qz < 4; ++qz)
					for (std::size_t qy = 0; qy < 4; ++qy)
						for (std::size_t qx = 0; qx < 4; ++qx) {
							CheckAdd(result.diagnostics.samples, 1, std::numeric_limits<std::size_t>::max(), "cut quadrature sample count overflows");
							std::array<double, 3> reference{};
							for (std::size_t axis = 0; axis < 3; ++axis) {
								const std::size_t q = axis == 0 ? qx : (axis == 1 ? qy : qz);
								reference[axis] = node.reference_lower[axis]+(node.reference_upper[axis]-node.reference_lower[axis])*kGaussFourPoints[q];
							}
							const auto location = domain.SurfaceIndex().LocatePoint(EvaluateElementGeometry(element, reference).physical);
							if (location == PointLocation::Inside) {
								CheckAdd(result.diagnostics.output_points, 1, options_.max_points, "octree cut quadrature point cap reached");
								const double weight = Volume(node.reference_lower, node.reference_upper)*kGaussFourWeights[qx]*kGaussFourWeights[qy]*kGaussFourWeights[qz]/8.0;
								points.push_back({reference, weight});
							} else if (location == PointLocation::Boundary) ++result.diagnostics.boundary_samples;
							else if (location == PointLocation::Ambiguous) { ++result.diagnostics.ambiguous_samples; ++result.diagnostics.predicate_ambiguities; result.usable = false; }
						}
				continue;
			}
			const std::size_t available_nodes = options_.max_nodes-result.diagnostics.nodes;
			if (stack.size() > available_nodes || 8 > available_nodes-stack.size())
				throw std::runtime_error("octree cut quadrature node cap reached");
			// LIFO reverse order gives Morton x-fast visitation 0..7.
			for (int child = 7; child >= 0; --child) {
				Node next; next.depth = node.depth+1;
				for (std::size_t axis = 0; axis < 3; ++axis) {
					const bool high = (child & (1 << axis)) != 0;
					next.reference_lower[axis] = high ? ref_mid[axis] : node.reference_lower[axis];
					next.reference_upper[axis] = high ? node.reference_upper[axis] : ref_mid[axis];
					next.physical_bounds.minimum[axis] = high ? physical_mid[axis] : node.physical_bounds.minimum[axis];
					next.physical_bounds.maximum[axis] = high ? node.physical_bounds.maximum[axis] : physical_mid[axis];
				}
				stack.push_back(next);
			}
		}
		result.rule = VolumeQuadratureRule(std::move(points));
		// Certified leaves are represented by their scaled Gauss weights; only
		// unresolved leaves need pointwise sampling, so this sum is exactly the
		// certified contribution plus accepted terminal samples.
		result.diagnostics.estimated_reference_volume = 0.0;
		for (const auto& point : result.rule.Points()) AddFinite(result.diagnostics.estimated_reference_volume, point.weight, "cut quadrature weight sum overflows");
		result.diagnostics.estimated_physical_volume = PhysicalRuleVolume(element, result.rule);
		result.diagnostics.lower_reference_volume = result.diagnostics.certified_reference_volume;
		result.diagnostics.upper_reference_volume = result.diagnostics.certified_reference_volume+result.diagnostics.unresolved_reference_volume;
		result.diagnostics.lower_physical_volume = result.diagnostics.certified_physical_volume;
		result.diagnostics.upper_physical_volume = result.diagnostics.certified_physical_volume+result.diagnostics.unresolved_physical_volume;
		if (result.usable && !result.rule.Points().empty()) ValidateVolumeQuadratureRule(element, result.rule);
		ValidateStoredRule(result);
		ValidateDiagnostics(result.diagnostics);
	}

	CubicCartesianGridSpec grid_spec_{};
	std::string surface_canonical_hash_;
	OctreeCutQuadratureOptions options_{};
	std::vector<CutCellVolumeQuadratureCell> cells_;
	CutCellVolumeQuadratureDiagnostics diagnostics_;
};

} // namespace iga

#endif
