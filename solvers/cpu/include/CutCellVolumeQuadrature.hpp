#ifndef IGA_CUT_CELL_VOLUME_QUADRATURE_HPP
#define IGA_CUT_CELL_VOLUME_QUADRATURE_HPP

#include "CartesianDomainClassification.hpp"
#include "Quadrature.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
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
	// Kept at the end so the original four-value aggregate initializers remain
	// source compatible.  max_points continues to cap expanded point storage.
	// Compact rules apply max_logical_points only after compression: this is the
	// number of points they would emit, not the number of retained records.
	// max_records bounds monotone append attempts, including rollback work.
	std::size_t max_records = 1000000;
	std::size_t max_retained_bytes = 256u*1024u*1024u;
	std::size_t max_logical_points = 64u*1024u*1024u;
	// Disabled by default.  A positive value retries an otherwise unresolved
	// empty Cut rule at successive absolute depths through this bounded cap.
	std::uint32_t empty_rule_rescue_max_depth = 0;
};

enum class CutCellVolumeQuadratureStorageMode { Expanded, Compact };
using CutCellVolumeStorageMode = CutCellVolumeQuadratureStorageMode;

// Coordinates are integers on the 2^max_depth reference lattice.  Keeping
// this representation through construction makes coalescing exact and avoids
// making a large transient vector of quadrature points.
struct CompactCutCellVolumeBlock {
	std::array<std::uint32_t, 3> lower{};
	std::array<std::uint32_t, 3> upper{};
};

struct CompactCutCellVolumeSampleLeaf {
	std::array<std::uint32_t, 3> key{};
	std::uint32_t depth = 0;
	std::uint64_t inside_mask = 0;
};

struct CompactCutCellVolumeRule {
	std::uint32_t max_depth = 0;
	std::vector<CompactCutCellVolumeBlock> certified_blocks;
	std::vector<CompactCutCellVolumeSampleLeaf> sample_leaves;
	// A fitted rule is already small and has non-lattice nodes. It exclusively
	// uses explicit positive points, preserving their expanded order and bits.
	std::vector<VolumeQuadraturePoint> fitted_points;
};

inline void ValidateCompactFittedVolumePoints(const CompactCutCellVolumeRule& rule)
{
	if(rule.fitted_points.empty())return;
	if(!rule.certified_blocks.empty()||!rule.sample_leaves.empty())
		throw std::invalid_argument("compact fitted points cannot mix with octree records");
	for(std::size_t i=0;i<rule.fitted_points.size();++i) {
		const auto& point=rule.fitted_points[i];
		if(!std::isfinite(point.weight)||!(point.weight>0))throw std::invalid_argument("compact fitted point weight is invalid");
		for(auto coordinate:point.parametric)if(!std::isfinite(coordinate)||coordinate<0||coordinate>1)
			throw std::invalid_argument("compact fitted point coordinate is invalid");
		if(i&&!(rule.fitted_points[i-1].parametric<point.parametric))
			throw std::invalid_argument("compact fitted points are not strictly ordered");
	}
}

inline std::size_t CompactCutCellVolumeLogicalPointCount(const CompactCutCellVolumeRule& rule)
{
	ValidateCompactFittedVolumePoints(rule);
	if(!rule.fitted_points.empty())return rule.fitted_points.size();
	if (rule.certified_blocks.size() > std::numeric_limits<std::size_t>::max()/64)
		throw std::overflow_error("compact cut-cell logical point count overflows");
	std::size_t result = 64*rule.certified_blocks.size();
	for (const auto& leaf : rule.sample_leaves) {
		const std::size_t count = static_cast<std::size_t>(__builtin_popcountll(leaf.inside_mask));
		if (count > std::numeric_limits<std::size_t>::max()-result)
			throw std::overflow_error("compact cut-cell logical point count overflows");
		result += count;
	}
	return result;
}

template <class Callback> void ForEachVolumePoint(const CompactCutCellVolumeRule& rule, Callback&& callback)
{
	if (rule.max_depth > 20) throw std::invalid_argument("compact cut-cell depth exceeds supported cap");
	ValidateCompactFittedVolumePoints(rule);
	if(!rule.fitted_points.empty()) {
		for(const auto& point:rule.fitted_points)callback(point);
		return;
	}
	const double lattice = static_cast<double>(std::uint64_t(1) << rule.max_depth);
	const auto emit = [&](const std::array<std::uint32_t, 3>& lower, const std::array<std::uint32_t, 3>& upper,
		std::uint64_t mask) {
		const double volume = (static_cast<double>(upper[0]-lower[0])/lattice)
			*(static_cast<double>(upper[1]-lower[1])/lattice)*(static_cast<double>(upper[2]-lower[2])/lattice);
		for (std::size_t qz = 0; qz < 4; ++qz) for (std::size_t qy = 0; qy < 4; ++qy) for (std::size_t qx = 0; qx < 4; ++qx) {
			const std::size_t bit = qx+4*qy+16*qz;
			if ((mask & (std::uint64_t(1) << bit)) == 0) continue;
			std::array<double, 3> point{};
			for (std::size_t axis = 0; axis < 3; ++axis) {
				const std::size_t q = axis == 0 ? qx : (axis == 1 ? qy : qz);
				point[axis] = static_cast<double>(lower[axis])/lattice
					+static_cast<double>(upper[axis]-lower[axis])/lattice*kGaussFourPoints[q];
			}
			callback(VolumeQuadraturePoint{point, volume*kGaussFourWeights[qx]*kGaussFourWeights[qy]*kGaussFourWeights[qz]/8.0});
		}
	};
	for (const auto& block : rule.certified_blocks) emit(block.lower, block.upper, ~std::uint64_t(0));
	for (const auto& leaf : rule.sample_leaves) {
		if (leaf.depth > rule.max_depth) throw std::invalid_argument("compact cut-cell sample leaf depth is invalid");
		const std::uint32_t width = std::uint32_t(1) << (rule.max_depth-leaf.depth);
		std::array<std::uint32_t, 3> lower{}, upper{};
		for (std::size_t axis = 0; axis < 3; ++axis) { lower[axis] = leaf.key[axis]*width; upper[axis] = lower[axis]+width; }
		emit(lower, upper, leaf.inside_mask);
	}
}

inline double CompactCutCellVolumeWeightSum(const CompactCutCellVolumeRule& rule)
{
	double result = 0.0;
	ForEachVolumePoint(rule, [&](const VolumeQuadraturePoint& point) { result += point.weight; });
	return result;
}

inline double CompactCutCellVolumeMoment(const CompactCutCellVolumeRule& rule, int x, int y, int z)
{
	if (x < 0 || y < 0 || z < 0) throw std::invalid_argument("compact cut-cell moment exponent is negative");
	double result = 0.0;
	ForEachVolumePoint(rule, [&](const VolumeQuadraturePoint& point) {
		result += point.weight*std::pow(point.parametric[0], x)*std::pow(point.parametric[1], y)*std::pow(point.parametric[2], z);
	});
	return result;
}

inline bool CompactCutCellVolumeBlocksOverlap(const CompactCutCellVolumeBlock& left,
	const CompactCutCellVolumeBlock& right)
{
	for (std::size_t axis = 0; axis < 3; ++axis)
		if (left.lower[axis] >= right.upper[axis] || right.lower[axis] >= left.upper[axis]) return false;
	return true;
}

inline std::uint64_t CompactMortonLower(const std::array<std::uint32_t, 3>& lower, std::uint32_t depth)
{
	std::uint64_t result = 0;
	for (std::uint32_t level = 0; level < depth; ++level) {
		const std::uint32_t bit = depth-1-level;
		result = (result << 3) | ((std::uint64_t((lower[0] >> bit) & 1u))
			| (std::uint64_t((lower[1] >> bit) & 1u) << 1) | (std::uint64_t((lower[2] >> bit) & 1u) << 2));
	}
	return result;
}

struct CompactCutCellVolumeBvhNode {
	CompactCutCellVolumeBlock bounds{};
	std::size_t begin = 0, end = 0;
	int left = -1, right = -1;
};

inline int BuildCompactCutCellVolumeBvh(const std::vector<CompactCutCellVolumeBlock>& boxes,
	std::vector<std::size_t>& order, std::vector<CompactCutCellVolumeBvhNode>& nodes, std::size_t begin, std::size_t end)
{
	CompactCutCellVolumeBvhNode node; node.begin = begin; node.end = end;
	for (std::size_t axis = 0; axis < 3; ++axis) { node.bounds.lower[axis] = std::numeric_limits<std::uint32_t>::max(); node.bounds.upper[axis] = 0; }
	for (std::size_t i = begin; i < end; ++i) for (std::size_t axis = 0; axis < 3; ++axis) {
		node.bounds.lower[axis] = std::min(node.bounds.lower[axis], boxes[order[i]].lower[axis]);
		node.bounds.upper[axis] = std::max(node.bounds.upper[axis], boxes[order[i]].upper[axis]);
	}
	const int index = static_cast<int>(nodes.size()); nodes.push_back(node);
	if (end-begin <= 8) return index;
	std::size_t axis = 0;
	for (std::size_t candidate = 1; candidate < 3; ++candidate)
		if (node.bounds.upper[candidate]-node.bounds.lower[candidate] > node.bounds.upper[axis]-node.bounds.lower[axis]) axis = candidate;
	const std::size_t middle = begin+(end-begin)/2;
	std::nth_element(order.begin()+static_cast<std::ptrdiff_t>(begin), order.begin()+static_cast<std::ptrdiff_t>(middle), order.begin()+static_cast<std::ptrdiff_t>(end),
		[&](std::size_t left, std::size_t right) {
			const std::uint64_t left_center = std::uint64_t(boxes[left].lower[axis])+boxes[left].upper[axis];
			const std::uint64_t right_center = std::uint64_t(boxes[right].lower[axis])+boxes[right].upper[axis];
			return left_center == right_center ? left < right : left_center < right_center;
		});
	nodes[index].left = BuildCompactCutCellVolumeBvh(boxes, order, nodes, begin, middle);
	nodes[index].right = BuildCompactCutCellVolumeBvh(boxes, order, nodes, middle, end);
	return index;
}

inline bool CompactCutCellVolumeBvhOverlaps(const std::vector<CompactCutCellVolumeBlock>& boxes,
	const std::vector<std::size_t>& order, const std::vector<CompactCutCellVolumeBvhNode>& nodes, int node_index,
	const CompactCutCellVolumeBlock& query, std::size_t skip = std::numeric_limits<std::size_t>::max())
{
	const auto& node = nodes[static_cast<std::size_t>(node_index)];
	if (!CompactCutCellVolumeBlocksOverlap(node.bounds, query)) return false;
	if (node.left < 0) {
		for (std::size_t i = node.begin; i < node.end; ++i)
			if (order[i] != skip && CompactCutCellVolumeBlocksOverlap(boxes[order[i]], query)) return true;
		return false;
	}
	return CompactCutCellVolumeBvhOverlaps(boxes, order, nodes, node.left, query, skip)
		|| CompactCutCellVolumeBvhOverlaps(boxes, order, nodes, node.right, query, skip);
}

inline void ValidateCompactCutCellVolumeRule(const CompactCutCellVolumeRule& rule)
{
	if (rule.max_depth > 20) throw std::invalid_argument("compact cut-cell depth exceeds supported cap");
	ValidateCompactFittedVolumePoints(rule);
	if(!rule.fitted_points.empty())return;
	const std::uint64_t lattice = std::uint64_t(1) << rule.max_depth;
	const auto valid_box = [&](const std::array<std::uint32_t, 3>& lower, const std::array<std::uint32_t, 3>& upper) {
		for (std::size_t axis = 0; axis < 3; ++axis)
			if (lower[axis] >= upper[axis] || upper[axis] > lattice) throw std::invalid_argument("compact cut-cell block is invalid");
	};
	for (std::size_t i = 0; i < rule.certified_blocks.size(); ++i) {
		valid_box(rule.certified_blocks[i].lower, rule.certified_blocks[i].upper);
		if (i && !(rule.certified_blocks[i-1].lower < rule.certified_blocks[i].lower
			|| (rule.certified_blocks[i-1].lower == rule.certified_blocks[i].lower
				&& rule.certified_blocks[i-1].upper < rule.certified_blocks[i].upper)))
			throw std::invalid_argument("compact cut-cell blocks are not strictly ordered");
	}
	std::vector<std::pair<std::uint64_t, CompactCutCellVolumeBlock>> sample_boxes;
	sample_boxes.reserve(rule.sample_leaves.size());
	for (const auto& leaf : rule.sample_leaves) {
		if (leaf.depth > rule.max_depth) throw std::invalid_argument("compact cut-cell sample leaf depth is invalid");
		const std::uint32_t span = std::uint32_t(1) << leaf.depth;
		const std::uint32_t width = std::uint32_t(1) << (rule.max_depth-leaf.depth);
		CompactCutCellVolumeBlock box;
		for (std::size_t axis = 0; axis < 3; ++axis) {
			if (leaf.key[axis] >= span) throw std::invalid_argument("compact cut-cell sample leaf key is invalid");
			box.lower[axis] = leaf.key[axis]*width; box.upper[axis] = box.lower[axis]+width;
		}
		sample_boxes.push_back({CompactMortonLower(box.lower, rule.max_depth), box});
	}
	for (std::size_t i = 1; i < sample_boxes.size(); ++i)
		if (sample_boxes[i-1].first >= sample_boxes[i].first)
			throw std::invalid_argument("compact cut-cell sample leaves are not strictly Morton ordered");
	// Terminal sample leaves are octree cells.  Morton intervals are therefore
	// disjoint iff their starts are strictly increasing in the construction's
	// x-fast visitation; this avoids pairwise scans at depth eight.
	std::vector<std::pair<std::uint64_t, std::uint64_t>> intervals;
	intervals.reserve(sample_boxes.size());
	for (const auto& entry : sample_boxes) {
		const auto& box = entry.second;
		const std::uint32_t width = box.upper[0]-box.lower[0];
		if (box.upper[1]-box.lower[1] != width || box.upper[2]-box.lower[2] != width)
			throw std::invalid_argument("compact cut-cell sample leaf is not cubic");
		const std::uint32_t leaf_depth = rule.max_depth-static_cast<std::uint32_t>(std::log2(width));
		const std::uint64_t length = std::uint64_t(1) << (3*(rule.max_depth-leaf_depth));
		intervals.push_back({entry.first, entry.first+length});
	}
	for (std::size_t i = 1; i < intervals.size(); ++i)
		if (intervals[i-1].second > intervals[i].first) throw std::invalid_argument("compact cut-cell sample leaves overlap");
	if (!rule.certified_blocks.empty()) {
		std::vector<std::size_t> order(rule.certified_blocks.size());
		for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
		std::vector<CompactCutCellVolumeBvhNode> nodes; nodes.reserve(2*order.size());
		const int root = BuildCompactCutCellVolumeBvh(rule.certified_blocks, order, nodes, 0, order.size());
		for (std::size_t i = 0; i < rule.certified_blocks.size(); ++i)
			if (CompactCutCellVolumeBvhOverlaps(rule.certified_blocks, order, nodes, root, rule.certified_blocks[i], i))
				throw std::invalid_argument("compact cut-cell blocks overlap");
		for (const auto& sample : sample_boxes)
			if (CompactCutCellVolumeBvhOverlaps(rule.certified_blocks, order, nodes, root, sample.second))
				throw std::invalid_argument("compact cut-cell block and sample leaf overlap");
	}
	ForEachVolumePoint(rule, [](const VolumeQuadraturePoint& point) {
		if (!QuadratureFinite(point.parametric) || point.parametric[0] < 0.0 || point.parametric[0] > 1.0
			|| point.parametric[1] < 0.0 || point.parametric[1] > 1.0 || point.parametric[2] < 0.0 || point.parametric[2] > 1.0
			|| !std::isfinite(point.weight) || !(point.weight > 0.0))
			throw std::invalid_argument("compact cut-cell point is invalid");
	});
}

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
	std::size_t certified_blocks = 0;
	std::size_t sample_leaves = 0;
	// Compact-only monotone append count.  This includes records discarded by
	// uniform-subtree rollback, so max_records cannot be bypassed by collapse.
	std::size_t record_attempts = 0;
	// Compact-only number of speculative records discarded by recursive uniform
	// subtree rollback.  This excludes fixed-point block coalescing.
	std::size_t rolled_back_records = 0;
	// Compact-only deterministic compact-record accounting peak during this
	// cell's construction, including discarded empty-rule rescue attempts.  It
	// counts the planned 1.5x capacities and both old and replacement plans
	// while reserve may reallocate.  This is the portable quantity constrained
	// by max_retained_bytes; it is not an allocator-byte peak because
	// std::vector::reserve may overallocate.
	std::size_t retained_bytes = 0;
	// Compact-only observed final std::vector capacity bytes for the accepted
	// rule.  This diagnostic is intentionally not capped and does not report a
	// transient allocator peak or discarded rescue-attempt capacity.
	std::size_t observed_retained_bytes = 0;
	std::size_t logical_output_points = 0;
	std::uint32_t reached_depth = 0;
	// Rescue attempts include the nominal build.  The work counters are
	// monotone across discarded retries; retained rule counters above describe
	// only the accepted rule.
	std::size_t rescue_attempts = 0;
	std::size_t attempted_nodes = 0;
	std::size_t attempted_leaves = 0;
	std::size_t attempted_output_points = 0;
	std::size_t attempted_samples = 0;
	std::size_t attempted_record_attempts = 0;
	std::size_t attempted_logical_output_points = 0;
	std::uint32_t rescue_effective_depth = 0;
};

struct CutCellVolumeQuadratureCell {
	std::uint64_t id = 0;
	CellClassification classification = CellClassification::Outside;
	VolumeQuadratureRule rule;
	CompactCutCellVolumeRule compact_rule;
	CutCellVolumeQuadratureDiagnostics diagnostics;
	bool usable = true;
};

class CutCellVolumeQuadratureCatalog {
public:
	using StorageModeType = CutCellVolumeQuadratureStorageMode;
	CutCellVolumeQuadratureCatalog(const CutCellVolumeQuadratureCatalog&) = delete;
	CutCellVolumeQuadratureCatalog& operator=(const CutCellVolumeQuadratureCatalog&) = delete;

	explicit CutCellVolumeQuadratureCatalog(const CartesianDomainClassification& domain,
		OctreeCutQuadratureOptions options = {},
		CutCellVolumeQuadratureStorageMode storage_mode = CutCellVolumeQuadratureStorageMode::Expanded)
		: grid_spec_(domain.Background().Spec()), surface_canonical_hash_(domain.SurfaceCanonicalHash()),
		options_(ValidateOptions(options)), storage_mode_(storage_mode)
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
			if (storage_mode_ == CutCellVolumeQuadratureStorageMode::Expanded) BuildCellWithRescue(domain, cell, result);
			else BuildCompactCellWithRescue(domain, cell, result);
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
	CutCellVolumeQuadratureStorageMode StorageMode() const noexcept { return storage_mode_; }
	const std::vector<CutCellVolumeQuadratureCell>& Cells() const noexcept { return cells_; }
	const CutCellVolumeQuadratureDiagnostics& Diagnostics() const noexcept { return diagnostics_; }

	const CutCellVolumeQuadratureCell& Cell(std::uint64_t id) const
	{
		if (id >= cells_.size()) throw std::out_of_range("cut-cell quadrature id is out of range");
		return cells_[static_cast<std::size_t>(id)];
	}

#ifdef IGA_EXACT_DYADIC_TESTING
	// Focused test-only probe for the same recursive rollback helper used below.
	// Eight speculative child blocks collapse before the parent record is
	// appended, so the cap must still include all nine append attempts.
	static CutCellVolumeQuadratureDiagnostics CompactRollbackProbe(std::size_t max_records)
	{
		CutCellVolumeQuadratureDiagnostics diagnostics;
		std::vector<CompactCutCellVolumeBlock> blocks;
		std::vector<CompactCutCellVolumeSampleLeaf> samples;
		for (std::size_t child = 0; child < 8; ++child) {
			CheckAdd(diagnostics.record_attempts, 1, max_records, "compact cut-cell record cap reached");
			blocks.push_back({{{static_cast<std::uint32_t>(child),0,0}}, {{static_cast<std::uint32_t>(child+1),1,1}}});
		}
		RollbackCompactRecords(blocks, samples, 0, 0, diagnostics, max_records);
		CheckAdd(diagnostics.record_attempts, 1, max_records, "compact cut-cell record cap reached");
		blocks.push_back({{{0,0,0}}, {{8,1,1}}});
		diagnostics.certified_blocks = blocks.size();
		return diagnostics;
	}
#endif

	// Rules are bound to the exact classified surface and Cartesian grid that
	// produced them; a matching id by itself is not a safe assembly binding.
	const VolumeQuadratureRule& UsableRule(const CartesianDomainClassification& domain,
		std::uint64_t id) const
	{
		ValidateBinding(domain, id);
		if (storage_mode_ != CutCellVolumeQuadratureStorageMode::Expanded)
			throw std::runtime_error("expanded cut-cell rule requested from compact catalog");
		const auto& cell = Cell(id);
		if (!cell.usable) throw std::runtime_error("cut-cell quadrature is unusable after predicate ambiguity");
		if (cell.rule.Points().empty()) throw std::runtime_error("cut-cell quadrature is certified empty");
		return cell.rule;
	}

	const CompactCutCellVolumeRule& UsableCompactRule(const CartesianDomainClassification& domain,
		std::uint64_t id) const
	{
		ValidateBinding(domain, id);
		if (storage_mode_ != CutCellVolumeQuadratureStorageMode::Compact)
			throw std::runtime_error("compact cut-cell rule requested from expanded catalog");
		const auto& cell = Cell(id);
		if (!cell.usable) throw std::runtime_error("cut-cell quadrature is unusable after predicate ambiguity");
		if (CompactCutCellVolumeLogicalPointCount(cell.compact_rule) == 0)
			throw std::runtime_error("cut-cell quadrature is certified empty");
		return cell.compact_rule;
	}

	void ValidateUsableRule(const CartesianDomainClassification& domain, std::uint64_t id) const
	{
		ValidateBinding(domain, id);
		if (storage_mode_ != CutCellVolumeQuadratureStorageMode::Expanded)
			throw std::runtime_error("expanded cut-cell rule requested from compact catalog");
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

	void ValidateUsableCompactRule(const CartesianDomainClassification& domain, std::uint64_t id) const
	{
		ValidateBinding(domain, id);
		if (storage_mode_ != CutCellVolumeQuadratureStorageMode::Compact)
			throw std::runtime_error("compact cut-cell rule requested from expanded catalog");
		const auto& cell = Cell(id);
		if (!cell.usable) throw std::runtime_error("cut-cell quadrature is unusable after predicate ambiguity");
		const Element element = domain.Background().MaterializeElement(id);
		ValidateCompactStoredRule(cell, &element);
	}

private:
	struct Node {
		std::array<double, 3> reference_lower{};
		std::array<double, 3> reference_upper{};
		SurfaceAabb physical_bounds;
		std::uint32_t depth = 0;
		std::array<std::uint32_t, 3> key{};
	};

	enum class CompactNodeState { UniformFull, UniformEmpty, Mixed };
	static OctreeCutQuadratureOptions ValidateOptions(OctreeCutQuadratureOptions options)
	{
		if (options.max_depth > 20) throw std::invalid_argument("octree cut quadrature depth exceeds supported cap");
		if (options.empty_rule_rescue_max_depth
			&& (options.empty_rule_rescue_max_depth <= options.max_depth || options.empty_rule_rescue_max_depth > 20))
			throw std::invalid_argument("octree empty-rule rescue depth must exceed max_depth and not exceed supported cap");
		if (!options.max_nodes || !options.max_leaves || !options.max_points
			|| !options.max_records || !options.max_retained_bytes || !options.max_logical_points)
			throw std::invalid_argument("octree cut quadrature caps must be positive");
		return options;
	}
	static void CheckAdd(std::size_t& value, std::size_t add, std::size_t cap, const char* message)
	{
		if (value > cap || add > cap-value) throw std::runtime_error(message);
		value += add;
	}
	static void RollbackCompactRecords(std::vector<CompactCutCellVolumeBlock>& blocks,
		std::vector<CompactCutCellVolumeSampleLeaf>& samples, std::size_t block_start, std::size_t sample_start,
		CutCellVolumeQuadratureDiagnostics& diagnostics, std::size_t max_records)
	{
		if (block_start > blocks.size() || sample_start > samples.size())
			throw std::logic_error("compact cut-cell rollback range is invalid");
		const std::size_t block_records = blocks.size()-block_start;
		const std::size_t sample_records = samples.size()-sample_start;
		if (sample_records > std::numeric_limits<std::size_t>::max()-block_records)
			throw std::overflow_error("compact cut-cell record rollback overflows");
		CheckAdd(diagnostics.rolled_back_records, block_records+sample_records, max_records,
			"compact cut-cell record rollback overflows");
		blocks.resize(block_start); samples.resize(sample_start);
	}
	static void AddFinite(double& total, double value, const char* message)
	{
		if (!std::isfinite(value) || value < 0.0 || !std::isfinite(total)
			|| total > std::numeric_limits<double>::max()-value)
			throw std::runtime_error(message);
		total += value;
	}
	static void AddCompensatedNonnegative(double& total, double& compensation, double value, const char* message)
	{
		if (!std::isfinite(value) || value < 0.0 || !std::isfinite(total) || !std::isfinite(compensation))
			throw std::runtime_error(message);
		const double corrected = value-compensation;
		const double next = total+corrected;
		if (!std::isfinite(next)) throw std::runtime_error(message);
		compensation = (next-total)-corrected;
		total = next;
	}
	static void ValidateReferenceEstimate(double estimate)
	{
		if (!std::isfinite(estimate) || estimate < 0.0 || estimate > 1.0)
			throw std::runtime_error("cut quadrature reference estimate is outside unit interval");
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
		add_count(target.certified_blocks, source.certified_blocks); add_count(target.sample_leaves, source.sample_leaves);
		add_count(target.record_attempts, source.record_attempts); add_count(target.rolled_back_records, source.rolled_back_records);
		add_count(target.retained_bytes, source.retained_bytes); add_count(target.observed_retained_bytes, source.observed_retained_bytes);
		add_count(target.logical_output_points, source.logical_output_points);
		add_count(target.rescue_attempts, source.rescue_attempts);
		add_count(target.attempted_nodes, source.attempted_nodes); add_count(target.attempted_leaves, source.attempted_leaves);
		add_count(target.attempted_output_points, source.attempted_output_points); add_count(target.attempted_samples, source.attempted_samples);
		add_count(target.attempted_record_attempts, source.attempted_record_attempts);
		add_count(target.attempted_logical_output_points, source.attempted_logical_output_points);
		target.reached_depth = std::max(target.reached_depth, source.reached_depth);
		target.rescue_effective_depth = std::max(target.rescue_effective_depth, source.rescue_effective_depth);
	}
	static void ValidateDiagnostics(const CutCellVolumeQuadratureDiagnostics& diagnostics)
	{
		const double operations = 1.0+std::sqrt(static_cast<double>(diagnostics.samples)+1.0);
		const auto close = [operations](double a, double b) {
			const double scale = std::max(std::numeric_limits<double>::min(), std::abs(a)+std::abs(b));
			return std::abs(a-b) <= 64.0*std::numeric_limits<double>::epsilon()*operations*scale;
		};
		const auto slack = [operations](double a, double b) {
			return 64.0*std::numeric_limits<double>::epsilon()*operations
				*std::max(std::numeric_limits<double>::min(), std::abs(a)+std::abs(b));
		};
		if (!close(diagnostics.lower_reference_volume, diagnostics.certified_reference_volume)
			|| !close(diagnostics.lower_physical_volume, diagnostics.certified_physical_volume)
			|| diagnostics.upper_reference_volume+slack(diagnostics.upper_reference_volume, diagnostics.lower_reference_volume) < diagnostics.lower_reference_volume
			|| diagnostics.upper_physical_volume+slack(diagnostics.upper_physical_volume, diagnostics.lower_physical_volume) < diagnostics.lower_physical_volume
			|| !close(diagnostics.upper_reference_volume, diagnostics.certified_reference_volume+diagnostics.unresolved_reference_volume)
			|| !close(diagnostics.upper_physical_volume, diagnostics.certified_physical_volume+diagnostics.unresolved_physical_volume))
			throw std::runtime_error("cut quadrature diagnostic bound sum is inconsistent");
		if (diagnostics.estimated_reference_volume+slack(diagnostics.estimated_reference_volume, diagnostics.lower_reference_volume) < diagnostics.lower_reference_volume
			|| diagnostics.estimated_reference_volume > diagnostics.upper_reference_volume+slack(diagnostics.estimated_reference_volume, diagnostics.upper_reference_volume))
			throw std::runtime_error("cut quadrature reference estimate is outside diagnostic bounds");
		if (diagnostics.estimated_physical_volume+slack(diagnostics.estimated_physical_volume, diagnostics.lower_physical_volume) < diagnostics.lower_physical_volume
			|| diagnostics.estimated_physical_volume > diagnostics.upper_physical_volume+slack(diagnostics.estimated_physical_volume, diagnostics.upper_physical_volume))
			throw std::runtime_error("cut quadrature physical estimate is outside diagnostic bounds");
	}
	static void ValidateStoredRule(const CutCellVolumeQuadratureCell& cell)
	{
		double weight_sum = 0.0, weight_compensation = 0.0;
		for (const auto& point : cell.rule.Points()) {
			if (!QuadratureFinite(point.parametric) || point.parametric[0] < 0.0 || point.parametric[0] > 1.0
				|| point.parametric[1] < 0.0 || point.parametric[1] > 1.0
				|| point.parametric[2] < 0.0 || point.parametric[2] > 1.0
				|| !std::isfinite(point.weight) || !(point.weight > 0.0))
				throw std::runtime_error("stored cut-cell quadrature point is invalid");
			AddCompensatedNonnegative(weight_sum, weight_compensation, point.weight,
				"stored cut-cell quadrature weight sum overflows");
		}
		if (cell.diagnostics.output_points != cell.rule.Points().size()
			|| cell.diagnostics.logical_output_points != cell.rule.Points().size()
			|| !QuadratureClose(weight_sum, cell.diagnostics.estimated_reference_volume))
			throw std::runtime_error("stored cut-cell quadrature weights are inconsistent");
	}
	static void ValidateCompactStoredRule(const CutCellVolumeQuadratureCell& cell, const Element* element = nullptr,
		bool permit_provisional_empty = false)
	{
		const auto& rule = cell.compact_rule;
		try { ValidateCompactCutCellVolumeRule(rule); }
		catch (const std::exception&) { throw std::runtime_error("stored compact cut-cell rule is invalid"); }
		double weight_sum = 0.0, physical_sum = 0.0, weight_compensation = 0.0, physical_compensation = 0.0;
		ForEachVolumePoint(rule, [&](const VolumeQuadraturePoint& point) {
			if (!QuadratureFinite(point.parametric) || point.parametric[0] < 0.0 || point.parametric[0] > 1.0
				|| point.parametric[1] < 0.0 || point.parametric[1] > 1.0 || point.parametric[2] < 0.0 || point.parametric[2] > 1.0
				|| !std::isfinite(point.weight) || !(point.weight > 0.0))
				throw std::runtime_error("stored compact cut-cell point is invalid");
			AddCompensatedNonnegative(weight_sum, weight_compensation, point.weight,
				"stored compact cut-cell weight sum overflows");
			if (element) {
				const auto geometry = EvaluateElementGeometry(*element, point.parametric);
				if (!std::isfinite(geometry.raw_determinant) || !(geometry.raw_determinant > 0.0)
					|| point.weight > std::numeric_limits<double>::max()/geometry.raw_determinant)
					throw std::runtime_error("stored compact cut-cell physical Jacobian is invalid");
				AddCompensatedNonnegative(physical_sum, physical_compensation, point.weight*geometry.raw_determinant,
					"stored compact cut-cell physical volume overflows");
			}
		});
		if (!std::isfinite(weight_sum) || !QuadratureClose(weight_sum, cell.diagnostics.estimated_reference_volume)
			|| CompactCutCellVolumeLogicalPointCount(rule) != cell.diagnostics.logical_output_points)
			throw std::runtime_error("stored compact cut-cell weights are inconsistent");
		if (!permit_provisional_empty && CompactCutCellVolumeLogicalPointCount(rule) == 0
			&& (cell.diagnostics.certified_reference_volume != 0.0
				|| cell.diagnostics.unresolved_reference_volume != 0.0))
			throw std::runtime_error("empty compact cut-cell rule is not certified empty");
		if (rule.certified_blocks.size() > std::numeric_limits<std::size_t>::max()-rule.sample_leaves.size()
			|| cell.diagnostics.record_attempts < rule.certified_blocks.size()+rule.sample_leaves.size())
			throw std::runtime_error("stored compact cut-cell record diagnostics are inconsistent");
		if (cell.diagnostics.rolled_back_records > cell.diagnostics.record_attempts)
			throw std::runtime_error("stored compact cut-cell rollback diagnostics are inconsistent");
		const auto capacity_bytes = [](std::size_t block_capacity, std::size_t sample_capacity) {
			if (block_capacity > std::numeric_limits<std::size_t>::max()/sizeof(CompactCutCellVolumeBlock)
				|| sample_capacity > std::numeric_limits<std::size_t>::max()/sizeof(CompactCutCellVolumeSampleLeaf))
				throw std::runtime_error("stored compact cut-cell capacity overflows");
			const std::size_t block_bytes = block_capacity*sizeof(CompactCutCellVolumeBlock);
			const std::size_t sample_bytes = sample_capacity*sizeof(CompactCutCellVolumeSampleLeaf);
			if (sample_bytes > std::numeric_limits<std::size_t>::max()-block_bytes)
				throw std::runtime_error("stored compact cut-cell capacity overflows");
			return block_bytes+sample_bytes;
		};
		if (cell.diagnostics.observed_retained_bytes
			< capacity_bytes(rule.certified_blocks.capacity(), rule.sample_leaves.capacity()))
			throw std::runtime_error("stored compact cut-cell observed retained byte diagnostics are inconsistent");
		if (element && (!std::isfinite(physical_sum) || !QuadratureClose(physical_sum, cell.diagnostics.estimated_physical_volume)))
			throw std::runtime_error("stored compact cut-cell physical volume is inconsistent");
		ValidateDiagnostics(cell.diagnostics);
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
		double result = 0.0, compensation = 0.0;
		for (const auto& point : rule.Points()) {
			const auto geometry = EvaluateElementGeometry(element, point.parametric);
			if (!std::isfinite(geometry.raw_determinant) || !(geometry.raw_determinant > 0.0)
				|| !std::isfinite(point.weight) || point.weight > std::numeric_limits<double>::max()/geometry.raw_determinant)
				throw std::runtime_error("cut-cell quadrature has invalid physical Jacobian");
			AddCompensatedNonnegative(result, compensation, point.weight*geometry.raw_determinant,
				"cut quadrature physical estimate overflows");
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
	static std::size_t RemainingCap(std::size_t cap, std::size_t used, const char* message)
	{
		if (used > cap) throw std::runtime_error(message);
		return cap-used;
	}
	static void AddAttemptWork(CutCellVolumeQuadratureDiagnostics& target,
		const CutCellVolumeQuadratureDiagnostics& attempt, std::size_t max_retained_bytes)
	{
		const auto add = [](std::size_t& left, std::size_t right) {
			if (right > std::numeric_limits<std::size_t>::max()-left)
				throw std::overflow_error("cut quadrature rescue work count overflows");
			left += right;
		};
		if (target.rescue_attempts == std::numeric_limits<std::size_t>::max())
			throw std::overflow_error("cut quadrature rescue work count overflows");
		++target.rescue_attempts;
		add(target.attempted_nodes, attempt.nodes); add(target.attempted_leaves, attempt.leaves);
		add(target.attempted_output_points, attempt.output_points); add(target.attempted_samples, attempt.samples);
		add(target.attempted_record_attempts, attempt.record_attempts);
		add(target.attempted_logical_output_points, attempt.logical_output_points);
		if (attempt.retained_bytes > max_retained_bytes)
			throw std::runtime_error("compact cut-cell retained byte cap reached");
		target.retained_bytes = std::max(target.retained_bytes, attempt.retained_bytes);
	}
	OctreeCutQuadratureOptions AttemptOptions(std::uint32_t depth,
		const CutCellVolumeQuadratureDiagnostics& work, bool compact) const
	{
		OctreeCutQuadratureOptions result = options_;
		result.max_depth = depth;
		result.max_nodes = RemainingCap(options_.max_nodes, work.attempted_nodes,
			"octree cut quadrature node cap reached");
		result.max_leaves = RemainingCap(options_.max_leaves, work.attempted_leaves,
			"octree cut quadrature leaf cap reached");
		result.max_points = RemainingCap(options_.max_points, work.attempted_output_points,
			"octree cut quadrature point cap reached");
		if (compact) {
			result.max_records = RemainingCap(options_.max_records, work.attempted_record_attempts,
				"compact cut-cell record cap reached");
			result.max_logical_points = RemainingCap(options_.max_logical_points, work.attempted_logical_output_points,
				"compact cut-cell logical point cap reached");
		}
		return result;
	}
	static bool IsExactCertifiedEmpty(const CutCellVolumeQuadratureCell& cell, bool compact)
	{
		const std::size_t logical_points = compact ? CompactCutCellVolumeLogicalPointCount(cell.compact_rule)
			: cell.rule.Points().size();
		return logical_points == 0 && cell.diagnostics.certified_reference_volume == 0.0
			&& cell.diagnostics.unresolved_reference_volume == 0.0;
	}
	static bool HasLogicalSupport(const CutCellVolumeQuadratureCell& cell, bool compact)
	{
		return compact ? CompactCutCellVolumeLogicalPointCount(cell.compact_rule) != 0 : !cell.rule.Points().empty();
	}
	void PublishAttempt(CutCellVolumeQuadratureCell& result, CutCellVolumeQuadratureCell&& attempt,
		const CutCellVolumeQuadratureDiagnostics& work, std::uint32_t effective_depth) const
	{
		result = std::move(attempt);
		result.diagnostics.retained_bytes = work.retained_bytes;
		result.diagnostics.rescue_attempts = work.rescue_attempts;
		result.diagnostics.attempted_nodes = work.attempted_nodes;
		result.diagnostics.attempted_leaves = work.attempted_leaves;
		result.diagnostics.attempted_output_points = work.attempted_output_points;
		result.diagnostics.attempted_samples = work.attempted_samples;
		result.diagnostics.attempted_record_attempts = work.attempted_record_attempts;
		result.diagnostics.attempted_logical_output_points = work.attempted_logical_output_points;
		result.diagnostics.rescue_effective_depth = effective_depth;
	}
	void BuildCellWithRescue(const CartesianDomainClassification& domain, const CartesianDomainCell& source,
		CutCellVolumeQuadratureCell& result) const
	{
		BuildWithRescue(domain, source, result, false);
	}
	void BuildCompactCellWithRescue(const CartesianDomainClassification& domain, const CartesianDomainCell& source,
		CutCellVolumeQuadratureCell& result) const
	{
		BuildWithRescue(domain, source, result, true);
	}
	void BuildWithRescue(const CartesianDomainClassification& domain, const CartesianDomainCell& source,
		CutCellVolumeQuadratureCell& result, bool compact) const
	{
		CutCellVolumeQuadratureDiagnostics work;
		for (std::uint32_t depth = options_.max_depth;; ++depth) {
			CutCellVolumeQuadratureCell attempt; attempt.id = source.id; attempt.classification = source.classification;
			const auto attempt_options = AttemptOptions(depth, work, compact);
			if (compact) BuildCompactCellAtDepth(domain, source, attempt, attempt_options);
			else BuildCellAtDepth(domain, source, attempt, attempt_options);
			AddAttemptWork(work, attempt.diagnostics, options_.max_retained_bytes);
			const bool positive = HasLogicalSupport(attempt, compact);
			const bool certified_empty = IsExactCertifiedEmpty(attempt, compact);
			const bool provisional_empty = !positive && !certified_empty;
			// A rescue retry cannot repair a midpoint that no longer subdivides.
			// Do not publish a usable, zero-support unresolved rule in this case:
			// it would conceal the precision failure in expanded storage while
			// compact consumer validation correctly rejects the same rule.
			if (provisional_empty && attempt.usable
				&& attempt.diagnostics.unresolved_reference_volume > 0.0
				&& options_.empty_rule_rescue_max_depth
				&& attempt.diagnostics.precision_limited_leaves != 0)
				throw std::runtime_error("cut-cell empty-rule rescue stopped by precision-limited subdivision for cell "
					+ std::to_string(source.id) + " at depth " + std::to_string(depth));
			if (positive || certified_empty || source.classification != CellClassification::Cut
				|| !options_.empty_rule_rescue_max_depth || source.ambiguous
				|| !attempt.usable || attempt.diagnostics.precision_limited_leaves != 0) {
				// With rescue disabled, retain a structurally valid zero-support Cut
				// attempt in either storage mode for backward-compatible construction.
				// Consumer-facing validation still rejects it as not certified empty.
				if (compact && attempt.usable) ValidateCompactStoredRule(attempt, nullptr,
					provisional_empty && !options_.empty_rule_rescue_max_depth);
				PublishAttempt(result, std::move(attempt), work, depth);
				return;
			}
			if (depth == options_.empty_rule_rescue_max_depth) {
				throw std::runtime_error("cut-cell empty-rule rescue exhausted for cell " + std::to_string(source.id)
					+ " at depth " + std::to_string(depth) + " with unresolved reference volume "
					+ std::to_string(attempt.diagnostics.unresolved_reference_volume));
			}
		}
	}

	void BuildCellAtDepth(const CartesianDomainClassification& domain, const CartesianDomainCell& source,
		CutCellVolumeQuadratureCell& result, const OctreeCutQuadratureOptions& options) const
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
			CheckAdd(result.diagnostics.nodes, 1, options.max_nodes, "octree cut quadrature node cap reached");
			CheckAdd(result.diagnostics.leaves, 1, options.max_leaves, "octree cut quadrature leaf cap reached");
			CheckAdd(result.diagnostics.output_points, result.rule.Points().size(), options.max_points, "octree cut quadrature point cap reached");
			CheckAdd(result.diagnostics.samples, result.rule.Points().size(), std::numeric_limits<std::size_t>::max(), "cut quadrature sample count overflows");
			result.diagnostics.logical_output_points = result.diagnostics.output_points;
			ValidateVolumeQuadratureRule(element, result.rule);
			ValidateReferenceEstimate(result.diagnostics.estimated_reference_volume);
			ValidateStoredRule(result);
			ValidateDiagnostics(result.diagnostics);
			return;
		}
		if (source.classification == CellClassification::Outside) {
			CheckAdd(result.diagnostics.nodes, 1, options.max_nodes, "octree cut quadrature node cap reached");
			CheckAdd(result.diagnostics.leaves, 1, options.max_leaves, "octree cut quadrature leaf cap reached");
			ValidateReferenceEstimate(result.diagnostics.estimated_reference_volume);
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
			CheckAdd(result.diagnostics.nodes, 1, options.max_nodes, "octree cut quadrature node cap reached");
			result.diagnostics.reached_depth = std::max(result.diagnostics.reached_depth, node.depth);
			bool ambiguous = false;
			const auto contact = domain.SurfaceIndex().IntersectBox(node.physical_bounds);
			if (contact == BoxContact::Ambiguous) ambiguous = true;
			else if (contact == BoxContact::None || contact == BoxContact::BoundaryOnly) {
				std::array<double, 3> center{};
				for (std::size_t axis = 0; axis < 3; ++axis) center[axis] = node.reference_lower[axis]+(node.reference_upper[axis]-node.reference_lower[axis])/2.0;
				const auto location = domain.SurfaceIndex().LocatePoint(EvaluateElementGeometry(element, center).physical);
				if (location == PointLocation::Inside) {
					CheckAdd(result.diagnostics.leaves, 1, options.max_leaves, "octree cut quadrature leaf cap reached");
					AddCertified(node, result.diagnostics, physical_volume);
					AppendScaledGauss(node, points, result.diagnostics, options);
					continue;
				}
				if (location == PointLocation::Outside) {
					CheckAdd(result.diagnostics.leaves, 1, options.max_leaves, "octree cut quadrature leaf cap reached");
					continue;
				}
				ambiguous = location == PointLocation::Ambiguous;
			}
			if (ambiguous) { result.usable = false; ++result.diagnostics.predicate_ambiguities; }
			std::array<double, 3> ref_mid{}, physical_mid{};
			bool midpoint_ok = node.depth < options.max_depth;
			for (std::size_t axis = 0; axis < 3 && midpoint_ok; ++axis)
				midpoint_ok = Midpoint(node.reference_lower[axis], node.reference_upper[axis], ref_mid[axis])
					&& Midpoint(node.physical_bounds.minimum[axis], node.physical_bounds.maximum[axis], physical_mid[axis]);
			if (!midpoint_ok) {
				CheckAdd(result.diagnostics.leaves, 1, options.max_leaves, "octree cut quadrature leaf cap reached");
				if (node.depth < options.max_depth) ++result.diagnostics.precision_limited_leaves;
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
								CheckAdd(result.diagnostics.output_points, 1, options.max_points, "octree cut quadrature point cap reached");
								const double weight = Volume(node.reference_lower, node.reference_upper)*kGaussFourWeights[qx]*kGaussFourWeights[qy]*kGaussFourWeights[qz]/8.0;
								points.push_back({reference, weight});
							} else if (location == PointLocation::Boundary) ++result.diagnostics.boundary_samples;
							else if (location == PointLocation::Ambiguous) { ++result.diagnostics.ambiguous_samples; ++result.diagnostics.predicate_ambiguities; result.usable = false; }
						}
				continue;
			}
			const std::size_t available_nodes = options.max_nodes-result.diagnostics.nodes;
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
		double reference_estimate_compensation = 0.0;
		for (const auto& point : result.rule.Points())
			AddCompensatedNonnegative(result.diagnostics.estimated_reference_volume, reference_estimate_compensation,
				point.weight, "cut quadrature weight sum overflows");
		ValidateReferenceEstimate(result.diagnostics.estimated_reference_volume);
		result.diagnostics.estimated_physical_volume = PhysicalRuleVolume(element, result.rule);
		result.diagnostics.logical_output_points = result.diagnostics.output_points;
		result.diagnostics.lower_reference_volume = result.diagnostics.certified_reference_volume;
		result.diagnostics.upper_reference_volume = result.diagnostics.certified_reference_volume+result.diagnostics.unresolved_reference_volume;
		result.diagnostics.lower_physical_volume = result.diagnostics.certified_physical_volume;
		result.diagnostics.upper_physical_volume = result.diagnostics.certified_physical_volume+result.diagnostics.unresolved_physical_volume;
		if (result.usable && !result.rule.Points().empty()) ValidateVolumeQuadratureRule(element, result.rule);
		ValidateStoredRule(result);
		ValidateDiagnostics(result.diagnostics);
	}

	static CompactCutCellVolumeBlock CompactBlock(const Node& node, std::uint32_t max_depth)
	{
		if (node.depth > max_depth) throw std::runtime_error("compact cut-cell depth is invalid");
		const std::uint32_t width = std::uint32_t(1) << (max_depth-node.depth);
		CompactCutCellVolumeBlock result;
		for (std::size_t axis = 0; axis < 3; ++axis) {
			result.lower[axis] = node.key[axis]*width;
			result.upper[axis] = result.lower[axis]+width;
		}
		return result;
	}
	static bool CompactBlockLess(const CompactCutCellVolumeBlock& left, const CompactCutCellVolumeBlock& right)
	{
		for (std::size_t axis = 0; axis < 3; ++axis) if (left.lower[axis] != right.lower[axis]) return left.lower[axis] < right.lower[axis];
		for (std::size_t axis = 0; axis < 3; ++axis) if (left.upper[axis] != right.upper[axis]) return left.upper[axis] < right.upper[axis];
		return false;
	}
	static void CoalesceCompactBlocks(std::vector<CompactCutCellVolumeBlock>& blocks)
	{
		bool changed = true;
		while (changed) {
			changed = false;
			for (std::size_t axis = 0; axis < 3; ++axis) {
				std::sort(blocks.begin(), blocks.end(), [axis](const CompactCutCellVolumeBlock& a, const CompactCutCellVolumeBlock& b) {
					for (std::size_t other = 0; other < 3; ++other) if (other != axis) {
						if (a.lower[other] != b.lower[other]) return a.lower[other] < b.lower[other];
						if (a.upper[other] != b.upper[other]) return a.upper[other] < b.upper[other];
					}
					if (a.lower[axis] != b.lower[axis]) return a.lower[axis] < b.lower[axis];
					return a.upper[axis] < b.upper[axis];
				});
				std::size_t retained = 0;
				for (const auto& block : blocks) {
					if (retained) {
						auto& previous = blocks[retained-1]; bool same = previous.upper[axis] == block.lower[axis];
						for (std::size_t other = 0; other < 3; ++other) if (other != axis)
							same = same && previous.lower[other] == block.lower[other] && previous.upper[other] == block.upper[other];
						if (same) { previous.upper[axis] = block.upper[axis]; changed = true; continue; }
					}
					blocks[retained++] = block;
				}
				blocks.resize(retained);
			}
		}
		std::sort(blocks.begin(), blocks.end(), CompactBlockLess);
	}

	void BuildCompactCellAtDepth(const CartesianDomainClassification& domain, const CartesianDomainCell& source,
		CutCellVolumeQuadratureCell& result, const OctreeCutQuadratureOptions& options) const
	{
		const auto cell = domain.Background().Cell(source.id);
		const Element element = domain.Background().MaterializeElement(source.id);
		const double physical_volume = Volume(cell.lower_m, cell.upper_m);
		result.compact_rule.max_depth = options.max_depth;
		double reference_estimate_compensation = 0.0, physical_estimate_compensation = 0.0;
		double certified_reference_compensation = 0.0, certified_physical_compensation = 0.0;
		double unresolved_reference_compensation = 0.0, unresolved_physical_compensation = 0.0;
		const auto add_certified = [&](const Node& node) {
			const double volume = Volume(node.reference_lower, node.reference_upper);
			AddCompensatedNonnegative(result.diagnostics.certified_reference_volume, certified_reference_compensation, volume, "cut quadrature certified volume overflows");
			AddCompensatedNonnegative(result.diagnostics.certified_physical_volume, certified_physical_compensation, volume*physical_volume, "cut quadrature physical volume overflows");
		};
		const auto add_unresolved = [&](const Node& node) {
			const double volume = Volume(node.reference_lower, node.reference_upper);
			AddCompensatedNonnegative(result.diagnostics.unresolved_reference_volume, unresolved_reference_compensation, volume, "cut quadrature unresolved volume overflows");
			AddCompensatedNonnegative(result.diagnostics.unresolved_physical_volume, unresolved_physical_compensation, volume*physical_volume, "cut quadrature physical volume overflows");
		};
		const auto add_node_estimate = [&](const Node& node, std::uint64_t mask) {
			const double volume = Volume(node.reference_lower, node.reference_upper);
			for (std::size_t qz = 0; qz < 4; ++qz) for (std::size_t qy = 0; qy < 4; ++qy) for (std::size_t qx = 0; qx < 4; ++qx) {
				const std::size_t bit = qx+4*qy+16*qz;
				if ((mask & (std::uint64_t(1) << bit)) == 0) continue;
				std::array<double, 3> reference{};
				for (std::size_t axis = 0; axis < 3; ++axis) { const std::size_t q = axis == 0 ? qx : (axis == 1 ? qy : qz); reference[axis] = node.reference_lower[axis]+(node.reference_upper[axis]-node.reference_lower[axis])*kGaussFourPoints[q]; }
				const double weight = volume*kGaussFourWeights[qx]*kGaussFourWeights[qy]*kGaussFourWeights[qz]/8.0;
				AddCompensatedNonnegative(result.diagnostics.estimated_reference_volume, reference_estimate_compensation, weight, "cut quadrature weight sum overflows");
				const auto geometry = EvaluateElementGeometry(element, reference);
				if (!std::isfinite(geometry.raw_determinant) || !(geometry.raw_determinant > 0.0)
					|| weight > std::numeric_limits<double>::max()/geometry.raw_determinant)
					throw std::runtime_error("cut-cell quadrature has invalid physical Jacobian");
				AddCompensatedNonnegative(result.diagnostics.estimated_physical_volume, physical_estimate_compensation, weight*geometry.raw_determinant, "cut quadrature physical estimate overflows");
			}
		};
		// Scratch vectors are shared by the whole recursive walk.  Record attempts
		// are deliberately monotone: collapsing an all-full/all-empty subtree
		// releases storage but cannot evade the compact work caps.
		std::vector<CompactCutCellVolumeBlock> blocks;
		std::vector<CompactCutCellVolumeSampleLeaf> samples;
		std::size_t record_attempts = 0;
		std::size_t block_accounted_capacity = 0, sample_accounted_capacity = 0;
		std::size_t peak_retained_bytes = 0;
		const auto checked_bytes = [&](std::size_t block_count, std::size_t sample_count) {
			if (block_count > std::numeric_limits<std::size_t>::max()/sizeof(CompactCutCellVolumeBlock)
				|| sample_count > std::numeric_limits<std::size_t>::max()/sizeof(CompactCutCellVolumeSampleLeaf))
				throw std::overflow_error("compact cut-cell retained byte count overflows");
			const std::size_t block_bytes = block_count*sizeof(CompactCutCellVolumeBlock);
			const std::size_t sample_bytes = sample_count*sizeof(CompactCutCellVolumeSampleLeaf);
			if (sample_bytes > std::numeric_limits<std::size_t>::max()-block_bytes)
				throw std::overflow_error("compact cut-cell retained byte count overflows");
			return block_bytes+sample_bytes;
		};
		const auto checked_count_sum = [](std::size_t left, std::size_t right) {
			if (right > std::numeric_limits<std::size_t>::max()-left)
				throw std::overflow_error("compact cut-cell retained byte count overflows");
			return left+right;
		};
		const auto checked_peak_bytes = [&](std::size_t old_block_capacity, std::size_t new_block_capacity,
			std::size_t old_sample_capacity, std::size_t new_sample_capacity) {
			return checked_bytes(checked_count_sum(old_block_capacity, new_block_capacity),
				checked_count_sum(old_sample_capacity, new_sample_capacity));
		};
		const auto grow_capacity = [](std::size_t capacity, std::size_t required) {
			if (capacity >= required) return capacity;
			const std::size_t increment = capacity/2+capacity%2;
			if (increment > std::numeric_limits<std::size_t>::max()-capacity)
				return std::numeric_limits<std::size_t>::max();
			return std::max(required, capacity+increment);
		};
		const auto observe_accounted_bytes = [&](std::size_t bytes) {
			if (bytes > options.max_retained_bytes)
				throw std::runtime_error("compact cut-cell retained byte cap reached");
			peak_retained_bytes = std::max(peak_retained_bytes, bytes);
		};
		const auto reserve_block = [&] {
			if (blocks.size() == block_accounted_capacity) {
				const std::size_t required = checked_count_sum(blocks.size(), 1);
				const std::size_t next_capacity = grow_capacity(block_accounted_capacity, required);
				observe_accounted_bytes(checked_peak_bytes(block_accounted_capacity, next_capacity, sample_accounted_capacity, 0));
				blocks.reserve(next_capacity);
				block_accounted_capacity = next_capacity;
			}
		};
		const auto reserve_sample = [&] {
			if (samples.size() == sample_accounted_capacity) {
				const std::size_t required = checked_count_sum(samples.size(), 1);
				const std::size_t next_capacity = grow_capacity(sample_accounted_capacity, required);
				observe_accounted_bytes(checked_peak_bytes(block_accounted_capacity, 0, sample_accounted_capacity, next_capacity));
				samples.reserve(next_capacity);
				sample_accounted_capacity = next_capacity;
			}
		};
		const auto append_block = [&](const CompactCutCellVolumeBlock& block) {
			CheckAdd(record_attempts, 1, options.max_records, "compact cut-cell record cap reached");
			reserve_block(); blocks.push_back(block);
		};
		const auto append_sample = [&](const CompactCutCellVolumeSampleLeaf& sample) {
			CheckAdd(record_attempts, 1, options.max_records, "compact cut-cell record cap reached");
			reserve_sample(); samples.push_back(sample);
		};
		const auto finish = [&] {
			CoalesceCompactBlocks(blocks);
			const std::size_t records = checked_count_sum(blocks.size(), samples.size());
			if (records > options.max_records) throw std::runtime_error("compact cut-cell record cap reached");
			CompactCutCellVolumeRule published; published.max_depth = options.max_depth;
			published.certified_blocks = std::move(blocks); published.sample_leaves = std::move(samples);
			const std::size_t logical_points = CompactCutCellVolumeLogicalPointCount(published);
			if (logical_points > options.max_logical_points)
				throw std::runtime_error("compact cut-cell logical point cap reached");
			result.compact_rule = std::move(published);
			result.diagnostics.certified_blocks = result.compact_rule.certified_blocks.size();
			result.diagnostics.sample_leaves = result.compact_rule.sample_leaves.size();
			result.diagnostics.record_attempts = record_attempts;
			result.diagnostics.retained_bytes = peak_retained_bytes;
			result.diagnostics.observed_retained_bytes = checked_bytes(result.compact_rule.certified_blocks.capacity(),
				result.compact_rule.sample_leaves.capacity());
			result.diagnostics.logical_output_points = logical_points;
			result.diagnostics.lower_reference_volume = result.diagnostics.certified_reference_volume;
			result.diagnostics.upper_reference_volume = result.diagnostics.certified_reference_volume+result.diagnostics.unresolved_reference_volume;
			result.diagnostics.lower_physical_volume = result.diagnostics.certified_physical_volume;
			result.diagnostics.upper_physical_volume = result.diagnostics.certified_physical_volume+result.diagnostics.unresolved_physical_volume;
			ValidateReferenceEstimate(result.diagnostics.estimated_reference_volume);
			// Empty unresolved attempts are structurally valid provisional rescue
			// candidates.  BuildWithRescue applies the publishability predicate.
			ValidateCompactStoredRule(result, &element, true); ValidateDiagnostics(result.diagnostics);
		};
		if (source.classification == CellClassification::Inside) {
			CheckAdd(result.diagnostics.nodes, 1, options.max_nodes, "octree cut quadrature node cap reached");
			CheckAdd(result.diagnostics.leaves, 1, options.max_leaves, "octree cut quadrature leaf cap reached");
			CheckAdd(result.diagnostics.samples, 64, std::numeric_limits<std::size_t>::max(), "cut quadrature sample count overflows");
			append_block(CompactBlock(Node{{{0,0,0}},{{1,1,1}},source.bounds,0,{{0,0,0}}}, options.max_depth));
			result.diagnostics.certified_reference_volume = 1.0; result.diagnostics.certified_physical_volume = physical_volume;
			add_node_estimate(Node{{{0,0,0}},{{1,1,1}},source.bounds,0,{{0,0,0}}}, ~std::uint64_t(0));
			finish(); return;
		}
		if (source.classification == CellClassification::Outside) {
			CheckAdd(result.diagnostics.nodes, 1, options.max_nodes, "octree cut quadrature node cap reached");
			CheckAdd(result.diagnostics.leaves, 1, options.max_leaves, "octree cut quadrature leaf cap reached");
			finish(); return;
		}
		if (source.ambiguous) { result.usable = false; ++result.diagnostics.predicate_ambiguities; }
		const auto recurse = [&](auto&& self, const Node& node) -> CompactNodeState {
			CheckAdd(result.diagnostics.nodes, 1, options.max_nodes, "octree cut quadrature node cap reached");
			result.diagnostics.reached_depth = std::max(result.diagnostics.reached_depth, node.depth);
			bool ambiguous = false;
			const auto contact = domain.SurfaceIndex().IntersectBox(node.physical_bounds);
			if (contact == BoxContact::Ambiguous) ambiguous = true;
			else if (contact == BoxContact::None || contact == BoxContact::BoundaryOnly) {
				std::array<double, 3> center{};
				for (std::size_t axis = 0; axis < 3; ++axis) center[axis] = node.reference_lower[axis]+(node.reference_upper[axis]-node.reference_lower[axis])/2.0;
				const auto location = domain.SurfaceIndex().LocatePoint(EvaluateElementGeometry(element, center).physical);
				if (location == PointLocation::Inside) {
					CheckAdd(result.diagnostics.leaves, 1, options.max_leaves, "octree cut quadrature leaf cap reached");
					CheckAdd(result.diagnostics.samples, 64, std::numeric_limits<std::size_t>::max(), "cut quadrature sample count overflows");
					add_certified(node); add_node_estimate(node, ~std::uint64_t(0)); return CompactNodeState::UniformFull;
				}
				if (location == PointLocation::Outside) { CheckAdd(result.diagnostics.leaves, 1, options.max_leaves, "octree cut quadrature leaf cap reached"); return CompactNodeState::UniformEmpty; }
				ambiguous = location == PointLocation::Ambiguous;
			}
			if (ambiguous) { result.usable = false; ++result.diagnostics.predicate_ambiguities; }
			std::array<double, 3> ref_mid{}, physical_mid{};
			bool midpoint_ok = node.depth < options.max_depth;
			for (std::size_t axis = 0; axis < 3 && midpoint_ok; ++axis)
				midpoint_ok = Midpoint(node.reference_lower[axis], node.reference_upper[axis], ref_mid[axis])
					&& Midpoint(node.physical_bounds.minimum[axis], node.physical_bounds.maximum[axis], physical_mid[axis]);
			if (!midpoint_ok) {
				CheckAdd(result.diagnostics.leaves, 1, options.max_leaves, "octree cut quadrature leaf cap reached");
				if (node.depth < options.max_depth) ++result.diagnostics.precision_limited_leaves;
				add_unresolved(node);
				CompactCutCellVolumeSampleLeaf leaf; leaf.key = node.key; leaf.depth = node.depth;
				for (std::size_t qz = 0; qz < 4; ++qz) for (std::size_t qy = 0; qy < 4; ++qy) for (std::size_t qx = 0; qx < 4; ++qx) {
					CheckAdd(result.diagnostics.samples, 1, std::numeric_limits<std::size_t>::max(), "cut quadrature sample count overflows");
					std::array<double, 3> reference{};
					for (std::size_t axis = 0; axis < 3; ++axis) { const std::size_t q = axis == 0 ? qx : (axis == 1 ? qy : qz); reference[axis] = node.reference_lower[axis]+(node.reference_upper[axis]-node.reference_lower[axis])*kGaussFourPoints[q]; }
					const auto location = domain.SurfaceIndex().LocatePoint(EvaluateElementGeometry(element, reference).physical);
					if (location == PointLocation::Inside) leaf.inside_mask |= std::uint64_t(1) << (qx+4*qy+16*qz);
					else if (location == PointLocation::Boundary) ++result.diagnostics.boundary_samples;
					else if (location == PointLocation::Ambiguous) { ++result.diagnostics.ambiguous_samples; ++result.diagnostics.predicate_ambiguities; result.usable = false; }
				}
				add_node_estimate(node, leaf.inside_mask);
				append_sample(leaf); return CompactNodeState::Mixed;
			}
			std::array<Node, 8> children{};
			for (int child = 0; child < 8; ++child) {
				auto& next = children[static_cast<std::size_t>(child)]; next.depth = node.depth+1;
				for (std::size_t axis = 0; axis < 3; ++axis) { const bool high = (child & (1 << axis)) != 0;
					next.reference_lower[axis] = high ? ref_mid[axis] : node.reference_lower[axis]; next.reference_upper[axis] = high ? node.reference_upper[axis] : ref_mid[axis];
					next.physical_bounds.minimum[axis] = high ? physical_mid[axis] : node.physical_bounds.minimum[axis]; next.physical_bounds.maximum[axis] = high ? node.physical_bounds.maximum[axis] : physical_mid[axis];
					next.key[axis] = 2*node.key[axis]+(high ? 1u : 0u); }
			}
			const std::size_t block_start = blocks.size(), sample_start = samples.size();
			bool all_full = true, all_empty = true;
			for (std::size_t child = 0; child < 8; ++child) {
				const CompactNodeState state = self(self, children[child]);
				all_full = all_full && state == CompactNodeState::UniformFull;
				all_empty = all_empty && state == CompactNodeState::UniformEmpty;
				if (state == CompactNodeState::UniformFull) append_block(CompactBlock(children[child], options.max_depth));
			}
			if (all_full || all_empty) {
				RollbackCompactRecords(blocks, samples, block_start, sample_start, result.diagnostics, options.max_records);
				return all_full ? CompactNodeState::UniformFull : CompactNodeState::UniformEmpty;
			}
			return CompactNodeState::Mixed;
		};
		Node root_node{{{0,0,0}}, {{1,1,1}}, source.bounds, 0, {{0,0,0}}};
		const CompactNodeState compact = recurse(recurse, root_node);
		if (compact == CompactNodeState::UniformFull) append_block(CompactBlock(root_node, options.max_depth));
		finish();
	}

	CubicCartesianGridSpec grid_spec_{};
	std::string surface_canonical_hash_;
	OctreeCutQuadratureOptions options_{};
	CutCellVolumeQuadratureStorageMode storage_mode_ = CutCellVolumeQuadratureStorageMode::Expanded;
	std::vector<CutCellVolumeQuadratureCell> cells_;
	CutCellVolumeQuadratureDiagnostics diagnostics_;
};

} // namespace iga

#endif
