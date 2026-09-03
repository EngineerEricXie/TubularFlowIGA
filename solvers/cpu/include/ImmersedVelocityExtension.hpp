#ifndef IGA_IMMERSED_VELOCITY_EXTENSION_HPP
#define IGA_IMMERSED_VELOCITY_EXTENSION_HPP

// Deterministic, standalone continuation of a committed immersed velocity to
// newly-active Cartesian spline nodes.  This intentionally has no runtime or
// PETSc ownership: callers either receive a complete immutable result or an
// exception before anything can be published.
#include "CubicCartesianSplineFace.hpp"
#include "ImmersedTransientState.hpp"
#include "MovingCutGeometry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace iga {

struct ImmersedVelocityExtensionOptions {
	std::size_t max_unknowns = 2048;
	std::size_t max_dense_bytes = 64u*1024u*1024u;
	std::size_t max_factor_work = 3000000000ull;
	std::size_t max_band_cells = 250000;
	std::size_t max_faces = 750000;
	std::size_t max_nodes = 250000;
	std::size_t max_trace_entries = 20000000;
	std::size_t max_reverse_visited_cells = 250000;
	double residual_tolerance = 1.e-11;
};

struct ImmersedVelocityExtensionDiagnostics {
	std::size_t old_positive_cells = 0, new_positive_cells = 0, band_cells = 0;
	std::size_t old_nodes = 0, band_nodes = 0, anchors = 0, unknowns = 0, faces = 0, trace_entries = 0;
	std::uint32_t maximum_forward_distance = 0, maximum_reverse_distance = 0;
	double minimum_diagonal = 0.0, maximum_diagonal = 0.0;
	double minimum_pivot = 0.0, maximum_pivot = 0.0, pivot_threshold = 0.0;
	std::array<double, 3> residual{{0.0,0.0,0.0}};
};

struct ImmersedScalarExtension {
	std::vector<std::int32_t> band_node_ids, target_node_ids;
	std::vector<double> band_values, target_values;
	std::string hash_sha256;
};

class ImmersedVelocityExtension {
public:
	ImmersedVelocityExtension(const ImmersedVelocityExtension&) = delete;
	ImmersedVelocityExtension& operator=(const ImmersedVelocityExtension&) = delete;
	ImmersedVelocityExtension(ImmersedVelocityExtension&&) noexcept = default;
	ImmersedVelocityExtension& operator=(ImmersedVelocityExtension&&) noexcept = default;

	static ImmersedVelocityExtension Build(const MovingCutGeometry& old_geometry,
		const ImmersedActiveLayout& old_layout, const ImmersedGlobalFlowState& old_state,
		const MovingCutGeometry& new_geometry, const ImmersedActiveLayout& new_layout,
		std::uint32_t layers, ImmersedVelocityExtensionOptions options = {})
	{
		ImmersedVelocityExtension result;
		result.options_ = ValidateOptions(options); result.layers_ = layers;
#ifdef IGA_IMMERSED_VELOCITY_EXTENSION_TESTING
		TestingResetCapCounters();
#endif
		result.old_geometry_identity_ = old_geometry.GeometryIdentitySha256();
		result.new_geometry_identity_ = new_geometry.GeometryIdentitySha256();
		result.old_layout_identity_ = old_layout.HashSha256(); result.new_layout_identity_ = new_layout.HashSha256();
		result.old_state_identity_ = old_state.HashSha256();
		result.old_time_s_ = old_state.TimeS(); result.new_time_s_ = new_geometry.Evaluation().EvaluatedTimeS();
		if (!old_layout.Valid() || !new_layout.Valid() || !old_state.Valid()) throw std::invalid_argument("immersed velocity extension received invalid state or layout");
		if (old_layout.GeometryIdentity() != result.old_geometry_identity_ || new_layout.GeometryIdentity() != result.new_geometry_identity_
			|| old_state.GeometryIdentity() != old_layout.GeometryIdentity() || old_state.NodeIds() != old_layout.NodeIds()
			|| old_state.PortIds() != old_layout.PortIds() || old_state.HasGaugeMultiplier() != old_layout.HasGaugeRow()
			|| old_state.TimeS() != old_geometry.Evaluation().EvaluatedTimeS())
			throw std::invalid_argument("immersed velocity extension geometry, layout, and state identities do not match");
		if (!(result.new_time_s_ > result.old_time_s_) || !std::isfinite(result.new_time_s_))
			throw std::invalid_argument("immersed velocity extension requires strictly increasing geometry time");
		if (!SameGrid(old_geometry.Domain().Background().Spec(), new_geometry.Domain().Background().Spec()))
			throw std::invalid_argument("immersed velocity extension requires one fixed Cartesian grid");
		result.background_ = old_geometry.Domain().Background(); result.grid_ = result.background_.Spec();
		auto old_cells = PositiveCells(old_geometry), new_cells = PositiveCells(new_geometry);
#ifdef IGA_IMMERSED_VELOCITY_EXTENSION_TESTING
		if (TestingFaultActive(TestingFault::EmptyOldCells)) old_cells.clear();
		if (TestingFaultActive(TestingFault::EmptyNewCells)) new_cells.clear();
#endif
		if (old_cells.empty() || new_cells.empty()) throw std::invalid_argument("immersed velocity extension positive-cell set is empty");
		result.old_cells_ = old_cells; result.new_cells_ = new_cells;
		if (old_cells.size() > options.max_band_cells) throw std::length_error("immersed velocity extension old positive-cell set exceeds band cap");
		if (new_cells.size() > options.max_band_cells) throw std::length_error("immersed velocity extension new positive-cell set exceeds band cap");
		RequireLayoutNodes(old_geometry, old_cells, old_layout, options.max_nodes, "old");
		RequireLayoutNodes(new_geometry, new_cells, new_layout, options.max_nodes, "new");
		result.BuildBand(old_geometry.Domain().Background(), old_cells);
		if (result.band_cells_.size() > options.max_band_cells) throw std::length_error("immersed velocity extension band cell cap exceeded");
		result.CheckNewCells();
		result.BuildNodes(old_layout, old_state);
		result.BuildFaces(old_geometry.Domain().Background());
		result.PreflightTraceGraph();
		result.BuildMatrixAndFactor();
		result.ExtendVelocity(old_state);
		result.BuildHistory(new_layout);
		result.operator_hash_ = result.HashOperator(); result.reduced_hash_ = result.HashReduced();
		result.extension_hash_ = result.HashExtension();
		return result;
	}

	const std::vector<std::uint64_t>& BandCells() const noexcept { return band_cells_; }
	const std::vector<std::uint32_t>& BandDistances() const noexcept { return band_distances_; }
	const std::vector<std::int32_t>& BandNodeIds() const noexcept { return band_nodes_; }
	const std::vector<std::int32_t>& AnchorNodeIds() const noexcept { return anchor_nodes_; }
	const std::vector<std::int32_t>& UnknownNodeIds() const noexcept { return unknown_nodes_; }
	const ImmersedVelocityHistory& TargetHistory() const noexcept { return history_; }
	const ImmersedVelocityExtensionDiagnostics& Diagnostics() const noexcept { return diagnostics_; }
	const std::string& OperatorHashSha256() const noexcept { return operator_hash_; }
	const std::string& ReducedHashSha256() const noexcept { return reduced_hash_; }
	const std::string& HashSha256() const noexcept { return extension_hash_; }

#ifdef IGA_IMMERSED_VELOCITY_EXTENSION_TESTING
	// These read-only seams and deliberate gate probes are compiled only into the
	// focused evidence executable.  They neither alter construction nor expose
	// mutable storage to a production caller.
	struct TestingFaceTrace {
		std::uint64_t minus_cell = 0, plus_cell = 0;
		std::uint8_t axis = 0;
		double h = 0.0, area = 0.0;
		std::vector<std::int32_t> nodes;
		std::vector<double> alpha;
		std::vector<std::vector<double>> jumps;
	};
	const std::vector<double>& TestingUpper() const noexcept { return upper_; }
	const std::vector<double>& TestingDiagonal() const noexcept { return diagonal_; }
	bool TestingFactorized() const noexcept { return !factor_.empty(); }
	std::size_t TestingPivotCount() const noexcept { return factor_.empty() ? 0u : diagonal_.size(); }
	const std::vector<std::array<double,3>>& TestingBandVelocities() const noexcept { return band_velocities_; }
	const std::vector<ImmersedVelocityHistoryProvenance>& TestingBandProvenance() const noexcept { return band_provenance_; }
	const std::vector<std::uint64_t>& TestingReverseCells() const noexcept { return reverse_cells_; }
	const std::vector<std::uint32_t>& TestingReverseDistances() const noexcept { return reverse_distances_; }
	static std::size_t TestingForwardVisitedCount() noexcept { return TestingCapCountersState().forward_visited; }
	static std::size_t TestingForwardInsertionAttempts() noexcept { return TestingCapCountersState().forward_insertion_attempts; }
	static std::size_t TestingReverseVisitedCount() noexcept { return TestingCapCountersState().reverse_visited; }
	static std::size_t TestingReverseInsertionAttempts() noexcept { return TestingCapCountersState().reverse_insertion_attempts; }
	static std::size_t TestingLayoutNodeCount() noexcept { return TestingCapCountersState().layout_nodes; }
	static std::size_t TestingBandNodeCount() noexcept { return TestingCapCountersState().band_nodes; }
	static std::size_t TestingUnknownNodeCount() noexcept { return TestingCapCountersState().unknown_nodes; }
	std::vector<double> TestingRhs(const std::vector<double>& full) const { return BuildRhs(full); }
	double TestingScalarResidual(const std::vector<double>& anchor_values_old_layout_order) const
	{
		if (anchor_values_old_layout_order.size() != old_layout_nodes_.size()) throw std::invalid_argument("immersed scalar extension does not match old layout order");
		std::vector<double> full(band_nodes_.size());
		for (std::size_t a=0; a<anchor_nodes_.size(); ++a) full[anchor_band_positions_[a]] = anchor_values_old_layout_order.at(old_layout_positions_[a]);
		const auto x = SolveComponent(full); for (std::size_t u=0; u<x.size(); ++u) full[unknown_band_positions_[u]]=x[u];
		return Residual(full, x, BuildRhs(full));
	}
	std::vector<TestingFaceTrace> TestingFaceTraces() const
	{
		std::vector<TestingFaceTrace> result; result.reserve(faces_.size());
		for (const auto& face : faces_)
			result.push_back({face.minus_cell, face.plus_cell, face.axis, face.h, face.area, face.nodes, face.alpha, face.jumps});
		return result;
	}
	static double TestingScaledCoefficient(double value, double diagonal_i, double diagonal_j)
	{ return ScaleCoefficient(value, std::sqrt(diagonal_i), std::sqrt(diagonal_j)); }
	enum class TestingFault : std::uint8_t { Diagonal, TraceAnchoring, Pivot, TargetCoverage, VelocityResidual, ScalarResidual, ReverseEmptyFrontier, EmptyOldCells, EmptyNewCells };
	class TestingFaultScope {
	public:
		explicit TestingFaultScope(TestingFault fault) : fault_(fault) { TestingFaults()[static_cast<std::size_t>(fault_)] = true; }
		~TestingFaultScope() { TestingFaults().fill(false); }
		TestingFaultScope(const TestingFaultScope&) = delete;
		TestingFaultScope& operator=(const TestingFaultScope&) = delete;
	private:
		TestingFault fault_;
	};
	static std::uint32_t TestingNarrowDistance(std::uint64_t distance)
	{ return NarrowDistance(distance, "immersed velocity extension test distance exceeds storage"); }
#endif

	ImmersedScalarExtension ExtendScalar(const std::vector<double>& anchor_values_old_layout_order) const
	{
		if (anchor_values_old_layout_order.size() != old_layout_nodes_.size()) throw std::invalid_argument("immersed scalar extension does not match old layout order");
		for (double v : anchor_values_old_layout_order) if (!std::isfinite(v)) throw std::invalid_argument("immersed scalar extension anchor is nonfinite");
		std::vector<double> full(band_nodes_.size());
		for (std::size_t a=0; a<anchor_nodes_.size(); ++a) full[anchor_band_positions_[a]] = anchor_values_old_layout_order.at(old_layout_positions_[a]);
		const auto x = SolveComponent(full); for (std::size_t u=0;u<x.size();++u) full[unknown_band_positions_[u]]=x[u];
		double residual = Residual(full, x, BuildRhs(full));
#ifdef IGA_IMMERSED_VELOCITY_EXTENSION_TESTING
		if (TestingFaultActive(TestingFault::ScalarResidual)) residual = std::numeric_limits<double>::quiet_NaN();
#endif
		if (!(std::isfinite(residual) && residual <= options_.residual_tolerance)) throw std::runtime_error("immersed scalar extension residual gate failed");
		ImmersedScalarExtension answer; answer.band_node_ids=band_nodes_; answer.band_values=std::move(full);
		answer.target_node_ids=target_node_ids_; answer.target_values.reserve(target_node_ids_.size());
		for (const auto id:target_node_ids_) answer.target_values.push_back(answer.band_values.at(BandPosition(id)));
		Sha256 h; String(h,"ImmersedScalarExtension/v1"); String(h,extension_hash_); Ids(h,answer.band_node_ids); for(double v:answer.band_values)h.AppendNormalizedDouble(v); Ids(h,answer.target_node_ids); for(double v:answer.target_values)h.AppendNormalizedDouble(v); answer.hash_sha256=h.Hex(); return answer;
	}

private:
	ImmersedVelocityExtension() = default;
	struct FaceTrace { std::uint64_t minus_cell=0, plus_cell=0; std::uint8_t axis=0; double h=0.0, area=0.0; std::vector<std::int32_t> nodes; std::vector<double> alpha; std::vector<std::vector<double>> jumps; };
#ifdef IGA_IMMERSED_VELOCITY_EXTENSION_TESTING
	struct TestingCapCounters {
		std::size_t forward_visited = 0, forward_insertion_attempts = 0;
		std::size_t reverse_visited = 0, reverse_insertion_attempts = 0;
		std::size_t layout_nodes = 0, band_nodes = 0, unknown_nodes = 0;
	};
	static TestingCapCounters& TestingCapCountersState()
	{
		static TestingCapCounters counters;
		return counters;
	}
	static void TestingResetCapCounters() { TestingCapCountersState() = {}; }
	static std::array<bool,9>& TestingFaults()
	{
		static std::array<bool,9> faults{{false,false,false,false,false,false,false,false,false}};
		return faults;
	}
	static bool TestingFaultActive(TestingFault fault)
	{ return TestingFaults()[static_cast<std::size_t>(fault)]; }
#endif
	static void String(Sha256& h,const std::string& v){ h.AppendLittleEndian64(v.size());h.Append(v.data(),v.size()); }
	static void Ids(Sha256& h,const std::vector<std::int32_t>& v){h.AppendLittleEndian64(v.size());for(auto x:v)h.AppendLittleEndian32(static_cast<std::uint32_t>(x));}
	static bool SameGrid(const CubicCartesianGridSpec&a,const CubicCartesianGridSpec&b){return a.lower_m==b.lower_m&&a.upper_m==b.upper_m&&a.cells==b.cells;}
	static std::size_t Mul(std::size_t a,std::size_t b,const char* m){if(a&&b>std::numeric_limits<std::size_t>::max()/a)throw std::overflow_error(m);return a*b;}
	static std::size_t Add(std::size_t a,std::size_t b,const char* m){if(b>std::numeric_limits<std::size_t>::max()-a)throw std::overflow_error(m);return a+b;}
	static double FiniteDouble(long double value,const char* message)
	{ if(!std::isfinite(value)||value>static_cast<long double>(std::numeric_limits<double>::max())||value<-static_cast<long double>(std::numeric_limits<double>::max()))throw std::runtime_error(message); return static_cast<double>(value); }
	static double AddProduct(double current,double left,double middle,double right,const char* message)
	{ if(!std::isfinite(current)||!std::isfinite(left)||!std::isfinite(middle)||!std::isfinite(right))throw std::runtime_error(message);return FiniteDouble(static_cast<long double>(current)+static_cast<long double>(left)*static_cast<long double>(middle)*static_cast<long double>(right),message); }
	static double ScaleCoefficient(double value,double diagonal_i,double diagonal_j)
	{ if(!std::isfinite(value)||!std::isfinite(diagonal_i)||!std::isfinite(diagonal_j)||!(diagonal_i>0.0)||!(diagonal_j>0.0))throw std::runtime_error("immersed velocity extension Jacobi scale is invalid"); const long double scaled=static_cast<long double>(value)/static_cast<long double>(diagonal_i)/static_cast<long double>(diagonal_j); return FiniteDouble(scaled,"immersed velocity extension scaled coefficient is not representable"); }
	static double ScaleValue(double value,double square_root_diagonal,const char* message)
	{ if(!std::isfinite(value)||!std::isfinite(square_root_diagonal)||!(square_root_diagonal>0.0))throw std::runtime_error(message);return FiniteDouble(static_cast<long double>(value)/square_root_diagonal,message); }
	static std::uint32_t NarrowDistance(std::uint64_t distance, const char* message)
	{
		if (distance > static_cast<std::uint64_t>(std::numeric_limits<std::uint32_t>::max())) throw std::overflow_error(message);
		return static_cast<std::uint32_t>(distance);
	}
	static ImmersedVelocityExtensionOptions ValidateOptions(ImmersedVelocityExtensionOptions o){if(!o.max_unknowns||!o.max_dense_bytes||!o.max_factor_work||!o.max_band_cells||!o.max_faces||!o.max_nodes||!o.max_trace_entries||!o.max_reverse_visited_cells||!std::isfinite(o.residual_tolerance)||!(o.residual_tolerance>0.0))throw std::invalid_argument("immersed velocity extension options are invalid");return o;}
	static std::vector<std::uint64_t> PositiveCells(const MovingCutGeometry& g)
	{
		const auto& d=g.Domain();const auto& v=g.Volume(); if(d.Cells().size()!=v.Cells().size()||d.Cells().size()!=d.Background().ElementCount()||!SameGrid(d.Background().Spec(),v.GridSpec())||d.SurfaceCanonicalHash()!=v.SurfaceCanonicalHash())throw std::invalid_argument("immersed velocity extension geometry catalog is incomplete or mismatched");
		std::vector<std::uint64_t> r; for(std::size_t i=0;i<d.Cells().size();++i){const auto& c=d.Cells()[i];const auto&q=v.Cells()[i];if(c.id!=i||q.id!=i||c.classification!=q.classification||c.ambiguous||!q.usable)throw std::invalid_argument("immersed velocity extension catalog cell is ambiguous or unusable");const double f=q.diagnostics.estimated_reference_volume;if(!std::isfinite(f)||f<0.0||f>1.0)throw std::invalid_argument("immersed velocity extension volume estimate is invalid");if(c.classification==CellClassification::Inside||(c.classification==CellClassification::Cut&&f>0.0))r.push_back(i);}return r;
	}
	static std::vector<std::int32_t> NodesFor(const CubicCartesianBackground& bg,const std::vector<std::uint64_t>& cells,std::size_t maximum_nodes)
	{
		std::set<std::int32_t> unique_nodes;
		for (const auto cell_id : cells) {
			const auto element = bg.MaterializeElement(cell_id);
			if (element.connectivity.size() != 64) throw std::runtime_error("immersed velocity extension requires cubic elements");
			for (const auto node_id : element.connectivity) {
				if (unique_nodes.find(node_id) != unique_nodes.end()) continue;
				if (unique_nodes.size() >= maximum_nodes) throw std::length_error("immersed velocity extension node cap exceeded");
				unique_nodes.insert(node_id);
#ifdef IGA_IMMERSED_VELOCITY_EXTENSION_TESTING
				TestingCapCountersState().layout_nodes = unique_nodes.size();
#endif
			}
		}
		return {unique_nodes.begin(), unique_nodes.end()};
	}
	static void RequireLayoutNodes(const MovingCutGeometry& g,const std::vector<std::uint64_t>& cells,const ImmersedActiveLayout& l,std::size_t maximum_nodes,const char* label)
	{if(l.NodeIds()!=NodesFor(g.Domain().Background(),cells,maximum_nodes))throw std::invalid_argument(std::string("immersed velocity extension ")+label+" layout is not the exact positive-cell node union");}
	std::vector<std::int32_t> BoundedBandNodes(const CubicCartesianBackground& bg) const
	{
		std::set<std::int32_t> unique_nodes;
		for (const auto cell_id : band_cells_) {
			const auto element = bg.MaterializeElement(cell_id);
			if (element.connectivity.size() != 64) throw std::runtime_error("immersed velocity extension requires cubic elements");
			for (const auto node_id : element.connectivity) {
				if (unique_nodes.find(node_id) != unique_nodes.end()) continue;
				if (unique_nodes.size() >= options_.max_nodes) throw std::length_error("immersed velocity extension node cap exceeded");
				unique_nodes.insert(node_id);
#ifdef IGA_IMMERSED_VELOCITY_EXTENSION_TESTING
				TestingCapCountersState().band_nodes = unique_nodes.size();
#endif
			}
		}
		return {unique_nodes.begin(), unique_nodes.end()};
	}
	void BuildBand(const CubicCartesianBackground& bg,const std::vector<std::uint64_t>& old)
	{
		if (old.size() > options_.max_band_cells) throw std::length_error("immersed velocity extension old positive-cell set exceeds band cap");
		std::map<std::uint64_t,std::uint32_t> distances;
		std::vector<std::uint64_t> frontier; frontier.reserve(old.size());
		for (const auto cell_id : old) {
			if (distances.find(cell_id) != distances.end()) continue;
			if (distances.size() >= options_.max_band_cells) throw std::length_error("immersed velocity extension band cell cap exceeded");
			distances.emplace(cell_id, 0); frontier.push_back(cell_id);
#ifdef IGA_IMMERSED_VELOCITY_EXTENSION_TESTING
			TestingCapCountersState().forward_visited = distances.size();
#endif
		}
		const std::array<std::size_t,6> order{{CubicCartesianBackground::XMinus,CubicCartesianBackground::XPlus,CubicCartesianBackground::YMinus,CubicCartesianBackground::YPlus,CubicCartesianBackground::ZMinus,CubicCartesianBackground::ZPlus}};
		for(std::uint64_t d=1; d<=static_cast<std::uint64_t>(layers_); ++d){
			if(frontier.empty()) break;
			const auto stored_distance=NarrowDistance(d,"immersed velocity extension band distance exceeds storage");
			std::vector<std::uint64_t> next;
			for(auto c:frontier){const auto cell=bg.Cell(c);for(auto slot:order){const auto q=cell.neighbor[slot];if(q==CubicCartesianBackground::kNoNeighbor||distances.find(q)!=distances.end())continue;
#ifdef IGA_IMMERSED_VELOCITY_EXTENSION_TESTING
				++TestingCapCountersState().forward_insertion_attempts;
#endif
				if(distances.size()>=options_.max_band_cells)throw std::length_error("immersed velocity extension band cell cap exceeded");
				distances.emplace(q,stored_distance); next.push_back(q);
#ifdef IGA_IMMERSED_VELOCITY_EXTENSION_TESTING
				TestingCapCountersState().forward_visited = distances.size();
#endif
			}}
			std::sort(next.begin(),next.end());
			frontier=std::move(next);
		}
		band_cells_.clear(); band_distances_.clear(); band_cells_.reserve(distances.size()); band_distances_.reserve(distances.size());
		for (const auto& value : distances) { band_cells_.push_back(value.first); band_distances_.push_back(value.second); }
	}
	void CheckNewCells()
	{
		diagnostics_.old_positive_cells=old_cells_.size();diagnostics_.new_positive_cells=new_cells_.size();diagnostics_.band_cells=band_cells_.size();std::unordered_map<std::uint64_t,std::uint32_t> d;for(std::size_t i=0;i<band_cells_.size();++i)d.emplace(band_cells_[i],band_distances_[i]);
		for(auto c:new_cells_){const auto p=d.find(c);if(p==d.end()||p->second>layers_)throw std::invalid_argument("immersed velocity extension new positive cell lies outside forward band");diagnostics_.maximum_forward_distance=std::max(diagnostics_.maximum_forward_distance,p->second);}
		if (new_cells_.size() > options_.max_reverse_visited_cells) throw std::length_error("immersed velocity extension reverse diagnostic cell cap exceeded");
		std::map<std::uint64_t,std::uint32_t> reverse;
		std::vector<std::uint64_t> frontier; frontier.reserve(new_cells_.size()); std::set<std::uint64_t> unresolved(old_cells_.begin(),old_cells_.end());
		for(const auto cell_id:new_cells_) {
			if (reverse.find(cell_id) != reverse.end()) continue;
#ifdef IGA_IMMERSED_VELOCITY_EXTENSION_TESTING
			++TestingCapCountersState().reverse_insertion_attempts;
#endif
			if (reverse.size() >= options_.max_reverse_visited_cells) throw std::length_error("immersed velocity extension reverse diagnostic cell cap exceeded");
			reverse.emplace(cell_id,0); frontier.push_back(cell_id); unresolved.erase(cell_id);
#ifdef IGA_IMMERSED_VELOCITY_EXTENSION_TESTING
			TestingCapCountersState().reverse_visited = reverse.size();
#endif
		}
#ifdef IGA_IMMERSED_VELOCITY_EXTENSION_TESTING
		if (TestingFaultActive(TestingFault::ReverseEmptyFrontier)) frontier.clear();
#endif
		const std::array<std::size_t,6> order{{CubicCartesianBackground::XMinus,CubicCartesianBackground::XPlus,CubicCartesianBackground::YMinus,CubicCartesianBackground::YPlus,CubicCartesianBackground::ZMinus,CubicCartesianBackground::ZPlus}};
		for(std::uint64_t distance=1;!frontier.empty()&&!unresolved.empty();++distance){const auto stored_distance=NarrowDistance(distance,"immersed velocity extension reverse distance exceeds storage");std::vector<std::uint64_t> next;for(auto c:frontier)for(auto slot:order){const auto n=background_.Cell(c).neighbor[slot];if(n==CubicCartesianBackground::kNoNeighbor||reverse.find(n)!=reverse.end())continue;
#ifdef IGA_IMMERSED_VELOCITY_EXTENSION_TESTING
			++TestingCapCountersState().reverse_insertion_attempts;
#endif
			if(reverse.size()>=options_.max_reverse_visited_cells)throw std::length_error("immersed velocity extension reverse diagnostic cell cap exceeded");
			reverse.emplace(n,stored_distance);unresolved.erase(n);next.push_back(n);
#ifdef IGA_IMMERSED_VELOCITY_EXTENSION_TESTING
			TestingCapCountersState().reverse_visited = reverse.size();
#endif
		}std::sort(next.begin(),next.end());frontier=std::move(next);}
		if(!unresolved.empty()) throw std::runtime_error("immersed velocity extension reverse diagnostic cannot reach old positive cell");
		for(auto c:old_cells_) diagnostics_.maximum_reverse_distance=std::max(diagnostics_.maximum_reverse_distance,reverse.at(c));

#ifdef IGA_IMMERSED_VELOCITY_EXTENSION_TESTING
		reverse_cells_.clear(); reverse_distances_.clear(); reverse_cells_.reserve(reverse.size()); reverse_distances_.reserve(reverse.size()); for(const auto& value:reverse){reverse_cells_.push_back(value.first);reverse_distances_.push_back(value.second);}
#endif
	}
	void BuildNodes(const ImmersedActiveLayout& old_layout,const ImmersedGlobalFlowState& old_state)
	{old_layout_nodes_=old_layout.NodeIds();band_nodes_=BoundedBandNodes(background_); anchor_nodes_=old_layout.NodeIds();unknown_nodes_.clear();unknown_nodes_.reserve(std::min(options_.max_unknowns,band_nodes_.size()));for(auto n:band_nodes_)if(!std::binary_search(anchor_nodes_.begin(),anchor_nodes_.end(),n)){
			if(unknown_nodes_.size()>=options_.max_unknowns)throw std::length_error("immersed velocity extension unknown cap exceeded");
			unknown_nodes_.push_back(n);
#ifdef IGA_IMMERSED_VELOCITY_EXTENSION_TESTING
			TestingCapCountersState().unknown_nodes=unknown_nodes_.size();
#endif
		}anchor_band_positions_.reserve(anchor_nodes_.size());old_layout_positions_.reserve(anchor_nodes_.size());for(std::size_t i=0;i<anchor_nodes_.size();++i){anchor_band_positions_.push_back(BandPosition(anchor_nodes_[i]));old_layout_positions_.push_back(i);}
		unknown_band_positions_.reserve(unknown_nodes_.size());for(auto n:unknown_nodes_)unknown_band_positions_.push_back(BandPosition(n));diagnostics_.old_nodes=old_layout.NodeIds().size();diagnostics_.band_nodes=band_nodes_.size();diagnostics_.anchors=anchor_nodes_.size();diagnostics_.unknowns=unknown_nodes_.size();(void)old_state;}
	std::size_t BandPosition(std::int32_t n) const {const auto it=std::lower_bound(band_nodes_.begin(),band_nodes_.end(),n);if(it==band_nodes_.end()||*it!=n)throw std::out_of_range("immersed velocity extension node is absent from band");return static_cast<std::size_t>(it-band_nodes_.begin());}
	void BuildFaces(const CubicCartesianBackground& bg)
	{
		background_=bg;grid_=bg.Spec(); trace_entries_=0;
		std::unordered_map<std::uint64_t,unsigned char> in; in.reserve(band_cells_.size());for(auto c:band_cells_)in.emplace(c,1);
		for(auto m:band_cells_)for(std::uint8_t axis=0;axis<3;++axis){
			const auto cm=bg.Cell(m); const std::size_t slot=axis==0?CubicCartesianBackground::XPlus:axis==1?CubicCartesianBackground::YPlus:CubicCartesianBackground::ZPlus;const auto p=cm.neighbor[slot];if(p==CubicCartesianBackground::kNoNeighbor||in.find(p)==in.end())continue;
			if(faces_.size()>=options_.max_faces)throw std::length_error("immersed velocity extension face cap exceeded");
			FaceTrace f;f.minus_cell=m;f.plus_cell=p;f.axis=axis;f.nodes=CubicCartesianSplineFaceConnectivityFor(bg,m,p);if(f.nodes.empty()||f.nodes.size()>options_.max_nodes)throw std::length_error("immersed velocity extension face node cap exceeded");if(f.nodes.size()>options_.max_trace_entries||trace_entries_>options_.max_trace_entries-f.nodes.size())throw std::length_error("immersed velocity extension trace entry cap exceeded"); f.h=cm.upper_m[axis]-cm.lower_m[axis];const int a=(axis+1)%3,b=(axis+2)%3;f.area=(cm.upper_m[a]-cm.lower_m[a])*(cm.upper_m[b]-cm.lower_m[b]);const auto minus=bg.MaterializeElement(m),plus=bg.MaterializeElement(p);const auto& quadrature=CubicCartesianSplineFaceTangentialQuadrature();f.alpha.reserve(quadrature.size());f.jumps.reserve(quadrature.size());
			for(const auto&q:quadrature){if(trace_entries_>options_.max_trace_entries-f.nodes.size())throw std::length_error("immersed velocity extension trace entry cap exceeded");const double alpha=std::pow(f.h,5)*q.weight*f.area/4.0;auto j=CubicCartesianSplineFaceNormalDerivativeJump(minus,plus,axis,3,q.tangential_0,q.tangential_1,f.h,f.nodes);if(!std::isfinite(alpha)||!(alpha>0.0)||j.size()!=f.nodes.size())throw std::runtime_error("immersed velocity extension face coefficient is invalid");for(double x:j)if(!std::isfinite(x))throw std::runtime_error("immersed velocity extension trace jump is nonfinite");f.alpha.push_back(alpha);f.jumps.push_back(std::move(j));trace_entries_=Add(trace_entries_,f.nodes.size(),"immersed velocity extension trace entry count overflows");}
			faces_.push_back(std::move(f));
		}
		diagnostics_.faces=faces_.size(); diagnostics_.trace_entries=trace_entries_;
	}
	static std::vector<std::int32_t> CubicCartesianSplineFaceConnectivityFor(const CubicCartesianBackground& bg,std::uint64_t m,std::uint64_t p){auto r=bg.MaterializeElement(m).connectivity;const auto q=bg.MaterializeElement(p).connectivity;r.insert(r.end(),q.begin(),q.end());std::sort(r.begin(),r.end());r.erase(std::unique(r.begin(),r.end()),r.end());return r;}
	void PreflightTraceGraph() const
	{
		if (unknown_nodes_.empty()) return;
		const std::size_t count=band_nodes_.size(); if(count>options_.max_nodes)throw std::length_error("immersed velocity extension trace graph node cap exceeded");
		std::vector<std::size_t> parent(count), rank(count); std::vector<unsigned char> anchored(count);
		for(std::size_t i=0;i<count;++i){parent[i]=i;anchored[i]=std::binary_search(anchor_nodes_.begin(),anchor_nodes_.end(),band_nodes_[i])?1u:0u;}
		auto root=[&](std::size_t value){while(parent[value]!=value){parent[value]=parent[parent[value]];value=parent[value];}return value;};
		auto join=[&](std::size_t left,std::size_t right){left=root(left);right=root(right);if(left==right)return;if(rank[left]<rank[right])std::swap(left,right);parent[right]=left;anchored[left]=static_cast<unsigned char>(anchored[left]||anchored[right]);if(rank[left]==rank[right])++rank[left];};
		std::size_t unions=0;
		for(const auto& face:faces_) for(const auto& jump:face.jumps) {
			if(jump.size()!=face.nodes.size())throw std::runtime_error("immersed velocity extension trace graph row is malformed");
			std::size_t first=count;
			for(std::size_t i=0;i<jump.size();++i){if(!std::isfinite(jump[i]))throw std::runtime_error("immersed velocity extension trace graph is nonfinite");if(jump[i]==0.0)continue;const auto position=BandPosition(face.nodes[i]);if(first==count){first=position;continue;}if(unions>=options_.max_trace_entries)throw std::length_error("immersed velocity extension trace graph union cap exceeded");++unions;join(first,position);}
		}
		for(const auto node:unknown_nodes_){const auto component=root(BandPosition(node));
#ifdef IGA_IMMERSED_VELOCITY_EXTENSION_TESTING
			if (TestingFaultActive(TestingFault::TraceAnchoring)) anchored[component]=0;
#endif
			if(!anchored[component])throw std::runtime_error("immersed velocity extension trace graph is not anchored");}
	}
	void BuildMatrixAndFactor()
	{
		const auto n=unknown_nodes_.size();if(n>options_.max_unknowns)throw std::length_error("immersed velocity extension unknown cap exceeded");if(!n)return;const auto entries=Mul(n,n,"immersed velocity extension dense matrix overflows");const auto matrix_bytes=Mul(entries,sizeof(double),"immersed velocity extension dense matrix bytes overflow");const auto retained_bytes=Add(Mul(matrix_bytes,2,"immersed velocity extension retained dense bytes overflow"),Mul(n,2*sizeof(double),"immersed velocity extension diagonal bytes overflow"),"immersed velocity extension retained dense bytes overflow");const auto work=Mul(entries,n,"immersed velocity extension work overflows");if(retained_bytes>options_.max_dense_bytes||work>options_.max_factor_work)throw std::length_error("immersed velocity extension dense factor cap exceeded");
		upper_.assign(entries,0.0);std::unordered_map<std::int32_t,std::size_t> u;u.reserve(n);for(std::size_t i=0;i<n;++i)u.emplace(unknown_nodes_[i],i);for(const auto& f:faces_){if(f.alpha.size()!=f.jumps.size())throw std::runtime_error("immersed velocity extension face trace framing is inconsistent");for(std::size_t q=0;q<f.alpha.size();++q){if(!std::isfinite(f.alpha[q])||f.jumps[q].size()!=f.nodes.size())throw std::runtime_error("immersed velocity extension matrix trace is invalid");for(std::size_t a=0;a<f.nodes.size();++a){const auto ia=u.find(f.nodes[a]);if(ia==u.end())continue;for(std::size_t b=a;b<f.nodes.size();++b){const auto ib=u.find(f.nodes[b]);if(ib==u.end())continue;std::size_t i=ia->second,j=ib->second;if(i>j)std::swap(i,j);upper_[i*n+j]=AddProduct(upper_[i*n+j],f.alpha[q],f.jumps[q][a],f.jumps[q][b],"immersed velocity extension matrix accumulation is nonfinite");}}}}
		diagonal_.resize(n);for(std::size_t i=0;i<n;++i){diagonal_[i]=upper_[i*n+i];
#ifdef IGA_IMMERSED_VELOCITY_EXTENSION_TESTING
			if (TestingFaultActive(TestingFault::Diagonal)) diagonal_[i]=0.0;
#endif
			if(!std::isfinite(diagonal_[i])||!(diagonal_[i]>0.0))throw std::runtime_error("immersed velocity extension unknown has nonpositive finite diagonal");}diagnostics_.minimum_diagonal=*std::min_element(diagonal_.begin(),diagonal_.end());diagnostics_.maximum_diagonal=*std::max_element(diagonal_.begin(),diagonal_.end());
		sqrt_diagonal_.resize(n);for(std::size_t i=0;i<n;++i){sqrt_diagonal_[i]=std::sqrt(diagonal_[i]);if(!std::isfinite(sqrt_diagonal_[i])||!(sqrt_diagonal_[i]>0.0))throw std::runtime_error("immersed velocity extension square-root diagonal is invalid");}
		factor_.assign(entries,0.0);for(std::size_t i=0;i<n;++i)for(std::size_t j=i;j<n;++j){const double k=upper_[i*n+j];factor_[i*n+j]=factor_[j*n+i]=ScaleCoefficient(k,sqrt_diagonal_[i],sqrt_diagonal_[j]);}
		const double threshold=256.0*std::numeric_limits<double>::epsilon()*n;
		diagnostics_.pivot_threshold=threshold;
		for(std::size_t i=0;i<n;++i){
			double p=factor_[i*n+i];
			for(std::size_t k=0;k<i;++k) p=FiniteDouble(static_cast<long double>(p)-static_cast<long double>(factor_[i*n+k])*factor_[i*n+k],"immersed velocity extension Cholesky pivot is nonfinite");
#ifdef IGA_IMMERSED_VELOCITY_EXTENSION_TESTING
			if (TestingFaultActive(TestingFault::Pivot)) p=0.0;
#endif
			if(!std::isfinite(p)||!(p>threshold)) throw std::runtime_error("immersed velocity extension scaled Cholesky pivot gate failed at unknown index "+std::to_string(i));
			factor_[i*n+i]=std::sqrt(p);
			if(!std::isfinite(factor_[i*n+i])) throw std::runtime_error("immersed velocity extension Cholesky pivot square root is nonfinite");
			for(std::size_t j=i+1;j<n;++j){
				double v=factor_[j*n+i];
				for(std::size_t k=0;k<i;++k) v=FiniteDouble(static_cast<long double>(v)-static_cast<long double>(factor_[j*n+k])*factor_[i*n+k],"immersed velocity extension Cholesky entry is nonfinite");
				factor_[j*n+i]=FiniteDouble(static_cast<long double>(v)/factor_[i*n+i],"immersed velocity extension Cholesky entry is not representable");
			}
			diagnostics_.minimum_pivot=i?std::min(diagnostics_.minimum_pivot,p):p;
			diagnostics_.maximum_pivot=std::max(diagnostics_.maximum_pivot,p);
		}
	}
	std::vector<double> BuildRhs(const std::vector<double>& full) const {if(full.size()!=band_nodes_.size())throw std::invalid_argument("immersed velocity extension RHS coverage is invalid");for(double value:full)if(!std::isfinite(value))throw std::runtime_error("immersed velocity extension RHS state is nonfinite");std::vector<double>b(unknown_nodes_.size());std::unordered_map<std::int32_t,std::size_t>u;u.reserve(unknown_nodes_.size());for(std::size_t i=0;i<unknown_nodes_.size();++i)u.emplace(unknown_nodes_[i],i);for(const auto&f:faces_)for(std::size_t q=0;q<f.alpha.size();++q){long double anchor=0;for(std::size_t a=0;a<f.nodes.size();++a)if(std::binary_search(anchor_nodes_.begin(),anchor_nodes_.end(),f.nodes[a]))anchor+=static_cast<long double>(f.jumps[q][a])*full[BandPosition(f.nodes[a])];const double finite_anchor=FiniteDouble(anchor,"immersed velocity extension RHS anchor is not representable");for(std::size_t i=0;i<f.nodes.size();++i){const auto it=u.find(f.nodes[i]);if(it!=u.end())b[it->second]=AddProduct(b[it->second],-f.alpha[q],f.jumps[q][i],finite_anchor,"immersed velocity extension RHS accumulation is nonfinite");}}return b;}
	std::vector<double> SolveComponent(const std::vector<double>& full) const {const auto n=unknown_nodes_.size();if(!n)return{};auto b=BuildRhs(full);for(std::size_t i=0;i<n;++i)b[i]=ScaleValue(b[i],sqrt_diagonal_[i],"immersed velocity extension scaled RHS is not representable");for(std::size_t i=0;i<n;++i){for(std::size_t k=0;k<i;++k)b[i]=FiniteDouble(static_cast<long double>(b[i])-static_cast<long double>(factor_[i*n+k])*b[k],"immersed velocity extension forward solve is nonfinite");b[i]=FiniteDouble(static_cast<long double>(b[i])/factor_[i*n+i],"immersed velocity extension forward solve is not representable");}for(std::size_t z=n;z-- >0;){for(std::size_t k=z+1;k<n;++k)b[z]=FiniteDouble(static_cast<long double>(b[z])-static_cast<long double>(factor_[k*n+z])*b[k],"immersed velocity extension backward solve is nonfinite");b[z]=FiniteDouble(static_cast<long double>(b[z])/factor_[z*n+z],"immersed velocity extension backward solve is not representable");}for(std::size_t i=0;i<n;++i)b[i]=ScaleValue(b[i],sqrt_diagonal_[i],"immersed velocity extension unscaled solution is not representable");return b;}
	static long double StableNorm(const std::vector<double>& values,const char* message){long double scale=0.0L,sum=1.0L;for(double value:values){if(!std::isfinite(value))throw std::runtime_error(message);const long double magnitude=std::abs(static_cast<long double>(value));if(magnitude==0.0L)continue;if(scale<magnitude){const long double ratio=scale/magnitude;sum=1.0L+sum*ratio*ratio;scale=magnitude;}else{const long double ratio=magnitude/scale;sum+=ratio*ratio;}}return scale==0.0L?0.0L:scale*std::sqrt(sum);}
	double Residual(const std::vector<double>& full,const std::vector<double>& x,const std::vector<double>& b) const {const auto n=x.size();if(!n)return 0.0;if(full.size()!=band_nodes_.size()||b.size()!=n)throw std::invalid_argument("immersed velocity extension residual coverage is invalid");std::vector<double> r(n);for(const auto&f:faces_)for(std::size_t q=0;q<f.alpha.size();++q){long double jv=0;for(std::size_t a=0;a<f.nodes.size();++a){const double value=full[BandPosition(f.nodes[a])];if(!std::isfinite(value)||!std::isfinite(f.jumps[q][a])||!std::isfinite(f.alpha[q]))throw std::runtime_error("immersed velocity extension residual contribution is nonfinite");jv+=static_cast<long double>(f.jumps[q][a])*value;}const double finite_jv=FiniteDouble(jv,"immersed velocity extension residual jump is not representable");for(std::size_t a=0;a<f.nodes.size();++a){const auto it=std::lower_bound(unknown_nodes_.begin(),unknown_nodes_.end(),f.nodes[a]);if(it!=unknown_nodes_.end()&&*it==f.nodes[a])r[static_cast<std::size_t>(it-unknown_nodes_.begin())]=AddProduct(r[static_cast<std::size_t>(it-unknown_nodes_.begin())],f.alpha[q],f.jumps[q][a],finite_jv,"immersed velocity extension residual accumulation is nonfinite");}}std::vector<double> rows(n);for(std::size_t i=0;i<n;++i){long double row=0;for(std::size_t j=0;j<n;++j){const double entry=i<=j?upper_[i*n+j]:upper_[j*n+i];if(!std::isfinite(entry))throw std::runtime_error("immersed velocity extension matrix norm entry is nonfinite");row+=std::abs(static_cast<long double>(entry));}rows[i]=FiniteDouble(row,"immersed velocity extension matrix norm is not representable");}const long double numerator=StableNorm(r,"immersed velocity extension residual norm is nonfinite");const long double matrix_norm=*std::max_element(rows.begin(),rows.end());const long double denominator=matrix_norm*StableNorm(x,"immersed velocity extension solution norm is nonfinite")+StableNorm(b,"immersed velocity extension RHS norm is nonfinite");if(!std::isfinite(numerator)||!std::isfinite(matrix_norm)||!std::isfinite(denominator))throw std::runtime_error("immersed velocity extension residual norm is nonfinite");if(denominator==0.0L)return numerator==0.0L?0.0:std::numeric_limits<double>::infinity();return FiniteDouble(numerator/denominator,"immersed velocity extension residual is not representable");}
	void ExtendVelocity(const ImmersedGlobalFlowState& state){std::array<std::vector<double>,3> values;for(int c=0;c<3;++c){values[c].resize(band_nodes_.size());for(std::size_t a=0;a<anchor_nodes_.size();++a)values[c][anchor_band_positions_[a]]=state.Coefficients()[old_layout_positions_[a]][c];const auto x=SolveComponent(values[c]);for(std::size_t u=0;u<x.size();++u)values[c][unknown_band_positions_[u]]=x[u];diagnostics_.residual[c]=Residual(values[c],x,BuildRhs(values[c]));
#ifdef IGA_IMMERSED_VELOCITY_EXTENSION_TESTING
			if (TestingFaultActive(TestingFault::VelocityResidual)) diagnostics_.residual[c]=std::numeric_limits<double>::quiet_NaN();
#endif
			if(!(std::isfinite(diagnostics_.residual[c])&&diagnostics_.residual[c]<=options_.residual_tolerance))throw std::runtime_error("immersed velocity extension residual gate failed");}band_velocities_.resize(band_nodes_.size());band_provenance_.resize(band_nodes_.size(),ImmersedVelocityHistoryProvenance::Extended);for(std::size_t i=0;i<band_nodes_.size();++i)for(int c=0;c<3;++c)band_velocities_[i][c]=values[c][i];for(std::size_t a=0;a<anchor_nodes_.size();++a){const auto p=anchor_band_positions_[a];band_provenance_[p]=ImmersedVelocityHistoryProvenance::Committed;for(int c=0;c<3;++c)if(std::memcmp(&band_velocities_[p][c],&state.Coefficients()[old_layout_positions_[a]][c],sizeof(double))!=0)throw std::runtime_error("immersed velocity extension anchor bitwise audit failed");}}
	void BuildHistory(const ImmersedActiveLayout& target){target_node_ids_=target.NodeIds();
#ifdef IGA_IMMERSED_VELOCITY_EXTENSION_TESTING
		if (TestingFaultActive(TestingFault::TargetCoverage) && !target_node_ids_.empty()) target_node_ids_.pop_back();
#endif
		std::vector<std::array<double,3>>v;std::vector<ImmersedVelocityHistoryProvenance>p;v.reserve(target_node_ids_.size());p.reserve(target_node_ids_.size());for(auto n:target_node_ids_){const auto i=BandPosition(n);v.push_back(band_velocities_[i]);p.push_back(band_provenance_[i]);}history_=ImmersedVelocityHistory(old_time_s_,new_time_s_,old_geometry_identity_,new_geometry_identity_,target_node_ids_,std::move(v),std::move(p));history_.ValidateCoverage(target);}
	std::string HashOperator()const{Sha256 h;String(h,"ImmersedVelocityExtension/operator/v2");for(double x:grid_.lower_m)h.AppendNormalizedDouble(x);for(double x:grid_.upper_m)h.AppendNormalizedDouble(x);for(auto x:grid_.cells)h.AppendLittleEndian32(x);h.AppendLittleEndian32(layers_);h.AppendLittleEndian64(band_cells_.size());for(std::size_t i=0;i<band_cells_.size();++i){h.AppendLittleEndian64(band_cells_[i]);h.AppendLittleEndian32(band_distances_[i]);}Ids(h,band_nodes_);h.AppendLittleEndian64(faces_.size());for(const auto&f:faces_){h.AppendLittleEndian64(f.minus_cell);h.AppendLittleEndian64(f.plus_cell);h.AppendLittleEndian32(f.axis);h.AppendNormalizedDouble(f.h);h.AppendNormalizedDouble(f.area);h.AppendLittleEndian64(f.nodes.size());Ids(h,f.nodes);h.AppendLittleEndian64(f.alpha.size());h.AppendLittleEndian64(f.jumps.size());for(std::size_t q=0;q<f.alpha.size();++q){h.AppendNormalizedDouble(f.alpha[q]);h.AppendLittleEndian64(f.jumps[q].size());for(std::size_t j=0;j<f.jumps[q].size();++j){h.AppendLittleEndian32(static_cast<std::uint32_t>(f.nodes[j]));h.AppendNormalizedDouble(f.jumps[q][j]);}}}return h.Hex();}
	std::string HashReduced()const{Sha256 h;String(h,"ImmersedVelocityExtension/reduced/v2");String(h,operator_hash_);Ids(h,unknown_nodes_);h.AppendLittleEndian64(upper_.size());for(double x:upper_)h.AppendNormalizedDouble(x);h.AppendLittleEndian64(diagonal_.size());for(double x:diagonal_)h.AppendNormalizedDouble(x);return h.Hex();}
	// Identity includes all public deterministic diagnostics that describe the
	// accepted topology/solve.  Temporary work buffers and recomputable factor
	// entries are intentionally omitted; reduced_hash_ binds their input matrix.
	void HashDiagnostics(Sha256& h)const{for(auto x:{diagnostics_.old_positive_cells,diagnostics_.new_positive_cells,diagnostics_.band_cells,diagnostics_.old_nodes,diagnostics_.band_nodes,diagnostics_.anchors,diagnostics_.unknowns,diagnostics_.faces,diagnostics_.trace_entries})h.AppendLittleEndian64(x);h.AppendLittleEndian32(diagnostics_.maximum_forward_distance);h.AppendLittleEndian32(diagnostics_.maximum_reverse_distance);for(double x:{diagnostics_.minimum_diagonal,diagnostics_.maximum_diagonal,diagnostics_.minimum_pivot,diagnostics_.maximum_pivot,diagnostics_.pivot_threshold,diagnostics_.residual[0],diagnostics_.residual[1],diagnostics_.residual[2]})h.AppendNormalizedDouble(x);}
	std::string HashExtension()const{Sha256 h;String(h,"ImmersedVelocityExtension/v2");for(const auto&s:{old_geometry_identity_,new_geometry_identity_,old_layout_identity_,new_layout_identity_,old_state_identity_,operator_hash_,reduced_hash_})String(h,s);h.AppendNormalizedDouble(old_time_s_);h.AppendNormalizedDouble(new_time_s_);h.AppendLittleEndian32(layers_);for(auto x:{options_.max_unknowns,options_.max_dense_bytes,options_.max_factor_work,options_.max_band_cells,options_.max_faces,options_.max_nodes,options_.max_trace_entries,options_.max_reverse_visited_cells})h.AppendLittleEndian64(x);h.AppendNormalizedDouble(options_.residual_tolerance);HashDiagnostics(h);Ids(h,band_nodes_);h.AppendLittleEndian64(band_velocities_.size());h.AppendLittleEndian64(band_provenance_.size());for(std::size_t i=0;i<band_nodes_.size();++i){for(double x:band_velocities_[i])h.AppendNormalizedDouble(x);h.AppendLittleEndian32(static_cast<std::uint32_t>(band_provenance_[i]));}String(h,history_.HashSha256());return h.Hex();}
	ImmersedVelocityExtensionOptions options_{};ImmersedVelocityExtensionDiagnostics diagnostics_{};std::uint32_t layers_=0;CubicCartesianGridSpec grid_{};CubicCartesianBackground background_{CubicCartesianGridSpec{{0,0,0},{1,1,1},{1,1,1}}};double old_time_s_=0,new_time_s_=0;std::string old_geometry_identity_,new_geometry_identity_,old_layout_identity_,new_layout_identity_,old_state_identity_,operator_hash_,reduced_hash_,extension_hash_;std::vector<std::uint64_t>old_cells_,new_cells_,band_cells_;std::vector<std::uint32_t>band_distances_;
#ifdef IGA_IMMERSED_VELOCITY_EXTENSION_TESTING
	std::vector<std::uint64_t>reverse_cells_;
	std::vector<std::uint32_t>reverse_distances_;
#endif
	std::vector<std::int32_t>old_layout_nodes_,band_nodes_,anchor_nodes_,unknown_nodes_,target_node_ids_;std::vector<std::size_t>anchor_band_positions_,old_layout_positions_,unknown_band_positions_;std::vector<FaceTrace>faces_;std::size_t trace_entries_=0;std::vector<double>upper_,diagonal_,sqrt_diagonal_,factor_;std::vector<std::array<double,3>>band_velocities_;std::vector<ImmersedVelocityHistoryProvenance>band_provenance_;ImmersedVelocityHistory history_{0,0,"x","x",{}, {}, {}};
};

} // namespace iga
#endif
