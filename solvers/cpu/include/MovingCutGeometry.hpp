#ifndef IGA_MOVING_CUT_GEOMETRY_HPP
#define IGA_MOVING_CUT_GEOMETRY_HPP

// Immutable correctness-oracle geometry for a material moving immersed
// boundary.  There is deliberately no reuse path: every evaluated time builds
// its own complete dependency chain before this object is made visible.
#include "MaterialSurfaceKinematics.hpp"
#include "CutCellGhostPenalty.hpp"
#include "ImmersedSurfaceQuadrature.hpp"
#include "Sha256.hpp"
#include "PhaseProfile.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

struct MovingCutGeometryOptions {
	OctreeCutQuadratureOptions volume;
	CutCellVolumeQuadratureStorageMode volume_storage = CutCellVolumeQuadratureStorageMode::Expanded;
	ImmersedSurfaceQuadratureOptions surface;
	CutCellGhostPenaltyOptions ghost;
};

struct MovingCutCellTransition {
	std::uint64_t id = 0;
	CellClassification old_classification = CellClassification::Outside;
	CellClassification new_classification = CellClassification::Outside;
};

struct MovingCutGeometryDiagnostics {
	double time_s = 0.0;
	std::string evaluated_motion_hash, surface_hash, geometry_identity_sha256, publication_identity_sha256;
	std::size_t outside_cells = 0, inside_cells = 0, cut_cells = 0, active_cells = 0;
	double surface_area_m2 = 0.0, source_surface_area_m2 = 0.0, surface_area_residual_m2 = 0.0;
	double catalog_lower_physical_volume_m3 = 0.0, catalog_estimated_physical_volume_m3 = 0.0;
	double catalog_upper_physical_volume_m3 = 0.0, closed_surface_physical_volume_m3 = 0.0, volume_residual_m3 = 0.0;
	CutCellVolumeQuadratureDiagnostics volume;
	ImmersedSurfaceQuadratureDiagnostics surface;
	CutCellGhostPenaltyDiagnostics ghost;
	std::size_t active_nodes = 0, ghost_faces = 0;
	std::size_t unchanged_cells = 0;
	std::map<std::pair<CellClassification, CellClassification>, std::size_t> transition_counts;
	std::vector<MovingCutCellTransition> transitions;
};

class MovingCutGeometry {
public:
	MovingCutGeometry(const MovingCutGeometry&) = delete;
	MovingCutGeometry& operator=(const MovingCutGeometry&) = delete;
	MovingCutGeometry(MovingCutGeometry&&) = delete;
	MovingCutGeometry& operator=(MovingCutGeometry&&) = delete;

	// Construct into a temporary and swap/replace the caller's committed owner
	// only after this returns.  A failed constructor cannot publish a partial
	// geometry state.
	static std::unique_ptr<MovingCutGeometry> Build(CubicCartesianGridSpec grid,
		MaterialSurfaceKinematics kinematics,
		MovingCutGeometryOptions options = {}, const MovingCutGeometry* previous = nullptr)
	{
		PhaseScope geometry_phase(ProfilePhase::Geometry);
		kinematics.Validate();
		if (previous) ValidatePreviousInput(grid, kinematics, *previous);
		return std::unique_ptr<MovingCutGeometry>(new MovingCutGeometry(std::move(grid), std::move(kinematics),
			std::move(options), previous));
	}

	const MaterialSurfaceKinematics& Kinematics() const noexcept { return kinematics_; }
	// Compatibility accessor for Phase 7 callers; its value is the neutral
	// material payload, not a prescribed-motion-specific representation.
	const MaterialSurfaceKinematics& Evaluation() const noexcept { return kinematics_; }
	const CartesianDomainClassification& Domain() const noexcept { return domain_; }
	const CutCellVolumeQuadratureCatalog& Volume() const noexcept { return volume_; }
	const ImmersedSurfaceQuadratureCatalog& Surface() const noexcept { return surface_; }
	const CutCellGhostPenaltyCatalog& Ghost() const noexcept { return ghost_; }
	const MovingCutGeometryOptions& Options() const noexcept { return options_; }
	const MovingCutGeometryDiagnostics& Diagnostics() const noexcept { return diagnostics_; }
	// Current geometry excludes predecessor context; publication identity adds
	// the retry-relevant transition context.
	const std::string& GeometryIdentitySha256() const noexcept { return diagnostics_.geometry_identity_sha256; }
	const std::string& PublicationIdentitySha256() const noexcept { return diagnostics_.publication_identity_sha256; }
	const std::string& IdentitySha256() const noexcept { return GeometryIdentitySha256(); }

private:
	MovingCutGeometry(CubicCartesianGridSpec grid, MaterialSurfaceKinematics kinematics,
		MovingCutGeometryOptions options, const MovingCutGeometry* previous)
		: kinematics_(std::move(kinematics)), domain_(CubicCartesianBackground(grid),
			SurfaceSpatialIndex(kinematics_.Surface())), volume_(domain_, options.volume, options.volume_storage),
		  surface_(domain_, options.surface), ghost_(domain_, volume_, options.ghost), options_(std::move(options))
	{
		kinematics_.Validate();
		ValidateCatalogs();
		diagnostics_.time_s = kinematics_.EvaluatedTimeS();
		diagnostics_.evaluated_motion_hash = kinematics_.ContentIdentitySha256();
		diagnostics_.surface_hash = domain_.SurfaceCanonicalHash();
		diagnostics_.outside_cells = domain_.Diagnostics().outside_count;
		diagnostics_.inside_cells = domain_.Diagnostics().inside_count;
		diagnostics_.cut_cells = domain_.Diagnostics().cut_count;
		diagnostics_.active_cells = domain_.ActiveCellIds().size();
		diagnostics_.volume = volume_.Diagnostics(); diagnostics_.surface = surface_.Diagnostics(); diagnostics_.ghost = ghost_.Diagnostics();
		diagnostics_.surface_area_m2 = surface_.Diagnostics().total_area_m2;
		diagnostics_.source_surface_area_m2 = surface_.Diagnostics().source_total_area_m2;
		diagnostics_.surface_area_residual_m2 = surface_.Diagnostics().total_area_residual_m2;
		diagnostics_.active_nodes = ghost_.Diagnostics().active_nodes; diagnostics_.ghost_faces = ghost_.Faces().size();
		diagnostics_.catalog_lower_physical_volume_m3 = volume_.Diagnostics().lower_physical_volume;
		diagnostics_.catalog_estimated_physical_volume_m3 = volume_.Diagnostics().estimated_physical_volume;
		diagnostics_.catalog_upper_physical_volume_m3 = volume_.Diagnostics().upper_physical_volume;
		diagnostics_.closed_surface_physical_volume_m3 = ClosedSurfaceVolume();
		diagnostics_.volume_residual_m3 = diagnostics_.catalog_estimated_physical_volume_m3-diagnostics_.closed_surface_physical_volume_m3;
		ValidateGlobalVolumeAudit();
		if (previous) ComparePrevious(*previous);
		diagnostics_.geometry_identity_sha256 = HashGeometryState();
		diagnostics_.publication_identity_sha256 = HashPublicationState(previous);
	}

	static bool SameGrid(const CubicCartesianGridSpec& left, const CubicCartesianGridSpec& right) noexcept
	{ return left.lower_m == right.lower_m && left.upper_m == right.upper_m && left.cells == right.cells; }
	static void ValidatePreviousInput(const CubicCartesianGridSpec& grid, const MaterialSurfaceKinematics& kinematics,
		const MovingCutGeometry& previous)
	{
		previous.kinematics_.Validate();
		if (!SameGrid(grid, previous.domain_.Background().Spec()))
			throw std::invalid_argument("moving cut geometry previous state has a different fixed grid");
		if (kinematics.SourceVerticesM().size() != previous.kinematics_.SourceVerticesM().size()
			|| kinematics.MaterialIdentitySha256() != previous.kinematics_.MaterialIdentitySha256()
			|| kinematics.TopologyIdentitySha256() != previous.kinematics_.TopologyIdentitySha256())
			throw std::invalid_argument("moving cut geometry previous state has different motion topology or labels");
	}
	double ClosedSurfaceVolume() const
	{
		const auto& vertices=kinematics_.Surface().Vertices(); const auto& triangles=kinematics_.Surface().Triangles();
		if (vertices.empty()) throw std::runtime_error("moving cut geometry surface has no vertices");
		const auto origin=vertices.front(); long double sum=0.0L;
		for (const auto& triangle : triangles) { std::array<long double,3> a{},b{},c{}; for (std::size_t axis=0;axis<3;++axis) { a[axis]=static_cast<long double>(vertices[triangle.indices[0]][axis])-origin[axis]; b[axis]=static_cast<long double>(vertices[triangle.indices[1]][axis])-origin[axis]; c[axis]=static_cast<long double>(vertices[triangle.indices[2]][axis])-origin[axis]; } sum += a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0]); }
		const double result=static_cast<double>(std::abs(sum)/6.0L); if (!std::isfinite(result)) throw std::runtime_error("moving cut geometry closed-surface volume is nonfinite"); return result;
	}
	void ValidateGlobalVolumeAudit() const
	{
		const double lower=diagnostics_.catalog_lower_physical_volume_m3, upper=diagnostics_.catalog_upper_physical_volume_m3, surface=diagnostics_.closed_surface_physical_volume_m3;
		const double tolerance=128.0*std::numeric_limits<double>::epsilon()*std::max({1.0,std::abs(lower),std::abs(upper),std::abs(surface)});
		if (!std::isfinite(lower)||!std::isfinite(upper)||!std::isfinite(surface)||surface < lower-tolerance||surface > upper+tolerance||std::abs(diagnostics_.volume_residual_m3) > (upper-lower)+tolerance) throw std::runtime_error("moving cut geometry global volume audit failed");
	}
	void ValidateCatalogs() const
	{
		if (kinematics_.Surface().CanonicalSha256() != domain_.SurfaceCanonicalHash()
			|| !SameGrid(domain_.Background().Spec(), volume_.GridSpec()) || !SameGrid(domain_.Background().Spec(), surface_.GridSpec())
			|| volume_.SurfaceCanonicalHash() != domain_.SurfaceCanonicalHash() || surface_.SurfaceCanonicalHash() != domain_.SurfaceCanonicalHash())
			throw std::runtime_error("moving cut geometry catalog surface/grid binding is inconsistent");
		if (!surface_.Usable() || !ghost_.Usable() || !ghost_.Matches(domain_, volume_))
			throw std::runtime_error("moving cut geometry catalog is unusable");
		for (std::uint64_t id = 0; id < domain_.Cells().size(); ++id) {
			const auto& cell = volume_.Cell(id);
			if (!cell.usable || cell.classification != domain_.Cells()[static_cast<std::size_t>(id)].classification)
				throw std::runtime_error("moving cut geometry volume catalog cell is unusable or mismatched");
			if (volume_.StorageMode() == CutCellVolumeQuadratureStorageMode::Expanded) volume_.ValidateUsableRule(domain_, id);
			else volume_.ValidateUsableCompactRule(domain_, id);
			if (domain_.Cells()[static_cast<std::size_t>(id)].classification == CellClassification::Cut) {
				surface_.ValidateUsableRule(domain_, id);
				const auto& rule = surface_.UsableRule(domain_, id); const auto& provenance = surface_.UsableProvenance(domain_, id);
				for (std::size_t point = 0; point < rule.Points().size(); ++point) {
					const auto& material = provenance[point]; const auto& triangle = kinematics_.Surface().Triangles().at(material.canonical_triangle);
					if (rule.Points()[point].boundary_id != static_cast<int>(triangle.boundary_id)) throw std::runtime_error("moving cut geometry surface label provenance is inconsistent");
					const auto& bounds = domain_.Cells()[id].bounds;
					std::array<double, 3> reconstruction{{0.0, 0.0, 0.0}}, local_point{{0.0, 0.0, 0.0}};
					double normal_dot = 0.0, residual = 0.0, triangle_scale = 0.0, cell_scale = 0.0;
					for (std::size_t corner = 0; corner < 3; ++corner) for (std::size_t axis = 0; axis < 3; ++axis) {
						const double local = kinematics_.Surface().Vertices()[triangle.indices[corner]][axis]-bounds.minimum[axis];
						reconstruction[axis] += material.canonical_barycentric[corner]*local;
						triangle_scale = std::max(triangle_scale, std::abs(local));
					}
					for (std::size_t axis = 0; axis < 3; ++axis) {
						const double width = bounds.maximum[axis]-bounds.minimum[axis];
						// The parametric coordinate is authoritative: physical coordinates
						// can lose cell-local bits after adding a large cell origin.
						local_point[axis] = rule.Points()[point].parametric[axis]*width;
						cell_scale = std::max(cell_scale, std::abs(width));
						residual = std::max(residual, std::abs(reconstruction[axis]-local_point[axis]));
						normal_dot += rule.Points()[point].normal[axis]*triangle.outward_unit_normal[axis];
					}
					const double provenance_scale = std::max(triangle_scale, cell_scale);
					if (residual > 2.0e-10*provenance_scale
						|| normal_dot < 1.0-128.0*std::numeric_limits<double>::epsilon())
						throw std::runtime_error("moving cut geometry surface physical or normal provenance is inconsistent");
					(void)kinematics_.WallVelocity(material.canonical_triangle, material.canonical_barycentric);
				}
			}
		}
		if (!std::isfinite(surface_.Diagnostics().total_area_m2)
			|| surface_.Diagnostics().total_area_absolute_residual_m2 > 1.0e-8*std::max(1.0, surface_.Diagnostics().source_total_area_m2))
			throw std::runtime_error("moving cut geometry surface-area audit failed");
		for (const auto& label : surface_.Diagnostics().area_absolute_residual_by_boundary_id)
			if (label.second > 1.0e-8*std::max(1.0, surface_.Diagnostics().source_area_by_boundary_id.at(label.first)))
				throw std::runtime_error("moving cut geometry label-area audit failed");
	}
	void ComparePrevious(const MovingCutGeometry& previous)
	{
		const auto& old_cells = previous.domain_.Cells(); const auto& new_cells = domain_.Cells();
		if (old_cells.size() != new_cells.size()) throw std::invalid_argument("moving cut geometry previous state has different cell count");
		for (std::size_t id = 0; id < new_cells.size(); ++id) {
			if (old_cells[id].id != id || new_cells[id].id != id) throw std::runtime_error("moving cut geometry cell ordering is invalid");
			const auto key = std::make_pair(old_cells[id].classification, new_cells[id].classification);
			if (key.first == key.second) {
				++diagnostics_.unchanged_cells;
				continue;
			}
			++diagnostics_.transition_counts[key];
			diagnostics_.transitions.push_back({static_cast<std::uint64_t>(id), key.first, key.second});
		}
	}
	static void AppendCount(Sha256& hash, std::size_t value) { hash.AppendLittleEndian64(value); }
	static void AppendString(Sha256& hash, const std::string& value) { AppendCount(hash, value.size()); hash.Append(value.data(), value.size()); }
	static void AppendClassification(Sha256& hash, CellClassification value) { hash.AppendLittleEndian32(static_cast<std::uint32_t>(value)); }
	static void AppendVolumeDiagnostics(Sha256& hash, const CutCellVolumeQuadratureDiagnostics& value)
	{
		for (double item : {value.certified_reference_volume, value.estimated_reference_volume, value.unresolved_reference_volume,
			value.lower_reference_volume, value.upper_reference_volume, value.certified_physical_volume,
			value.estimated_physical_volume, value.unresolved_physical_volume, value.lower_physical_volume,
			value.upper_physical_volume}) hash.AppendNormalizedDouble(item);
		for (std::size_t item : {value.nodes, value.leaves, value.output_points, value.samples, value.boundary_samples,
			value.ambiguous_samples, value.predicate_ambiguities, value.precision_limited_leaves, value.certified_blocks,
			value.sample_leaves, value.record_attempts, value.rolled_back_records, value.retained_bytes,
			value.observed_retained_bytes, value.logical_output_points, value.rescue_attempts, value.attempted_nodes,
			value.attempted_leaves, value.attempted_output_points, value.attempted_samples,
			value.attempted_record_attempts, value.attempted_logical_output_points}) AppendCount(hash, item);
		hash.AppendLittleEndian32(value.reached_depth); hash.AppendLittleEndian32(value.rescue_effective_depth);
	}
	static void AppendGhostDiagnostics(Sha256& hash, const CutCellGhostPenaltyDiagnostics& value)
	{
		for (std::size_t item : {value.active_cells, value.full_cells, value.cut_cells, value.active_nodes,
			value.candidate_faces, value.selected_faces, value.cut_inside_faces, value.cut_cut_faces,
			value.uncovered_cells, value.uncovered_components, value.quadrature_points}) AppendCount(hash, item);
		for (std::size_t item : value.selected_by_axis) AppendCount(hash, item);
		for (double item : {value.minimum_h_normal_m, value.maximum_h_normal_m, value.minimum_face_area_m2,
			value.maximum_face_area_m2, value.maximum_abs_trace_jump, value.maximum_abs_coefficient}) hash.AppendNormalizedDouble(item);
	}
	void AppendVolumeCell(Sha256& hash, const CutCellVolumeQuadratureCell& cell) const
	{
		hash.AppendLittleEndian64(cell.id); AppendClassification(hash, cell.classification); hash.AppendLittleEndian32(cell.usable ? 1u : 0u);
		AppendVolumeDiagnostics(hash, cell.diagnostics);
		AppendCount(hash, cell.rule.Points().size());
		for (const auto& point : cell.rule.Points()) { for (double value : point.parametric) hash.AppendNormalizedDouble(value); hash.AppendNormalizedDouble(point.weight); }
		AppendCompactCutCellVolumeRuleHash(hash,cell.compact_rule);
	}
	std::string HashGeometryState() const
	{
		Sha256 hash; AppendString(hash, "MovingCutGeometry/current/v5");
		// The content digest binds current coordinates and wall velocities even
		// when an adapter retains a compatibility epoch identifier.
		AppendString(hash, kinematics_.ContentIdentitySha256());
		for (double value : domain_.Background().Spec().lower_m) hash.AppendNormalizedDouble(value);
		for (double value : domain_.Background().Spec().upper_m) hash.AppendNormalizedDouble(value);
		for (auto value : domain_.Background().Spec().cells) hash.AppendLittleEndian32(value);
		hash.AppendLittleEndian32(options_.volume.max_depth); hash.AppendLittleEndian64(options_.volume.max_nodes); hash.AppendLittleEndian64(options_.volume.max_leaves); hash.AppendLittleEndian64(options_.volume.max_points); hash.AppendLittleEndian64(options_.volume.max_records); hash.AppendLittleEndian64(options_.volume.max_retained_bytes); hash.AppendLittleEndian64(options_.volume.max_logical_points); hash.AppendLittleEndian32(options_.volume.empty_rule_rescue_max_depth); hash.AppendLittleEndian32(static_cast<std::uint32_t>(options_.volume_storage));
		hash.AppendLittleEndian64(options_.surface.max_candidates); hash.AppendLittleEndian64(options_.surface.max_fragments); hash.AppendLittleEndian64(options_.surface.max_points); hash.AppendLittleEndian64(options_.surface.max_exact_limbs);
		hash.AppendNormalizedDouble(options_.ghost.gamma_u); hash.AppendNormalizedDouble(options_.ghost.gamma_p); hash.AppendLittleEndian64(options_.ghost.max_faces); hash.AppendLittleEndian64(options_.ghost.max_quadrature_points); hash.AppendLittleEndian64(options_.ghost.max_trace_entries);
		AppendCount(hash, domain_.Cells().size());
		for (const auto& cell : domain_.Cells()) {
			hash.AppendLittleEndian64(cell.id); for (auto value : cell.index) hash.AppendLittleEndian32(value);
			for (double value : cell.bounds.minimum) hash.AppendNormalizedDouble(value);
			for (double value : cell.bounds.maximum) hash.AppendNormalizedDouble(value);
			AppendClassification(hash, cell.classification); AppendCount(hash, cell.triangle_ids.size()); for (auto value : cell.triangle_ids) AppendCount(hash, value);
			hash.AppendLittleEndian32(cell.boundary_only_contact ? 1u : 0u); hash.AppendLittleEndian32(cell.ambiguous ? 1u : 0u);
		}
		AppendCount(hash, volume_.Cells().size()); for (const auto& cell : volume_.Cells()) AppendVolumeCell(hash, cell);
		AppendCount(hash, surface_.Cells().size());
		for (std::uint64_t id = 0; id < surface_.Cells().size(); ++id) {
			const auto& cell = surface_.Cell(id); hash.AppendLittleEndian64(cell.id); hash.AppendLittleEndian32(cell.usable ? 1u : 0u); hash.AppendLittleEndian32(cell.ambiguous ? 1u : 0u); hash.AppendNormalizedDouble(cell.area_m2); AppendCount(hash, cell.point_count);
			const auto& rule = surface_.UsableRule(domain_, id); const auto& provenance = surface_.UsableProvenance(domain_, id); AppendCount(hash, rule.Points().size());
			for (std::size_t point = 0; point < rule.Points().size(); ++point) { const auto& value=rule.Points()[point]; hash.AppendLittleEndian32(provenance[point].canonical_triangle); for (double coordinate : provenance[point].canonical_barycentric) hash.AppendNormalizedDouble(coordinate); for (double coordinate : value.parametric) hash.AppendNormalizedDouble(coordinate); for (double coordinate : value.physical) hash.AppendNormalizedDouble(coordinate); for (double coordinate : value.normal) hash.AppendNormalizedDouble(coordinate); hash.AppendNormalizedDouble(value.weight); hash.AppendLittleEndian32(static_cast<std::uint32_t>(value.boundary_id)); }
		}
		const auto& ghost_options = ghost_.Options(); hash.AppendNormalizedDouble(ghost_options.gamma_u); hash.AppendNormalizedDouble(ghost_options.gamma_p); AppendCount(hash, ghost_options.max_faces); AppendCount(hash, ghost_options.max_quadrature_points); AppendCount(hash, ghost_options.max_trace_entries);
		AppendGhostDiagnostics(hash, ghost_.Diagnostics()); hash.AppendLittleEndian32(ghost_.Usable() ? 1u : 0u); hash.AppendLittleEndian32(ghost_.Version());
		AppendCount(hash, ghost_.Faces().size()); for (const auto& face : ghost_.Faces()) { hash.AppendLittleEndian64(face.minus_cell); hash.AppendLittleEndian64(face.plus_cell); hash.AppendLittleEndian32(face.axis); hash.AppendNormalizedDouble(face.h_normal_m); hash.AppendNormalizedDouble(face.area_m2); }
		AppendCount(hash, ghost_.ActiveNodes().size()); for (auto node : ghost_.ActiveNodes()) hash.AppendLittleEndian32(static_cast<std::uint32_t>(node));
		AppendCount(hash, ghost_.UncoveredCells().size()); for (auto id : ghost_.UncoveredCells()) hash.AppendLittleEndian64(id);
		AppendCount(hash, domain_.Cells().size()); for (std::uint64_t id = 0; id < domain_.Cells().size(); ++id) hash.AppendLittleEndian32(ghost_.Covered(id) ? 1u : 0u);
		return hash.Hex();
	}
	std::string HashPublicationState(const MovingCutGeometry* previous) const
	{
		Sha256 hash; AppendString(hash,"MovingCutGeometry/publication/v3"); AppendString(hash,diagnostics_.geometry_identity_sha256); AppendString(hash,previous ? previous->GeometryIdentitySha256() : std::string{}); hash.AppendLittleEndian64(diagnostics_.unchanged_cells); hash.AppendLittleEndian64(diagnostics_.transitions.size()); for (const auto& item:diagnostics_.transition_counts) { AppendClassification(hash,item.first.first); AppendClassification(hash,item.first.second); hash.AppendLittleEndian64(item.second); } for (const auto& transition:diagnostics_.transitions) { hash.AppendLittleEndian64(transition.id); AppendClassification(hash,transition.old_classification); AppendClassification(hash,transition.new_classification); } return hash.Hex();
	}

	// Declaration order is the required lifetime order: evaluation first, then
	// the domain (which owns its spatial index), then its dependent catalogs.
	MaterialSurfaceKinematics kinematics_;
	CartesianDomainClassification domain_;
	CutCellVolumeQuadratureCatalog volume_;
	ImmersedSurfaceQuadratureCatalog surface_;
	CutCellGhostPenaltyCatalog ghost_;
	MovingCutGeometryOptions options_;
	MovingCutGeometryDiagnostics diagnostics_;
};

} // namespace iga

#endif
