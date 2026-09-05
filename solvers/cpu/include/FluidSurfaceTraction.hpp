#ifndef IGA_FLUID_SURFACE_TRACTION_HPP
#define IGA_FLUID_SURFACE_TRACTION_HPP

// Rank-local extraction of fluid-on-structure Cauchy traction.  This is a
// producer-neutral kernel: it reads an IGA state and immutable material
// geometry but does not publish it through a runtime or mutate either input.
#include "DistributedSurfaceInterface.hpp"
#include "FluidCauchyStress.hpp"
#include "ImmersedSurfaceQuadrature.hpp"
#include "MaterialSurfacePatchMap.hpp"
#include "MaterialSurfaceKinematics.hpp"
#include "NavierStokesElement.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

struct FluidSurfaceElementState {
	std::uint64_t cell_id = 0;
	// Cartesian velocity followed by pressure, in the exact local IGA order.
	std::vector<std::array<double, 4>> nodal_state;
};

struct FluidSurfaceTractionProjectionOptions {
	// This intentionally rejects an unimplemented distributed/dense global
	// assembly path.  The first immersed FSI slice is rank-local only.
	std::size_t maximum_nodes = 4096;
	// SI tolerances for the force [N] and moment [N m] conservation audits.
	// The relative term is scaled componentwise by accumulated absolute
	// contributions, so cancellation never makes an audit artificially strict.
	double conservation_absolute_force_tolerance_n = 1.0e-12;
	double conservation_absolute_moment_tolerance_n_m = 1.0e-12;
	double conservation_relative_tolerance = 2.0e-10;
};

struct FluidSurfaceTractionDiagnostics {
	std::array<double, 3> quadrature_resultant_n{};
	std::array<double, 3> nodal_resultant_n{};
	std::array<double, 3> quadrature_moment_n_m{};
	std::array<double, 3> nodal_moment_n_m{};
	std::size_t retained_quadrature_points = 0;
};

struct FluidSurfaceTractionResult {
	SurfaceTraction traction;
	FluidSurfaceTractionDiagnostics diagnostics;
};

namespace fluid_surface_traction_detail {

inline void AppendString(Sha256& hash, const std::string& value)
{ distributed_surface_detail::AppendString(hash, value); }

inline void AppendVector(Sha256& hash, const std::array<double, 3>& value)
{ distributed_surface_detail::AppendVector(hash, value); }

inline std::array<double, 3> Cross(const std::array<double, 3>& left,
	const std::array<double, 3>& right)
{
	return {{left[1]*right[2]-left[2]*right[1], left[2]*right[0]-left[0]*right[2],
		left[0]*right[1]-left[1]*right[0]}};
}

struct CompensatedSum {
	long double sum = 0.0L;
	long double correction = 0.0L;
	void Add(double value, const char* what)
	{
		if (!std::isfinite(value)) throw std::overflow_error(std::string(what)+" is nonfinite");
		const long double promoted = static_cast<long double>(value);
		const long double next = sum+promoted;
		if (std::abs(sum) >= std::abs(promoted)) correction += (sum-next)+promoted;
		else correction += (promoted-next)+sum;
		sum = next;
		if (!std::isfinite(sum) || !std::isfinite(correction))
			throw std::overflow_error(std::string(what)+" overflows");
	}
	double Value(const char* what) const
	{
		const long double value = sum+correction;
		if (!std::isfinite(value) || !std::isfinite(static_cast<double>(value)))
			throw std::overflow_error(std::string(what)+" overflows");
		return static_cast<double>(value);
	}
};

struct CompensatedVector {
	std::array<CompensatedSum, 3> signed_sum{};
	std::array<CompensatedSum, 3> absolute_sum{};
	void AddScaled(const std::array<double, 3>& value, double scale, const char* what)
	{
		if (!std::isfinite(scale)) throw std::overflow_error(std::string(what)+" scale is nonfinite");
		for (int component = 0; component < 3; ++component) {
			const double addition = value[component]*scale;
			if (!std::isfinite(addition)) throw std::overflow_error(std::string(what)+" is nonfinite");
			signed_sum[component].Add(addition, what);
			absolute_sum[component].Add(std::abs(addition), what);
		}
	}
	std::array<double, 3> Value(const char* what) const
	{
		std::array<double, 3> result{};
		for (int component = 0; component < 3; ++component)
			result[component] = signed_sum[component].Value(what);
		return result;
	}
};

inline std::uint64_t DoubleBits(double value)
{
	std::uint64_t result = 0;
	static_assert(sizeof(result) == sizeof(value), "double bit width is unsupported");
	std::memcpy(&result, &value, sizeof(result));
	return result;
}

inline std::array<std::uint64_t, 3> SortedTriangle(std::array<std::uint64_t, 3> value)
{ std::sort(value.begin(), value.end()); return value; }

inline void AppendCartesianDiscretization(Sha256& hash,
	const CartesianDomainClassification& domain, const ImmersedSurfaceQuadratureCatalog& catalog)
{
	if (!catalog.Usable() || catalog.SurfaceCanonicalHash() != domain.SurfaceCanonicalHash()
		|| catalog.GridSpec().lower_m != domain.Background().Spec().lower_m
		|| catalog.GridSpec().upper_m != domain.Background().Spec().upper_m
		|| catalog.GridSpec().cells != domain.Background().Spec().cells)
		throw std::invalid_argument("fluid surface traction quadrature catalog does not bind the Cartesian domain");
	AppendString(hash, "FluidSurfaceTraction/cartesian-discretization/v1");
	const auto& spec = domain.Background().Spec();
	for (const double value : spec.lower_m) hash.AppendNormalizedDouble(value);
	for (const double value : spec.upper_m) hash.AppendNormalizedDouble(value);
	for (const auto value : spec.cells) hash.AppendLittleEndian32(value);
	AppendString(hash, domain.SurfaceCanonicalHash());
	hash.AppendLittleEndian64(domain.Background().ElementCount());
	for (std::uint64_t id = 0; id < domain.Background().ElementCount(); ++id) {
		const auto& cell = domain.Cells().at(static_cast<std::size_t>(id));
		if (cell.id != id) throw std::runtime_error("fluid surface traction Cartesian domain cell ids are not canonical");
		hash.AppendLittleEndian32(static_cast<std::uint32_t>(cell.classification));
		for (const auto value : cell.index) hash.AppendLittleEndian32(value);
		for (const double value : cell.bounds.minimum) hash.AppendNormalizedDouble(value);
		for (const double value : cell.bounds.maximum) hash.AppendNormalizedDouble(value);
		hash.AppendLittleEndian32(cell.boundary_only_contact ? 1U : 0U);
		hash.AppendLittleEndian32(cell.ambiguous ? 1U : 0U);
		hash.AppendLittleEndian64(static_cast<std::uint64_t>(cell.triangle_ids.size()));
		for (const auto triangle : cell.triangle_ids) hash.AppendLittleEndian64(static_cast<std::uint64_t>(triangle));
		const auto element = domain.Background().MaterializeElement(id);
		for (const auto connectivity : element.connectivity) hash.AppendLittleEndian32(static_cast<std::uint32_t>(connectivity));
		for (const auto& row : element.extraction) for (const double value : row) hash.AppendNormalizedDouble(value);
		for (const auto& point : element.bezier_points) for (const double value : point) hash.AppendNormalizedDouble(value);
		const auto& quadrature_cell = catalog.Cell(id);
		hash.AppendLittleEndian32(quadrature_cell.usable ? 1U : 0U);
		hash.AppendLittleEndian32(quadrature_cell.ambiguous ? 1U : 0U);
		hash.AppendNormalizedDouble(quadrature_cell.area_m2);
		hash.AppendLittleEndian64(static_cast<std::uint64_t>(quadrature_cell.point_count));
		const auto& rule = catalog.UsableRule(domain, id); const auto& provenance = catalog.UsableProvenance(domain, id);
		for (std::size_t point = 0; point < rule.Points().size(); ++point) {
			for (const double value : rule.Points()[point].parametric) hash.AppendNormalizedDouble(value);
			for (const double value : rule.Points()[point].physical) hash.AppendNormalizedDouble(value);
			for (const double value : rule.Points()[point].normal) hash.AppendNormalizedDouble(value);
			hash.AppendNormalizedDouble(rule.Points()[point].weight);
			hash.AppendLittleEndian32(rule.Points()[point].boundary_id);
			hash.AppendLittleEndian32(provenance[point].canonical_triangle);
			for (const double value : provenance[point].canonical_barycentric) hash.AppendNormalizedDouble(value);
		}
	}
}

inline std::string BuildStateIdentity(const CartesianDomainClassification& domain,
	const ImmersedSurfaceQuadratureCatalog& catalog, const MaterialSurfaceKinematics& material,
	const MaterialSurfacePatchMap& patch_map,
	double viscosity_pa_s, const std::vector<FluidSurfaceElementState>& state)
{
	Sha256 hash; AppendString(hash, "FluidSurfaceTraction/producer-state/v3");
	AppendCartesianDiscretization(hash, domain, catalog);
	AppendString(hash, material.ContentIdentitySha256());
	AppendString(hash, material.MaterialIdentitySha256());
	AppendString(hash, material.TopologyIdentitySha256());
	AppendString(hash, patch_map.IdentitySha256());
	AppendString(hash, patch_map.ReferenceIdentitySha256());
	hash.AppendLittleEndian32(patch_map.PatchLabel());
	hash.AppendLittleEndian64(static_cast<std::uint64_t>(patch_map.LayoutTriangleToSourceTriangles().size()));
	for (const auto source_triangle : patch_map.LayoutTriangleToSourceTriangles()) hash.AppendLittleEndian32(source_triangle);
	hash.AppendNormalizedDouble(viscosity_pa_s);
	hash.AppendLittleEndian64(static_cast<std::uint64_t>(state.size()));
	for (const auto& element : state) {
		hash.AppendLittleEndian64(element.cell_id);
		hash.AppendLittleEndian64(static_cast<std::uint64_t>(element.nodal_state.size()));
		for (const auto& node : element.nodal_state)
			for (const double value : node) hash.AppendNormalizedDouble(value);
	}
	return hash.Hex();
}

inline std::vector<double> SolveConsistentMass(std::vector<double> matrix,
	const std::vector<double>& right_hand_side, std::size_t node_count)
{
	if (right_hand_side.size() != 3*node_count || matrix.size() != node_count*node_count)
		throw std::invalid_argument("consistent surface mass dimensions are invalid");
	double scale = 0.0;
	for (std::size_t i = 0; i < node_count; ++i) scale = std::max(scale, std::abs(matrix[i*node_count+i]));
	if (!(scale > 0.0) || !std::isfinite(scale)) throw std::runtime_error("consistent surface mass has no selected area");
	for (std::size_t row = 0; row < node_count; ++row) {
		for (std::size_t column = 0; column <= row; ++column) {
			double value = matrix[row*node_count+column];
			for (std::size_t k = 0; k < column; ++k) value -= matrix[row*node_count+k]*matrix[column*node_count+k];
			if (row == column) {
				if (!(value > 256.0*std::numeric_limits<double>::epsilon()*scale) || !std::isfinite(value))
					throw std::runtime_error("consistent surface mass is singular");
				matrix[row*node_count+column] = std::sqrt(value);
			} else matrix[row*node_count+column] = value/matrix[column*node_count+column];
		}
	}
	std::vector<double> result(3*node_count, 0.0);
	for (int component = 0; component < 3; ++component) {
		for (std::size_t row = 0; row < node_count; ++row) {
			double value = right_hand_side[3*row+component];
			for (std::size_t k = 0; k < row; ++k) value -= matrix[row*node_count+k]*result[3*k+component];
			result[3*row+component] = value/matrix[row*node_count+row];
		}
		for (std::size_t reverse = node_count; reverse-- > 0;) {
			double value = result[3*reverse+component];
			for (std::size_t k = reverse+1; k < node_count; ++k) value -= matrix[k*node_count+reverse]*result[3*k+component];
			result[3*reverse+component] = value/matrix[reverse*node_count+reverse];
			if (!std::isfinite(result[3*reverse+component])) throw std::overflow_error("consistent surface mass solution is nonfinite");
		}
	}
	return result;
}

} // namespace fluid_surface_traction_detail

inline std::string BuildFluidSurfaceTractionStateIdentitySha256(
	const CartesianDomainClassification& domain, const ImmersedSurfaceQuadratureCatalog& catalog,
	const MaterialSurfaceKinematics& material, const MaterialSurfacePatchMap& patch_map,
	double dynamic_viscosity_pa_s,
	const std::vector<FluidSurfaceElementState>& state)
{
	return fluid_surface_traction_detail::BuildStateIdentity(domain, catalog, material, patch_map, dynamic_viscosity_pa_s, state);
}

inline FluidSurfaceTractionResult BuildFluidSurfaceTraction(
	const CartesianDomainClassification& domain, const ImmersedSurfaceQuadratureCatalog& catalog,
	const MaterialSurfaceKinematics& material, const DistributedSurfaceInterface& surface,
	const DistributedSurfaceLayout& layout, const MaterialSurfacePatchMap& patch_map,
	const SurfaceFieldStamp& stamp,
	double dynamic_viscosity_pa_s, const std::vector<FluidSurfaceElementState>& state,
	FluidSurfaceTractionProjectionOptions options = {})
{
	using namespace fluid_surface_traction_detail;
	material.Validate(); ValidateDistributedSurfaceInterface(surface); ValidateDistributedSurfaceLayout(layout);
	if (!std::isfinite(dynamic_viscosity_pa_s) || !(dynamic_viscosity_pa_s > 0.0)
		|| options.maximum_nodes == 0 || !std::isfinite(options.conservation_absolute_force_tolerance_n)
		|| !(options.conservation_absolute_force_tolerance_n >= 0.0)
		|| !std::isfinite(options.conservation_absolute_moment_tolerance_n_m)
		|| !(options.conservation_absolute_moment_tolerance_n_m >= 0.0)
		|| !std::isfinite(options.conservation_relative_tolerance)
		|| !(options.conservation_relative_tolerance >= 0.0)
		|| (!(options.conservation_absolute_force_tolerance_n > 0.0)
			&& !(options.conservation_absolute_moment_tolerance_n_m > 0.0)
			&& !(options.conservation_relative_tolerance > 0.0)))
		throw std::invalid_argument("fluid surface traction parameters are invalid");
	if (surface.provides != std::vector<SurfaceFieldQuantity>{SurfaceFieldQuantity::TractionOnStructure}
		|| surface.requires != std::vector<SurfaceFieldQuantity>{
			SurfaceFieldQuantity::Displacement, SurfaceFieldQuantity::Velocity})
		throw std::invalid_argument("fluid surface traction requires the exact traction/displacement/velocity role");
	if (layout.partition_count != 1 || layout.partition_rank != 0)
		throw std::invalid_argument("fluid surface traction supports only the complete single-rank partition");
	if (layout.owned_global_node_ids.size() > options.maximum_nodes)
		throw std::invalid_argument("fluid surface traction node count exceeds the bounded consistent projection");
	if (material.MaterialIdentitySha256() != patch_map.FullReference().MaterialIdentitySha256()
		|| material.TopologyIdentitySha256() != patch_map.FullReference().TopologyIdentitySha256()
		|| material.Surface().CanonicalSha256() != domain.SurfaceCanonicalHash()
		|| catalog.SurfaceCanonicalHash() != domain.SurfaceCanonicalHash())
		throw std::invalid_argument("fluid surface traction material, patch map, or quadrature identity is inconsistent");
	if (surface.reference_mesh_identity_sha256 != patch_map.ReferenceIdentitySha256()
		|| layout.reference_mesh_identity_sha256 != patch_map.ReferenceIdentitySha256()
		|| surface.boundary_labels != patch_map.Interface().boundary_labels
		|| surface.boundary_labels != std::vector<std::int64_t>{static_cast<std::int64_t>(patch_map.PatchLabel())}
		|| layout.layout_identity_sha256 != patch_map.Layout().layout_identity_sha256
		|| BuildDistributedSurfacePartitionIdentitySha256(layout)
			!= BuildDistributedSurfacePartitionIdentitySha256(patch_map.Layout()))
		throw std::invalid_argument("fluid surface traction interface or layout is not the mapped patch authority");
	if (stamp.time_s != material.EvaluatedTimeS())
		throw std::invalid_argument("fluid surface traction stamp time does not match material state");
	ValidateSurfaceFieldStamp(stamp, layout);

	const std::size_t nodes = layout.owned_global_node_ids.size();
	std::vector<std::uint32_t> source_for_local(nodes);
	for (std::size_t local = 0; local < nodes; ++local)
		source_for_local[local] = patch_map.SourceVertexForGlobalNode(layout.owned_global_node_ids[local]);
	std::vector<std::size_t> layout_for_canonical(material.CanonicalTriangleProvenance().size(),
		std::numeric_limits<std::size_t>::max());
	for (std::size_t triangle = 0; triangle < patch_map.LayoutTriangleToSourceTriangles().size(); ++triangle) {
		const auto canonical = patch_map.CanonicalTriangleForLayoutTriangle(triangle);
		if (canonical >= layout_for_canonical.size() || layout_for_canonical[canonical] != std::numeric_limits<std::size_t>::max())
			throw std::invalid_argument("fluid surface traction patch map canonical membership is ambiguous");
		layout_for_canonical[canonical] = triangle;
	}

	std::vector<std::uint64_t> required_cells;
	for (std::uint64_t id = 0; id < domain.Cells().size(); ++id) {
		if (domain.Cells()[static_cast<std::size_t>(id)].classification != CellClassification::Cut) continue;
		const auto& provenance = catalog.UsableProvenance(domain, id);
		if (std::any_of(provenance.begin(), provenance.end(), [&layout_for_canonical](const ImmersedSurfaceQuadraturePointProvenance& point) {
			return point.canonical_triangle < layout_for_canonical.size()
				&& layout_for_canonical[point.canonical_triangle] != std::numeric_limits<std::size_t>::max();
		})) required_cells.push_back(id);
	}
	if (state.size() != required_cells.size()) throw std::invalid_argument("fluid surface traction state does not cover exactly the retained interface cells");
	std::vector<std::array<double, 4>> global_coefficients(domain.Background().NodeCount());
	std::vector<bool> global_coefficient_present(domain.Background().NodeCount(), false);
	for (std::size_t item = 0; item < state.size(); ++item) {
		if (state[item].cell_id != required_cells[item]
			|| (item && state[item-1].cell_id >= state[item].cell_id))
			throw std::invalid_argument("fluid surface traction state cell ownership is ambiguous");
		const auto element = domain.Background().MaterializeElement(state[item].cell_id);
		if (state[item].nodal_state.size() != element.connectivity.size())
			throw std::invalid_argument("fluid surface traction element state does not match IGA element");
		for (std::size_t local = 0; local < state[item].nodal_state.size(); ++local) {
			const auto& value = state[item].nodal_state[local];
			for (double component : value)
				if (!std::isfinite(component)) throw std::invalid_argument("fluid surface traction state is nonfinite");
			const auto global = element.connectivity[local];
			if (global < 0 || static_cast<std::uint64_t>(global) >= domain.Background().NodeCount())
				throw std::runtime_error("fluid surface traction IGA connectivity is invalid");
			const auto node = static_cast<std::size_t>(global);
			if (global_coefficient_present[node]) {
				for (int component = 0; component < 4; ++component)
					if (DoubleBits(global_coefficients[node][component]) != DoubleBits(value[component]))
						throw std::invalid_argument("fluid surface traction shared IGA node coefficients are inconsistent");
			} else { global_coefficients[node] = value; global_coefficient_present[node] = true; }
		}
	}
	if (stamp.producer_state_identity_sha256 != BuildStateIdentity(domain, catalog, material, patch_map, dynamic_viscosity_pa_s, state))
		throw std::invalid_argument("fluid surface traction producer state identity is not derived from the IGA state");

	std::vector<double> mass(nodes*nodes, 0.0); std::vector<CompensatedSum> force_accumulator(3*nodes);
	CompensatedVector quadrature_resultant, quadrature_moment;
	FluidSurfaceTractionDiagnostics diagnostics;
	std::vector<std::array<double, 3>> current_positions(nodes);
	for (std::size_t local = 0; local < nodes; ++local) current_positions[local] = material.SourceVerticesM()[source_for_local[local]];
	for (std::size_t item = 0; item < state.size(); ++item) {
		const auto element = domain.Background().MaterializeElement(state[item].cell_id);
		const auto& rule = catalog.UsableRule(domain, state[item].cell_id);
		const auto& provenance = catalog.UsableProvenance(domain, state[item].cell_id);
		for (std::size_t point_index = 0; point_index < rule.Points().size(); ++point_index) {
			const auto& point = rule.Points()[point_index];
			const auto& material_point = provenance[point_index];
			if (material_point.canonical_triangle >= layout_for_canonical.size()
				|| layout_for_canonical[material_point.canonical_triangle] == std::numeric_limits<std::size_t>::max()) continue;
			const auto layout_triangle = layout_for_canonical[material_point.canonical_triangle];
			const auto& triangle_provenance = material.CanonicalTriangleProvenance().at(material_point.canonical_triangle);
			const auto source_triangle = patch_map.SourceTriangleForLayoutTriangle(layout_triangle);
			if (triangle_provenance.source_triangle != source_triangle
				|| source_triangle >= material.SourceTriangles().size()
				|| point.boundary_id != static_cast<std::int32_t>(patch_map.PatchLabel())
				|| triangle_provenance.boundary_id != patch_map.PatchLabel()
				|| material.SourceTriangles()[source_triangle].boundary_id != patch_map.PatchLabel())
				throw std::invalid_argument("fluid surface traction point label or source provenance disagrees with patch map");
			std::array<std::size_t, 3> local_nodes{};
			for (std::size_t corner = 0; corner < 3; ++corner) {
				const auto source_vertex = triangle_provenance.source_vertex_indices[
					triangle_provenance.canonical_corner_to_source_corner[corner]];
				const auto& map_triangle = layout.reference_triangles[layout_triangle];
				std::size_t layout_corner = 0;
				while (layout_corner < 3 && patch_map.SourceVertexForGlobalNode(map_triangle[layout_corner]) != source_vertex) ++layout_corner;
				if (layout_corner == 3) throw std::runtime_error("fluid surface traction mapped triangle omits canonical source vertex");
				const auto global = map_triangle[layout_corner];
				const auto local = std::lower_bound(layout.owned_global_node_ids.begin(), layout.owned_global_node_ids.end(), global);
				if (local == layout.owned_global_node_ids.end() || *local != global)
					throw std::runtime_error("fluid surface traction mapped triangle has no owned local node");
				local_nodes[corner] = static_cast<std::size_t>(local-layout.owned_global_node_ids.begin());
			}
			const auto basis = EvaluateBasis(element, point.parametric[0], point.parametric[1], point.parametric[2], false);
			double pressure = 0.0; std::array<std::array<double, 3>, 3> gradient{};
			for (std::size_t a = 0; a < state[item].nodal_state.size(); ++a) {
				pressure += state[item].nodal_state[a][3]*basis.value[a];
				for (int velocity = 0; velocity < 3; ++velocity)
					for (int direction = 0; direction < 3; ++direction)
						gradient[velocity][direction] += state[item].nodal_state[a][velocity]*basis.gradient[a][direction];
			}
			const auto traction = FluidOnStructureCauchyTraction(pressure, gradient, dynamic_viscosity_pa_s, point.normal);
			quadrature_resultant.AddScaled(traction, point.weight, "fluid surface traction resultant");
			quadrature_moment.AddScaled(Cross(point.physical, traction), point.weight, "fluid surface traction moment");
			for (std::size_t left = 0; left < 3; ++left) {
				const double shape_left = material_point.canonical_barycentric[left];
				if (!std::isfinite(shape_left)) throw std::runtime_error("fluid surface traction barycentric provenance is nonfinite");
				for (int component = 0; component < 3; ++component) {
					const double addition = traction[component]*shape_left*point.weight;
					if (!std::isfinite(addition)) throw std::overflow_error("fluid surface traction force is nonfinite");
					force_accumulator[3*local_nodes[left]+component].Add(addition, "fluid surface traction force");
				}
				for (std::size_t right = 0; right < 3; ++right)
					mass[local_nodes[left]*nodes+local_nodes[right]] += shape_left*material_point.canonical_barycentric[right]*point.weight;
			}
			++diagnostics.retained_quadrature_points;
		}
	}
	if (!diagnostics.retained_quadrature_points) throw std::runtime_error("fluid surface traction has no retained quadrature points");
	diagnostics.quadrature_resultant_n = quadrature_resultant.Value("fluid surface traction resultant");
	diagnostics.quadrature_moment_n_m = quadrature_moment.Value("fluid surface traction moment");
	std::vector<double> force(3*nodes, 0.0);
	for (std::size_t index = 0; index < force.size(); ++index)
		force[index] = force_accumulator[index].Value("fluid surface traction force");
	const auto nodal_traction = SolveConsistentMass(mass, force, nodes);
	FluidSurfaceTractionResult result;
	result.traction.interface = surface.id; result.traction.stamp = stamp;
	result.diagnostics = diagnostics;
	result.traction.traction_on_structure_pa.resize(nodes); result.traction.consistent_nodal_force_n.resize(nodes);
	CompensatedVector nodal_resultant, nodal_moment;
	for (std::size_t node = 0; node < nodes; ++node) {
		for (int component = 0; component < 3; ++component) {
			result.traction.traction_on_structure_pa[node][component] = nodal_traction[3*node+component];
			result.traction.consistent_nodal_force_n[node][component] = force[3*node+component];
		}
		nodal_resultant.AddScaled(result.traction.consistent_nodal_force_n[node], 1.0, "fluid surface traction nodal resultant");
		nodal_moment.AddScaled(Cross(current_positions[node], result.traction.consistent_nodal_force_n[node]), 1.0, "fluid surface traction nodal moment");
	}
	result.diagnostics.nodal_resultant_n = nodal_resultant.Value("fluid surface traction nodal resultant");
	result.diagnostics.nodal_moment_n_m = nodal_moment.Value("fluid surface traction nodal moment");
	for (int component = 0; component < 3; ++component)
		if (std::abs(diagnostics.quadrature_resultant_n[component]-result.diagnostics.nodal_resultant_n[component]) > options.conservation_absolute_force_tolerance_n
			+options.conservation_relative_tolerance*(quadrature_resultant.absolute_sum[component].Value("fluid surface traction resultant scale")
				+nodal_resultant.absolute_sum[component].Value("fluid surface traction nodal resultant scale"))
			|| std::abs(diagnostics.quadrature_moment_n_m[component]-result.diagnostics.nodal_moment_n_m[component]) > options.conservation_absolute_moment_tolerance_n_m
			+options.conservation_relative_tolerance*(quadrature_moment.absolute_sum[component].Value("fluid surface traction moment scale")
				+nodal_moment.absolute_sum[component].Value("fluid surface traction nodal moment scale")))
			throw std::runtime_error("fluid surface traction consistent projection does not conserve resultant or moment");
	Sha256 identity; AppendString(identity, "FluidSurfaceTractionProjection/consistent-p1/v2");
	AppendString(identity, BuildDistributedSurfaceInterfaceIdentitySha256(surface));
	AppendString(identity, BuildSurfaceFieldStampIdentitySha256(stamp, layout));
	AppendString(identity, patch_map.IdentitySha256()); AppendString(identity, patch_map.ReferenceIdentitySha256());
	AppendString(identity, material.MaterialIdentitySha256());
	AppendString(identity, material.TopologyIdentitySha256()); AppendString(identity, material.ContentIdentitySha256());
	AppendString(identity, catalog.SurfaceCanonicalHash());
	AppendString(identity, BuildStateIdentity(domain, catalog, material, patch_map, dynamic_viscosity_pa_s, state));
	identity.AppendNormalizedDouble(dynamic_viscosity_pa_s); identity.AppendLittleEndian64(static_cast<std::uint64_t>(options.maximum_nodes));
	identity.AppendNormalizedDouble(options.conservation_absolute_force_tolerance_n);
	identity.AppendNormalizedDouble(options.conservation_absolute_moment_tolerance_n_m);
	identity.AppendNormalizedDouble(options.conservation_relative_tolerance);
	identity.AppendLittleEndian32(patch_map.PatchLabel());
	for (const auto source_triangle : patch_map.LayoutTriangleToSourceTriangles()) identity.AppendLittleEndian32(source_triangle);
	for (std::size_t item = 0; item < state.size(); ++item) {
		const auto& rule = catalog.UsableRule(domain, state[item].cell_id); const auto& provenance = catalog.UsableProvenance(domain, state[item].cell_id);
		for (std::size_t point = 0; point < rule.Points().size(); ++point) if (provenance[point].canonical_triangle < layout_for_canonical.size()
			&& layout_for_canonical[provenance[point].canonical_triangle] != std::numeric_limits<std::size_t>::max()) {
			identity.AppendLittleEndian64(state[item].cell_id); identity.AppendLittleEndian32(provenance[point].canonical_triangle);
			for (double value : provenance[point].canonical_barycentric) identity.AppendNormalizedDouble(value);
			for (double value : rule.Points()[point].physical) identity.AppendNormalizedDouble(value);
			for (double value : rule.Points()[point].normal) identity.AppendNormalizedDouble(value);
			identity.AppendNormalizedDouble(rule.Points()[point].weight);
		}
	}
	for (const auto& value : result.traction.consistent_nodal_force_n) AppendVector(identity, value);
	result.traction.projection_identity_sha256 = identity.Hex();
	ValidateSurfaceTraction(result.traction, layout);
	return result;
}

} // namespace iga

#endif
