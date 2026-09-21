#ifndef IGA_GEOMETRY_CONTRACT_MANIFEST_HPP
#define IGA_GEOMETRY_CONTRACT_MANIFEST_HPP

#include "MovingCutGeometry.hpp"

#include <algorithm>
#include <array>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

namespace iga {

inline bool GeometryContractRegionRole(const std::string& role)
{
	return role == "fluid" || role == "solid" || role == "shell"
		|| role == "porous" || role == "network" || role == "lumped";
}

inline bool GeometryContractSameAsReference(const MaterialSurfaceKinematics& surface)
{
	if (surface.ReferenceMaterialVerticesM() != surface.SourceVerticesM()) return false;
	return std::all_of(surface.SourceVertexVelocitiesMPerS().begin(),
		surface.SourceVertexVelocitiesMPerS().end(), [](const std::array<double, 3>& value) {
			return value == std::array<double, 3>{{0.0, 0.0, 0.0}};
		});
}

inline std::string ImmersedGeometryContractJson(const MovingCutGeometry& geometry,
	const std::string& region_role)
{
	if (!GeometryContractRegionRole(region_role))
		throw std::invalid_argument("region role is not part of the geometry contract");
	const auto& kinematics = geometry.Kinematics();
	kinematics.Validate();
	const auto& grid = geometry.Domain().Background().Spec();
	const bool same_as_reference = GeometryContractSameAsReference(kinematics);
	std::ostringstream output;
	output << std::setprecision(17)
		<< "{\n  \"schema_version\": 1,\n"
		<< "  \"route\": \"surface_to_immersed_background\",\n"
		<< "  \"coordinate_system\": \"cartesian\",\n"
		<< "  \"length_unit\": \"m\",\n"
		<< "  \"reference_geometry\": {\"identity_sha256\": \""
		<< kinematics.MaterialIdentitySha256() << "\"},\n"
		<< "  \"current_geometry\": {\"kind\": \""
		<< (same_as_reference ? "same_as_reference" : "evaluated")
		<< "\", \"identity_sha256\": \""
		<< (same_as_reference ? kinematics.MaterialIdentitySha256()
			: kinematics.ContentIdentitySha256())
		<< "\", \"reference_identity_sha256\": \""
		<< kinematics.MaterialIdentitySha256() << "\", \"time_s\": "
		<< kinematics.EvaluatedTimeS() << ", \"canonical_surface_sha256\": \""
		<< kinematics.Surface().CanonicalSha256()
		<< "\", \"runtime_geometry_identity_sha256\": \""
		<< geometry.GeometryIdentitySha256() << "\"},\n"
		<< "  \"background\": {\"lower_m\": [" << grid.lower_m[0] << ", "
		<< grid.lower_m[1] << ", " << grid.lower_m[2] << "], \"upper_m\": ["
		<< grid.upper_m[0] << ", " << grid.upper_m[1] << ", " << grid.upper_m[2]
		<< "], \"cells\": [" << grid.cells[0] << ", " << grid.cells[1] << ", "
		<< grid.cells[2] << "]},\n"
		<< "  \"stable_ids\": {\"surface_vertices\": \"material_vertex_index\", "
		<< "\"surface_triangles\": \"canonical_triangle_index\", "
		<< "\"background_nodes\": \"cartesian_global_node_id\", "
		<< "\"background_cells\": \"cartesian_global_cell_id\"},\n"
		<< "  \"regions\": [{\"id\": \"immersed_volume\", \"role\": \""
		<< region_role << "\", \"dimension\": 3}],\n  \"boundaries\": [";
	std::size_t index = 0;
	for (const auto& entry : geometry.Surface().Diagnostics().source_area_by_boundary_id) {
		if (index++) output << ',';
		output << "{\"label\": " << entry.first << ", \"area_m2\": " << entry.second
			<< ", \"semantic_role\": \"explicitly_configured_elsewhere\"}";
	}
	output << "]\n}\n";
	return output.str();
}

} // namespace iga

#endif
