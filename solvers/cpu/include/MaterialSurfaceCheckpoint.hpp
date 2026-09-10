#ifndef IGA_MATERIAL_SURFACE_CHECKPOINT_HPP
#define IGA_MATERIAL_SURFACE_CHECKPOINT_HPP

#include "MaterialSurfaceKinematics.hpp"
#include "CheckpointMetadataCodec.hpp"

namespace iga {

// Bounded replicated geometry metadata only. Distributed fluid fields belong
// in separate owned shards, and atomic bundle publication is a caller concern.
struct MaterialSurfaceCheckpointLimits {
	std::size_t maximum_vertices=100000,maximum_triangles=200000;
};

inline std::string SerializeMaterialSurfaceCheckpoint(const MaterialSurfaceKinematics& material,
	MaterialSurfaceCheckpointLimits limits={})
{
	using checkpoint_metadata::Require;
	material.Validate();
	Require(material.SourceVerticesM().size()<=limits.maximum_vertices
		&& material.SourceTriangles().size()<=limits.maximum_triangles,"material geometry exceeds checkpoint limits");
	checkpoint_metadata::Writer output;
	output.Text("IGA_MATERIAL_SURFACE/1");output.Text(material.ContentIdentitySha256());
	output.Real(material.EvaluatedTimeS());output.Real(material.StepStartS());output.Real(material.StepEndS());
	output.Unsigned(material.SourceVerticesM().size());output.Unsigned(material.SourceTriangles().size());
	for(const auto* vertices:{&material.ReferenceMaterialVerticesM(),&material.SourceVerticesM(),&material.SourceVertexVelocitiesMPerS()})
		for(const auto& vertex:*vertices)for(double value:vertex)output.Real(value);
	for(const auto& triangle:material.SourceTriangles()) {
		for(auto node:triangle.source_vertex_indices)output.Unsigned(node);
		output.Unsigned(triangle.boundary_id);
	}
	return output.Bytes();
}

inline MaterialSurfaceKinematics ParseMaterialSurfaceCheckpoint(std::string_view bytes,
	const std::string& expected_content_identity,MaterialSurfaceCheckpointLimits limits={})
{
	using checkpoint_metadata::Require;
	checkpoint_metadata::Reader input(bytes);
	Require(input.Text()=="IGA_MATERIAL_SURFACE/1","unsupported material checkpoint version");
	const auto identity=input.Text();
	Require(identity==expected_content_identity,"material checkpoint differs from expected epoch");
	const double time=input.Real(),start=input.Real(),end=input.Real();
	const auto vertices=input.Unsigned(),triangles=input.Unsigned();
	Require(vertices>0&&triangles>0&&vertices<=limits.maximum_vertices&&triangles<=limits.maximum_triangles
		&&vertices<=UINT32_MAX&&triangles<=UINT32_MAX,"invalid material checkpoint geometry counts");
	// Check a lower bound before allocation; a tiny corrupt stream cannot
	// cause allocations proportional to a forged count, even with high caps.
	Require(vertices<=bytes.size()/72&&triangles<=(bytes.size()-72*vertices)/32,
		"truncated material geometry payload");
	std::vector<std::array<double,3>> reference(vertices),positions(vertices),velocities(vertices);
	for(auto* values:{&reference,&positions,&velocities})
		for(auto& vertex:*values)for(auto& value:vertex)value=input.Real();
	std::vector<RawSurfaceTriangle> topology(triangles);
	for(auto& triangle:topology) {
		for(auto& index:triangle.indices) {
			const auto node=input.Unsigned();Require(node<vertices,"material checkpoint connectivity is outside vertex array");
			index=static_cast<std::int64_t>(node);
		}
		const auto label=input.Unsigned();Require(label<=UINT32_MAX,"material checkpoint label overflow");
		triangle.boundary_id=static_cast<std::int64_t>(label);
	}
	input.Finish();
	auto result=MaterialSurfaceKinematics::CreateFromSourceTopology(std::move(reference),std::move(positions),
		std::move(velocities),std::move(topology),time,start,end);
	Require(result.ContentIdentitySha256()==identity,"material checkpoint content identity mismatch");
	return result;
}

} // namespace iga
#endif
