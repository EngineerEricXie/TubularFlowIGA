#ifndef IGA_MOVING_FSI_PUBLICATION_CHECKPOINT_HPP
#define IGA_MOVING_FSI_PUBLICATION_CHECKPOINT_HPP

#include "CheckpointMetadataCodec.hpp"
#include "FsiDomainRuntime.hpp"

namespace iga {

struct MovingFsiPublicationCheckpoint {
	FsiTrialContext context;
	SurfaceTraction traction;
	std::string composition_identity,material_identity,geometry_identity;
};

struct MovingFsiPublicationRecords {
	MovingFsiPublicationCheckpoint state;
	std::vector<std::uint64_t> node_ids;
	std::string reference_identity,layout_identity,partition_identity,producer_identity,projection_identity;
};

// One owned surface slice. The outer bundle authenticates bytes and binds the
// physical configuration. Repartitioning must explicitly reconstruct slices;
// the parser never relabels an old publication with a new partition identity.
inline std::string SerializeMovingFsiPublicationCheckpoint(const MovingFsiPublicationCheckpoint& state,
	const DistributedSurfaceLayout& layout,const std::string& configuration,std::size_t maximum_nodes=4096)
{
	using checkpoint_metadata::Require;
	Require(maximum_nodes>0&&layout.global_node_count<=maximum_nodes,"FSI surface exceeds checkpoint bound");
	ValidateFsiTrialContext(state.context);ValidateSurfaceTraction(state.traction,layout);
	Require(state.context.step>0&&state.context.start_time_s>=0,"invalid accepted FSI context");
	const auto& stamp=state.traction.stamp;
	Require(stamp.step==state.context.step&&stamp.time_s==state.context.EndTime()
		&&stamp.coupling_iteration==state.context.coupling_iteration,"FSI publication context differs");
	checkpoint_metadata::Writer out;out.Text("IGA_MOVING_FSI_PUBLICATION/1");
	for(const auto* id:{&configuration,&state.composition_identity,&state.material_identity,&state.geometry_identity}) {
		Require(IsLowercaseSha256(*id),"invalid FSI checkpoint identity");out.Text(*id);
	}
	out.Unsigned(state.context.step);out.Real(state.context.start_time_s);out.Real(state.context.dt_s);
	out.Unsigned(state.context.coupling_iteration);
	for(const auto* text:{&state.traction.interface.domain_id,&state.traction.interface.subsystem_id,&state.traction.interface.interface_id,
		&stamp.reference_mesh_identity_sha256,&stamp.layout_identity_sha256,&stamp.partition_identity_sha256,
		&stamp.producer_state_identity_sha256,&state.traction.projection_identity_sha256})out.Text(*text);
	out.Unsigned(layout.owned_global_node_ids.size());
	for(std::size_t row=0;row<layout.owned_global_node_ids.size();++row) {
		out.Unsigned(layout.owned_global_node_ids[row]);
		for(double value:state.traction.traction_on_structure_pa[row])out.Real(value);
		for(double value:state.traction.consistent_nodal_force_n[row])out.Real(value);
	}
	return out.Bytes();
}

inline MovingFsiPublicationCheckpoint ParseMovingFsiPublicationCheckpoint(std::string_view bytes,
	const std::string& expected_payload_identity,const DistributedSurfaceLayout& layout,
	const SurfaceInterfaceRef& interface,const std::string& configuration,
	const std::string& material_identity,const std::string& geometry_identity,
	std::uint64_t step,double time_s,std::size_t maximum_nodes=4096)
{
	using checkpoint_metadata::Require;
	Require(bytes.size()<=checkpoint_metadata::maximum_bytes&&IsLowercaseSha256(expected_payload_identity),"invalid FSI payload authority");
	Sha256 hash;hash.Append(bytes.data(),bytes.size());
	Require(hash.Hex()==expected_payload_identity,"FSI checkpoint payload hash differs");
	Require(maximum_nodes>0&&layout.global_node_count<=maximum_nodes,"FSI surface exceeds checkpoint bound");
	checkpoint_metadata::Reader in(bytes);MovingFsiPublicationCheckpoint result;
	Require(in.Text()=="IGA_MOVING_FSI_PUBLICATION/1","unsupported FSI publication checkpoint");
	Require(in.Text()==configuration,"FSI configuration differs");
	result.composition_identity=in.Text();result.material_identity=in.Text();result.geometry_identity=in.Text();
	Require(result.material_identity==material_identity&&result.geometry_identity==geometry_identity,"FSI checkpoint geometry differs");
	result.context.step=in.Unsigned();result.context.start_time_s=in.Real();result.context.dt_s=in.Real();
	result.context.coupling_iteration=in.Unsigned();
	Require(result.context.step==step&&result.context.EndTime()==time_s,"FSI accepted clock differs");
	auto& value=result.traction;auto& stamp=value.stamp;
	value.interface.domain_id=in.Text();value.interface.subsystem_id=in.Text();value.interface.interface_id=in.Text();
	Require(value.interface==interface,"FSI checkpoint interface differs");
	stamp.reference_mesh_identity_sha256=in.Text();stamp.layout_identity_sha256=in.Text();stamp.partition_identity_sha256=in.Text();
	stamp.producer_state_identity_sha256=in.Text();value.projection_identity_sha256=in.Text();
	stamp.step=result.context.step;stamp.time_s=result.context.EndTime();stamp.coupling_iteration=result.context.coupling_iteration;
	const auto count=in.Unsigned();Require(count==layout.owned_global_node_ids.size()&&count<=maximum_nodes,"FSI owned slice differs");
	value.traction_on_structure_pa.resize(count);value.consistent_nodal_force_n.resize(count);
	for(std::size_t row=0;row<count;++row) {
		Require(in.Unsigned()==layout.owned_global_node_ids[row],"FSI checkpoint owned node differs");
		for(double& component:value.traction_on_structure_pa[row])component=in.Real();
		for(double& component:value.consistent_nodal_force_n[row])component=in.Real();
	}
	in.Finish();
	Require(SerializeMovingFsiPublicationCheckpoint(result,layout,configuration,maximum_nodes)==bytes,"noncanonical FSI publication checkpoint");
	return result;
}

// Decode an authenticated source slice without pretending that its saved
// partition stamp belongs to the destination layout. Global coverage and the
// new publication provenance are established by the repartitioning factory.
inline MovingFsiPublicationRecords DecodeMovingFsiPublicationRecords(std::string_view bytes,
	const std::string& expected_payload_identity,const SurfaceInterfaceRef& interface,
	const std::string& configuration,const std::string& material_identity,
	const std::string& geometry_identity,std::uint64_t step,double time_s,std::size_t maximum_nodes=4096)
{
	using checkpoint_metadata::Require;
	Require(bytes.size()<=checkpoint_metadata::maximum_bytes&&IsLowercaseSha256(expected_payload_identity),"invalid FSI payload authority");
	Sha256 hash;hash.Append(bytes.data(),bytes.size());Require(hash.Hex()==expected_payload_identity,"FSI checkpoint payload hash differs");
	checkpoint_metadata::Reader in(bytes);MovingFsiPublicationRecords records;auto& result=records.state;
	Require(in.Text()=="IGA_MOVING_FSI_PUBLICATION/1","unsupported FSI publication checkpoint");
	Require(in.Text()==configuration,"FSI configuration differs");
	result.composition_identity=in.Text();result.material_identity=in.Text();result.geometry_identity=in.Text();
	Require(IsLowercaseSha256(result.composition_identity)&&result.material_identity==material_identity
		&&result.geometry_identity==geometry_identity,"FSI checkpoint geometry differs");
	result.context.step=in.Unsigned();result.context.start_time_s=in.Real();result.context.dt_s=in.Real();
	result.context.coupling_iteration=in.Unsigned();ValidateFsiTrialContext(result.context);
	Require(result.context.step==step&&result.context.EndTime()==time_s,"FSI accepted clock differs");
	auto& value=result.traction;value.interface.domain_id=in.Text();value.interface.subsystem_id=in.Text();value.interface.interface_id=in.Text();
	Require(value.interface==interface,"FSI checkpoint interface differs");
	records.reference_identity=in.Text();records.layout_identity=in.Text();records.partition_identity=in.Text();
	records.producer_identity=in.Text();records.projection_identity=in.Text();
	for(const auto* identity:{&records.reference_identity,&records.layout_identity,&records.partition_identity,
		&records.producer_identity,&records.projection_identity})Require(IsLowercaseSha256(*identity),"invalid FSI source publication identity");
	const auto count=in.Unsigned();Require(count<=maximum_nodes,"FSI owned slice exceeds bound");
	records.node_ids.resize(count);value.traction_on_structure_pa.resize(count);value.consistent_nodal_force_n.resize(count);
	for(std::size_t row=0;row<count;++row) {
		records.node_ids[row]=in.Unsigned();
		Require(row==0||records.node_ids[row-1]<records.node_ids[row],"FSI source node IDs must be sorted and unique");
		for(double& component:value.traction_on_structure_pa[row])component=in.Real();
		for(double& component:value.consistent_nodal_force_n[row])component=in.Real();
	}
	in.Finish();return records;
}

} // namespace iga
#endif
