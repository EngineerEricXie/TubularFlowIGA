#include "MovingFsiPublicationCheckpoint.hpp"
#define main MembraneRegressionMain
#include "test_pretensioned_membrane.cpp"
#undef main

int main()
{
	try {
		const auto layout=SquareLayout(2);const auto configuration=HashText("configuration");
		iga::MovingFsiPublicationCheckpoint state;state.context={3,.2,.1,7};
		state.material_identity=HashText("material");state.geometry_identity=HashText("geometry");
		state.composition_identity=HashText("composition");
		auto& value=state.traction;value.interface=FluidInterface();
		value.stamp={state.context.EndTime(),3,7,layout.reference_mesh_identity_sha256,layout.layout_identity_sha256,
			iga::BuildDistributedSurfacePartitionIdentitySha256(layout),HashText("producer")};
		value.projection_identity_sha256=HashText("projection");
		for(auto id:layout.owned_global_node_ids) {
			value.traction_on_structure_pa.push_back({{double(id),-double(id),.3}});
			value.consistent_nodal_force_n.push_back({{.5,-.2,double(id)*.01}});
		}
		const auto hash=[](std::string_view bytes) { iga::Sha256 digest;digest.Append(bytes.data(),bytes.size());return digest.Hex(); };
		const auto bytes=iga::SerializeMovingFsiPublicationCheckpoint(state,layout,configuration);
		const auto parse=[&](std::string_view input,const std::string& authority) {
			return iga::ParseMovingFsiPublicationCheckpoint(input,authority,layout,FluidInterface(),configuration,
				state.material_identity,state.geometry_identity,3,state.context.EndTime());
		};
		const auto restored=parse(bytes,hash(bytes));
		assert(iga::SerializeMovingFsiPublicationCheckpoint(restored,layout,configuration)==bytes);
		assert(iga::BuildSurfaceTractionIdentitySha256(restored.traction,layout)==iga::BuildSurfaceTractionIdentitySha256(value,layout));
		auto corrupt=bytes;corrupt.back()^=1;Reject([&] { (void)parse(corrupt,hash(bytes)); });
		for(const auto& malformed:{bytes.substr(0,bytes.size()-1),bytes+"x"})Reject([&] { (void)parse(malformed,hash(malformed)); });
		for(int mismatch=0;mismatch<5;++mismatch) {
			auto changed=state;
			if(mismatch==0)changed.material_identity=HashText("wrong");
			if(mismatch==1)changed.geometry_identity=HashText("wrong");
			if(mismatch==2)changed.traction.interface.domain_id="wrong";
			if(mismatch==3) { ++changed.context.step;++changed.traction.stamp.step; }
			const auto wrong=iga::SerializeMovingFsiPublicationCheckpoint(changed,layout,mismatch==4?HashText("wrong"):configuration);
			Reject([&] { (void)parse(wrong,hash(wrong)); });
		}
		auto changed=state;changed.traction.traction_on_structure_pa[0][0]=std::numeric_limits<double>::infinity();
		Reject([&] { (void)iga::SerializeMovingFsiPublicationCheckpoint(changed,layout,configuration); });
		Reject([&] { (void)iga::SerializeMovingFsiPublicationCheckpoint(state,layout,configuration,1); });
		Reject([&] { (void)iga::ParseMovingFsiPublicationCheckpoint(bytes,hash(bytes),layout,FluidInterface(),configuration,
			state.material_identity,state.geometry_identity,3,state.context.EndTime(),1); });
		auto different_partition=layout;different_partition.partition_count=2;different_partition.partition_rank=1;
		different_partition.layout_identity_sha256=iga::BuildDistributedSurfaceLayoutIdentitySha256(different_partition);
		Reject([&] { (void)iga::ParseMovingFsiPublicationCheckpoint(bytes,hash(bytes),different_partition,FluidInterface(),configuration,
			state.material_identity,state.geometry_identity,3,state.context.EndTime()); });
		auto nonfinite=bytes;const double infinity=std::numeric_limits<double>::infinity();std::uint64_t bits;
		std::memcpy(&bits,&infinity,sizeof(bits));
		for(unsigned byte=0;byte<8;++byte)nonfinite[nonfinite.size()-8+byte]=static_cast<char>((bits>>(8*byte))&255);
		Reject([&] { (void)parse(nonfinite,hash(nonfinite)); });
		// Empty owned slices are legitimate and still carry authenticated stamps.
		auto empty=layout;empty.partition_count=2;empty.partition_rank=0;
		empty.owned_global_node_ids.clear();empty.owned_reference_lumped_areas_m2.clear();
		empty.layout_identity_sha256=iga::BuildDistributedSurfaceLayoutIdentitySha256(empty);
		changed=state;changed.traction.traction_on_structure_pa.clear();changed.traction.consistent_nodal_force_n.clear();
		changed.traction.stamp.layout_identity_sha256=empty.layout_identity_sha256;
		changed.traction.stamp.partition_identity_sha256=iga::BuildDistributedSurfacePartitionIdentitySha256(empty);
		const auto empty_bytes=iga::SerializeMovingFsiPublicationCheckpoint(changed,empty,configuration);
		const auto empty_state=iga::ParseMovingFsiPublicationCheckpoint(empty_bytes,hash(empty_bytes),empty,FluidInterface(),configuration,
			state.material_identity,state.geometry_identity,3,state.context.EndTime());
		assert(empty_state.traction.traction_on_structure_pa.empty());
		Reject([&] { (void)parse(empty_bytes,hash(empty_bytes)); });
		std::cout<<"moving FSI publication checkpoint guards passed\n";return 0;
	} catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
}
