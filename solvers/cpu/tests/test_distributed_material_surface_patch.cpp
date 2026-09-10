#include "DistributedMaterialSurfacePatchKinematics.hpp"
#define main SerialMaterialPatchRegression
#include "test_material_surface_patch_kinematics.cpp"
#undef main
#include <iostream>

namespace {
template<class Function> void RejectCollectively(MPI_Comm comm,Function&& function)
{
	bool rejected=false;try { function(); } catch(const std::exception&) { rejected=true; }
	iga::CollectiveLocalStage(comm,"expected distributed patch rejection",[&] {
		if(!rejected)throw std::runtime_error("missing distributed patch rejection");
	});
}

void Run(MPI_Comm comm)
{
	int rank=0,size=0;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
	const auto full=Full();const auto serial_map=Map(full);auto layout=Layout(full);
	layout.partition_count=size;layout.partition_rank=rank;
	const auto all=layout.owned_global_node_ids;
	layout.owned_global_node_ids.clear();layout.owned_reference_lumped_areas_m2.clear();
	for(std::size_t i=0;i<all.size();++i)if((size==1 ? 0 : 1+static_cast<int>(i)%(size-1))==rank) {
		layout.owned_global_node_ids.push_back(all[i]);layout.owned_reference_lumped_areas_m2.push_back(.1);
	}
	layout.layout_identity_sha256=iga::BuildDistributedSurfaceLayoutIdentitySha256(layout);
	const auto map=iga::MaterialSurfacePatchMap::Create(Interface(layout.reference_mesh_identity_sha256),layout,
		full,7,MappingData(),TriangleData(),Clamps());
	iga::FsiTrialContext context;context.step=1;context.dt_s=.5;context.start_time_s=0;context.coupling_iteration=2;
	const auto serial_input=Trial(serial_map,context);auto input=Trial(map,context);
	input.displacement_m.clear();input.velocity_m_per_s.clear();
	for(auto id:layout.owned_global_node_ids) {
		const auto row=std::lower_bound(all.begin(),all.end(),id)-all.begin();
		input.displacement_m.push_back(serial_input.displacement_m[row]);input.velocity_m_per_s.push_back(serial_input.velocity_m_per_s[row]);
	}
	iga::Sha256 producer;producer.AppendLittleEndian64(rank);input.stamp.producer_state_identity_sha256=producer.Hex();
	const auto expected=input.stamp;const auto input_hash=iga::BuildSurfaceKinematicsIdentitySha256(input,layout);
	const auto first=iga::ComposeDistributedMaterialSurfacePatch(comm,map,full,input,expected,full.ContentIdentitySha256(),context);
	const auto serial_first=iga::MaterialSurfacePatchKinematics::ComposeTarget(serial_map,full,serial_input,context);
	assert(first.target.ContentIdentitySha256()==serial_first.target.ContentIdentitySha256());
	assert(first.target.IdentitySha256()==serial_first.target.IdentitySha256());
	iga::RequireCollectiveSameText(comm,"patch composition result identity",first.composition_identity_sha256);
	assert(input_hash==iga::BuildSurfaceKinematicsIdentitySha256(input,layout));
	RejectCollectively(comm,[&] {
		auto wrong=expected;if(rank==0)wrong.time_s+=.01;
		(void)iga::ComposeDistributedMaterialSurfacePatch(comm,map,full,input,wrong,full.ContentIdentitySha256(),context);
	});
	RejectCollectively(comm,[&] {
		(void)iga::ComposeDistributedMaterialSurfacePatch(comm,map,full,input,expected,
			rank==0 ? std::string(64,'0') : full.ContentIdentitySha256(),context);
	});
	for(bool seam:{false,true})RejectCollectively(comm,[&] {
		auto wrong=input;
		for(std::size_t row=0;row<layout.owned_global_node_ids.size();++row) {
			if(seam && layout.owned_global_node_ids[row]==10)wrong.displacement_m[row][2]=.01;
			if(!seam && layout.owned_global_node_ids[row]==14)wrong.velocity_m_per_s[row][2]+=.01;
		}
		(void)iga::ComposeDistributedMaterialSurfacePatch(comm,map,full,wrong,expected,full.ContentIdentitySha256(),context);
	});
	RejectCollectively(comm,[&] {
		(void)iga::ComposeDistributedMaterialSurfacePatch(comm,map,full,input,expected,full.ContentIdentitySha256(),context,rank==0 ? 3 : 4096);
	});
	RejectCollectively(comm,[&] {
		auto wrong=input;
		for(std::size_t row=0;row<layout.owned_global_node_ids.size();++row)
			if(layout.owned_global_node_ids[row]==14)wrong.velocity_m_per_s[row][2]=std::numeric_limits<double>::quiet_NaN();
		(void)iga::ComposeDistributedMaterialSurfacePatch(comm,map,full,wrong,expected,full.ContentIdentitySha256(),context);
	});
	RejectCollectively(comm,[&] {
		std::optional<iga::MaterialSurfacePatchMap> missing_map;auto missing_input=input;
		iga::CollectiveLocalStage(comm,"missing patch owner fixture",[&] {
			auto missing_layout=layout;
			if(rank==(size==1 ? 0 : 1)) {
				missing_layout.owned_global_node_ids.erase(missing_layout.owned_global_node_ids.begin());
				missing_layout.owned_reference_lumped_areas_m2.erase(missing_layout.owned_reference_lumped_areas_m2.begin());
				missing_input.displacement_m.erase(missing_input.displacement_m.begin());
				missing_input.velocity_m_per_s.erase(missing_input.velocity_m_per_s.begin());
			}
			missing_map.emplace(iga::MaterialSurfacePatchMap::Create(Interface(missing_layout.reference_mesh_identity_sha256),
				missing_layout,full,7,MappingData(),TriangleData(),Clamps()));
			missing_input.stamp.partition_identity_sha256=iga::BuildDistributedSurfacePartitionIdentitySha256(missing_layout);
		});
		(void)iga::ComposeDistributedMaterialSurfacePatch(comm,*missing_map,full,missing_input,missing_input.stamp,full.ContentIdentitySha256(),context);
	});
	const auto retry=iga::ComposeDistributedMaterialSurfacePatch(comm,map,full,input,expected,full.ContentIdentitySha256(),context);
	assert(retry.composition_identity_sha256==first.composition_identity_sha256);
	context.step=2;context.start_time_s=.5;context.coupling_iteration=0;
	auto second_input=input;second_input.stamp.time_s=context.EndTime();second_input.stamp.step=2;second_input.stamp.coupling_iteration=0;
	for(std::size_t row=0;row<layout.owned_global_node_ids.size();++row)if(layout.owned_global_node_ids[row]==14) {
		second_input.displacement_m[row][2]=.15;second_input.velocity_m_per_s[row][2]=.1;
	}
	auto serial_second_input=Trial(serial_map,context);serial_second_input.displacement_m[3][2]=.15;serial_second_input.velocity_m_per_s[3][2]=.1;
	const auto second=iga::ComposeDistributedMaterialSurfacePatch(comm,map,first.target,second_input,second_input.stamp,
		first.target.ContentIdentitySha256(),context);
	const auto serial_second=iga::MaterialSurfacePatchKinematics::ComposeTarget(serial_map,serial_first.target,serial_second_input,context);
	assert(second.target.ContentIdentitySha256()==serial_second.target.ContentIdentitySha256());
	assert(second.target.MaterialIdentitySha256()==full.MaterialIdentitySha256());
	assert(second.target.TopologyIdentitySha256()==full.TopologyIdentitySha256());
	std::cout<<"distributed_material_patch rank="<<rank<<" ranks="<<size<<" owned="<<layout.owned_global_node_ids.size()
		<<" steps=2 serial_identity_equal=1 unchanged_input=1 rejection_retry=1 passed\n";
}
}

int main(int argc,char** argv)
{
	MPI_Init(&argc,&argv);int status=0;
	try { Run(MPI_COMM_WORLD); } catch(const std::exception& error) { std::cerr<<error.what()<<'\n';status=1; }
	MPI_Finalize();return status;
}
