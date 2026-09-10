#include "MovingCheckpointBundle.hpp"
#include "CompliantChannelFsiFixture.hpp"
#include <fstream>
#define main MovingMetadataFixtureMain
#include "test_moving_checkpoint_metadata.cpp"
#undef main

int main(int argc,char** argv)
{
	try {
		Require(argc==2);const std::filesystem::path root(argv[1]);Require(std::filesystem::create_directory(root));
		auto state=State();state.time_s=2.;state.conservation.source_time_s=1.;state.conservation.target_time_s=2.;state.conservation.dt_s=1.;
		const auto previous=iga::compliant_channel_fixture::InitialMaterial();
		std::vector<iga::RawSurfaceTriangle> topology;
		for(const auto& item:previous.SourceTriangles()) {
			iga::RawSurfaceTriangle triangle;triangle.boundary_id=item.boundary_id;
			for(int i=0;i<3;++i)triangle.indices[i]=item.source_vertex_indices[i];
			topology.push_back(triangle);
		}
		const auto current=iga::MaterialSurfaceKinematics::CreateFromSourceTopology(previous.ReferenceMaterialVerticesM(),
			previous.SourceVerticesM(),previous.SourceVertexVelocitiesMPerS(),topology,2.,1.,2.);
		state.previous_material.emplace(previous);state.current_material.emplace(current);state.material_identity_sha256=current.ContentIdentitySha256();
		const std::string config(64,'1');const auto hash=[](const std::string& text) { return iga::coupled_checkpoint_detail::Hash(text); };
		const iga::CoupledCheckpointEpoch epoch{hash("moving-bundle"),{}, {hash("case"),config,hash("mapping"),3},2,2.,1.};
		const iga::MovingCheckpointLayout layout{"flow",{2,0,1}};
		iga::CreateCoupledCheckpointEpoch(root,epoch);std::vector<iga::CoupledCheckpointShard> receipts;
		for(std::size_t rank=0;rank<3;++rank) {
			auto source=state;source.field={3,rank,rank+1,{rank+1.25}};
			auto written=iga::WriteMovingCheckpointShards(root,epoch,layout,rank,source,config);
			receipts.insert(receipts.end(),written.begin(),written.end());
		}
		std::sort(receipts.begin(),receipts.end(),[](const auto& a,const auto& b) { return a.spec.id<b.spec.id; });
		const auto manifest=iga::PublishCoupledCheckpoint(root,epoch,iga::MovingCheckpointCatalog(layout,3),receipts);
		const auto materials=iga::LoadMovingCheckpointMaterials(root,manifest,layout,config);
		Require(iga::LoadMovingCheckpointPortControls(root,manifest,layout,config)==state.port_control_values);
		Require(materials.previous.ContentIdentitySha256()==previous.ContentIdentitySha256());
		Require(materials.current.ContentIdentitySha256()==current.ContentIdentitySha256());
		for(std::uint64_t ranks:{1u,2u,4u})for(std::uint64_t rank=0;rank<ranks;++rank) {
			auto target=state;const auto begin=3*rank/ranks,end=3*(rank+1)/ranks;
			target.field={3,begin,end,std::vector<double>(end-begin,-1.)};target.conservation={};
			iga::LoadMovingCheckpointFields(root,manifest,layout,config,target);
			for(std::size_t i=0;i<target.field.values.size();++i)Require(target.field.values[i]==begin+i+1.25);
			Require(iga::SerializeMovingConservationCheckpoint(target.conservation)==iga::SerializeMovingConservationCheckpoint(state.conservation));
		}
		Reject([&] { (void)iga::LoadMovingCheckpointMaterials(root,manifest,layout,std::string(64,'2')); });
		auto missing=manifest;missing.shards.pop_back();
		Reject([&] { (void)iga::LoadMovingCheckpointMaterials(root,missing,layout,config); });
		const auto field=root/epoch.id/"flow.field.rank-1.shard";
		std::ifstream input(field,std::ios::binary);const std::string original((std::istreambuf_iterator<char>(input)),{});input.close();
		auto damaged=original;damaged.back()^=1;
		const auto replace=[&](const std::string& bytes) { std::ofstream output(field,std::ios::binary|std::ios::trunc);output.write(bytes.data(),bytes.size());output.close();Require(bool(output)); };
		replace(damaged);
		Reject([&] { auto target=state;target.field={3,0,1,{-1.}};iga::LoadMovingCheckpointFields(root,manifest,layout,config,target); });
		replace(original);
		auto retry=state;retry.field={3,0,1,{-1.}};iga::LoadMovingCheckpointFields(root,manifest,layout,config,retry);Require(retry.field.values[0]==1.25);
		std::cout<<"moving checkpoint real-file bundle: materials, configuration, source mapping, 3-to-1/2/4 repartition, corruption and retry passed\n";
	} catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
	return 0;
}
