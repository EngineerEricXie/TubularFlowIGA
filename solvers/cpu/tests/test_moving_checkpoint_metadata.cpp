#include "MovingCheckpointRestore.hpp"
#include <iostream>

namespace {
void Require(bool condition)
{
	if(!condition)throw std::runtime_error("moving checkpoint metadata test failed");
}
template<class F> void Reject(F action)
{
	bool rejected=false;try { action(); }catch(const std::exception&) { rejected=true; }
	Require(rejected);
}
iga::ImmersedMovingAcceptedCheckpoint State()
{
	iga::ImmersedMovingAcceptedCheckpoint state;
	state.geometry_identity_sha256=std::string(64,'a');state.publication_identity_sha256=std::string(64,'b');
	state.layout_identity_sha256=std::string(64,'c');state.material_identity_sha256=std::string(64,'d');
	state.step=2;state.time_s=.25;state.field={3,1,3,{2.3,-4.}};
	state.port_control_values={{"inlet",.1},{"outlet",0.}};
	auto& c=state.conservation;c.source_index=1;c.target_index=2;c.source_time_s=.125;c.target_time_s=.25;c.dt_s=.125;
	c.source_geometry_identity_sha256=std::string(64,'e');c.target_geometry_identity_sha256=state.geometry_identity_sha256;
	c.source_publication_identity_sha256=std::string(64,'f');c.target_publication_identity_sha256=state.publication_identity_sha256;
	c.source_audited_volume_m3=c.target_audited_volume_m3=1.;
	c.endpoint.surface_flow_by_boundary_label_m3_s={{7,-.2},{8,.2}};
	return state;
}
}
int main()
{
	try {
		const auto state=State();const std::string config(64,'1');
		const auto bytes=iga::SerializeMovingAcceptedCheckpointMetadata(state,config);
		const auto authority=iga::ParseMovingAcceptedCheckpointMetadata(bytes,state,config);
		const auto conservation=iga::SerializeMovingConservationCheckpoint(state.conservation);
		const auto restored=iga::ParseMovingConservationCheckpoint(conservation,authority,state.geometry_identity_sha256,
			state.publication_identity_sha256,state.step,state.time_s);
		Require(iga::SerializeMovingConservationCheckpoint(restored)==conservation);
		Reject([&] { (void)iga::ParseMovingAcceptedCheckpointMetadata(bytes,state,std::string(64,'2')); });
		for(int mode=0;mode<8;++mode) {
			auto wrong=state;
			if(mode==0)wrong.geometry_identity_sha256=std::string(64,'0');
			if(mode==1)wrong.publication_identity_sha256=std::string(64,'0');
			if(mode==2)wrong.layout_identity_sha256=std::string(64,'0');
			if(mode==3)wrong.material_identity_sha256=std::string(64,'0');
			if(mode==4)wrong.step++;
			if(mode==5)wrong.time_s+=.125;
			if(mode==6)wrong.field.global_rows++;
			if(mode==7)wrong.port_control_values["inlet"]+=.1;
			Reject([&] { (void)iga::ParseMovingAcceptedCheckpointMetadata(bytes,wrong,config); });
		}
		for(const auto& damaged:{bytes.substr(0,bytes.size()-1),bytes+"x"})
			Reject([&] { (void)iga::ParseMovingAcceptedCheckpointMetadata(damaged,state,config); });
		iga::ImmersedTransientFlowOptions options;
		iga::ImmersedFlowPortDefinition inlet,outlet;inlet.id="inlet";inlet.boundary_label=8;inlet.value=-1.;
		outlet.id="outlet";outlet.boundary_label=9;outlet.value=4.;options.ports={inlet,outlet};
		iga::ApplyMovingCheckpointPortControls(options,state.port_control_values);
		Require(options.ports[0].value==.1&&options.ports[1].value==0.&&options.ports[0].boundary_label==8&&options.ports[1].boundary_label==9);
		for(int mode=0;mode<3;++mode) {
			auto invalid=state.port_control_values;
			if(mode==0)invalid.erase("outlet");
			if(mode==1) { invalid.erase("outlet");invalid["unknown"]=1.; }
			if(mode==2)invalid["outlet"]=std::numeric_limits<double>::infinity();
			Reject([&] { iga::ApplyMovingCheckpointPortControls(options,invalid); });
			Require(options.ports[0].value==.1&&options.ports[1].value==0.);
		}
		iga::checkpoint_metadata::Writer controls;controls.Reals(state.port_control_values);
		auto legacy=bytes.substr(0,bytes.size()-controls.Bytes().size());legacy[8+std::string("IGA_MOVING_ACCEPTED/2").size()-1]='1';
		Reject([&] { (void)iga::ParseMovingAcceptedCheckpointMetadata(legacy,state,config); });
		auto port_free=state;port_free.port_control_values.clear();
		Require(iga::ParseMovingAcceptedCheckpointMetadata(legacy,port_free,config)==authority);
		auto invalid_controls=state;invalid_controls.port_control_values["inlet"]=std::numeric_limits<double>::infinity();
		Reject([&] { (void)iga::SerializeMovingAcceptedCheckpointMetadata(invalid_controls,config); });
		auto another=state;another.field={3,0,1,{9.}};
		Require(iga::SerializeMovingAcceptedCheckpointMetadata(another,config)==bytes);
		std::string field;
		iga::WriteOwnedCheckpointField(state.field,[&](const void* data,std::size_t count) { field.append(static_cast<const char*>(data),count); });
		Require(field.size()==iga::OwnedCheckpointFieldBytes(state.field));
		auto candidate=state.field;candidate.values.assign(2,0.);
		auto input=iga::ReadOwnedCheckpointField(candidate);
		for(char byte:field)input.Consume(&byte,1);
		input.Finish();Require(candidate.values==state.field.values);
		for(int mode=0;mode<4;++mode) {
			auto target=state.field;
			if(mode==0)target.row_begin=0;
			auto reader=iga::ReadOwnedCheckpointField(target);
			Reject([&] {
				if(mode==1) { reader.Consume(field.data(),field.size()-1);reader.Finish(); }
				else if(mode==2) { reader.Consume(field.data(),field.size());reader.Consume("x",1); }
				else if(mode==3) {
					auto invalid=field;
					for(int i=0;i<8;++i)invalid[24+i]=0;
					invalid[30]=static_cast<char>(0xf0);invalid[31]=static_cast<char>(0x7f);
					reader.Consume(invalid.data(),invalid.size());
				} else reader.Consume(field.data(),field.size());
			});
		}
		const auto encode=[](const iga::OwnedCheckpointVector& value) {
			std::string result;iga::WriteOwnedCheckpointField(value,[&](const void* data,std::size_t size) { result.append(static_cast<const char*>(data),size); });return result;
		};
		const std::vector<std::string> shards={encode({8,0,3,{1.,2.,3.}}),encode({8,3,6,{4.,5.,6.}}),encode({8,6,8,{7.,8.}}),encode({8,8,8,{}})};
		for(std::uint64_t ranks:{1u,2u,4u,9u})for(std::uint64_t rank=0;rank<ranks;++rank) {
			const auto begin=8*rank/ranks,end=8*(rank+1)/ranks;
			iga::OwnedCheckpointVector target{8,begin,end,std::vector<double>(end-begin,-1.)};
			iga::OwnedCheckpointRepartitionReader reader(target);
			for(auto it=shards.rbegin();it!=shards.rend();++it)reader.ReadShard(it->size(),[&](auto consume) {
				for(char byte:*it)consume(&byte,1);
			});
			reader.Finish();for(std::size_t i=0;i<target.values.size();++i)Require(target.values[i]==begin+i+1.);
		}
		for(int mode=0;mode<4;++mode) {
			iga::OwnedCheckpointVector target{8,0,1,{-1.}};iga::OwnedCheckpointRepartitionReader reader(target);
			const auto feed=[&](const std::string& data,std::size_t trim=0) {
				reader.ReadShard(data.size(),[&](auto consume) { consume(data.data(),data.size()-trim); });
			};
			Reject([&] {
				feed(shards[0]);
				if(mode==0)feed(shards[0]); // overlapping source ranges
				if(mode==1) { feed(shards[2]);reader.Finish(); } // global gap outside target
				if(mode==2)feed(shards[2],1); // truncated source outside target
				if(mode==3)feed(encode({8,2,4,{3.,4.}})); // partially overlapping range
			});
			Reject([&] { reader.Finish(); });
		}
		std::cout<<"moving checkpoint metadata configuration/geometry/epoch binding, conservation authority and owned stream roundtrip/rejection passed\n";
	} catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
	return 0;
}
