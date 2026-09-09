#include "BodyFittedCheckpointBundle.hpp"
#include "CollectiveCheckpointReceipts.hpp"
#define IGA_BODY_FITTED_CHECKPOINT_FIXTURE_ONLY
#include "test_body_fitted_accepted_checkpoint.cpp"
#undef IGA_BODY_FITTED_CHECKPOINT_FIXTURE_ONLY

namespace {
using namespace iga;
int codec_rejections = 0;
std::string Hash(const std::string& value) { return coupled_checkpoint_detail::Hash(value); }
template<class Work> void CodecReject(Work work)
{
	try { work(); } catch (const std::exception&) { ++codec_rejections; return; }
	throw std::runtime_error("expected codec rejection");
}
template<class Write> std::string Encoded(Write write)
{
	std::string bytes; write([&](const void* data, std::size_t size) { bytes.append(static_cast<const char*>(data), size); }); return bytes;
}
std::string FlowImage(const FlowAcceptedCheckpointState& flow)
{
	return SerializeBodyFittedFlowMetadata(flow)+Encoded([&](auto sink) { WriteOwnedCheckpointField(flow.field, sink); })
		+Encoded([&](auto sink) { WriteBodyFittedBoundaries(flow, sink); });
}
std::string TransportImage(const TransportAcceptedCheckpointState& transport)
{
	return SerializeBodyFittedTransportMetadata(transport,kDt)+Encoded([&](auto sink) { WriteOwnedCheckpointField(transport.field,sink); });
}
void CodecTests()
{
	FlowAcceptedCheckpointState state; state.configuration_identity_sha256 = std::string(64,'a');
	state.accepted_steps=3; state.accepted_time_s=0.03; state.macro_dt_s=0.01; state.total_linear_iterations=7;
	state.field={40,8,24,std::vector<double>(16,1.25)}; state.field.values[2]=-0.0;
	state.boundaries.velocity.assign(10,{1,2,3}); state.boundaries.pressure.assign(10,-2.5);
	state.pressure_tractions={{2,3.5}}; OutletModelState outlet; outlet.label=2; outlet.kind=FieldBoundaryKind::WindkesselRC;
	outlet.resistance=1; outlet.capacitance=2; outlet.capacitor_pressure=3; outlet.pressure=4; outlet.flow=-0.2; state.outlets={outlet};
	const CoupledCheckpointEpoch epoch{Hash("codec"),{}, {Hash("case"),Hash("config"),Hash("execution"),3},3,0.03,0.01};
	const auto metadata=SerializeBodyFittedFlowMetadata(state), field=Encoded([&](auto sink){WriteOwnedCheckpointField(state.field,sink);});
	const auto boundary=Encoded([&](auto sink){WriteBodyFittedBoundaries(state,sink);});
	auto target=[&] {auto value=state;value.accepted_steps=0;value.accepted_time_s=0;return value;};
	for(std::size_t size=0;size<metadata.size();++size)CodecReject([&]{auto value=target();ParseBodyFittedFlowMetadata(std::string_view(metadata).substr(0,size),value,epoch);});
	CodecReject([&]{auto value=target();ParseBodyFittedFlowMetadata(metadata+"x",value,epoch);});
	for(int fault=0;fault<5;++fault)CodecReject([&]{
		auto value=target();auto wrong=epoch;
		if(fault==0)++wrong.accepted_steps;
		if(fault==1)wrong.time_s+=1;
		if(fault==2)wrong.dt_s*=2;
		if(fault==3)value.configuration_identity_sha256[0]='b';
		if(fault==4)value.boundaries.velocity.pop_back();
		ParseBodyFittedFlowMetadata(metadata,value,wrong);
	});
	for(std::size_t size=0;size<field.size();++size)CodecReject([&]{auto value=target();auto reader=ReadOwnedCheckpointField(value.field);reader.Consume(field.data(),size);reader.Finish();});
	for(std::size_t size=0;size<boundary.size();++size)CodecReject([&]{auto value=target();auto reader=ReadBodyFittedBoundaries(value);reader.Consume(boundary.data(),size);reader.Finish();});
	for(int fault=0;fault<3;++fault) {
		auto bad=field;bad[8*fault]^=1;auto value=target();auto reader=ReadOwnedCheckpointField(value.field);
		CodecReject([&]{reader.Consume(bad.data(),bad.size());});CodecReject([&]{reader.Finish();});
	}
	{
		auto bad=field;const std::uint64_t infinity=0x7ff0000000000000ULL;
		for(unsigned i=0;i<8;++i)bad[24+i]=static_cast<char>(infinity>>(8*i));
		auto value=target();auto reader=ReadOwnedCheckpointField(value.field);
		CodecReject([&]{reader.Consume(bad.data(),bad.size());});CodecReject([&]{reader.Consume(field.data(),field.size());});
	}
	{
		auto value=target();ParseBodyFittedFlowMetadata(metadata,value,epoch);auto owned=ReadOwnedCheckpointField(value.field);auto bounds=ReadBodyFittedBoundaries(value);
		for(std::size_t i=0;i<field.size();++i)owned.Consume(field.data()+i,1);
		for(std::size_t i=0;i<boundary.size();++i)bounds.Consume(boundary.data()+i,1);
		owned.Finish();bounds.Finish();Require(FlowImage(value)==FlowImage(state),"fragmented flow codec differs");
		CodecReject([&]{owned.Consume("x",1);});
	}
	TransportAcceptedCheckpointState transport{std::string(64,'b'),3,{10,2,6,{2,3,4,5}}};
	const auto scalar=SerializeBodyFittedTransportMetadata(transport,0.01);
	for(std::size_t size=0;size<scalar.size();++size)CodecReject([&]{auto value=transport;value.accepted_steps=0;ParseBodyFittedTransportMetadata(std::string_view(scalar).substr(0,size),value,0.01,epoch);});
	CodecReject([&]{auto value=transport;value.accepted_steps=0;auto wrong=epoch;wrong.time_s+=1;ParseBodyFittedTransportMetadata(scalar,value,0.01,wrong);});
	CodecReject([&]{BodyFittedCheckpointCatalog({"bad/path",{0},true},3);});
	CodecReject([&]{BodyFittedCheckpointCatalog({"body",{0,0},true},3);});
	CodecReject([&]{BodyFittedCheckpointCatalog({"body",{3},true},3);});
	{
		checkpoint_stream::Writer writer([](const void*,std::size_t){throw std::runtime_error("writer failure");});
		CodecReject([&]{for(int i=0;i<8192;++i)writer.Word(0);});CodecReject([&]{writer.Word(0);});
	}
	for(std::size_t count:{std::size_t{0},std::size_t{3*1024*1024}}) {
		OwnedCheckpointVector source{count,0,count,std::vector<double>(count,1.25)},copy{count,0,count,std::vector<double>(count,0)};
		auto reader=ReadOwnedCheckpointField(copy);std::size_t maximum=0;
		WriteOwnedCheckpointField(source,[&](const void* data,std::size_t size){maximum=std::max(maximum,size);reader.Consume(data,size);});reader.Finish();
		Require(copy.values==source.values && maximum<=65536,"large or empty owned stream differs");
		std::cout<<"body_checkpoint_codec bytes="<<OwnedCheckpointFieldBytes(source)<<" maximum_chunk="<<maximum<<'\n';
	}
	std::cout<<"body_checkpoint_codec status=passed rejections="<<codec_rejections<<'\n';
}

struct Fixture {
	MPI_Comm comm;
	std::unique_ptr<Database> database;
	SimulationConfiguration configuration,scalar;
	std::unique_ptr<TransientFlowRuntime> flow;
	std::unique_ptr<TransientTransportRuntime> transport;
	BodyFittedCheckpointLayout layout;
	std::string identity;
	int rank=0,world_rank=0,mode=0;
	Fixture(const fs::path& root,MPI_Comm communicator,int selected):comm(communicator),mode(selected)
	{
		MPI_Comm_rank(comm,&rank);MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);int size=1;MPI_Comm_size(comm,&size);
		layout.domain_id="body";layout.transport=true;layout.world_ranks.resize(size);
		const auto global=static_cast<std::uint32_t>(world_rank);MPI_Allgather(&global,1,MPI_UINT32_T,layout.world_ranks.data(),1,MPI_UINT32_T,comm);
		LabeledHexMesh mesh;std::vector<std::array<double,3>> velocity;std::vector<OutletModelState> outlets;
		ResolvedBoundaryConditions boundaries;std::set<std::int32_t> wall;
		CollectiveLocalStage(comm,"bundle fixture input",[&]{
			database=std::make_unique<Database>((root/"fixture/group.ntiga").string());
			mesh=ReadLabeledHexMesh((root/"fixture/controlmesh.vtk").string(),64,1);
			velocity=ReadVelocity((root/"fixture/initial_velocityfield.txt").string(),64);
			configuration=ReadSimulationConfiguration((root/"fixture/simulation_config.json").string());configuration.time.steps=6;
			if(mode==0)configuration.equation_systems.front().time_integration="steady";
			if(mode==2){auto& outlet=configuration.boundaries[2].conditions[0];outlet.kind=FieldBoundaryKind::WindkesselRC;outlet.resistance=1e-3;outlet.capacitance=1.0;}
			scalar=ScalarConfiguration();const auto& definition=FirstNavierStokesSystem(configuration);outlets=InitializeOutletModels(configuration,definition);
			boundaries=ResolveFlowBoundaries(MaterializeOutletPressures(configuration,outlets),definition,mesh.labels,velocity);wall=WallTraceBasis(*database,mesh,0);
			Sha256 hash;const auto recipe=std::string("body-bundle-fixture-v1-mode-")+std::to_string(mode);hash.Append(recipe.data(),recipe.size());
			for(const auto* name:{"group.ntiga","controlmesh.vtk","initial_velocityfield.txt","simulation_config.json"}) {
				std::ifstream file(root/"fixture"/name,std::ios::binary);Require(bool(file),"missing identity input");
				const std::string bytes((std::istreambuf_iterator<char>(file)),{});hash.AppendLittleEndian64(bytes.size());hash.Append(bytes.data(),bytes.size());
			}
			identity=hash.Hex();
		});
		const auto& definition=FirstNavierStokesSystem(configuration);
		flow=std::make_unique<TransientFlowRuntime>(*database,comm,true,mode!=0,NavierStokesParameters{definition.density,definition.viscosity,mode? kDt:0.0},
			boundaries,mesh.labels,velocity,wall,outlets,std::set<std::string>{},identity,kDt);flow->InitializeState(configuration);
		CompiledLinearSystem system;CollectiveLocalStage(comm,"bundle scalar compile",[&]{system=CompileLinearSystem(scalar,"transport");});
		transport=std::make_unique<TransientTransportRuntime>(*database,comm,scalar,system,mesh.labels,std::map<std::uint64_t,VolumeQuadratureRule>{},std::set<std::string>{},Hash(identity+"transport-v1"));
	}
	void Advance(int step)
	{
		auto changed=configuration;changed.boundaries[1].conditions[0].scale=1+0.05*step;
		flow->Advance(changed,step,flow->AcceptedTime()+kDt,12,1e-8,1e-12,1e-6);
		transport->Advance(scalar,flow->RequiredNodes(),flow->GatherRequiredVelocity());
	}
	CoupledCheckpointEpoch Epoch(int steps)const
	{
		std::string mapping;for(auto rank:layout.world_ranks)mapping+=std::to_string(rank)+",";
		return {Hash(identity+"-epoch-"+std::to_string(steps)),steps==3?"":Hash(identity+"-epoch-3"),{Hash(identity+"case"),identity,Hash(mapping),3},
			static_cast<std::uint64_t>(steps),steps==3?0.03:0.04,kDt};
	}
	std::string Fingerprint()
	{
		const auto f=flow->CaptureCheckpointState();const auto t=transport->CaptureCheckpointState();
		const auto mass=transport->TotalMass();std::string result;
		CollectiveLocalStage(comm,"bundle reference fingerprint",[&]{checkpoint_metadata::Writer extra;extra.Reals(mass);result=Hash(FlowImage(f)+TransportImage(t)+extra.Bytes());});return result;
	}
	void Close(){transport->Close();flow->Close();}
};

std::string BroadcastManifest(MPI_Comm comm,const std::string& root_bytes)
{
	int rank=0;MPI_Comm_rank(comm,&rank);std::uint64_t size=rank==0?root_bytes.size():0;MPI_Bcast(&size,1,MPI_UINT64_T,0,comm);
	std::string bytes;CollectiveLocalStage(comm,"bundle manifest receive",[&]{Require(size<=coupled_checkpoint_detail::maximum_manifest_bytes,"manifest too large");bytes=rank==0?root_bytes:std::string(size,'\0');});
	MPI_Bcast(bytes.data(),static_cast<int>(size),MPI_CHAR,0,comm);return bytes;
}
CoupledCheckpointManifest LoadLatest(Fixture& fixture,const fs::path& root,bool reject_expected)
{
	std::string bytes;
	CollectiveLocalStage(fixture.comm,"bundle latest verified manifest",[&]{if(fixture.rank==0){
		const auto expected=fixture.Epoch(3);const auto discovery=FindLatestCoupledCheckpoint(root,expected.compatibility,BodyFittedCheckpointCatalog(fixture.layout,3));
		Require(discovery.latest && discovery.latest->epoch.id==expected.id,"last complete epoch differs");
		if(reject_expected)Require(!discovery.rejected.empty(),"interrupted epoch was not reported");
		std::cout<<"body_checkpoint_discovery mode="<<fixture.mode<<" rejected="<<discovery.rejected.size()<<'\n';bytes=SerializeCoupledCheckpointManifest(*discovery.latest);
	}});
	bytes=BroadcastManifest(fixture.comm,bytes);CoupledCheckpointManifest manifest;
	CollectiveLocalStage(fixture.comm,"bundle manifest parse",[&]{manifest=ParseCoupledCheckpointManifest(bytes);});return manifest;
}
void Restore(Fixture& fixture,const fs::path& root,const CoupledCheckpointManifest& manifest,bool negatives)
{
	auto flow=fixture.flow->CreateCheckpointRestoreCandidate();auto transport=fixture.transport->CreateCheckpointRestoreCandidate();
	if(negatives){
		Reject(fixture.comm,"bundle candidate load",[&]{CollectiveLocalStage(fixture.comm,"bundle candidate load",[&]{
			LoadBodyFittedCheckpointShards(fixture.rank==1?root/"missing":root,manifest,fixture.layout,fixture.rank,flow,&transport,kDt);
		});});
		Require(fixture.flow->AcceptedSteps()==0 && fixture.transport->Steps()==0,"failed local load published state");
		flow=fixture.flow->CreateCheckpointRestoreCandidate();transport=fixture.transport->CreateCheckpointRestoreCandidate();
	}
	CollectiveLocalStage(fixture.comm,"bundle candidate load",[&]{LoadBodyFittedCheckpointShards(root,manifest,fixture.layout,fixture.rank,flow,&transport,kDt);});
	fixture.flow->RestoreCheckpointState(flow);fixture.transport->RestoreCheckpointState(transport);
}
void Save(Fixture& fixture,const fs::path& root)
{
	for(int step=0;step<3;++step)fixture.Advance(step);
	const auto epoch=fixture.Epoch(3);const auto f=fixture.flow->CaptureCheckpointState();const auto t=fixture.transport->CaptureCheckpointState();
	CollectiveLocalStage(fixture.comm,"bundle epoch create",[&]{if(fixture.rank==0)CreateCoupledCheckpointEpoch(root,epoch);});
	std::vector<CoupledCheckpointShard> local;
	CollectiveLocalStage(fixture.comm,"bundle shard write",[&]{local=WriteBodyFittedCheckpointShards(root,epoch,fixture.layout,fixture.rank,f,&t);});
	// Exercise empty producers and pre-publication receipt agreement failures.
	const CoupledCheckpointShard probe{{"probe","probe-v1",fixture.layout.world_ranks[0]},0,Hash("probe")};
	std::vector<CoupledCheckpointShard> empty_probe;
	if(fixture.rank==0)empty_probe.push_back(probe);
	const auto gathered=GatherCheckpointReceipts(fixture.comm,epoch,empty_probe);
	CollectiveLocalStage(fixture.comm,"empty producer check",[&]{Require(gathered.size()==(fixture.rank==0?1u:0u),"empty producer receipt gather differs");});
	auto wrong_epoch=epoch;if(fixture.rank==2)wrong_epoch.id=Hash("wrong epoch");
	Reject(fixture.comm,"checkpoint receipt epoch agreement",[&]{GatherCheckpointReceipts(fixture.comm,wrong_epoch,{});});
	Reject(fixture.comm,"checkpoint receipt catalog validation",[&]{GatherCheckpointReceipts(fixture.comm,epoch,{probe});});
	const auto receipts=GatherCheckpointReceipts(fixture.comm,epoch,local);
	CollectiveLocalStage(fixture.comm,"bundle publish",[&]{if(fixture.rank==0)PublishCoupledCheckpoint(root,epoch,BodyFittedCheckpointCatalog(fixture.layout,3),receipts);});
	std::ofstream reference;CollectiveLocalStage(fixture.comm,"reference open",[&]{reference.open(root/("reference-"+std::to_string(fixture.world_rank)+".txt"));Require(bool(reference),"reference open failed");});
	for(int step=2;step<6;++step){if(step>2)fixture.Advance(step);const auto value=fixture.Fingerprint();CollectiveLocalStage(fixture.comm,"reference write",[&]{reference<<value<<'\n';Require(bool(reference),"reference write failed");});}
	reference.close();CollectiveLocalStage(fixture.comm,"reference close",[&]{Require(bool(reference),"reference close failed");});
}
void Resume(Fixture& fixture,const fs::path& root,bool interrupted)
{
	const auto manifest=LoadLatest(fixture,root,interrupted);Restore(fixture,root,manifest,true);
	std::ifstream reference;CollectiveLocalStage(fixture.comm,"reference read open",[&]{reference.open(root/("reference-"+std::to_string(fixture.world_rank)+".txt"));Require(bool(reference),"reference missing");});
	for(int step=2;step<6;++step){if(step>2)fixture.Advance(step);const auto actual=fixture.Fingerprint();CollectiveLocalStage(fixture.comm,"new-job exact comparison",[&]{std::string expected;Require(bool(reference>>expected) && actual==expected,"new MPI job differs from uninterrupted reference");});}
	if(fixture.rank==0)std::cout<<"body_checkpoint_restart mode="<<fixture.mode<<" steps=3-6 exact=1 status=passed\n";
}
void Crash(Fixture& fixture,const fs::path& root,bool in_stream)
{
	const auto manifest=LoadLatest(fixture,root,false);Restore(fixture,root,manifest,false);fixture.Advance(3);
	const auto epoch=fixture.Epoch(4);const auto f=fixture.flow->CaptureCheckpointState();const auto t=fixture.transport->CaptureCheckpointState();
	CollectiveLocalStage(fixture.comm,"crash epoch",[&]{if(fixture.rank==0)CreateCoupledCheckpointEpoch(root,epoch);});
	CollectiveLocalStage(fixture.comm,"interrupted shard write",[&]{
		if(fixture.rank==0 && in_stream){
			const CoupledCheckpointShardSpec spec{"body.flow.rank-"+std::to_string(fixture.world_rank),"owned-real-field-v1",static_cast<std::uint32_t>(fixture.world_rank)};
			WriteCoupledCheckpointShard(root,epoch,spec,OwnedCheckpointFieldBytes(f.field),[&](auto& output){const char partial[5]={1,2,3,4,5};output.Write(partial,sizeof(partial));std::cerr<<"checkpoint_crash stage=field-stream bytes=5 epoch="<<epoch.id<<std::endl;::_exit(86);});
		}else WriteBodyFittedCheckpointShards(root,epoch,fixture.layout,fixture.rank,f,&t);
	});
	if(fixture.rank==0){std::cerr<<"checkpoint_crash stage=before-manifest epoch="<<epoch.id<<std::endl;::_exit(87);}
	MPI_Barrier(fixture.comm);throw std::runtime_error("interrupted rank unexpectedly returned");
}
} // namespace

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);int rank=0,size=1;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
	try{
		Require(argc==3,"usage: body_fitted_checkpoint_bundle_test --codec|--save|--resume|--crash-stream|--crash-manifest|--resume-interrupted ROOT");
		const std::string action=argv[1];const fs::path root=argv[2];
		if(action=="--codec"){CodecTests();PetscFinalize();return 0;}
		Require(size==3,"bundle test requires three ranks");
		for(int mode=0;mode<3;++mode){
			if(action=="--crash-stream" && mode!=0)continue;
			if(action=="--crash-manifest" && mode!=1)continue;
			MPI_Comm comm;MPI_Comm_split(MPI_COMM_WORLD,0,mode==2?(rank+1)%3:rank,&comm);int local_rank=0;MPI_Comm_rank(comm,&local_rank);
			const auto directory=root/("mode-"+std::to_string(mode));
			if(action=="--save")CollectiveLocalStage(comm,"fixture creation",[&]{if(local_rank==0){Require(fs::create_directories(directory/"fixture"),"fixture exists");WriteThreeDCase(directory/"fixture");WriteDatabase(directory/"fixture/group.ntiga",3);}});
			Fixture fixture(directory,comm,mode);
			if(action=="--save")Save(fixture,directory);
			else if(action=="--resume" || action=="--resume-interrupted")Resume(fixture,directory,action=="--resume-interrupted" && mode<2);
			else if(action=="--crash-stream" || action=="--crash-manifest")Crash(fixture,directory,action=="--crash-stream");
			else throw std::runtime_error("unknown action");
			fixture.Close();MPI_Comm_free(&comm);
		}
		if(rank==0)std::cout<<"body_checkpoint_bundle action="<<action<<" status=passed\n";
	}catch(const std::exception& error){std::cerr<<"rank "<<rank<<": "<<error.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);return 1;}
	PetscFinalize();return 0;
}
