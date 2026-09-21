#include "CoupledCheckpointBundle.hpp"
#include "ZeroDSourceReservoirSpeciesDomainRuntime.hpp"
#include "ZeroDTerminalRcrSpeciesDomainRuntime.hpp"

#include <cassert>
#include <filesystem>
#include <functional>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {

using namespace iga;
namespace fs=std::filesystem;

std::string Digest(std::string_view value)
{
	Sha256 hash;hash.Append(value.data(),value.size());return hash.Hex();
}

ZeroDSourceReservoirSpeciesDomainRuntime Source()
{
	ZeroDFlowModel hydraulic;
	hydraulic.role=ZeroDFlowRole::SourceReservoir;
	hydraulic.source={0.01,0.1,0.1};
	ZeroDSpeciesReservoirModel species;
	species.port_ids={"graph","pump"};species.species_ids={"tracer"};
	return {"source",hydraulic,{0.},species,{1.,{{"tracer",2.}}},{{"tracer",2.}}};
}

ZeroDTerminalRcrSpeciesDomainRuntime Terminal()
{
	ZeroDFlowModel hydraulic;
	hydraulic.role=ZeroDFlowRole::TerminalRcr;
	hydraulic.terminal={1.,10.,1.,0.};
	ZeroDSpeciesReservoirModel species;
	species.port_ids={"graph","distal"};species.species_ids={"tracer"};
	return {"terminal",hydraulic,{0.},species,{1.,{{"tracer",2.}}}};
}

void Advance(ZeroDSourceReservoirSpeciesDomainRuntime& source,
	ZeroDTerminalRcrSpeciesDomainRuntime& terminal,int index)
{
	const DomainStepContext step{index,index*0.05,0.05};
	source.BeginStep(step);
	PortBoundaryData pressure;pressure.time_s=step.EndTime();pressure.mean_pressure_pa=0.;
	source.SetPortInput("port",pressure);
	source.SolveTrial();source.PrepareCommitStep();source.FinalizeCommitStep();
	terminal.BeginStep(step);
	PortBoundaryData flow;flow.time_s=step.EndTime();flow.outward_flow_m3_s=-0.1;
	terminal.SetPortInput("port",flow);
	terminal.SolveHydraulicTrial();
	terminal.SetTransportConcentration("port",step.EndTime(),{{"tracer",2.}});
	terminal.SolveTransportTrial();terminal.PrepareCommitStep();terminal.FinalizeCommitStep();
}

CoupledCheckpointEpoch Epoch()
{
	return {Digest("zero-d-species-step-0"),{},
		{Digest("zero-d-species-case"),Digest("zero-d-species-config"),
			Digest("zero-d-species-execution"),1},1,0.05,0.05};
}

std::vector<CoupledCheckpointShardSpec> Catalog()
{
	return {{"source","zero-d-flow-species-v1",0},
		{"terminal","zero-d-flow-species-v1",0}};
}

void RequireRejected(const std::function<void()>& work)
{
	bool rejected=false;
	try{work();}catch(const std::exception&){rejected=true;}
	if(!rejected)throw std::runtime_error("invalid 0D species checkpoint was accepted");
}

void Save(const fs::path& root)
{
	if(!fs::create_directory(root))throw std::runtime_error("checkpoint test directory exists");
	auto source=Source();auto terminal=Terminal();Advance(source,terminal,0);
	const auto epoch=Epoch();const auto catalog=Catalog();
	CreateCoupledCheckpointEpoch(root,epoch);
	std::vector<CoupledCheckpointShard> receipts;
	for(std::size_t index=0;index<catalog.size();++index){
		const auto state=index==0?source.CaptureCheckpointState()
			:terminal.CaptureCheckpointState();
		const auto bytes=SerializeZeroDFlowSpeciesCheckpoint(state);
		receipts.push_back(WriteCoupledCheckpointShard(root,epoch,catalog[index],
			bytes.size(),[&](auto& output){output.Write(bytes.data(),bytes.size());}));
	}
	PublishCoupledCheckpoint(root,epoch,catalog,receipts);
}

void Resume(const fs::path& root)
{
	const auto epoch=Epoch();const auto catalog=Catalog();
	const auto manifest=LoadCoupledCheckpoint(root,epoch.id,epoch.compatibility,catalog);
	std::vector<ZeroDFlowSpeciesCheckpointState> restored;
	for(std::size_t index=0;index<catalog.size();++index){
		std::string bytes;
		ReadCoupledCheckpointShard(root,manifest,index,[&](const void* data,std::size_t count){
			bytes.append(static_cast<const char*>(data),count);
		});
		const auto value=ParseZeroDFlowSpeciesCheckpoint(bytes,manifest.epoch);
		if(SerializeZeroDFlowSpeciesCheckpoint(value)!=bytes)
			throw std::runtime_error("0D species checkpoint codec changed accepted bytes");
		RequireRejected([&]{ParseZeroDFlowSpeciesCheckpoint(
			std::string_view(bytes).substr(0,bytes.size()-1));});
		auto corrupted=bytes;corrupted[corrupted.size()/2]^=1;
		RequireRejected([&]{ParseZeroDFlowSpeciesCheckpoint(corrupted);});
		auto wrong_epoch=manifest.epoch;wrong_epoch.accepted_steps++;
		RequireRejected([&]{ParseZeroDFlowSpeciesCheckpoint(bytes,wrong_epoch);});
		restored.push_back(value);
	}
	auto source=Source();auto terminal=Terminal();
	RequireRejected([&]{source.RestoreCheckpointState(restored[1]);});
	RequireRejected([&]{terminal.RestoreCheckpointState(restored[0]);});
	auto wrong=restored[0];wrong.species.amount_mol["unknown"]=1.;
	RequireRejected([&]{source.RestoreCheckpointState(wrong);});
	wrong=restored[0];wrong.model_identity_sha256=Digest("different-model");
	RequireRejected([&]{source.RestoreCheckpointState(wrong);});
	if(source.CommittedStepIndex()!=-1||terminal.CommittedStepIndex()!=-1)
		throw std::runtime_error("rejected restore mutated 0D state");
	source.RestoreCheckpointState(restored[0]);
	terminal.RestoreCheckpointState(restored[1]);
	RequireRejected([&]{source.RestoreCheckpointState(restored[0]);});
	Advance(source,terminal,1);
	auto reference_source=Source();auto reference_terminal=Terminal();
	Advance(reference_source,reference_terminal,0);
	Advance(reference_source,reference_terminal,1);
	if(SerializeZeroDFlowSpeciesCheckpoint(source.CaptureCheckpointState())
		!=SerializeZeroDFlowSpeciesCheckpoint(reference_source.CaptureCheckpointState())
		||SerializeZeroDFlowSpeciesCheckpoint(terminal.CaptureCheckpointState())
		!=SerializeZeroDFlowSpeciesCheckpoint(reference_terminal.CaptureCheckpointState()))
		throw std::runtime_error("cross-process 0D species continuation differs");
}

} // namespace

int main(int argc,char** argv)
{
	if(argc!=3)throw std::runtime_error("usage: checkpoint-test save|resume ROOT");
	const std::string mode=argv[1];
	if(mode=="save")Save(argv[2]);
	else if(mode=="resume")Resume(argv[2]);
	else throw std::runtime_error("unknown checkpoint test mode");
	return 0;
}
