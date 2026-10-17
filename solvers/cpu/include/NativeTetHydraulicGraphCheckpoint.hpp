#ifndef IGA_NATIVE_TET_HYDRAULIC_GRAPH_CHECKPOINT_HPP
#define IGA_NATIVE_TET_HYDRAULIC_GRAPH_CHECKPOINT_HPP

#include "CollectiveCheckpointReceipts.hpp"
#include "CoupledCheckpointBundle.hpp"
#include "NativeTetAleFlowDomainAdapter.hpp"
#include "PressureFlowComponentExecutor.hpp"
#include "SimulationGraph.hpp"
#include "ZeroDFlowCheckpoint.hpp"

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <map>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <string>
#include <utility>
#include <vector>

namespace iga {

// One contiguous range of accepted steps that used the same 0D boundary
// values (identified by the hash of every accepted 0D model).
struct NativeTetHydraulicBoundaryHistoryEntry
{
	std::uint64_t first_step=0,last_step=0;
	std::string values_sha256;
};

// Opt-in policy for controlled restarts that may change declared 0D boundary
// values. The configuration identity then uses each 0D owner's structural
// identity (its model without the declared values) instead of the full model
// identity, and an extra shard carries the exact accepted 0D models plus the
// boundary-value history, so every 0D state is still validated against the
// model that produced it.
struct NativeTetHydraulicBoundaryPolicy
{
	std::map<std::string,std::string> zero_d_structural_identity;
};

struct NativeTetHydraulicGraphImage
{
	CoupledCheckpointEpoch epoch;
	NativeTetAleFlowCheckpoint native;
	std::map<std::string,ZeroDFlowCheckpointState> zero_d;
	std::map<std::string,ZeroDFlowModel> accepted_zero_d_models;
	std::vector<NativeTetHydraulicBoundaryHistoryEntry> boundary_history;
};

namespace native_tet_hydraulic_checkpoint_detail {

inline void WriteModel(checkpoint_metadata::Writer& output,const ZeroDFlowModel& model)
{
	ValidateZeroDFlowModel(model);
	output.Unsigned(static_cast<std::uint64_t>(model.role));
	if(model.role==ZeroDFlowRole::SourceReservoir)
		for(const double value:{model.source.capacitance_m3_pa,model.source.resistance_pa_s_m3,
			model.source.prescribed_flow_m3_s})output.Real(value);
	else for(const double value:{model.terminal.proximal_resistance_pa_s_m3,
		model.terminal.distal_resistance_pa_s_m3,model.terminal.capacitance_m3_pa,
		model.terminal.distal_pressure_pa})output.Real(value);
}

inline ZeroDFlowModel ReadModel(checkpoint_metadata::Reader& input)
{
	ZeroDFlowModel model;
	const auto role=input.Unsigned();
	checkpoint_metadata::Require(role<=1,"unknown 0D role in boundary values");
	model.role=static_cast<ZeroDFlowRole>(role);
	if(model.role==ZeroDFlowRole::SourceReservoir){
		model.source.capacitance_m3_pa=input.Real();model.source.resistance_pa_s_m3=input.Real();
		model.source.prescribed_flow_m3_s=input.Real();
	}else{
		model.terminal.proximal_resistance_pa_s_m3=input.Real();
		model.terminal.distal_resistance_pa_s_m3=input.Real();
		model.terminal.capacitance_m3_pa=input.Real();model.terminal.distal_pressure_pa=input.Real();
	}
	ValidateZeroDFlowModel(model);
	return model;
}

// Hash of every 0D model, keyed and ordered by domain id.
inline std::string BoundaryValuesSha256(const std::map<std::string,ZeroDFlowModel>& models)
{
	Sha256 hash;
	const auto append=[&](const std::string& value){
		hash.AppendLittleEndian64(value.size());hash.Append(value.data(),value.size());
	};
	append("native-tet-hydraulic-boundary-values/v1");
	hash.AppendLittleEndian64(models.size());
	for(const auto& item:models){append(item.first);append(BuildZeroDFlowModelIdentitySha256(item.second));}
	return hash.Hex();
}

inline std::string SerializeBoundaryValues(const std::map<std::string,ZeroDFlowModel>& models,
	const std::vector<NativeTetHydraulicBoundaryHistoryEntry>& history)
{
	checkpoint_metadata::Writer output;
	output.Text("IGA_NATIVE_TET_HYDRAULIC_BOUNDARY_VALUES");output.Unsigned(1);
	output.Unsigned(models.size());
	for(const auto& item:models){output.Text(item.first);WriteModel(output,item.second);}
	output.Unsigned(history.size());
	for(const auto& entry:history){
		output.Unsigned(entry.first_step);output.Unsigned(entry.last_step);output.Text(entry.values_sha256);
	}
	return output.Bytes();
}

inline void ParseBoundaryValues(std::string_view bytes,NativeTetHydraulicGraphImage& image)
{
	checkpoint_metadata::Reader input(bytes);
	checkpoint_metadata::Require(input.Text()=="IGA_NATIVE_TET_HYDRAULIC_BOUNDARY_VALUES"
		&&input.Unsigned()==1,"unsupported boundary-values schema");
	const auto count=input.Unsigned();
	checkpoint_metadata::Require(count<=checkpoint_metadata::maximum_map_entries,"too many 0D models");
	for(std::uint64_t index=0;index<count;++index){
		auto domain=input.Text();
		checkpoint_metadata::Require(image.accepted_zero_d_models.emplace(domain,ReadModel(input)).second,
			"duplicate 0D model in boundary values");
	}
	const auto entries=input.Unsigned();
	checkpoint_metadata::Require(entries>=1&&entries<=checkpoint_metadata::maximum_map_entries,
		"boundary-value history length is invalid");
	std::uint64_t next=1;
	for(std::uint64_t index=0;index<entries;++index){
		NativeTetHydraulicBoundaryHistoryEntry entry;
		entry.first_step=input.Unsigned();entry.last_step=input.Unsigned();entry.values_sha256=input.Text();
		checkpoint_metadata::Require(entry.first_step==next&&entry.last_step>=entry.first_step
			&&entry.values_sha256.size()==64,"boundary-value history is not contiguous");
		next=entry.last_step+1;
		image.boundary_history.push_back(std::move(entry));
	}
	input.Finish();
}

} // namespace native_tet_hydraulic_checkpoint_detail

// A file-backed accepted publication for a hydraulic native tetra graph.
// The caller supplies case and execution identities; the configuration identity
// is derived from every owner model. Use a fresh, unpublished graph for restore.
class NativeTetHydraulicGraphCheckpoint
{
public:
	NativeTetHydraulicGraphCheckpoint(MPI_Comm comm,const SimulationGraph& graph,
		const PressureFlowExecutionControls& controls,
		NativeTetAleFlowDomainAdapter& native,
		std::map<std::string,ZeroDFlowDomainRuntime*> zero_d,
		std::optional<NativeTetHydraulicBoundaryPolicy> policy=std::nullopt)
		:comm_(comm),native_(native),zero_d_(std::move(zero_d)),policy_(std::move(policy))
	{
		MPI_Comm_rank(comm_,&rank_);MPI_Comm_size(comm_,&ranks_);
		ValidatePressureFlowExecutionControls(controls);
		if(zero_d_.empty())
			throw std::invalid_argument("native tetra hydraulic checkpoint needs 0D owners");
		if(graph.Domains().size()!=zero_d_.size()+1||graph.Edges().empty()
			||!graph.Species().empty()||!graph.FsiEdges().empty())
			throw std::invalid_argument("native tetra hydraulic graph topology differs");
		const auto prefix=[](const std::string& id){
			return coupled_checkpoint_detail::Hash("native-tet-hydraulic/v1:"+id);
		};
		catalog_.push_back({prefix(native_.DomainId())+".native",
			"native-tet-ale-flow-v1",0});
		for(const auto& item:zero_d_){
			if(item.first.empty()||item.first==native_.DomainId()||!item.second
				||item.second->DomainId()!=item.first)
				throw std::invalid_argument("native tetra hydraulic checkpoint owner differs");
			catalog_.push_back({prefix(item.first)+".zero-d","zero-d-accepted-v1",0});
			if(policy_&&!policy_->zero_d_structural_identity.count(item.first))
				throw std::invalid_argument("native tetra hydraulic boundary policy misses a 0D owner");
		}
		if(policy_){
			if(policy_->zero_d_structural_identity.size()!=zero_d_.size())
				throw std::invalid_argument("native tetra hydraulic boundary policy has extra owners");
			catalog_.push_back({BoundaryShardId(),"native-tet-hydraulic-boundary-values-v1",0});
		}
		std::sort(catalog_.begin(),catalog_.end(),[](const auto& a,const auto& b){
			return a.id<b.id;
		});
		coupled_checkpoint_detail::ValidateCatalog(catalog_,ranks_);
		Sha256 hash;
		const auto append=[&](const std::string& value){
			hash.AppendLittleEndian64(value.size());hash.Append(value.data(),value.size());
		};
		append("native-tet-hydraulic-models/v1");
		append(native_.DomainId());append(native_.ModelIdentitySha256());
		if(policy_)append("native-tet-hydraulic-boundary-policy/v1");
		for(const auto& item:zero_d_){
			append(item.first);append(policy_?policy_->zero_d_structural_identity.at(item.first)
				:item.second->ModelIdentitySha256());
		}
		hash.AppendLittleEndian64(graph.Domains().size());
		for(const auto& item:graph.Domains()){
			const auto& domain=item.second;
			if((item.first==native_.DomainId()
				&&domain.kind!=DomainKind::ThreeDBodyFittedFlow)
				|| (item.first!=native_.DomainId()
					&&(!zero_d_.count(item.first)||domain.kind!=DomainKind::ZeroDFlow))
				||!domain.species_bindings.empty()
				||!domain.surface_interfaces.empty()||!domain.surface_layouts.empty())
				throw std::invalid_argument("native tetra hydraulic graph owners differ");
			append(item.first);
			hash.AppendLittleEndian32(static_cast<std::uint32_t>(domain.kind));
			hash.AppendLittleEndian64(domain.ports.size());
			for(const auto& port:domain.ports){
				append(port.id);append(port.subsystem_id);
				append(port.locator_kind);append(port.locator);
				hash.AppendLittleEndian32(static_cast<std::uint32_t>(
					port.orientation.native_to_outward_sign));
				hash.AppendLittleEndian64(port.provides.size());
				for(const auto quantity:port.provides)
					hash.AppendLittleEndian32(static_cast<std::uint32_t>(quantity));
				hash.AppendLittleEndian64(port.requires.size());
				for(const auto quantity:port.requires)
					hash.AppendLittleEndian32(static_cast<std::uint32_t>(quantity));
				hash.AppendLittleEndian64(port.species.size());
				for(const auto& species:port.species)append(species);
			}
		}
		hash.AppendLittleEndian64(graph.Edges().size());
		for(const auto& edge:graph.Edges()){
			if(edge.law!=CouplingLaw::PressureFlow||!edge.species.empty())
				throw std::invalid_argument("native tetra hydraulic graph edge differs");
			append(edge.id);append(edge.first.domain_id);append(edge.first.port_id);
			append(edge.second.domain_id);append(edge.second.port_id);
			hash.AppendLittleEndian32(static_cast<std::uint32_t>(edge.law));
		}
		hash.AppendLittleEndian32(static_cast<std::uint32_t>(controls.method));
		hash.AppendLittleEndian32(static_cast<std::uint32_t>(controls.maximum_iterations));
		for(const double value:{controls.pressure_relative_tolerance,
			controls.pressure_reference_pa,controls.flow_relative_tolerance,
			controls.relaxation_factor,controls.minimum_relaxation,
			controls.maximum_relaxation})hash.AppendNormalizedDouble(value);
		configuration_sha256_=hash.Hex();
	}

	CoupledCheckpointCompatibility Compatibility(const std::string& case_sha256,
		const std::string& execution_sha256)const
	{
		CoupledCheckpointCompatibility value{case_sha256,configuration_sha256_,
			execution_sha256,static_cast<std::uint32_t>(ranks_)};
		coupled_checkpoint_detail::Validate(value);return value;
	}

	const std::vector<CoupledCheckpointShardSpec>& Catalog()const noexcept
	{return catalog_;}

	void Save(const std::filesystem::path& root,const CoupledCheckpointEpoch& epoch)
	{
		std::vector<std::pair<CoupledCheckpointShardSpec,std::string>> images;
		CollectiveLocalStage(comm_,"native hydraulic checkpoint capture",[&]{
			coupled_checkpoint_detail::Validate(epoch);
			coupled_checkpoint_detail::Require(
				epoch.compatibility.configuration_sha256==configuration_sha256_
				&&epoch.compatibility.ranks==static_cast<std::uint32_t>(ranks_)
				&&native_.AcceptedSteps()==epoch.accepted_steps
				&&native_.CommittedTime()==epoch.time_s,
				"native hydraulic checkpoint epoch or model differs");
			auto native=SerializeNativeTetAleFlowCheckpoint(native_.CaptureCheckpoint());
			images.emplace_back(Spec(native_.DomainId(),".native"),std::move(native));
			for(const auto& item:zero_d_){
				auto value=item.second->CaptureCheckpointState();
				checkpoint_metadata::RequireAcceptedStep(value.accepted_step,
					epoch.accepted_steps,epoch.time_s,epoch.dt_s);
				images.emplace_back(Spec(item.first,".zero-d"),
					SerializeZeroDFlowCheckpoint(value));
			}
			if(policy_){
				const auto models=AcceptedModels();
				const auto values=native_tet_hydraulic_checkpoint_detail::BoundaryValuesSha256(models);
				auto history=history_;
				if(!history.empty()&&history.back().values_sha256==values
					&&history.back().last_step+1==epoch.accepted_steps)
					history.back().last_step=epoch.accepted_steps;
				else{
					coupled_checkpoint_detail::Require(history.empty()
						?epoch.accepted_steps==1:history.back().last_step+1==epoch.accepted_steps,
						"native hydraulic boundary-value history is not contiguous");
					history.push_back({epoch.accepted_steps,epoch.accepted_steps,values});
				}
				images.emplace_back(BoundarySpec(),
					native_tet_hydraulic_checkpoint_detail::SerializeBoundaryValues(models,history));
				pending_history_=std::move(history);
			}
			std::sort(images.begin(),images.end(),[](const auto& a,const auto& b){
				return a.first.id<b.first.id;
			});
		});
		RequireCollectiveSameText(comm_,"native hydraulic epoch agreement",
			CheckpointEpochAgreementBytes(epoch));
		for(const auto& item:images)
			RequireCollectiveSameText(comm_,"native hydraulic replicated shard",
				coupled_checkpoint_detail::Hash(item.second));
		CollectiveLocalStage(comm_,"native hydraulic epoch create",[&]{
			if(rank_==0)CreateCoupledCheckpointEpoch(root,epoch);
		});
		std::vector<CoupledCheckpointShard> local;
		CollectiveLocalStage(comm_,"native hydraulic shard write",[&]{
			if(rank_!=0)return;
			for(const auto& item:images)
				local.push_back(WriteCoupledCheckpointShard(root,epoch,item.first,
					item.second.size(),[&](auto& output){
						output.Write(item.second.data(),item.second.size());
					}));
		});
		const auto receipts=GatherCheckpointReceipts(comm_,epoch,local);
		CollectiveLocalStage(comm_,"native hydraulic manifest publish",[&]{
			if(rank_==0)PublishCoupledCheckpoint(root,epoch,catalog_,receipts);
		});
		history_.swap(pending_history_);
		pending_history_.clear();
	}

	// Boundary-value history up to the latest published or restored epoch.
	const std::vector<NativeTetHydraulicBoundaryHistoryEntry>& BoundaryHistory()const noexcept
	{return history_;}

	// Accepted models of every 0D owner, keyed by domain id.
	std::map<std::string,ZeroDFlowModel> AcceptedModels()const
	{
		std::map<std::string,ZeroDFlowModel> models;
		for(const auto& item:zero_d_)models.emplace(item.first,item.second->AcceptedModel());
		return models;
	}

	NativeTetHydraulicGraphImage LoadLatest(const std::filesystem::path& root,
		const CoupledCheckpointCompatibility& expected)const
	{
		std::string manifest_bytes;
		CollectiveLocalStage(comm_,"native hydraulic restore compatibility",[&]{
			coupled_checkpoint_detail::Validate(expected);
			coupled_checkpoint_detail::Require(
				expected.configuration_sha256==configuration_sha256_
				&&expected.ranks==static_cast<std::uint32_t>(ranks_),
				"native hydraulic restore model or rank count differs");
		});
		RequireCollectiveSameText(comm_,"native hydraulic restore identity",
			expected.case_sha256+expected.configuration_sha256
				+expected.execution_sha256);
		CollectiveLocalStage(comm_,"native hydraulic checkpoint discovery",[&]{
			if(rank_!=0)return;
			const auto found=FindLatestCoupledCheckpoint(root,expected,catalog_);
			coupled_checkpoint_detail::Require(bool(found.latest),
				"no compatible complete native hydraulic checkpoint");
			manifest_bytes=SerializeCoupledCheckpointManifest(*found.latest);
		});
		std::uint64_t length=manifest_bytes.size();
		MPI_Bcast(&length,1,MPI_UINT64_T,0,comm_);
		CollectiveLocalStage(comm_,"native hydraulic manifest allocation",[&]{
			coupled_checkpoint_detail::Require(
				length<=coupled_checkpoint_detail::maximum_manifest_bytes,
				"native hydraulic manifest exceeds limit");
			manifest_bytes.resize(static_cast<std::size_t>(length));
		});
		MPI_Bcast(manifest_bytes.data(),static_cast<int>(length),MPI_CHAR,0,comm_);
		NativeTetHydraulicGraphImage image;
		CollectiveLocalStage(comm_,"native hydraulic candidate read",[&]{
			const auto manifest=ParseCoupledCheckpointManifest(manifest_bytes);
			coupled_checkpoint_detail::RequireCompatible(manifest.epoch.compatibility,
				expected);
			image.epoch=manifest.epoch;
			const auto read=[&](const CoupledCheckpointShardSpec& spec){
				const auto found=std::find_if(manifest.shards.begin(),manifest.shards.end(),
					[&](const auto& shard){return shard.spec.id==spec.id;});
				coupled_checkpoint_detail::Require(found!=manifest.shards.end()
					&&coupled_checkpoint_detail::Same(found->spec,spec)
					&&found->payload_bytes<=checkpoint_metadata::maximum_bytes,
					"native hydraulic shard catalog or size differs");
				std::string bytes;
				ReadCoupledCheckpointShard(root,manifest,
					static_cast<std::size_t>(found-manifest.shards.begin()),
					[&](const void* data,std::size_t count){
						checkpoint_metadata::Require(count<=checkpoint_metadata::maximum_bytes
							-bytes.size(),"native hydraulic shard exceeds limit");
						bytes.append(static_cast<const char*>(data),count);
					});
				return bytes;
			};
			image.native=ParseNativeTetAleFlowCheckpoint(
				read(Spec(native_.DomainId(),".native")));
			coupled_checkpoint_detail::Require(
				image.native.model_identity_sha256==native_.ModelIdentitySha256()
				&&image.native.accepted_steps==image.epoch.accepted_steps
				&&image.native.time_s==image.epoch.time_s,
				"native hydraulic state, model, or epoch differs");
			for(const auto& item:zero_d_)
				image.zero_d.emplace(item.first,ParseZeroDFlowCheckpoint(
					read(Spec(item.first,".zero-d")),image.epoch));
			if(policy_){
				native_tet_hydraulic_checkpoint_detail::ParseBoundaryValues(read(BoundarySpec()),image);
				coupled_checkpoint_detail::Require(image.accepted_zero_d_models.size()==zero_d_.size()
					&&image.boundary_history.back().last_step==image.epoch.accepted_steps
					&&image.boundary_history.back().values_sha256==native_tet_hydraulic_checkpoint_detail::
						BoundaryValuesSha256(image.accepted_zero_d_models),
					"native hydraulic boundary values differ from the epoch");
				for(const auto& item:zero_d_){
					const auto found=image.accepted_zero_d_models.find(item.first);
					coupled_checkpoint_detail::Require(found!=image.accepted_zero_d_models.end()
						&&BuildZeroDFlowModelIdentitySha256(found->second)
							==image.zero_d.at(item.first).model_identity_sha256,
						"native hydraulic 0D state and accepted boundary values differ");
				}
			}
		});
		return image;
	}

	// The graph owner must discard every candidate runtime if this fails.
	void RestoreFresh(const NativeTetHydraulicGraphImage& image)
	{
		CollectiveLocalStage(comm_,"native hydraulic fresh restore",[&]{
			coupled_checkpoint_detail::Require(
				image.epoch.compatibility.configuration_sha256==configuration_sha256_
				&&image.epoch.compatibility.ranks==static_cast<std::uint32_t>(ranks_)
				&&image.native.accepted_steps==image.epoch.accepted_steps
				&&image.native.time_s==image.epoch.time_s
				&&image.zero_d.size()==zero_d_.size(),
				"native hydraulic 0D catalog differs");
			native_.RestoreCheckpoint(image.native);
			for(const auto& item:zero_d_)
				if(policy_)item.second->RestoreCheckpointStateUnderModel(image.zero_d.at(item.first),
					image.accepted_zero_d_models.at(item.first));
				else item.second->RestoreCheckpointState(image.zero_d.at(item.first));
			history_=image.boundary_history;
		});
	}

private:
	static std::string BoundaryShardId()
	{
		return coupled_checkpoint_detail::Hash("native-tet-hydraulic/v1:boundary-values")
			+".boundary-values";
	}
	const CoupledCheckpointShardSpec& BoundarySpec()const
	{
		const auto id=BoundaryShardId();
		const auto found=std::find_if(catalog_.begin(),catalog_.end(),
			[&](const auto& item){return item.id==id;});
		coupled_checkpoint_detail::Require(found!=catalog_.end(),"unknown native hydraulic shard");
		return *found;
	}
	const CoupledCheckpointShardSpec& Spec(const std::string& domain,
		const std::string& suffix)const
	{
		const auto id=coupled_checkpoint_detail::Hash(
			"native-tet-hydraulic/v1:"+domain)+suffix;
		const auto found=std::find_if(catalog_.begin(),catalog_.end(),
			[&](const auto& item){return item.id==id;});
		coupled_checkpoint_detail::Require(found!=catalog_.end(),"unknown native hydraulic shard");
		return *found;
	}
	MPI_Comm comm_;
	NativeTetAleFlowDomainAdapter& native_;
	std::map<std::string,ZeroDFlowDomainRuntime*> zero_d_;
	std::optional<NativeTetHydraulicBoundaryPolicy> policy_;
	std::vector<NativeTetHydraulicBoundaryHistoryEntry> history_,pending_history_;
	int rank_=0,ranks_=1;
	std::string configuration_sha256_;
	std::vector<CoupledCheckpointShardSpec> catalog_;
};

} // namespace iga

#endif
