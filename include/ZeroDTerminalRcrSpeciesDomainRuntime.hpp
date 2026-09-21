#ifndef IGA_ZERO_D_TERMINAL_RCR_SPECIES_DOMAIN_RUNTIME_HPP
#define IGA_ZERO_D_TERMINAL_RCR_SPECIES_DOMAIN_RUNTIME_HPP

#include "ZeroDTerminalRcrSpecies.hpp"
#include "ZeroDFlowSpeciesCheckpoint.hpp"

#include <cmath>
#include <exception>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

inline CouplingPort MakeZeroDTerminalRcrSpeciesPort(const std::string& domain_id,
	const std::set<std::string>& species_ids)
{
	if(species_ids.empty())throw std::invalid_argument("0D terminal species set is empty");
	auto port=MakeZeroDFlowPort(domain_id,ZeroDFlowRole::TerminalRcr);
	port.provides.insert(PortQuantity::SpeciesConcentration);
	port.provides.insert(PortQuantity::SpeciesFlux);
	port.requires.insert(PortQuantity::SpeciesConcentration);
	port.requires.insert(PortQuantity::SpeciesFlux);
	port.species=species_ids;
	return port;
}

class ZeroDTerminalRcrSpeciesDomainRuntime final : public CoupledDomainRuntime,
	public StagedFlowTransportDomainRuntime
{
public:
	ZeroDTerminalRcrSpeciesDomainRuntime(std::string domain_id,
		ZeroDFlowModel hydraulic_model,ZeroDFlowState initial_hydraulic,
		ZeroDSpeciesReservoirModel species_model,
		ZeroDSpeciesReservoirState initial_species,
		std::map<std::string,double> distal_donor_concentration_mol_m3={},
		double initial_time_s=0.0)
		:domain_id_(std::move(domain_id)),hydraulic_model_(std::move(hydraulic_model)),
		 species_model_(std::move(species_model)),hydraulic_committed_(initial_hydraulic),
		 species_committed_(std::move(initial_species)),
		 distal_donor_(std::move(distal_donor_concentration_mol_m3)),
		 committed_time_s_(initial_time_s)
	{
		if(domain_id_.empty()||!std::isfinite(initial_time_s))
			throw std::invalid_argument("0D terminal species domain id or time is invalid");
		if(hydraulic_model_.role!=ZeroDFlowRole::TerminalRcr
			||species_model_.port_ids!=std::set<std::string>{"graph","distal"})
			throw std::invalid_argument("0D terminal species requires RCR graph/distal ports");
		ValidateZeroDFlowModel(hydraulic_model_);
		ValidateZeroDFlowState(hydraulic_committed_);
		ValidateZeroDSpeciesReservoirState(species_model_,species_committed_);
		if(!distal_donor_.empty()){
			if(distal_donor_.size()!=species_model_.species_ids.size())
				throw std::invalid_argument("0D terminal distal donor species set is incomplete");
			for(const auto& species:species_model_.species_ids){
				const auto value=distal_donor_.find(species);
				if(value==distal_donor_.end()||!(value->second>=0.0)
					||!std::isfinite(value->second))
					throw std::invalid_argument("0D terminal distal donor concentration is invalid");
			}
		}
		ports_.push_back(MakeZeroDTerminalRcrSpeciesPort(domain_id_,
			species_model_.species_ids));
	}

	const std::string& DomainId()const noexcept override{return domain_id_;}
	DomainKind Kind()const noexcept override{return DomainKind::ZeroDFlow;}
	const std::vector<CouplingPort>& Ports()const noexcept override{return ports_;}
	const ZeroDFlowState& CommittedHydraulicState()const noexcept{return hydraulic_committed_;}
	const ZeroDSpeciesReservoirState& CommittedSpeciesState()const noexcept
	{return species_committed_;}
	int CommittedStepIndex()const noexcept{return committed_step_index_;}
	std::string CheckpointModelIdentitySha256()const{return ModelIdentitySha256();}
	ZeroDFlowSpeciesCheckpointState CaptureCheckpointState()const
	{
		RequirePhase(Phase::Committed,"capture checkpoint");
		if(committed_step_index_<0)
			throw std::runtime_error("0D terminal checkpoint requires an accepted step");
		return {domain_id_,ModelIdentitySha256(),committed_step_,
			hydraulic_committed_,species_committed_};
	}
	void RestoreCheckpointState(ZeroDFlowSpeciesCheckpointState value)
	{
		RequirePhase(Phase::Committed,"restore checkpoint");
		if(committed_step_index_>=0||value.domain_id!=domain_id_
			||value.model_identity_sha256!=ModelIdentitySha256())
			throw std::runtime_error("0D terminal checkpoint identity or fresh-state mismatch");
		ValidateZeroDFlowSpeciesCheckpointState(value);
		ValidateZeroDSpeciesReservoirState(species_model_,value.species);
		committed_step_=value.accepted_step;
		committed_time_s_=committed_step_.EndTime();
		committed_step_index_=committed_step_.step_index;
		hydraulic_committed_=value.hydraulic;
		species_committed_.volume_m3=value.species.volume_m3;
		species_committed_.amount_mol.swap(value.species.amount_mol);
	}

	void BeginStep(const DomainStepContext& step)override
	{
		RequirePhase(Phase::Committed,"begin step");step.Validate();
		if(committed_step_index_==std::numeric_limits<int>::max()
			||step.start_time_s!=committed_time_s_
			||step.step_index!=committed_step_index_+1)
			throw std::runtime_error("0D terminal species step differs from committed clock");
		step_=step;ClearTrial();phase_=Phase::TrialReady;
	}

	void SetPortInput(const std::string& port_id,const PortBoundaryData& input)override
	{
		RequirePhase(Phase::TrialReady,"set hydraulic input");RequirePort(port_id);
		ValidatePortBoundaryData(input);
		if(input.time_s!=step_.EndTime()||!input.outward_flow_m3_s
			||input.mean_pressure_pa||input.mean_normal_traction_pa||input.total_pressure_pa
			||!input.concentration.empty()||!input.outward_species_flux.empty()
			||flow_input_)
			throw std::invalid_argument("0D terminal species requires one flow-only input");
		flow_input_=*input.outward_flow_m3_s;
	}

	void SolveHydraulicTrial()override
	{
		RequirePhase(Phase::TrialReady,"solve hydraulic trial");
		if(!flow_input_)throw std::runtime_error("0D terminal species is missing flow input");
		hydraulic_trial_=EvaluateZeroDFlowTrial(hydraulic_model_,hydraulic_committed_,
			*flow_input_,step_.dt_s,step_.start_time_s);
		phase_=Phase::HydraulicSolved;
	}
	void SolveTrial()override
	{
		SolveHydraulicTrial();SolveTransportTrial();
	}
	PortState GetHydraulicPortState(const std::string& port_id)const override
	{
		RequirePort(port_id);
		if(phase_!=Phase::HydraulicSolved&&phase_!=Phase::TransportSolved
			&&phase_!=Phase::Prepared)
			throw std::runtime_error("0D terminal hydraulic trial state is unavailable");
		return hydraulic_trial_->port;
	}
	PortState GetPortState(const std::string& port_id)const override
	{
		return phase_==Phase::TransportSolved||phase_==Phase::Prepared
			?GetTransportPortState(port_id):GetHydraulicPortState(port_id);
	}
	void RollbackHydraulicTrial()override
	{
		RequirePhase(Phase::HydraulicSolved,"rollback hydraulic trial");
		ClearTrial();phase_=Phase::TrialReady;
	}
	void RollbackTrial()override
	{
		if(phase_==Phase::TransportSolved)RollbackTransportTrial();
		RollbackHydraulicTrial();
	}

	void SetTransportConcentration(const std::string& port_id,double time_s,
		const std::map<std::string,double>& concentration)override
	{
		RequirePhase(Phase::HydraulicSolved,"set donor concentration");RequirePort(port_id);
		if(time_s!=step_.EndTime()||*flow_input_>=0.0||graph_donor_)
			throw std::invalid_argument("0D terminal graph donor is out of phase or direction");
		if(concentration.size()!=species_model_.species_ids.size())
			throw std::invalid_argument("0D terminal graph donor species set is incomplete");
		for(const auto& species:species_model_.species_ids){
			const auto value=concentration.find(species);
			if(value==concentration.end()||!(value->second>=0.0)
				||!std::isfinite(value->second))
				throw std::invalid_argument("0D terminal graph donor concentration is invalid");
		}
		graph_donor_=concentration;
	}

	void SolveTransportTrial()override
	{
		RequirePhase(Phase::HydraulicSolved,"solve transport trial");
		if(*flow_input_<0.0&&!graph_donor_)
			throw std::runtime_error("0D terminal graph inflow needs donor concentration");
		auto combined=CombinedTrial();
		transport_trial_=std::move(combined.transport);
		phase_=Phase::TransportSolved;
	}
	PortState GetTransportPortState(const std::string& port_id)const override
	{
		RequirePort(port_id);
		if(phase_==Phase::TransportSolved||phase_==Phase::Prepared)
			return MergePort(transport_trial_->ports.at("graph"));
		if(phase_!=Phase::HydraulicSolved||*flow_input_<0.0)
			throw std::runtime_error("0D terminal outgoing donor state is unavailable");
		return MergePort(CombinedTrial().transport.ports.at("graph"));
	}
	void RollbackTransportTrial()override
	{
		RequirePhase(Phase::TransportSolved,"rollback transport trial");
		transport_trial_.reset();graph_donor_.reset();phase_=Phase::HydraulicSolved;
	}
	std::map<std::string,SpeciesStepAccounting> GetSpeciesStepAccounting()const override
	{
		if(phase_!=Phase::TransportSolved&&phase_!=Phase::Prepared)
			throw std::runtime_error("0D terminal species accounting is unavailable");
		auto result=transport_trial_->balance.species;
		for(auto& species:result){
			auto node=species.second.outward_port_amount.extract("graph");
			if(node.empty())throw std::runtime_error("0D terminal graph amount is missing");
			node.key()="port";
			species.second.outward_port_amount.insert(std::move(node));
		}
		return result;
	}

	void AbortStep()override
	{
		if(phase_==Phase::Committed)
			throw std::runtime_error("0D terminal abort requires an active step");
		ClearTrial();phase_=Phase::Committed;
	}
	void PrepareCommitStep()override
	{
		RequirePhase(Phase::TransportSolved,"prepare commit");
		ValidateZeroDFlowState(hydraulic_trial_->state);
		ValidateZeroDSpeciesReservoirState(species_model_,transport_trial_->state);
		prepared_=Prepared{hydraulic_trial_->state,transport_trial_->state};
		phase_=Phase::Prepared;
	}
	void FinalizeCommitStep()noexcept override
	{
		if(phase_!=Phase::Prepared||!prepared_)std::terminate();
		hydraulic_committed_=prepared_->hydraulic;
		species_committed_.volume_m3=prepared_->species.volume_m3;
		species_committed_.amount_mol.swap(prepared_->species.amount_mol);
		committed_time_s_=step_.EndTime();committed_step_index_=step_.step_index;
		committed_step_=step_;
		ClearTrial();phase_=Phase::Committed;
	}

private:
	enum class Phase{Committed,TrialReady,HydraulicSolved,TransportSolved,Prepared};
	struct Prepared
	{
		ZeroDFlowState hydraulic;
		ZeroDSpeciesReservoirState species;
	};
	void RequirePhase(Phase expected,const char* operation)const
	{
		if(phase_!=expected)
			throw std::runtime_error(std::string("0D terminal ")+operation+" is out of phase");
	}
	void RequirePort(const std::string& port_id)const
	{
		if(port_id!="port")throw std::invalid_argument("0D terminal species port is unknown");
	}
	std::string ModelIdentitySha256()const
	{
		return BuildZeroDFlowSpeciesModelIdentitySha256(hydraulic_model_,
			species_model_,distal_donor_);
	}
	ZeroDTerminalRcrSpeciesTrial CombinedTrial()const
	{
		const double distal_flow=hydraulic_trial_->storage.distal_sink_amount_m3/step_.dt_s;
		const std::map<std::string,double> empty;
		return EvaluateZeroDTerminalRcrSpeciesTrial(hydraulic_model_,
			hydraulic_committed_,species_model_,species_committed_,*flow_input_,
			graph_donor_?*graph_donor_:empty,distal_flow<0.0?distal_donor_:empty,
			step_.dt_s,step_.start_time_s);
	}
	PortState MergePort(PortState state)const
	{
		state.mean_pressure_pa=hydraulic_trial_->port.mean_pressure_pa;
		return state;
	}
	void ClearTrial()noexcept
	{
		flow_input_.reset();graph_donor_.reset();hydraulic_trial_.reset();
		transport_trial_.reset();prepared_.reset();
	}
	std::string domain_id_;
	ZeroDFlowModel hydraulic_model_;
	ZeroDSpeciesReservoirModel species_model_;
	ZeroDFlowState hydraulic_committed_;
	ZeroDSpeciesReservoirState species_committed_;
	std::map<std::string,double> distal_donor_;
	std::vector<CouplingPort> ports_;
	DomainStepContext step_;
	DomainStepContext committed_step_;
	double committed_time_s_=0.0;
	int committed_step_index_=-1;
	Phase phase_=Phase::Committed;
	std::optional<double> flow_input_;
	std::optional<std::map<std::string,double>> graph_donor_;
	std::optional<ZeroDFlowTrial> hydraulic_trial_;
	std::optional<ZeroDSpeciesReservoirTrial> transport_trial_;
	std::optional<Prepared> prepared_;
};

} // namespace iga

#endif
