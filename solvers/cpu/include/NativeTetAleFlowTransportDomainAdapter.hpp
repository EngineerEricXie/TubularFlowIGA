#ifndef IGA_NATIVE_TET_ALE_FLOW_TRANSPORT_DOMAIN_ADAPTER_HPP
#define IGA_NATIVE_TET_ALE_FLOW_TRANSPORT_DOMAIN_ADAPTER_HPP

#include "NativeTetAleGraphPorts.hpp"
#include "NativeTetAleFlowCheckpoint.hpp"
#include "NativeTetAlePetscRuntime.hpp"
#include "NativeTetMovingSpeciesCheckpoint.hpp"
#include "NativeTetMovingSpeciesPorts.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

struct NativeTetAleFlowSpeciesCheckpoint
{
	NativeTetAleFlowCheckpoint flow;
	NativeTetMovingSpeciesCheckpoint species;
};

// Project-owned ALE P2/P1 hydraulics + conservative P1 scalar behind the
// graph's staged transaction API. PETSc is used only for algebra/MPI.
class NativeTetAleFlowTransportDomainAdapter final : public CoupledDomainRuntime,
	public StagedFlowTransportDomainRuntime
{
public:
	using MeshProvider=std::function<NativeTetMesh(const DomainStepContext&,
		const NativeTetMesh&)>;

	NativeTetAleFlowTransportDomainAdapter(std::string domain_id,
		NativeTetMesh reference_mesh,std::vector<CouplingPort> ports,
		std::string species_id,NativeNavierStokesParameters parameters,
		std::vector<double> initial_flow_state,
		std::vector<double> initial_concentration_mol_m3,
		double diffusivity_m2_s=0.,double source_mol_m3_s=0.,
		MeshProvider mesh_provider={},std::set<int> no_slip_wall_labels={},
		double nonlinear_tolerance=1e-8,std::size_t maximum_iterations=12,
		std::string motion_model_id={},bool monotone_species=false,
		bool zero_flow_concentration_fallback=false)
		:domain_id_(std::move(domain_id)),ports_(std::move(ports)),
		 species_id_(std::move(species_id)),parameters_(parameters),
		 committed_flow_state_(std::move(initial_flow_state)),
		 species_state_(InitializeNativeTetMovingSpeciesState(reference_mesh,
			std::move(initial_concentration_mol_m3),diffusivity_m2_s)),
		 hydraulic_ports_(domain_id_,reference_mesh,ports_),
		 species_ports_(domain_id_,species_id_,reference_mesh,ports_),
		 mesh_provider_(std::move(mesh_provider)),wall_labels_(std::move(no_slip_wall_labels)),
		 source_mol_m3_s_(source_mol_m3_s),nonlinear_tolerance_(nonlinear_tolerance),
		 maximum_iterations_(maximum_iterations),motion_model_id_(std::move(motion_model_id)),
		 monotone_species_(monotone_species),
		 zero_flow_concentration_fallback_(zero_flow_concentration_fallback)
	{
		const auto topology=BuildNativeTaylorHoodTopology(reference_mesh);
		const std::size_t velocity_nodes=reference_mesh.points.size()+topology.edges.size();
		if(committed_flow_state_.size()!=3*velocity_nodes+reference_mesh.points.size()
			||!(parameters_.density>0.)||!(parameters_.dynamic_viscosity>0.)
			||!std::isfinite(parameters_.density)
			||!std::isfinite(parameters_.dynamic_viscosity)
			||!std::isfinite(source_mol_m3_s_)
			||!(nonlinear_tolerance_>0.)||!std::isfinite(nonlinear_tolerance_)
			||maximum_iterations_==0)
			throw std::invalid_argument("native tetra staged model is invalid");
		for(const double value:committed_flow_state_)
			if(!std::isfinite(value))
				throw std::invalid_argument("native tetra staged initial flow is nonfinite");
		for(const int label:wall_labels_){
			if(topology.boundary_velocity_nodes.count(label)==0)
				throw std::invalid_argument("native tetra staged wall label is absent");
			for(const auto& port:ports_)
				if(ParseThreeDFlowBoundaryLabel(port)==label)
					throw std::invalid_argument("native tetra staged wall is a graph port");
		}
		const auto flow_id=NativeTetAleFlowModelIdentitySha256(domain_id_,
			reference_mesh,ports_,parameters_,wall_labels_,motion_model_id_,
			nonlinear_tolerance_,maximum_iterations_);
		Sha256 hash;
		const auto append_text=[&hash](const std::string& value){
			hash.AppendLittleEndian64(value.size());
			hash.Append(value.data(),value.size());
		};
		append_text("NativeTetAleFlowSpeciesModel/v1");
		append_text(flow_id);
		append_text(species_state_.model_identity_sha256);
		append_text(species_id_);
		hash.AppendNormalizedDouble(source_mol_m3_s_);
		// Preserve the existing Galerkin identity while separating opt-in
		// low-order checkpoint/model compatibility.
		if(monotone_species_)append_text("monotone-graph-diffusion/v1");
		model_identity_sha256_=hash.Hex();
	}

	const std::string& DomainId()const noexcept override{return domain_id_;}
	DomainKind Kind()const noexcept override{return DomainKind::ThreeDBodyFittedFlow;}
	const std::vector<CouplingPort>& Ports()const noexcept override{return ports_;}
	const std::vector<double>& CommittedFlowState()const noexcept{return committed_flow_state_;}
	const NativeTetMovingSpeciesCommittedState& CommittedSpeciesState()const noexcept
	{return species_state_;}
	const std::string& CheckpointModelIdentitySha256()const noexcept
	{return model_identity_sha256_;}
	NativeTetAleFlowSpeciesCheckpoint CaptureCheckpoints()const
	{
		RequirePhase(Phase::Committed,"capture checkpoints");
		if(motion_model_id_.empty()||species_state_.accepted_steps==0)
			throw std::runtime_error("native tetra species checkpoint needs accepted state and motion id");
		NativeTetAleFlowSpeciesCheckpoint result;
		result.species=CaptureNativeTetMovingSpeciesCheckpoint(species_state_);
		result.flow.model_identity_sha256=model_identity_sha256_;
		result.flow.accepted_steps=species_state_.accepted_steps;
		result.flow.time_s=species_state_.time_s;
		result.flow.current_points_m=species_state_.current_mesh.points;
		result.flow.flow_state=committed_flow_state_;
		result.flow.state_identity_sha256=
			native_tet_ale_flow_checkpoint_detail::StateIdentity(result.flow);
		native_tet_ale_flow_checkpoint_detail::Validate(result.flow);
		return result;
	}
	void RestoreCheckpoints(const NativeTetAleFlowSpeciesCheckpoint& value)
	{
		RequirePhase(Phase::Committed,"restore checkpoints");
		if(motion_model_id_.empty()||species_state_.accepted_steps!=0)
			throw std::runtime_error("native tetra species restore needs fresh state and motion id");
		native_tet_ale_flow_checkpoint_detail::Validate(value.flow);
		native_tet_moving_species_checkpoint_detail::Validate(value.species);
		if(value.flow.model_identity_sha256!=model_identity_sha256_
			||value.flow.accepted_steps!=value.species.accepted_steps
			||value.flow.time_s!=value.species.time_s
			||value.flow.current_points_m!=value.species.current_points_m
			||value.flow.flow_state.size()!=committed_flow_state_.size())
			throw std::invalid_argument("native tetra flow/species checkpoint pair differs");
		auto restored=RestoreNativeTetMovingSpeciesCheckpoint(value.species,
			species_state_.reference_mesh,species_state_.diffusivity_m2_s);
		auto flow=value.flow.flow_state;
		using std::swap;
		swap(species_state_,restored);
		committed_flow_state_.swap(flow);
	}

	void BeginStep(const DomainStepContext& step)override
	{
		RequirePhase(Phase::Committed,"begin step");step.Validate();
		if(species_state_.accepted_steps>static_cast<std::uint64_t>(
			std::numeric_limits<int>::max())
			||step.step_index!=static_cast<int>(species_state_.accepted_steps)
			||std::abs(step.start_time_s-species_state_.time_s)>
				1e-12*std::max({1.,std::abs(step.start_time_s),
					std::abs(species_state_.time_s)}))
			throw std::invalid_argument("native tetra staged step clock differs from committed state");
		auto trial_mesh=mesh_provider_
			?mesh_provider_(step,species_state_.current_mesh)
			:species_state_.current_mesh;
		if(!native_tet_moving_species_checkpoint_detail::SameTopology(
			trial_mesh,species_state_.reference_mesh))
			throw std::invalid_argument("native tetra staged trial topology changed");
		std::vector<std::array<double,3>> velocity(trial_mesh.points.size());
		for(std::size_t node=0;node<velocity.size();++node)
			for(int axis=0;axis<3;++axis){
				velocity[node][axis]=(trial_mesh.points[node][axis]
					-species_state_.current_mesh.points[node][axis])/step.dt_s;
				if(!std::isfinite(velocity[node][axis]))
					throw std::invalid_argument("native tetra staged mesh velocity is nonfinite");
			}
		for(const auto& cell:trial_mesh.cells)EvaluateNativeTetGeometry(trial_mesh,cell);
		trial_mesh_=std::move(trial_mesh);
		mesh_velocity_=std::move(velocity);
		step_=step;hydraulic_inputs_.clear();concentration_inputs_.clear();
		phase_=Phase::TrialReady;
	}

	void SetPortInput(const std::string& port_id,const PortBoundaryData& input)override
	{
		RequirePhase(Phase::TrialReady,"set hydraulic input");
		if(!HasPort(port_id)||hydraulic_inputs_.count(port_id))
			throw std::invalid_argument("native tetra staged hydraulic port is unknown or duplicate");
		ValidatePortBoundaryData(input);
		if(!input.concentration.empty()||!input.outward_species_flux.empty())
			throw std::invalid_argument("native tetra staged hydraulic input contains species");
		hydraulic_inputs_.emplace(port_id,input);
	}

	void SolveTrial()override{SolveHydraulicTrial();}
	void SolveHydraulicTrial()override
	{
		RequirePhase(Phase::TrialReady,"solve hydraulic trial");
		const auto conditions=hydraulic_ports_.BoundaryConditions(trial_mesh_,
			mesh_velocity_,hydraulic_inputs_,step_.EndTime());
		const auto prescribed=WallVelocity();
		const auto result=SolveNativeTetAlePetscTransient(trial_mesh_,mesh_velocity_,
			committed_flow_state_,committed_flow_state_,prescribed,
			std::numeric_limits<std::uint32_t>::max(),parameters_,step_.dt_s,
			nonlinear_tolerance_,maximum_iterations_,conditions);
		auto states=hydraulic_ports_.Observe(trial_mesh_,mesh_velocity_,
			result.replicated_state,step_.EndTime());
		hydraulic_result_=result;hydraulic_states_=std::move(states);
		phase_=Phase::HydraulicSolved;
	}

	PortState GetPortState(const std::string& port_id)const override
	{
		if(phase_==Phase::TransportSolved||phase_==Phase::Prepared)
			return transport_states_.at(port_id);
		return GetHydraulicPortState(port_id);
	}
	PortState GetHydraulicPortState(const std::string& port_id)const override
	{
		if(phase_!=Phase::HydraulicSolved&&phase_!=Phase::TransportSolved
			&&phase_!=Phase::Prepared)
			throw std::runtime_error("native tetra staged hydraulic state is unavailable");
		return hydraulic_states_.at(port_id);
	}

	void RollbackTrial()override{RollbackHydraulicTrial();}
	void RollbackHydraulicTrial()override
	{
		RequirePhase(Phase::HydraulicSolved,"rollback hydraulic trial");
		hydraulic_result_.reset();hydraulic_states_.clear();hydraulic_inputs_.clear();
		phase_=Phase::TrialReady;
	}

	void SetTransportConcentration(const std::string& port_id,double time_s,
		const std::map<std::string,double>& concentration)override
	{
		RequirePhase(Phase::HydraulicSolved,"set transport concentration");
		if(!HasPort(port_id)||concentration_inputs_.count(port_id)
			||concentration.size()!=1||concentration.count(species_id_)!=1
			||!std::isfinite(concentration.at(species_id_))
			||concentration.at(species_id_)<0.
			||std::abs(time_s-step_.EndTime())>
				1e-12*std::max({1.,std::abs(time_s),std::abs(step_.EndTime())}))
			throw std::invalid_argument("native tetra staged species input is invalid");
		PortBoundaryData input;input.time_s=time_s;input.concentration=concentration;
		concentration_inputs_.emplace(port_id,std::move(input));
	}

	void SolveTransportTrial()override
	{
		RequirePhase(Phase::HydraulicSolved,"solve transport trial");
		std::map<std::string,double> flows;
		for(const auto& port:ports_)
			flows.emplace(port.id,*hydraulic_states_.at(port.id).outward_flow_m3_s);
		auto inlet=species_ports_.InflowByLabel(flows,
			concentration_inputs_,step_.EndTime(),flow_epsilon_m3_s_);
		if(zero_flow_concentration_fallback_){
			const double zero_flow_concentration=std::accumulate(
				species_state_.concentration_mol_m3.begin(),
				species_state_.concentration_mol_m3.end(),0.)
				/species_state_.concentration_mol_m3.size();
			for(const auto& port:ports_)
				inlet.try_emplace(ParseThreeDFlowBoundaryLabel(port),zero_flow_concentration);
		}
		const auto topology=BuildNativeTaylorHoodTopology(trial_mesh_);
		const std::size_t velocity_nodes=trial_mesh_.points.size()+topology.edges.size();
		std::vector<std::array<double,3>> velocity(velocity_nodes);
		for(std::size_t node=0;node<velocity_nodes;++node)
			for(int axis=0;axis<3;++axis)
				velocity[node][axis]=hydraulic_result_->replicated_state[3*node+axis];
		const auto result=SolveNativeTetMovingSpeciesPetscStep(
			species_state_.current_mesh,trial_mesh_,velocity,mesh_velocity_,
			species_state_.concentration_mol_m3,inlet,
			species_state_.diffusivity_m2_s,source_mol_m3_s_,step_.dt_s,
			monotone_species_);
		auto scalar_states=species_ports_.Observe(trial_mesh_,result.step,
			step_.EndTime(),flow_epsilon_m3_s_);
		for(const auto& port:ports_){
			auto& state=scalar_states.at(port.id);
			const auto& hydraulic=hydraulic_states_.at(port.id);
			const double flow=*hydraulic.outward_flow_m3_s;
			if(std::abs(*state.outward_flow_m3_s-flow)>
				1e-10*std::max({1.,std::abs(flow),
					std::abs(*state.outward_flow_m3_s)}))
				throw std::runtime_error("native tetra staged hydraulic/species port flow mismatch");
			state.area_m2=hydraulic.area_m2;
			state.mean_pressure_pa=hydraulic.mean_pressure_pa;
		}
		const auto accounting=species_ports_.Accounting(result.step,step_.dt_s);
		const double scale=std::max({1.,std::abs(accounting.initial_mass),
			std::abs(accounting.final_mass),std::abs(accounting.source_amount)});
		if(!std::isfinite(accounting.residual)
			||std::abs(accounting.residual)>1e-10*scale)
			throw std::runtime_error("native tetra staged species amount balance failed");
		transport_result_=result;transport_states_=std::move(scalar_states);
		accounting_=accounting;phase_=Phase::TransportSolved;
	}

	PortState GetTransportPortState(const std::string& port_id)const override
	{
		if(phase_!=Phase::TransportSolved&&phase_!=Phase::Prepared)
			throw std::runtime_error("native tetra staged transport state is unavailable");
		return transport_states_.at(port_id);
	}
	void RollbackTransportTrial()override
	{
		RequirePhase(Phase::TransportSolved,"rollback transport trial");
		transport_result_.reset();transport_states_.clear();accounting_.reset();
		concentration_inputs_.clear();phase_=Phase::HydraulicSolved;
	}
	std::map<std::string,SpeciesStepAccounting> GetSpeciesStepAccounting()const override
	{
		if((phase_!=Phase::TransportSolved&&phase_!=Phase::Prepared)||!accounting_)
			throw std::runtime_error("native tetra staged species accounting is unavailable");
		return {{species_id_,*accounting_}};
	}

	void AbortStep()override
	{
		if(phase_==Phase::Committed)
			throw std::runtime_error("native tetra staged abort needs an active step");
		ClearTrial();phase_=Phase::Committed;
	}
	void PrepareCommitStep()override
	{
		RequirePhase(Phase::TransportSolved,"prepare commit");
		Prepared next;
		next.mesh=trial_mesh_;
		next.flow=hydraulic_result_->replicated_state;
		next.concentration=transport_result_->step.concentration_mol_m3;
		next.time_s=step_.EndTime();
		next.accepted_steps=species_state_.accepted_steps+1;
		prepared_.emplace(std::move(next));phase_=Phase::Prepared;
	}
	void FinalizeCommitStep()noexcept override
	{
		if(phase_!=Phase::Prepared||!prepared_)std::terminate();
		using std::swap;
		swap(species_state_.current_mesh,prepared_->mesh);
		committed_flow_state_.swap(prepared_->flow);
		species_state_.concentration_mol_m3.swap(prepared_->concentration);
		species_state_.time_s=prepared_->time_s;
		species_state_.accepted_steps=prepared_->accepted_steps;
		ClearTrial();phase_=Phase::Committed;
	}

private:
	enum class Phase{Committed,TrialReady,HydraulicSolved,TransportSolved,Prepared};
	struct Prepared
	{
		NativeTetMesh mesh;
		std::vector<double> flow,concentration;
		double time_s=0.;
		std::uint64_t accepted_steps=0;
	};
	void RequirePhase(Phase expected,const char* operation)const
	{
		if(phase_!=expected)
			throw std::runtime_error(std::string("native tetra staged ")+operation
				+" is out of phase");
	}
	bool HasPort(const std::string& id)const
	{
		for(const auto& port:ports_)if(port.id==id)return true;
		return false;
	}
	std::map<std::uint32_t,std::array<double,3>> WallVelocity()const
	{
		std::map<std::uint32_t,std::array<double,3>> result;
		const auto topology=BuildNativeTaylorHoodTopology(trial_mesh_);
		std::vector<std::array<double,3>> p2=mesh_velocity_;
		for(const auto& edge:topology.edges){
			std::array<double,3> average{};
			for(int axis=0;axis<3;++axis)
				average[axis]=0.5*(mesh_velocity_[edge[0]][axis]
					+mesh_velocity_[edge[1]][axis]);
			p2.push_back(average);
		}
		for(const int label:wall_labels_)
			for(const auto node:topology.boundary_velocity_nodes.at(label))
				result.emplace(node,p2.at(node));
		return result;
	}
	void ClearTrial()noexcept
	{
		hydraulic_inputs_.clear();concentration_inputs_.clear();
		hydraulic_states_.clear();transport_states_.clear();
		hydraulic_result_.reset();transport_result_.reset();accounting_.reset();
		prepared_.reset();trial_mesh_=NativeTetMesh{};mesh_velocity_.clear();
	}
	std::string domain_id_;
	std::vector<CouplingPort> ports_;
	std::string species_id_;
	NativeNavierStokesParameters parameters_;
	std::vector<double> committed_flow_state_;
	NativeTetMovingSpeciesCommittedState species_state_;
	NativeTetAleGraphPorts hydraulic_ports_;
	NativeTetMovingSpeciesPorts species_ports_;
	MeshProvider mesh_provider_;
	std::set<int> wall_labels_;
	double source_mol_m3_s_=0.,nonlinear_tolerance_=1e-8;
	std::size_t maximum_iterations_=12;
	std::string motion_model_id_,model_identity_sha256_;
	bool monotone_species_=false;
	bool zero_flow_concentration_fallback_=false;
	double flow_epsilon_m3_s_=1e-10;
	Phase phase_=Phase::Committed;
	DomainStepContext step_;
	NativeTetMesh trial_mesh_;
	std::vector<std::array<double,3>> mesh_velocity_;
	std::map<std::string,PortBoundaryData> hydraulic_inputs_,concentration_inputs_;
	std::optional<NativeTetAlePetscSolveResult> hydraulic_result_;
	std::optional<NativeTetMovingSpeciesPetscResult> transport_result_;
	std::map<std::string,PortState> hydraulic_states_,transport_states_;
	std::optional<SpeciesStepAccounting> accounting_;
	std::optional<Prepared> prepared_;
};

} // namespace iga

#endif
