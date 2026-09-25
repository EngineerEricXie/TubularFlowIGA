#ifndef IGA_NATIVE_TET_ALE_FLOW_DOMAIN_ADAPTER_HPP
#define IGA_NATIVE_TET_ALE_FLOW_DOMAIN_ADAPTER_HPP

#include "NativeTetAleGraphPorts.hpp"
#include "NativeTetAleFlowCheckpoint.hpp"
#include "NativeTetAlePetscRuntime.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

// Hydraulic-only graph owner. FEM basis, assembly and nonlinear residuals
// remain in the native tetra implementation; PETSc supplies algebra/MPI.
class NativeTetAleFlowDomainAdapter final : public CoupledDomainRuntime
{
public:
	using MeshProvider=std::function<NativeTetMesh(const DomainStepContext&,
		const NativeTetMesh&)>;

	NativeTetAleFlowDomainAdapter(std::string domain_id,NativeTetMesh reference_mesh,
		std::vector<CouplingPort> ports,NativeNavierStokesParameters parameters,
		std::vector<double> initial_flow_state,MeshProvider mesh_provider={},
		std::string motion_model_id={},
		std::set<int> no_slip_wall_labels={},double nonlinear_tolerance=1e-8,
		std::size_t maximum_iterations=12)
		:domain_id_(std::move(domain_id)),ports_(std::move(ports)),
		 parameters_(parameters),reference_mesh_(reference_mesh),
		 committed_mesh_(std::move(reference_mesh)),
		 committed_flow_state_(std::move(initial_flow_state)),
		 graph_ports_(domain_id_,committed_mesh_,ports_),
		 mesh_provider_(std::move(mesh_provider)),
		 motion_model_id_(std::move(motion_model_id)),
		 wall_labels_(std::move(no_slip_wall_labels)),
		 nonlinear_tolerance_(nonlinear_tolerance),maximum_iterations_(maximum_iterations)
	{
		const auto topology=BuildNativeTaylorHoodTopology(committed_mesh_);
		const std::size_t velocity_nodes=committed_mesh_.points.size()+topology.edges.size();
		if(committed_flow_state_.size()!=3*velocity_nodes+committed_mesh_.points.size()
			||!(parameters_.density>0.)||!(parameters_.dynamic_viscosity>0.)
			||!std::isfinite(parameters_.density)
			||!std::isfinite(parameters_.dynamic_viscosity)
			||!(nonlinear_tolerance_>0.)||!std::isfinite(nonlinear_tolerance_)
			||maximum_iterations_==0)
			throw std::invalid_argument("native tetra flow graph model is invalid");
		for(const double value:committed_flow_state_)
			if(!std::isfinite(value))
				throw std::invalid_argument("native tetra flow graph initial state is nonfinite");
		for(const int label:wall_labels_){
			if(topology.boundary_velocity_nodes.count(label)==0)
				throw std::invalid_argument("native tetra flow graph wall label is absent");
			for(const auto& port:ports_)
				if(ParseThreeDFlowBoundaryLabel(port)==label)
					throw std::invalid_argument("native tetra flow graph wall is a port");
		}
		if(mesh_provider_&&motion_model_id_.empty())
			throw std::invalid_argument("native tetra flow graph motion model id is required");
		model_identity_sha256_=NativeTetAleFlowModelIdentitySha256(domain_id_,
			reference_mesh_,ports_,parameters_,wall_labels_,motion_model_id_,
			nonlinear_tolerance_,maximum_iterations_);
	}

	const std::string& DomainId()const noexcept override{return domain_id_;}
	DomainKind Kind()const noexcept override{return DomainKind::ThreeDBodyFittedFlow;}
	const std::vector<CouplingPort>& Ports()const noexcept override{return ports_;}
	const NativeTetMesh& ReferenceMesh()const noexcept{return reference_mesh_;}
	const NativeTetMesh& CommittedMesh()const noexcept{return committed_mesh_;}
	const std::vector<double>& CommittedFlowState()const noexcept{return committed_flow_state_;}
	double CommittedTime()const noexcept{return committed_time_s_;}
	std::uint64_t AcceptedSteps()const noexcept{return accepted_steps_;}
	const std::string& ModelIdentitySha256()const noexcept{return model_identity_sha256_;}

	NativeTetAleFlowCheckpoint CaptureCheckpoint()const
	{
		RequirePhase(Phase::Committed,"capture checkpoint");
		if(accepted_steps_==0||!SameTopology(committed_mesh_,reference_mesh_))
			throw std::runtime_error("native tetra flow graph needs an accepted state");
		NativeTetAleFlowCheckpoint value;
		value.model_identity_sha256=model_identity_sha256_;
		value.accepted_steps=accepted_steps_;value.time_s=committed_time_s_;
		value.current_points_m=committed_mesh_.points;
		value.flow_state=committed_flow_state_;
		value.state_identity_sha256=
			native_tet_ale_flow_checkpoint_detail::StateIdentity(value);
		native_tet_ale_flow_checkpoint_detail::Validate(value);
		return value;
	}

	void RestoreCheckpoint(const NativeTetAleFlowCheckpoint& value)
	{
		RequirePhase(Phase::Committed,"restore checkpoint");
		if(accepted_steps_!=0)
			throw std::runtime_error("native tetra flow graph restore requires a fresh runtime");
		native_tet_ale_flow_checkpoint_detail::Validate(value);
		if(value.model_identity_sha256!=model_identity_sha256_
			||value.current_points_m.size()!=reference_mesh_.points.size()
			||value.flow_state.size()!=committed_flow_state_.size())
			throw std::invalid_argument("native tetra flow graph checkpoint model mismatch");
		auto mesh=reference_mesh_;
		mesh.points=value.current_points_m;
		for(const auto& cell:mesh.cells)EvaluateNativeTetGeometry(mesh,cell);
		auto flow=value.flow_state;
		using std::swap;
		swap(committed_mesh_,mesh);
		committed_flow_state_.swap(flow);
		committed_time_s_=value.time_s;
		accepted_steps_=value.accepted_steps;
	}

	void BeginStep(const DomainStepContext& step)override
	{
		RequirePhase(Phase::Committed,"begin step");step.Validate();
		if(accepted_steps_>static_cast<std::uint64_t>(std::numeric_limits<int>::max())
			||step.step_index!=static_cast<int>(accepted_steps_)
			||std::abs(step.start_time_s-committed_time_s_)>
				1e-12*std::max({1.,std::abs(step.start_time_s),
					std::abs(committed_time_s_)}))
			throw std::invalid_argument("native tetra flow graph step clock differs");
		auto mesh=mesh_provider_?mesh_provider_(step,committed_mesh_):committed_mesh_;
		if(!SameTopology(mesh,committed_mesh_))
			throw std::invalid_argument("native tetra flow graph trial topology changed");
		std::vector<std::array<double,3>> velocity(mesh.points.size());
		for(std::size_t node=0;node<mesh.points.size();++node)
			for(int axis=0;axis<3;++axis){
				velocity[node][axis]=(mesh.points[node][axis]
					-committed_mesh_.points[node][axis])/step.dt_s;
				if(!std::isfinite(velocity[node][axis]))
					throw std::invalid_argument("native tetra flow graph mesh velocity is nonfinite");
			}
		for(const auto& cell:mesh.cells)EvaluateNativeTetGeometry(mesh,cell);
		trial_mesh_=std::move(mesh);mesh_velocity_=std::move(velocity);
		step_=step;inputs_.clear();phase_=Phase::Open;
	}

	void SetPortInput(const std::string& port_id,const PortBoundaryData& input)override
	{
		RequirePhase(Phase::Open,"set input");
		if(!HasPort(port_id)||inputs_.count(port_id))
			throw std::invalid_argument("native tetra flow graph input port is unknown or duplicate");
		ValidatePortBoundaryData(input);
		if(!input.concentration.empty()||!input.outward_species_flux.empty())
			throw std::invalid_argument("native tetra flow graph input contains species");
		inputs_.emplace(port_id,input);
	}

	void SolveTrial()override
	{
		RequirePhase(Phase::Open,"solve trial");
		const auto conditions=graph_ports_.BoundaryConditions(trial_mesh_,
			mesh_velocity_,inputs_,step_.EndTime());
		const auto result=SolveNativeTetAlePetscTransient(trial_mesh_,mesh_velocity_,
			committed_flow_state_,committed_flow_state_,WallVelocity(),
			std::numeric_limits<std::uint32_t>::max(),parameters_,step_.dt_s,
			nonlinear_tolerance_,maximum_iterations_,conditions);
		auto states=graph_ports_.Observe(trial_mesh_,mesh_velocity_,
			result.replicated_state,step_.EndTime());
		trial_result_=result;trial_states_=std::move(states);phase_=Phase::Solved;
	}

	PortState GetPortState(const std::string& port_id)const override
	{
		if(phase_!=Phase::Solved&&phase_!=Phase::Prepared)
			throw std::runtime_error("native tetra flow graph trial state is unavailable");
		return trial_states_.at(port_id);
	}

	void RollbackTrial()override
	{
		RequirePhase(Phase::Solved,"rollback trial");
		trial_result_.reset();trial_states_.clear();inputs_.clear();phase_=Phase::Open;
	}

	void AbortStep()override
	{
		if(phase_==Phase::Committed)
			throw std::runtime_error("native tetra flow graph abort needs an active step");
		ClearTrial();phase_=Phase::Committed;
	}

	void PrepareCommitStep()override
	{
		RequirePhase(Phase::Solved,"prepare commit");
		Prepared prepared;
		prepared.mesh=trial_mesh_;
		prepared.flow=trial_result_->replicated_state;
		prepared.time_s=step_.EndTime();
		prepared.accepted_steps=accepted_steps_+1;
		prepared_=std::move(prepared);phase_=Phase::Prepared;
	}

	void FinalizeCommitStep()noexcept override
	{
		if(phase_!=Phase::Prepared||!prepared_)std::terminate();
		using std::swap;
		swap(committed_mesh_,prepared_->mesh);
		committed_flow_state_.swap(prepared_->flow);
		committed_time_s_=prepared_->time_s;
		accepted_steps_=prepared_->accepted_steps;
		ClearTrial();phase_=Phase::Committed;
	}

private:
	enum class Phase{Committed,Open,Solved,Prepared};
	struct Prepared
	{
		NativeTetMesh mesh;
		std::vector<double> flow;
		double time_s=0.;
		std::uint64_t accepted_steps=0;
	};
	static bool SameTopology(const NativeTetMesh& a,const NativeTetMesh& b)
	{
		if(a.points.size()!=b.points.size()||a.cells.size()!=b.cells.size()
			||a.boundary_triangles.size()!=b.boundary_triangles.size())return false;
		for(std::size_t i=0;i<a.cells.size();++i)
			if(a.cells[i].id!=b.cells[i].id||a.cells[i].nodes!=b.cells[i].nodes)
				return false;
		for(std::size_t i=0;i<a.boundary_triangles.size();++i)
			if(a.boundary_triangles[i].id!=b.boundary_triangles[i].id
				||a.boundary_triangles[i].nodes!=b.boundary_triangles[i].nodes
				||a.boundary_triangles[i].boundary_label
					!=b.boundary_triangles[i].boundary_label)return false;
		return true;
	}
	void RequirePhase(Phase phase,const char* operation)const
	{
		if(phase_!=phase)
			throw std::runtime_error(std::string("native tetra flow graph ")
				+operation+" is out of phase");
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
		inputs_.clear();trial_states_.clear();trial_result_.reset();prepared_.reset();
		trial_mesh_=NativeTetMesh{};mesh_velocity_.clear();
	}
	std::string domain_id_;
	std::vector<CouplingPort> ports_;
	NativeNavierStokesParameters parameters_;
	NativeTetMesh reference_mesh_,committed_mesh_,trial_mesh_;
	std::vector<double> committed_flow_state_;
	NativeTetAleGraphPorts graph_ports_;
	MeshProvider mesh_provider_;
	std::string motion_model_id_,model_identity_sha256_;
	std::set<int> wall_labels_;
	double nonlinear_tolerance_;
	std::size_t maximum_iterations_;
	double committed_time_s_=0.;
	std::uint64_t accepted_steps_=0;
	DomainStepContext step_;
	std::vector<std::array<double,3>> mesh_velocity_;
	std::map<std::string,PortBoundaryData> inputs_;
	std::map<std::string,PortState> trial_states_;
	std::optional<NativeTetAlePetscSolveResult> trial_result_;
	std::optional<Prepared> prepared_;
	Phase phase_=Phase::Committed;
};

} // namespace iga

#endif
