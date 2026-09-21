#ifndef IGA_NATIVE_TET_ALE_FSI_RUNTIME_HPP
#define IGA_NATIVE_TET_ALE_FSI_RUNTIME_HPP

// Single-partition reference adapter connecting native ALE mesh motion,
// native P2/P1 flow, matching traction transfer, and the shared FSI lifecycle.
#include "FsiDomainRuntime.hpp"
#include "NativeTetAleDenseRuntime.hpp"
#include "NativeTetAleKinematics.hpp"
#include "NativeTetAleMeshMotion.hpp"
#include "NativeTetMatchingFsiInterface.hpp"
#include "StrongFluidStructureCoupling.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

class NativeTetAleFsiRuntime final : public FsiFluidDomainRuntime {
public:
	NativeTetAleFsiRuntime(const NativeTetAleFsiRuntime&)=delete;
	NativeTetAleFsiRuntime& operator=(const NativeTetAleFsiRuntime&)=delete;
	NativeTetAleFsiRuntime(NativeTetAleFsiRuntime&&)=delete;
	NativeTetAleFsiRuntime& operator=(NativeTetAleFsiRuntime&&)=delete;

	NativeTetAleFsiRuntime(std::string domain_id,std::string subsystem_id,FsiCouplingEdge edge,
		DistributedSurfaceInterface fluid_surface,DistributedSurfaceInterface structure_surface,
		DistributedSurfaceLayout fluid_layout,DistributedSurfaceLayout structure_layout,
		NativeTetMesh reference_mesh,int interface_label,NativeNavierStokesParameters parameters,
		std::uint32_t pressure_gauge_node,std::vector<double> initial_fluid_state={},
		double initial_time_s=0.0,double nonlinear_tolerance=1e-8,std::size_t maximum_iterations=12,
		std::map<std::uint32_t,std::array<double,3>> prescribed_noninterface_velocity={},
		std::set<int> natural_velocity_labels={},
		std::vector<std::array<double,3>> initial_ale_displacement_m={},
		std::vector<std::array<double,3>> initial_interface_displacement_m={})
		:domain_id_(std::move(domain_id)),subsystem_id_(std::move(subsystem_id)),edge_(std::move(edge)),
		 catalog_{std::move(fluid_surface)},structure_surface_(std::move(structure_surface)),
		 fluid_layout_(std::move(fluid_layout)),structure_layout_(std::move(structure_layout)),
		 reference_mesh_(std::move(reference_mesh)),topology_(BuildNativeTaylorHoodTopology(reference_mesh_)),
		 interface_label_(interface_label),parameters_(parameters),pressure_gauge_node_(pressure_gauge_node),
		 nonlinear_tolerance_(nonlinear_tolerance),maximum_iterations_(maximum_iterations),
		 prescribed_noninterface_velocity_(std::move(prescribed_noninterface_velocity)),
		 natural_velocity_labels_(std::move(natural_velocity_labels)),
		 ale_(reference_mesh_.points,reference_mesh_.cells,initial_time_s,std::move(initial_ale_displacement_m)),
		 lifecycle_(domain_id_,subsystem_id_,edge_,catalog_.front(),structure_surface_,fluid_layout_,structure_layout_)
	{
		if(!(catalog_.front().id==edge_.fluid)||!(structure_surface_.id==edge_.structure)
			||catalog_.front().id.domain_id!=domain_id_||catalog_.front().id.subsystem_id!=subsystem_id_)
			throw std::runtime_error("native tetrahedral ALE FSI endpoints are inconsistent");
		if(fluid_layout_.partition_count!=1||fluid_layout_.partition_rank!=0)
			throw std::runtime_error("native tetrahedral ALE FSI reference adapter is single-partition");
		if(catalog_.front().boundary_labels!=std::vector<std::int64_t>{interface_label_})
			throw std::runtime_error("native tetrahedral ALE FSI label differs from its surface contract");
		if(!(parameters_.density>0.0)||!(parameters_.dynamic_viscosity>0.0)
			||!(nonlinear_tolerance_>0.0)||maximum_iterations_==0
			||pressure_gauge_node_>=reference_mesh_.points.size())
			throw std::runtime_error("native tetrahedral ALE FSI numerical options are invalid");
		if(natural_velocity_labels_.count(interface_label_))
			throw std::runtime_error("native tetrahedral ALE FSI interface cannot be a natural-velocity boundary");
		std::set<int> boundary_labels;
		for(const auto& triangle:reference_mesh_.boundary_triangles)
			boundary_labels.insert(triangle.boundary_label);
		for(int label:natural_velocity_labels_)
			if(label<0||!boundary_labels.count(label))
				throw std::runtime_error("native tetrahedral ALE FSI natural-velocity label is absent");
		const std::size_t velocity_nodes=reference_mesh_.points.size()+topology_.edges.size();
		for(const auto& prescribed:prescribed_noninterface_velocity_){
			if(prescribed.first>=velocity_nodes)throw std::runtime_error("native tetrahedral ALE FSI prescribed velocity node is invalid");
			for(double value:prescribed.second)if(!std::isfinite(value))
				throw std::runtime_error("native tetrahedral ALE FSI prescribed velocity is nonfinite");
		}
		ValidateInterfaceTopology();
		const std::size_t dofs=3*velocity_nodes+reference_mesh_.points.size();
		if(initial_fluid_state.empty())initial_fluid_state.assign(dofs,0.0);
		if(initial_fluid_state.size()!=dofs)throw std::runtime_error("native tetrahedral ALE FSI state size mismatch");
		for(double value:initial_fluid_state)if(!std::isfinite(value))
			throw std::runtime_error("native tetrahedral ALE FSI initial state is nonfinite");
		committed_fluid_state_=std::move(initial_fluid_state);
		if(initial_interface_displacement_m.empty())
			initial_interface_displacement_m.assign(fluid_layout_.owned_global_node_ids.size(),{{0,0,0}});
		if(initial_interface_displacement_m.size()!=fluid_layout_.owned_global_node_ids.size())
			throw std::runtime_error("native tetrahedral ALE FSI initial interface displacement size mismatch");
		for(std::size_t row=0;row<initial_interface_displacement_m.size();++row){
			const auto node=static_cast<std::size_t>(fluid_layout_.owned_global_node_ids[row]);
			for(int component=0;component<3;++component)
				if(!std::isfinite(initial_interface_displacement_m[row][component])
					||initial_interface_displacement_m[row][component]!=ale_.CommittedDisplacement()[node][component])
					throw std::runtime_error("native tetrahedral ALE FSI restored interface and ALE states differ");
		}
		committed_interface_displacement_=std::move(initial_interface_displacement_m);
		model_identity_sha256_=BuildModelIdentity();
		committed_composition_identity_sha256_=BuildCompositionIdentity(committed_fluid_state_,
			committed_interface_displacement_,initial_time_s);
	}

	const std::vector<DistributedSurfaceInterface>& SurfaceInterfaces()const noexcept override{return catalog_;}
	const FsiTrialLifecycle& Lifecycle()const noexcept{return lifecycle_;}
	const FsiCouplingEdge& StrongCouplingEdge()const noexcept{return edge_;}
	const DistributedSurfaceLayout& StrongCouplingLayout()const noexcept{return fluid_layout_;}
	const std::vector<std::array<double,3>>& StrongCouplingCommittedDisplacementM()const noexcept
	{return committed_interface_displacement_;}
	const std::string& CommittedCompositionIdentitySha256()const noexcept
	{return committed_composition_identity_sha256_;}
	const std::string& ModelIdentitySha256()const noexcept{return model_identity_sha256_;}
	const std::vector<double>& CommittedFluidState()const noexcept{return committed_fluid_state_;}
	const std::vector<std::array<double,3>>& CommittedAleDisplacementM()const noexcept
	{return ale_.CommittedDisplacement();}
	const std::vector<std::array<double,3>>& CommittedInterfaceDisplacementM()const noexcept
	{return committed_interface_displacement_;}
	double CommittedAleTimeS()const noexcept{return ale_.CommittedTime();}

	void BeginMacroStep(const DomainStepContext& step)
	{
		if(trial_||prepared_)throw std::runtime_error("native tetrahedral ALE FSI retains stale trial state");
		if(step.start_time_s!=ale_.CommittedTime())
			throw std::runtime_error("native tetrahedral ALE FSI step time differs from committed ALE time");
		lifecycle_.BeginStep(step);
	}
	void BeginCouplingIteration(std::uint64_t iteration,const SurfaceFieldStamp& expected_kinematics,
		const SurfaceFieldStampEnvelope& output_envelope)
	{
		if(trial_||prepared_||ale_.HasTrial())throw std::runtime_error("native tetrahedral ALE FSI retains stale iteration state");
		SurfaceFieldStamp input_copy=expected_kinematics;SurfaceFieldStampEnvelope output_copy=output_envelope;
		lifecycle_.BeginIteration(iteration,input_copy,output_copy);using std::swap;
		swap(expected_kinematics_,input_copy);swap(output_envelope_,output_copy);
	}
	void SetSurfaceKinematics(const std::string& interface_id,const SurfaceKinematics& kinematics)override
	{
		if(interface_id!=edge_.fluid.interface_id)throw std::runtime_error("unknown native tetrahedral ALE FSI interface");
		SurfaceKinematics candidate=kinematics;
		ValidateFsiFluidKinematicsInput(edge_,candidate,fluid_layout_,expected_kinematics_);
		lifecycle_.MarkInput(candidate.interface,candidate.stamp);kinematics_=std::move(candidate);
	}

	void SolveFluidTrial()
	{
		lifecycle_.RequireSolveAllowed();if(!kinematics_)throw std::runtime_error("native tetrahedral ALE FSI input is unavailable");
		const auto boundary_displacement=BuildBoundaryDisplacement(*kinematics_);
		const auto motion=SolveNativeTetAleHarmonicMotion(reference_mesh_,boundary_displacement);
		const auto& geometry=ale_.BeginTrial(motion.displacement_m,lifecycle_.Context().EndTime());
		try{
			const auto current_mesh=BuildNativeTetAleCurrentMesh(reference_mesh_,geometry);
			const auto prescribed_velocity=BuildFluidBoundaryVelocity(geometry);
			NativeTetAleFsiTrial candidate;
			candidate.solve=SolveNativeTetAleDenseTransient(current_mesh,geometry.mesh_velocity_m_s,
				committed_fluid_state_,committed_fluid_state_,prescribed_velocity,pressure_gauge_node_,
				parameters_,geometry.dt_s,nonlinear_tolerance_,maximum_iterations_);
			const auto triangle_traction=EvaluateNativeTetFluidTractionOnMatchingStructure(current_mesh,
				topology_,candidate.solve.state,parameters_.dynamic_viscosity,interface_label_);
			std::vector<double> displacement(3*reference_mesh_.points.size()),velocity(displacement.size());
			for(std::size_t node=0;node<reference_mesh_.points.size();++node)for(int component=0;component<3;++component){
				displacement[3*node+component]=geometry.displacement_m[node][component];
				velocity[3*node+component]=geometry.mesh_velocity_m_s[node][component];}
			candidate.transfer=BuildNativeTetMatchingFsiTransfer(reference_mesh_,topology_,interface_label_,
				displacement,velocity,triangle_traction);
			candidate.traction=BuildTraction(candidate,triangle_traction);
			lifecycle_.MarkSolved(candidate.traction.interface,candidate.traction.stamp);
			trial_.emplace(std::move(candidate));
		}catch(...){ale_.RejectTrial();throw;}
	}
	SurfaceTraction GetSurfaceTraction(const std::string& interface_id)const override
	{
		if(interface_id!=edge_.fluid.interface_id||!trial_)
			throw std::runtime_error("native tetrahedral ALE FSI traction is unavailable");
		lifecycle_.RequireTrialOutput(trial_->traction.interface,trial_->traction.stamp);return trial_->traction;
	}
	SurfaceTraction GetCommittedSurfaceTraction(const std::string& interface_id,
		const SurfaceFieldStamp& stamp)const
	{
		if(interface_id!=edge_.fluid.interface_id)throw std::runtime_error("unknown native tetrahedral ALE FSI interface");
		lifecycle_.RequireCommittedOutput(catalog_.front().id,stamp);return committed_traction_;
	}
	StrongCouplingAcceptanceDiagnostics TrialCouplingAcceptanceDiagnostics()const
	{
		if(!trial_)throw std::runtime_error("native tetrahedral ALE FSI diagnostics require a trial");
		return {trial_->solve.final_free_residual_l2,
			std::abs(trial_->transfer.nodal_power_w-trial_->transfer.quadrature_power_w)};
	}
	void RejectCouplingIteration()
	{
		if(!lifecycle_.HasActiveStep())throw std::runtime_error("native tetrahedral ALE FSI reject requires a step");
		DiscardTrialNoexcept();lifecycle_.RejectIteration();
	}
	void PrepareCommitStep()
	{
		if(!trial_)throw std::runtime_error("native tetrahedral ALE FSI prepare requires a trial");
		NativeTetAleFsiPrepared candidate;candidate.fluid_state=trial_->solve.state;
		candidate.traction=trial_->traction;candidate.interface_displacement=kinematics_->displacement_m;
		candidate.composition_identity_sha256=BuildCompositionIdentity(candidate.fluid_state,
			candidate.interface_displacement,lifecycle_.Context().EndTime());
		try{ale_.PrepareTrialCommit();lifecycle_.PrepareCommit();prepared_.emplace(std::move(candidate));}
		catch(...){
			ale_.RejectTrial();trial_.reset();kinematics_.reset();prepared_.reset();
			if(lifecycle_.HasActiveStep())lifecycle_.RejectIteration();
			throw;
		}
	}
	void AbortStep()noexcept{DiscardTrialNoexcept();lifecycle_.AbortStep();}
	std::string CoordinatorPreparedCommittedTractionIdentitySha256()const
	{if(!prepared_)throw std::runtime_error("native tetrahedral ALE FSI prepared traction unavailable");return BuildSurfaceTractionIdentitySha256(prepared_->traction,fluid_layout_);}
	std::string CoordinatorPreparedCommittedCompositionIdentitySha256()const
	{if(!prepared_)throw std::runtime_error("native tetrahedral ALE FSI prepared state unavailable");return prepared_->composition_identity_sha256;}

private:
	friend class StrongFluidStructureCouplingAccess;
	struct NativeTetAleFsiTrial{NativeTetAleDenseSolveResult solve;NativeTetMatchingFsiTransfer transfer;SurfaceTraction traction;};
	struct NativeTetAleFsiPrepared{std::vector<double> fluid_state;std::vector<std::array<double,3>> interface_displacement;
		SurfaceTraction traction;std::string composition_identity_sha256;};
	static_assert(std::is_nothrow_move_constructible<NativeTetAleFsiPrepared>::value,
		"prepared native ALE FSI state must move without throwing");
	void CoordinatorRequireFinalizeAllowed()const
	{if(!trial_||!prepared_)throw std::runtime_error("native tetrahedral ALE FSI finalize requires prepared state");
		ale_.RequireFinalizePreparedTrialAllowed();lifecycle_.RequireFinalizeAllowed();}
	void CoordinatorFinalizeCommitNoexcept()noexcept
	{
		using std::swap;swap(committed_fluid_state_,prepared_->fluid_state);
		swap(committed_interface_displacement_,prepared_->interface_displacement);
		swap(committed_traction_,prepared_->traction);
		swap(committed_composition_identity_sha256_,prepared_->composition_identity_sha256);
		ale_.FinalizePreparedTrialNoexcept();lifecycle_.FinalizePreparedCommitNoexcept();
		trial_.reset();prepared_.reset();kinematics_.reset();
	}
	void DiscardTrialNoexcept()noexcept
	{ale_.RejectTrial();trial_.reset();prepared_.reset();kinematics_.reset();}

	void ValidateInterfaceTopology()const
	{
		std::set<std::array<std::uint64_t,3>> mesh_faces,layout_faces;
		for(const auto& triangle:reference_mesh_.boundary_triangles)if(triangle.boundary_label==interface_label_){
			std::array<std::uint64_t,3> face{{triangle.nodes[0],triangle.nodes[1],triangle.nodes[2]}};std::sort(face.begin(),face.end());
			if(!mesh_faces.insert(face).second)throw std::runtime_error("native tetrahedral ALE FSI duplicate interface face");}
		for(auto face:fluid_layout_.reference_triangles){std::sort(face.begin(),face.end());
			if(!layout_faces.insert(face).second)throw std::runtime_error("native tetrahedral ALE FSI duplicate layout face");}
		if(mesh_faces.empty()||mesh_faces!=layout_faces)throw std::runtime_error("native tetrahedral ALE FSI layout does not match labelled faces");
		for(const auto& position:fluid_layout_.reference_positions)
			if(position.global_node_id>=reference_mesh_.points.size()
				||reference_mesh_.points[position.global_node_id]!=position.position_m)
				throw std::runtime_error("native tetrahedral ALE FSI reference positions differ from volume mesh");
	}
	std::map<std::uint32_t,std::array<double,3>> BuildBoundaryDisplacement(const SurfaceKinematics& value)const
	{
		std::map<std::uint32_t,std::array<double,3>> result;std::set<std::uint32_t> interface_nodes;
		for(std::size_t row=0;row<fluid_layout_.owned_global_node_ids.size();++row){
			const auto node=static_cast<std::uint32_t>(fluid_layout_.owned_global_node_ids[row]);
			result[node]=value.displacement_m[row];interface_nodes.insert(node);}
		for(const auto& triangle:reference_mesh_.boundary_triangles)if(triangle.boundary_label!=interface_label_)
			for(auto node:triangle.nodes){const auto found=result.find(node);if(found==result.end())result[node]={{0,0,0}};
				else if(interface_nodes.count(node)&&NativeTetInterfaceNorm(found->second)>1e-12)
					throw std::runtime_error("native tetrahedral ALE FSI moving/fixed seam displacement conflicts");}
		return result;
	}
	std::map<std::uint32_t,std::array<double,3>> BuildFluidBoundaryVelocity(const NativeTetAleTrial& geometry)const
	{
		std::vector<double> displacement(3*reference_mesh_.points.size()),velocity(displacement.size());
		for(std::size_t node=0;node<reference_mesh_.points.size();++node)for(int component=0;component<3;++component){
			displacement[3*node+component]=geometry.displacement_m[node][component];
			velocity[3*node+component]=geometry.mesh_velocity_m_s[node][component];}
		std::vector<NativeTetInterfaceTriangleTraction> zero;
		for(const auto& triangle:reference_mesh_.boundary_triangles)if(triangle.boundary_label==interface_label_)
			zero.push_back({triangle.id,{{0,0,0}}});
		auto result=BuildNativeTetMatchingFsiTransfer(reference_mesh_,topology_,interface_label_,
			displacement,velocity,zero).fluid_no_slip_velocity_m_s;
		std::set<std::uint32_t> consumed;
		for(const auto& by_label:topology_.boundary_velocity_nodes)if(by_label.first!=interface_label_
			&&!natural_velocity_labels_.count(by_label.first))
			for(auto node:by_label.second){
				const auto prescribed=prescribed_noninterface_velocity_.find(node);
				if(prescribed==prescribed_noninterface_velocity_.end())
					throw std::runtime_error("native tetrahedral ALE FSI essential boundary node lacks velocity");
				consumed.insert(node);const auto found=result.find(node);
				if(found==result.end())result[node]=prescribed->second;
				else for(int component=0;component<3;++component)
					if(std::abs(found->second[component]-prescribed->second[component])>1e-12)
						throw std::runtime_error("native tetrahedral ALE FSI interface/essential velocity conflicts");
			}
		for(const auto& prescribed:prescribed_noninterface_velocity_)
			if(!consumed.count(prescribed.first))
				throw std::runtime_error("native tetrahedral ALE FSI prescribed velocity is extraneous");
		return result;
	}
	SurfaceTraction BuildTraction(const NativeTetAleFsiTrial& candidate,
		const std::vector<NativeTetInterfaceTriangleTraction>& triangles)const
	{
		SurfaceTraction result;result.interface=catalog_.front().id;
		result.stamp={output_envelope_.time_s,output_envelope_.step,output_envelope_.coupling_iteration,
			output_envelope_.reference_mesh_identity_sha256,output_envelope_.layout_identity_sha256,
			output_envelope_.partition_identity_sha256,{}};
		Sha256 projection;distributed_surface_detail::AppendString(projection,"NativeTetAleFsiProjection/v1");
		for(const auto& triangle:triangles){projection.AppendLittleEndian64(triangle.triangle_id);
			for(double value:triangle.traction_on_structure_pa)projection.AppendNormalizedDouble(value);}
		result.projection_identity_sha256=projection.Hex();
		for(std::size_t row=0;row<fluid_layout_.owned_global_node_ids.size();++row){
			const std::size_t node=fluid_layout_.owned_global_node_ids[row];std::array<double,3> force{},traction{};
			for(int component=0;component<3;++component){force[component]=candidate.transfer.solid_nodal_force_n[3*node+component];
				traction[component]=force[component]/fluid_layout_.owned_reference_lumped_areas_m2[row];}
			result.consistent_nodal_force_n.push_back(force);result.traction_on_structure_pa.push_back(traction);}
		Sha256 producer;distributed_surface_detail::AppendString(producer,"NativeTetAleFsiTrial/v1");
		distributed_surface_detail::AppendString(producer,model_identity_sha256_);
		distributed_surface_detail::AppendString(producer,committed_composition_identity_sha256_);
		distributed_surface_detail::AppendString(producer,result.projection_identity_sha256);
		for(double value:candidate.solve.state)producer.AppendNormalizedDouble(value);
		result.stamp.producer_state_identity_sha256=producer.Hex();ValidateSurfaceTraction(result,fluid_layout_);return result;
	}
	std::string BuildModelIdentity()const
	{
		Sha256 hash;distributed_surface_detail::AppendString(hash,"NativeTetAleFsiModel/v1");
		distributed_surface_detail::AppendString(hash,BuildFsiCouplingEdgeIdentitySha256(edge_));
		distributed_surface_detail::AppendString(hash,fluid_layout_.layout_identity_sha256);
		for(const auto& point:reference_mesh_.points)for(double value:point)hash.AppendNormalizedDouble(value);
		for(const auto& cell:reference_mesh_.cells)for(auto node:cell.nodes)hash.AppendLittleEndian32(node);
		hash.AppendLittleEndian32(static_cast<std::uint32_t>(interface_label_));
		hash.AppendNormalizedDouble(parameters_.density);hash.AppendNormalizedDouble(parameters_.dynamic_viscosity);
		for(const auto& component:parameters_.body_acceleration_m_s2)for(double value:component)
			hash.AppendNormalizedDouble(value);
		hash.AppendLittleEndian32(pressure_gauge_node_);hash.AppendNormalizedDouble(nonlinear_tolerance_);
		hash.AppendLittleEndian64(maximum_iterations_);
		for(const auto& prescribed:prescribed_noninterface_velocity_){hash.AppendLittleEndian32(prescribed.first);
			for(double value:prescribed.second)hash.AppendNormalizedDouble(value);}
		for(int label:natural_velocity_labels_)hash.AppendLittleEndian32(static_cast<std::uint32_t>(label));
		return hash.Hex();
	}
	std::string BuildCompositionIdentity(const std::vector<double>& state,
		const std::vector<std::array<double,3>>& displacement,double time)const
	{
		Sha256 hash;distributed_surface_detail::AppendString(hash,"NativeTetAleFsiComposition/v1");
		distributed_surface_detail::AppendString(hash,model_identity_sha256_.empty()?BuildModelIdentity():model_identity_sha256_);
		hash.AppendNormalizedDouble(time);for(double value:state)hash.AppendNormalizedDouble(value);
		for(const auto& vector:displacement)for(double value:vector)hash.AppendNormalizedDouble(value);
		return hash.Hex();
	}

	std::string domain_id_,subsystem_id_;FsiCouplingEdge edge_;std::vector<DistributedSurfaceInterface> catalog_;
	DistributedSurfaceInterface structure_surface_;DistributedSurfaceLayout fluid_layout_,structure_layout_;
	NativeTetMesh reference_mesh_;NativeTaylorHoodTopology topology_;int interface_label_=-1;
	NativeNavierStokesParameters parameters_;std::uint32_t pressure_gauge_node_=0;
	double nonlinear_tolerance_=1e-8;std::size_t maximum_iterations_=12;
	std::map<std::uint32_t,std::array<double,3>> prescribed_noninterface_velocity_;
	std::set<int> natural_velocity_labels_;
	NativeTetAleKinematics ale_;
	FsiTrialLifecycle lifecycle_;std::vector<double> committed_fluid_state_;
	std::vector<std::array<double,3>> committed_interface_displacement_;
	std::string model_identity_sha256_,committed_composition_identity_sha256_;
	SurfaceFieldStamp expected_kinematics_;SurfaceFieldStampEnvelope output_envelope_;
	std::optional<SurfaceKinematics> kinematics_;std::optional<NativeTetAleFsiTrial> trial_;
	std::optional<NativeTetAleFsiPrepared> prepared_;SurfaceTraction committed_traction_;
};

} // namespace iga
#endif
