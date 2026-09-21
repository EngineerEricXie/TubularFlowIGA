#include "NativeTetSolidFsiRuntime.hpp"

#include <cassert>
#include <map>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>

namespace {
template<class Function> void Reject(Function&& function)
{ bool rejected=false;try{function();}catch(const std::exception&){rejected=true;}assert(rejected); }
std::string Hash(const char* text)
{ iga::Sha256 hash;hash.Append(text,std::char_traits<char>::length(text));return hash.Hex(); }
iga::NativeTetMesh Mesh()
{
	iga::NativeTetMesh mesh;mesh.points={{{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}}}};
	mesh.cells={iga::NativeTetCell{1,{{0,1,2,3}}}};
	mesh.boundary_triangles={iga::NativeTetTriangle{10,{{0,1,2}},7},
		iga::NativeTetTriangle{11,{{0,3,1}},0},iga::NativeTetTriangle{12,{{0,2,3}},0},
		iga::NativeTetTriangle{13,{{1,3,2}},0}};return mesh;
}
iga::DistributedSurfaceLayout Layout()
{
	iga::DistributedSurfaceLayout value;value.reference_mesh_identity_sha256=Hash("native-tet-fsi-mesh");
	value.global_node_count=3;value.partition_count=1;value.partition_rank=0;
	value.owned_global_node_ids={0,1,2};value.reference_positions={{0,{{0,0,0}}},{1,{{1,0,0}}},{2,{{0,1,0}}}};
	value.reference_triangles={{{0,1,2}}};value.owned_reference_lumped_areas_m2={1./6.,1./6.,1./6.};
	value.layout_identity_sha256=iga::BuildDistributedSurfaceLayoutIdentitySha256(value);return value;
}
iga::DistributedSurfaceInterface Structure(const std::string& mesh)
{
	iga::DistributedSurfaceInterface value;value.id={"solid","native-tet-solid","wall"};value.subsystem_id="native-tet-solid";
	value.boundary_labels={7};value.reference_mesh_identity_sha256=mesh;
	value.provides={iga::SurfaceFieldQuantity::Displacement,iga::SurfaceFieldQuantity::Velocity};
	value.requires={iga::SurfaceFieldQuantity::TractionOnStructure};return value;
}
iga::DistributedSurfaceInterface Fluid(const std::string& mesh)
{
	iga::DistributedSurfaceInterface value;value.id={"fluid","native-tet-ale","wall"};value.subsystem_id="native-tet-ale";
	value.boundary_labels={7};value.reference_mesh_identity_sha256=mesh;
	value.provides={iga::SurfaceFieldQuantity::TractionOnStructure};
	value.requires={iga::SurfaceFieldQuantity::Displacement,iga::SurfaceFieldQuantity::Velocity};return value;
}
iga::SurfaceFieldStamp Stamp(const iga::DistributedSurfaceLayout& layout,double time,
	std::uint64_t step,std::uint64_t iteration,const char* state)
{return {time,step,iteration,layout.reference_mesh_identity_sha256,layout.layout_identity_sha256,
	iga::BuildDistributedSurfacePartitionIdentitySha256(layout),Hash(state)};}
iga::SurfaceFieldStampEnvelope Envelope(const iga::DistributedSurfaceLayout& layout,double time,
	std::uint64_t step,std::uint64_t iteration)
{return {time,step,iteration,layout.reference_mesh_identity_sha256,layout.layout_identity_sha256,
	iga::BuildDistributedSurfacePartitionIdentitySha256(layout)};}
iga::SurfaceTraction Traction(const iga::DistributedSurfaceLayout& layout,
	const iga::SurfaceInterfaceRef& fluid,const iga::SurfaceFieldStamp& stamp,double force)
{
	iga::SurfaceTraction value;value.interface=fluid;value.stamp=stamp;
	value.traction_on_structure_pa.assign(layout.owned_global_node_ids.size(),{{0,0,0}});
	value.consistent_nodal_force_n.assign(layout.owned_global_node_ids.size(),{{0,0,0}});
	value.consistent_nodal_force_n[2][2]=force;value.projection_identity_sha256=Hash("native-projection");return value;
}

class ConstantLoadFluid {
public:
	ConstantLoadFluid(iga::FsiCouplingEdge edge,iga::DistributedSurfaceLayout layout,double force)
		:edge_(std::move(edge)),layout_(std::move(layout)),force_(force){}
	const iga::FsiCouplingEdge& StrongCouplingEdge()const noexcept{return edge_;}
	const iga::DistributedSurfaceLayout& StrongCouplingLayout()const noexcept{return layout_;}
	void BeginMacroStep(const iga::DomainStepContext& step){if(active_)throw std::runtime_error("fluid step active");step_=step;active_=true;}
	void BeginCouplingIteration(std::uint64_t iteration,const iga::SurfaceFieldStamp& expected,
		const iga::SurfaceFieldStampEnvelope& envelope)
	{
		if(!active_||trial_)throw std::runtime_error("fluid iteration state invalid");
		if(iteration!=envelope.coupling_iteration)throw std::runtime_error("fluid iteration mismatch");
		expected_=expected;envelope_=envelope;has_input_=false;
	}
	void SetSurfaceKinematics(const std::string& interface_id,const iga::SurfaceKinematics& value)
	{
		if(interface_id!=edge_.fluid.interface_id)throw std::runtime_error("unknown fluid interface");
		iga::ValidateFsiFluidKinematicsInput(edge_,value,layout_,expected_);input_=value;has_input_=true;
	}
	void SolveFluidTrial()
	{
		if(!has_input_)throw std::runtime_error("fluid input missing");
		iga::SurfaceTraction value;value.interface=edge_.fluid;
		value.stamp={envelope_.time_s,envelope_.step,envelope_.coupling_iteration,
			envelope_.reference_mesh_identity_sha256,envelope_.layout_identity_sha256,
			envelope_.partition_identity_sha256,Hash("constant-load-fluid-state")};
		value.traction_on_structure_pa.assign(layout_.owned_global_node_ids.size(),{{0,0,0}});
		value.consistent_nodal_force_n.assign(layout_.owned_global_node_ids.size(),{{0,0,0}});
		value.consistent_nodal_force_n[2][2]=force_;value.projection_identity_sha256=Hash("constant-load-projection");
		trial_=std::move(value);
	}
	iga::SurfaceTraction GetSurfaceTraction(const std::string& interface_id)const
	{if(interface_id!=edge_.fluid.interface_id||!trial_)throw std::runtime_error("fluid trial unavailable");return *trial_;}
	void RejectCouplingIteration(){if(!active_||!trial_)throw std::runtime_error("fluid reject invalid");trial_.reset();has_input_=false;}
	void PrepareCommitStep(){if(!trial_)throw std::runtime_error("fluid prepare invalid");prepared_=*trial_;}
	std::string CoordinatorPreparedCommittedTractionIdentitySha256()const
	{if(!prepared_)throw std::runtime_error("fluid prepared traction unavailable");return iga::BuildSurfaceTractionIdentitySha256(*prepared_,layout_);}
	std::string CoordinatorPreparedCommittedCompositionIdentitySha256()const
	{if(!prepared_)throw std::runtime_error("fluid prepared state unavailable");return Hash("constant-load-composition");}
	void AbortStep()noexcept{trial_.reset();prepared_.reset();has_input_=false;active_=false;}
private:
	friend class iga::StrongFluidStructureCouplingAccess;
	void CoordinatorRequireFinalizeAllowed()const{if(!active_||!trial_||!prepared_)throw std::runtime_error("fluid finalize invalid");}
	void CoordinatorFinalizeCommitNoexcept()noexcept{committed_=std::move(prepared_);trial_.reset();has_input_=false;active_=false;}
	iga::FsiCouplingEdge edge_;iga::DistributedSurfaceLayout layout_;double force_=0.;
	iga::DomainStepContext step_;iga::SurfaceFieldStamp expected_;iga::SurfaceFieldStampEnvelope envelope_;
	iga::SurfaceKinematics input_;std::optional<iga::SurfaceTraction> trial_,prepared_,committed_;
	bool active_=false,has_input_=false;
};
}

static_assert(!std::is_copy_constructible<iga::NativeTetSolidFsiRuntime>::value,"runtime must have one owner");
static_assert(!std::is_move_constructible<iga::NativeTetSolidFsiRuntime>::value,"runtime capabilities must not move");

int main()
{
	const auto layout=Layout();const auto structure=Structure(layout.reference_mesh_identity_sha256);
	const auto fluid=Fluid(layout.reference_mesh_identity_sha256);
	const iga::FsiCouplingEdge edge("native-wall",fluid.id,structure.id,
		iga::FsiCouplingLaw::FluidStructureTractionKinematics);
	// Nodes 0,1,3 are fixed. Node 2 may move only in z: one stable nonlinear DOF.
	const std::map<std::size_t,double> constraints{{0,0},{1,0},{2,0},{3,0},{4,0},{5,0},
		{6,0},{7,0},{9,0},{10,0},{11,0}};
	iga::NativeTetSolidStaticOptions options;options.relative_tolerance=1e-7;options.absolute_tolerance_n=1e-9;
	iga::NativeTetSolidFsiRuntime runtime("solid","native-tet-solid",edge,structure,fluid,
		layout,layout,Mesh(),7,{1500.,.3,980.},constraints,options);
	const iga::DomainStepContext step{2,.1,.05};runtime.BeginMacroStep(step);
	const auto stamp0=Stamp(layout,.1+.05,2,0,"fluid-0");const auto envelope0=Envelope(layout,.1+.05,2,0);
	auto traction0=Traction(layout,fluid.id,stamp0,.02);auto stale=traction0;stale.stamp.coupling_iteration=4;
	runtime.BeginCouplingIteration(0,stamp0,envelope0);Reject([&]{runtime.SetSurfaceTraction("wall",stale);});
	runtime.SetSurfaceTraction("wall",traction0);runtime.SolveSolidTrial();
	const auto rejected_trial=runtime.GetSurfaceKinematics("wall");assert(rejected_trial.displacement_m[2][2]>0.0);
	const auto committed_before=runtime.CommittedStateIdentitySha256();runtime.RejectCouplingIteration();
	assert(runtime.CommittedStateIdentitySha256()==committed_before);

	const auto stamp1=Stamp(layout,.1+.05,2,1,"fluid-1");const auto envelope1=Envelope(layout,.1+.05,2,1);
	runtime.BeginCouplingIteration(1,stamp1,envelope1);
	runtime.SetSurfaceTraction("wall",Traction(layout,fluid.id,stamp1,.03));runtime.SolveMembraneTrial();
	const auto prepared_stamp=runtime.GetSurfaceKinematics("wall").stamp;runtime.PrepareCommitStep();runtime.AbortStep();
	assert(runtime.CommittedStateIdentitySha256()==committed_before);
	Reject([&]{runtime.GetCommittedSurfaceKinematics("wall",prepared_stamp);});

	runtime.BeginStep(step);const auto stamp2=Stamp(layout,.1+.05,2,2,"fluid-2");const auto envelope2=Envelope(layout,.1+.05,2,2);
	runtime.BeginIteration(2,stamp2,envelope2);
	runtime.SetSurfaceTraction("wall",Traction(layout,fluid.id,stamp2,.04));runtime.SolveTrial();
	const auto accepted=runtime.GetSurfaceKinematics("wall");runtime.PrepareCommit();
	assert(runtime.CoordinatorPreparedCommittedKinematicsIdentitySha256()==
		iga::BuildSurfaceKinematicsIdentitySha256(accepted,layout));
	runtime.FinalizeCommit();assert(runtime.CommittedStateIdentitySha256()!=committed_before);
	const auto committed=runtime.GetCommittedSurfaceKinematics("wall",accepted.stamp);
	assert(committed.displacement_m==accepted.displacement_m);
	assert(runtime.CommittedVelocityMPerS()[8]==runtime.CommittedDisplacementM()[8]/step.dt_s);
	Reject([&]{runtime.GetCommittedSurfaceKinematics("wall",rejected_trial.stamp);});

	// Instantiate the existing generic coordinator with the real native solid
	// adapter. Unit relaxation makes the constant-load fixed point converge on
	// the second iteration, proving reject/Aitken/retry/paired commit wiring.
	iga::NativeTetSolidFsiRuntime coupled_solid("solid","native-tet-solid",edge,structure,fluid,
		layout,layout,Mesh(),7,{1500.,.3,980.},constraints,options);
	ConstantLoadFluid coupled_fluid(edge,layout,.04);
	iga::StrongFluidStructureCouplingOptions coupling_options;coupling_options.maximum_iterations=4;
	coupling_options.absolute_displacement_tolerance_m=1e-12;
	coupling_options.relative_displacement_tolerance=1e-10;
	coupling_options.aitken_controls.initial_relaxation=1.0;
	coupling_options.aitken_controls.minimum_relaxation=1.0;
	coupling_options.aitken_controls.maximum_relaxation=1.0;
	iga::StrongFluidStructureCoupling<ConstantLoadFluid,iga::NativeTetSolidFsiRuntime> coordinator(
		coupled_fluid,coupled_solid,edge,layout,coupling_options);
	const auto coupled_result=coordinator.Execute({3,.15,.05});
	assert(coupled_result.converged&&coupled_result.iterations==2);
	assert(coupled_result.history.front().aitken_proposal_applied);
	assert(coupled_solid.CommittedDisplacementM()[8]>0.0);

	iga::NativeTetSolidFsiRuntime rollback_solid("solid","native-tet-solid",edge,structure,fluid,
		layout,layout,Mesh(),7,{1500.,.3,980.},constraints,options);
	ConstantLoadFluid rollback_fluid(edge,layout,.04);
	auto rollback_options=coupling_options;rollback_options.fail_before_finalize_for_testing=true;
	iga::StrongFluidStructureCoupling<ConstantLoadFluid,iga::NativeTetSolidFsiRuntime> rollback_coordinator(
		rollback_fluid,rollback_solid,edge,layout,rollback_options);
	const auto rollback_identity=rollback_solid.CommittedStateIdentitySha256();
	Reject([&]{(void)rollback_coordinator.Execute({4,.2,.05});});
	assert(rollback_solid.CommittedStateIdentitySha256()==rollback_identity);
	return 0;
}
