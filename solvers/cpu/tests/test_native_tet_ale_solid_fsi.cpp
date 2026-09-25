#include "NativeTetAleFsiRuntime.hpp"
#include "NativeTetSolidFsiRuntime.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

namespace {
template<class Function>void Reject(Function&& function)
{bool rejected=false;try{function();}catch(const std::exception&){rejected=true;}assert(rejected);}
std::string Hash(const char* text)
{iga::Sha256 hash;hash.Append(text,std::char_traits<char>::length(text));return hash.Hex();}
iga::NativeTetMesh Mesh()
{
	iga::NativeTetMesh mesh;mesh.points={{{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}},{{.25,.25,.25}}}};
	mesh.cells={iga::NativeTetCell{1,{{4,1,2,3}}},iga::NativeTetCell{2,{{0,4,2,3}}},
		iga::NativeTetCell{3,{{0,1,4,3}}},iga::NativeTetCell{4,{{0,1,2,4}}}};
	mesh.boundary_triangles={iga::NativeTetTriangle{1,{{1,2,3}},7},
		iga::NativeTetTriangle{2,{{0,3,2}},7},iga::NativeTetTriangle{3,{{0,1,3}},7},
		iga::NativeTetTriangle{4,{{0,2,1}},7}};return mesh;
}
iga::DistributedSurfaceLayout Layout()
{
	iga::DistributedSurfaceLayout value;value.reference_mesh_identity_sha256=Hash("native-ale-solid-fsi-mesh");
	value.global_node_count=4;value.partition_count=1;value.partition_rank=0;
	value.owned_global_node_ids={0,1,2,3};value.reference_positions={{0,{{0,0,0}}},{1,{{1,0,0}}},
		{2,{{0,1,0}}},{3,{{0,0,1}}}};
	value.reference_triangles={{{1,2,3}},{{0,3,2}},{{0,1,3}},{{0,2,1}}};
	const double corner=1./3.+std::sqrt(3.)/6.;
	value.owned_reference_lumped_areas_m2={.5,corner,corner,corner};
	value.layout_identity_sha256=iga::BuildDistributedSurfaceLayoutIdentitySha256(value);return value;
}
iga::DistributedSurfaceInterface Fluid(const std::string& mesh)
{
	iga::DistributedSurfaceInterface value;value.id={"fluid","native-ale","wall"};value.subsystem_id="native-ale";
	value.boundary_labels={7};value.reference_mesh_identity_sha256=mesh;
	value.provides={iga::SurfaceFieldQuantity::TractionOnStructure};
	value.requires={iga::SurfaceFieldQuantity::Displacement,iga::SurfaceFieldQuantity::Velocity};return value;
}
iga::DistributedSurfaceInterface Solid(const std::string& mesh)
{
	iga::DistributedSurfaceInterface value;value.id={"solid","native-solid","wall"};value.subsystem_id="native-solid";
	value.boundary_labels={7};value.reference_mesh_identity_sha256=mesh;
	value.provides={iga::SurfaceFieldQuantity::Displacement,iga::SurfaceFieldQuantity::Velocity};
	value.requires={iga::SurfaceFieldQuantity::TractionOnStructure};return value;
}
std::vector<double> DecayingFlowState(const iga::NativeTetMesh& mesh,double speed)
{
	const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
	const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();
	std::vector<double> state(3*velocity_nodes+mesh.points.size(),0.0);
	for(std::size_t edge=0;edge<topology.edges.size();++edge)
		if(topology.edges[edge][0]==4||topology.edges[edge][1]==4){
			const std::size_t node=mesh.points.size()+edge;
			state[3*node]=speed;
			state[3*node+1]=-0.35*speed;
		}
	return state;
}
}

int main()
{
	const auto layout=Layout();const auto fluid_surface=Fluid(layout.reference_mesh_identity_sha256);
	const auto solid_surface=Solid(layout.reference_mesh_identity_sha256);
	const iga::FsiCouplingEdge edge("native-closed-wall",fluid_surface.id,solid_surface.id,
		iga::FsiCouplingLaw::FluidStructureTractionKinematics);
	iga::NativeTetAleFsiRuntime fluid("fluid","native-ale",edge,fluid_surface,solid_surface,
		layout,layout,Mesh(),7,{1000.,.004},0,{},0.,1e-8,8);
	std::map<std::size_t,double> constraints;
	for(std::size_t node=0;node<4;++node)for(int component=0;component<3;++component)
		constraints.emplace(3*node+component,0.0);
	iga::NativeTetSolidStaticOptions solid_options;solid_options.relative_tolerance=1e-6;
	solid_options.absolute_tolerance_n=1e-9;
	iga::NativeTetSolidFsiRuntime solid("solid","native-solid",edge,solid_surface,fluid_surface,
		layout,layout,Mesh(),7,{1500.,.3,980.},constraints,solid_options);
	iga::StrongFluidStructureCouplingOptions options;options.maximum_iterations=3;
	options.absolute_displacement_tolerance_m=1e-12;options.relative_displacement_tolerance=1e-8;
	options.absolute_nodal_force_residual_tolerance_n=1e-12;
	options.absolute_traction_residual_tolerance_pa=1e-12;
	options.fluid_residual_tolerance=1e-8;
	options.structure_residual_tolerance_n=1e-9;
	options.interface_power_defect_tolerance_w=1e-12;
	iga::StrongFluidStructureCoupling<iga::NativeTetAleFsiRuntime,iga::NativeTetSolidFsiRuntime>
		coordinator(fluid,solid,edge,layout,options);
	const auto result=coordinator.Execute({1,0.,.05});
	assert(result.converged&&result.iterations==2);
	assert(!result.history[0].traction_gate_passed&&result.history[0].displacement_gate_passed);
	assert(result.history[1].traction_gate_passed&&result.history[1].nodal_force_gate_passed
		&&result.history[1].fluid_gate_passed
		&&result.history[1].structure_gate_passed&&result.history[1].interface_power_gate_passed);
	assert(fluid.CommittedAleTimeS()==.05);
	for(double value:solid.CommittedDisplacementM())assert(value==0.0);
	assert(!result.final_committed_fluid_traction_identity_sha256.empty());
	assert(!result.final_committed_structure_state_identity_sha256.empty());

	iga::NativeTetAleFsiRuntime rollback_fluid("fluid","native-ale",edge,fluid_surface,solid_surface,
		layout,layout,Mesh(),7,{1000.,.004},0,{},0.,1e-8,8);
	iga::NativeTetSolidFsiRuntime rollback_solid("solid","native-solid",edge,solid_surface,fluid_surface,
		layout,layout,Mesh(),7,{1500.,.3,980.},constraints,solid_options);
	auto rollback_options=options;rollback_options.fail_before_finalize_for_testing=true;
	iga::StrongFluidStructureCoupling<iga::NativeTetAleFsiRuntime,iga::NativeTetSolidFsiRuntime>
		rollback(rollback_fluid,rollback_solid,edge,layout,rollback_options);
	const auto solid_identity=rollback_solid.CommittedStateIdentitySha256();
	const auto fluid_identity=rollback_fluid.CommittedCompositionIdentitySha256();
	Reject([&]{(void)rollback.Execute({1,0.,.05});});
	assert(rollback_fluid.CommittedAleTimeS()==0.0);
	assert(rollback_fluid.CommittedCompositionIdentitySha256()==fluid_identity);
	assert(rollback_solid.CommittedStateIdentitySha256()==solid_identity);

	// Nonzero function-first case: an interior native P2 velocity pulse decays
	// against the closed moving wall, its computed Cauchy traction loads the
	// real hyperelastic solid, and the coupled state commits only after all
	// independent gates pass.
	const auto nonzero_mesh=Mesh();
	iga::NativeTetAleFsiRuntime loaded_fluid("fluid","native-ale",edge,fluid_surface,solid_surface,
		layout,layout,nonzero_mesh,7,{1000.,.004},0,DecayingFlowState(nonzero_mesh,1e-5),0.,1e-6,30);
	const std::map<std::size_t,double> supported{{0,0},{1,0},{2,0},{4,0},{5,0},{8,0}};
	iga::NativeTetSolidStaticOptions loaded_solid_options;loaded_solid_options.relative_tolerance=1e-5;
	loaded_solid_options.absolute_tolerance_n=1e-8;loaded_solid_options.maximum_iterations=40;
	iga::NativeTetSolidFsiRuntime loaded_solid("solid","native-solid",edge,solid_surface,fluid_surface,
		layout,layout,Mesh(),7,{1.5e6,.3,980.},supported,loaded_solid_options);
	iga::StrongFluidStructureCouplingOptions loaded_options;
	loaded_options.maximum_iterations=12;
	loaded_options.absolute_displacement_tolerance_m=2e-5;
	loaded_options.relative_displacement_tolerance=2e-2;
	loaded_options.absolute_nodal_force_residual_tolerance_n=2e-2;
	loaded_options.absolute_traction_residual_tolerance_pa=1e-1;
	loaded_options.fluid_residual_tolerance=1e-6;
	loaded_options.structure_residual_tolerance_n=1e-7;
	loaded_options.interface_power_defect_tolerance_w=1e-10;
	loaded_options.aitken_controls.initial_relaxation=.8;
	iga::StrongFluidStructureCoupling<iga::NativeTetAleFsiRuntime,iga::NativeTetSolidFsiRuntime>
		loaded(loaded_fluid,loaded_solid,edge,layout,loaded_options);
	const auto loaded_result=loaded.Execute({2,0.,.02});
	assert(loaded_result.converged&&loaded_result.iterations>=2);
	double displacement_norm=0.0;
	for(double value:loaded_solid.CommittedDisplacementM())displacement_norm+=value*value;
	assert(displacement_norm>1e-20);
	const auto& accepted_loaded=loaded_result.history.back();
	assert(accepted_loaded.displacement_gate_passed&&accepted_loaded.traction_gate_passed
		&&accepted_loaded.nodal_force_gate_passed&&accepted_loaded.fluid_gate_passed
		&&accepted_loaded.structure_gate_passed&&accepted_loaded.interface_power_gate_passed);
	std::cout<<"native ALE-solid nonzero FSI functional test passed iterations="
		<<loaded_result.iterations<<" displacement_norm_m="<<std::sqrt(displacement_norm)
		<<" displacement_residual_m="<<accepted_loaded.area_weighted_rms_residual_m
		<<" traction_residual_pa="<<*accepted_loaded.traction_residual_pa
		<<" fluid_residual="<<*accepted_loaded.fluid_acceptance_residual
		<<" solid_residual_n="<<*accepted_loaded.structure_acceptance_residual_n
		<<" power_defect_w="<<*accepted_loaded.interface_power_defect_w<<'\n';
	return 0;
}
