#include "ImmersedMovingDistributedFsiRuntime.hpp"
#include "DistributedStrongFsiStep.hpp"
#include "CompliantChannelFsiFixture.hpp"
#include "DistributedFsiCommitCoordinator.hpp"
#define IGA_SINGLE_OWNER_MEMBRANE_TESTING
#include "SingleOwnerMembraneRuntime.hpp"
#define main MovingRuntimeRegressionMain
#include "test_immersed_moving_distributed_runtime.cpp"
#undef main

namespace {
void RunPaired(MPI_Comm comm,bool deforming,bool strong,bool exhaustion)
{
	int rank=0,size=0;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
	auto soup=Cube();for(auto& vertex:soup.vertices)for(auto& value:vertex)value+=.1;
	iga::PrescribedSurfaceMotion motion({{0,soup},{1,soup}});
	const auto serial_map=[&] {
		if(!deforming)return BottomPatch(motion.Evaluate(0,0,.125));
		const auto old=iga::compliant_channel_fixture::InitialMaterial();
		std::vector<iga::RawSurfaceTriangle> triangles;
		for(const auto& source:old.SourceTriangles()) {
			iga::RawSurfaceTriangle triangle;triangle.boundary_id=source.boundary_id;
			for(int axis=0;axis<3;++axis)triangle.indices[axis]=source.source_vertex_indices[axis];
			triangles.push_back(triangle);
		}
		const auto initial=iga::MaterialSurfaceKinematics::CreateFromSourceTopology(old.ReferenceMaterialVerticesM(),
			old.SourceVerticesM(),old.SourceVertexVelocitiesMPerS(),triangles,0.,-1.,0.);
		return iga::compliant_channel_fixture::PatchMap(initial);
	}();
	auto layout=serial_map.Layout();const auto ids=layout.owned_global_node_ids;
	const auto areas=layout.owned_reference_lumped_areas_m2;
	layout.partition_rank=rank;layout.partition_count=size;
	layout.owned_global_node_ids.clear();layout.owned_reference_lumped_areas_m2.clear();
	for(std::size_t i=0;i<ids.size();++i)if((size==1?0:1+static_cast<int>(i)%(size-1))==rank) {
		layout.owned_global_node_ids.push_back(ids[i]);layout.owned_reference_lumped_areas_m2.push_back(areas[i]);
	}
	layout.layout_identity_sha256=iga::BuildDistributedSurfaceLayoutIdentitySha256(layout);
	std::vector<iga::MaterialSurfacePatchMap::GlobalToSourceVertex> vertices;
	for(std::size_t i=0;i<ids.size();++i)vertices.emplace_back(ids[i],i);
	std::vector<std::uint32_t> triangles;for(std::size_t i=0;i<layout.reference_triangles.size();++i)triangles.push_back(i);
	const auto& clamps=serial_map.ConfiguredClampedGlobalNodeIds();
	const auto map=iga::MaterialSurfacePatchMap::Create(serial_map.Interface(),layout,serial_map.FullReference(),
		deforming?7:8,vertices,triangles,clamps);
	auto interface=map.Interface();interface.id={"fluid","flow","bottom"};
	interface.subsystem_id="flow";
	interface.provides={iga::SurfaceFieldQuantity::TractionOnStructure};
	interface.requires={iga::SurfaceFieldQuantity::Displacement,iga::SurfaceFieldQuantity::Velocity};
	iga::CubicCartesianGridSpec grid{{{0,0,0}},{{1.2,1.2,1.2}},{{4,3,3}}};
	iga::MovingCutGeometryOptions geometry_options;geometry_options.volume.max_depth=2;
	geometry_options.volume_storage=iga::CutCellVolumeQuadratureStorageMode::Compact;
	const auto channel=iga::compliant_channel_fixture::FlowOptions();
	if(deforming) { grid=channel.grid;geometry_options=channel.geometry; }
	std::unique_ptr<iga::MovingCutGeometry> initial;
	iga::CollectiveLocalStage(comm,"paired initial geometry",[&] {
		initial=iga::MovingCutGeometry::Build(grid,map.FullReference(),geometry_options);
	});
	iga::ImmersedTransientFlowOptions options;options.parameters={1,.1,0};options.wall_labels={7,8,9};
	if(deforming)options=channel.flow;
	if(strong&&deforming) {
		// Fixed dimensional scale from the channel data, not from a measured
		// residual: Q_ref = width * height^3 * delta_p / (12 * mu * length).
		// This parallel-plate upper bound neglects the square duct sidewalls.
		const double width=.6,height=.6,length=.6,pressure_drop=.1;
		options.flow_controller_reference_flow_m3_s=width*height*height*height*pressure_drop
			/(12.*options.parameters.dynamic_viscosity*length);
	}
	iga::ImmersedMovingTransientDistributedRuntime runtime(comm,std::move(initial),options);
	iga::ImmersedMovingDistributedFsiRuntime::Policy policy{1e-8,1e-8,1e-8,1e-8,1e-8};
	if(deforming)policy={.03,1e-8,.03,.03,1e-8};
	iga::ImmersedMovingDistributedFsiRuntime fluid(runtime,map,interface,
		[&](const iga::MovingCutGeometry&,const iga::MaterialSurfaceKinematics& target) {
			return iga::MovingCutGeometry::Build(grid,target,geometry_options);
		},size-1,policy);
	iga::SingleOwnerMembraneRuntime structure(comm,layout,size-1,map.Interface(),interface,
		deforming?iga::compliant_channel_fixture::MembraneMaterial():iga::PretensionedMembraneMaterial{1,.1,1,1},clamps);
	const auto original=Owned(runtime.CommittedState());
	const auto material=runtime.CommittedGeometry().Evaluation().ContentIdentitySha256();
	iga::FsiTrialContext context;context.step=1;context.dt_s=deforming?1.:.125;context.coupling_iteration=2;
	iga::SurfaceKinematics input;input.interface=map.Interface().id;
	input.stamp.time_s=context.EndTime();input.stamp.step=context.step;input.stamp.coupling_iteration=context.coupling_iteration;
	input.stamp.reference_mesh_identity_sha256=layout.reference_mesh_identity_sha256;
	input.stamp.layout_identity_sha256=layout.layout_identity_sha256;
	input.stamp.partition_identity_sha256=iga::BuildDistributedSurfacePartitionIdentitySha256(layout);
	iga::Sha256 hash;hash.AppendLittleEndian64(rank);input.stamp.producer_state_identity_sha256=hash.Hex();
	input.displacement_m.assign(layout.owned_global_node_ids.size(),{{0,0,0}});
	input.velocity_m_per_s=input.displacement_m;
	if(deforming)for(std::size_t i=0;i<layout.owned_global_node_ids.size();++i)if(layout.owned_global_node_ids[i]==15) {
		input.displacement_m[i][2]=.006;input.velocity_m_per_s[i][2]=.006;
	}
	if(strong) {
		if(exhaustion) {
			context.coupling_iteration=0;
			auto controls=iga::compliant_channel_fixture::CouplingOptions();controls.maximum_iterations=1;
			std::string failure;
			try { (void)iga::SolveDistributedStrongFsiStep(comm,fluid,structure,map,
				runtime.CommittedGeometry().Evaluation(),context,controls); }
			catch(const std::exception& error) { failure=error.what(); }
			iga::CollectiveLocalStage(comm,"strong FSI iteration exhaustion",[&] {
				Require(failure=="distributed strong FSI iteration did not converge","expected iteration exhaustion, not an unrelated failure");
				Require(runtime.Clock().index==0 && runtime.Clock().time_s==0 && !runtime.Clock().trial_active,"iteration exhaustion changed accepted clock");
				Require(Owned(runtime.CommittedState())==original,"iteration exhaustion changed accepted field");
				Require(runtime.CommittedGeometry().Evaluation().ContentIdentitySha256()==material,"iteration exhaustion changed material");
			});
			Reject(comm,[&] { (void)fluid.TrialTraction(); });
			Reject(comm,[&] { (void)structure.TrialKinematics(); });
			Reject(comm,[&] { (void)fluid.CommittedTraction(); });
			Reject(comm,[&] { (void)structure.CommittedKinematics(); });
			std::cout<<"strong_fsi_exhaustion rank="<<rank<<" exact_rollback=1 passed"<<std::endl;
		}
		if(!deforming) {
			context.coupling_iteration=0;
			auto controls=iga::compliant_channel_fixture::CouplingOptions();
			const auto check_abort=[&] {
				iga::CollectiveLocalStage(comm,"strong FSI failure isolation",[&] {
					Require(runtime.Clock().index==0 && runtime.Clock().time_s==0 && !runtime.Clock().trial_active,"strong failure changed accepted clock");
					Require(Owned(runtime.CommittedState())==original,"strong failure changed accepted field");
					Require(runtime.CommittedGeometry().Evaluation().ContentIdentitySha256()==material,"strong failure changed accepted material");
				});
				Reject(comm,[&] { (void)fluid.TrialTraction(); });
				Reject(comm,[&] { (void)structure.TrialKinematics(); });
				Reject(comm,[&] { (void)fluid.CommittedTraction(); });
				Reject(comm,[&] { (void)structure.CommittedKinematics(); });
			};
			if(size>1) {
				Reject(comm,[&] { (void)iga::SolveDistributedStrongFsiStep(MPI_COMM_SELF,fluid,structure,map,
					runtime.CommittedGeometry().Evaluation(),context,controls); });
				check_abort();
			}
			controls.fail_before_finalize_for_testing=true;
			Reject(comm,[&] { (void)iga::SolveDistributedStrongFsiStep(comm,fluid,structure,map,
				runtime.CommittedGeometry().Evaluation(),context,controls); });
			check_abort();controls.fail_before_finalize_for_testing=false;
			structure.SetFailureForTesting(iga::SingleOwnerMembraneRuntime::FailurePoint::AfterPrepare,0);
			Reject(comm,[&] { (void)iga::SolveDistributedStrongFsiStep(comm,fluid,structure,map,
				runtime.CommittedGeometry().Evaluation(),context,controls); });
			check_abort();structure.SetFailureForTesting(iga::SingleOwnerMembraneRuntime::FailurePoint::None,-1);
			if(rank==0)controls.maximum_iterations=0;
			Reject(comm,[&] { (void)iga::SolveDistributedStrongFsiStep(comm,fluid,structure,map,
				runtime.CommittedGeometry().Evaluation(),context,controls); });
			check_abort();
			std::cout<<"strong_fsi_abort rank="<<rank<<" communicator_guard="<<(size>1)<<" precommit_failure=1 structure_prepare_failure=1 invalid_controls=1 exact_rollback=1 passed"<<std::endl;
		}
		for(int step=0;step<2;++step) {
			context.step=runtime.Clock().index+1;context.start_time_s=runtime.Clock().time_s;context.coupling_iteration=0;
			const auto result=iga::SolveDistributedStrongFsiStep(comm,fluid,structure,map,
				runtime.CommittedGeometry().Evaluation(),context,iga::compliant_channel_fixture::CouplingOptions());
			iga::CollectiveLocalStage(comm,"strong FSI accepted step",[&] {
				Require(result.converged&&!result.history.empty(),"strong FSI step did not converge");
				Require(runtime.Clock().index==context.step && runtime.Clock().time_s==context.EndTime(),"strong FSI accepted clock mismatch");
				Require(fluid.CommittedTraction().stamp.step==context.step && structure.CommittedKinematics().stamp.step==context.step,"strong FSI publication step mismatch");
			});
			const auto& last=result.history.back();
			std::cout<<std::setprecision(17)<<"strong_fsi_step rank="<<rank<<" ranks="<<size<<" step="<<context.step
				<<" iterations="<<result.history.size()<<" rms="<<last.displacement.area_weighted_rms_residual_m
				<<" velocity_residual_times_dt="<<last.maximum_velocity_residual_times_dt_m
				<<" threshold="<<last.displacement.convergence_threshold_m<<" passed"<<std::endl;
		}
		runtime.Close();return;
	}
	const auto solve=[&] {
		fluid.SolveTrial(context,input,input.stamp,material);
		const auto& traction=fluid.TrialTraction();
		structure.SolveTrial(context,traction,traction.stamp,traction.projection_identity_sha256);
		if(deforming) {
			double local_force=0,local_displacement=0;
			for(const auto& force:traction.consistent_nodal_force_n)for(double value:force)local_force+=value*value;
			for(const auto& displacement:structure.TrialKinematics().displacement_m)
				for(double value:displacement)local_displacement+=value*value;
			double force=0,displacement=0;MPI_Allreduce(&local_force,&force,1,MPI_DOUBLE,MPI_SUM,comm);
			MPI_Allreduce(&local_displacement,&displacement,1,MPI_DOUBLE,MPI_SUM,comm);
			iga::CollectiveLocalStage(comm,"deforming paired nontrivial response",[&] {
				Require(std::isfinite(force)&&force>1e-12,"deforming patch force is trivial");
				Require(std::isfinite(displacement)&&displacement>1e-16,"membrane response is trivial");
				Require(runtime.TrialGeometry().Evaluation().SourceVerticesM()[4][2]>
					runtime.CommittedGeometry().Evaluation().SourceVerticesM()[4][2],"fluid material patch did not deform");
			});
		}

	};
	const auto unchanged=[&] {
		iga::CollectiveLocalStage(comm,"paired accepted state isolation",[&] {
			Require(runtime.Clock().index==0 && runtime.Clock().time_s==0 && !runtime.Clock().trial_active,"paired abort changed clock");
			Require(Owned(runtime.CommittedState())==original,"paired abort changed owned field");
			Require(runtime.CommittedGeometry().Evaluation().ContentIdentitySha256()==material,"paired abort changed material");
		});
		Reject(comm,[&] { (void)fluid.TrialTraction(); });
		Reject(comm,[&] { (void)structure.TrialKinematics(); });
		Reject(comm,[&] { (void)fluid.CommittedTraction(); });
		Reject(comm,[&] { (void)structure.CommittedKinematics(); });
	};
	solve();
	structure.SetFailureForTesting(iga::SingleOwnerMembraneRuntime::FailurePoint::AfterPrepare,0);
	Reject(comm,[&] { iga::DistributedFsiCommitCoordinator::Commit(comm,context,fluid,structure); });
	unchanged();structure.SetFailureForTesting(iga::SingleOwnerMembraneRuntime::FailurePoint::None,-1);
	solve();auto wrong=context;if(rank==0)++wrong.coupling_iteration;
	Reject(comm,[&] { iga::DistributedFsiCommitCoordinator::Commit(comm,wrong,fluid,structure); });
	unchanged();solve();
	const auto expected=iga::BuildSurfaceTractionIdentitySha256(fluid.TrialTraction(),layout);
	const auto expected_structure=iga::BuildSurfaceKinematicsIdentitySha256(structure.TrialKinematics(),layout);
	iga::DistributedFsiCommitCoordinator::Commit(comm,context,fluid,structure);
	iga::CollectiveLocalStage(comm,"paired accepted publications",[&] {
		Require(runtime.Clock().index==1 && runtime.Clock().time_s==context.EndTime() && !runtime.Clock().trial_active,"paired commit clock mismatch");
		Require(iga::BuildSurfaceTractionIdentitySha256(fluid.CommittedTraction(),layout)==expected,"paired fluid publication mismatch");
		Require(iga::BuildSurfaceKinematicsIdentitySha256(structure.CommittedKinematics(),layout)==expected_structure,"paired structure publication mismatch");
	});
	if(!deforming)for(int step=0;step<3;++step) {
		context.start_time_s=runtime.Clock().time_s;context.dt_s=.1;++context.step;
		input.stamp.step=context.step;input.stamp.time_s=context.EndTime();
		const auto source=runtime.CommittedGeometry().Evaluation().ContentIdentitySha256();
		fluid.SolveTrial(context,input,input.stamp,source);
		const auto& traction=fluid.TrialTraction();
		structure.SolveTrial(context,traction,traction.stamp,traction.projection_identity_sha256);
		iga::DistributedFsiCommitCoordinator::Commit(comm,context,fluid,structure);
		iga::CollectiveLocalStage(comm,"paired decimal step clock",[&] {
			Require(runtime.Clock().index==context.step && runtime.Clock().time_s==context.EndTime(),"decimal paired clock mismatch");
			Require(fluid.CommittedTraction().stamp.time_s==context.EndTime()
				&& structure.CommittedKinematics().stamp.time_s==context.EndTime(),"decimal publication time mismatch");
		});
	}
	// Test-only owned output: the comparison tool can reconstruct the small
	// fixture without adding a production gather or publication path.
	if(deforming) {
		std::cout<<std::setprecision(17);
		std::cout<<"paired_layout "<<runtime.CommittedLayout().NodeFieldRows()<<' '
			<<runtime.CommittedLayout().Rows()<<' '<<layout.global_node_count<<'\n';
		const auto values=Owned(runtime.CommittedState());
		for(std::size_t row=0;row<values.size();++row)
			std::cout<<"paired_field "<<runtime.CommittedRowBegin()+static_cast<PetscInt>(row)<<' '<<values[row]<<'\n';
		const auto& traction=fluid.CommittedTraction();const auto& kinematics=structure.CommittedKinematics();
		for(std::size_t row=0;row<layout.owned_global_node_ids.size();++row) {
			std::cout<<"paired_surface "<<layout.owned_global_node_ids[row];
			for(const auto* value:{&kinematics.displacement_m[row],&kinematics.velocity_m_per_s[row],
				&traction.traction_on_structure_pa[row],&traction.consistent_nodal_force_n[row]})
				for(double component:*value)std::cout<<' '<<component;
			std::cout<<'\n';
		}
	}
	runtime.Close();
	std::cout<<"moving_fsi_paired rank="<<rank<<" ranks="<<size<<" owned="<<layout.owned_global_node_ids.size()
		<<" deforming="<<deforming<<" prepare_failure=1 context_failure=1 exact_rollback=1 paired_commit=1 decimal_steps="<<(!deforming?3:0)<<" transaction_fixture=1 passed\n";
}
}
int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);int status=0;
	try {
		const std::string mode=argc>1?argv[1]:"";
		const bool exhaustion=mode=="strong-exhaustion";
		RunPaired(PETSC_COMM_WORLD,exhaustion||mode=="strong"||mode=="deforming",
			exhaustion||mode=="strong"||mode=="strong-zero",exhaustion);
	} catch(const std::exception& error) { std::cerr<<error.what()<<'\n';status=1; }
	int global=0;MPI_Allreduce(&status,&global,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);PetscFinalize();return global;
}
