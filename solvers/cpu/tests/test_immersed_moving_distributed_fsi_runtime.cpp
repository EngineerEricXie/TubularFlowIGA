#include "MovingCheckpointRestore.hpp"
#include "CollectiveCheckpointReceipts.hpp"
#include "ImmersedMovingDistributedFsiRuntime.hpp"
#include "DistributedStrongFsiStep.hpp"
#include "CompliantChannelFsiFixture.hpp"
#include "DistributedFsiCommitCoordinator.hpp"
#define IGA_SINGLE_OWNER_MEMBRANE_TESTING
#include "SingleOwnerMembraneRuntime.hpp"
#include "MovingFsiCheckpointBundle.hpp"
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
		char paired_root[4096]={};
		Require(!PetscOptionsGetString(nullptr,nullptr,"-moving_fsi_checkpoint_root",paired_root,sizeof(paired_root),nullptr),"invalid paired checkpoint option");
		PetscBool paired_read_only=PETSC_FALSE;
		Require(!PetscOptionsGetBool(nullptr,nullptr,"-moving_fsi_checkpoint_read_only",&paired_read_only,nullptr),"invalid paired checkpoint read mode");
		PetscInt paired_source_ranks=size;
		Require(!PetscOptionsGetInt(nullptr,nullptr,"-moving_fsi_checkpoint_source_ranks",&paired_source_ranks,nullptr)
			&&paired_source_ranks>0,"invalid paired checkpoint source ranks");
		Require(!paired_read_only||paired_root[0],"paired read-only mode requires checkpoint root");
		Require(paired_read_only||paired_source_ranks==size,"paired writer source ranks must equal communicator size");
		std::unique_ptr<iga::RestoredMovingFsiPair> restored_pair;
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
			if(deforming) {
				const auto conservation=runtime.ConservationDiagnostics();
				iga::CollectiveLocalStage(comm,"strong FSI owned validation output",[&] {
					const auto& active=runtime.CommittedLayout();const auto values=Owned(runtime.CommittedState());
					std::cout<<"strong_layout "<<context.step<<' '<<active.NodeFieldRows()<<' '<<active.Rows()<<' '<<layout.global_node_count<<'\n';
					for(std::size_t local=0;local<values.size();++local) {
						const auto row=runtime.CommittedRowBegin()+static_cast<PetscInt>(local);
						const bool node=static_cast<std::size_t>(row)<active.NodeFieldRows();
						std::cout<<"strong_field "<<context.step<<' '<<row<<' '
							<<(node?active.NodeIds()[row/4]:UINT64_MAX)<<' '<<(node?row%4:4)<<' '<<values[local]<<'\n';
					}
					const auto& traction=fluid.CommittedTraction();const auto& kinematics=structure.CommittedKinematics();
					for(std::size_t local=0;local<layout.owned_global_node_ids.size();++local) {
						std::cout<<"strong_surface "<<context.step<<' '<<layout.owned_global_node_ids[local];
						for(const auto* value:{&kinematics.displacement_m[local],&kinematics.velocity_m_per_s[local],
							&traction.traction_on_structure_pa[local],&traction.consistent_nodal_force_n[local]})
							for(double component:*value)std::cout<<' '<<component;
						std::cout<<'\n';
					}
					for(std::size_t iteration=0;iteration<result.history.size();++iteration) {
						const auto& value=result.history[iteration];const auto& displacement=value.displacement;
						std::cout<<"strong_iteration "<<context.step<<' '<<iteration<<' '
							<<displacement.area_weighted_rms_residual_m<<' '<<displacement.max_residual_m<<' '
							<<displacement.displacement_scale_m<<' '<<displacement.convergence_threshold_m<<' '
							<<value.maximum_velocity_residual_times_dt_m<<' '<<value.relaxation<<'\n';
					}
					for(const auto& port:runtime.Diagnostics().ports)
						std::cout<<"strong_port "<<context.step<<' '<<port.id<<' '<<port.measurement.area_m2<<' '
							<<port.measurement.outward_flow_m3_s<<' '<<port.measurement.mean_pressure_pa<<' '
							<<port.measurement.mean_normal_traction_pa<<'\n';
					std::cout<<"strong_conservation "<<context.step;
					for(double value:{conservation.normalized_divergence_theorem_defect,conservation.normalized_reynolds_defect,
						conservation.normalized_moving_mass_defect,conservation.normalized_wall_relative_leakage,
						conservation.endpoint.normalized_discrete_moving_wall_continuity_defect,
						policy.divergence_theorem,policy.reynolds,policy.moving_mass,policy.wall_relative_leakage,policy.discrete_continuity})
						std::cout<<' '<<value;
					std::cout<<std::endl;
				});
			}
			if(paired_root[0]&&step==0) {
				const auto hash=[](const std::string& bytes) { return iga::moving_fsi_checkpoint_detail::Hash(bytes); };
				const auto configuration=hash("paired-strong-checkpoint/v1");
				const iga::CoupledCheckpointEpoch epoch{hash(paired_root),{},
					{configuration,configuration,hash(std::to_string(paired_source_ranks)),static_cast<std::uint32_t>(paired_source_ranks)},context.step,context.EndTime(),context.dt_s};
				iga::MovingCheckpointLayout source;source.domain_id="flow";
				for(PetscInt member=0;member<paired_source_ranks;++member)source.world_ranks.push_back(member);
				const auto catalog=iga::MovingFsiCheckpointCatalog(source,paired_source_ranks,paired_source_ranks-1);
				if(!paired_read_only) {
					iga::CollectiveLocalStage(comm,"paired file fixture create",[&] {
						if(rank==0) { Require(std::filesystem::create_directory(paired_root),"paired fixture directory exists");iga::CreateCoupledCheckpointEpoch(paired_root,epoch); }
					});
					const auto local=iga::WriteMovingFsiCheckpointShards(comm,paired_root,epoch,source,size-1,configuration,runtime,fluid,structure);
					const auto receipts=iga::GatherCheckpointReceipts(comm,epoch,local);
					iga::CollectiveLocalStage(comm,"paired file fixture publish",[&] {
						if(rank==0)(void)iga::PublishCoupledCheckpoint(paired_root,epoch,catalog,receipts);
					});
				}
				std::optional<iga::CoupledCheckpointManifest> manifest;
				iga::CollectiveLocalStage(comm,"paired file fixture discover",[&] {
					const auto found=iga::FindLatestCoupledCheckpoint(paired_root,epoch.compatibility,catalog);
					Require(found.latest.has_value()&&found.rejected.empty(),"paired checkpoint not found");manifest=found.latest;
				});
				const auto restore=[&](const iga::CoupledCheckpointManifest& saved) {
					return iga::RestoreMovingFsiCheckpoint(comm,paired_root,saved,source,paired_source_ranks-1,configuration,grid,geometry_options,options,map,interface,
						[&](const iga::MovingCutGeometry&,const iga::MaterialSurfaceKinematics& target) { return iga::MovingCutGeometry::Build(grid,target,geometry_options); },
						size-1,policy,deforming?iga::compliant_channel_fixture::MembraneMaterial():iga::PretensionedMembraneMaterial{1,.1,1,1});
				};
				const auto original_field=Owned(runtime.CommittedState());const auto original_publication=fluid.CaptureCheckpoint();
				const auto original_membrane=structure.CaptureCheckpoint();
				auto corrupt=*manifest;
				for(auto& shard:corrupt.shards)if(shard.spec.id=="flow.membrane.numerical")shard.sha256=hash("wrong payload hash");
				Reject(comm,[&] { auto rejected=restore(corrupt);rejected->Close(); });
				restored_pair=restore(*manifest);
				const auto restored_publication=restored_pair->fluid->CaptureCheckpoint();
				const auto restored_membrane=restored_pair->membrane->CaptureCheckpoint();
				double local_restore_error[2]={},global_restore_error[2]={};
				iga::CollectiveLocalStage(comm,"paired file fixture exact restore",[&] {
					const auto near=[](double expected,double obtained) { Require(std::isfinite(expected)&&std::isfinite(obtained)
						&&std::abs(expected-obtained)<=1e-12+1e-6*std::abs(expected),"paired accepted surface differs"); };
					const auto actual=Owned(restored_pair->flow->CommittedState());Require(actual.size()==original_field.size(),"paired restored field shape differs");
					const auto& restored_layout=restored_pair->flow->CommittedLayout();const auto& baseline_layout=runtime.CommittedLayout();
					Require(restored_layout.NodeIds()==baseline_layout.NodeIds()&&restored_layout.PortIds()==baseline_layout.PortIds()
						&&restored_layout.HasGaugeRow()==baseline_layout.HasGaugeRow(),"paired restored active row catalog differs");
					if(paired_source_ranks==size)Require(restored_layout.HashSha256()==baseline_layout.HashSha256(),"paired restored active layout differs");
					Require(original_membrane.has_value()==restored_membrane.has_value(),"paired restored membrane owner differs");
					for(std::size_t row=0;row<actual.size();++row) { const double delta=actual[row]-original_field[row];local_restore_error[0]+=delta*delta;local_restore_error[1]+=original_field[row]*original_field[row]; }
					if(!paired_read_only) {
						Require(actual==original_field,"paired restored field differs");
						Require(restored_publication==original_publication,"paired restored publication differs");
						if(original_membrane)Require(original_membrane->metadata==restored_membrane->metadata&&original_membrane->numerical==restored_membrane->numerical,"paired restored membrane differs");
					}
					const auto& expected_traction=fluid.CommittedTraction();const auto& obtained_traction=restored_pair->fluid->CommittedTraction();
					const auto& expected_kinematics=structure.CommittedKinematics();const auto& obtained_kinematics=restored_pair->membrane->CommittedKinematics();
					Require(expected_traction.traction_on_structure_pa.size()==obtained_traction.traction_on_structure_pa.size()
						&&expected_kinematics.displacement_m.size()==obtained_kinematics.displacement_m.size(),"paired accepted surface shape differs");
					for(std::size_t row=0;row<expected_traction.traction_on_structure_pa.size();++row)for(int axis=0;axis<3;++axis) {
						near(expected_traction.traction_on_structure_pa[row][axis],obtained_traction.traction_on_structure_pa[row][axis]);
						near(expected_traction.consistent_nodal_force_n[row][axis],obtained_traction.consistent_nodal_force_n[row][axis]);
						near(expected_kinematics.displacement_m[row][axis],obtained_kinematics.displacement_m[row][axis]);
						near(expected_kinematics.velocity_m_per_s[row][axis],obtained_kinematics.velocity_m_per_s[row][axis]);
					}
					const auto indices=iga::moving_fsi_checkpoint_detail::Indices(*manifest,source,paired_source_ranks-1);
					if(paired_source_ranks==size)Require(restored_publication==iga::moving_checkpoint_bundle_detail::Metadata(paired_root,*manifest,
						indices.at("fsi.rank-"+std::to_string(source.world_ranks.at(rank)))),"paired publication differs from saved file");
					if(restored_membrane) {
						Require(restored_membrane->metadata==iga::moving_checkpoint_bundle_detail::Metadata(paired_root,*manifest,indices.at("membrane.metadata")),"paired membrane metadata differs from file");
						Require(restored_membrane->numerical==iga::moving_checkpoint_bundle_detail::Metadata(paired_root,*manifest,indices.at("membrane.numerical")),"paired membrane state differs from file");
					}
					Require(Owned(runtime.CommittedState())==original_field,"paired failed candidate changed original fluid");
				});
				MPI_Allreduce(local_restore_error,global_restore_error,2,MPI_DOUBLE,MPI_SUM,comm);
				const double restore_error=std::sqrt(global_restore_error[0])/std::max(1.,std::sqrt(global_restore_error[1]));Require(restore_error<1e-8,"paired restored field differs from independent baseline");
				std::cout<<"moving_fsi_file_checkpoint rank="<<rank<<" ranks="<<size<<" source_ranks="<<paired_source_ranks
					<<" read_only="<<paired_read_only<<" exact_file_publications="<<(paired_source_ranks==size)
					<<" surface_repartition="<<(paired_source_ranks!=size)<<" accepted_field_scaled_l2="<<restore_error<<" checksum_cleanup_retry=1 passed"<<std::endl;
			}
			if(restored_pair&&step==1) {
				const auto continued=iga::SolveDistributedStrongFsiStep(comm,*restored_pair->fluid,*restored_pair->membrane,map,
					restored_pair->flow->CommittedGeometry().Evaluation(),context,iga::compliant_channel_fixture::CouplingOptions());
				const auto reference=Owned(runtime.CommittedState()),actual=Owned(restored_pair->flow->CommittedState());
				const auto reference_conservation=runtime.ConservationDiagnostics(),actual_conservation=restored_pair->flow->ConservationDiagnostics();
				double local_error[2]={},global_error[2]={};
				iga::CollectiveLocalStage(comm,"paired continuation comparison",[&] {
					const auto near=[](double expected,double obtained) {
						Require(std::isfinite(expected)&&std::isfinite(obtained)&&std::abs(expected-obtained)<=1e-12+1e-6*std::abs(expected),"paired continued diagnostic differs");
					};
					Require(actual.size()==reference.size()&&continued.history.size()==result.history.size(),"paired continuation shape or iterations differ");
					for(std::size_t row=0;row<actual.size();++row) { const double delta=actual[row]-reference[row];local_error[0]+=delta*delta;local_error[1]+=reference[row]*reference[row]; }
					const auto& expected=structure.CommittedKinematics();const auto& obtained=restored_pair->membrane->CommittedKinematics();
					const auto& expected_traction=fluid.CommittedTraction();const auto& obtained_traction=restored_pair->fluid->CommittedTraction();
					Require(expected.displacement_m.size()==obtained.displacement_m.size(),"paired continued surface shape differs");
					for(std::size_t row=0;row<expected.displacement_m.size();++row)for(int axis=0;axis<3;++axis) {
						Require(std::abs(expected.displacement_m[row][axis]-obtained.displacement_m[row][axis])<1e-12,"paired continued displacement differs");
						Require(std::abs(expected.velocity_m_per_s[row][axis]-obtained.velocity_m_per_s[row][axis])<1e-12,"paired continued velocity differs");
						near(expected_traction.traction_on_structure_pa[row][axis],obtained_traction.traction_on_structure_pa[row][axis]);
						near(expected_traction.consistent_nodal_force_n[row][axis],obtained_traction.consistent_nodal_force_n[row][axis]);
					}
					for(std::size_t iteration=0;iteration<result.history.size();++iteration) {
						const auto& a=result.history[iteration];const auto& b=continued.history[iteration];
						near(a.displacement.area_weighted_rms_residual_m,b.displacement.area_weighted_rms_residual_m);
						near(a.displacement.max_residual_m,b.displacement.max_residual_m);
						near(a.displacement.displacement_scale_m,b.displacement.displacement_scale_m);
						near(a.displacement.convergence_threshold_m,b.displacement.convergence_threshold_m);
						near(a.maximum_velocity_residual_times_dt_m,b.maximum_velocity_residual_times_dt_m);near(a.relaxation,b.relaxation);
					}
					const auto& ports=runtime.Diagnostics().ports;const auto& restored_ports=restored_pair->flow->Diagnostics().ports;
					Require(ports.size()==restored_ports.size(),"paired continued port count differs");
					for(std::size_t port=0;port<ports.size();++port) {
						Require(ports[port].id==restored_ports[port].id,"paired continued port ID differs");
						const auto& a=ports[port].measurement;const auto& b=restored_ports[port].measurement;
						near(a.area_m2,b.area_m2);near(a.outward_flow_m3_s,b.outward_flow_m3_s);
						near(a.mean_pressure_pa,b.mean_pressure_pa);near(a.mean_normal_traction_pa,b.mean_normal_traction_pa);
					}
					const auto& a=reference_conservation;const auto& b=actual_conservation;
					near(a.normalized_divergence_theorem_defect,b.normalized_divergence_theorem_defect);
					near(a.normalized_reynolds_defect,b.normalized_reynolds_defect);near(a.normalized_moving_mass_defect,b.normalized_moving_mass_defect);
					near(a.normalized_wall_relative_leakage,b.normalized_wall_relative_leakage);
					near(a.endpoint.normalized_discrete_moving_wall_continuity_defect,b.endpoint.normalized_discrete_moving_wall_continuity_defect);
				});
				MPI_Allreduce(local_error,global_error,2,MPI_DOUBLE,MPI_SUM,comm);
				const double error=std::sqrt(global_error[0])/std::max(1.,std::sqrt(global_error[1]));Require(error<1e-8,"paired continued field differs");
				std::cout<<"moving_fsi_file_continuation rank="<<rank<<" ranks="<<size<<" field_scaled_l2="<<error<<" iterations="<<continued.history.size()<<" surface_history_ports_conservation=1 passed"<<std::endl;
			}
		}
		if(restored_pair)restored_pair->Close();
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
	const auto publication_checkpoint=fluid.CaptureCheckpoint();
	const auto publication_fluid_values=Owned(runtime.CommittedState());
	const auto publication_hash=[](const std::string& bytes) { iga::Sha256 hash;hash.Append(bytes.data(),bytes.size());return hash.Hex(); };
	{
		iga::ImmersedMovingDistributedFsiRuntime fresh(runtime,map,interface,
			[&](const iga::MovingCutGeometry&,const iga::MaterialSurfaceKinematics& target) {
				return iga::MovingCutGeometry::Build(grid,target,geometry_options);
			},size-1,policy);
		auto corrupt=publication_checkpoint;if(rank==0)corrupt.back()^=1;
		Reject(comm,[&] { fresh.RestoreCheckpoint(corrupt,publication_hash(publication_checkpoint)); });
		Reject(comm,[&] { (void)fresh.CommittedTraction(); });
		auto invalid_state=iga::ParseMovingFsiPublicationCheckpoint(publication_checkpoint,publication_hash(publication_checkpoint),layout,
			interface.id,fresh.CheckpointConfigurationIdentity(),runtime.CommittedGeometry().Evaluation().ContentIdentitySha256(),
			runtime.CommittedGeometry().GeometryIdentitySha256(),runtime.Clock().index,runtime.Clock().time_s);
		invalid_state.composition_identity="invalid";
		Reject(comm,[&] { fresh.RestoreCheckpointState(invalid_state); });
		Reject(comm,[&] { (void)fresh.CommittedTraction(); });
		if(size>1) {
			// Individually valid payloads can still disagree on the shared epoch.
			// Change the saved coupling iteration on one rank and authenticate the
			// changed bytes to exercise the collective gate after local parsing.
			auto inconsistent=publication_checkpoint;
			const std::size_t iteration_offset=8+std::string("IGA_MOVING_FSI_PUBLICATION/1").size()+4*(8+64)+3*8;
			if(rank==0)inconsistent.at(iteration_offset)^=1;
			Reject(comm,[&] { fresh.RestoreCheckpoint(inconsistent,publication_hash(inconsistent)); });
			Reject(comm,[&] { (void)fresh.CommittedTraction(); });
		}
		fresh.RestoreCheckpoint(publication_checkpoint,publication_hash(publication_checkpoint));
		iga::CollectiveLocalStage(comm,"FSI restored publication exact",[&] {
			Require(iga::BuildSurfaceTractionIdentitySha256(fresh.CommittedTraction(),layout)==expected,"FSI restored traction differs");
		});
		Require(fresh.CaptureCheckpoint()==publication_checkpoint,"FSI restored publication bytes differ");
		Reject(comm,[&] { fresh.RestoreCheckpoint(publication_checkpoint,publication_hash(publication_checkpoint)); });
		Require(Owned(runtime.CommittedState())==publication_fluid_values,"FSI publication restore changed fluid");
		std::cout<<"moving_fsi_publication_checkpoint rank="<<rank<<" exact_roundtrip=1 corrupt_rejected=1 state_guard=1 overwrite_rejected=1 passed\n";
	}
	char checkpoint_root[4096]={};
	Require(!PetscOptionsGetString(nullptr,nullptr,"-moving_port_checkpoint_root",checkpoint_root,sizeof(checkpoint_root),nullptr),"invalid port checkpoint option");
	if(checkpoint_root[0]) {
		Require(deforming,"port checkpoint fixture requires deforming channel");
		const auto accepted_values=Owned(runtime.CommittedState());
		const auto original_controls=runtime.CaptureAcceptedCheckpoint().port_control_values;
		runtime.SetPortControlValue("inlet",.2);
		const auto saved=runtime.CaptureAcceptedCheckpoint();Require(saved.port_control_values.at("inlet")==.2,"queued control not captured");
		Require(runtime.CommittedGeometry().PreviousMaterialIdentitySha256().empty(),"fixture must exercise independent target publication");
		const auto hash=[](const std::string& value) { return iga::coupled_checkpoint_detail::Hash(value); };
		const std::string configuration=hash("queued-pressure-port-checkpoint/v1");
		const iga::CoupledCheckpointEpoch epoch{hash(checkpoint_root),{},
			{configuration,configuration,hash(std::to_string(size)),static_cast<std::uint32_t>(size)},1,1.,1.};
		iga::MovingCheckpointLayout source;source.domain_id="flow";
		for(int member=0;member<size;++member)source.world_ranks.push_back(member);
		iga::CollectiveLocalStage(comm,"port checkpoint create",[&] {
			if(rank==0) { Require(std::filesystem::create_directory(checkpoint_root),"port checkpoint directory exists");iga::CreateCoupledCheckpointEpoch(checkpoint_root,epoch); }
		});
		std::vector<iga::CoupledCheckpointShard> local;
		iga::CollectiveLocalStage(comm,"port checkpoint write",[&] { local=iga::WriteMovingCheckpointShards(checkpoint_root,epoch,source,rank,saved,configuration); });
		const auto receipts=iga::GatherCheckpointReceipts(comm,epoch,local);
		iga::CollectiveLocalStage(comm,"port checkpoint publish",[&] {
			if(rank==0)(void)iga::PublishCoupledCheckpoint(checkpoint_root,epoch,iga::MovingCheckpointCatalog(source,size),receipts);
		});
		std::optional<iga::CoupledCheckpointManifest> manifest;
		iga::CollectiveLocalStage(comm,"port checkpoint discover",[&] {
			const auto found=iga::FindLatestCoupledCheckpoint(checkpoint_root,epoch.compatibility,iga::MovingCheckpointCatalog(source,size));
			Require(found.latest.has_value()&&found.rejected.empty(),"port checkpoint not found");manifest=found.latest;
		});
		auto restored=iga::RestoreMovingCheckpointRuntime(comm,checkpoint_root,*manifest,source,configuration,grid,geometry_options,options);
		Require(restored->CaptureAcceptedCheckpoint().port_control_values==saved.port_control_values,"restored port controls differ");
		Require(Owned(restored->CommittedState())==accepted_values,"queued control changed saved field");
		Require(restored->CommittedGeometry().PublicationIdentitySha256()==runtime.CommittedGeometry().PublicationIdentitySha256(),"independent publication restore differs");
		const auto next=[&](const iga::MovingCutGeometry& previous) {
			std::unique_ptr<iga::MovingCutGeometry> target;
			iga::CollectiveLocalStage(comm,"port checkpoint next geometry",[&] {
				const auto& old=previous.Evaluation();std::vector<iga::RawSurfaceTriangle> triangles;
				for(const auto& item:old.SourceTriangles()) {
					iga::RawSurfaceTriangle triangle;triangle.boundary_id=item.boundary_id;
					for(int axis=0;axis<3;++axis)triangle.indices[axis]=item.source_vertex_indices[axis];
					triangles.push_back(triangle);
				}
				const auto material=iga::MaterialSurfaceKinematics::CreateFromSourceTopology(old.ReferenceMaterialVerticesM(),old.SourceVerticesM(),
					std::vector<std::array<double,3>>(old.SourceVerticesM().size(),{{0,0,0}}),triangles,2.,1.,2.);
				target=iga::MovingCutGeometry::Build(grid,material,geometry_options,&previous);
			});return target;
		};
		runtime.BeginTrial(next(runtime.CommittedGeometry()),2,2);Require(runtime.SolveTrial(),"queued-pressure reference step failed");
		restored->BeginTrial(next(restored->CommittedGeometry()),2,2);Require(restored->SolveTrial(),"queued-pressure restored step failed");
		const auto reference=Owned(runtime.TrialState()),actual=Owned(restored->TrialState());Require(reference.size()==actual.size(),"port continuation shape differs");
		double local_error[2]={},global_error[2]={};
		for(std::size_t row=0;row<actual.size();++row) { const double delta=actual[row]-reference[row];local_error[0]+=delta*delta;local_error[1]+=reference[row]*reference[row]; }
		MPI_Allreduce(local_error,global_error,2,MPI_DOUBLE,MPI_SUM,comm);
		const double error=std::sqrt(global_error[0])/std::max(1.,std::sqrt(global_error[1]));Require(std::isfinite(error)&&error<1e-8,"queued-pressure continuation differs");
		restored->Commit();restored->Close();runtime.AbortTrial();
		for(const auto& value:original_controls)runtime.SetPortControlValue(value.first,value.second);
		Require(Owned(runtime.CommittedState())==accepted_values,"port fixture changed original accepted state");
		std::cout<<"moving_port_checkpoint rank="<<rank<<" ranks="<<size<<" queued_inlet=0.2 independent_publication=1 next_step_scaled_l2="<<error<<" passed\n";
	}
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
