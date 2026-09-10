#include "ThreeDImmersedMovingDistributedFlowDomain.hpp"
#include "MovingImmersedTransientFlowRuntime.hpp"
#include "PrescribedSurfaceMotion.hpp"
#include <cstring>
#include <iostream>

namespace {
iga::RawSurfaceSoup Soup()
{
	iga::RawSurfaceSoup soup;
	soup.vertices = {{{{0,0,0}},{{1,0,0}},{{1,1,0}},{{0,1,0}},{{0,0,1}},{{1,0,1}},{{1,1,1}},{{0,1,1}}}};
	const std::array<std::array<std::int64_t,3>,12> triangles{{{{0,2,1}},{{0,3,2}},{{4,5,6}},{{4,6,7}},
		{{0,1,5}},{{0,5,4}},{{1,2,6}},{{1,6,5}},{{2,3,7}},{{2,7,6}},{{3,0,4}},{{3,4,7}}}};
	for (std::size_t i = 0; i < triangles.size(); ++i) {
		iga::RawSurfaceTriangle triangle; triangle.indices = triangles[i]; triangle.boundary_id = i < 2 ? 8 : i < 4 ? 9 : 7;
		soup.triangles.push_back(triangle);
	}
	return soup;
}
void Require(bool value,const char* message)
{
	if(!value)throw std::runtime_error(message);
}
template<class Function> void Reject(MPI_Comm comm,Function&& function)
{
	bool rejected=false;try { function(); } catch(const std::exception&) { rejected=true; }
	iga::CollectiveLocalStage(comm,"moving graph expected rejection",[&] { Require(rejected,"missing moving graph rejection"); });
}
iga::PortBoundaryData Input(double time,double flow)
{
	iga::PortBoundaryData result;result.time_s=time;result.outward_flow_m3_s=flow;return result;
}
void Run(MPI_Comm comm,bool reject_mass)
{
	int rank=0,size=0;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
	iga::RawSurfaceSoup source=Soup();for(auto& x:source.vertices)for(double& value:x)value+=.1;
	auto target=source;for(auto& x:target.vertices)x[0]+=.74;
	iga::PrescribedSurfaceMotion motion({{0,source},{.25,target}});
	iga::DomainStepContext decimal_step;decimal_step.start_time_s=.1;decimal_step.dt_s=.02;
	const auto decimal_geometry=motion.Evaluate(decimal_step.EndTime(),decimal_step.start_time_s,decimal_step.EndTime());
	Require(decimal_geometry.DtS()!=decimal_step.dt_s,"decimal interval fixture does not exercise subtraction roundoff");
	iga::ValidateImmersedMovingGraphInterval(decimal_geometry,decimal_step);
	auto wrong_decimal=decimal_step;wrong_decimal.start_time_s+=.001;
	Reject(comm,[&] { iga::ValidateImmersedMovingGraphInterval(decimal_geometry,wrong_decimal); });
	auto shrunk=target;for(auto& x:shrunk.vertices)x[0]-=.06*(x[0]-.84);
	iga::PrescribedSurfaceMotion contraction({{0,source},{.25,shrunk}});
	const iga::CubicCartesianGridSpec grid{{{0,0,0}},{{2.4,1.2,1.2}},{{8,3,3}}};
	iga::MovingCutGeometryOptions geometry_options;geometry_options.volume.max_depth=2;
	geometry_options.volume_storage=iga::CutCellVolumeQuadratureStorageMode::Compact;
	geometry_options.volume_fitting.emplace();
	geometry_options.volume_fitting->support_expansion=3.;
	bool fail_geometry=false,fail_force=false,bad_time=false;
	iga::ImmersedTransientFlowOptions options;options.parameters={1,.1,0};options.wall_labels={7};
	options.flow_controller_reference_flow_m3_s=1.;
	options.ports={{"inlet",8,iga::ImmersedFlowPortControlMode::FlowRate,0},{"outlet",9,iga::ImmersedFlowPortControlMode::FlowRate,0}};
	options.body_force=[&](const std::array<double,3>&) {
		if(fail_force&&rank==0)throw std::runtime_error("injected graph force failure");
		return std::array<double,3>{};
	};
	std::unique_ptr<iga::MovingCutGeometry> initial;
	iga::CollectiveLocalStage(comm,"moving graph initial geometry",[&] {
		initial=iga::MovingCutGeometry::Build(grid,motion.Evaluate(0,0,.125),geometry_options);
	});
	iga::ImmersedMovingTransientDistributedRuntime runtime(comm,std::move(initial),options);
	std::vector<PetscScalar> seed(runtime.CommittedRowEnd()-runtime.CommittedRowBegin(),0);
	for(PetscInt row=runtime.CommittedRowBegin();row<runtime.CommittedRowEnd();++row)
		if(static_cast<std::size_t>(row)<runtime.CommittedLayout().NodeFieldRows()&&row%4==0)seed[row-runtime.CommittedRowBegin()]=2.96;
	runtime.SetCommittedOwnedState(seed);
	iga::ImmersedMovingDistributedGraphBackend::GeometryProvider provider=[&](const iga::MovingCutGeometry& accepted,const iga::DomainStepContext& step) {
		if(fail_geometry&&rank==0)throw std::runtime_error("injected moving graph geometry failure");
		const double end=step.EndTime()+(bad_time&&rank==0?.001:0);
		return iga::MovingCutGeometry::Build(grid,(reject_mass?contraction:motion).Evaluate(end,step.start_time_s,end),geometry_options,&accepted);
	};
	const iga::ImmersedMovingConservationLimits limits{1e-2,1e-8,1e-2,1e-2,1e-8};
	if(size>1) {
		auto different=limits;if(rank==0)different.reynolds*=2;
		Reject(comm,[&] { iga::ImmersedMovingDistributedGraphBackend rejected(runtime,provider,2,different); });
	}
	iga::ImmersedMovingDistributedGraphBackend backend(runtime,provider,2,limits);
	std::vector<iga::CouplingPort> ports;
	for(const auto& definition:options.ports) {
		iga::CouplingPort port;port.id=definition.id;port.subsystem_id="immersed";
		port.locator_kind="boundary_label";port.locator=std::to_string(definition.boundary_label);
		port.requires={iga::PortQuantity::FlowRate};
		port.provides={iga::PortQuantity::Area,iga::PortQuantity::FlowRate,iga::PortQuantity::MeanPressure,iga::PortQuantity::MeanNormalTraction};
		ports.push_back(port);
	}
	iga::ThreeDImmersedMovingDistributedFlowDomain domain("immersed",backend,ports);
	if(reject_mass) {
		iga::DomainStepContext step;step.start_time_s=0;step.dt_s=.125;step.step_index=0;
		domain.BeginStep(step);domain.SetPortInput("inlet",Input(.125,-1e-4));domain.SetPortInput("outlet",Input(.125,1e-4));
		std::string message;try { domain.SolveTrial(); } catch(const std::exception& error) { message=error.what(); }
		iga::CollectiveLocalStage(comm,"moving graph physical rejection evidence",[&] {
			Require(message.find("moving graph physical conservation")!=std::string::npos,"contracting graph did not fail its physical gate");
			const auto actual=domain.CommittedOwnedBackendState();
			Require(actual==seed&&!runtime.Clock().trial_active&&runtime.Clock().index==0
				&&domain.Diagnostics().committed_steps==0&&domain.CommittedPortStates().empty(),"physical failure published a trial or changed accepted state");
			for(const auto& port:runtime.PortDefinitions())Require(port.value==0,"physical failure did not restore controls");
		});
		reject_mass=false;domain.SolveTrial();domain.PrepareCommitStep();domain.FinalizeCommitStep();
		Require(runtime.Clock().index==1&&domain.Diagnostics().committed_steps==1,"physical failure retry did not commit");
		runtime.Close();
		std::cout<<"immersed_moving_domain_mass_rejection rank="<<rank<<" reason="<<message<<" accepted_unchanged=1 healthy_retry=1 passed\n";return;
	}
	std::unique_ptr<iga::MovingImmersedTransientFlowRuntime> oracle;
	iga::CollectiveLocalStage(comm,"moving graph serial oracle",[&] {
		if(rank)return;
		iga::MovingImmersedTransientFlowOptions settings;settings.grid=grid;settings.geometry=geometry_options;settings.extension_layers=2;
		settings.flow=options;settings.flow.solver_options_prefix="moving_graph_oracle_";
		settings.flow.ports[0].value=-1e-4;settings.flow.ports[1].value=1e-4;
		oracle=std::make_unique<iga::MovingImmersedTransientFlowRuntime>(motion.Evaluate(0,0,.125),settings);
		const auto& layout=oracle->CommittedLayout();std::vector<std::array<double,4>> fields(layout.NodeIds().size(),{{2.96,0,0,0}});
		oracle->InitializeCommittedGlobalState(iga::ImmersedGlobalFlowState(0,0,layout,fields,{0,0},true,0));
	});
	double maximum_field_error=0,maximum_relative_field_error=0,maximum_port_error=0;
	for(int step=0;step<2;++step) {
		iga::DomainStepContext context;context.start_time_s=.125*step;context.dt_s=.125;context.step_index=step;
		const auto accepted=domain.CommittedOwnedBackendState();
		const auto accepted_geometry=runtime.CommittedGeometry().PublicationIdentitySha256();
		const auto unchanged=[&] {
			iga::CollectiveLocalStage(comm,"moving graph accepted isolation",[&] {
				const auto actual=domain.CommittedOwnedBackendState();
				Require(actual.size()==accepted.size()&&std::memcmp(actual.data(),accepted.data(),actual.size()*sizeof(PetscScalar))==0,"graph changed accepted field or ownership");
				Require(runtime.Clock().index==static_cast<std::uint64_t>(step)&&runtime.Clock().time_s==context.start_time_s
					&&runtime.CommittedGeometry().PublicationIdentitySha256()==accepted_geometry,"graph changed accepted geometry or clock");
			});
		};
		const auto inputs=[&] { domain.SetPortInput("inlet",Input(context.EndTime(),-1e-4));domain.SetPortInput("outlet",Input(context.EndTime(),1e-4)); };
		domain.BeginStep(context);inputs();
		if(!step) {
			fail_geometry=true;Reject(comm,[&] { domain.SolveTrial(); });fail_geometry=false;unchanged();
			bad_time=true;Reject(comm,[&] { domain.SolveTrial(); });bad_time=false;unchanged();
			fail_force=true;Reject(comm,[&] { domain.SolveTrial(); });fail_force=false;unchanged();
		}
		domain.SolveTrial();unchanged();
		Require(runtime.TrialLayout().NodeIds()!=runtime.CommittedLayout().NodeIds(),"graph fixture did not change active nodes");
		if(!step) {
			if(rank==size-1)runtime.FailNextPrepareForTesting();
			Reject(comm,[&] { domain.PrepareCommitStep(); });domain.RollbackTrial();
		} else {
			domain.PrepareCommitStep();domain.AbortStep();domain.FinalizeCommitStep();domain.BeginStep(context);
		}
		unchanged();
		for(const auto& port:runtime.PortDefinitions())Require(port.value==(step?(port.id=="inlet"?-1e-4:1e-4):0),"graph rollback did not restore accepted controls");
		inputs();domain.SolveTrial();unchanged();
		domain.PrepareCommitStep();unchanged();
		const auto record=runtime.ConservationDiagnostics();
		domain.FinalizeCommitStep();domain.FinalizeCommitStep();
		Require(domain.Diagnostics().committed_steps==static_cast<std::size_t>(step+1)&&runtime.Clock().index==static_cast<std::uint64_t>(step+1),"graph commit clocks differ");
		Require(runtime.ConservationDiagnostics().target_publication_identity_sha256==record.target_publication_identity_sha256,"graph did not retain committed conservation");
		std::vector<PetscScalar> expected;
		iga::CollectiveLocalStage(comm,"moving graph oracle solve",[&] {
			if(rank)return;
			oracle->BeginTrial(motion.Evaluate(context.EndTime(),context.start_time_s,context.EndTime()),step+1,context.dt_s);
			Require(oracle->SolveTrial(),"moving graph oracle failed");oracle->Commit();expected=oracle->CommittedState();
			for(const auto& entry:domain.CommittedPortStates()) {
				const auto& list=oracle->CommittedDiagnostics().ports;
				const auto found=std::find_if(list.begin(),list.end(),[&](const auto& port) { return port.id==entry.first; });
				Require(found!=list.end(),"oracle port missing");const auto& a=entry.second;const auto& b=found->measurement;
				for(double difference:{*a.area_m2-b.area_m2,*a.outward_flow_m3_s-b.outward_flow_m3_s,*a.mean_pressure_pa-b.mean_pressure_pa,*a.mean_normal_traction_pa-b.mean_normal_traction_pa})
					maximum_port_error=std::max(maximum_port_error,std::abs(difference));
				Require(maximum_port_error<1e-8,"graph port differs from moving oracle");
			}
		});
		iga::CollectiveLocalStage(comm,"moving graph expected field storage",[&] { if(rank)expected.resize(runtime.CommittedLayout().Rows()); });
		MPI_Bcast(expected.data(),static_cast<int>(expected.size()),MPI_DOUBLE,0,comm);
		std::array<double,12> local{},global{};
		iga::CollectiveLocalStage(comm,"moving graph owned field comparison",[&] {
			const auto actual=domain.CommittedOwnedBackendState();
			for(PetscInt row=domain.RowBegin();row<domain.RowEnd();++row) {
				const auto field=static_cast<std::size_t>(row)<runtime.CommittedLayout().NodeFieldRows()?row%4:row+1==static_cast<PetscInt>(runtime.CommittedLayout().Rows())?5:4;
				const double delta=actual[row-domain.RowBegin()]-expected[row];
				local[2*field]+=delta*delta;local[2*field+1]+=expected[row]*expected[row];
			}
		});
		MPI_Allreduce(local.data(),global.data(),12,MPI_DOUBLE,MPI_SUM,comm);
		for(int field=0;field<6;++field) {
			const double error=std::sqrt(global[2*field]),norm=std::sqrt(global[2*field+1]);
			maximum_field_error=std::max(maximum_field_error,error/std::max(1.,norm));
			if(norm>1e-10) {
				const double relative=error/norm;maximum_relative_field_error=std::max(maximum_relative_field_error,relative);
				Require(std::isfinite(relative)&&relative<1e-6,"moving graph relative field error exceeds gate");
			} else Require(error<1e-10,"moving graph near-zero field exceeds absolute gate");
		}
		Require(std::isfinite(maximum_field_error)&&maximum_field_error<1e-8,"moving graph field differs from serial");
		if(!rank)std::cout<<"moving_graph_step="<<step<<" field_scaled_l2="<<maximum_field_error
			<<" reynolds="<<record.normalized_reynolds_defect<<" mass="<<record.normalized_moving_mass_defect
			<<" leakage="<<record.normalized_wall_relative_leakage<<" continuity="<<record.endpoint.normalized_discrete_moving_wall_continuity_defect<<std::endl;
	}
	runtime.Close();oracle.reset();MPI_Bcast(&maximum_port_error,1,MPI_DOUBLE,0,comm);
	std::cout<<"immersed_moving_domain_mpi rank="<<rank<<" ranks="<<size<<" steps="<<domain.Diagnostics().committed_steps
		<<" field_scaled_l2="<<maximum_field_error<<" field_relative_l2="<<maximum_relative_field_error<<" port_max_abs="<<maximum_port_error<<" physical_gates=1 active_set_changes=2 passed\n";
}
}
int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);int status=0;
	try { Run(PETSC_COMM_WORLD,argc>1&&std::string(argv[1])=="reject-mass"); } catch(const std::exception& error) { std::cerr<<error.what()<<'\n';status=1; }
	int global=0;MPI_Allreduce(&status,&global,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);PetscFinalize();return global;
}
