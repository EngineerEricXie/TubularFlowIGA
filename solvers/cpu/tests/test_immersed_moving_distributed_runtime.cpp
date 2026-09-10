#include "ImmersedMovingTransientDistributedRuntime.hpp"
#include "PrescribedSurfaceMotion.hpp"
#include "MovingImmersedTransientFlowRuntime.hpp"
#include <iostream>
#include <cstring>
#include <optional>

namespace {
iga::RawSurfaceSoup Cube()
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
void Require(bool condition, const char* message)
{
	if (!condition) throw std::runtime_error(message);
}
template<class Function> void Reject(MPI_Comm comm, Function&& function)
{
	std::string message;
	try { function(); } catch (const std::exception& error) { message = error.what(); }
	iga::CollectiveLocalStage(comm,"history test rejection",[&] { Require(!message.empty(),"missing history rejection"); });
	iga::RequireCollectiveSameText(comm,"history rejection agreement",message);
}
class ScopedGlobalOption
{
public:
	ScopedGlobalOption(const char* key,const char* value) : key_(key)
	{
		const char* previous=nullptr;PetscBool found=PETSC_FALSE;
		Require(!PetscOptionsFindPair(nullptr,nullptr,key,&previous,&found),"cannot read global option");
		had_key_=found;had_value_=previous!=nullptr;if(previous)previous_=previous;
		Require(!PetscOptionsSetValue(nullptr,key,value),"cannot replace global option");
	}
	~ScopedGlobalOption()
	{
		if(had_key_)PetscOptionsSetValue(nullptr,key_,had_value_?previous_.c_str():nullptr);
		else PetscOptionsClearValue(nullptr,key_);
	}
	ScopedGlobalOption(const ScopedGlobalOption&)=delete;
	ScopedGlobalOption& operator=(const ScopedGlobalOption&)=delete;
private:
	const char* key_;std::string previous_;bool had_key_=false,had_value_=false;
};
std::vector<PetscScalar> Owned(Vec vector)
{
	iga::PetscReadArray array;array.Acquire(vector);
	PetscInt count=0;if(VecGetLocalSize(vector,&count))throw std::runtime_error("local size failed");
	std::vector<PetscScalar> result(array.Data(),array.Data()+count);array.Restore();return result;
}
void Run(MPI_Comm comm,bool changing)
{
	int rank=0;MPI_Comm_rank(comm,&rank);
	auto soup=Cube();for(auto& x:soup.vertices)for(double& c:x)c+=.1;
	auto end=soup;for(auto& x:end.vertices)x[0]+=changing ? .74 : .04;
	iga::PrescribedSurfaceMotion motion({{0,soup},{.25,end}});
	iga::MovingCutGeometryOptions geometry_options;geometry_options.volume.max_depth=2;
	geometry_options.volume_storage=iga::CutCellVolumeQuadratureStorageMode::Compact;
	const iga::CubicCartesianGridSpec grid{{{0,0,0}},{{changing ? 2.4 : 1.2,1.2,1.2}},{{changing ? 8u : 4u,3,3}}};
	const double speed=changing ? 2.96 : .16;
	const std::uint32_t layers=changing ? 2u : 1u;
	std::unique_ptr<iga::MovingCutGeometry> initial;
	iga::CollectiveLocalStage(comm,"moving runtime initial geometry",[&] {
		initial=iga::MovingCutGeometry::Build(grid,motion.Evaluate(0,0,.125),geometry_options);
	});
	iga::ImmersedTransientFlowOptions options;options.parameters={1,.1,0};options.wall_labels={7,8,9};
	bool fail_force=false;
	options.body_force=[&](const std::array<double,3>&) {
		if(fail_force&&rank==0)throw std::runtime_error("injected moving force freeze failure");
		return std::array<double,3>{{0,0,0}};
	};
	iga::ImmersedMovingTransientDistributedRuntime runtime(comm,std::move(initial),options);
	// Root-only serial oracle and its broadcast are confined to this small
	// regression; the production runtime never constructs a global field.
	std::unique_ptr<iga::MovingImmersedTransientFlowRuntime> oracle;
	iga::CollectiveLocalStage(comm,"moving serial oracle initial state",[&] {
		if(rank)return;
		iga::MovingImmersedTransientFlowOptions reference;
		reference.grid=grid;reference.geometry=geometry_options;reference.flow=options;reference.extension_layers=layers;
		reference.flow.solver_options_prefix="moving_oracle_";
		oracle=std::make_unique<iga::MovingImmersedTransientFlowRuntime>(motion.Evaluate(0,0,.125),reference);
		const auto& layout=oracle->CommittedLayout();
		std::vector<std::array<double,4>> fields(layout.NodeIds().size(),{{speed,0,0,0}});
		oracle->InitializeCommittedGlobalState(iga::ImmersedGlobalFlowState(0,0,layout,fields,{},true,0));
	});
	double maximum_field_error=0;PetscInt total_newton=0,total_ksp=0;
	const auto compare=[&](double start,std::uint64_t index) {
		std::vector<PetscScalar> expected;
		iga::CollectiveLocalStage(comm,"moving serial oracle solve",[&] {
			if(rank)return;
			oracle->BeginTrial(motion.Evaluate(start+.125,start,start+.125),index,.125);
			Require(oracle->SolveTrial(),"serial moving oracle did not converge");
			expected=oracle->TrialState();
			Require(oracle->TrialLayout().HashSha256()==runtime.TrialLayout().HashSha256(),"serial moving layout differs");
		});
		iga::CollectiveLocalStage(comm,"moving oracle receive storage",[&] {
			if(rank)expected.resize(runtime.TrialLayout().Rows());
			Require(expected.size()<=static_cast<std::size_t>(std::numeric_limits<int>::max()),"oracle MPI count overflow");
		});
		MPI_Bcast(expected.data(),static_cast<int>(expected.size()),MPI_DOUBLE,0,comm);
		std::array<double,10> local{},global{};
		iga::CollectiveLocalStage(comm,"moving solved field comparison",[&] {
			PetscInt begin=0,end=0;Require(!VecGetOwnershipRange(runtime.TrialState(),&begin,&end),"trial ownership lookup failed");
			const auto actual=Owned(runtime.TrialState());
			for(PetscInt row=begin;row<end;++row) {
				const auto field=static_cast<std::size_t>(row)<runtime.TrialLayout().NodeFieldRows()?row%4:4;
				const double error=actual[row-begin]-expected[row];
				Require(std::isfinite(error),"nonfinite moving solved field error");
				local[2*field]+=error*error;local[2*field+1]+=expected[row]*expected[row];
			}
		});
		MPI_Allreduce(local.data(),global.data(),10,MPI_DOUBLE,MPI_SUM,comm);
		for(int field=0;field<5;++field)maximum_field_error=std::max(maximum_field_error,std::sqrt(global[2*field])/std::max(1.,std::sqrt(global[2*field+1])));
		Require(std::isfinite(maximum_field_error)&&maximum_field_error<1e-8,"moving solved field differs from serial oracle");
		total_newton+=runtime.TrialDiagnostics().nonlinear_iterations;total_ksp+=runtime.TrialDiagnostics().ksp_iterations;
	};
	PetscInt first=0,last=0;VecGetOwnershipRange(runtime.CommittedState(),&first,&last);
	std::vector<PetscScalar> seed(last-first,0.);
	for(PetscInt row=first;row<last;++row)if(static_cast<std::size_t>(row)<runtime.CommittedLayout().NodeFieldRows()&&row%4==0)seed[row-first]=speed;
	runtime.SetCommittedOwnedState(seed);
	const auto original_geometry=runtime.CommittedGeometry().GeometryIdentitySha256();
	const auto build=[&](double start) {
		std::unique_ptr<iga::MovingCutGeometry> result;
		iga::CollectiveLocalStage(comm,"moving runtime target geometry",[&] {
			result=iga::MovingCutGeometry::Build(grid,motion.Evaluate(start+.125,start,start+.125),geometry_options,&runtime.CommittedGeometry());
		});
		return result;
	};
	const auto unchanged=[&](const std::vector<PetscScalar>& expected,double time,std::uint64_t index) {
		iga::CollectiveLocalStage(comm,"moving runtime committed isolation",[&] {
			const auto actual=Owned(runtime.CommittedState());
			Require(actual.size()==expected.size()&&std::memcmp(actual.data(),expected.data(),actual.size()*sizeof(PetscScalar))==0,"accepted vector changed during trial");
			Require(runtime.Clock().time_s==time&&runtime.Clock().index==index,"accepted clock changed during trial");
		});
	};
	Reject(comm,[&] { runtime.BeginTrial(build(0),2,layers); });unchanged(seed,0,0);
	fail_force=true;Reject(comm,[&] { runtime.BeginTrial(build(0),1,layers); });fail_force=false;unchanged(seed,0,0);
	Require(!runtime.Clock().trial_active,"failed history freeze published a trial");
	runtime.BeginTrial(build(0),1,layers);
	if(changing)Require(runtime.TrialLayout().NodeIds()!=runtime.CommittedLayout().NodeIds(),"first moving step did not change active nodes");
	Require(runtime.TrialHistory().HasMappedSource(),"moving runtime did not map source history");
	Reject(comm,[&] { runtime.PrepareCommit(); });unchanged(seed,0,0);
	Require(runtime.SolveTrial(),"moving runtime did not converge");unchanged(seed,0,0);
	if(rank==0)runtime.FailNextPrepareForTesting();
	Reject(comm,[&] { runtime.PrepareCommit(); });unchanged(seed,0,0);
	runtime.AbortTrial();runtime.AbortTrial();unchanged(seed,0,0);
	Require(runtime.CommittedGeometry().GeometryIdentitySha256()==original_geometry,"abort changed accepted geometry");
	runtime.BeginTrial(build(0),1,layers);Require(runtime.SolveTrial(),"moving retry did not converge");
	compare(0,1);
	const auto first_trial=Owned(runtime.TrialState());
	const auto first_geometry=runtime.TrialGeometry().GeometryIdentitySha256();
	const auto first_publication=runtime.TrialGeometry().PublicationIdentitySha256();
	const auto first_layout=runtime.TrialLayout().HashSha256();
	runtime.PrepareCommit();unchanged(seed,0,0);runtime.FinalizeCommit();runtime.FinalizeCommit();
	unchanged(first_trial,.125,1);
	iga::CollectiveLocalStage(comm,"first moving epoch publication",[&] {
		Require(runtime.CommittedGeometry().GeometryIdentitySha256()==first_geometry
			&&runtime.CommittedGeometry().PublicationIdentitySha256()==first_publication
			&&runtime.CommittedLayout().HashSha256()==first_layout,"first committed epoch differs from prepared epoch");
	});
	iga::CollectiveLocalStage(comm,"moving oracle first commit",[&] { if(!rank)oracle->Commit(); });
	Require(runtime.Clock().time_s==.125&&runtime.Clock().index==1&&!runtime.Clock().trial_active,"first publication clock differs");
	const auto accepted=Owned(runtime.CommittedState());
	if(changing) {
		Reject(comm,[&] { runtime.BeginTrial(build(.125),2,1); });unchanged(accepted,.125,1);
		Require(!runtime.Clock().trial_active,"insufficient extension band published a trial");
	}
	{
		std::optional<ScopedGlobalOption> changed_flow,changed_extension;
		iga::CollectiveLocalStage(comm,"moving global options mutation",[&] {
			changed_flow.emplace("-immersed_moving_ksp_type",rank==0?"unavailable_first":"unavailable_other");
			changed_extension.emplace("-immersed_extension_ksp_type",rank==0?"unavailable_first":"unavailable_other");
		});
		// Both the new Newton solver and extension must use the constructor
		// snapshot, even while invalid replacements exist in the global DB.
		runtime.BeginTrial(build(.125),2,layers);
	}
	if(changing)Require(runtime.TrialLayout().NodeIds()!=runtime.CommittedLayout().NodeIds(),"second moving step did not change active nodes");
	Require(runtime.SolveTrial(),"second moving step did not converge");
	runtime.PrepareCommit();runtime.AbortPrepared();runtime.Rollback();unchanged(accepted,.125,1);
	Require(runtime.SolveTrial(),"prepared rollback retry did not converge");
	compare(.125,2);
	const auto second_trial=Owned(runtime.TrialState());
	const auto second_geometry=runtime.TrialGeometry().GeometryIdentitySha256();
	const auto second_publication=runtime.TrialGeometry().PublicationIdentitySha256();
	const auto second_layout=runtime.TrialLayout().HashSha256();
	runtime.Commit();unchanged(second_trial,.25,2);
	iga::CollectiveLocalStage(comm,"second moving epoch publication",[&] {
		Require(runtime.CommittedGeometry().GeometryIdentitySha256()==second_geometry
			&&runtime.CommittedGeometry().PublicationIdentitySha256()==second_publication
			&&runtime.CommittedLayout().HashSha256()==second_layout&&!runtime.Clock().trial_active,"second committed epoch differs from prepared epoch");
	});
	runtime.Close();runtime.Close();Reject(comm,[&] { runtime.AbortTrial(); });
	std::cout<<"immersed_moving_runtime_mpi rank="<<rank<<" changing="<<changing<<" extension_layers="<<layers<<" maximum_field_scaled_l2="<<maximum_field_error
		<<" newton_iterations="<<total_newton<<" ksp_iterations="<<total_ksp
		<<" accepted_steps=2 freeze_failure=1 prepare_failure=1 abort_retry=1 exact_publication=1 global_options_frozen=1 passed\n";
}
}
int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);int status=0,rank=0;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	try { Run(PETSC_COMM_WORLD,argc>1&&std::string(argv[1])=="changing"); } catch(const std::exception& error) { std::cerr<<"rank "<<rank<<": "<<error.what()<<'\n';status=1; }
	int global=0;MPI_Allreduce(&status,&global,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);PetscFinalize();return global;
}
