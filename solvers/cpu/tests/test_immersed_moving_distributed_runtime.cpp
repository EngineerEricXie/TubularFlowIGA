#include "ImmersedMovingTransientDistributedRuntime.hpp"
#include "MaterialSurfaceCheckpoint.hpp"
#include "MovingCheckpointRestore.hpp"
#include "CollectiveCheckpointReceipts.hpp"
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
iga::MaterialSurfacePatchMap BottomPatch(const iga::MaterialSurfaceKinematics& input)
{
	std::vector<iga::RawSurfaceTriangle> topology;
	for (const auto& source : input.SourceTriangles()) {
		iga::RawSurfaceTriangle triangle;triangle.boundary_id=source.boundary_id;
		for (int corner=0;corner<3;++corner)triangle.indices[corner]=source.source_vertex_indices[corner];
		topology.push_back(triangle);
	}
	const auto& reference=input.ReferenceMaterialVerticesM();
	const auto material=iga::MaterialSurfaceKinematics::CreateFromSourceTopology(reference,reference,
		std::vector<std::array<double,3>>(reference.size(),{{0,0,0}}),topology,0.,-.125,0.);
	iga::DistributedSurfaceLayout layout;
	layout.global_node_count=4;layout.partition_count=1;layout.partition_rank=0;
	layout.owned_global_node_ids={10,11,12,13};
	std::vector<iga::MaterialSurfacePatchMap::GlobalToSourceVertex> vertices;
	for (std::uint32_t node=0;node<4;++node) {
		vertices.emplace_back(10+node,node);
		layout.reference_positions.push_back({10+node,material.ReferenceMaterialVerticesM()[node]});
	}
	layout.reference_triangles={{{10,12,11}},{{10,13,12}}};
	layout.owned_reference_lumped_areas_m2={1./3.,1./6.,1./3.,1./6.};
	layout.reference_mesh_identity_sha256=iga::MaterialSurfacePatchMap::BuildReferenceIdentitySha256(
		material,8,vertices,layout.reference_triangles,{0,1},layout.owned_global_node_ids);
	layout.layout_identity_sha256=iga::BuildDistributedSurfaceLayoutIdentitySha256(layout);
	iga::DistributedSurfaceInterface interface;interface.id={"structure","membrane","bottom"};interface.subsystem_id="membrane";
	interface.boundary_labels={8};interface.reference_mesh_identity_sha256=layout.reference_mesh_identity_sha256;
	interface.provides={iga::SurfaceFieldQuantity::Displacement,iga::SurfaceFieldQuantity::Velocity};
	interface.requires={iga::SurfaceFieldQuantity::TractionOnStructure};
	return iga::MaterialSurfacePatchMap::Create(interface,layout,material,8,vertices,{0,1},layout.owned_global_node_ids);
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
	double maximum_field_error=0,maximum_conservation_error=0;PetscInt total_newton=0,total_ksp=0;
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
		const auto measured=runtime.ConservationDiagnostics();
		iga::CollectiveLocalStage(comm,"moving transition conservation comparison",[&] {
			Require(measured.source_index+1==index&&measured.target_index==index&&measured.source_time_s==start
				&&measured.target_time_s==start+.125&&measured.dt_s==.125,"moving conservation clock differs");
			Require(measured.source_geometry_identity_sha256==runtime.CommittedGeometry().GeometryIdentitySha256()
				&&measured.target_geometry_identity_sha256==runtime.TrialGeometry().GeometryIdentitySha256()
				&&measured.source_publication_identity_sha256==runtime.CommittedGeometry().PublicationIdentitySha256()
				&&measured.target_publication_identity_sha256==runtime.TrialGeometry().PublicationIdentitySha256(),"moving conservation geometry differs");
			if(rank)return;
			const auto expected=oracle->ConservationDiagnostics();
			const std::array<double,12> actual{{measured.source_audited_volume_m3,measured.target_audited_volume_m3,
				measured.backward_euler_volume_rate_m3_s,measured.reynolds_defect_m3_s,measured.moving_mass_defect_m3_s,
				measured.normalization_scale_m3_s,measured.normalized_divergence_theorem_defect,measured.normalized_reynolds_defect,
				measured.normalized_moving_mass_defect,measured.normalized_wall_relative_leakage,
				measured.endpoint.total_material_surface_outward_flow_m3_s,measured.endpoint.discrete_moving_wall_continuity_defect_m3_s}};
			const std::array<double,12> reference{{expected.source_audited_volume_m3,expected.target_audited_volume_m3,
				expected.backward_euler_volume_rate_m3_s,expected.reynolds_defect_m3_s,expected.moving_mass_defect_m3_s,
				expected.normalization_scale_m3_s,expected.normalized_divergence_theorem_defect,expected.normalized_reynolds_defect,
				expected.normalized_moving_mass_defect,expected.normalized_wall_relative_leakage,
				expected.total_material_surface_outward_flow_m3_s,expected.discrete_moving_wall_continuity_defect_m3_s}};
			for(std::size_t i=0;i<actual.size();++i)maximum_conservation_error=std::max(maximum_conservation_error,std::abs(actual[i]-reference[i])/std::max(1.,std::abs(reference[i])));
			Require(std::isfinite(maximum_conservation_error)&&maximum_conservation_error<1e-8,"moving conservation differs from serial");
		});
		total_newton+=runtime.TrialDiagnostics().nonlinear_iterations;total_ksp+=runtime.TrialDiagnostics().ksp_iterations;
	};
	PetscInt first=0,last=0;VecGetOwnershipRange(runtime.CommittedState(),&first,&last);
	std::vector<PetscScalar> seed(last-first,0.);
	for(PetscInt row=first;row<last;++row)if(static_cast<std::size_t>(row)<runtime.CommittedLayout().NodeFieldRows()&&row%4==0)seed[row-first]=speed;
	runtime.SetCommittedOwnedState(seed);
	Reject(comm,[&] { (void)runtime.ConservationDiagnostics(); });
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
	const auto patch=BottomPatch(runtime.CommittedGeometry().Evaluation());
	int ranks=1;MPI_Comm_size(comm,&ranks);
	auto patch_layout=patch.Layout();patch_layout.partition_count=ranks;patch_layout.partition_rank=rank;
	patch_layout.owned_global_node_ids.clear();patch_layout.owned_reference_lumped_areas_m2.clear();
	for(std::size_t node=0;node<4;++node)if(static_cast<int>(node%ranks)==rank) {
		patch_layout.owned_global_node_ids.push_back(patch.Layout().owned_global_node_ids[node]);
		patch_layout.owned_reference_lumped_areas_m2.push_back(patch.Layout().owned_reference_lumped_areas_m2[node]);
	}
	patch_layout.layout_identity_sha256=iga::BuildDistributedSurfaceLayoutIdentitySha256(patch_layout);
	const auto distributed_patch=iga::MaterialSurfacePatchMap::Create(patch.Interface(),patch_layout,patch.FullReference(),8,
		patch.GlobalToSourceVertices(),patch.LayoutTriangleToSourceTriangles(),patch.ConfiguredClampedGlobalNodeIds());
	auto fluid_interface=patch.Interface();fluid_interface.id={"fluid","moving","bottom"};fluid_interface.subsystem_id="moving";
	fluid_interface.provides={iga::SurfaceFieldQuantity::TractionOnStructure};
	fluid_interface.requires={iga::SurfaceFieldQuantity::Displacement,iga::SurfaceFieldQuantity::Velocity};
	const auto check_traction=[&](double start,std::uint64_t index) {
		iga::FsiTrialContext context;context.step=index;context.start_time_s=start;context.dt_s=.125;context.coupling_iteration=2;
		const auto& geometry=runtime.TrialGeometry();
		const auto expected_material=geometry.Evaluation().ContentIdentitySha256();
		for(int mode=0;mode<2;++mode) {
			auto expected=context;auto identity=expected_material;
			if(rank==0&&mode==0)++expected.step;
			if(rank==0&&mode==1)identity=std::string(64,'0');
			Reject(comm,[&] { (void)runtime.BuildTrialSurfaceTraction(distributed_patch,fluid_interface,expected,identity,ranks-1); });
		}
		if(rank==0)std::cerr<<"moving_stage traction_begin step="<<index<<"\n";
		const auto publication=runtime.BuildTrialSurfaceTraction(distributed_patch,fluid_interface,context,expected_material,ranks-1);
		if(rank==0)std::cerr<<"moving_stage traction_complete step="<<index<<"\n";
		iga::CollectiveLocalStage(comm,"moving traction publication stamp",[&] {
		Require(publication.stamp.step==index&&publication.stamp.time_s==context.EndTime()
			&&publication.stamp.coupling_iteration==2,"moving traction publication stamp differs");
		});
		std::array<double,24> expected{};
		iga::CollectiveLocalStage(comm,"moving traction serial oracle",[&] {
			if(rank)return;
			const auto fields=oracle->TrialState();const auto& layout=oracle->TrialLayout();
			std::vector<iga::FluidSurfaceElementState> state;
			for(std::uint64_t cell=0;cell<geometry.Domain().Cells().size();++cell) {
				if(geometry.Domain().Cells()[cell].classification!=iga::CellClassification::Cut)continue;
				const auto& points=geometry.Surface().UsableRule(geometry.Domain(),cell).Points();
				if(!std::any_of(points.begin(),points.end(),[](const auto& p){return p.boundary_id==8;}))continue;
				const auto element=geometry.Domain().Background().MaterializeElement(cell);
				iga::FluidSurfaceElementState item;item.cell_id=cell;item.nodal_state.resize(element.connectivity.size());
				for(std::size_t node=0;node<element.connectivity.size();++node)for(int field=0;field<4;++field)
					item.nodal_state[node][field]=PetscRealPart(fields[4*layout.LocalNode(element.connectivity[node])+field]);
				state.push_back(std::move(item));
			}
			iga::SurfaceFieldStamp stamp;stamp.time_s=context.EndTime();stamp.step=index;stamp.coupling_iteration=2;
			stamp.reference_mesh_identity_sha256=patch.ReferenceIdentitySha256();stamp.layout_identity_sha256=patch.Layout().layout_identity_sha256;
			stamp.partition_identity_sha256=iga::BuildDistributedSurfacePartitionIdentitySha256(patch.Layout());
			stamp.producer_state_identity_sha256=iga::BuildFluidSurfaceTractionStateIdentitySha256(geometry.Domain(),geometry.Surface(),geometry.Evaluation(),patch,.1,state);
			const auto serial=iga::BuildFluidSurfaceTraction(geometry.Domain(),geometry.Surface(),geometry.Evaluation(),
				fluid_interface,patch.Layout(),patch,stamp,.1,state);
			for(std::size_t node=0;node<4;++node)for(int axis=0;axis<3;++axis) {
				expected[6*node+axis]=serial.traction.consistent_nodal_force_n[node][axis];
				expected[6*node+3+axis]=serial.traction.traction_on_structure_pa[node][axis];
			}
		});
		MPI_Bcast(expected.data(),expected.size(),MPI_DOUBLE,0,comm);
		iga::CollectiveLocalStage(comm,"moving traction oracle comparison",[&] {
			for(std::size_t local=0;local<patch_layout.owned_global_node_ids.size();++local) {
				const auto node=patch_layout.owned_global_node_ids[local]-10;
				for(int axis=0;axis<3;++axis) {
					Require(std::abs(publication.consistent_nodal_force_n[local][axis]-expected[6*node+axis])<1e-8,"moving nodal force differs from serial oracle");
					Require(std::abs(publication.traction_on_structure_pa[local][axis]-expected[6*node+3+axis])<1e-8,"moving traction differs from serial oracle");
				}
			}
		});
	};

	Reject(comm,[&] { (void)runtime.CaptureTrialPatchState(patch); });
	Reject(comm,[&] { runtime.BeginTrial(build(0),2,layers); });unchanged(seed,0,0);
	fail_force=true;Reject(comm,[&] { runtime.BeginTrial(build(0),1,layers); });fail_force=false;unchanged(seed,0,0);
	Require(!runtime.Clock().trial_active,"failed history freeze published a trial");
	runtime.BeginTrial(build(0),1,layers);
	if(changing)Require(runtime.TrialLayout().NodeIds()!=runtime.CommittedLayout().NodeIds(),"first moving step did not change active nodes");
	Reject(comm,[&] { (void)runtime.CaptureTrialPatchState(patch); });
	Require(runtime.TrialHistory().HasMappedSource(),"moving runtime did not map source history");
	Reject(comm,[&] { runtime.PrepareCommit(); });unchanged(seed,0,0);
	Require(runtime.SolveTrial(),"moving runtime did not converge");unchanged(seed,0,0);
	if(rank==0)std::cerr<<"moving_stage first_solve_complete\n";
	const auto first_patch=runtime.CaptureTrialPatchState(patch);
	const auto repeated_patch=runtime.CaptureTrialPatchState(patch);
	iga::CollectiveLocalStage(comm,"moving identical trial capture",[&] {
		Require(first_patch.size()==repeated_patch.size(),"repeated capture coverage differs");
		for(std::size_t i=0;i<first_patch.size();++i)
			Require(first_patch[i].cell_id==repeated_patch[i].cell_id && first_patch[i].nodal_state==repeated_patch[i].nodal_state,"same trial capture changed coefficients");
	});
	std::uint64_t local_patch_cells=first_patch.size(),global_patch_cells=0;
	MPI_Allreduce(&local_patch_cells,&global_patch_cells,1,MPI_UINT64_T,MPI_SUM,comm);
	Require(global_patch_cells>0,"moving trial patch capture was empty globally");
	// Structure preparation can fail after the fluid is already prepared.
	// Revoke that publication locally and reuse the accepted epoch safely.
	runtime.PrepareCommit();
	runtime.AbortTrialDeferredNoexcept();unchanged(seed,0,0);
	Reject(comm,[&] { (void)runtime.ConservationDiagnostics(); });
	runtime.BeginTrial(build(0),1,layers);
	Require(runtime.SolveTrial(),"prepared deferred-abort retry did not converge");
	if(rank==0)runtime.FailNextPrepareForTesting();
	Reject(comm,[&] { runtime.PrepareCommit(); });unchanged(seed,0,0);
	static_assert(noexcept(runtime.AbortTrialDeferredNoexcept()),"paired FSI abort must be noexcept");
	runtime.AbortTrialDeferredNoexcept();runtime.AbortTrialDeferredNoexcept();unchanged(seed,0,0);
	Require(!runtime.Clock().trial_active && runtime.Clock().target_time_s==0 && runtime.Clock().target_index==0,"deferred abort retained a trial clock");
	Reject(comm,[&] { (void)runtime.CaptureTrialPatchState(patch); });
	Require(runtime.CommittedGeometry().GeometryIdentitySha256()==original_geometry,"abort changed accepted geometry");
	runtime.BeginTrial(build(0),1,layers);Require(runtime.SolveTrial(),"moving retry did not converge");
	const auto retry_patch=runtime.CaptureTrialPatchState(patch);
	// Accepted state and repeat extraction are bitwise checks. A fresh MPI
	// nonlinear solve uses the same per-field numerical gate as the oracle.
	std::array<double,8> patch_error{},global_patch_error{};
	double maximum_patch_retry_scaled_l2=0;
	iga::CollectiveLocalStage(comm,"moving patch retry comparison",[&] {
		Require(retry_patch.size()==first_patch.size(),"moving patch retry coverage differs");
		for(std::size_t i=0;i<first_patch.size();++i) {
			Require(retry_patch[i].cell_id==first_patch[i].cell_id && retry_patch[i].nodal_state.size()==first_patch[i].nodal_state.size(),"moving patch retry topology differs");
			for(std::size_t node=0;node<first_patch[i].nodal_state.size();++node)for(int field=0;field<4;++field) {
				const double before=first_patch[i].nodal_state[node][field],after=retry_patch[i].nodal_state[node][field];
				Require(std::isfinite(before)&&std::isfinite(after),"nonfinite patch retry coefficient");
				patch_error[2*field]+=(before-after)*(before-after);patch_error[2*field+1]+=before*before;
			}
		}
	});
	MPI_Allreduce(patch_error.data(),global_patch_error.data(),8,MPI_DOUBLE,MPI_SUM,comm);
	iga::CollectiveLocalStage(comm,"moving patch retry numerical gate",[&] {
		for(int field=0;field<4;++field) {
			const double error=std::sqrt(global_patch_error[2*field])/std::max(1.,std::sqrt(global_patch_error[2*field+1]));
			Require(std::isfinite(error)&&error<1e-8,"moving patch retry differs beyond existing field gate");
			maximum_patch_retry_scaled_l2=std::max(maximum_patch_retry_scaled_l2,error);
		}
	});
	if(rank==0)std::cerr<<"moving_stage first_oracle_begin\n";
	compare(0,1);check_traction(0,1);
	const auto first_trial=Owned(runtime.TrialState());
	const auto first_geometry=runtime.TrialGeometry().GeometryIdentitySha256();
	const auto first_publication=runtime.TrialGeometry().PublicationIdentitySha256();
	const auto first_layout=runtime.TrialLayout().HashSha256();
	runtime.PrepareCommit();unchanged(seed,0,0);
	Reject(comm,[&] { (void)runtime.CaptureTrialPatchState(patch); });
	const auto prepared_conservation=runtime.ConservationDiagnostics();
	runtime.FinalizeCommit();runtime.FinalizeCommit();
	Reject(comm,[&] { (void)runtime.CaptureTrialPatchState(patch); });
	const auto committed_conservation=runtime.ConservationDiagnostics();
	iga::CollectiveLocalStage(comm,"moving conservation checkpoint roundtrip",[&] {
		const auto bytes=iga::SerializeMovingConservationCheckpoint(committed_conservation);
		iga::Sha256 hash;hash.Append(bytes.data(),bytes.size());const auto identity=hash.Hex();
		const auto parse=[&](std::string_view input,const std::string& authority,const std::string& geometry,std::uint64_t step) {
			return iga::ParseMovingConservationCheckpoint(input,authority,geometry,
				committed_conservation.target_publication_identity_sha256,step,runtime.Clock().time_s);
		};
		const auto restored=parse(bytes,identity,runtime.CommittedGeometry().GeometryIdentitySha256(),runtime.Clock().index);
		Require(iga::SerializeMovingConservationCheckpoint(restored)==bytes,"conservation roundtrip differs");
		const auto rejects=[&](auto action) {
			bool failed=false;try { action(); }catch(const std::exception&) { failed=true; }
			Require(failed,"invalid conservation checkpoint accepted");
		};
		rejects([&] { (void)parse(bytes,std::string(64,'0'),restored.target_geometry_identity_sha256,restored.target_index); });
		rejects([&] { (void)parse(bytes,identity,std::string(64,'0'),restored.target_index); });
		rejects([&] { (void)parse(bytes,identity,restored.target_geometry_identity_sha256,restored.target_index+1); });
		for(auto damaged:{bytes.substr(0,bytes.size()-1),bytes+"x"}) {
			iga::Sha256 changed;changed.Append(damaged.data(),damaged.size());
			rejects([&] { (void)parse(damaged,changed.Hex(),restored.target_geometry_identity_sha256,restored.target_index); });
		}
	});

	Require(prepared_conservation.target_publication_identity_sha256==committed_conservation.target_publication_identity_sha256
		&&prepared_conservation.source_publication_identity_sha256==committed_conservation.source_publication_identity_sha256
		&&prepared_conservation.moving_mass_defect_m3_s==committed_conservation.moving_mass_defect_m3_s
		&&prepared_conservation.endpoint.material_surface_outward_flow_by_boundary_label_m3_s==committed_conservation.endpoint.material_surface_outward_flow_by_boundary_label_m3_s,
		"publication changed retained conservation");
	unchanged(first_trial,.125,1);
	iga::CollectiveLocalStage(comm,"first moving epoch publication",[&] {
		Require(runtime.CommittedGeometry().GeometryIdentitySha256()==first_geometry
			&&runtime.CommittedGeometry().PublicationIdentitySha256()==first_publication
			&&runtime.CommittedLayout().HashSha256()==first_layout,"first committed epoch differs from prepared epoch");
	});
	iga::CollectiveLocalStage(comm,"moving oracle first commit",[&] { if(!rank)oracle->Commit(); });
	Require(runtime.Clock().time_s==.125&&runtime.Clock().index==1&&!runtime.Clock().trial_active,"first publication clock differs");
	const auto accepted=Owned(runtime.CommittedState());
	const auto checkpoint=runtime.CaptureAcceptedCheckpoint();
	std::unique_ptr<iga::MovingCutGeometry> restored_geometry;
	iga::CollectiveLocalStage(comm,"moving checkpoint geometry reconstruction",[&] {
		Require(checkpoint.previous_material.has_value()&&checkpoint.current_material.has_value(),"missing checkpoint material history");
		const auto previous_bytes=iga::SerializeMaterialSurfaceCheckpoint(*checkpoint.previous_material);
		const auto current_bytes=iga::SerializeMaterialSurfaceCheckpoint(*checkpoint.current_material);
		const auto previous=iga::ParseMaterialSurfaceCheckpoint(previous_bytes,checkpoint.previous_material->ContentIdentitySha256());
		const auto current=iga::ParseMaterialSurfaceCheckpoint(current_bytes,checkpoint.current_material->ContentIdentitySha256());
		auto predecessor=iga::MovingCutGeometry::Build(grid,previous,geometry_options);
		restored_geometry=iga::MovingCutGeometry::Build(grid,current,geometry_options,predecessor.get());
	});
	auto resumed=iga::AllocateCollectiveRuntime<iga::ImmersedMovingTransientDistributedRuntime>(comm,
		comm,std::move(restored_geometry),options,1);
	for(int mode=0;mode<9;++mode) {
		auto invalid=checkpoint;
		if(rank==0) {
			if(mode==0)invalid.layout_identity_sha256=std::string(64,'0');
			if(mode==1)invalid.field.row_end++;
			if(mode==2)invalid.time_s+=.125;
			if(mode==3)invalid.field.values.front()=std::numeric_limits<double>::infinity();
			if(mode==4)invalid.conservation.normalized_wall_relative_leakage=std::numeric_limits<double>::quiet_NaN();
			if(mode==5)invalid.conservation.source_index++;
			if(mode==6) { if(ranks>1)invalid.conservation.normalized_wall_relative_leakage+=1.;else invalid.current_material.reset(); }
			if(mode==7)invalid.previous_material.reset();
			if(mode==8)invalid.previous_material.emplace(*invalid.current_material);
		}
		Reject(comm,[&] { resumed->RestoreAcceptedCheckpoint(invalid); });
		Reject(comm,[&] { (void)resumed->CaptureAcceptedCheckpoint(); });
	}
	resumed->RestoreAcceptedCheckpoint(checkpoint);
	Require(Owned(resumed->CommittedState())==accepted,"restored owned field differs");
	Require(iga::SerializeMovingConservationCheckpoint(resumed->ConservationDiagnostics())
		==iga::SerializeMovingConservationCheckpoint(committed_conservation),"restored conservation differs");
	Reject(comm,[&] { resumed->RestoreAcceptedCheckpoint(checkpoint); });
	char bundle_root[4096]={};
	Require(!PetscOptionsGetString(nullptr,nullptr,"-moving_checkpoint_bundle_root",bundle_root,sizeof(bundle_root),nullptr),"cannot read bundle test path");
	if(bundle_root[0]) {
		PetscBool read_only=PETSC_FALSE;PetscInt source_ranks=ranks;
		Require(!PetscOptionsGetBool(nullptr,nullptr,"-moving_checkpoint_read_only",&read_only,nullptr),"cannot read bundle mode");
		Require(!PetscOptionsGetInt(nullptr,nullptr,"-moving_checkpoint_source_ranks",&source_ranks,nullptr),"cannot read source rank count");
		Require(source_ranks>0&&source_ranks<=4096,"invalid source rank count");
		iga::RequireCollectiveSameInt(comm,"moving bundle read mode",read_only);
		iga::RequireCollectiveSameInt(comm,"moving bundle source ranks",static_cast<int>(source_ranks));
		Require(read_only||source_ranks==ranks,"writer source rank count differs");
		const auto hash=[](const std::string& value) { return iga::coupled_checkpoint_detail::Hash(value); };
		const std::string configuration=hash(changing?"moving-changing-fixture-v1":"moving-fixed-fixture-v1");
		const iga::CoupledCheckpointEpoch epoch{hash(std::string(bundle_root)),{},
			{configuration,configuration,hash(std::to_string(source_ranks)),static_cast<std::uint32_t>(source_ranks)},1,.125,.125};
		iga::MovingCheckpointLayout layout;layout.domain_id="flow";
		for(int member=0;member<source_ranks;++member)layout.world_ranks.push_back(member);
		if(!read_only) {
			iga::CollectiveLocalStage(comm,"moving numerical bundle create",[&] {
				if(rank==0) {
					Require(std::filesystem::create_directory(bundle_root),"moving bundle fixture exists");
					iga::CreateCoupledCheckpointEpoch(bundle_root,epoch);
				}
			});
			std::vector<iga::CoupledCheckpointShard> local_receipts;
			iga::CollectiveLocalStage(comm,"moving numerical bundle write",[&] {
				local_receipts=iga::WriteMovingCheckpointShards(bundle_root,epoch,layout,rank,checkpoint,configuration);
			});
			const auto receipts=iga::GatherCheckpointReceipts(comm,epoch,local_receipts);
			iga::CollectiveLocalStage(comm,"moving numerical bundle publish",[&] {
				if(rank==0)(void)iga::PublishCoupledCheckpoint(bundle_root,epoch,iga::MovingCheckpointCatalog(layout,ranks),receipts);
			});
		}
		std::optional<iga::CoupledCheckpointManifest> manifest;
		iga::CollectiveLocalStage(comm,"moving numerical bundle discovery",[&] {
			const auto found=iga::FindLatestCoupledCheckpoint(bundle_root,epoch.compatibility,iga::MovingCheckpointCatalog(layout,static_cast<std::uint32_t>(source_ranks)));
			Require(found.latest.has_value()&&found.rejected.empty(),"moving bundle discovery failed");manifest=found.latest;
		});
		{
			auto bad=*manifest;
			const auto field=std::find_if(bad.shards.begin(),bad.shards.end(),[](const auto& shard) { return shard.spec.format=="owned-real-field-v1"; });
			Require(field!=bad.shards.end(),"missing fixture shard");field->sha256[0]=field->sha256[0]=='a'?'b':'a';
			Reject(comm,[&] { (void)iga::RestoreMovingCheckpointRuntime(comm,bundle_root,bad,layout,configuration,grid,geometry_options,options); });
			unchanged(accepted,.125,1);
			auto wrong_configuration=configuration;if(rank==0)wrong_configuration[0]=wrong_configuration[0]=='a'?'b':'a';
			Reject(comm,[&] { (void)iga::RestoreMovingCheckpointRuntime(comm,bundle_root,*manifest,layout,wrong_configuration,grid,geometry_options,options); });
			unchanged(accepted,.125,1);
		}
		resumed->Close();
		resumed=iga::RestoreMovingCheckpointRuntime(comm,bundle_root,*manifest,layout,configuration,grid,geometry_options,options);
		const auto file_values=Owned(resumed->CommittedState());double local_error[2]={},global_error[2]={};
		for(std::size_t i=0;i<accepted.size();++i) {
			const double delta=PetscRealPart(file_values[i]-accepted[i]);local_error[0]+=delta*delta;
			local_error[1]+=PetscRealPart(accepted[i])*PetscRealPart(accepted[i]);
		}
		MPI_Allreduce(local_error,global_error,2,MPI_DOUBLE,MPI_SUM,comm);
		const double file_error=std::sqrt(global_error[0])/std::max(1.,std::sqrt(global_error[1]));
		Require(std::isfinite(file_error)&&file_error<1e-8,"file restored field differs beyond gate");
		if(!read_only) {
			Require(file_values==accepted,"same-run file restored owned field differs");
			Require(iga::SerializeMovingConservationCheckpoint(resumed->ConservationDiagnostics())
				==iga::SerializeMovingConservationCheckpoint(committed_conservation),"file restored conservation differs");
		}
		std::cout<<"moving_checkpoint_files rank="<<rank<<" ranks="<<ranks<<" source_ranks="<<source_ranks
			<<" owned_restore_scaled_l2="<<file_error<<" passed\n";
	}


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
	if(rank==0)std::cerr<<"moving_stage second_oracle_begin\n";
	compare(.125,2);check_traction(.125,2);
	const auto second_trial=Owned(runtime.TrialState());
	std::unique_ptr<iga::MovingCutGeometry> resumed_target;
	iga::CollectiveLocalStage(comm,"moving checkpoint continuation geometry",[&] {
		resumed_target=iga::MovingCutGeometry::Build(grid,motion.Evaluate(.25,.125,.25),geometry_options,&resumed->CommittedGeometry());
	});
	resumed->BeginTrial(std::move(resumed_target),2,layers);
	Require(resumed->SolveTrial(),"restored moving next step did not converge");
	const auto continued=Owned(resumed->TrialState());double local[2]={},global[2]={};
	for(std::size_t i=0;i<continued.size();++i) {
		const double delta=PetscRealPart(continued[i]-second_trial[i]);
		local[0]+=delta*delta;local[1]+=PetscRealPart(second_trial[i])*PetscRealPart(second_trial[i]);
	}
	MPI_Allreduce(local,global,2,MPI_DOUBLE,MPI_SUM,comm);
	const double restart_error=std::sqrt(global[0])/std::max(1.,std::sqrt(global[1]));
	Require(std::isfinite(restart_error)&&restart_error<1e-8,"restored next step field differs");
	resumed->Commit();resumed->Close();
	std::cout<<"moving_checkpoint rank="<<rank<<" changing="<<changing<<" owned_roundtrip_exact=1 next_step_scaled_l2="<<restart_error<<" passed\n";

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
	Reject(comm,[&] { (void)runtime.ConservationDiagnostics(); });
	MPI_Bcast(&maximum_conservation_error,1,MPI_DOUBLE,0,comm);
	std::cout<<"immersed_moving_runtime_mpi rank="<<rank<<" changing="<<changing<<" extension_layers="<<layers<<" maximum_conservation_scaled_error="<<maximum_conservation_error<<" maximum_field_scaled_l2="<<maximum_field_error
		<<" patch_retry_scaled_l2="<<maximum_patch_retry_scaled_l2<<" same_trial_capture_exact=1 deferred_abort=1"
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
