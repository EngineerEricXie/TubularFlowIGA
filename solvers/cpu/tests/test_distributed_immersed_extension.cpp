#include "PrescribedSurfaceMotion.hpp"
#include "DistributedImmersedVelocityExtension.hpp"
#include "ImmersedDistributedTransientVolume.hpp"
#include <iostream>
#include <cstring>

namespace {
iga::RawSurfaceSoup Surface(double x)
{
	iga::RawSurfaceSoup s;s.vertices={{{{x,.35,.35}},{{x+.72,.35,.35}},{{x,1.12,.35}},{{x,.35,1.12}}}};
	for(const auto& indices:{std::array<std::int64_t,3>{{0,2,1}},{{0,1,3}},{{0,3,2}},{{1,2,3}}}) {
		iga::RawSurfaceTriangle t;t.indices=indices;t.boundary_id=1;s.triangles.push_back(t);
	}
	return s;
}
void Require(bool good,const char* message) { if(!good)throw std::runtime_error(message); }
}
int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int rank=0,ranks=1;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	int status=0;
	try {
		iga::CubicCartesianGridSpec grid{{{0,0,0}},{{2,2,2}},{{6,6,6}}};
		iga::MovingCutGeometryOptions options;options.volume.max_depth=3;
		options.volume.max_nodes=options.volume.max_leaves=options.volume.max_points=options.volume.max_records=options.volume.max_logical_points=1000000;
		options.volume.max_retained_bytes=100000000;
		iga::PrescribedSurfaceMotion motion({{0.,Surface(.18)},{1.,Surface(.40)}});
		auto old=iga::MovingCutGeometry::Build(grid,motion.Evaluate(0,0,1),options);
		auto target=iga::MovingCutGeometry::Build(grid,motion.Evaluate(1,0,1),options,old.get());
		const auto old_layout=iga::ImmersedActiveLayout::Build(old->Domain(),old->Volume(),old->GeometryIdentitySha256());
		const auto new_layout=iga::ImmersedActiveLayout::Build(target->Domain(),target->Volume(),target->GeometryIdentitySha256());
		Require(old_layout.NodeIds()!=new_layout.NodeIds(),"fixture did not change active nodes");
		std::vector<std::uint64_t> owned,queries,reference_ids;std::vector<double> values,reference_values;
		std::vector<std::array<double,4>> full;
		for(std::size_t row=0;row<old_layout.NodeIds().size();++row) {
			const auto id=old_layout.NodeIds()[row];
			std::array<double,4> value{{std::sin(.07*id),std::cos(.03*id),-0.,std::sin(.11*id)}};
			if(rank==0)full.push_back(value);
			if(row%ranks==static_cast<std::size_t>(rank)) { owned.push_back(id);values.insert(values.end(),value.begin(),value.end()); }
		}
		for(std::size_t row=0;row<new_layout.NodeIds().size();++row)
			if(row%ranks==static_cast<std::size_t>(rank))queries.push_back(new_layout.NodeIds()[row]);
		iga::CollectiveLocalStage(PETSC_COMM_WORLD,"serial extension oracle",[&] {
			if(rank!=0)return;
			iga::ImmersedGlobalFlowState state(0.,0,old_layout,full);
			const auto serial=iga::ImmersedVelocityExtension::Build(*old,old_layout,state,*target,new_layout,3);
			std::vector<double> pressure;for(const auto& value:full)pressure.push_back(value[3]);
			const auto scalar=serial.ExtendScalar(pressure);
			for(std::size_t row=0;row<new_layout.NodeIds().size();++row) {
				reference_ids.push_back(new_layout.NodeIds()[row]);
				for(double value:serial.TargetHistory().Velocities()[row])reference_values.push_back(value);
				reference_values.push_back(scalar.target_values[row]);
			}
		});
		const auto expected=iga::FetchOwnedPointValues(PETSC_COMM_WORLD,reference_ids,reference_values,queries,4);
		iga::DistributedImmersedVelocityExtension extension(PETSC_COMM_WORLD,*old,old_layout,*target,new_layout,3);
		const auto result=extension.Extend(owned,values,queries);
		long double local[2]{},global[2]{};double local_max=0.,maximum=0.;
		for(std::size_t i=0;i<result.size();++i) {
			local[0]+=static_cast<long double>(result[i]-expected[i])*(result[i]-expected[i]);
			local[1]+=static_cast<long double>(expected[i])*expected[i];
			local_max=std::max(local_max,std::abs(result[i]-expected[i]));
			if(std::binary_search(old_layout.NodeIds().begin(),old_layout.NodeIds().end(),queries[i/4]))
				Require(std::memcmp(&result[i],&expected[i],sizeof(double))==0,"anchor bits changed");
		}
		MPI_Allreduce(local,global,2,MPI_LONG_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
		MPI_Allreduce(&local_max,&maximum,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);
		const double relative=std::sqrt(global[0]/global[1]);Require(relative<2.e-10,"sparse extension differs from serial");
		// Runtime state uses contiguous whole-node ownership, unlike the
		// round-robin tuple input above. No full field scatter is constructed.
		const auto node_count=old_layout.NodeIds().size();
		const auto first_node=node_count*static_cast<std::size_t>(rank)/ranks;
		const auto last_node=node_count*static_cast<std::size_t>(rank+1)/ranks;
		Vec committed=nullptr;
		iga::RequireCollectivePetscSuccess(PETSC_COMM_WORLD,"extension fixture vector",VecCreateMPI(PETSC_COMM_WORLD,
			4*(last_node-first_node),old_layout.Rows(),&committed));
		iga::CollectiveLocalStage(PETSC_COMM_WORLD,"extension fixture vector fill",[&] {
			PetscScalar* data=nullptr;Require(VecGetArray(committed,&data)==0,"fixture vector view failed");
			for(std::size_t i=first_node;i<last_node;++i) {
				const auto id=old_layout.NodeIds()[i];
				const std::array<double,4> value{{std::sin(.07*id),std::cos(.03*id),-0.,std::sin(.11*id)}};
				for(int c=0;c<4;++c)data[4*(i-first_node)+c]=value[c];
			}
			Require(VecRestoreArray(committed,&data)==0,"fixture vector restore failed");
		});
		const auto from_vec=extension.ExtendCommitted(committed,old_layout,0.,0,1.,1,1.,queries);
		long double vector_local[2]{},vector_global[2]{};
		for(std::size_t i=0;i<from_vec.size();++i) {
			vector_local[0]+=static_cast<long double>(from_vec[i]-expected[i])*(from_vec[i]-expected[i]);
			vector_local[1]+=static_cast<long double>(expected[i])*expected[i];
			if(std::binary_search(old_layout.NodeIds().begin(),old_layout.NodeIds().end(),queries[i/4]))
				Require(std::memcmp(&from_vec[i],&expected[i],sizeof(double))==0,"Vec source anchor bits changed");
		}
		MPI_Allreduce(vector_local,vector_global,2,MPI_LONG_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
		const double vector_relative=std::sqrt(vector_global[0]/vector_global[1]);
		Require(vector_relative<2.e-10,"Vec extension differs from serial");
		auto reject_source=[&](auto&& operation) {
			int rejected=0,total=0;try{operation();}catch(const std::runtime_error&){rejected=1;}
			MPI_Allreduce(&rejected,&total,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
			Require(total==ranks,"Vec extension rejection was not collective");
		};
		reject_source([&] { (void)extension.ExtendCommitted(rank==0?nullptr:committed,old_layout,0.,0,1.,1,1.,queries); });
		reject_source([&] { (void)extension.ExtendCommitted(committed,new_layout,0.,0,1.,1,1.,queries); });
		reject_source([&] { (void)extension.ExtendCommitted(committed,old_layout,rank==0?.25:0.,0,1.,1,1.,queries); });
		if(ranks>1) {
			reject_source([&] { (void)extension.ExtendCommitted(committed,old_layout,0.,rank==0?1:0,1.,rank==0?2:1,1.,queries); });
			Vec wrong_comm=nullptr,split_node=nullptr;
			iga::RequireCollectivePetscSuccess(PETSC_COMM_WORLD,"extension wrong communicator fixture",VecCreateSeq(PETSC_COMM_SELF,old_layout.Rows(),&wrong_comm));
			reject_source([&] { (void)extension.ExtendCommitted(wrong_comm,old_layout,0.,0,1.,1,1.,queries); });
			iga::RequireCollectivePetscSuccess(PETSC_COMM_WORLD,"extension wrong communicator destroy",VecDestroy(&wrong_comm));
			const PetscInt split_rows=4*(last_node-first_node)+(rank==0?1:0)-(rank==ranks-1?1:0);
			iga::RequireCollectivePetscSuccess(PETSC_COMM_WORLD,"extension split node fixture",VecCreateMPI(PETSC_COMM_WORLD,split_rows,old_layout.Rows(),&split_node));
			reject_source([&] { (void)extension.ExtendCommitted(split_node,old_layout,0.,0,1.,1,1.,queries); });
			iga::RequireCollectivePetscSuccess(PETSC_COMM_WORLD,"extension split node destroy",VecDestroy(&split_node));
		}
		(void)extension.ExtendCommitted(committed,old_layout,0.,0,1.,1,1.,queries);
		iga::CollectiveLocalStage(PETSC_COMM_WORLD,"extension source immutability",[&] {
			iga::PetscReadArray view;view.Acquire(committed);
			for(std::size_t i=first_node;i<last_node;++i) {
				const auto id=old_layout.NodeIds()[i];
				const std::array<double,4> value{{std::sin(.07*id),std::cos(.03*id),-0.,std::sin(.11*id)}};
				for(int c=0;c<4;++c)Require(std::memcmp(&view.Data()[4*(i-first_node)+c],&value[c],sizeof(double))==0,"extension mutated committed Vec");
			}
			view.Restore();
		});
		Vec target_template=nullptr;
		const auto target_nodes=new_layout.NodeIds().size();
		const auto target_first=target_nodes*static_cast<std::size_t>(rank)/ranks;
		const auto target_last=target_nodes*static_cast<std::size_t>(rank+1)/ranks;
		iga::RequireCollectivePetscSuccess(PETSC_COMM_WORLD,"mapped history target vector",VecCreateMPI(PETSC_COMM_WORLD,4*(target_last-target_first),new_layout.Rows(),&target_template));
		iga::RequireCollectivePetscSuccess(PETSC_COMM_WORLD,"target template sentinel",VecSet(target_template,3.));
		std::vector<std::uint64_t> seed_queries;
		for(std::size_t i=target_first;i<target_last;++i)seed_queries.push_back(new_layout.NodeIds()[i]);
		const auto expected_seed=iga::FetchOwnedPointValues(PETSC_COMM_WORLD,reference_ids,reference_values,seed_queries,4);
		const auto seed=extension.BuildOwnedTargetSeed(committed,old_layout,target_template,new_layout,0.,0,1.,1,1.);
		Require(seed.size()==expected_seed.size(),"moving target seed row count differs");
		double seed_local_error=0.,seed_error=0.;
		for(std::size_t i=0;i<seed.size();++i)seed_local_error=std::max(seed_local_error,std::abs(seed[i]-expected_seed[i]));
		MPI_Allreduce(&seed_local_error,&seed_error,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);
		Require(seed_error<2.e-10,"moving target seed differs from serial");
		iga::CollectiveLocalStage(PETSC_COMM_WORLD,"target template unchanged",[&] {
			iga::PetscReadArray view;view.Acquire(target_template);
			for(std::size_t i=0;i<seed.size();++i)Require(view.Data()[i]==3.,"seed construction wrote target Vec");
			view.Restore();
		});
		if(rank==0)std::cout<<"moving_target_seed=passed ranks="<<ranks<<" max_abs="<<seed_error<<'\n';
		std::vector<std::uint64_t> cells;int force_failure_rank=-1;
		for(const auto& cell:target->Domain().Cells()) {
			const bool positive=cell.classification==iga::CellClassification::Inside
				||(cell.classification==iga::CellClassification::Cut&&target->Volume().Cell(cell.id).diagnostics.estimated_reference_volume>0.);
			if(!positive)continue;
			if(force_failure_rank<0)force_failure_rank=cell.id%ranks;
			if(cell.id%ranks==static_cast<std::uint64_t>(rank))cells.push_back(cell.id);
		}
		iga::ImmersedDistributedTransientVolume inputs(PETSC_COMM_WORLD,target->Domain(),target->Volume(),new_layout,target_template,cells);
		const auto& required=inputs.History().RequiredNodeIds();
		const std::vector<std::uint64_t> required_queries(required.begin(),required.end());
		const auto expected_halo=iga::FetchOwnedPointValues(PETSC_COMM_WORLD,reference_ids,reference_values,required_queries,4);
		iga::NavierStokesParameters parameters;parameters.density=1.;parameters.dynamic_viscosity=.1;parameters.dt=1.;
		const iga::NavierStokesBodyForceEvaluator force=[](const std::array<double,3>&) { return std::array<double,3>{{0.,0.,0.}}; };
		auto freeze=[&] { inputs.Freeze(committed,old_layout,0.,0,1.,1,parameters,force,&extension); };
		freeze();
		Require(inputs.History().HasMappedSource() && inputs.History().SourceGeometryIdentity()==old_layout.GeometryIdentity(),"mapped source metadata lost");
		Require(inputs.History().SourceTimeS()==0. && inputs.History().TargetTimeS()==1. && inputs.History().SourceIndex()==0 && inputs.History().TargetIndex()==1,"mapped history clock changed");
		for(std::size_t i=0;i<required.size();++i) {
			const auto expected_provenance=std::binary_search(old_layout.NodeIds().begin(),old_layout.NodeIds().end(),required[i])
				?iga::ImmersedVelocityHistoryProvenance::Committed:iga::ImmersedVelocityHistoryProvenance::Extended;
			Require(inputs.History().RequiredProvenance()[i]==expected_provenance,"mapped provenance differs");
		}
		reject_source(freeze);Require(inputs.History().Active(),"active history destroyed by duplicate freeze");
		inputs.ReleaseTrial();
		const iga::NavierStokesBodyForceEvaluator bad_force=[&](const std::array<double,3>& x) {
			if(rank==force_failure_rank)throw std::runtime_error("injected mapped force failure");
			return force(x);
		};
		reject_source([&] { inputs.Freeze(committed,old_layout,0.,0,1.,1,parameters,bad_force,&extension); });
		Require(!inputs.History().Active() && inputs.FrozenPointCount()==0,"failed mapped freeze published inputs");
		freeze();
		// Later source mutation cannot change the accepted local halo snapshot.
		iga::RequireCollectivePetscSuccess(PETSC_COMM_WORLD,"mapped source mutation fixture",VecSet(committed,0.));
		double assembly_error=0.,global_assembly_error=0.;
		iga::CollectiveLocalStage(PETSC_COMM_WORLD,"mapped volume oracle",[&] {
			for(auto cell:cells) {
				const auto element=target->Domain().Background().MaterializeElement(cell);
				const auto actual_history=inputs.History().Localize(element,1.);
				iga::NavierStokesVelocityHistory history(element.connectivity.size());
				for(std::size_t i=0;i<element.connectivity.size();++i) {
					const auto position=std::lower_bound(required.begin(),required.end(),element.connectivity[i])-required.begin();
					for(int c=0;c<3;++c) {
						history[i][c]=expected_halo[4*position+c];
						Require(std::abs(actual_history[i][c]-history[i][c])<2.e-10,"mapped halo differs after source mutation");
					}
				}
				const std::vector<std::array<double,4>> current(element.connectivity.size());
				const auto actual=inputs.BuildVolume(cell,current);
				const auto expected_system=iga::BuildNavierStokesElementFromLocalVelocityHistory(element,current,history,parameters,
					target->Volume().UsableRule(target->Domain(),cell),force,iga::NavierStokesResolvedMixedForm::Conservative);
				for(std::size_t i=0;i<actual.negative_residual.size();++i)assembly_error=std::max(assembly_error,std::abs(actual.negative_residual[i]-expected_system.negative_residual[i]));
				for(std::size_t i=0;i<actual.jacobian.size();++i)assembly_error=std::max(assembly_error,std::abs(actual.jacobian[i]-expected_system.jacobian[i]));
			}
		});
		MPI_Allreduce(&assembly_error,&global_assembly_error,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);
		Require(global_assembly_error<2.e-10,"mapped volume differs from serial history assembly");
		inputs.Close();
		iga::RequireCollectivePetscSuccess(PETSC_COMM_WORLD,"mapped history target destroy",VecDestroy(&target_template));
		if(rank==0)std::cout<<"mapped_history_volume=passed ranks="<<ranks<<" assembly_max_abs="<<global_assembly_error<<" force_failure_retry=passed\n";
		iga::RequireCollectivePetscSuccess(PETSC_COMM_WORLD,"extension fixture vector destroy",VecDestroy(&committed));
		if(rank==0)std::cout<<"committed_vec_extension=passed ranks="<<ranks<<" relative_l2="<<vector_relative<<" source_immutable=1\n";
		for(int mode=0;mode<3;++mode) {
			auto bad_ids=owned,bad_queries=queries;auto bad_values=values;
			if(rank==0) {
				if(mode==0){bad_ids.erase(bad_ids.begin());bad_values.erase(bad_values.begin(),bad_values.begin()+4);}
				if(mode==1)bad_values[0]=std::numeric_limits<double>::infinity();
				if(mode==2)bad_queries.push_back(std::numeric_limits<std::uint64_t>::max());
			}
			int rejected=0,total=0;
			try{(void)extension.Extend(bad_ids,bad_values,bad_queries);}catch(const std::runtime_error&){rejected=1;}
			MPI_Allreduce(&rejected,&total,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);Require(total==ranks,"extension failure not collective");
		}
		const auto retry=extension.Extend(owned,values,queries);
		double retry_local=0.,retry_max=0.,scale_local=1.,scale=1.;
		for(std::size_t i=0;i<retry.size();++i) {
			retry_local=std::max(retry_local,std::abs(retry[i]-result[i]));scale_local=std::max(scale_local,std::abs(result[i]));
			if(std::binary_search(old_layout.NodeIds().begin(),old_layout.NodeIds().end(),queries[i/4]))
				Require(std::memcmp(&retry[i],&result[i],sizeof(double))==0,"retry anchor bits changed");
		}
		MPI_Allreduce(&retry_local,&retry_max,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);
		MPI_Allreduce(&scale_local,&scale,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);
		Require(retry_max/scale<1.e-12,"extension retry numerical mismatch");
		if(rank==0)std::cout<<"distributed_extension=passed ranks="<<ranks<<" relative_l2="<<relative<<" max_abs="<<maximum<<" retry_relative_max="<<retry_max/scale<<" residual="<<*std::max_element(extension.Residuals().begin(),extension.Residuals().end())<<'\n';
		extension.Close();
		int closed_rejected=0,closed_total=0;
		try{(void)extension.Extend(owned,values,queries);}catch(const std::runtime_error&){closed_rejected=1;}
		MPI_Allreduce(&closed_rejected,&closed_total,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
		Require(closed_total==ranks,"closed extension remained usable");
		iga::PrescribedSurfaceMotion stationary_motion({{0.,Surface(.18)},{1.,Surface(.18)}});
		auto stationary=iga::MovingCutGeometry::Build(grid,stationary_motion.Evaluate(1,0,1),options,old.get());
		const auto stationary_layout=iga::ImmersedActiveLayout::Build(stationary->Domain(),stationary->Volume(),stationary->GeometryIdentitySha256());
		Require(stationary_layout.NodeIds()==old_layout.NodeIds(),"stationary fixture changed active nodes");
		iga::DistributedImmersedVelocityExtension identity_extension(PETSC_COMM_WORLD,*old,old_layout,*stationary,stationary_layout,0);
		Require(identity_extension.OwnedRows()==0,"stationary extension introduced unknowns");
		std::vector<std::uint64_t> root_ids;std::vector<double> root_values;
		if(rank==0) {
			for(auto id:old_layout.NodeIds())root_ids.push_back(id);
			for(const auto& value:full)root_values.insert(root_values.end(),value.begin(),value.end());
		}
		const auto identity_values=identity_extension.Extend(root_ids,root_values,owned);
		Require(identity_values.size()==values.size() && std::memcmp(identity_values.data(),values.data(),values.size()*sizeof(double))==0,"stationary identity extension changed bits");
		identity_extension.Close();
		// The last rank may own only controller/gauge rows and no nodes.
		const auto controlled_layout=iga::ImmersedActiveLayout::Build(old->Domain(),old->Volume(),old->GeometryIdentitySha256(),{17},true);
		iga::DistributedImmersedVelocityExtension controlled(PETSC_COMM_WORLD,*old,controlled_layout,*stationary,stationary_layout,0);
		const PetscInt controlled_rows=(rank==0?controlled_layout.NodeFieldRows():0)+(rank==ranks-1?2:0);
		Vec controlled_state=nullptr;
		iga::RequireCollectivePetscSuccess(PETSC_COMM_WORLD,"controlled extension vector",VecCreateMPI(PETSC_COMM_WORLD,controlled_rows,controlled_layout.Rows(),&controlled_state));
		iga::CollectiveLocalStage(PETSC_COMM_WORLD,"controlled extension fill",[&] {
			PetscScalar* data=nullptr;Require(VecGetArray(controlled_state,&data)==0,"controlled view failed");
			if(rank==0)for(std::size_t i=0;i<root_values.size();++i)data[i]=root_values[i];
			if(rank==ranks-1){data[controlled_rows-2]=1234.5;data[controlled_rows-1]=-5678.9;}
			Require(VecRestoreArray(controlled_state,&data)==0,"controlled restore failed");
		});
		const auto controlled_values=controlled.ExtendCommitted(controlled_state,controlled_layout,0.,0,1.,1,1.,owned);
		Require(controlled_values.size()==values.size() && std::memcmp(controlled_values.data(),values.data(),values.size()*sizeof(double))==0,"controller rows entered node extension");
		iga::CollectiveLocalStage(PETSC_COMM_WORLD,"controller rows unchanged",[&] {
			iga::PetscReadArray view;view.Acquire(controlled_state);
			if(rank==ranks-1)Require(view.Data()[controlled_rows-2]==1234.5 && view.Data()[controlled_rows-1]==-5678.9,"extension mutated controller rows");
			view.Restore();
		});
		const auto seed_layout=iga::ImmersedActiveLayout::Build(stationary->Domain(),stationary->Volume(),stationary->GeometryIdentitySha256(),{17},true);
		iga::DistributedImmersedVelocityExtension seed_extension(PETSC_COMM_WORLD,*old,controlled_layout,*stationary,seed_layout,0);
		Vec seed_template=nullptr;
		iga::RequireCollectivePetscSuccess(PETSC_COMM_WORLD,"controller seed template",VecCreateMPI(PETSC_COMM_WORLD,rank==0?seed_layout.Rows():0,seed_layout.Rows(),&seed_template));
		iga::RequireCollectivePetscSuccess(PETSC_COMM_WORLD,"controller seed sentinel",VecSet(seed_template,-99.));
		auto build_seed=[&] { return seed_extension.BuildOwnedTargetSeed(controlled_state,controlled_layout,seed_template,seed_layout,0.,0,1.,1,1.); };
		const auto controller_seed=build_seed();
		if(rank==0) {
			Require(controller_seed.size()==seed_layout.Rows(),"controller target seed coverage differs");
			for(std::size_t i=0;i<root_values.size();++i)Require(std::memcmp(&controller_seed[i],&root_values[i],sizeof(double))==0,"stationary seed node bits changed");
			Require(controller_seed[seed_layout.ControllerRow(17)]==1234.5,"controller stable ID transfer failed");
			Require(controller_seed[seed_layout.GaugeRow()]==0. && !std::signbit(controller_seed[seed_layout.GaugeRow()]),"target gauge was not canonical zero");
		} else Require(controller_seed.empty(),"empty target owner received seed rows");
		iga::CollectiveLocalStage(PETSC_COMM_WORLD,"inject nonfinite source controller",[&] {
			if(rank!=ranks-1)return;
			PetscScalar* data=nullptr;Require(VecGetArray(controlled_state,&data)==0,"controller write view failed");
			data[controlled_rows-2]=std::numeric_limits<double>::infinity();
			Require(VecRestoreArray(controlled_state,&data)==0,"controller write restore failed");
		});
		reject_source([&] { (void)build_seed(); });
		iga::CollectiveLocalStage(PETSC_COMM_WORLD,"restore source controller",[&] {
			if(rank!=ranks-1)return;
			PetscScalar* data=nullptr;Require(VecGetArray(controlled_state,&data)==0,"controller retry view failed");data[controlled_rows-2]=1234.5;
			Require(VecRestoreArray(controlled_state,&data)==0,"controller retry restore failed");
		});
		Require(build_seed()==controller_seed,"controller seed retry changed values");
		iga::CollectiveLocalStage(PETSC_COMM_WORLD,"controller seed template unchanged",[&] {
			iga::PetscReadArray view;view.Acquire(seed_template);
			for(std::size_t i=0;i<controller_seed.size();++i)Require(view.Data()[i]==-99.,"controller seed wrote target Vec");
			view.Restore();
		});
		seed_extension.Close();
		iga::RequireCollectivePetscSuccess(PETSC_COMM_WORLD,"controller seed template destroy",VecDestroy(&seed_template));
		if(rank==0)std::cout<<"controller_seed=passed ranks="<<ranks<<" owner_transfer=passed gauge_reset=passed retry=passed\n";
		controlled.Close();
		iga::RequireCollectivePetscSuccess(PETSC_COMM_WORLD,"controlled extension destroy",VecDestroy(&controlled_state));
		if(rank==0)std::cout<<"controller_rows_excluded=passed ranks="<<ranks<<'\n';
		if(rank==0)std::cout<<"stationary_extension=passed no_unknowns empty_source_ranks="<<ranks-1<<'\n';
	} catch(const std::exception& error) { std::cerr<<error.what()<<'\n';status=1; }
	int global_status=0;MPI_Allreduce(&status,&global_status,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);
	PetscFinalize();return global_status;
}
