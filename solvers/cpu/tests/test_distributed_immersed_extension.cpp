#include "PrescribedSurfaceMotion.hpp"
#include "DistributedImmersedVelocityExtension.hpp"
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
		if(rank==0)std::cout<<"stationary_extension=passed no_unknowns empty_source_ranks="<<ranks-1<<'\n';
	} catch(const std::exception& error) { std::cerr<<error.what()<<'\n';status=1; }
	int global_status=0;MPI_Allreduce(&status,&global_status,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);
	PetscFinalize();return global_status;
}
