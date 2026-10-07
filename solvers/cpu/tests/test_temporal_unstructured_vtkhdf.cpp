#include "NativeTetHydraulicVisualization.hpp"
#include "NativeTetSpeciesVisualization.hpp"
#include "TemporalUnstructuredVtkHdf.hpp"
#include <iostream>

static bool fail_flush=false;
extern "C" herr_t __real_H5Fflush(hid_t,H5F_scope_t);
extern "C" herr_t __wrap_H5Fflush(hid_t file,H5F_scope_t scope)
{
	return fail_flush?-1:__real_H5Fflush(file,scope);
}
template<class Function> void Reject(Function&& operation)
{
	bool rejected=false;
	try{operation();}catch(const std::exception&){rejected=true;}
	if(!rejected)throw std::runtime_error("invalid VTKHDF operation was accepted");
}
int main(int argc,char** argv)
{
	if(argc!=2)throw std::invalid_argument("test requires a fresh output directory");
	const std::filesystem::path root=argv[1];
	std::filesystem::create_directories(root);
	iga::NativeTetMesh mesh;
	mesh.points={{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}}};
	mesh.cells={{9007199254740993ULL,{{0,1,2,3}}}};
	std::vector<double> state(34,0.);
	for(std::size_t index=0;index<state.size();++index)state[index]=index*.125;
	auto initial=iga::BuildNativeTetHydraulicVtkPartition(mesh,mesh,state,.0035,0,1);
	for(bool moving:{false,true}){
		const auto stem=moving?"moving":"fixed";
		const auto file=root/(std::string(stem)+".vtkhdf");
		{
			iga::TemporalUnstructuredVtkHdfWriter writer(file,initial,moving);
			writer.Append(0.,initial);
			Reject([&]{writer.Append(0.,initial);});
			auto bad=initial;bad.point_arrays[0].values[0]=std::numeric_limits<double>::quiet_NaN();
			Reject([&]{writer.Append(.1,bad);});
			bad=initial;bad.cell_ids[0]+=1;
			Reject([&]{writer.Append(.1,bad);});
			bad=initial;bad.cell_arrays[0].components=3;
			Reject([&]{writer.Append(.1,bad);});
			if(!moving){bad=initial;bad.grid.points[0]+=.01;Reject([&]{writer.Append(.1,bad);});}
			writer.Close();writer.Close();
			Reject([&]{writer.Append(.1,initial);});
		}
		if(!moving){
			// The original zero-based writer did not record an offset.
			auto legacy=iga::hdf_detail::RequireHandle(H5Fopen(file.string().c_str(),
				H5F_ACC_RDWR,H5P_DEFAULT),H5Fclose,"open legacy fixture");
			auto metadata=iga::hdf_detail::OpenGroup(legacy.get(),"TubularFlowIGA");
			iga::hdf_detail::Require(H5Adelete(metadata.get(),"StepOffset"),"remove offset");
		}
		Reject([&]{iga::TemporalUnstructuredVtkHdfWriter duplicate(file,initial,moving);});
		Reject([&]{iga::TemporalUnstructuredVtkHdfWriter mismatch(file,initial,!moving,4,1);});
		Reject([&]{iga::TemporalUnstructuredVtkHdfWriter mismatch(file,initial,moving,4,2);});
		iga::TemporalUnstructuredVtkHdfWriter resumed(file,initial,moving,4,1);
		iga::WriteVtuPartition(root/(std::string(stem)+"_0.vtu"),initial,0.);
		for(int step=1;step<=2;++step){
			auto current=mesh;
			if(moving)for(auto& point:current.points)point[0]+=.1*step;
			auto values=state;for(auto& value:values)value+=step;
			const auto piece=iga::BuildNativeTetHydraulicVtkPartition(mesh,current,values,.0035,0,1);
			resumed.Append(.1*step,piece);
			iga::WriteVtuPartition(root/(std::string(stem)+"_"+std::to_string(step)+".vtu"),piece,.1*step);
		}
		resumed.Close();
		const auto fresh=root/(std::string(stem)+"_fresh.vtkhdf");
		{
			iga::TemporalUnstructuredVtkHdfWriter writer(fresh,initial,moving,4,0,1);
			writer.Append(.2,initial);writer.Close();
		}
		Reject([&]{iga::TemporalUnstructuredVtkHdfWriter stale(fresh,initial,moving,4,1);});
		Reject([&]{iga::TemporalUnstructuredVtkHdfWriter ahead(fresh,initial,moving,4,3);});
		{
			iga::TemporalUnstructuredVtkHdfWriter writer(fresh,initial,moving,4,2);
			writer.Append(.3,initial);writer.Close();
		}
		iga::TemporalUnstructuredVtkHdfWriter twice(fresh,initial,moving,4,3);
		if(twice.LastTime()!=.3)throw std::runtime_error("fresh series lost resumed time");
		twice.Close();
	}
	const auto species=iga::BuildNativeTetSpeciesVtkPartition(mesh,mesh,{1.,2.,3.,4.},0,1);
	iga::TemporalUnstructuredVtkHdfWriter concentration(root/"species.vtkhdf",species);
	concentration.Append(0.,species);concentration.Close();
	iga::TemporalUnstructuredVtkHdfWriter failure(root/"failure.vtkhdf",initial);
	failure.Append(0.,initial);
	fail_flush=true;Reject([&]{failure.Append(.1,initial);});
	Reject([&]{failure.Close();});fail_flush=false;
	Reject([&]{failure.Append(.1,initial);});failure.Close();
	Reject([&]{iga::TemporalUnstructuredVtkHdfWriter partial(root/"failure.vtkhdf",initial,false,4,1);});
	std::cout<<"temporal_unstructured_vtkhdf_test: PASS\n";
}
