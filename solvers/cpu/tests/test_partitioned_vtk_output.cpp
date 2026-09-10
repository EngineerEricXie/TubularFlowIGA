#include "PartitionedVtkOutput.hpp"
#include <iostream>

namespace {
void Require(bool value,const char* message)
{
	if(!value)throw std::runtime_error(message);
}
template<class F> void Reject(F function)
{
	bool rejected=false;try { function(); } catch(const std::invalid_argument&) { rejected=true; }
	Require(rejected,"invalid partition was accepted");
}
iga::VtkPartition Piece(int rank,double time)
{
	iga::VtkPartition piece;
	piece.point_arrays={{"velocity",3,{}},{"temperature & tracer",1,{}}};
	piece.cell_arrays={{"owner",1,{}}};
	if(rank==2)return piece;
	piece.grid.points={0,0,0,1,0,0,0,1,0,0,0,rank==0?1.:-1.};
	piece.grid.connectivity=rank==0?std::vector<std::int64_t>{0,1,2,3}:std::vector<std::int64_t>{0,2,1,3};
	piece.grid.offsets={4};piece.grid.types={10};
	constexpr std::int64_t base=9007199254740993LL;
	piece.point_ids={base,base+1,base+2,base+3+rank};piece.cell_ids={100+rank};
	for(std::size_t point=0;point<4;++point) {
		double sum=time;
		for(unsigned q=0;q<3;++q) {
			const double value=piece.grid.points[3*point+q];sum+=value;
			piece.point_arrays[0].values.push_back(value+time);
		}
		piece.point_arrays[1].values.push_back(sum);
	}
	piece.cell_arrays[0].values={static_cast<double>(rank)};return piece;
}
}
int main(int argc,char** argv)
{
	try {
		if(argc!=2)throw std::invalid_argument("usage: partitioned_vtk_output_test OUTPUT");
		const std::filesystem::path root(argv[1]);
		std::vector<std::pair<double,std::filesystem::path>> snapshots;
		for(int step=0;step<2;++step) {
			const double time=.25*(step+1);std::vector<std::filesystem::path> names;
			for(int rank=0;rank<3;++rank) {
				const auto name="step"+std::to_string(step)+".rank"+std::to_string(rank)+".vtu";
				iga::WriteVtuPartition(root/name,Piece(rank,time),time);names.emplace_back(name);
			}
			const auto index=root/("step"+std::to_string(step)+".pvtu");
			iga::WritePvtu(index,names,{{"velocity",3},{"temperature & tracer",1}},{{"owner",1}});
			snapshots.push_back({time,index});
		}
		iga::WritePvd(root/"series.pvd",snapshots);
		const auto valid=Piece(0,.25);
		Reject([&]{auto bad=valid;bad.point_ids[1]=bad.point_ids[0];iga::ValidateVtkPartition(bad,0);});
		Reject([&]{auto bad=valid;bad.cell_ids[0]=-1;iga::ValidateVtkPartition(bad,0);});
		Reject([&]{auto bad=valid;bad.grid.connectivity[0]=4;iga::ValidateVtkPartition(bad,0);});
		Reject([&]{auto bad=valid;bad.grid.offsets[0]=3;iga::ValidateVtkPartition(bad,0);});
		Reject([&]{auto bad=valid;bad.grid.types[0]=256;iga::ValidateVtkPartition(bad,0);});
		Reject([&]{auto bad=valid;bad.point_arrays[0].values.pop_back();iga::ValidateVtkPartition(bad,0);});
		Reject([&]{auto bad=valid;bad.point_arrays[0].name="GlobalPointIds";iga::ValidateVtkPartition(bad,0);});
		Reject([&]{auto bad=valid;bad.point_arrays[1].name="velocity";iga::ValidateVtkPartition(bad,0);});
		Reject([&]{iga::ValidateVtkPartition(valid,std::numeric_limits<double>::infinity());});
		Reject([&]{iga::WritePvtu(root/"bad.pvtu",{"../escape.vtu"},{},{});});
		Reject([&]{iga::WritePvtu(root/"bad.pvtu",{"same.vtu","same.vtu"},{},{});});
		Require(!std::filesystem::exists(root/"bad.pvtu"),"invalid index was opened before validation");
		std::cout<<"partitioned_vtk_output_test: PASS two steps, shared Int64 ids, empty piece, 11 rejections\n";
	} catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
}
