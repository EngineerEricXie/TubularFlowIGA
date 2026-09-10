#include "ParallelVtkOutput.hpp"
#include <sys/resource.h>
#include <csignal>
#include <iostream>

namespace {
void Require(bool value,const char* message)
{
	if(!value)throw std::runtime_error(message);
}
std::string Read(const std::filesystem::path& path)
{
	std::ifstream input(path);Require(bool(input),"cannot read snapshot evidence");
	return std::string((std::istreambuf_iterator<char>(input)),{});
}
iga::VtkPartition Piece(int rank,int ranks,double seed)
{
	iga::VtkPartition piece;piece.point_arrays={{"value",1,{}}};
	if(rank==ranks-1&&ranks>1)return piece;
	piece.grid.points={static_cast<double>(rank),seed,0};piece.grid.connectivity={0};piece.grid.offsets={1};piece.grid.types={1};
	piece.point_ids={rank};piece.cell_ids={rank};piece.point_arrays[0].values={seed+rank};return piece;
}
void Run(MPI_Comm comm,const std::filesystem::path& root,double seed)
{
	int rank=0,ranks=0;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&ranks);
	const auto piece=Piece(rank,ranks,seed);
	for(const auto& mode:std::vector<std::string>{"schema","time","local","write"}) {
		auto candidate=piece;double time=.25;
		if(rank==ranks-1) {
			if(mode=="schema")candidate.point_arrays[0].name="different";
			if(mode=="time")time=.5;
			if(mode=="local")candidate.point_arrays[0].values.push_back(0);
		}
		// With one rank, a different but valid schema/time is not a mismatch.
		if(ranks==1&&(mode=="schema"||mode=="time"))continue;
		struct rlimit saved{};using Handler=void(*)(int);Handler handler=SIG_DFL;
		iga::CollectiveLocalStage(comm,"test file limit setup",[&] {
			if(mode=="write"&&rank==ranks-1) {
				Require(getrlimit(RLIMIT_FSIZE,&saved)==0,"getrlimit failed");
				handler=std::signal(SIGXFSZ,SIG_IGN);Require(handler!=SIG_ERR,"signal setup failed");
				auto limit=saved;limit.rlim_cur=1;Require(setrlimit(RLIMIT_FSIZE,&limit)==0,"setrlimit failed");
			}
		});
		bool rejected=false;std::string message;
		try { iga::WriteParallelVtkSnapshot(comm,root/mode,candidate,time); }
		catch(const std::exception& error) { rejected=true;message=error.what(); }
		iga::CollectiveLocalStage(comm,"test file limit restore",[&] {
			if(mode=="write"&&rank==ranks-1) {
				Require(setrlimit(RLIMIT_FSIZE,&saved)==0,"restore limit failed");
				Require(std::signal(SIGXFSZ,handler)!=SIG_ERR,"restore signal failed");
			}
		});
		iga::CollectiveLocalStage(comm,"test failure agreement",[&] {
			Require(rejected,"expected parallel VTK failure was accepted");
			Require(!std::filesystem::exists(root/mode/"snapshot.pvtu"),"failed snapshot published an index");
			if(mode=="write")Require(message.find("piece write")!=std::string::npos,"did not fail in piece write");
			else Require(!std::filesystem::exists(root/mode),"preflight failure created output");
		});
		iga::WriteParallelVtkSnapshot(comm,root/(mode+"-retry"),piece,.25);
	}
	const auto healthy=root/"healthy";
	iga::WriteParallelVtkSnapshot(comm,healthy,piece,.5);
	std::string saved_index,saved_piece;
	iga::CollectiveLocalStage(comm,"test snapshot capture",[&] {
		saved_index=Read(healthy/"snapshot.pvtu");
		saved_piece=Read(healthy/("rank"+std::to_string(rank)+".vtu"));
	});
	bool rejected=false;
	try { iga::WriteParallelVtkSnapshot(comm,healthy,piece,.75); } catch(const std::exception&) { rejected=true; }
	iga::CollectiveLocalStage(comm,"test immutable snapshot",[&] {
		Require(rejected,"existing snapshot was overwritten");
		Require(Read(healthy/"snapshot.pvtu")==saved_index,"existing snapshot index changed");
		Require(Read(healthy/("rank"+std::to_string(rank)+".vtu"))==saved_piece,"existing snapshot piece changed");
	});
	if(rank==0)std::cout<<"parallel_vtk ranks="<<ranks<<" failures="<<(ranks==1?3:5)<<" retries=passed\n";
}
}
int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int rank=0,ranks=0;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	try {
		Require(argc==2&&ranks==3,"require three ranks and a fresh output directory");
		Run(PETSC_COMM_WORLD,std::filesystem::path(argv[1])/"world",10);
		MPI_Comm group=MPI_COMM_NULL;MPI_Comm_split(PETSC_COMM_WORLD,rank==0?0:1,rank,&group);
		Run(group,std::filesystem::path(argv[1])/(rank==0?"self":"pair"),rank==0?20:30);
		MPI_Comm_free(&group);
	} catch(const std::exception& error) { std::cerr<<"rank "<<rank<<": "<<error.what()<<'\n';MPI_Abort(PETSC_COMM_WORLD,1);return 1; }
	PetscFinalize();return 0;
}
