// Strongly coupled native tetrahedral FSI from one combined case file.
#include "NativeTetFsiCase.hpp"

#include <petscksp.h>

#include <filesystem>
#include <iostream>
#include <stdexcept>

int main(int argc,char** argv)
{
	if(PetscInitialize(&argc,&argv,nullptr,nullptr))return 1;
	int rank=0,ranks=1,status=0;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	try{
		// PETSc options follow the two positional arguments.
		if(argc<3||argv[1][0]=='-'||argv[2][0]=='-')throw std::invalid_argument(
			"usage: native_tet_fsi_steps case.json output_dir [PETSc options]");
		const auto options=iga::ReadNativeTetFsiSolverOptions(ranks);
		iga::NativeTetFsiSession session(iga::ReadNativeTetFsiCase(argv[1]),options);
		const std::filesystem::path output(argv[2]);
		const int steps=session.Scenario().steps;
		for(int step=1;step<=steps;++step)session.RunStep(step,output,rank,ranks);
		session.WriteSeries(output);
		if(rank==0)std::cout<<"native_fsi_complete accepted_steps="<<steps<<'\n';
	}catch(const std::exception& error){
		if(rank==0)std::cerr<<"native_tet_fsi_steps: "<<error.what()<<'\n';
		status=1;
	}
	PetscFinalize();return status;
}
