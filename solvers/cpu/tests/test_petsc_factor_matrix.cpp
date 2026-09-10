#include "PetscSolverOptions.hpp"
#include "ExecutionResources.hpp"
#include <iostream>

namespace {
void Check(PetscErrorCode code)
{
	if (code) throw std::runtime_error("PETSc factor matrix test operation failed");
}
struct Options {
	PetscOptions value = nullptr;
	Options()
	{
		Check(PetscOptionsCreate(&value));
	}
	~Options()
	{
		PetscOptionsDestroy(&value);
	}
	void Set(const char* key, const char* value_text)
	{
		Check(PetscOptionsSetValue(value,key,value_text));
	}
};
struct Objects {
	Mat matrix = nullptr; KSP solver = nullptr; Vec rhs = nullptr, solution = nullptr;
	~Objects()
	{
		KSPDestroy(&solver); VecDestroy(&rhs); VecDestroy(&solution); MatDestroy(&matrix);
	}
};
// A known SPD problem tests registered AIJ factor interfaces and their guards.
void Run(MPI_Comm comm)
{
	int rank = 0, size = 0; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&size);
	for (const char* backend : {"petsc","mumps","superlu_dist","superlu","umfpack","unavailable_backend"}) {
		for (const char* type : {PCLU,PCCHOLESKY}) {
			Options source; source.Set("-ksp_type","preonly"); source.Set("-pc_type",type);
			source.Set("-pc_factor_mat_solver_type",backend); source.Set("-ksp_view",nullptr);
			iga::PetscSolverOptions options(comm,"factor_",source.value);
			Objects objects;
			Check(MatCreateAIJ(comm,2,2,2*size,2*size,2,nullptr,2,nullptr,&objects.matrix));
			PetscInt begin = 0, end = 0; Check(MatGetOwnershipRange(objects.matrix,&begin,&end));
			for (PetscInt row = begin; row < end; ++row) {
				Check(MatSetValue(objects.matrix,row,row,4.,INSERT_VALUES));
				if (row) Check(MatSetValue(objects.matrix,row,row-1,-1.,INSERT_VALUES));
				if (row+1<2*size) Check(MatSetValue(objects.matrix,row,row+1,-1.,INSERT_VALUES));
			}
			Check(MatAssemblyBegin(objects.matrix,MAT_FINAL_ASSEMBLY)); Check(MatAssemblyEnd(objects.matrix,MAT_FINAL_ASSEMBLY));
			Check(MatSetOption(objects.matrix,MAT_SYMMETRIC,PETSC_TRUE)); Check(MatSetOption(objects.matrix,MAT_SPD,PETSC_TRUE));
			PetscBool available = PETSC_FALSE;
			options.Call("factor capability",[&] { return MatGetFactorAvailable(objects.matrix,backend,std::string(type)==PCLU ? MAT_FACTOR_LU : MAT_FACTOR_CHOLESKY,&available); });
			iga::RequireCollectiveSameText(comm,"factor capability agreement",available ? "available" : "unavailable");
			Check(KSPCreate(comm,&objects.solver)); options.Attach(objects.solver);
			bool rejected = false;
			try {
				options.Call("factor options",[&] { return KSPSetFromOptions(objects.solver); });
				iga::RequireKspFactorBackend(objects.solver,objects.matrix,comm);
			} catch (const std::exception&) { rejected = true; }
			if (rejected == static_cast<bool>(available)) throw std::runtime_error("capability and preflight disagree");
			PetscReal norm = 0.; PetscInt iterations = 0; int reason = 0;
			if (available) {
				Check(MatCreateVecs(objects.matrix,&objects.solution,&objects.rhs)); Check(VecSet(objects.solution,1.));
				Check(MatMult(objects.matrix,objects.solution,objects.rhs)); Check(VecSet(objects.solution,0.));
				Check(KSPSetOperators(objects.solver,objects.matrix,objects.matrix));
				options.Call("factor solve",[&] { return KSPSolve(objects.solver,objects.rhs,objects.solution); });
				const auto configuration = iga::CaptureKspConfiguration(objects.solver);
				iterations = configuration.last_iterations; reason = static_cast<int>(configuration.last_reason);
				if (configuration.last_reason<=0 || configuration.factor_backend!=backend) throw std::runtime_error("factor solve configuration differs");
				Check(VecShift(objects.solution,-1.)); Check(VecNorm(objects.solution,NORM_2,&norm));
				if (norm>1e-10) throw std::runtime_error("factor field differs from exact solution");
			}
			options.RecordUsed();
			if (rank==0) { std::cout << "factor_matrix ranks=" << size << " backend=" << backend << " pc=" << type << " available=" << int(available) << " rejected=" << rejected << " iterations=" << iterations << " reason=" << reason
				<< " error=";
				if (available) std::cout << norm; else std::cout << "not_solved";
				std::cout << std::endl;
			}
		}
	}
}
}
int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr); int status = 0,rank = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank); MPI_Comm group = MPI_COMM_NULL;
	try { Run(PETSC_COMM_WORLD); MPI_Comm_split(PETSC_COMM_WORLD,rank==0 ? 0 : 1,rank,&group); Run(group); }
	catch (const std::exception& error) { std::cerr << error.what() << '\n'; status = 1; }
	if (group!=MPI_COMM_NULL) MPI_Comm_free(&group);
	int global = 0; MPI_Allreduce(&status,&global,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);
	PetscFinalize(); return global;
}
