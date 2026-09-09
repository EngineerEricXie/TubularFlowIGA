#include "ExecutionResources.hpp"
#include <iostream>
#include <vector>

namespace {
void Check(bool value, const char* message)
{
	if (!value) throw std::runtime_error(message);
}

template<class Work>
void Reject(MPI_Comm comm, const char* stage, Work&& work)
{
	std::string message;
	try { work(); } catch (const std::exception& error) { message = error.what(); }
	iga::CollectiveLocalStage(comm,"resource test diagnostic",[&] {
		Check(message.find(stage) != std::string::npos,"missing or wrong resource diagnostic");
	});
	iga::RequireCollectiveSameText(comm,"resource common diagnostic",message);
}

struct Environment {
	std::string name, value;
	bool existed;
	explicit Environment(const char* key) : name(key), existed(std::getenv(key) != nullptr)
	{
		if (existed) value = std::getenv(key);
	}
	~Environment() { if (existed) ::setenv(name.c_str(),value.c_str(),1); else ::unsetenv(name.c_str()); }
};

void Run(MPI_Comm comm)
{
	int rank = 0, ranks = 1, cases = 0;
	MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&ranks);
	std::ostringstream report;
	iga::RequireExecutionResources(comm,&report);
	iga::CollectiveLocalStage(comm,"resource report",[&] {
		if (rank == 0) Check(report.str().find("execution_resources ranks="+std::to_string(ranks)) != std::string::npos,"missing resource summary");
	});
	const std::vector<std::pair<const char*,const char*>> invalid{
		{"OMP_NUM_THREADS","0"},{"OMP_NUM_THREADS","2,"},{"OMP_NUM_THREADS","1,x"},
		{"OMP_NUM_THREADS","2147483648"},{"OMP_NUM_THREADS",""},{"OMP_THREAD_LIMIT","-1"},
		{"OPENBLAS_NUM_THREADS","2x"},{"MKL_NUM_THREADS","-4"},{"BLIS_NUM_THREADS","1,2"}};
	for (const auto& setting : invalid) {
		Environment saved(setting.first);
		if (rank == ranks-1) Check(::setenv(setting.first,setting.second,1)==0,"setenv failed");
		Reject(comm,"execution resource preflight",[&] { iga::RequireExecutionResources(comm); }); ++cases;
	}
	{
		Environment omp("OMP_NUM_THREADS"), blas("OPENBLAS_NUM_THREADS");
		::setenv("OMP_NUM_THREADS"," 2, 3 ",1);
		::setenv("OPENBLAS_NUM_THREADS",rank == 0 ? "0" : "3",1);
		iga::RequireExecutionResources(comm); // valid nested list and heterogeneous BLAS request
	}
	for (int mode = 0; mode < 3; ++mode) {
		Reject(comm,"packed execution test",[&] {
			iga::CollectiveLocalStage(comm,"packed execution test",[&] {
				auto partitions = static_cast<std::uint64_t>(ranks);
				std::uint64_t nodes = 64, fields = 4;
				if (rank == ranks-1) {
					if (mode == 0) ++partitions;
					if (mode == 1) fields = 0;
					if (mode == 2) nodes = std::numeric_limits<std::uint64_t>::max();
				}
				iga::ValidatePackedExecution(partitions,nodes,fields,ranks);
			});
		}); ++cases;
	}
	std::ostringstream broken;
	if (rank == 0) broken.setstate(std::ios::badbit);
	Reject(comm,"execution resource report",[&] { iga::RequireExecutionResources(comm,&broken); }); ++cases;
	Mat matrix = nullptr;
	iga::RequireCollectivePetscSuccess(comm,"resource matrix",MatCreateAIJ(comm,PETSC_DECIDE,PETSC_DECIDE,8,8,1,nullptr,0,nullptr,&matrix));
	Reject(comm,"factor backend availability",[&] { iga::RequireFactorBackend(matrix,comm,"not_a_backend"); }); ++cases;
	Reject(comm,"factor backend preparation",[&] { iga::RequireFactorBackend(matrix,comm,rank == ranks-1 ? "" : "mumps"); }); ++cases;
	if (ranks > 1) {
		Reject(comm,"factor backend agreement",[&] { iga::RequireFactorBackend(matrix,comm,rank == ranks-1 ? "petsc" : "mumps"); }); ++cases;
	}
	iga::RequireFactorBackend(matrix,comm,"mumps");
	KSP solver = nullptr; PC pc = nullptr;
	KSPCreate(comm,&solver); KSPGetPC(solver,&pc); PCSetType(pc,PCLU);
	PCFactorSetMatSolverType(pc,"not_a_backend");
	Reject(comm,"factor backend availability",[&] { iga::RequireKspFactorBackend(solver,matrix,comm); }); ++cases;
	PCFactorSetMatSolverType(pc,"mumps");
	iga::RequireKspFactorBackend(solver,matrix,comm);
	PCSetType(pc,PCJACOBI); iga::RequireKspFactorBackend(solver,matrix,comm);
	KSPDestroy(&solver); MatDestroy(&matrix);
	iga::RequireExecutionResources(comm);
	if (rank == 0) std::cout << "execution_resources ranks=" << ranks << " cases=" << cases << " passed\n";
}
} // namespace

int main(int argc, char** argv)
{
	const bool single = argc == 2 && std::string(argv[1]) == "single";
	if (single) PETSC_MPI_THREAD_REQUIRED = MPI_THREAD_SINGLE;
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int rank = 0, ranks = 1;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank); MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	try {
		Check(ranks == 3,"resource test requires 3 MPI ranks");
		if (single) {
			int provided = -1; MPI_Query_thread(&provided);
			Check(provided == MPI_THREAD_SINGLE,"test MPI did not provide requested SINGLE mode");
#ifdef _OPENMP
			if (omp_get_max_threads() > 1) {
				Reject(PETSC_COMM_WORLD,"execution resource preflight",[&] { iga::RequireExecutionResources(PETSC_COMM_WORLD); });
				if (rank == 0) std::cout << "execution_resources MPI_THREAD_SINGLE rejection passed\n";
			} else
#endif
			{
				iga::RequireExecutionResources(PETSC_COMM_WORLD);
				if (rank == 0) std::cout << "execution_resources MPI_THREAD_SINGLE serial passed\n";
			}
		} else {
			Run(PETSC_COMM_WORLD);
			MPI_Comm group = MPI_COMM_NULL;
			MPI_Comm_split(PETSC_COMM_WORLD,rank == 0 ? 0 : 1,rank,&group);
			Run(group); MPI_Comm_free(&group);
			if (rank == 0) std::cout << "execution_resources all comparisons passed\n";
		}
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		MPI_Abort(PETSC_COMM_WORLD,1); return 1;
	}
	PetscFinalize(); return 0;
}
