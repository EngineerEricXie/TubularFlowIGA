// Exercise the native implementation with injected returned errors and real
// PETSc objects. The test retains one reference to inspect cleanup after each run.
#define main LegacyTransportCliMain
#include "../src/iga_transport.cpp"
#undef main

#include <array>
#include <iterator>
#include <sstream>

namespace {

enum class Fault {
	None, Allocation, Insertion, MatrixAssembly, Rows, InitialValues, BoundaryValues,
	Multiply, Copy, Options, Solve, Norm, Scatter, Read, Restore, Destroy
};
Fault fault = Fault::None;
int remaining = 0;
bool observing = false;
std::array<PetscObject, 32> retained{};
std::size_t retained_count = 0;

void RequireLegacy(bool condition, const std::string& message)
{
	if (!condition) throw std::runtime_error(message);
}

bool Hit(Fault candidate)
{
	if (fault != candidate || --remaining > 0) return false;
	fault = Fault::None;
	return true;
}

void Retain(PetscObject object)
{
	if (!observing || !object) return;
	RequireLegacy(retained_count < retained.size(), "too many observed PETSc objects");
	RequireLegacy(PetscObjectReference(object) == 0, "cannot retain object");
	retained[retained_count++] = object;
}

void Release(MPI_Comm communicator)
{
	iga::RequireCollectiveSameInt(communicator, "legacy test observed objects", static_cast<int>(retained_count));
	// KSP/scatter were retained after their data: release parents first so
	// their owned data references are gone when those data objects are checked.
	while (retained_count) {
		const auto object = retained[retained_count-1];
		iga::CollectiveLocalStage(communicator, "legacy test ownership", [&] {
			PetscInt references = 0;
			RequireLegacy(PetscObjectGetReference(object, &references) == 0 && references == 1,
				"native legacy run leaked an ownership reference");
		});
		iga::RequireCollectivePetscSuccess(communicator, "legacy test release", PetscObjectDereference(object));
		retained[--retained_count] = nullptr;
	}
}

std::string ReadLegacy(const fs::path& path)
{
	std::ifstream input(path);
	RequireLegacy(bool(input), "cannot read legacy output");
	return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
}

} // namespace

extern "C" {
PetscErrorCode __real_MatCreateAIJ(MPI_Comm c, PetscInt m, PetscInt n, PetscInt M, PetscInt N, PetscInt d, const PetscInt* dn, PetscInt o, const PetscInt* on, Mat* a);
PetscErrorCode __wrap_MatCreateAIJ(MPI_Comm c, PetscInt m, PetscInt n, PetscInt M, PetscInt N, PetscInt d, const PetscInt* dn, PetscInt o, const PetscInt* on, Mat* a)
{
	const auto status = __real_MatCreateAIJ(c,m,n,M,N,d,dn,o,on,a);
	if (status) return status;
	Retain(reinterpret_cast<PetscObject>(*a));
	return status;
}
PetscErrorCode __real_VecCreateMPI(MPI_Comm c, PetscInt n, PetscInt N, Vec* v);
PetscErrorCode __wrap_VecCreateMPI(MPI_Comm c, PetscInt n, PetscInt N, Vec* v)
{
	const auto status = __real_VecCreateMPI(c,n,N,v);
	if (status) return status;
	Retain(reinterpret_cast<PetscObject>(*v));
	return status;
}
PetscErrorCode __real_KSPCreate(MPI_Comm c, KSP* k);
PetscErrorCode __wrap_KSPCreate(MPI_Comm c, KSP* k)
{
	const auto status = __real_KSPCreate(c,k);
	if (status) return status;
	Retain(reinterpret_cast<PetscObject>(*k));
	return status;
}
PetscErrorCode __real_VecScatterCreateToZero(Vec x, VecScatter* s, Vec* y);
PetscErrorCode __wrap_VecScatterCreateToZero(Vec x, VecScatter* s, Vec* y)
{
	const auto status = __real_VecScatterCreateToZero(x,s,y);
	if (status) return status;
	Retain(reinterpret_cast<PetscObject>(*y)); Retain(reinterpret_cast<PetscObject>(*s));
	return status;
}
PetscErrorCode __real_MatAssemblyEnd(Mat a, MatAssemblyType t);
PetscErrorCode __wrap_MatAssemblyEnd(Mat a, MatAssemblyType t)
{
	const auto status = __real_MatAssemblyEnd(a,t);
	if (status) return status;
	return Hit(Fault::MatrixAssembly) ? PETSC_ERR_USER : status;
}
PetscErrorCode __real_MatZeroRows(Mat a, PetscInt n, const PetscInt* rows, PetscScalar d, Vec x, Vec b);
PetscErrorCode __wrap_MatZeroRows(Mat a, PetscInt n, const PetscInt* rows, PetscScalar d, Vec x, Vec b)
{
	const auto status = __real_MatZeroRows(a,n,rows,d,x,b);
	if (status) return status;
	return Hit(Fault::Rows) ? PETSC_ERR_USER : status;
}
PetscErrorCode __real_MatMult(Mat a, Vec x, Vec y);
PetscErrorCode __wrap_MatMult(Mat a, Vec x, Vec y)
{
	const auto status = __real_MatMult(a,x,y);
	if (status) return status;
	return Hit(Fault::Multiply) ? PETSC_ERR_USER : status;
}
PetscErrorCode __real_VecCopy(Vec x, Vec y);
PetscErrorCode __wrap_VecCopy(Vec x, Vec y)
{
	const auto status = __real_VecCopy(x,y);
	if (status) return status;
	return Hit(Fault::Copy) ? PETSC_ERR_USER : status;
}
PetscErrorCode __real_KSPSetFromOptions(KSP k);
PetscErrorCode __wrap_KSPSetFromOptions(KSP k)
{
	const auto status = __real_KSPSetFromOptions(k);
	if (status) return status;
	return Hit(Fault::Options) ? PETSC_ERR_USER : status;
}
PetscErrorCode __real_KSPSolve(KSP k, Vec b, Vec x);
PetscErrorCode __wrap_KSPSolve(KSP k, Vec b, Vec x)
{
	const auto status = __real_KSPSolve(k,b,x);
	if (status) return status;
	return Hit(Fault::Solve) ? PETSC_ERR_USER : status;
}
PetscErrorCode __real_VecNorm(Vec v, NormType t, PetscReal* n);
PetscErrorCode __wrap_VecNorm(Vec v, NormType t, PetscReal* n)
{
	const auto status = __real_VecNorm(v,t,n);
	if (status) return status;
	return Hit(Fault::Norm) ? PETSC_ERR_USER : status;
}
PetscErrorCode __real_VecScatterEnd(VecScatter s, Vec x, Vec y, InsertMode m, ScatterMode d);
PetscErrorCode __wrap_VecScatterEnd(VecScatter s, Vec x, Vec y, InsertMode m, ScatterMode d)
{
	const auto status = __real_VecScatterEnd(s,x,y,m,d);
	if (status) return status;
	return Hit(Fault::Scatter) ? PETSC_ERR_USER : status;
}
PetscErrorCode __real_KSPDestroy(KSP* k);
PetscErrorCode __wrap_KSPDestroy(KSP* k)
{
	const auto status = __real_KSPDestroy(k);
	if (status) return status;
	return Hit(Fault::Destroy) ? PETSC_ERR_USER : status;
}
PetscErrorCode __real_MatSetValues(Mat a, PetscInt m, const PetscInt* rows, PetscInt n, const PetscInt* columns, const PetscScalar* values, InsertMode mode);
PetscErrorCode __wrap_MatSetValues(Mat a, PetscInt m, const PetscInt* rows, PetscInt n, const PetscInt* columns, const PetscScalar* values, InsertMode mode)
{
	if (Hit(Fault::Allocation)) throw std::bad_alloc();
	if (Hit(Fault::Insertion)) return PETSC_ERR_USER;
	return __real_MatSetValues(a,m,rows,n,columns,values,mode);
}
PetscErrorCode __real_VecSetValues(Vec v, PetscInt n, const PetscInt* rows, const PetscScalar* values, InsertMode mode);
PetscErrorCode __wrap_VecSetValues(Vec v, PetscInt n, const PetscInt* rows, const PetscScalar* values, InsertMode mode)
{
	if (Hit(Fault::InitialValues) || Hit(Fault::BoundaryValues)) return PETSC_ERR_USER;
	return __real_VecSetValues(v,n,rows,values,mode);
}
PetscErrorCode __real_VecGetArrayRead(Vec v, const PetscScalar** a);
PetscErrorCode __wrap_VecGetArrayRead(Vec v, const PetscScalar** a)
{
	if (Hit(Fault::Read)) return PETSC_ERR_USER;
	return __real_VecGetArrayRead(v,a);
}
PetscErrorCode __real_VecRestoreArrayRead(Vec v, const PetscScalar** a);
PetscErrorCode __wrap_VecRestoreArrayRead(Vec v, const PetscScalar** a)
{
	if (Hit(Fault::Restore)) return PETSC_ERR_USER;
	return __real_VecRestoreArrayRead(v,a);
}
} // extern "C"

namespace {

void Exercise(MPI_Comm communicator, const fs::path& source, const fs::path& directory, int seed)
{
	int rank = 0, ranks = 1;
	MPI_Comm_rank(communicator, &rank);
	MPI_Comm_size(communicator, &ranks);
	iga::CollectiveLocalStage(communicator, "legacy test fixture", [&] {
		if (rank != 0) return;
		RequireLegacy(fs::create_directories(directory), "test directory must be new");
		for (const auto* file : {"controlmesh.vtk", "initial_velocityfield.txt", "simulation_parameter.txt"})
			fs::copy_file(source/file, directory/file);
		fs::copy_file(source/(ranks == 1 ? "serial.ntiga" : "group.ntiga"), directory/"database.ntiga");
		if (seed) {
			auto text = ReadLegacy(directory/"simulation_parameter.txt");
			const auto index = text.find("N0bc 1.25");
			RequireLegacy(index != std::string::npos, "unexpected fixture inlet");
			text.replace(index, 9, "N0bc "+std::to_string(1.25+seed));
			std::ofstream output(directory/"simulation_parameter.txt");
			output << text;
			output.close();
			RequireLegacy(bool(output), "cannot change test inlet");
		}
	});
	struct Case { Fault type; const char* name; const char* stage; int occurrence = 1; };
	const std::vector<Case> cases{
		{Fault::Allocation, "allocation", "legacy transport element assembly"},
		{Fault::Insertion, "insertion", "legacy transport element assembly"},
		{Fault::MatrixAssembly, "matrix-assembly", "matrix assembly end"},
		{Fault::Rows, "rows", "legacy transport MatZeroRows left"},
		{Fault::InitialValues, "initial-values", "legacy transport initial values"},
		{Fault::BoundaryValues, "boundary-values", "legacy transport boundary values", 2},
		{Fault::Multiply, "multiply", "legacy transport MatMult"},
		{Fault::Copy, "copy", "legacy transport VecCopy"},
		{Fault::Options, "options", "legacy transport KSPSetFromOptions"},
		{Fault::Solve, "solve", "legacy transport linear solve"},
		{Fault::Norm, "norm", "legacy transport VecNorm"},
		{Fault::Scatter, "scatter", "legacy output scatter end"},
		{Fault::Read, "read", "legacy transport output"},
		{Fault::Restore, "restore", "legacy transport output"},
		{Fault::Destroy, "destroy", "legacy transport KSPDestroy"}
	};
	auto run = [&](const fs::path& path) {
		std::vector<std::string> arguments{"iga_transport", (directory/"database.ntiga").string(),
			directory.string(), "2", path.string()};
		std::vector<char*> argv;
		for (auto& argument : arguments) argv.push_back(argument.data());
		std::ostringstream log;
		const auto original = std::cout.rdbuf(log.rdbuf());
		std::string error;
		observing = true;
		try { RunLegacyTransport(static_cast<int>(argv.size()), argv.data(), communicator); }
		catch (const std::exception& failure) { error = failure.what(); }
		observing = false;
		std::cout.rdbuf(original);
		Release(communicator);
		return std::make_pair(error, log.str());
	};
	const auto baseline = directory/"baseline.txt";
	const auto healthy = run(baseline);
	iga::CollectiveLocalStage(communicator, "legacy test baseline", [&] {
		RequireLegacy(healthy.first.empty(), "baseline failed: "+healthy.first);
	});
	// Different numbers of calls in disjoint groups expose accidental world collectives.
	if (seed == 1) {
		const auto extra = run(directory/"extra.txt");
		RequireLegacy(extra.first.empty(), "independent extra run failed");
	}
	for (const auto& test : cases) {
		const auto path = directory/(std::string(test.name)+".txt");
		const int failed_rank = test.type == Fault::Read || test.type == Fault::Restore ? 0 : ranks-1;
		fault = rank == failed_rank ? test.type : Fault::None;
		remaining = test.occurrence;
		const auto failure = run(path);
		iga::CollectiveLocalStage(communicator, "legacy test rejection", [&] {
			RequireLegacy(fault == Fault::None, "fault was not reached");
			RequireLegacy(failure.first.find(std::string(test.stage)+": rank "+std::to_string(failed_rank)+":") == 0,
				"wrong failure: "+failure.first);
			RequireLegacy(failure.second.find("transport_v2 nodes=") == std::string::npos, "false success summary");
		});
		iga::RequireCollectiveSameText(communicator, "legacy test common diagnostic", failure.first);
		const auto retry = run(path);
		iga::CollectiveLocalStage(communicator, "legacy test retry", [&] {
			RequireLegacy(retry.first.empty(), "retry failed: "+retry.first);
			if (rank == 0) RequireLegacy(ReadLegacy(path) == ReadLegacy(baseline), "retry changed the field");
		});
	}
	// Modify the actual options database after PetscInitialize; some PETSc
	// installations synchronize environment-supplied options during startup.
	if (ranks > 1) {
		if (rank == ranks-1) PetscOptionsSetValue(nullptr, "-legacy_test_option", "different");
		const auto rejected = run(directory/"different-options.txt");
		PetscOptionsClearValue(nullptr, "-legacy_test_option");
		iga::CollectiveLocalStage(communicator, "legacy test option agreement", [&] {
			RequireLegacy(rejected.first.find("PETSc option agreement: rank ") == 0, "different effective options accepted");
		});
		const auto retry = run(directory/"same-options.txt");
		iga::CollectiveLocalStage(communicator, "legacy test options retry", [&] {
			RequireLegacy(retry.first.empty(), "retry after option agreement failed");
			if (rank == 0) RequireLegacy(ReadLegacy(directory/"same-options.txt") == ReadLegacy(baseline),
				"options retry changed the field");
		});
	}
	if (rank == 0) std::cout << "legacy faults=" << cases.size()+(ranks > 1) << " retries=" << cases.size()+(ranks > 1)
		<< " ranks=" << ranks << " seed=" << seed << " ownership=PASS\n";
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, "legacy transport failure test\n");
	int rank = 0, ranks = 1;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
	int status = 0;
	try {
		RequireLegacy(argc == 3 && ranks == 2, "usage: mpiexec -np 2 legacy_transport_failure_test FIXTURE NEW_OUTPUT");
		PetscOptionsSetValue(nullptr, "-ksp_type", "gmres");
		PetscOptionsSetValue(nullptr, "-pc_type", "lu");
		PetscOptionsSetValue(nullptr, "-pc_factor_mat_solver_type", "mumps");
		Exercise(PETSC_COMM_WORLD, argv[1], fs::path(argv[2])/"world", 0);
		MPI_Comm group = MPI_COMM_NULL;
		MPI_Comm_split(PETSC_COMM_WORLD, rank, 0, &group);
		Exercise(group, argv[1], fs::path(argv[2])/("group-"+std::to_string(rank)), rank+1);
		MPI_Comm_free(&group);
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		status = 1;
	}
	int global = 0;
	MPI_Allreduce(&status, &global, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
	PetscFinalize();
	return global;
}
