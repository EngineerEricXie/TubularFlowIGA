#include "CollectivePetscOptions.hpp"
#include <iostream>
#include <vector>

namespace {

struct Options {
	Options() { iga::petsc_options_detail::Check(PetscOptionsCreate(&value), "PetscOptionsCreate"); }
	~Options() { PetscOptionsDestroy(&value); }
	Options(const Options&) = delete;
	Options& operator=(const Options&) = delete;
	void Set(const char* name, const char* text)
	{
		iga::petsc_options_detail::Check(PetscOptionsSetValue(value, name, text), "PetscOptionsSetValue");
	}
	PetscOptions value = nullptr;
};

void Require(bool condition, const char* message)
{
	if (!condition) throw std::runtime_error(message);
}

PetscInt Unused(PetscOptions options)
{
	PetscInt count = 0;
	iga::petsc_options_detail::Check(PetscOptionsAllUsed(options, &count), "PetscOptionsAllUsed");
	return count;
}

void Reject(PetscOptions options, MPI_Comm communicator)
{
	int rejected = 0, all_rejected = 0;
	try { iga::RequireCollectivePetscOptions(communicator, options); }
	catch (const std::runtime_error& error) {
		rejected = std::string(error.what()).find("PETSc option agreement:") == 0;
	}
	MPI_Allreduce(&rejected, &all_rejected, 1, MPI_INT, MPI_MIN, communicator);
	Require(all_rejected == 1, "option disagreement did not reach all ranks");
}

std::string Display(PetscOptions options)
{
	iga::petsc_options_detail::Buffers buffers(options);
	iga::petsc_options_detail::Check(PetscOptionsGetAll(options, &buffers.all), "PetscOptionsGetAll");
	return buffers.all ? buffers.all : "";
}

void CheckOptions(MPI_Comm communicator)
{
	int rank = 0, ranks = 1;
	MPI_Comm_rank(communicator, &rank);
	MPI_Comm_size(communicator, &ranks);
	Options options;
	const std::vector<std::pair<const char*, const char*>> entries{
		{"-ksp_type", "gmres"}, {"-hpc_text", "spaces \"quotes\"\n-ksp_type cg "},
		{"-hpc_empty", ""}, {"-hpc_flag", nullptr}, {"-hpc_negative", "-1e-5"}};
	if (rank%2) for (auto entry = entries.rbegin(); entry != entries.rend(); ++entry)
		options.Set(entry->first, entry->second);
	else for (const auto& entry : entries) options.Set(entry.first, entry.second);
	if (rank == 0) {
		char value[128]{};
		PetscOptionsGetString(options.value, nullptr, "-ksp_type", value, sizeof(value), nullptr);
		PetscOptionsGetString(options.value, nullptr, "-hpc_text", value, sizeof(value), nullptr);
		PetscBool flag = PETSC_FALSE;
		PetscOptionsGetBool(options.value, nullptr, "-hpc_flag", &flag, nullptr);
	}
	const auto before = Unused(options.value);
	const auto display = Display(options.value);
	iga::RequireCollectivePetscOptions(communicator, options.value);
	Require(Unused(options.value) == before && Display(options.value) == display,
		"option inspection changed values or unused tracking");
	if (ranks > 1) {
		if (rank == ranks-1) options.Set("-ksp_type", "cg");
		Reject(options.value, communicator);
		options.Set("-ksp_type", "gmres");
		if (rank == ranks-1) options.Set("-fieldsplit_0_pc_type", "hypre");
		Reject(options.value, communicator);
		PetscOptionsClearValue(options.value, "-fieldsplit_0_pc_type");
	}
	const auto local_path = "rank "+std::to_string(rank)+" /path -ksp_type cg";
	options.Set("--graph-case", local_path.c_str());
	iga::RequireCollectivePetscOptions(communicator, options.value, {"--graph-case"});
	Require(Display(options.value).find(local_path) != std::string::npos,
		"excluding an application argument mutated the real options");
	Options one, two;
	one.Set("-aaa", "value -zzz another");
	two.Set("-aaa", "value"); two.Set("-zzz", "another");
	Require(Display(one.value) == Display(two.value), "fixture no longer produces a display collision");
	Require(iga::CapturePetscOptions(one.value) != iga::CapturePetscOptions(two.value),
		"different options with identical display strings compared equal");
	if (ranks > 1) Reject(rank == ranks-1 ? one.value : two.value, communicator);
	iga::RequireCollectivePetscOptions(communicator, options.value, {"--graph-case"});
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int rank = 0, ranks = 1;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
	try {
		Require(ranks == 3, "collective_petsc_options_test requires three ranks");
		CheckOptions(PETSC_COMM_WORLD);
		MPI_Comm group = MPI_COMM_NULL;
		MPI_Comm_split(PETSC_COMM_WORLD, rank == 0 ? 0 : 1, rank, &group);
		CheckOptions(group);
		MPI_Comm_free(&group);
		std::cout << "PETSc options rank=" << rank << " passed\n";
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		MPI_Abort(PETSC_COMM_WORLD, 2);
	}
	PetscFinalize();
	return 0;
}
