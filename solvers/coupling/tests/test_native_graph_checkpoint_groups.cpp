#include "MultidomainRunner.hpp"
#include "CollectiveFailure.hpp"
#include "CheckedText.hpp"
#include <petscsys.h>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv)
{
	int petsc_argc = 1; char* petsc_arguments[] = {argv[0], nullptr}; char** petsc_argv = petsc_arguments;
	PetscInitialize(&petsc_argc, &petsc_argv, nullptr, nullptr);
	iga::CurrentPhaseProfile().EnableFromEnvironment();
	int rank = 0, size = 1; MPI_Comm_rank(PETSC_COMM_WORLD, &rank); MPI_Comm_size(PETSC_COMM_WORLD, &size);
	int status = 0; MPI_Comm group = MPI_COMM_NULL;
	try {
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "checkpoint groups arguments", [&] {
			if (size != 3 || argc != 3) throw std::runtime_error("usage: mpiexec -np 3 native_graph_checkpoint_groups_test ROOT full|save|resume");
			if (std::string(argv[2]) != "full" && std::string(argv[2]) != "save" && std::string(argv[2]) != "resume") throw std::runtime_error("invalid group action");
		});
		const int color = rank == 0 ? 0 : 1;
		MPI_Comm_split(PETSC_COMM_WORLD, color, rank, &group);
		PetscOptionsSetValue(nullptr, "-ksp_rtol", color == 0 ? "1e-9" : "1e-10");
		std::vector<std::string> arguments; std::vector<char*> native_argv;
		iga::CollectiveLocalStage(group, "checkpoint group command", [&] {
			const auto root = std::filesystem::path(argv[1])/("group-"+std::to_string(color)); const std::string action = argv[2];
			const auto bundle = root/(action == "full" ? "full-bundle" : "split-bundle");
			arguments = {"embedded-native-graph", "--graph-case", (root/"fixture").string(), "--output-dir", (root/(action+"-output")).string(), "--checkpoint-dir", bundle.string()};
			if (action == "save") { arguments.push_back("--stop-after-step"); arguments.push_back("3"); }
			if (action == "resume") { arguments.push_back("--restart-dir"); arguments.push_back(bundle.string()); }
			for (auto& argument : arguments) native_argv.push_back(argument.data());
		});
		status = iga::RunMultidomainFlow(static_cast<int>(native_argv.size()), native_argv.data(), group);
		MPI_Comm_free(&group);
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "checkpoint group profiles", [&] {
			iga::CurrentPhaseProfile().Write(std::cout, rank, size, status); iga::FlushCheckedText(std::cout);
		});
	} catch (const std::exception& error) { std::cerr << "rank " << rank << ": " << error.what() << '\n'; status = 1; }
	int global = 0; MPI_Allreduce(&status, &global, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
	if (group != MPI_COMM_NULL) MPI_Comm_free(&group);
	PetscFinalize(); return global;
}
