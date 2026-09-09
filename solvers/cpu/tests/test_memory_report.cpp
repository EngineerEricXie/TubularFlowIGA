#include "MemoryReport.hpp"

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>

namespace {
template <class Function>
void Reject(MPI_Comm communicator, Function&& work, const char* diagnostic)
{
	std::string message;
	try { work(); }
	catch (const std::exception& error) { message = error.what(); }
	iga::CollectiveLocalStage(communicator, "memory report test rejection", [&] {
		if (message.find(diagnostic) == std::string::npos)
			throw std::runtime_error("missing coordinated rejection: "+message);
	});
	iga::RequireCollectiveSameText(communicator, "memory report test diagnostic agreement", message);
}

void Exercise(MPI_Comm communicator, const std::filesystem::path& path)
{
	int rank = 0;
	MPI_Comm_rank(communicator, &rank);
	std::unique_ptr<iga::DistributedMemoryRecorder> recorder;
	const auto create = [&](const std::filesystem::path& output) {
		iga::CollectiveLocalStage(communicator, "memory report test setup", [&] {
			recorder = std::make_unique<iga::DistributedMemoryRecorder>(communicator, output);
		});
	};
	create({});
	recorder->Record("disabled");
	recorder->Close();
	recorder->Close();
	create(path.string()+".healthy");
	recorder->Record("before_close");
	recorder->Close();
	recorder->Close();
	Reject(communicator, [&] { recorder->Record("after_close"); }, "memory report is closed");
	create(path);
	recorder->Record("before_failed_close");
	iga::CollectiveLocalStage(communicator, "memory report test injection setup", [&] {
		if (setenv("TUBULARFLOWIGA_TEST_TEXT_CLOSE", path.c_str(), 1))
			throw std::runtime_error("cannot select close injection");
	});
	Reject(communicator, [&] { recorder->Close(); }, "memory report close: rank 0: cannot close memory report");
	Reject(communicator, [&] { recorder->Close(); }, "memory report close: rank 0: cannot close memory report");
	Reject(communicator, [&] { recorder->Record("after_failed_close"); }, "memory report is closed");
	iga::CollectiveLocalStage(communicator, "memory report test injection reset", [&] {
		if (unsetenv("TUBULARFLOWIGA_TEST_TEXT_CLOSE")) throw std::runtime_error("cannot reset injection");
	});
	create(path.string()+".retry");
	recorder->Record("new_recorder");
	recorder->Close();
	if (rank == 0) std::cout << "memory report close/repeat/reject/retry passed: " << path << '\n';
}
}

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int rank = 0, size = 0, status = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	MPI_Comm communicator = MPI_COMM_NULL;
	try {
		if (argc != 3 || size != 3) throw std::runtime_error("use three ranks: test_memory_report OUTPUT_DIR world|split");
		const int group = std::string(argv[2]) == "split" && rank > 0 ? 1 : 0;
		MPI_Comm_split(PETSC_COMM_WORLD, group, rank, &communicator);
		std::filesystem::path path;
		iga::CollectiveLocalStage(communicator, "memory report test path", [&] {
			path = std::filesystem::absolute(argv[1])/std::to_string(group);
			std::filesystem::create_directories(path.parent_path());
		});
		Exercise(communicator, path);
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		status = 1;
	}
	if (communicator != MPI_COMM_NULL) MPI_Comm_free(&communicator);
	int global = 0;
	MPI_Allreduce(&status, &global, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
	PetscFinalize();
	return global;
}
