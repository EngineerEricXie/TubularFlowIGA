#include "BoundarySupport.hpp"
#include "TransientFlowRuntime.hpp"
#include <atomic>
#include <filesystem>
#include <thread>

namespace {
void Require(bool value, const char* message)
{
	iga::CollectiveLocalStage(PETSC_COMM_WORLD, "parallel flow test", [&] {
		if (!value) throw std::runtime_error(message);
	});
}

template<class Work> void Reject(Work work, const char* expected)
{
	std::string diagnostic;
	try { work(); } catch (const std::exception& error) { diagnostic = error.what(); }
	Require(diagnostic.find(expected) != std::string::npos, "expected collective failure was absent or incorrect");
	iga::RequireCollectiveSameText(PETSC_COMM_WORLD, "parallel flow common failure", diagnostic);
}

std::vector<double> Gather(Vec vector)
{
	Vec all = nullptr;
	VecScatter scatter = nullptr;
	VecScatterCreateToAll(vector, &scatter, &all);
	VecScatterBegin(scatter, vector, all, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(scatter, vector, all, INSERT_VALUES, SCATTER_FORWARD);
	PetscInt size = 0;
	VecGetSize(all, &size);
	iga::PetscReadArray view;
	view.Acquire(all);
	std::vector<double> result(static_cast<std::size_t>(size));
	for (PetscInt i = 0; i < size; ++i) result[i] = PetscRealPart(view.Data()[i]);
	view.Restore();
	VecScatterDestroy(&scatter);
	VecDestroy(&all);
	return result;
}

double Compare(const std::vector<double>& reference, const std::vector<double>& candidate, bool pressure)
{
	Require(reference.size() == candidate.size(), "state dimensions differ");
	double error = 0, scale = 0;
	for (std::size_t i = 0; i < reference.size(); ++i) if ((i % 4 == 3) == pressure) {
		error = std::hypot(error, candidate[i] - reference[i]);
		scale = std::hypot(scale, reference[i]);
	}
	Require(std::isfinite(error) && std::isfinite(scale) && (scale > 0 ? error / scale <= 1e-6 : error <= 1e-12),
		"thread configurations differ beyond the predeclared field gate");
	return scale > 0 ? error / scale : error;
}
}

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int status = 0, rank = 0, ranks = 1;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
	iga::CurrentPhaseProfile().EnableFromEnvironment();
	std::cout << std::unitbuf;
	try {
		Require(argc >= 3 && ranks >= 2, "require DATABASE CASE and at least two ranks");
		const std::filesystem::path case_dir = argv[2];
		iga::Database database(argv[1]);
		const auto mesh = iga::ReadLabeledHexMesh((case_dir / "controlmesh.vtk").string(), database.header().nodes, database.header().elements);
		const auto velocity = iga::ReadVelocity((case_dir / "initial_velocityfield.txt").string(), database.header().nodes);
		const auto configuration = iga::ReadSimulationConfiguration((case_dir / "simulation_config.json").string());
		const auto& definition = iga::FirstNavierStokesSystem(configuration);
		const auto boundaries = iga::ResolveFlowBoundaries(configuration, definition, mesh.labels, velocity);
		const auto wall = iga::WallTraceBasis(database, mesh);
		const auto caller = std::this_thread::get_id();
		std::cout << "parallel_flow_stage=resources rank=" << rank << '\n';
		::setenv("IGA_ASSEMBLY_THREADS", rank == 1 ? "0" : "1", 1);
		Reject([&] {
			iga::TransientFlowRuntime invalid(database, PETSC_COMM_WORLD, true, false,
				{definition.density, definition.viscosity, 0}, boundaries, mesh.labels, velocity, wall, {});
		}, "flow assembly resources: rank 1:");
		std::vector<double> reference;
		for (const int threads : {1, 2}) {
			::setenv("IGA_ASSEMBLY_THREADS", std::to_string(threads).c_str(), 1);
			::setenv("IGA_ASSEMBLY_BATCH_SIZE", "8", 1);
			iga::TransientFlowRuntime flow(database, PETSC_COMM_WORLD, true, false,
				{definition.density, definition.viscosity, 0}, boundaries, mesh.labels, velocity, wall, {});
			Require(flow.Elements().size() > 13, "fixture needs at least fourteen required elements per rank");
			flow.Assembler().ValidateOwnership();
			std::cout << "parallel_flow_stage=initialize rank=" << rank << " threads=" << threads << '\n';
			flow.InitializeState(configuration);
			const auto initial = Gather(flow.State());
			std::atomic<std::size_t> worker_calls{0}, completed{0};
			auto observe = [&](std::size_t) {
				++completed;
				if (std::this_thread::get_id() != caller) ++worker_calls;
			};
			const auto fail = [&](std::size_t index) {
				observe(index);
				if (rank == 1 && (index == 10 || index == 12))
					throw std::runtime_error("volume index " + std::to_string(index));
			};
			// Fail before a complete sparse pattern has ever been assembled.
			std::cout << "parallel_flow_stage=first_injection rank=" << rank << " threads=" << threads << '\n';
			flow.SetVolumeProbeForTesting(fail);
			flow.BeginStep(0, 0.0, 12, 1e-8, 1e-12, 1e-6);
			flow.SetTrialBoundaryConfiguration(configuration);
			Reject([&] { flow.SolveTrial(); }, "flow element assembly: rank 1: volume index 10");
			Require(rank != 1 || completed == 16, "failed worker batch did not join or stopped at wrong boundary");
			Reject([&] { flow.PrepareCommitStep(); }, "requires a successful 3D flow trial solve");
			flow.RollbackTrial();
			Require(Gather(flow.State()) == initial, "worker failure changed committed state");
			completed = 0; worker_calls = 0;
			flow.SetVolumeProbeForTesting(observe);
			flow.SetTrialBoundaryConfiguration(configuration);
			flow.SolveTrial(); flow.CommitStep();
			const auto& batch = flow.LastVolumeBatchStatistics();
			Require(batch.items == flow.Elements().size() && batch.maximum_resident_items <= 8 && batch.maximum_team_size == threads,
				"actual team, required-element coverage or batch bound differs");
			Require(threads == 1 ? worker_calls == 0 : worker_calls > 0, "actual volume worker participation differs");
			const auto solved = Gather(flow.State());
			if (threads == 1) reference = solved;
			const double velocity_error = Compare(reference, solved, false), pressure_error = Compare(reference, solved, true);
			// Also exercise a later failure after the full pattern exists.
			std::cout << "parallel_flow_stage=later_injection rank=" << rank << " threads=" << threads << '\n';
			flow.SetVolumeProbeForTesting(fail);
			flow.BeginStep(1, 0.0, 12, 1e-8, 1e-12, 1e-6);
			flow.SetTrialBoundaryConfiguration(configuration);
			Reject([&] { flow.SolveTrial(); }, "flow element assembly: rank 1: volume index 10");
			flow.AbortStep();
			Require(Gather(flow.State()) == solved, "abort changed committed state");
			flow.SetVolumeProbeForTesting({});
			flow.Advance(configuration, 1, 0.0, 12, 1e-8, 1e-12, 1e-6);
			const auto retried = Gather(flow.State());
			Compare(solved, retried, false); Compare(solved, retried, true);
			if (rank == 0) std::cout << "parallel_body_fitted_flow ranks=" << ranks << " threads=" << threads
				<< " first_assembly_failure_retry=1 abort_new_step=1 velocity_l2=" << velocity_error
				<< " pressure_l2=" << pressure_error << " passed\n";
		}
	} catch (const std::exception& error) { std::cerr << error.what() << '\n'; status = 1; }
	iga::CurrentPhaseProfile().Write(std::cout, rank, ranks, status);
	PetscFinalize();
	return status;
}
