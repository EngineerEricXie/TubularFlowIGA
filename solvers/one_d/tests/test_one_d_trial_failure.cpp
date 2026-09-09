#include "OneDRuntime.hpp"
#include "OneDCheckpoint.hpp"
#include "CollectiveFailure.hpp"

#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>

namespace {

void Require(bool condition, const char* message)
{
	if (!condition) throw std::runtime_error(message);
}

std::vector<double> Packed(const iga::OneDFlowRuntime& runtime)
{
	return iga::PackOneDCheckpointState(runtime.FlowState(), runtime.Transports(), runtime.Network());
}

void CheckGroup(MPI_Comm group, const std::filesystem::path& directory, int repetitions)
{
	int rank = 0, ranks = 1;
	MPI_Comm_rank(group, &rank);
	MPI_Comm_size(group, &ranks);
	const int target = ranks > 1 ? 1 : 0;
	std::ifstream file(directory/"simulation_config.json");
	std::ostringstream text;
	text << file.rdbuf();
	auto configuration = iga::ParseOneDConfiguration(text.str());
	const auto inlet = iga::ResolveOneDInlet(configuration);
	const auto flow = configuration.flow_systems.front();
	const auto network = iga::ReadOneDNetwork(directory/configuration.geometry.file,
		configuration.geometry.length_scale_to_m, flow.discretization.cells_per_segment,
		flow.dynamic_viscosity, configuration.geometry.root_node_id);
	const double seed = iga::EvaluateOneDInlet(configuration, inlet, directory, 0.0, network.segments.front().area0);
	const double dt = 3*configuration.time.dt;
	for (int repetition = 0; repetition < repetitions; ++repetition) {
		for (const std::string mode : {"1d trial preparation", "1d trial accounting",
			"1d substep preparation", "1d substep transport", "1d trial completion",
			"1d local flow solve", "missing-inlet", "missing-advance"}) {
			bool armed = true;
			int advances = 0;
			auto selected = flow;
			if (mode == "1d local flow solve") {
				selected.scheme = iga::OneDFlowScheme::SteadyPoiseuille;
				selected.model = iga::OneDFlowModel::Rigid;
			}
			iga::OneDFlowRuntime::ImplicitAdvance advance;
			if (mode != "missing-advance" || rank != target)
				advance = [&](const iga::OneDNetwork& n, const iga::OneDFlowSystemDefinition& f,
					iga::OneDFlowState& state, double q, double step) {
					++advances;
					iga::AdvanceImplicitOneD(n, f, state, q, step, group);
				};
			iga::OneDFlowRuntime runtime(configuration, selected, network, inlet, directory, advance,
				[&](const char* stage, std::exception_ptr error) {
					iga::CollectiveLocalStage(group, stage, [&] {
						if (error) std::rethrow_exception(error);
						if (armed && rank == target && mode == stage)
							throw std::runtime_error("injected local trial phase failure");
					});
				});
			runtime.InitializeOpenLoop(seed);
			const auto before = Packed(runtime);
			runtime.BeginStep(0.0, dt);
			if (mode != "missing-inlet" || rank != target)
				runtime.SetOpenLoopInlet(runtime.OpenLoopInlet(dt, seed*(1.0+0.1*repetition)));
			bool rejected = false;
			try { runtime.SolveTrial(); }
			catch (const std::runtime_error& error) {
				const auto stage = mode == "missing-inlet" ? "1d trial preparation" :
					mode == "missing-advance" ? "1d substep preparation" : mode;
				rejected = std::string(error.what()).find(stage+": rank "+std::to_string(target)+":") != std::string::npos;
			}
			Require(rejected, "trial failure was not coordinated on every group rank");
			const int expected_advances = mode == "1d substep transport" ? 1 : mode == "1d trial completion" ? 3 : 0;
			Require(advances == expected_advances, "trial entered an unexpected implicit collective after failure");
			bool commit_rejected = false;
			try { runtime.PrepareCommitStep(); }
			catch (const std::runtime_error&) { commit_rejected = true; }
			Require(commit_rejected, "failed trial could be committed");
			runtime.AbortStep();
			Require(runtime.CurrentPhase() == iga::OneDFlowRuntime::Phase::Ready
				&& runtime.FlowState().completed_step == 0 && runtime.FlowState().physical_time == 0.0
				&& Packed(runtime) == before, "abort did not restore flow, species, outlets and dynamic radius");
			if (mode == "missing-advance") continue;
			armed = false;
			runtime.BeginStep(0.0, dt);
			runtime.SetOpenLoopInlet(runtime.OpenLoopInlet(dt, seed*(1.0+0.1*repetition)));
			runtime.SolveTrial();
			runtime.CommitStep();
			iga::OneDFlowRuntime reference(configuration, selected, network, inlet, directory,
				[](const iga::OneDNetwork& n, const iga::OneDFlowSystemDefinition& f,
					iga::OneDFlowState& state, double q, double step) {
					iga::AdvanceImplicitOneD(n, f, state, q, step, PETSC_COMM_SELF);
				});
			reference.InitializeOpenLoop(seed);
			reference.BeginStep(0.0, dt);
			reference.SetOpenLoopInlet(reference.OpenLoopInlet(dt, seed*(1.0+0.1*repetition)));
			reference.SolveTrial();
			reference.CommitStep();
			const auto expected = Packed(reference), current = Packed(runtime);
			Require(expected.size() == current.size(), "retry field size differs");
			double norm = 0.0, error = 0.0;
			for (std::size_t i = 0; i < expected.size(); ++i) {
				Require(std::isfinite(expected[i]) && std::isfinite(current[i]), "nonfinite retry field");
				norm = std::hypot(norm, expected[i]);
				error = std::hypot(error, expected[i]-current[i]);
			}
			Require(norm == 0.0 ? error <= 1e-12 : error/norm <= 1e-6, "retry differs from serial reference");
			Require(runtime.FlowState().completed_step == 3, "retry completed wrong substep count");
		}
	}
}

void CompareStaged(const iga::OneDFlowRuntime& reference, const iga::OneDFlowRuntime& current)
{
	auto compare = [](const std::vector<double>& a, const std::vector<double>& b) {
		Require(a.size() == b.size(), "staged field size differs");
		double norm = 0.0, error = 0.0;
		for (std::size_t i = 0; i < a.size(); ++i) {
			Require(std::isfinite(a[i]) && std::isfinite(b[i]), "nonfinite staged field");
			norm = std::hypot(norm, a[i]);
			error = std::hypot(error, a[i]-b[i]);
		}
		Require(norm == 0.0 ? error <= 1e-12 : error/norm <= 1e-6, "staged retry differs from serial reference");
	};
	const auto& a = reference.FlowState();
	const auto& b = current.FlowState();
	compare(a.area, b.area); compare(a.flow, b.flow); compare(a.pressure, b.pressure);
	compare(a.node_pressure, b.node_pressure); compare(a.segment_flow, b.segment_flow);
	Require(a.outlets.size() == b.outlets.size(), "staged outlet count differs");
	for (std::size_t i = 0; i < a.outlets.size(); ++i) {
		compare({a.outlets[i].pressure}, {b.outlets[i].pressure});
		compare({a.outlets[i].capacitor_pressure}, {b.outlets[i].capacitor_pressure});
		compare({a.outlets[i].flow}, {b.outlets[i].flow});
	}
	Require(reference.Transports().size() == current.Transports().size(), "transport count differs");
	for (std::size_t i = 0; i < reference.Transports().size(); ++i) {
		const auto& first = reference.Transports()[i];
		const auto& second = current.Transports()[i];
		Require(first.species.size() == second.species.size(), "species count differs");
		for (std::size_t j = 0; j < first.species.size(); ++j) {
			Require(first.species[j].definition.field == second.species[j].definition.field, "species identity differs");
			compare(first.species[j].concentration, second.species[j].concentration);
		}
	}
}

void CheckStagedGroup(MPI_Comm group, const std::filesystem::path& directory, int repetitions)
{
	int rank = 0, ranks = 1;
	MPI_Comm_rank(group, &rank); MPI_Comm_size(group, &ranks);
	const int target = ranks > 1 ? 1 : 0;
	std::ifstream file(directory/"simulation_config.json");
	std::ostringstream text; text << file.rdbuf();
	auto configuration = iga::ParseOneDConfiguration(text.str());
	configuration.physiology.vasodilation = false; // Existing staged capability requirement.
	const auto inlet = iga::ResolveOneDInlet(configuration);
	const auto flow = configuration.flow_systems.front();
	const auto network = iga::ReadOneDNetwork(directory/configuration.geometry.file,
		configuration.geometry.length_scale_to_m, flow.discretization.cells_per_segment,
		flow.dynamic_viscosity, configuration.geometry.root_node_id);
	const double seed = iga::EvaluateOneDInlet(configuration, inlet, directory, 0.0, network.segments.front().area0);
	const double dt = 3*configuration.time.dt;
	for (int repeat = 0; repeat < repetitions; ++repeat) {
		for (const std::string mode : {"1d hydraulic preparation", "1d hydraulic substep preparation",
			"1d hydraulic frame capture", "1d hydraulic completion", "1d hydraulic local solve",
			"1d staged transport", "nonfinite-concentration"}) {
			bool armed = true;
			int advances = 0;
			auto selected = flow;
			if (mode == "1d hydraulic local solve") {
				selected.scheme = iga::OneDFlowScheme::SteadyPoiseuille;
				selected.model = iga::OneDFlowModel::Rigid;
			}
			iga::OneDFlowRuntime runtime(configuration, selected, network, inlet, directory,
				[&](const iga::OneDNetwork& n, const iga::OneDFlowSystemDefinition& f,
					iga::OneDFlowState& state, double q, double step) {
					++advances;
					iga::AdvanceImplicitOneD(n, f, state, q, step, group);
				}, [&](const char* stage, std::exception_ptr error) {
					iga::CollectiveLocalStage(group, stage, [&] {
						if (error) std::rethrow_exception(error);
						if (armed && rank == target && mode == stage)
							throw std::runtime_error("injected staged phase failure");
					});
				});
			runtime.InitializeOpenLoop(seed);
			const auto initial = Packed(runtime);
			const auto port = runtime.OpenLoopInlet(dt, seed*(1.0+0.1*repeat));
			runtime.BeginStep(0.0, dt); runtime.SetOpenLoopInlet(port);
			const bool scalar_fault = mode == "1d staged transport" || mode == "nonfinite-concentration";
			if (!scalar_fault) {
				bool rejected = false;
				try { runtime.SolveHydraulicTrial(); }
				catch (const std::runtime_error& error) {
					rejected = std::string(error.what()).find(mode+": rank "+std::to_string(target)+":") != std::string::npos;
				}
				Require(rejected, "hydraulic fault not coordinated");
				Require(advances == (mode == "1d hydraulic frame capture" ? 1 : mode == "1d hydraulic completion" ? 3 : 0),
					"unexpected solve after hydraulic fault");
				bool scalar_rejected = false;
				try { runtime.SolveStagedTransportTrial(port.species, {}); }
				catch (const std::runtime_error&) { scalar_rejected = true; }
				Require(scalar_rejected, "failed hydraulic trial allowed scalar replay");
				runtime.AbortStep();
				Require(Packed(runtime) == initial && runtime.FlowState().completed_step == 0,
					"hydraulic abort did not restore snapshot");
				armed = false;
				runtime.BeginStep(0.0, dt); runtime.SetOpenLoopInlet(port);
			}
			runtime.SolveHydraulicTrial();
			const auto hydraulics = Packed(runtime);
			const auto last_concentrations = runtime.GetPortState("root").concentration;
			const auto solve_count = advances;
			if (scalar_fault) {
				auto concentrations = port.species;
				if (mode == "nonfinite-concentration" && rank == target) concentrations.at("oxygen") = std::numeric_limits<double>::quiet_NaN();
				bool rejected = false;
				try { runtime.SolveStagedTransportTrial(concentrations, {}); }
				catch (const std::runtime_error& error) {
					rejected = std::string(error.what()).find("1d staged transport: rank "+std::to_string(target)+":") != std::string::npos;
				}
				Require(rejected && Packed(runtime) == hydraulics && advances == solve_count,
					"rejected scalar replay modified accepted hydraulic/scalar state");
				Require(runtime.CurrentPhase() == iga::OneDFlowRuntime::Phase::HydraulicSolved
					&& runtime.GetPortState("root").concentration == last_concentrations,
					"rejected scalar replay published inlet or phase");
				armed = false;
			}
			runtime.SolveStagedTransportTrial(port.species, {});
			const auto accepted = Packed(runtime);
			runtime.RollbackStagedTransportTrial();
			Require(Packed(runtime) == hydraulics, "scalar rollback changed hydraulic snapshot");
			runtime.SolveStagedTransportTrial(port.species, {});
			Require(Packed(runtime) == accepted && advances == solve_count, "scalar retry recomputed flow or changed values");
			runtime.CommitStep();
			iga::OneDFlowRuntime reference(configuration, selected, network, inlet, directory,
				[](const iga::OneDNetwork& n, const iga::OneDFlowSystemDefinition& f,
					iga::OneDFlowState& state, double q, double step) {
					iga::AdvanceImplicitOneD(n, f, state, q, step, PETSC_COMM_SELF);
				});
			reference.InitializeOpenLoop(seed);
			reference.BeginStep(0.0, dt); reference.SetOpenLoopInlet(port);
			reference.SolveHydraulicTrial(); reference.SolveStagedTransportTrial(port.species, {});
			reference.CommitStep();
			CompareStaged(reference, runtime);
			Require(runtime.FlowState().completed_step == 3, "staged retry has wrong step count");
		}
	}
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int rank = 0, ranks = 1;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
	try {
		Require(argc == 2 && ranks == 3, "require three ranks and a native multispecies case");
		int threads = 0;
		#ifdef _OPENMP
		#pragma omp parallel reduction(+:threads)
		#endif
		{ threads += 1; }
		CheckGroup(PETSC_COMM_WORLD, argv[1], 1);
		CheckStagedGroup(PETSC_COMM_WORLD, argv[1], 1);
		MPI_Comm group = MPI_COMM_NULL;
		MPI_Comm_split(PETSC_COMM_WORLD, rank == 0 ? 0 : 1, rank, &group);
		CheckGroup(group, argv[1], rank == 0 ? 1 : 2);
		CheckStagedGroup(group, argv[1], rank == 0 ? 1 : 2);
		MPI_Comm_free(&group);
		std::cout << "trial failure rank=" << rank << " passed threads=" << threads << '\n';
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n';
		MPI_Abort(PETSC_COMM_WORLD, 2);
	}
	PetscFinalize();
	return 0;
}
