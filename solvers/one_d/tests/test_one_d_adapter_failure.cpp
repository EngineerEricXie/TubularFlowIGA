#include "OneDFlowDomainAdapter.hpp"
#include "OneDCheckpoint.hpp"

#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>

namespace {

void Require(bool condition, const std::string& message)
{
	if (!condition) throw std::runtime_error(message);
}

std::vector<double> Packed(const iga::OneDFlowRuntime& runtime)
{
	return iga::PackOneDCheckpointState(runtime.FlowState(), runtime.Transports(), runtime.Network());
}

void Compare(const iga::OneDFlowRuntime& a, const iga::OneDFlowRuntime& b)
{
	const auto compare = [](const std::vector<double>& x, const std::vector<double>& y) {
		Require(x.size() == y.size(), "adapter field size differs");
		double norm = 0.0, error = 0.0;
		for (std::size_t i = 0; i < x.size(); ++i) {
			Require(std::isfinite(x[i]) && std::isfinite(y[i]), "nonfinite adapter field");
			norm = std::hypot(norm, x[i]); error = std::hypot(error, x[i]-y[i]);
		}
		Require(norm == 0.0 ? error <= 1e-12 : error/norm <= 1e-6, "adapter retry field differs");
	};
	compare(a.FlowState().area, b.FlowState().area); compare(a.FlowState().flow, b.FlowState().flow);
	compare(a.FlowState().pressure, b.FlowState().pressure);
	compare(a.FlowState().node_pressure, b.FlowState().node_pressure);
	compare(a.FlowState().segment_flow, b.FlowState().segment_flow);
	Require(a.FlowState().outlets.size() == b.FlowState().outlets.size(), "outlet count differs");
	for (std::size_t i = 0; i < a.FlowState().outlets.size(); ++i) {
		const auto& x = a.FlowState().outlets[i]; const auto& y = b.FlowState().outlets[i];
		compare({x.pressure}, {y.pressure}); compare({x.capacitor_pressure}, {y.capacitor_pressure});
		compare({x.flow}, {y.flow});
	}
	Require(a.Transports().size() == b.Transports().size(), "transport count differs");
	for (std::size_t i = 0; i < a.Transports().size(); ++i) {
		Require(a.Transports()[i].species.size() == b.Transports()[i].species.size(), "species count differs");
		for (std::size_t j = 0; j < a.Transports()[i].species.size(); ++j) {
			const auto& x = a.Transports()[i].species[j]; const auto& y = b.Transports()[i].species[j];
			Require(x.definition.field == y.definition.field, "species identity differs");
			compare(x.concentration, y.concentration);
		}
	}
}

void CheckGroup(MPI_Comm group, const std::filesystem::path& directory, int repetitions)
{
	int rank = 0, ranks = 1;
	MPI_Comm_rank(group, &rank); MPI_Comm_size(group, &ranks);
	const int target = ranks > 1 ? 1 : 0;
	std::ifstream input(directory/"simulation_config.json");
	std::ostringstream text; text << input.rdbuf();
	auto configuration = iga::ParseOneDConfiguration(text.str());
	configuration.physiology.vasodilation = false;
	const auto flow = configuration.flow_systems.front();
	const auto inlet = iga::ResolveOneDInlet(configuration);
	const auto network = iga::ReadOneDNetwork(directory/configuration.geometry.file,
		configuration.geometry.length_scale_to_m, flow.discretization.cells_per_segment,
		flow.dynamic_viscosity, configuration.geometry.root_node_id);
	const auto seed = iga::EvaluateOneDInlet(configuration, inlet, directory, 0.0, network.segments.front().area0);
	const double dt = 3*configuration.time.dt;
	int checked = 0;
	for (int repeat = 0; repeat < repetitions; ++repeat) for (bool staged : {false, true}) {
		const std::string prefix = staged ? "1d staged adapter " : "1d flow adapter ";
		std::vector<std::string> modes = {"begin step", "input", "prepare commit", "abort"};
		if (staged) modes.insert(modes.end(), {"hydraulic preparation", "concentration", "transport state", "transport preparation",
			"rollback preparation", "hydraulic rollback preparation", "transport rollback preparation",
			"bad-concentration", "missing-concentration", "missing-hydraulics"});
		else modes.insert(modes.end(), {"rollback", "bad-input-time"});
		for (const auto& mode : modes) {
			bool armed = true;
			int advances = 0;
			iga::OneDFlowRuntime runtime(configuration, flow, network, inlet, directory,
				[&](const iga::OneDNetwork& n, const iga::OneDFlowSystemDefinition& f,
					iga::OneDFlowState& state, double q, double step) {
					++advances; iga::AdvanceImplicitOneD(n, f, state, q, step, group);
				}, [&](const char* stage, std::exception_ptr error) {
					iga::CollectiveLocalStage(group, stage, [&] {
						if (error) std::rethrow_exception(error);
						if (armed && rank == target && stage == prefix+mode)
							throw std::runtime_error("injected adapter local-stage failure");
					});
				});
			runtime.InitializeOpenLoop(seed);
			const auto initial = Packed(runtime);
			const auto native_inlet = runtime.OpenLoopInlet(dt, seed*(1.0+0.05*repeat));
			iga::CouplingPort root;
			root.id = "logical-root"; root.subsystem_id = "domain"; root.locator_kind = "runtime_port"; root.locator = "root";
			root.provides = {iga::PortQuantity::FlowRate, iga::PortQuantity::MeanPressure, iga::PortQuantity::Area};
			root.requires = {iga::PortQuantity::FlowRate};
			std::map<std::string, std::string> bindings;
			std::map<std::string, double> concentrations;
			if (staged) {
				root.provides.insert(iga::PortQuantity::SpeciesConcentration);
				root.provides.insert(iga::PortQuantity::SpeciesFlux);
				root.requires.insert(iga::PortQuantity::SpeciesConcentration);
				for (const auto& species : native_inlet.species) {
					const auto logical = "logical_"+species.first;
					root.species.insert(logical); bindings.emplace(logical, species.first);
					concentrations.emplace(logical, species.second);
				}
			}
			std::unique_ptr<iga::CoupledDomainRuntime> adapter;
			if (staged) adapter = std::make_unique<iga::OneDFlowTransportDomainAdapter>(
				"domain", runtime, std::vector<iga::CouplingPort>{root}, iga::OneDInletPolicy::CoupledRoot, bindings);
			else adapter = std::make_unique<iga::OneDFlowDomainAdapter>(
				"domain", runtime, std::vector<iga::CouplingPort>{root}, iga::OneDInletPolicy::CoupledRoot);
			auto* scalar = dynamic_cast<iga::OneDFlowTransportDomainAdapter*>(adapter.get());
			iga::PortBoundaryData boundary;
			boundary.time_s = dt; boundary.outward_flow_m3_s = -native_inlet.flow_m3_s;
			if (!staged) boundary.concentration = native_inlet.species;
			if (staged && mode == "input") boundary.concentration = concentrations;
			std::vector<double> hydraulic_image;
			const auto solve = [&](bool faulty) {
				adapter->BeginStep({0, 0.0, dt});
				auto supplied = boundary;
				if (faulty && mode == "bad-input-time" && rank == target) supplied.time_s = dt*2;
				if (faulty && (mode == "missing-concentration" || mode == "missing-hydraulics")) {
					supplied.concentration = concentrations;
					if (rank == target) {
						if (mode == "missing-hydraulics") supplied.outward_flow_m3_s.reset();
						else supplied.concentration.clear();
					}
				}
				adapter->SetPortInput(root.id, supplied);
				if (scalar) {
					scalar->SolveHydraulicTrial();
					hydraulic_image = Packed(runtime);
					if (rank == 0) (void)scalar->GetHydraulicPortState(root.id);
					MPI_Barrier(group); // Port queries must remain local.
					if (!(faulty && (mode == "missing-concentration" || mode == "missing-hydraulics"))
						&& supplied.concentration.empty()) {
						auto values = concentrations;
						if (faulty && mode == "bad-concentration" && rank == target) values.erase(values.begin());
						scalar->SetTransportConcentration(root.id, dt, values);
					}
					scalar->SolveTransportTrial();
					if (!faulty) {
						const auto accepted = Packed(runtime);
						bool hydraulic_rejected = false, transport_rejected = false;
						try { scalar->SolveHydraulicTrial(); }
						catch (const std::runtime_error&) { hydraulic_rejected = true; }
						try { scalar->SolveTransportTrial(); }
						catch (const std::runtime_error&) { transport_rejected = true; }
						Require(hydraulic_rejected && transport_rejected && Packed(runtime) == accepted,
							"illegal repeated solve changed accepted state");
						// Port queries and PrepareCommit below must still accept the trial.
					}
				} else adapter->SolveTrial();
				if (rank == 0) {
					(void)adapter->GetPortState(root.id);
					if (scalar) (void)scalar->GetSpeciesStepAccounting();
				}
				MPI_Barrier(group);
				if (faulty && mode.find("rollback") != std::string::npos) adapter->RollbackTrial();
				else if (faulty && mode == "abort") adapter->AbortStep();
				else { adapter->PrepareCommitStep(); adapter->FinalizeCommitStep(); }
			};
			std::string stage = prefix+mode;
			if (mode == "bad-input-time") stage = prefix+"input";
			if (mode == "bad-concentration") stage = prefix+"concentration";
			if (mode == "missing-concentration") stage = prefix+"transport preparation";
			if (mode == "missing-hydraulics") stage = prefix+"hydraulic preparation";
			bool rejected = false;
			try { solve(true); }
			catch (const std::runtime_error& error) {
				rejected = std::string(error.what()).find(stage+": rank "+std::to_string(target)+":") != std::string::npos;
			}
			Require(rejected, "adapter error not coordinated: "+mode);
			const bool before_solve = mode == "begin step" || mode == "input"
				|| mode == "hydraulic preparation" || mode == "missing-hydraulics" || mode == "bad-input-time";
			Require(advances == (before_solve ? 0 : 3), "adapter failure entered an unexpected implicit solve");
			if (staged && (mode == "concentration" || mode == "bad-concentration"
				|| mode == "missing-concentration" || mode == "transport state" || mode == "transport preparation"))
				Require(Packed(runtime) == hydraulic_image, "adapter preparation changed accepted hydraulic/scalar image");
			armed = false;
			if (staged && mode == "transport state") {
				const auto count = advances;
				scalar->SolveTransportTrial(); // State rejection must preserve the supplied map.
				Require(advances == count, "transport state retry recomputed hydraulics");
			}
			if (staged && mode == "input") adapter->SetPortInput(root.id, boundary);
			if (staged && (mode == "concentration" || mode == "bad-concentration" || mode == "missing-concentration")) {
				const auto hydraulic = Packed(runtime); const auto count = advances;
				scalar->SetTransportConcentration(root.id, dt, concentrations);
				Require(Packed(runtime) == hydraulic, "concentration retry changed fields before replay");
				scalar->SolveTransportTrial();
				Require(advances == count, "concentration retry recomputed hydraulics");
			}
			adapter->AbortStep();
			Require(Packed(runtime) == initial && runtime.CurrentPhase() == iga::OneDFlowRuntime::Phase::Ready
				&& runtime.FlowState().completed_step == 0 && runtime.FlowState().physical_time == 0.0,
				"adapter abort failed to restore committed state");
			solve(false);
			iga::OneDFlowRuntime reference(configuration, flow, network, inlet, directory,
				[](const iga::OneDNetwork& n, const iga::OneDFlowSystemDefinition& f,
					iga::OneDFlowState& state, double q, double step) { iga::AdvanceImplicitOneD(n, f, state, q, step, PETSC_COMM_SELF); });
			reference.InitializeOpenLoop(seed); reference.BeginStep(0.0, dt); reference.SetOpenLoopInlet(native_inlet);
			if (staged) { reference.SolveHydraulicTrial(); reference.SolveStagedTransportTrial(native_inlet.species, {},
				iga::OneDStagedRootTransportOwnership::BoundaryConcentration); }
			else reference.SolveTrial();
			reference.CommitStep(); Compare(reference, runtime);
			Require(runtime.FlowState().completed_step == 3, "adapter retry has wrong step count");
			++checked;
		}
	}
	if (rank == 0) std::cout << "adapter failure group ranks=" << ranks << " cases=" << checked << " passed\n";
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int rank = 0, ranks = 1;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank); MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
	try {
		Require(argc == 2 && ranks == 3, "require three ranks and a native multispecies case");
		int threads = 0;
		#ifdef _OPENMP
		#pragma omp parallel reduction(+:threads)
		#endif
		{ threads += 1; }
		CheckGroup(PETSC_COMM_WORLD, argv[1], 1);
		MPI_Comm group = MPI_COMM_NULL;
		MPI_Comm_split(PETSC_COMM_WORLD, rank == 0 ? 0 : 1, rank, &group);
		CheckGroup(group, argv[1], rank == 0 ? 1 : 2);
		MPI_Comm_free(&group);
		std::cout << "adapter failure rank=" << rank << " passed threads=" << threads << '\n';
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n'; MPI_Abort(PETSC_COMM_WORLD, 2);
	}
	PetscFinalize();
	return 0;
}
