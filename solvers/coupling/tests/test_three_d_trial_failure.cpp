#include "BoundarySupport.hpp"
#include "ThreeDBodyFittedFlowDomainAdapter.hpp"
#include "TransientTransportRuntime.hpp"
#include <utility>
#include <streambuf>
#include "../../cpu/tests/StringStreamFailure.hpp"

#define main BifurcationFixtureMain
#include "test_bifurcation_coupling_smoke.cpp"
#undef main

namespace {

template <class Work>
void Failure(MPI_Comm comm, const std::string& stage, Work&& work)
{
	std::string diagnostic;
	try { work(); } catch (const std::exception& error) { diagnostic = error.what(); }
	iga::CollectiveLocalStage(comm, "test failure check", [&] {
		if (diagnostic.find(stage) == std::string::npos)
			throw std::runtime_error("missing/wrong failure: "+diagnostic);
	});
	iga::RequireCollectiveSameText(comm, "test failure diagnostic", diagnostic);
}

std::vector<double> OwnedState(Vec state)
{
	PetscInt count = 0;
	VecGetLocalSize(state, &count);
	iga::PetscReadArray view;
	view.Acquire(state);
	std::vector<double> result(static_cast<std::size_t>(count));
	for (PetscInt i = 0; i < count; ++i) result[i] = PetscRealPart(view.Data()[i]);
	return result;
}

void Require(bool condition, const char* message)
{
	if (!condition) throw std::runtime_error(message);
}

void Near(double expected, double actual)
{
	if (!std::isfinite(actual) || !std::isfinite(expected)
		|| (expected == 0.0 ? std::abs(actual) > 1e-12
			: std::abs((actual-expected)/expected) > 1e-6))
		throw std::runtime_error("trial quantity differs from reference");
}

class FailingOutput : private std::streambuf {
public:
	explicit FailingOutput(bool enabled) : enabled_(enabled)
	{
		if (!enabled_) return;
		mask_ = std::cout.exceptions(); state_ = std::cout.rdstate();
		previous_ = std::cout.rdbuf(this);
		std::cout.exceptions(std::ios::badbit | std::ios::failbit);
	}
	~FailingOutput()
	{
		if (!enabled_) return;
		std::cout.exceptions(std::ios::goodbit);
		std::cout.rdbuf(previous_); std::cout.clear(state_); std::cout.exceptions(mask_);
	}
private:
	int_type overflow(int_type) override { throw std::runtime_error("injected local output failure"); }
	std::streamsize xsputn(const char*, std::streamsize) override { throw std::runtime_error("injected local output failure"); }
	bool enabled_;
	std::streambuf* previous_ = nullptr;
	std::ios::iostate mask_ = std::ios::goodbit, state_ = std::ios::goodbit;
};

std::vector<double> Run(const fs::path& root, const fs::path& path, MPI_Comm comm, bool faults)
{
	int rank = 0, ranks = 1;
	MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &ranks);
	int cases = 0;
	iga::Database database(path.string());
	const auto mesh = iga::ReadLabeledHexMesh((root/"controlmesh.vtk").string(), 64, 1);
	const auto velocity = iga::ReadVelocity((root/"initial_velocityfield.txt").string(), 64);
	const auto configuration = iga::ReadSimulationConfiguration((root/"simulation_config.json").string());
	const auto& definition = iga::FirstNavierStokesSystem(configuration);
	const auto boundaries = iga::ResolveFlowBoundaries(configuration, definition, mesh.labels, velocity);
	iga::TransientFlowRuntime flow(database, comm, true, true,
		{definition.density, definition.viscosity, kDt}, boundaries, mesh.labels,
		velocity, iga::WallTraceBasis(database, mesh, 0), {});
	flow.InitializeState(configuration);
	const auto initial = OwnedState(flow.State());
	const auto begin = [&] { flow.BeginStep(0, kDt, 12, 1e-8, 1e-12, 1e-6); };
	if (faults) {
		iga_test::ExpectStringStreamFailure(comm, "flow initialization preparation", [&] { flow.InitializeState(configuration); });
		Require(OwnedState(flow.State()) == initial, "failed signature changed initial state");
		flow.InitializeState(configuration);
		iga_test::ExpectStringStreamFailure(comm, "flow begin step preparation", begin);
		Require(flow.Phase() == iga::FlowStepPhase::Committed && OwnedState(flow.State()) == initial,
			"failed signature opened or changed flow trial");
		cases += 6;
		for (int mode = 0; mode < 3; ++mode) {
			const bool configured = mode == 0, transient = mode == 2;
			iga::TransientFlowRuntime alternate(database,comm,configured,transient,
				{definition.density,definition.viscosity,transient ? kDt : 0.0},boundaries,mesh.labels,
				velocity,iga::WallTraceBasis(database,mesh,0),{});
			if (!transient) VecSet(alternate.State(),3.0);
			const auto before = OwnedState(alternate.State());
			alternate.InitializeState();
			Require(OwnedState(alternate.State()) == (transient ? initial : before),
				"unconfigured/steady initialization changed existing semantics");
			auto configured_initial = configuration;
			configured_initial.boundaries[2].conditions[0].value = {0.25};
			alternate.InitializeState(configured_initial);
			Require(OwnedState(alternate.State()) == (transient ? initial : before),
				"steady or ignored configuration changed field coefficients");
			if (configured) Near(0.25,alternate.PressureTractionValue(2).value());
		}
		const auto traction = flow.PressureTractionValue(2);
		for (int mode = 0; mode < 6; ++mode) {
			if (mode >= 3 && ranks == 1) continue;
			auto invalid = configuration;
			if (mode == 0) {
				invalid.boundaries[1].conditions[0].scale = 2.0;
				invalid.boundaries[2].conditions[0].value = {0.25};
			}
			if (rank == ranks-1) {
				if (mode == 0) invalid.boundaries[1].conditions[0].profile = "unsupported";
				if (mode == 1) invalid.boundaries[1].conditions[0].scale = std::numeric_limits<double>::quiet_NaN();
				if (mode == 2) invalid.boundaries[2].conditions[0].kind = iga::FieldBoundaryKind::Dirichlet;
				if (mode == 3) invalid.boundaries[1].conditions[0].scale = 2.0;
				if (mode == 4) invalid.boundaries[2].conditions[0].value = {0.25};
			}
			Failure(comm, mode < 3 ? "flow initialization preparation" : "flow initialization agreement", [&] {
				if (mode == 5 && rank == ranks-1) flow.InitializeState();
				else flow.InitializeState(invalid);
			});
			Require(OwnedState(flow.State()) == initial && flow.PressureTractionValue(2) == traction,
				"rejected initialization published state or traction");
			flow.InitializeState();
			Require(OwnedState(flow.State()) == initial && flow.PressureTractionValue(2) == traction,
				"rejected initialization leaked boundary configuration");
			++cases;
		}
		begin();
		Failure(comm,"flow initialization preparation",[&] { flow.InitializeState(configuration); });
		Require(flow.Phase() == iga::FlowStepPhase::TrialReady && OwnedState(flow.State()) == initial,
			"initialization during trial changed state or phase");
		flow.AbortStep(); ++cases;
		const auto handle = flow.State();
		auto scaled = configuration;
		scaled.boundaries[1].conditions[0].scale = 2.0;
		flow.InitializeState(scaled);
		flow.InitializeState(configuration);
		Require(flow.State() == handle && OwnedState(flow.State()) == initial,
			"initialization retry did not preserve handle or original field");
		Failure(comm, "flow begin step preparation", [&] {
			flow.BeginStep(0, kDt, rank == ranks-1 ? 0 : 12, 1e-8, 1e-12, 1e-6);
		});
		++cases;
		if (ranks > 1) {
			Failure(comm, "flow step controls agreement", [&] {
				flow.BeginStep(0, rank == ranks-1 ? 2*kDt : kDt, 12, 1e-8, 1e-12, 1e-6);
			});
			++cases;
		}
		Require(flow.Phase() == iga::FlowStepPhase::Committed && OwnedState(flow.State()) == initial,
			"failed flow begin changed committed state");
		auto invalid = configuration;
		if (rank == ranks-1) invalid.boundaries[1].conditions[0].profile = "unsupported";
		Failure(comm, "flow advance boundaries", [&] { flow.Advance(invalid, 0, kDt, 12, 1e-8, 1e-12, 1e-6); });
		Require(flow.Phase() == iga::FlowStepPhase::Committed && OwnedState(flow.State()) == initial,
			"failed Advance did not restore committed state");
		++cases;
		KSP local_solver = nullptr;
		KSPCreate(PETSC_COMM_SELF, &local_solver);
		PetscPushErrorHandler(PetscReturnErrorHandler, nullptr);
		const auto status = KSPSetTolerances(local_solver, rank == ranks-1 ? -1.0 : 1e-8,
			PETSC_DEFAULT, PETSC_DEFAULT, 10);
		PetscPopErrorHandler();
		Failure(comm, "returned status", [&] {
			iga::RequireCollectivePetscSuccess(comm, "returned status", status);
		});
		KSPDestroy(&local_solver);
		++cases;
		// Corrupt only a pressure coefficient; velocity history stays valid so
		// this exercises the post-assembly conservation diagnostics.
		PetscScalar* state_values = nullptr;
		PetscScalar pressure = 0.0;
		VecGetArray(flow.State(), &state_values);
		if (rank == ranks-1) {
			pressure = state_values[3]; state_values[3] = std::numeric_limits<double>::quiet_NaN();
		}
		VecRestoreArray(flow.State(), &state_values);
		Failure(comm,"flow conservation result",[&] {
			flow.Advance(configuration,0,kDt,12,1e-8,1e-12,1e-6);
		});
		Require(flow.Phase() == iga::FlowStepPhase::Committed,"invalid pressure left trial open");
		VecGetArray(flow.State(), &state_values);
		if (rank == ranks-1) state_values[3] = pressure;
		VecRestoreArray(flow.State(), &state_values);
		Require(OwnedState(flow.State()) == initial,"invalid pressure recovery changed other coefficients");
		++cases;
		// Exhaust the nonlinear budget after a genuine update, then verify the
		// convenience Advance path restores its committed snapshot.
		Failure(comm, "Navier-Stokes nonlinear solve reached MAX_NEWTON", [&] {
			flow.Advance(configuration, 0, kDt, 1, 1e-8, 1e-12, 1e-6);
		});
		Require(flow.Phase() == iga::FlowStepPhase::Committed && OwnedState(flow.State()) == initial,
			"nonlinear exhaustion failed to restore committed state"); ++cases;
		{
			FailingOutput output(rank == 0);
			Failure(comm, "flow iteration logging", [&] {
				flow.Advance(configuration, 0, kDt, 12, 1e-8, 1e-12, 1e-6);
			});
		}
		Require(flow.Phase() == iga::FlowStepPhase::Committed && OwnedState(flow.State()) == initial,
			"logging failure after update failed to restore state"); ++cases;
		for (int mode = 0; mode < 2; ++mode) {
			const std::vector<std::pair<std::string, std::string>> options{
				{"-ksp_type","richardson"}, {"-pc_type","none"}, {"-ksp_max_it",mode == 0 ? "1" : "2"},
				{"-ksp_rtol","1e-14"}, {"-ksp_atol","1e-50"},
				{"-ksp_divtol",mode == 0 ? "1e4" : "1.01"},
				{"-ksp_richardson_scale",mode == 0 ? "1" : "1e6"},
				{"-ksp_error_if_not_converged",mode == 0 ? "false" : "true"}};
			std::vector<std::pair<bool,std::string>> saved;
			for (const auto& option : options) {
				char value[1024]{}; PetscBool present = PETSC_FALSE;
				Require(PetscOptionsGetString(nullptr,nullptr,option.first.c_str(),value,sizeof(value),&present) == 0,
					"cannot read test solver option");
				saved.emplace_back(present,value);
				Require(PetscOptionsSetValue(nullptr,option.first.c_str(),option.second.c_str()) == 0,
					"cannot set test solver option");
			}
			iga::TransientFlowRuntime limited(database, comm, true, true,
				{definition.density, definition.viscosity, kDt}, boundaries, mesh.labels,
				velocity, iga::WallTraceBasis(database, mesh, 0), {});
			for (std::size_t i = 0; i < options.size(); ++i)
				Require((saved[i].first ? PetscOptionsSetValue(nullptr,options[i].first.c_str(),saved[i].second.c_str())
					: PetscOptionsClearValue(nullptr,options[i].first.c_str())) == 0,"cannot restore test solver option");
			limited.InitializeState(configuration);
			const auto snapshot = OwnedState(limited.State());
			limited.BeginStep(0,kDt,12,1e-8,1e-12,1e-6);
			limited.SetTrialBoundaryConfiguration(configuration);
			const auto stage = mode == 0 ? "flow linear convergence" : "flow linear solve";
			PetscPushErrorHandler(PetscReturnErrorHandler,nullptr);
			Failure(comm,stage,[&] { limited.SolveTrial(); });
			PetscPopErrorHandler();
			Failure(comm,"flow prepare commit",[&] { limited.PrepareCommitStep(); });
			Require(OwnedState(limited.State()) == snapshot && limited.TrialLinearIterations() == 0,
				"failed KSP published update or iterations");
			limited.AbortStep();
			Require(limited.Phase() == iga::FlowStepPhase::Committed && OwnedState(limited.State()) == snapshot,
				"failed KSP Abort did not restore state");
			PetscPushErrorHandler(PetscReturnErrorHandler,nullptr);
			Failure(comm,stage,[&] { limited.Advance(configuration,0,kDt,12,1e-8,1e-12,1e-6); });
			PetscPopErrorHandler();
			Require(limited.Phase() == iga::FlowStepPhase::Committed && OwnedState(limited.State()) == snapshot,
				"failed KSP Advance did not restore state");
			cases += 3;
		}
	}
	begin();
	if (faults) {
		if (rank != ranks-1) flow.SetTrialBoundaryConfiguration(configuration);
		Failure(comm, "flow solve preparation", [&] { flow.SolveTrial(); });
		Require(flow.Phase() == iga::FlowStepPhase::TrialReady && OwnedState(flow.State()) == initial,
			"flow solve preflight changed trial state");
		Failure(comm, "flow prepare commit", [&] { flow.PrepareCommitStep(); });
		Failure(comm, "flow rollback preparation", [&] { flow.RollbackTrial(); });
		cases += 3;
	}
	flow.SetTrialBoundaryConfiguration(configuration);
	flow.SolveTrial();
	const auto solved = OwnedState(flow.State());
	iga::CouplingPort port;
	port.id = "inlet"; port.subsystem_id = "flow";
	port.locator_kind = "boundary_label"; port.locator = "1";
	port.provides = {iga::PortQuantity::Area, iga::PortQuantity::FlowRate, iga::PortQuantity::MeanPressure};
	const auto measured = flow.GetPortState(port);
	Require(std::abs(measured.outward_flow_m3_s.value()) > 1e-4, "flow reference must be nonzero");
	flow.RollbackTrial();
	Require(OwnedState(flow.State()) == initial, "flow rollback snapshot differs");
	flow.SetTrialBoundaryConfiguration(configuration); flow.SolveTrial();
	const auto retry = OwnedState(flow.State());
	for (std::size_t i = 0; i < retry.size(); ++i) Near(solved[i], retry[i]);
	flow.PrepareCommitStep();
	if (faults && ranks > 1) {
		// Finalize is intentionally noexcept and local. A divergent caller must
		// be rejected by Abort's phase agreement before any VecCopy or return.
		if (rank == ranks-1) flow.FinalizeCommitStep();
		Failure(comm, "flow abort phase agreement", [&] { flow.AbortStep(); });
		if (rank != ranks-1) flow.FinalizeCommitStep();
		++cases;
	} else flow.FinalizeCommitStep();
	flow.AbortStep();
	Require(OwnedState(flow.State()) == retry, "committed Abort modified flow state");
	flow.CopyStateToPrevious();
	const auto required_velocity = flow.GatherRequiredVelocity();

	iga::SimulationConfiguration scalar;
	scalar.fields = {{"tracer", iga::FieldKind::Scalar, 2.0}};
	scalar.time = {kDt, 1};
	iga::EquationSystemDefinition system;
	system.name = "transport"; system.kind = iga::EquationKind::LinearTransport;
	system.unknowns = {"tracer"};
	system.terms = {{iga::TermKind::TimeDerivative, "tracer", "tracer", 1.0, ""},
		{iga::TermKind::VolumeSource, "tracer", "tracer", 1.0, ""}};
	scalar.equation_systems.push_back(system);
	scalar.velocity_sources.push_back({"prescribed", "prescribed", "", "", "error"});
	for (int label = 0; label < 4; ++label) {
		iga::FieldBoundaryCondition condition;
		condition.field = "tracer"; condition.kind = iga::FieldBoundaryKind::NoFlux;
		scalar.boundaries.push_back({label, "boundary"+std::to_string(label), {condition}});
	}
	const auto compiled = iga::CompileLinearSystem(scalar, "transport");
	iga::TransientTransportRuntime transport(database, comm, scalar, compiled, mesh.labels);
	const auto scalar_initial = transport.GatherRequiredState();
	transport.BeginStep();
	if (faults) {
		Failure(comm, "transport begin step preparation", [&] { transport.BeginStep(); });
		Failure(comm, "transport prepare commit", [&] { transport.PrepareCommitStep(); });
		Failure(comm, "transport rollback preparation", [&] { transport.RollbackTrial(); });
		auto invalid_nodes = transport.RequiredNodes();
		if (rank == ranks-1) invalid_nodes.pop_back();
		Failure(comm, "transport trial input", [&] { transport.SolveTrial(scalar, invalid_nodes, required_velocity); });
		Require(transport.Phase() == iga::TransportStepPhase::TrialOpen
			&& transport.GatherRequiredState() == scalar_initial, "transport preflight modified snapshot");
		cases += 4;
	}
	transport.SolveTrial(scalar, flow.RequiredNodes(), required_velocity);
	const auto scalar_solved = transport.GatherRequiredState();
	for (double value : scalar_solved) Near(2.0+kDt, value);
	transport.RollbackTrial();
	Require(transport.GatherRequiredState() == scalar_initial && transport.Steps() == 0,
		"transport rollback did not restore state and step");
	transport.SolveTrial(scalar, flow.RequiredNodes(), required_velocity);
	const auto scalar_retry = transport.GatherRequiredState();
	for (std::size_t i = 0; i < scalar_retry.size(); ++i) Near(scalar_solved[i], scalar_retry[i]);
	transport.PrepareCommitStep();
	if (faults && ranks > 1) {
		if (rank == ranks-1) transport.FinalizeCommitStep();
		Failure(comm, "transport abort phase agreement", [&] { transport.AbortStep(); });
		if (rank != ranks-1) transport.FinalizeCommitStep();
		++cases;
	} else transport.FinalizeCommitStep();
	transport.AbortStep();
	Require(transport.GatherRequiredState() == scalar_retry && transport.Steps() == 1,
		"committed transport Abort modified state");
	if (faults) {
		transport.BeginStep();
		auto invalid = scalar;
		if (rank == ranks-1) invalid.boundaries[1].conditions[0].kind = iga::FieldBoundaryKind::Resistance;
		Failure(comm, "transport trial boundaries", [&] { transport.SolveTrial(invalid, flow.RequiredNodes(), required_velocity); });
		Failure(comm, "transport prepare commit", [&] { transport.PrepareCommitStep(); });
		transport.AbortStep();
		Require(transport.GatherRequiredState() == scalar_retry && transport.Steps() == 1,
			"failed scalar trial Abort modified committed state");
		cases += 2;
	}
	// Exercise an actual post-solve local outlet-model exception. This test
	// owns the runtime; only the selected rank's model kind is made invalid.
	auto outlet_configuration = configuration;
	auto& condition = outlet_configuration.boundaries[2].conditions[0];
	condition.kind = iga::FieldBoundaryKind::WindkesselRC;
	condition.resistance = 1e-3; condition.capacitance = 1.0;
	const auto models = iga::InitializeOutletModels(outlet_configuration,
		iga::FirstNavierStokesSystem(outlet_configuration));
	const auto outlet_boundaries = iga::ResolveFlowBoundaries(
		iga::MaterializeOutletPressures(outlet_configuration, models), definition, mesh.labels, velocity);
	iga::TransientFlowRuntime outlet(database, comm, true, true,
		{definition.density, definition.viscosity, kDt}, outlet_boundaries, mesh.labels,
		velocity, iga::WallTraceBasis(database, mesh, 0), models);
	outlet.InitializeState(outlet_configuration);
	const auto outlet_initial = OwnedState(outlet.State());
	if (faults) {
		if (rank == ranks-1)
			outlet.OutletModels().at(0).kind = iga::FieldBoundaryKind::Dirichlet;
		outlet.BeginStep(0, kDt, 12, 1e-8, 1e-12, 1e-6);
		outlet.SetTrialBoundaryConfiguration(outlet_configuration);
		Failure(comm, "flow outlet evaluation", [&] { outlet.SolveTrial(); });
		const auto& rejected = std::as_const(outlet).OutletModels().at(0);
		Require(rejected.pressure == models.at(0).pressure
			&& rejected.flow == models.at(0).flow
			&& rejected.capacitor_pressure == models.at(0).capacitor_pressure,
			"rejected outlet evaluation published a partial candidate");
		Failure(comm, "flow prepare commit", [&] { outlet.PrepareCommitStep(); });
		outlet.AbortStep();
		Require(OwnedState(outlet.State()) == outlet_initial, "outlet failure did not restore flow snapshot");
		if (rank == ranks-1)
			outlet.OutletModels().at(0).kind = iga::FieldBoundaryKind::WindkesselRC;
		cases += 2;
	}
	outlet.Advance(outlet_configuration, 0, kDt, 12, 1e-8, 1e-12, 1e-6);
	if (rank == 0 && faults) std::cout << "three_d_trial_failure ranks=" << ranks << " cases=" << cases << " passed\n";
	return {measured.area_m2.value(), measured.outward_flow_m3_s.value(), measured.mean_pressure_pa.value(),
		transport.TotalMass().at("tracer"), outlet.OutletModels().at(0).flow,
		outlet.OutletModels().at(0).pressure, outlet.OutletModels().at(0).capacitor_pressure};
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int rank = 0, ranks = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank); MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
	MPI_Comm comm = MPI_COMM_NULL;
	try {
		Require(argc == 2 && ranks == 3, "usage: mpiexec -np 3 three_d_trial_failure_test FRESH_OUTPUT");
		for (int round = 0; round < 2; ++round) {
			const int color = round == 0 || rank == 0 ? 0 : 1;
			MPI_Comm_split(PETSC_COMM_WORLD, color, rank, &comm);
			int local_rank = 0, local_size = 1;
			MPI_Comm_rank(comm, &local_rank); MPI_Comm_size(comm, &local_size);
			const auto root = fs::path(argv[1])/(std::to_string(round)+"-"+std::to_string(color));
			iga::CollectiveLocalStage(comm, "fixture", [&] {
				if (local_rank == 0) {
					Require(fs::create_directories(root), "output already exists");
					WriteThreeDCase(root); WriteDatabase(root/"group.ntiga", local_size); WriteDatabase(root/"serial.ntiga", 1);
				}
			});
			std::vector<double> reference;
			iga::CollectiveLocalStage(comm, "serial fixture", [&] {
				if (local_rank == 0) reference = Run(root, root/"serial.ntiga", PETSC_COMM_SELF, false);
			});
			const auto actual = Run(root, root/"group.ntiga", comm, true);
			iga::CollectiveLocalStage(comm, "reference comparison", [&] {
				if (local_rank == 0) {
					Require(reference.size() == actual.size(), "reference quantity count");
					for (std::size_t i = 0; i < actual.size(); ++i) Near(reference[i], actual[i]);
				}
			});
			MPI_Comm_free(&comm);
		}
		if (rank == 0) std::cout << "three_d_trial_failure all comparisons passed\n";
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		MPI_Abort(PETSC_COMM_WORLD, 1);
		return 1;
	}
	PetscFinalize();
	return 0;
}
