#include "BoundarySupport.hpp"
#include "TransientFlowRuntime.hpp"
#include "TransientTransportRuntime.hpp"
#include <memory>

#define main BifurcationFixtureMain
#include "../../coupling/tests/test_bifurcation_coupling_smoke.cpp"
#undef main

namespace {
void Require(bool condition, const char* reason)
{
	if (!condition) throw std::runtime_error(reason);
}
template<class Work> void Reject(MPI_Comm comm, const char* stage, Work&& work)
{
	std::string message;
	try { work(); } catch (const std::exception& error) { message = error.what(); }
	iga::CollectiveLocalStage(comm, "checkpoint test expected rejection", [&] {
		if (message.find(stage) == std::string::npos) throw std::runtime_error("expected "+std::string(stage)+", got "+message);
	});
	iga::RequireCollectiveSameText(comm, "checkpoint test diagnostic agreement", message);
}
iga::SimulationConfiguration ScalarConfiguration()
{
	iga::SimulationConfiguration scalar;
	scalar.fields = {{"tracer", iga::FieldKind::Scalar, 2.0}}; scalar.time = {kDt, 8};
	iga::EquationSystemDefinition system; system.name = "transport";
	system.kind = iga::EquationKind::LinearTransport; system.unknowns = {"tracer"};
	system.terms = {{iga::TermKind::TimeDerivative, "tracer", "tracer", 1.0, ""},
		{iga::TermKind::VolumeSource, "tracer", "tracer", 1.0, ""},
		{iga::TermKind::Diffusion, "tracer", "tracer", 0.01, ""}};
	scalar.equation_systems.push_back(system);
	scalar.velocity_sources.push_back({"prescribed", "prescribed", "", "", "error"});
	for (int label = 0; label < 4; ++label) {
		iga::FieldBoundaryCondition condition; condition.field = "tracer"; condition.kind = iga::FieldBoundaryKind::NoFlux;
		scalar.boundaries.push_back({label, "boundary"+std::to_string(label), {condition}});
	}
	return scalar;
}

void Compare(MPI_Comm comm, const iga::FlowAcceptedCheckpointState& first, const iga::FlowAcceptedCheckpointState& second)
{
	iga::CollectiveLocalStage(comm, "checkpoint flow exact comparison", [&] {
		Require(first.configuration_identity_sha256 == second.configuration_identity_sha256
			&& first.accepted_steps == second.accepted_steps && first.accepted_time_s == second.accepted_time_s
			&& first.macro_dt_s == second.macro_dt_s && first.total_linear_iterations == second.total_linear_iterations && first.field.values == second.field.values
			&& first.field.row_begin == second.field.row_begin && first.field.row_end == second.field.row_end
			&& first.boundaries.velocity == second.boundaries.velocity && first.boundaries.pressure == second.boundaries.pressure
			&& first.pressure_tractions == second.pressure_tractions && first.outlets.size() == second.outlets.size(), "accepted flow states differ");
		for (std::size_t i = 0; i < first.outlets.size(); ++i)
			Require(first.outlets[i].pressure == second.outlets[i].pressure && first.outlets[i].flow == second.outlets[i].flow
				&& first.outlets[i].capacitor_pressure == second.outlets[i].capacitor_pressure, "accepted outlet states differ");
	});
}
void Compare(MPI_Comm comm, const iga::TransportAcceptedCheckpointState& first, const iga::TransportAcceptedCheckpointState& second)
{
	iga::CollectiveLocalStage(comm, "checkpoint transport exact comparison", [&] {
		Require(first.configuration_identity_sha256 == second.configuration_identity_sha256 && first.accepted_steps == second.accepted_steps
			&& first.field.values == second.field.values && first.field.row_begin == second.field.row_begin
			&& first.field.row_end == second.field.row_end, "accepted transport states differ");
	});
}

void Run(const fs::path& root, MPI_Comm comm, bool transient, bool outlet_model)
{
	int rank = 0, ranks = 1; MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &ranks);
	iga::Database database((root/"group.ntiga").string());
	const auto mesh = iga::ReadLabeledHexMesh((root/"controlmesh.vtk").string(), 64, 1);
	const auto velocity = iga::ReadVelocity((root/"initial_velocityfield.txt").string(), 64);
	auto configuration = iga::ReadSimulationConfiguration((root/"simulation_config.json").string());
	if (outlet_model) {
		auto& condition = configuration.boundaries[2].conditions[0]; condition.kind = iga::FieldBoundaryKind::WindkesselRC;
		condition.resistance = 1e-3; condition.capacitance = 1.0;
	}
	const auto& definition = iga::FirstNavierStokesSystem(configuration);
	const auto models = iga::InitializeOutletModels(configuration, definition);
	const auto boundaries = iga::ResolveFlowBoundaries(iga::MaterializeOutletPressures(configuration, models), definition, mesh.labels, velocity);
	const auto wall = iga::WallTraceBasis(database, mesh, 0);
	const auto scalar = ScalarConfiguration(); const auto system = iga::CompileLinearSystem(scalar, "transport");
	const std::string identity(64, static_cast<char>('a'+2*outlet_model+!transient+ranks-1));
	auto make_flow = [&] {
		auto runtime = std::make_unique<iga::TransientFlowRuntime>(database, comm, true, transient,
			iga::NavierStokesParameters{definition.density, definition.viscosity, transient ? kDt : 0.0}, boundaries, mesh.labels,
			velocity, wall, models, std::set<std::string>{}, identity, kDt);
		runtime->InitializeState(configuration); return runtime;
	};
	auto make_transport = [&] { return std::make_unique<iga::TransientTransportRuntime>(database, comm, scalar, system,
		mesh.labels, std::map<std::uint64_t, iga::VolumeQuadratureRule>{}, std::set<std::string>{}, identity); };
	auto original = make_flow(), restored = make_flow(); auto transport = make_transport(), restored_transport = make_transport();
	auto advance = [&](iga::TransientFlowRuntime& flow, iga::TransientTransportRuntime& species, int step) {
		auto changed = configuration; changed.boundaries[1].conditions[0].scale = 1.0+0.05*step;
		flow.Advance(changed, step, flow.AcceptedTime()+kDt, 12, 1e-8, 1e-12, 1e-6);
		const auto required = flow.GatherRequiredVelocity();
		species.Advance(scalar, flow.RequiredNodes(), required);
	};
	Reject(comm, "flow accepted checkpoint capture", [&] { original->CaptureCheckpointState(); });
	Reject(comm, "transport accepted checkpoint capture", [&] { transport->CaptureCheckpointState(); });
	for (int i = 0; i < 3; ++i) advance(*original, *transport, i);
	const auto flow_state = original->CaptureCheckpointState(); const auto transport_state = transport->CaptureCheckpointState();
	const auto initial = iga::CaptureOwnedCheckpointVector(restored->State()).values;
	const auto initial_transport = restored_transport->GatherRequiredState();
	for (int fault = 0; fault < 13; ++fault) {
		auto bad = flow_state;
		if (rank == ranks-1) {
			if (fault == 0) bad.configuration_identity_sha256[0] = identity[0] == 'f' ? 'e' : 'f';
			if (fault == 1) bad.accepted_steps = 0;
			if (fault == 2) bad.accepted_time_s += 1;
			if (fault == 3) ++bad.field.row_begin;
			if (fault == 4) bad.field.values.front() = NAN;
			if (fault == 5) bad.boundaries.velocity_constrained.front() ^= 1;
			if (fault == 6) bad.pressure_tractions[99] = 1;
			if (fault == 7) bad.total_linear_iterations = -1;
			if (fault == 8) bad.boundaries.pressure.front() = NAN;
			if (fault == 9) bad.field.values.pop_back();
			if (fault == 10) bad.boundaries.velocity_nodes += 1;
			if (fault == 11) bad.macro_dt_s *= 2;
			if (fault == 12) bad.macro_dt_s = NAN;
		}
		Reject(comm, "flow accepted checkpoint validation", [&] { restored->RestoreCheckpointState(bad); });
		iga::CollectiveLocalStage(comm, "checkpoint unchanged target", [&] {
			Require(restored->AcceptedSteps() == 0 && iga::CaptureOwnedCheckpointVector(restored->State()).values == initial, "bad restore changed target");
		});
	}
	for (int fault = 0; fault < 5; ++fault) {
		auto bad = transport_state;
		if (rank == ranks-1) {
			if (fault == 0) bad.configuration_identity_sha256[0] = identity[0] == 'f' ? 'e' : 'f';
			if (fault == 1) bad.accepted_steps = 0;
			if (fault == 2) ++bad.field.row_end;
			if (fault == 3) bad.field.values.front() = NAN;
			if (fault == 4) bad.field.values.pop_back();
		}
		Reject(comm, "transport accepted checkpoint validation", [&] { restored_transport->RestoreCheckpointState(bad); });
		const auto current = restored_transport->GatherRequiredState();
		iga::CollectiveLocalStage(comm, "checkpoint unchanged transport", [&] {
			Require(restored_transport->Steps() == 0 && current == initial_transport, "bad restore changed transport");
		});
	}
	if (ranks > 1) {
		auto bad = flow_state; if (rank == ranks-1) bad.boundaries.velocity.front()[0] += 1;
		Reject(comm, "flow checkpoint replicated state agreement", [&] { restored->RestoreCheckpointState(bad); });
		auto bad_transport = transport_state; if (rank == ranks-1) ++bad_transport.accepted_steps;
		Reject(comm, "transport checkpoint clock agreement", [&] { restored_transport->RestoreCheckpointState(bad_transport); });
	}
	const auto restore_start = std::chrono::steady_clock::now();
	restored->RestoreCheckpointState(flow_state); restored_transport->RestoreCheckpointState(transport_state);
	const double restore_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-restore_start).count();
	std::cout << "body_fitted_checkpoint_restore local_rank=" << rank << " ranks=" << ranks << " flow_rows=" << flow_state.field.values.size()
		<< " transport_rows=" << transport_state.field.values.size() << " seconds=" << restore_seconds << '\n';
	Compare(comm, original->CaptureCheckpointState(), restored->CaptureCheckpointState());
	Compare(comm, transport->CaptureCheckpointState(), restored_transport->CaptureCheckpointState());
	Reject(comm, "flow accepted checkpoint validation", [&] { restored->RestoreCheckpointState(flow_state); });
	Reject(comm, "transport accepted checkpoint validation", [&] { restored_transport->RestoreCheckpointState(transport_state); });
	original->BeginStep(3, original->AcceptedTime()+kDt, 12, 1e-8, 1e-12, 1e-6); transport->BeginStep();
	Reject(comm, "flow accepted checkpoint capture", [&] { original->CaptureCheckpointState(); });
	Reject(comm, "transport accepted checkpoint capture", [&] { transport->CaptureCheckpointState(); });
	original->SetTrialBoundaryConfiguration(configuration); original->SolveTrial();
	transport->SolveTrial(scalar, original->RequiredNodes(), original->GatherRequiredVelocity());
	original->PrepareCommitStep(); transport->PrepareCommitStep();
	Reject(comm, "flow accepted checkpoint capture", [&] { original->CaptureCheckpointState(); });
	original->AbortStep(); transport->AbortStep();
	Compare(comm, flow_state, original->CaptureCheckpointState()); Compare(comm, transport_state, transport->CaptureCheckpointState());
	for (int i = 3; i < 6; ++i) {
		advance(*original, *transport, i); advance(*restored, *restored_transport, i);
		Compare(comm, original->CaptureCheckpointState(), restored->CaptureCheckpointState());
		Compare(comm, transport->CaptureCheckpointState(), restored_transport->CaptureCheckpointState());
		const auto a = transport->TotalMass(), b = restored_transport->TotalMass();
		iga::CollectiveLocalStage(comm, "checkpoint mass comparison", [&] { Require(a == b, "restored mass differs"); });
	}
	auto exhausted = make_transport(); auto last = transport_state; last.accepted_steps = std::numeric_limits<int>::max(); exhausted->RestoreCheckpointState(last);
	Reject(comm, "transport begin step preparation", [&] { exhausted->BeginStep(); }); exhausted->Close();
	original->Close(); restored->Close(); transport->Close(); restored_transport->Close();
	Reject(comm, "flow accepted checkpoint capture", [&] { original->CaptureCheckpointState(); });
	Reject(comm, "transport accepted checkpoint capture", [&] { transport->CaptureCheckpointState(); });
	if (rank == 0) std::cout << "body_fitted_accepted_checkpoint ranks=" << ranks << " transient=" << transient
		<< " outlet=" << outlet_model << " restored_steps=3 continued_steps=3 exact_owned_fields=1 status=passed\n";
}
} // namespace

#ifdef IGA_BODY_FITTED_CHECKPOINT_FIXTURE_ONLY
int BodyFittedAcceptedCheckpointFixtureMain(int argc, char** argv)
#else
int main(int argc, char** argv)
#endif
{
	PetscInitialize(&argc, &argv, nullptr, nullptr); iga::CurrentPhaseProfile().EnableFromEnvironment(); int rank = 0, ranks = 1;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank); MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
	try {
		Require(argc == 2 && ranks == 3, "requires three MPI ranks and a fresh output directory");
		for (int round = 0; round < 2; ++round) {
			const int color = round == 0 || rank == 0 ? 0 : 1;
			MPI_Comm comm; MPI_Comm_split(PETSC_COMM_WORLD, color, rank, &comm);
			int local_rank = 0, local_size = 1; MPI_Comm_rank(comm, &local_rank); MPI_Comm_size(comm, &local_size);
			const auto root = fs::path(argv[1])/(std::to_string(round)+"-"+std::to_string(color));
			iga::CollectiveLocalStage(comm, "checkpoint fixture", [&] {
				if (local_rank == 0) { Require(fs::create_directories(root), "fixture exists"); WriteThreeDCase(root); WriteDatabase(root/"group.ntiga", local_size); }
			});
			for (bool transient : {false, true}) Run(root, comm, transient, false);
			Run(root, comm, true, true); MPI_Comm_free(&comm);
		}
		if (rank == 0) std::cout << "body_fitted_accepted_checkpoint all comparisons passed\n";
	} catch (const std::exception& error) { std::cerr << "rank " << rank << ": " << error.what() << '\n'; MPI_Abort(PETSC_COMM_WORLD, 1); return 1; }
	iga::CurrentPhaseProfile().Write(std::cout, rank, ranks, 0);
	PetscFinalize(); return 0;
}
