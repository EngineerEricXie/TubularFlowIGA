#include "BoundarySupport.hpp"
#include "TransientFlowRuntime.hpp"
#include "TransientTransportRuntime.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
namespace fs = std::filesystem;

double Compare(const std::vector<double>& reference, const std::vector<double>& current)
{
	if (reference.size() != current.size()) throw std::runtime_error("field sizes differ");
	double norm = 0.0, error = 0.0;
	for (std::size_t i = 0; i < reference.size(); ++i) {
		if (!std::isfinite(reference[i]) || !std::isfinite(current[i]))
			throw std::runtime_error("nonfinite field");
		norm = std::hypot(norm, reference[i]);
		error = std::hypot(error, current[i]-reference[i]);
	}
	if (norm == 0.0 ? error > 1e-12 : error/norm > 1e-6)
		throw std::runtime_error("subcommunicator differs from independent serial reference");
	return norm > 0.0 ? error/norm : error;
}

std::vector<double> Gather(Vec state)
{
	Vec all = nullptr;
	VecScatter scatter = nullptr;
	VecScatterCreateToAll(state, &scatter, &all);
	VecScatterBegin(scatter, state, all, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(scatter, state, all, INSERT_VALUES, SCATTER_FORWARD);
	PetscInt count = 0;
	VecGetSize(all, &count);
	const PetscScalar* values = nullptr;
	VecGetArrayRead(all, &values);
	std::vector<double> result(static_cast<std::size_t>(count));
	for (PetscInt i = 0; i < count; ++i) result[i] = PetscRealPart(values[i]);
	VecRestoreArrayRead(all, &values);
	VecScatterDestroy(&scatter); VecDestroy(&all);
	return result;
}

std::vector<double> Component(const std::vector<double>& state, bool pressure)
{
	std::vector<double> result;
	for (std::size_t i = 0; i < state.size(); ++i)
		if ((i%4 == 3) == pressure) result.push_back(state[i]);
	return result;
}

struct Result {
	std::vector<double> flow, flow_norms, reference_inlet, transport, mass, sources;
};

void CheckFlowHalo(iga::TransientFlowRuntime& flow)
{
	Vec saved = nullptr;
	VecDuplicate(flow.State(), &saved);
	VecCopy(flow.State(), saved);
	const auto original = Gather(flow.State());
	PetscInt begin = 0, end = 0;
	VecGetOwnershipRange(flow.State(), &begin, &end);
	for (int epoch = 1; epoch <= 2; ++epoch) {
		PetscScalar* values = nullptr;
		VecGetArray(flow.State(), &values);
		for (PetscInt row = begin; row < end; ++row)
			values[row-begin] = 100000.0*epoch+row+1;
		VecRestoreArray(flow.State(), &values);
		const auto halo = flow.GatherRequiredVelocity();
		for (std::size_t node = 0; node < flow.RequiredNodes().size(); ++node)
			for (int component = 0; component < 3; ++component)
				if (halo.at(node)[component] != 100000.0*epoch+4.0*flow.RequiredNodes()[node]+component+1)
					throw std::runtime_error("required flow halo has wrong ID or stale state");
	}
	VecCopy(saved, flow.State());
	VecDestroy(&saved);
	const auto restored = flow.GatherRequiredVelocity();
	for (std::size_t node = 0; node < flow.RequiredNodes().size(); ++node)
		for (int component = 0; component < 3; ++component)
			if (restored.at(node)[component] != original.at(4*flow.RequiredNodes()[node]+component))
				throw std::runtime_error("required flow halo did not refresh restored state");
}

void CheckTransportHalo(const iga::TransientTransportRuntime& transport)
{
	const auto global = transport.GatherState();
	const auto required = transport.GatherRequiredState();
	std::vector<double> reference;
	const auto fields = transport.System().fields.size();
	for (auto node : transport.RequiredNodes())
		for (std::size_t field = 0; field < fields; ++field)
			reference.push_back(global.at(static_cast<std::size_t>(node)*fields+field));
	Compare(reference, required);
}

Result Run(const fs::path& database_path, const fs::path& flow_case,
	const fs::path& transport_case, MPI_Comm communicator, int color,
	const fs::path& checkpoint)
{
	iga::Database database(database_path.string());
	const auto mesh = iga::ReadLabeledHexMesh((flow_case/"controlmesh.vtk").string(),
		database.header().nodes, database.header().elements);
	const auto velocity = iga::ReadVelocity((flow_case/"initial_velocityfield.txt").string(),
		database.header().nodes);
	auto configuration = iga::ReadSimulationConfiguration((flow_case/"simulation_config.json").string());
	for (auto& boundary : configuration.boundaries)
		for (auto& condition : boundary.conditions)
			if (boundary.label == 1 && condition.field == "velocity") condition.scale *= color+1;
	const auto& definition = iga::FirstNavierStokesSystem(configuration);
	const auto boundaries = iga::ResolveFlowBoundaries(configuration, definition, mesh.labels, velocity);
	Result result;
	{
		iga::TransientFlowRuntime flow(database, communicator, true, false,
			{definition.density, definition.viscosity, 0.0}, boundaries, mesh.labels,
			velocity, iga::WallTraceBasis(database, mesh), {});
		flow.Assembler().ValidateOwnership();
		CheckFlowHalo(flow);
		iga::RequireValidGeometry(flow.Elements(),
			[&](const iga::Element& element) { return flow.Assembler().OwnsElementByMinimumNode(element); },
			communicator);
		flow.Advance(configuration, 0, 0.0, 12, 1e-8, 1e-12, 1e-6);
		result.flow = Gather(flow.State());
		const auto summary = flow.Summary();
		result.flow_norms = {summary.velocity_l2, summary.pressure_l2, summary.state_l2};
		result.reference_inlet = {flow.ReferenceBoundaryFlow(1)};
	}
	auto transport_configuration = iga::ReadSimulationConfiguration(
		(transport_case/"simulation_config.json").string());
	for (auto& boundary : transport_configuration.boundaries)
		for (auto& condition : boundary.conditions)
			if (condition.kind == iga::FieldBoundaryKind::Dirichlet)
				for (auto& value : condition.value) value *= color+1;
	const auto system = iga::CompileLinearSystem(transport_configuration, "neuron_transport");
	const auto prescribed = iga::ReadVelocity((transport_case/"initial_velocityfield.txt").string(),
		database.header().nodes);
	{
		iga::TransientTransportRuntime transport(database, communicator, transport_configuration, system, mesh.labels);
		std::vector<std::array<double, 3>> required_velocity;
		for (auto node : transport.RequiredNodes()) required_velocity.push_back(prescribed.at(node));
		std::vector<double> first_step;
		for (int step = 0; step < color+1; ++step) {
			transport.BeginStep();
			transport.SolveTrial(transport_configuration, transport.RequiredNodes(), required_velocity);
			CheckTransportHalo(transport);
			const auto first_trial = transport.GatherState();
			transport.RollbackTrial();
			CheckTransportHalo(transport);
			transport.SolveTrial(transport_configuration, transport.RequiredNodes(), required_velocity);
			Compare(first_trial, transport.GatherState());
			transport.CommitStep();
			CheckTransportHalo(transport);
			if (step == 0) {
				first_step = transport.GatherState();
				if (!checkpoint.empty()) transport.WriteState(checkpoint);
			}
		}
		result.transport = transport.GatherState();
		for (const auto& mass : transport.TotalMass()) result.mass.push_back(mass.second);
		for (const auto& source : transport.SourceIntegrals()) result.sources.push_back(source.second);
		if (!checkpoint.empty()) {
			transport.ReadState(checkpoint);
			CheckTransportHalo(transport);
			Compare(first_step, transport.GatherState());
			if (transport.Steps() != 1) throw std::runtime_error("checkpoint step mismatch");
		}
	}
	return result;
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int rank = 0, ranks = 1;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank); MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
	MPI_Comm group = MPI_COMM_NULL;
	try {
		if (argc != 6 || ranks != 3)
			throw std::runtime_error("usage: mpiexec -np 3 body_fitted_subcommunicator_test DB1 DB2 FLOW_CASE TRANSPORT_CASE OUTPUT_DIR");
		const int color = rank == 0 ? 0 : 1;
		MPI_Comm_split(PETSC_COMM_WORLD, color, rank, &group);
		int local_rank = 0;
		MPI_Comm_rank(group, &local_rank);
		const auto reference = Run(argv[1], argv[3], argv[4], PETSC_COMM_SELF, color, {});
		const auto checkpoint = fs::path(argv[5])/("group-"+std::to_string(color)+".petsc");
		const auto current = Run(argv[color+1], argv[3], argv[4], group, color, checkpoint);
		const double velocity = Compare(Component(reference.flow, false), Component(current.flow, false));
		const double pressure = Compare(Component(reference.flow, true), Component(current.flow, true));
		const double transport = Compare(reference.transport, current.transport);
		const double mass = Compare(reference.mass, current.mass);
		Compare(reference.flow_norms, current.flow_norms);
		Compare(reference.reference_inlet, current.reference_inlet);
		Compare(reference.sources, current.sources);
		if (local_rank == 0) std::cout << "body_fitted_subcommunicator group=" << color
			<< " velocity_error=" << velocity << " pressure_error=" << pressure
			<< " transport_error=" << transport << " mass_error=" << mass << " passed\n";
		// All runtimes and their PETSc objects have died before the borrowed
		// communicator is freed, and this happens before PetscFinalize.
		MPI_Comm_free(&group);
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		MPI_Abort(PETSC_COMM_WORLD, 1);
		return 1;
	}
	PetscFinalize();
	return 0;
}
