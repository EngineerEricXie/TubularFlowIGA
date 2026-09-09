#include "OneDCheckpoint.hpp"

#include <iostream>

namespace {

void Require(bool condition, const std::string& message)
{
	if (!condition) throw std::runtime_error(message);
}

void Compare(const iga::OneDFlowState& a, const iga::OneDFlowState& b)
{
	const auto compare = [](const std::vector<double>& x, const std::vector<double>& y) {
		Require(x.size() == y.size(), "implicit field size differs");
		double norm = 0.0, error = 0.0;
		for (std::size_t i = 0; i < x.size(); ++i) {
			Require(std::isfinite(x[i]) && std::isfinite(y[i]), "nonfinite implicit field");
			norm = std::hypot(norm, x[i]); error = std::hypot(error, x[i]-y[i]);
		}
		Require(norm == 0.0 ? error <= 1e-12 : error/norm <= 1e-6, "implicit retry differs from serial field");
	};
	compare(a.area, b.area); compare(a.flow, b.flow); compare(a.pressure, b.pressure);
	compare(a.node_pressure, b.node_pressure); compare(a.segment_flow, b.segment_flow);
	Require(a.outlets.size() == b.outlets.size(), "outlet count differs");
	for (std::size_t i = 0; i < a.outlets.size(); ++i) {
		compare({a.outlets[i].pressure}, {b.outlets[i].pressure});
		compare({a.outlets[i].flow}, {b.outlets[i].flow});
		compare({a.outlets[i].capacitor_pressure}, {b.outlets[i].capacitor_pressure});
	}
}

void Option(const char* name, const char* value)
{
	iga::OneDPetscCheck(value ? PetscOptionsSetValue(nullptr, name, value)
		: PetscOptionsClearValue(nullptr, name), "test option");
}

void CheckCallbacks(MPI_Comm group, const iga::OneDNetwork& network,
	const iga::OneDFlowSystemDefinition& flow, iga::OneDFlowState state)
{
	int rank = 0, ranks = 1;
	MPI_Comm_rank(group, &rank); MPI_Comm_size(group, &ranks);
	const int target = ranks-1;
	const auto graph = iga::BuildOneDImplicitGraph(network, false);
	const PetscInt unknowns = graph.nodes+static_cast<int>(graph.edges.size());
	for (const bool jacobian_fault : {false, true}) {
		iga::OneDPetscObjects objects;
		iga::OneDCollectivePetscCheck(group, VecCreateMPI(group, PETSC_DECIDE, unknowns, &objects.solution), "test vector");
		iga::OneDCollectivePetscCheck(group, VecDuplicate(objects.solution, &objects.rhs), "test residual");
		iga::OneDCollectivePetscCheck(group, MatCreateAIJ(group, PETSC_DECIDE, PETSC_DECIDE,
			unknowns, unknowns, 10, nullptr, 10, nullptr, &objects.matrix), "test matrix");
		iga::OneDNonlinearContext healthy;
		healthy.network = &network; healthy.graph = &graph; healthy.flow = &flow; healthy.state = &state;
		healthy.dt = 1e-3; healthy.inlet_flow = 1e-9;
		healthy.old_pressure = state.node_pressure;
		healthy.old_flow.assign(graph.edges.size(), 0.0);
		healthy.compliance = iga::OneDNodeCompliance(graph, flow, healthy.old_pressure);
		auto broken = healthy;
		if (rank == target) broken.old_flow.clear();
		std::vector<double> initial(static_cast<std::size_t>(unknowns), 0.0);
		std::copy(healthy.old_pressure.begin(), healthy.old_pressure.end(), initial.begin());
		iga::OneDSetInitialVector(objects.solution, initial);
		iga::OneDCollectivePetscCheck(group, SNESCreate(group, &objects.snes), "test SNES");
		iga::OneDCollectivePetscCheck(group, SNESSetFunction(objects.snes, objects.rhs,
			iga::OneDNonlinearResidual, jacobian_fault ? &healthy : &broken), "test function");
		iga::OneDCollectivePetscCheck(group, SNESSetJacobian(objects.snes, objects.matrix,
			objects.matrix, iga::OneDNonlinearJacobian, jacobian_fault ? &broken : &healthy), "test jacobian");
		// Exercise PETSc's C callback stack, not a direct C++ invocation.
		const auto error = SNESSolve(objects.snes, nullptr, objects.solution);
		Require(error != 0 && broken.callback_error, "SNES callback error escaped or was lost");
		bool diagnostic = false;
		try { std::rethrow_exception(broken.callback_error); }
		catch (const std::runtime_error& failure) {
			const auto expected = std::string("1d nonlinear ")+(jacobian_fault ? "jacobian" : "residual")
				+" assembly: rank "+std::to_string(target)+":";
			diagnostic = std::string(failure.what()).find(expected) != std::string::npos;
		}
		Require(diagnostic, "SNES callback diagnostic is not common to the group");
		objects.Close(group);
	}
	iga::OneDPetscObjects objects;
	iga::OneDCollectivePetscCheck(group, VecCreateMPI(group, PETSC_DECIDE, unknowns, &objects.solution), "test initialization vector");
	std::vector<double> values(static_cast<std::size_t>(unknowns), 2.0);
	if (rank == target) values.pop_back();
	bool rejected = false;
	try { iga::OneDSetInitialVector(objects.solution, values); }
	catch (const std::runtime_error& error) {
		rejected = std::string(error.what()).find("1d implicit vector initialization: rank "
			+std::to_string(target)+":") != std::string::npos;
	}
	Require(rejected, "initial vector size error was not coordinated");
	values.assign(static_cast<std::size_t>(unknowns), 2.0);
	iga::OneDSetInitialVector(objects.solution, values);
	std::vector<double> restored;
	iga::OneDGetVectorAll(objects.solution, restored);
	Require(restored == values, "vector retry failed");
	objects.Close(group);
}

void CheckGroup(MPI_Comm group, const char* path, int repetitions)
{
	int rank = 0, ranks = 1;
	MPI_Comm_rank(group, &rank); MPI_Comm_size(group, &ranks);
	const int target = ranks-1;
	iga::OneDFlowSystemDefinition flow;
	flow.name = "failure-flow"; flow.model = iga::OneDFlowModel::Compliant;
	flow.scheme = iga::OneDFlowScheme::ImplicitPetsc;
	flow.dynamic_viscosity = 0.004; flow.density = 1060.0;
	const auto network = iga::ReadOneDNetwork(path, 1.0, 3, flow.dynamic_viscosity);
	iga::OneDFlowState initial;
	for (const auto node : network.outlet_nodes) {
		iga::OneDOutletState outlet; outlet.node = node;
		outlet.kind = iga::OneDOutletKind::Pressure; initial.outlets.push_back(outlet);
	}
	iga::InitializeCompliantOneDFromRigid(network, flow, initial, 1e-9, 1e-3);
	int failures = 0;
	for (int repeat = 0; repeat < repetitions; ++repeat) {
		for (const auto formulation : {iga::OneDImplicitFormulation::PressureNetwork,
			iga::OneDImplicitFormulation::LinearizedAQ, iga::OneDImplicitFormulation::NonlinearAQ,
			iga::OneDImplicitFormulation::ImplicitPde}) {
			flow.formulation = formulation;
			for (const std::string mode : {"step", "graph", "assembly", "publication", "layout", "method", "petsc-options", "convergence"}) {
				if (ranks == 1 && (mode == "layout" || mode == "method")) continue;
				auto local_flow = flow; auto local_network = network; auto current = initial;
				double dt = 1e-3;
				std::string expected;
				if (mode == "step") {
					if (rank == target) dt = 0.0;
					expected = "1d implicit graph preparation";
				} else if (mode == "graph") {
					if (rank == target) local_network.segments.front().parent = -1;
					expected = "1d implicit graph preparation";
				} else if (mode == "assembly") {
					if (rank == target) {
						local_flow.wall.model = iga::OneDWallModel::Olufsen;
						current.node_pressure.assign(network.nodes.size(), 1e12);
					}
					expected = formulation == iga::OneDImplicitFormulation::PressureNetwork ? "1d pressure assembly" : "1d linearized assembly";
				} else if (mode == "publication") {
					if (rank == target) current.outlets.front().node = -1;
					expected = "1d implicit state update";
				} else if (mode == "layout") {
					if (rank == target) local_network.nodes.push_back(local_network.nodes.back());
					expected = "1d implicit layout";
				} else if (mode == "method") {
					if (rank == target) local_flow.formulation = formulation == iga::OneDImplicitFormulation::PressureNetwork
						? iga::OneDImplicitFormulation::LinearizedAQ : iga::OneDImplicitFormulation::PressureNetwork;
					expected = "1d implicit formulation";
				} else if (mode == "petsc-options") {
					Option("-ksp_type", "no_such_iga_ksp"); expected = "1d solver options";
				} else {
					Option("-ksp_type", "gmres"); Option("-pc_type", "none"); Option("-ksp_max_it", "0");
					expected = "1d implicit convergence";
				}
				const auto before = iga::PackOneDCheckpointState(current, {}, local_network);
				const auto before_inlet = current.inlet_flow;
				const auto before_outlet_node = current.outlets.front().node;
				bool rejected = false;
				try { iga::AdvanceImplicitOneD(local_network, local_flow, current, 2e-9, dt, group); }
				catch (const std::runtime_error& error) {
					rejected = std::string(error.what()).find(expected+": rank ") != std::string::npos;
				}
				Require(rejected, "implicit failure not coordinated: "+mode);
				Require(iga::PackOneDCheckpointState(current, {}, local_network) == before, "failed implicit advance published state: "+mode);
				Require(current.inlet_flow == before_inlet && current.outlets.front().node == before_outlet_node
					&& current.completed_step == initial.completed_step && current.physical_time == initial.physical_time
					&& current.internal_substeps == initial.internal_substeps, "failed implicit advance changed metadata: "+mode);
				Option("-ksp_type", "preonly"); Option("-pc_type", "lu"); Option("-ksp_max_it", nullptr);
				current = initial;
				auto reference = initial;
				iga::AdvanceImplicitOneD(network, flow, current, 2e-9, 1e-3, group);
				iga::AdvanceImplicitOneD(network, flow, reference, 2e-9, 1e-3, PETSC_COMM_SELF);
				Compare(reference, current);
				++failures;
			}
		}
		CheckCallbacks(group, network, flow, initial);
	}
	if (rank == 0) std::cout << "implicit failure group ranks=" << ranks << " repetitions=" << repetitions
		<< " solve_failures=" << failures << " callback_failures=" << 2*repetitions
		<< " vector_failures=" << repetitions << " passed\n";
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	PetscPushErrorHandler(PetscReturnErrorHandler, nullptr);
	int rank = 0, ranks = 1;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank); MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
	try {
		Require(argc == 2 && ranks == 3, "require three ranks and a valid SWC");
		CheckGroup(PETSC_COMM_WORLD, argv[1], 1);
		MPI_Comm group = MPI_COMM_NULL;
		MPI_Comm_split(PETSC_COMM_WORLD, rank == 0 ? 0 : 1, rank, &group);
		CheckGroup(group, argv[1], rank == 0 ? 1 : 2);
		MPI_Comm_free(&group);
		std::cout << "implicit failure rank=" << rank << " passed\n";
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		MPI_Abort(PETSC_COMM_WORLD, 2);
	}
	PetscPopErrorHandler();
	PetscFinalize();
	return 0;
}
