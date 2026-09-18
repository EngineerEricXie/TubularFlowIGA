#ifndef IGA_ONE_D_IMPLICIT_HPP
#define IGA_ONE_D_IMPLICIT_HPP

#include "OneDFlow.hpp"
#include "OneDPetscSupport.hpp"

#include <petscksp.h>
#include <petscsnes.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace iga {

struct OneDImplicitEdge {
	int from = -1;
	int to = -1;
	int segment = -1;
	int cell = 0;
	double length = 0.0;
	double area0 = 0.0;
	double radius0 = 0.0;
};

struct OneDImplicitGraph {
	int original_nodes = 0;
	int root = -1;
	std::vector<OneDImplicitEdge> edges;
	std::vector<int> outlet_graph_nodes;
	std::vector<int> outlet_state_indices;
	int nodes = 0;
};

inline OneDImplicitGraph BuildOneDImplicitGraph(const OneDNetwork& network,
	bool expand_cells)
{
	const auto limit = static_cast<std::size_t>(std::numeric_limits<int>::max());
	if (network.nodes.empty() || network.nodes.size() > limit
		|| network.root < 0 || static_cast<std::size_t>(network.root) >= network.nodes.size())
		throw std::runtime_error("invalid implicit network node layout");
	std::size_t expanded_nodes = network.nodes.size(), edges = 0;
	for (const auto& segment : network.segments) {
		if (segment.parent < 0 || segment.child < 0
			|| static_cast<std::size_t>(segment.parent) >= network.nodes.size()
			|| static_cast<std::size_t>(segment.child) >= network.nodes.size()
			|| segment.cells <= 0)
			throw std::runtime_error("invalid implicit segment layout");
		const auto pieces = static_cast<std::size_t>(expand_cells ? segment.cells : 1);
		if (pieces > limit-edges || pieces-1 > limit-expanded_nodes)
			throw std::runtime_error("implicit network exceeds supported index range");
		edges += pieces;
		expanded_nodes += pieces-1;
	}
	if (edges > limit-expanded_nodes)
		throw std::runtime_error("implicit unknown count exceeds supported index range");
	OneDImplicitGraph graph;
	graph.original_nodes = static_cast<int>(network.nodes.size());
	graph.nodes = graph.original_nodes;
	graph.root = network.root;
	for (const auto& segment : network.segments) {
		const int pieces = expand_cells ? segment.cells : 1;
		int from = segment.parent;
		for (int cell = 0; cell < pieces; ++cell) {
			const int to = cell == pieces-1 ? segment.child : graph.nodes++;
			graph.edges.push_back({from, to, segment.index, cell,
				segment.length/pieces, segment.area0, segment.radius0});
			from = to;
		}
	}
	for (std::size_t i = 0; i < network.outlet_nodes.size(); ++i) {
		graph.outlet_graph_nodes.push_back(network.outlet_nodes[i]);
		graph.outlet_state_indices.push_back(static_cast<int>(i));
	}
	return graph;
}

inline std::vector<double> OneDNodeCompliance(const OneDImplicitGraph& graph,
	const OneDNetwork& network, const OneDFlowSystemDefinition& flow,
	const std::vector<double>& pressure)
{
	std::vector<double> result(static_cast<std::size_t>(graph.nodes), 0.0);
	if (flow.formulation == OneDImplicitFormulation::SteadyR) return result;
	for (const auto& edge : graph.edges) {
		const auto& segment = network.segments[static_cast<std::size_t>(edge.segment)];
		const int child_id = network.nodes[static_cast<std::size_t>(segment.child)].id;
		const auto configured = flow.lumped.segment_compliance.find(child_id);
		double compliance = 0.0;
		if (configured != flow.lumped.segment_compliance.end())
			compliance = configured->second*edge.length/segment.length;
		else {
			const double mean_pressure = 0.5*(pressure[static_cast<std::size_t>(edge.from)]
				+pressure[static_cast<std::size_t>(edge.to)]);
			compliance = flow.lumped.compliance_scale*OneDWallAreaDerivative(mean_pressure,
				edge.area0, edge.radius0, flow.wall)*edge.length;
		}
		result[static_cast<std::size_t>(edge.from)] += 0.5*compliance;
		result[static_cast<std::size_t>(edge.to)] += 0.5*compliance;
	}
	return result;
}

inline const OneDOutletState* OneDOutletAtGraphNode(const OneDImplicitGraph& graph,
	const std::vector<OneDOutletState>& outlets, int node)
{
	for (std::size_t i = 0; i < graph.outlet_graph_nodes.size(); ++i)
		if (graph.outlet_graph_nodes[i] == node)
			return &outlets.at(static_cast<std::size_t>(graph.outlet_state_indices[i]));
	return nullptr;
}

inline void OneDGetVectorAll(Vec vector, std::vector<double>& values)
{
	const auto communicator = PetscObjectComm(reinterpret_cast<PetscObject>(vector));
	OneDPetscObjects objects;
	auto& all = objects.solution;
	auto& scatter = objects.scatter;
	OneDCollectivePetscCheck(communicator, VecScatterCreateToAll(vector, &scatter, &all), "VecScatterCreateToAll");
	OneDCollectivePetscCheck(communicator, VecScatterBegin(scatter, vector, all, INSERT_VALUES, SCATTER_FORWARD), "VecScatterBegin");
	OneDCollectivePetscCheck(communicator, VecScatterEnd(scatter, vector, all, INSERT_VALUES, SCATTER_FORWARD), "VecScatterEnd");
	std::vector<double> candidate;
	CollectiveLocalStage(communicator, "1d implicit vector read", [&] {
		PetscInt size = 0;
		OneDPetscCheck(VecGetSize(all, &size), "VecGetSize");
		candidate.resize(static_cast<std::size_t>(size));
		PetscReadArray array;
		array.Acquire(all);
		for (PetscInt i = 0; i < size; ++i)
			candidate[static_cast<std::size_t>(i)] = PetscRealPart(array.Data()[i]);
		array.Restore();
	});
	objects.Close(communicator);
	values.swap(candidate);
}

inline void OneDSetInitialVector(Vec vector, const std::vector<double>& values)
{
	const auto communicator = PetscObjectComm(reinterpret_cast<PetscObject>(vector));
	CollectiveLocalStage(communicator, "1d implicit vector initialization", [&] {
		PetscInt first = 0, last = 0, size = 0;
		OneDPetscCheck(VecGetSize(vector, &size), "VecGetSize");
		if (values.size() != static_cast<std::size_t>(size))
			throw std::runtime_error("initial vector size does not match distributed vector");
		OneDPetscCheck(VecGetOwnershipRange(vector, &first, &last), "VecGetOwnershipRange");
		for (PetscInt i = first; i < last; ++i)
			OneDPetscCheck(VecSetValue(vector, i, values[static_cast<std::size_t>(i)], INSERT_VALUES), "VecSetValue");
	});
	OneDCollectivePetscCheck(communicator, VecAssemblyBegin(vector), "VecAssemblyBegin");
	OneDCollectivePetscCheck(communicator, VecAssemblyEnd(vector), "VecAssemblyEnd");
}

inline void OneDUpdateStateFromImplicitSolution(const OneDNetwork& network,
	const OneDImplicitGraph& graph, const OneDFlowSystemDefinition& flow,
	const std::vector<double>& pressure, const std::vector<double>& edge_flow,
	OneDFlowState& state, double dt)
{
	state.node_pressure.assign(pressure.begin(), pressure.begin()+graph.original_nodes);
	state.area.resize(static_cast<std::size_t>(network.cells));
	state.flow.resize(static_cast<std::size_t>(network.cells));
	state.pressure.resize(static_cast<std::size_t>(network.cells));
	state.segment_flow.assign(network.segments.size(), 0.0);
	std::vector<int> counts(network.segments.size(), 0);
	for (std::size_t edge_index = 0; edge_index < graph.edges.size(); ++edge_index) {
		const auto& edge = graph.edges[edge_index];
		const auto& segment = network.segments[static_cast<std::size_t>(edge.segment)];
		const double mean_pressure = 0.5*(pressure[static_cast<std::size_t>(edge.from)]
			+pressure[static_cast<std::size_t>(edge.to)]);
		if (graph.edges.size() == network.segments.size()) {
			for (int cell = 0; cell < segment.cells; ++cell) {
				const auto index = static_cast<std::size_t>(segment.cell_offset+cell);
				state.pressure[index] = mean_pressure;
				state.area[index] = OneDAreaFromPressure(mean_pressure, segment.area0, segment.radius0, flow.wall);
				state.flow[index] = edge_flow[edge_index];
			}
		} else {
			const auto index = static_cast<std::size_t>(segment.cell_offset+edge.cell);
			state.pressure[index] = mean_pressure;
			state.area[index] = OneDAreaFromPressure(mean_pressure, segment.area0, segment.radius0, flow.wall);
			state.flow[index] = edge_flow[edge_index];
		}
		state.segment_flow[static_cast<std::size_t>(edge.segment)] += edge_flow[edge_index];
		++counts[static_cast<std::size_t>(edge.segment)];
	}
	for (std::size_t i = 0; i < state.segment_flow.size(); ++i)
		state.segment_flow[i] /= std::max(counts[i], 1);
	for (auto& outlet : state.outlets) {
		const int incoming = OneDSegmentIntoNode(network, outlet.node);
		if (incoming < 0) throw std::runtime_error("implicit outlet has no incoming segment");
		AdvanceOutletState(outlet, state.segment_flow[static_cast<std::size_t>(incoming)], dt);
	}
}

// Implicit solver calls borrow communicator for the duration of the call.
// The caller owns/free it; these functions neither duplicate nor retain it.
// The default preserves existing CLI and injected callback behavior.
inline void SolveOneDPressureNetworkPetsc(const OneDNetwork& network,
	const OneDFlowSystemDefinition& flow, OneDFlowState& state,
	double inlet_flow, double dt, MPI_Comm communicator = PETSC_COMM_WORLD,
	OneDPetscSolverContext* solver_context = nullptr)
{
	std::optional<OneDPetscSolverContext> fallback;
	if (!solver_context) { fallback.emplace(communicator); solver_context = &*fallback; }
	OneDImplicitGraph graph;
	CollectiveLocalStage(communicator, "1d implicit graph preparation", [&] {
		graph = BuildOneDImplicitGraph(network, false);
		if (!std::isfinite(dt) || dt <= 0.0 || !std::isfinite(inlet_flow))
			throw std::runtime_error("implicit step requires finite flow and positive finite dt");
	});
	const PetscInt n = graph.nodes;
	OneDRequireSameInteger(communicator, n, "1d implicit layout");
	OneDPetscObjects objects;
	auto& matrix = objects.matrix;
	auto& rhs = objects.rhs;
	auto& solution = objects.solution;
	auto& solver = objects.ksp;
	OneDCollectivePetscCheck(communicator, MatCreateAIJ(communicator, PETSC_DECIDE, PETSC_DECIDE, n, n,
		8, nullptr, 8, nullptr, &matrix), "MatCreateAIJ");
	OneDCollectivePetscCheck(communicator, VecCreateMPI(communicator, PETSC_DECIDE, n, &rhs), "VecCreateMPI");
	OneDCollectivePetscCheck(communicator, VecDuplicate(rhs, &solution), "VecDuplicate");
	std::vector<double> previous;
	std::vector<double> compliance;
	CollectiveLocalStage(communicator, "1d pressure assembly", [&] {
		PetscInt first = 0, last = 0;
		OneDPetscCheck(MatGetOwnershipRange(matrix, &first, &last), "MatGetOwnershipRange");
		previous = state.node_pressure;
		if (previous.size() != static_cast<std::size_t>(n)) previous.assign(static_cast<std::size_t>(n), flow.wall.reference_pressure);
		compliance = OneDNodeCompliance(graph, network, flow, previous);
		for (PetscInt row = first; row < last; ++row) {
			const auto* outlet = OneDOutletAtGraphNode(graph, state.outlets, static_cast<int>(row));
			if (outlet && outlet->kind == OneDOutletKind::Pressure) {
				OneDPetscCheck(MatSetValue(matrix, row, row, 1.0, ADD_VALUES), "MatSetValue");
				OneDPetscCheck(VecSetValue(rhs, row, outlet->pressure, INSERT_VALUES), "VecSetValue");
				continue;
			}
			double diagonal = compliance[static_cast<std::size_t>(row)]/dt;
			double value = diagonal*previous[static_cast<std::size_t>(row)];
			for (const auto& edge : graph.edges) {
				if (edge.from != row && edge.to != row) continue;
				const auto& segment = network.segments[static_cast<std::size_t>(edge.segment)];
				const double resistance = flow.model == OneDFlowModel::Lumped
					? OneDLumpedSegmentResistance(network, flow, segment)
					: 8.0*OneDPi*flow.dynamic_viscosity*edge.length/(edge.area0*edge.area0);
				const int other = edge.from == row ? edge.to : edge.from;
				diagonal += 1.0/resistance;
				const auto* other_outlet = OneDOutletAtGraphNode(graph, state.outlets, other);
				if (other_outlet && other_outlet->kind == OneDOutletKind::Pressure)
					value += other_outlet->pressure/resistance;
				else OneDPetscCheck(MatSetValue(matrix, row, other,
					-1.0/resistance, ADD_VALUES), "MatSetValue");
			}
			if (row == graph.root) value += inlet_flow;
			if (outlet) {
				const double resistance = OutletEffectiveResistance(*outlet, dt);
				const double pressure = OutletEffectivePressure(*outlet, dt);
				diagonal += 1.0/resistance;
				value += pressure/resistance;
			}
			OneDPetscCheck(MatSetValue(matrix, row, row, diagonal, ADD_VALUES), "MatSetValue");
			OneDPetscCheck(VecSetValue(rhs, row, value, INSERT_VALUES), "VecSetValue");
		}
	});
	OneDCollectivePetscCheck(communicator, MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY), "MatAssemblyBegin");
	OneDCollectivePetscCheck(communicator, MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY), "MatAssemblyEnd");
	OneDCollectivePetscCheck(communicator, VecAssemblyBegin(rhs), "VecAssemblyBegin");
	OneDCollectivePetscCheck(communicator, VecAssemblyEnd(rhs), "VecAssemblyEnd");
	OneDCollectivePetscCheck(communicator, KSPCreate(communicator, &solver), "KSPCreate");
	OneDCollectivePetscCheck(communicator, KSPSetOperators(solver, matrix, matrix), "KSPSetOperators");
	solver_context->options.Attach(solver);
	solver_context->options.Call("1d solver options", [&] { return KSPSetFromOptions(solver); });
	RequireKspFactorBackend(solver, matrix, communicator);
	solver_context->options.Call("1d linear solve", [&] { return KSPSolve(solver, rhs, solution); });
	solver_context->Record(communicator, solver);
	CollectiveLocalStage(communicator, "1d implicit convergence", [&] {
		KSPConvergedReason reason;
		OneDPetscCheck(KSPGetConvergedReason(solver, &reason), "KSPGetConvergedReason");
		if (reason <= 0) throw std::runtime_error("PETSc pressure_network did not converge");
	});
	std::vector<double> pressure;
	OneDGetVectorAll(solution, pressure);
	OneDFlowState candidate;
	CollectiveLocalStage(communicator, "1d implicit state update", [&] {
		candidate = state;
		std::vector<double> edge_flow(graph.edges.size());
		for (std::size_t i = 0; i < graph.edges.size(); ++i) {
			const auto& edge = graph.edges[i];
			const auto& segment = network.segments[static_cast<std::size_t>(edge.segment)];
			const double resistance = flow.model == OneDFlowModel::Lumped
				? OneDLumpedSegmentResistance(network, flow, segment)
				: 8.0*OneDPi*flow.dynamic_viscosity*edge.length/(edge.area0*edge.area0);
			edge_flow[i] = (pressure[static_cast<std::size_t>(edge.from)]-pressure[static_cast<std::size_t>(edge.to)])/resistance;
		}
		OneDUpdateStateFromImplicitSolution(network, graph, flow, pressure, edge_flow, candidate, dt);
		double storage_rate = 0.0;
		for (int node = 0; node < graph.nodes; ++node) {
			const auto* outlet = OneDOutletAtGraphNode(graph, candidate.outlets, node);
			if (outlet && outlet->kind == OneDOutletKind::Pressure) continue;
			storage_rate += compliance[static_cast<std::size_t>(node)]
				*(pressure[static_cast<std::size_t>(node)]-previous[static_cast<std::size_t>(node)])/dt;
		}
		double outlet_flow = 0.0;
		for (const auto& outlet : candidate.outlets) outlet_flow += outlet.flow;
		const double residual = inlet_flow-outlet_flow-storage_rate;
		candidate.storage_rate = storage_rate;
		candidate.relative_continuity_residual = std::abs(residual)
			/std::max({std::abs(inlet_flow), std::abs(outlet_flow),
				std::abs(storage_rate), 1.0e-30});
		candidate.has_conservation_diagnostic = true;
	});
	objects.Close(communicator);
	static_assert(std::is_nothrow_move_assignable<OneDFlowState>::value, "implicit state publication must not throw");
	state = std::move(candidate);
}

inline void SolveOneDLinearizedAQPetsc(const OneDNetwork& network,
	const OneDFlowSystemDefinition& flow, OneDFlowState& state,
	double inlet_flow, double dt, bool expand_cells, MPI_Comm communicator = PETSC_COMM_WORLD,
	OneDPetscSolverContext* solver_context = nullptr)
{
	std::optional<OneDPetscSolverContext> fallback;
	if (!solver_context) { fallback.emplace(communicator); solver_context = &*fallback; }
	OneDImplicitGraph graph;
	CollectiveLocalStage(communicator, "1d implicit graph preparation", [&] {
		graph = BuildOneDImplicitGraph(network, expand_cells);
		if (!std::isfinite(dt) || dt <= 0.0 || !std::isfinite(inlet_flow))
			throw std::runtime_error("implicit step requires finite flow and positive finite dt");
	});
	const PetscInt unknowns = graph.nodes+static_cast<int>(graph.edges.size());
	OneDRequireSameInteger(communicator, unknowns, "1d implicit layout");
	OneDPetscObjects objects;
	auto& matrix = objects.matrix;
	auto& rhs = objects.rhs;
	auto& solution = objects.solution;
	auto& solver = objects.ksp;
	OneDCollectivePetscCheck(communicator, MatCreateAIJ(communicator, PETSC_DECIDE, PETSC_DECIDE, unknowns, unknowns,
		10, nullptr, 10, nullptr, &matrix), "MatCreateAIJ");
	OneDCollectivePetscCheck(communicator, VecCreateMPI(communicator, PETSC_DECIDE, unknowns, &rhs), "VecCreateMPI");
	OneDCollectivePetscCheck(communicator, VecDuplicate(rhs, &solution), "VecDuplicate");
	CollectiveLocalStage(communicator, "1d linearized assembly", [&] {
		std::vector<double> previous_pressure(static_cast<std::size_t>(graph.nodes), flow.wall.reference_pressure);
		for (int i = 0; i < graph.original_nodes && i < static_cast<int>(state.node_pressure.size()); ++i)
			previous_pressure[static_cast<std::size_t>(i)] = state.node_pressure[static_cast<std::size_t>(i)];
		const auto compliance = OneDNodeCompliance(graph, network, flow, previous_pressure);
		PetscInt first = 0, last = 0;
		OneDPetscCheck(MatGetOwnershipRange(matrix, &first, &last), "MatGetOwnershipRange");
		for (PetscInt row = first; row < last; ++row) {
			if (row < graph.nodes) {
				const auto* outlet = OneDOutletAtGraphNode(graph, state.outlets, static_cast<int>(row));
				if (outlet && outlet->kind == OneDOutletKind::Pressure) {
					OneDPetscCheck(MatSetValue(matrix, row, row, 1.0, ADD_VALUES), "MatSetValue");
					OneDPetscCheck(VecSetValue(rhs, row, outlet->pressure, INSERT_VALUES), "VecSetValue");
					continue;
				}
				double diagonal = compliance[static_cast<std::size_t>(row)]/dt;
				double value = diagonal*previous_pressure[static_cast<std::size_t>(row)];
				for (std::size_t edge_index = 0; edge_index < graph.edges.size(); ++edge_index) {
					const auto& edge = graph.edges[edge_index];
					if (edge.from == row) OneDPetscCheck(MatSetValue(matrix, row, graph.nodes+edge_index, 1.0, ADD_VALUES), "MatSetValue");
					if (edge.to == row) OneDPetscCheck(MatSetValue(matrix, row, graph.nodes+edge_index, -1.0, ADD_VALUES), "MatSetValue");
				}
				if (row == graph.root) value += inlet_flow;
				if (outlet) {
					const double resistance = OutletEffectiveResistance(*outlet, dt);
					const double pressure = OutletEffectivePressure(*outlet, dt);
					diagonal += 1.0/resistance;
					value += pressure/resistance;
				}
				OneDPetscCheck(MatSetValue(matrix, row, row, diagonal, ADD_VALUES), "MatSetValue");
				OneDPetscCheck(VecSetValue(rhs, row, value, INSERT_VALUES), "VecSetValue");
			} else {
				const int edge_index = static_cast<int>(row)-graph.nodes;
				const auto& edge = graph.edges[static_cast<std::size_t>(edge_index)];
				const double resistance = 8.0*OneDPi*flow.dynamic_viscosity*edge.length/(edge.area0*edge.area0);
				const double inertance = flow.density*edge.length/edge.area0;
				double previous_flow = 0.0;
				if (expand_cells) {
					const auto& segment = network.segments[static_cast<std::size_t>(edge.segment)];
					if (state.flow.size() == static_cast<std::size_t>(network.cells))
						previous_flow = state.flow[static_cast<std::size_t>(segment.cell_offset+edge.cell)];
				} else if (state.segment_flow.size() == network.segments.size())
					previous_flow = state.segment_flow[static_cast<std::size_t>(edge.segment)];
				OneDPetscCheck(MatSetValue(matrix, row, edge.from, -1.0, ADD_VALUES), "MatSetValue");
				OneDPetscCheck(MatSetValue(matrix, row, edge.to, 1.0, ADD_VALUES), "MatSetValue");
				OneDPetscCheck(MatSetValue(matrix, row, row, resistance+inertance/dt, ADD_VALUES), "MatSetValue");
				OneDPetscCheck(VecSetValue(rhs, row, inertance*previous_flow/dt, INSERT_VALUES), "VecSetValue");
			}
		}
	});
	OneDCollectivePetscCheck(communicator, MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY), "MatAssemblyBegin");
	OneDCollectivePetscCheck(communicator, MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY), "MatAssemblyEnd");
	OneDCollectivePetscCheck(communicator, VecAssemblyBegin(rhs), "VecAssemblyBegin");
	OneDCollectivePetscCheck(communicator, VecAssemblyEnd(rhs), "VecAssemblyEnd");
	OneDCollectivePetscCheck(communicator, KSPCreate(communicator, &solver), "KSPCreate");
	OneDCollectivePetscCheck(communicator, KSPSetOperators(solver, matrix, matrix), "KSPSetOperators");
	solver_context->options.Attach(solver);
	solver_context->options.Call("1d solver options", [&] { return KSPSetFromOptions(solver); });
	RequireKspFactorBackend(solver, matrix, communicator);
	solver_context->options.Call("1d linear solve", [&] { return KSPSolve(solver, rhs, solution); });
	solver_context->Record(communicator, solver);
	CollectiveLocalStage(communicator, "1d implicit convergence", [&] {
		KSPConvergedReason reason; OneDPetscCheck(KSPGetConvergedReason(solver, &reason), "KSPGetConvergedReason");
		if (reason <= 0) throw std::runtime_error("PETSc linearized_aq did not converge");
	});
	std::vector<double> values; OneDGetVectorAll(solution, values);
	OneDFlowState candidate;
	CollectiveLocalStage(communicator, "1d implicit state update", [&] {
		candidate = state;
		std::vector<double> pressure(values.begin(), values.begin()+graph.nodes);
		std::vector<double> edge_flow(values.begin()+graph.nodes, values.end());
		OneDUpdateStateFromImplicitSolution(network, graph, flow, pressure, edge_flow, candidate, dt);
	});
	objects.Close(communicator);
	static_assert(std::is_nothrow_move_assignable<OneDFlowState>::value, "implicit state publication must not throw");
	state = std::move(candidate);
}

struct OneDNonlinearContext {
	std::exception_ptr callback_error;
	const OneDNetwork* network = nullptr;
	const OneDImplicitGraph* graph = nullptr;
	const OneDFlowSystemDefinition* flow = nullptr;
	OneDFlowState* state = nullptr;
	double inlet_flow = 0.0;
	double dt = 0.0;
	std::vector<double> old_pressure;
	std::vector<double> old_flow;
	std::vector<double> compliance;
};

inline PetscErrorCode OneDNonlinearResidual(SNES, Vec input, Vec residual, void* raw) noexcept
{
	auto& context = *static_cast<OneDNonlinearContext*>(raw);
	try {
		const auto communicator = PetscObjectComm(reinterpret_cast<PetscObject>(input));
		std::vector<double> x;
		OneDGetVectorAll(input, x);
		CollectiveLocalStage(communicator, "1d nonlinear residual assembly", [&] {
			const auto& graph = *context.graph;
			std::vector<double> values(x.size(), 0.0);
			for (int node = 0; node < graph.nodes; ++node) {
				const auto* outlet = OneDOutletAtGraphNode(graph, context.state->outlets, node);
				if (outlet && outlet->kind == OneDOutletKind::Pressure) {
					values[static_cast<std::size_t>(node)] = x[static_cast<std::size_t>(node)]-outlet->pressure;
					continue;
				}
				double value = context.compliance[static_cast<std::size_t>(node)]
					*(x[static_cast<std::size_t>(node)]-context.old_pressure[static_cast<std::size_t>(node)])/context.dt;
				for (std::size_t edge = 0; edge < graph.edges.size(); ++edge) {
					if (graph.edges[edge].from == node) value += x[static_cast<std::size_t>(graph.nodes)+edge];
					if (graph.edges[edge].to == node) value -= x[static_cast<std::size_t>(graph.nodes)+edge];
				}
				if (node == graph.root) value -= context.inlet_flow;
				if (outlet) value += (x[static_cast<std::size_t>(node)]-OutletEffectivePressure(*outlet, context.dt))
					/OutletEffectiveResistance(*outlet, context.dt);
				values[static_cast<std::size_t>(node)] = value;
			}
			for (std::size_t edge_index = 0; edge_index < graph.edges.size(); ++edge_index) {
				const auto& edge = graph.edges[edge_index];
				const double p_from = x[static_cast<std::size_t>(edge.from)];
				const double p_to = x[static_cast<std::size_t>(edge.to)];
				const double pressure = 0.5*(p_from+p_to);
				const double area = OneDAreaFromPressure(pressure, edge.area0, edge.radius0, context.flow->wall);
				const double resistance = 8.0*OneDPi*context.flow->dynamic_viscosity*edge.length/(area*area);
				const double inertance = context.flow->density*edge.length/area;
				const double q = x[static_cast<std::size_t>(graph.nodes)+edge_index];
				double value = inertance*(q-context.old_flow.at(edge_index))/context.dt+resistance*q+p_to-p_from;
				if (edge.cell == 0) {
					const int incoming = OneDSegmentIntoNode(*context.network, edge.from);
					if (incoming >= 0 && context.network->nodes[static_cast<std::size_t>(edge.from)].children.size() > 1) {
						const auto& parent = context.network->segments[static_cast<std::size_t>(incoming)];
						const auto& child = context.network->segments[static_cast<std::size_t>(edge.segment)];
						const double k = OneDJunctionLossCoefficient(*context.flow, *context.network,
							edge.from, parent, child, 0.5);
						value += 0.5*context.flow->density*k*q*std::abs(q)/(area*area);
					}
				}
				values[static_cast<std::size_t>(graph.nodes)+edge_index] = value;
			}
			PetscInt first = 0, last = 0;
			OneDPetscCheck(VecGetOwnershipRange(residual, &first, &last), "VecGetOwnershipRange");
			for (PetscInt row = first; row < last; ++row) OneDPetscCheck(VecSetValue(residual, row, values[static_cast<std::size_t>(row)], INSERT_VALUES), "VecSetValue");
		});
		OneDCollectivePetscCheck(communicator, VecAssemblyBegin(residual), "VecAssemblyBegin");
		OneDCollectivePetscCheck(communicator, VecAssemblyEnd(residual), "VecAssemblyEnd");
		return 0;
	} catch (...) {
		context.callback_error = std::current_exception();
		return PETSC_ERR_USER;
	}
}

inline PetscErrorCode OneDNonlinearJacobian(SNES, Vec input, Mat jacobian, Mat, void* raw) noexcept
{
	auto& context = *static_cast<OneDNonlinearContext*>(raw);
	try {
		const auto communicator = PetscObjectComm(reinterpret_cast<PetscObject>(input));
		std::vector<double> x; OneDGetVectorAll(input, x);
		OneDCollectivePetscCheck(communicator, MatZeroEntries(jacobian), "MatZeroEntries");
		CollectiveLocalStage(communicator, "1d nonlinear jacobian assembly", [&] {
			const auto& graph = *context.graph;
			PetscInt first = 0, last = 0; OneDPetscCheck(MatGetOwnershipRange(jacobian, &first, &last), "MatGetOwnershipRange");
			for (PetscInt row = first; row < last; ++row) {
				if (row < graph.nodes) {
					const auto* outlet = OneDOutletAtGraphNode(graph, context.state->outlets, static_cast<int>(row));
					if (outlet && outlet->kind == OneDOutletKind::Pressure) {
						OneDPetscCheck(MatSetValue(jacobian, row, row, 1.0, INSERT_VALUES), "MatSetValue");
						continue;
					}
					double diagonal = context.compliance[static_cast<std::size_t>(row)]/context.dt;
					if (outlet) diagonal += 1.0/OutletEffectiveResistance(*outlet, context.dt);
					OneDPetscCheck(MatSetValue(jacobian, row, row, diagonal, INSERT_VALUES), "MatSetValue");
					for (std::size_t edge = 0; edge < graph.edges.size(); ++edge) {
						if (graph.edges[edge].from == row) OneDPetscCheck(MatSetValue(jacobian, row, graph.nodes+edge, 1.0, INSERT_VALUES), "MatSetValue");
						if (graph.edges[edge].to == row) OneDPetscCheck(MatSetValue(jacobian, row, graph.nodes+edge, -1.0, INSERT_VALUES), "MatSetValue");
					}
				} else {
					const std::size_t edge_index = static_cast<std::size_t>(row-graph.nodes);
					const auto& edge = graph.edges[edge_index];
					const double pressure = 0.5*(x[static_cast<std::size_t>(edge.from)]+x[static_cast<std::size_t>(edge.to)]);
					const double area = OneDAreaFromPressure(pressure, edge.area0, edge.radius0, context.flow->wall);
					const double da_dp = OneDWallAreaDerivative(pressure, edge.area0, edge.radius0, context.flow->wall);
					const double resistance = 8.0*OneDPi*context.flow->dynamic_viscosity*edge.length/(area*area);
					const double inertance = context.flow->density*edge.length/area;
					const double q = x[static_cast<std::size_t>(graph.nodes)+edge_index];
					const double dr_dp_endpoint = -resistance*da_dp/area;
					const double di_dp_endpoint = -0.5*inertance*da_dp/area;
					const double material = di_dp_endpoint*(q-context.old_flow.at(edge_index))/context.dt+dr_dp_endpoint*q;
					double dq = inertance/context.dt+resistance;
					double pressure_derivative = material;
					if (edge.cell == 0) {
						const int incoming = OneDSegmentIntoNode(*context.network, edge.from);
						if (incoming >= 0 && context.network->nodes[static_cast<std::size_t>(edge.from)].children.size() > 1) {
							const auto& parent = context.network->segments[static_cast<std::size_t>(incoming)];
							const auto& child = context.network->segments[static_cast<std::size_t>(edge.segment)];
							const double k = OneDJunctionLossCoefficient(*context.flow, *context.network, edge.from, parent, child, 0.5);
							const double coefficient = 0.5*context.flow->density*k;
							dq += 2.0*coefficient*std::abs(q)/(area*area);
							pressure_derivative -= coefficient*q*std::abs(q)*da_dp/(area*area*area);
						}
					}
					OneDPetscCheck(MatSetValue(jacobian, row, edge.from, pressure_derivative-1.0, INSERT_VALUES), "MatSetValue");
					OneDPetscCheck(MatSetValue(jacobian, row, edge.to, pressure_derivative+1.0, INSERT_VALUES), "MatSetValue");
					OneDPetscCheck(MatSetValue(jacobian, row, row, dq, INSERT_VALUES), "MatSetValue");
				}
			}
		});
		OneDCollectivePetscCheck(communicator, MatAssemblyBegin(jacobian, MAT_FINAL_ASSEMBLY), "MatAssemblyBegin");
		OneDCollectivePetscCheck(communicator, MatAssemblyEnd(jacobian, MAT_FINAL_ASSEMBLY), "MatAssemblyEnd");
		return 0;
	} catch (...) {
		context.callback_error = std::current_exception();
		return PETSC_ERR_USER;
	}
}

inline void SolveOneDNonlinearAQPetsc(const OneDNetwork& network,
	const OneDFlowSystemDefinition& flow, OneDFlowState& state,
	double inlet_flow, double dt, bool expand_cells, MPI_Comm communicator = PETSC_COMM_WORLD,
	OneDPetscSolverContext* solver_context = nullptr)
{
	std::optional<OneDPetscSolverContext> fallback;
	if (!solver_context) { fallback.emplace(communicator); solver_context = &*fallback; }
	OneDFlowState linear_guess;
	CollectiveLocalStage(communicator, "1d nonlinear initial state", [&] {
		linear_guess = state;
	});
	SolveOneDLinearizedAQPetsc(network, flow, linear_guess, inlet_flow, dt, expand_cells, communicator, solver_context);
	OneDImplicitGraph graph;
	CollectiveLocalStage(communicator, "1d nonlinear graph preparation", [&] {
		graph = BuildOneDImplicitGraph(network, expand_cells);
	});
	const PetscInt unknowns = graph.nodes+static_cast<int>(graph.edges.size());
	OneDRequireSameInteger(communicator, unknowns, "1d nonlinear layout");
	OneDPetscObjects objects;
	auto& solution = objects.solution;
	auto& residual = objects.rhs;
	auto& jacobian = objects.matrix;
	auto& solver = objects.snes;
	OneDCollectivePetscCheck(communicator, VecCreateMPI(communicator, PETSC_DECIDE, unknowns, &solution), "VecCreateMPI");
	OneDCollectivePetscCheck(communicator, VecDuplicate(solution, &residual), "VecDuplicate");
	OneDCollectivePetscCheck(communicator, MatCreateAIJ(communicator, PETSC_DECIDE, PETSC_DECIDE, unknowns, unknowns,
		10, nullptr, 10, nullptr, &jacobian), "MatCreateAIJ");
	OneDNonlinearContext context;
	std::vector<double> initial;
	CollectiveLocalStage(communicator, "1d nonlinear preparation", [&] {
		context.network = &network; context.graph = &graph; context.flow = &flow;
		context.state = &state; context.inlet_flow = inlet_flow; context.dt = dt;
		context.old_pressure.assign(static_cast<std::size_t>(graph.nodes), flow.wall.reference_pressure);
		for (int i = 0; i < graph.original_nodes && i < static_cast<int>(state.node_pressure.size()); ++i)
			context.old_pressure[static_cast<std::size_t>(i)] = state.node_pressure[static_cast<std::size_t>(i)];
		context.old_flow.assign(graph.edges.size(), 0.0);
		for (std::size_t i = 0; i < graph.edges.size(); ++i) {
			const auto& edge = graph.edges[i];
			if (expand_cells && state.flow.size() == static_cast<std::size_t>(network.cells)) {
				const auto& segment = network.segments[static_cast<std::size_t>(edge.segment)];
				context.old_flow[i] = state.flow[static_cast<std::size_t>(segment.cell_offset+edge.cell)];
			} else if (state.segment_flow.size() == network.segments.size())
				context.old_flow[i] = state.segment_flow[static_cast<std::size_t>(edge.segment)];
		}
		context.compliance = OneDNodeCompliance(graph, network, flow, context.old_pressure);
		initial.assign(static_cast<std::size_t>(unknowns), 0.0);
		std::copy(context.old_pressure.begin(), context.old_pressure.end(), initial.begin());
		for (int i = 0; i < graph.original_nodes && i < static_cast<int>(linear_guess.node_pressure.size()); ++i)
			initial[static_cast<std::size_t>(i)] = linear_guess.node_pressure[static_cast<std::size_t>(i)];
		for (std::size_t i = 0; i < graph.edges.size(); ++i) {
			const auto& edge = graph.edges[i];
			if (expand_cells && linear_guess.flow.size() == static_cast<std::size_t>(network.cells)) {
				const auto& segment = network.segments[static_cast<std::size_t>(edge.segment)];
				initial[static_cast<std::size_t>(graph.nodes)+i] =
					linear_guess.flow[static_cast<std::size_t>(segment.cell_offset+edge.cell)];
			} else if (linear_guess.segment_flow.size() == network.segments.size())
				initial[static_cast<std::size_t>(graph.nodes)+i] =
					linear_guess.segment_flow[static_cast<std::size_t>(edge.segment)];
		}
	});
	OneDSetInitialVector(solution, initial);
	OneDCollectivePetscCheck(communicator, SNESCreate(communicator, &solver), "SNESCreate");
	OneDCollectivePetscCheck(communicator, SNESSetFunction(solver, residual, OneDNonlinearResidual, &context), "SNESSetFunction");
	OneDCollectivePetscCheck(communicator, SNESSetJacobian(solver, jacobian, jacobian, OneDNonlinearJacobian, &context), "SNESSetJacobian");
	KSP nonlinear_ksp = nullptr;
	PC nonlinear_pc = nullptr;
	OneDCollectivePetscCheck(communicator, SNESGetKSP(solver, &nonlinear_ksp), "SNESGetKSP");
	OneDCollectivePetscCheck(communicator, KSPGetPC(nonlinear_ksp, &nonlinear_pc), "KSPGetPC");
	int mpi_size = 1;
	MPI_Comm_size(communicator, &mpi_size);
	if (mpi_size == 1) {
		OneDCollectivePetscCheck(communicator, KSPSetType(nonlinear_ksp, KSPPREONLY), "KSPSetType");
		OneDCollectivePetscCheck(communicator, PCSetType(nonlinear_pc, PCLU), "PCSetType");
	} else {
		OneDCollectivePetscCheck(communicator, KSPSetType(nonlinear_ksp, KSPPREONLY), "KSPSetType");
		OneDCollectivePetscCheck(communicator, PCSetType(nonlinear_pc, PCLU), "PCSetType");
		OneDCollectivePetscCheck(communicator, PCFactorSetMatSolverType(nonlinear_pc, MATSOLVERMUMPS), "PCFactorSetMatSolverType");
	}
	solver_context->options.Attach(solver);
	solver_context->options.Call("1d nonlinear options", [&] { return SNESSetFromOptions(solver); });
	RequireKspFactorBackend(nonlinear_ksp, jacobian, communicator);
	solver_context->options.Call("1d nonlinear solve", [&] {
		const auto code = SNESSolve(solver, nullptr, solution);
		if (context.callback_error) std::rethrow_exception(context.callback_error);
		return code;
	});
	solver_context->Record(communicator, nonlinear_ksp, solver);

	CollectiveLocalStage(communicator, "1d nonlinear convergence", [&] {
		SNESConvergedReason reason; OneDPetscCheck(SNESGetConvergedReason(solver, &reason), "SNESGetConvergedReason");
		if (reason <= 0) throw std::runtime_error("PETSc nonlinear 1d solve did not converge");
	});
	std::vector<double> values; OneDGetVectorAll(solution, values);
	OneDFlowState candidate;
	CollectiveLocalStage(communicator, "1d nonlinear state update", [&] {
		candidate = state;
		std::vector<double> pressure(values.begin(), values.begin()+graph.nodes);
		std::vector<double> edge_flow(values.begin()+graph.nodes, values.end());
		OneDUpdateStateFromImplicitSolution(network, graph, flow, pressure, edge_flow, candidate, dt);
	});
	objects.Close(communicator);
	state = std::move(candidate);
}

inline void AdvanceImplicitOneD(const OneDNetwork& network,
	const OneDFlowSystemDefinition& flow, OneDFlowState& state,
	double inlet_flow, double dt, MPI_Comm communicator = PETSC_COMM_WORLD,
	OneDPetscSolverContext* solver_context = nullptr)
{
	std::optional<OneDPetscSolverContext> fallback;
	if (!solver_context) { fallback.emplace(communicator); solver_context = &*fallback; }
	OneDRequireSameInteger(communicator, static_cast<int>(flow.formulation), "1d implicit formulation");
	OneDFlowState candidate;
	CollectiveLocalStage(communicator, "1d implicit advance preparation", [&] {
		candidate = state;
		candidate.inlet_flow = inlet_flow;
	});
	if (flow.formulation == OneDImplicitFormulation::SteadyR
		|| flow.formulation == OneDImplicitFormulation::PressureNetwork)
		SolveOneDPressureNetworkPetsc(network, flow, candidate, inlet_flow, dt, communicator, solver_context);
	else if (flow.formulation == OneDImplicitFormulation::LinearizedAQ)
		SolveOneDLinearizedAQPetsc(network, flow, candidate, inlet_flow, dt, false, communicator, solver_context);
	else if (flow.formulation == OneDImplicitFormulation::NonlinearAQ)
		SolveOneDNonlinearAQPetsc(network, flow, candidate, inlet_flow, dt, false, communicator, solver_context);
	else SolveOneDNonlinearAQPetsc(network, flow, candidate, inlet_flow, dt, true, communicator, solver_context);
	state = std::move(candidate);
}

} // namespace iga

#endif
