#ifndef IGA_ONE_D_IMPLICIT_HPP
#define IGA_ONE_D_IMPLICIT_HPP

#include "OneDFlow.hpp"

#include <petscksp.h>
#include <petscsnes.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

inline void OneDPetscCheck(PetscErrorCode error, const std::string& operation)
{
	if (error) throw std::runtime_error(operation + " failed with PETSc error " + std::to_string(error));
}

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
	std::vector<int> outlet_state_at_node;
	std::vector<int> node_edge_offsets;
	std::vector<int> node_edges;
	int nodes = 0;
};

inline OneDImplicitGraph BuildOneDImplicitGraph(const OneDNetwork& network,
	bool expand_cells)
{
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
	graph.outlet_state_at_node.assign(static_cast<std::size_t>(graph.nodes), -1);
	for (std::size_t i = 0; i < graph.outlet_graph_nodes.size(); ++i)
		graph.outlet_state_at_node[static_cast<std::size_t>(graph.outlet_graph_nodes[i])]
			= graph.outlet_state_indices[i];
	graph.node_edge_offsets.assign(static_cast<std::size_t>(graph.nodes)+1, 0);
	for (const auto& edge : graph.edges) {
		++graph.node_edge_offsets[static_cast<std::size_t>(edge.from)+1];
		++graph.node_edge_offsets[static_cast<std::size_t>(edge.to)+1];
	}
	for (int node = 0; node < graph.nodes; ++node)
		graph.node_edge_offsets[static_cast<std::size_t>(node)+1]
			+= graph.node_edge_offsets[static_cast<std::size_t>(node)];
	graph.node_edges.resize(2*graph.edges.size());
	auto next = graph.node_edge_offsets;
	for (std::size_t edge = 0; edge < graph.edges.size(); ++edge) {
		const auto& item = graph.edges[edge];
		graph.node_edges[static_cast<std::size_t>(next[static_cast<std::size_t>(item.from)]++)]
			= static_cast<int>(edge);
		graph.node_edges[static_cast<std::size_t>(next[static_cast<std::size_t>(item.to)]++)]
			= static_cast<int>(edge);
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
	if (node < 0 || node >= static_cast<int>(graph.outlet_state_at_node.size())) return nullptr;
	const int index = graph.outlet_state_at_node[static_cast<std::size_t>(node)];
	return index < 0 ? nullptr : &outlets[static_cast<std::size_t>(index)];
}

inline void OneDGetVectorAll(Vec vector, std::vector<double>& values)
{
	Vec all = nullptr;
	VecScatter scatter = nullptr;
	OneDPetscCheck(VecScatterCreateToAll(vector, &scatter, &all), "VecScatterCreateToAll");
	OneDPetscCheck(VecScatterBegin(scatter, vector, all, INSERT_VALUES, SCATTER_FORWARD), "VecScatterBegin");
	OneDPetscCheck(VecScatterEnd(scatter, vector, all, INSERT_VALUES, SCATTER_FORWARD), "VecScatterEnd");
	PetscInt size = 0;
	OneDPetscCheck(VecGetSize(all, &size), "VecGetSize");
	const PetscScalar* array = nullptr;
	OneDPetscCheck(VecGetArrayRead(all, &array), "VecGetArrayRead");
	values.resize(static_cast<std::size_t>(size));
	for (PetscInt i = 0; i < size; ++i) values[static_cast<std::size_t>(i)] = PetscRealPart(array[i]);
	OneDPetscCheck(VecRestoreArrayRead(all, &array), "VecRestoreArrayRead");
	VecScatterDestroy(&scatter);
	VecDestroy(&all);
}

inline void OneDSetInitialVector(Vec vector, const std::vector<double>& values)
{
	PetscInt first = 0, last = 0;
	OneDPetscCheck(VecGetOwnershipRange(vector, &first, &last), "VecGetOwnershipRange");
	for (PetscInt i = first; i < last; ++i)
		OneDPetscCheck(VecSetValue(vector, i, values[static_cast<std::size_t>(i)], INSERT_VALUES), "VecSetValue");
	OneDPetscCheck(VecAssemblyBegin(vector), "VecAssemblyBegin");
	OneDPetscCheck(VecAssemblyEnd(vector), "VecAssemblyEnd");
}

inline void OneDCreateMixedSparseSystem(const OneDImplicitGraph& graph,
	Mat* matrix, Vec* vector)
{
	PetscInt global = graph.nodes+static_cast<PetscInt>(graph.edges.size());
	PetscInt local = PETSC_DECIDE;
	OneDPetscCheck(PetscSplitOwnership(PETSC_COMM_WORLD, &local, &global),
		"PetscSplitOwnership mixed system");
	PetscInt first = 0;
	OneDPetscCheck(MPI_Exscan(&local, &first, 1, MPIU_INT, MPI_SUM,
		PETSC_COMM_WORLD), "MPI_Exscan mixed system");
	int rank = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	if (rank == 0) first = 0;
	const PetscInt last = first+local;
	std::vector<PetscInt> diagonal(static_cast<std::size_t>(local), 0);
	std::vector<PetscInt> off_diagonal(static_cast<std::size_t>(local), 0);
	auto count = [&](PetscInt row, PetscInt column) {
		auto& target = column >= first && column < last ? diagonal : off_diagonal;
		++target[static_cast<std::size_t>(row-first)];
	};
	for (PetscInt row = first; row < last; ++row) {
		count(row, row);
		if (row < graph.nodes) {
			for (int position = graph.node_edge_offsets[static_cast<std::size_t>(row)];
				position < graph.node_edge_offsets[static_cast<std::size_t>(row)+1]; ++position)
				count(row, graph.nodes+graph.node_edges[static_cast<std::size_t>(position)]);
		} else {
			const auto& edge = graph.edges[static_cast<std::size_t>(row-graph.nodes)];
			count(row, edge.from);
			count(row, edge.to);
		}
	}
	OneDPetscCheck(MatCreateAIJ(PETSC_COMM_WORLD, local, local, global, global,
		0, diagonal.data(), 0, off_diagonal.data(), matrix),
		"MatCreateAIJ mixed system");
	OneDPetscCheck(MatSetOption(*matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE),
		"MatSetOption mixed system");
	OneDPetscCheck(VecCreateMPI(PETSC_COMM_WORLD, local, global, vector),
		"VecCreateMPI mixed system");
}

struct OneDPressureNetworkWorkspace {
	OneDImplicitGraph graph;
	Mat matrix = nullptr;
	Vec rhs = nullptr;
	Vec solution = nullptr;
	KSP solver = nullptr;
	PetscInt global_nodes = 0;
	PetscInt local_nodes = 0;
	PetscInt first = 0;
	PetscInt last = 0;

	OneDPressureNetworkWorkspace() = default;
	OneDPressureNetworkWorkspace(const OneDPressureNetworkWorkspace&) = delete;
	OneDPressureNetworkWorkspace& operator=(const OneDPressureNetworkWorkspace&) = delete;

	~OneDPressureNetworkWorkspace()
	{
		KSPDestroy(&solver);
		VecDestroy(&solution);
		VecDestroy(&rhs);
		MatDestroy(&matrix);
	}

	void Initialize(const OneDNetwork& network)
	{
		if (matrix) return;
		graph = BuildOneDImplicitGraph(network, false);
		global_nodes = graph.nodes;
		local_nodes = PETSC_DECIDE;
		OneDPetscCheck(PetscSplitOwnership(PETSC_COMM_WORLD,
			&local_nodes, &global_nodes), "PetscSplitOwnership pressure network");
		OneDPetscCheck(MPI_Exscan(&local_nodes, &first, 1, MPIU_INT, MPI_SUM,
			PETSC_COMM_WORLD), "MPI_Exscan pressure network");
		int rank = 0;
		MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
		if (rank == 0) first = 0;
		last = first+local_nodes;
		std::vector<PetscInt> diagonal_nonzeros(
			static_cast<std::size_t>(local_nodes), 1);
		std::vector<PetscInt> off_diagonal_nonzeros(
			static_cast<std::size_t>(local_nodes), 0);
		for (PetscInt row = first; row < last; ++row)
			for (int position = graph.node_edge_offsets[static_cast<std::size_t>(row)];
				position < graph.node_edge_offsets[static_cast<std::size_t>(row)+1]; ++position) {
				const auto& edge = graph.edges[static_cast<std::size_t>(
					graph.node_edges[static_cast<std::size_t>(position)])];
				const int other = edge.from == row ? edge.to : edge.from;
				if (other >= first && other < last)
					++diagonal_nonzeros[static_cast<std::size_t>(row-first)];
				else ++off_diagonal_nonzeros[static_cast<std::size_t>(row-first)];
			}
		OneDPetscCheck(MatCreateAIJ(PETSC_COMM_WORLD, local_nodes, local_nodes,
			global_nodes, global_nodes, 0, diagonal_nonzeros.data(), 0,
			off_diagonal_nonzeros.data(), &matrix), "MatCreateAIJ pressure network");
		OneDPetscCheck(MatSetOption(matrix, MAT_NEW_NONZERO_ALLOCATION_ERR,
			PETSC_TRUE), "MatSetOption pressure network");
		OneDPetscCheck(VecCreateMPI(PETSC_COMM_WORLD, local_nodes,
			global_nodes, &rhs), "VecCreateMPI pressure network");
		OneDPetscCheck(VecDuplicate(rhs, &solution), "VecDuplicate pressure network");
		OneDPetscCheck(KSPCreate(PETSC_COMM_WORLD, &solver), "KSPCreate pressure network");
		OneDPetscCheck(KSPSetType(solver, KSPCG), "KSPSetType pressure network");
		PC preconditioner = nullptr;
		OneDPetscCheck(KSPGetPC(solver, &preconditioner), "KSPGetPC pressure network");
		OneDPetscCheck(PCSetType(preconditioner, PCJACOBI), "PCSetType pressure network");
		OneDPetscCheck(KSPSetTolerances(solver, 1.0e-12, PETSC_DEFAULT,
			PETSC_DEFAULT, PETSC_DEFAULT), "KSPSetTolerances pressure network");
		OneDPetscCheck(KSPSetFromOptions(solver), "KSPSetFromOptions pressure network");
	}
};

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
		double outlet_flow = state.segment_flow[static_cast<std::size_t>(incoming)];
		if (outlet.kind != OneDOutletKind::Pressure)
			outlet_flow = (pressure[static_cast<std::size_t>(outlet.node)]
				-OutletEffectivePressure(outlet, dt))/OutletEffectiveResistance(outlet, dt);
		AdvanceOutletState(outlet, outlet_flow, dt);
	}
}

inline void SolveOneDPressureNetworkPetsc(const OneDNetwork& network,
	const OneDFlowSystemDefinition& flow, OneDFlowState& state,
	double inlet_flow, double dt, OneDPressureNetworkWorkspace* persistent = nullptr)
{
	OneDPressureNetworkWorkspace local;
	auto& workspace = persistent ? *persistent : local;
	workspace.Initialize(network);
	const auto& graph = workspace.graph;
	const PetscInt n = workspace.global_nodes;
	const PetscInt local_n = workspace.local_nodes;
	const PetscInt first = workspace.first;
	const PetscInt last = workspace.last;
	Mat matrix = workspace.matrix;
	Vec rhs = workspace.rhs;
	Vec solution = workspace.solution;
	KSP solver = workspace.solver;
	OneDPetscCheck(MatZeroEntries(matrix), "MatZeroEntries pressure network");
	OneDPetscCheck(VecZeroEntries(rhs), "VecZeroEntries pressure network rhs");
	std::vector<double> previous = state.node_pressure;
	if (previous.size() != static_cast<std::size_t>(n)) previous.assign(static_cast<std::size_t>(n), flow.wall.reference_pressure);
	const auto compliance = OneDNodeCompliance(graph, network, flow, previous);
	std::vector<double> diagonal(static_cast<std::size_t>(local_n), 0.0);
	std::vector<double> right_hand_side(static_cast<std::size_t>(local_n), 0.0);
	for (PetscInt row = first; row < last; ++row) {
		const auto* outlet = OneDOutletAtGraphNode(graph, state.outlets, static_cast<int>(row));
		if (outlet && outlet->kind == OneDOutletKind::Pressure) {
			diagonal[static_cast<std::size_t>(row-first)] = 1.0;
			right_hand_side[static_cast<std::size_t>(row-first)] = outlet->pressure;
			continue;
		}
		auto& row_diagonal = diagonal[static_cast<std::size_t>(row-first)];
		auto& value = right_hand_side[static_cast<std::size_t>(row-first)];
		row_diagonal = compliance[static_cast<std::size_t>(row)]/dt;
		value = row_diagonal*previous[static_cast<std::size_t>(row)];
		if (row == graph.root) value += inlet_flow;
		if (outlet) {
			const double resistance = OutletEffectiveResistance(*outlet, dt);
			const double pressure = OutletEffectivePressure(*outlet, dt);
			row_diagonal += 1.0/resistance;
			value += pressure/resistance;
		}
	}
	for (const auto& edge : graph.edges) {
		const auto& segment = network.segments[static_cast<std::size_t>(edge.segment)];
		const double conductance = 1.0/OneDLumpedSegmentResistance(network, flow, segment);
		for (const auto& row_other : {std::pair<int, int>{edge.from, edge.to},
			std::pair<int, int>{edge.to, edge.from}}) {
			const int row = row_other.first;
			if (row < first || row >= last) continue;
			const auto* outlet = OneDOutletAtGraphNode(graph, state.outlets, row);
			if (outlet && outlet->kind == OneDOutletKind::Pressure) continue;
			diagonal[static_cast<std::size_t>(row-first)] += conductance;
			const auto* other_outlet = OneDOutletAtGraphNode(
				graph, state.outlets, row_other.second);
			if (other_outlet && other_outlet->kind == OneDOutletKind::Pressure) {
				right_hand_side[static_cast<std::size_t>(row-first)]
					+= conductance*other_outlet->pressure;
				continue;
			}
			OneDPetscCheck(MatSetValue(matrix, row, row_other.second,
				-conductance, INSERT_VALUES), "MatSetValue off-diagonal");
		}
	}
	for (PetscInt row = first; row < last; ++row) {
		OneDPetscCheck(MatSetValue(matrix, row, row,
			diagonal[static_cast<std::size_t>(row-first)], INSERT_VALUES),
			"MatSetValue diagonal");
		OneDPetscCheck(VecSetValue(rhs, row,
			right_hand_side[static_cast<std::size_t>(row-first)], INSERT_VALUES),
			"VecSetValue rhs");
	}
	OneDPetscCheck(MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY), "MatAssemblyBegin");
	OneDPetscCheck(MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY), "MatAssemblyEnd");
	OneDPetscCheck(VecAssemblyBegin(rhs), "VecAssemblyBegin");
	OneDPetscCheck(VecAssemblyEnd(rhs), "VecAssemblyEnd");
	OneDSetInitialVector(solution, previous);
	OneDPetscCheck(KSPSetOperators(solver, matrix, matrix), "KSPSetOperators");
	KSPType solver_type = nullptr;
	OneDPetscCheck(KSPGetType(solver, &solver_type), "KSPGetType pressure network");
	const PetscBool use_initial_guess = solver_type
		&& std::string(solver_type) == KSPPREONLY ? PETSC_FALSE : PETSC_TRUE;
	OneDPetscCheck(KSPSetInitialGuessNonzero(solver, use_initial_guess),
		"KSPSetInitialGuessNonzero");
	OneDPetscCheck(KSPSolve(solver, rhs, solution), "KSPSolve");
	KSPConvergedReason reason;
	KSPGetConvergedReason(solver, &reason);
	if (reason <= 0) throw std::runtime_error("PETSc 0D R/RC network did not converge");
	std::vector<double> pressure;
	OneDGetVectorAll(solution, pressure);
	std::vector<double> edge_flow(graph.edges.size());
	for (std::size_t i = 0; i < graph.edges.size(); ++i) {
		const auto& edge = graph.edges[i];
		const auto& segment = network.segments[static_cast<std::size_t>(edge.segment)];
		const double resistance = OneDLumpedSegmentResistance(network, flow, segment);
		edge_flow[i] = (pressure[static_cast<std::size_t>(edge.from)]-pressure[static_cast<std::size_t>(edge.to)])/resistance;
	}
	OneDUpdateStateFromImplicitSolution(network, graph, flow, pressure, edge_flow, state, dt);
	double storage_rate = 0.0;
	for (int node = 0; node < graph.nodes; ++node) {
		const auto* outlet = OneDOutletAtGraphNode(graph, state.outlets, node);
		if (outlet && outlet->kind == OneDOutletKind::Pressure) continue;
		storage_rate += compliance[static_cast<std::size_t>(node)]
			*(pressure[static_cast<std::size_t>(node)]-previous[static_cast<std::size_t>(node)])/dt;
	}
	double outlet_flow = 0.0;
	for (const auto& outlet : state.outlets) outlet_flow += outlet.flow;
	const double residual = inlet_flow-outlet_flow-storage_rate;
	state.storage_rate = storage_rate;
	state.relative_continuity_residual = std::abs(residual)
		/std::max({std::abs(inlet_flow), std::abs(outlet_flow),
			std::abs(storage_rate), 1.0e-30});
	state.has_conservation_diagnostic = true;
}

inline void SolveOneDLinearizedAQPetsc(const OneDNetwork& network,
	const OneDFlowSystemDefinition& flow, OneDFlowState& state,
	double inlet_flow, double dt, bool expand_cells)
{
	const auto graph = BuildOneDImplicitGraph(network, expand_cells);
	Mat matrix = nullptr;
	Vec rhs = nullptr, solution = nullptr;
	KSP solver = nullptr;
	OneDCreateMixedSparseSystem(graph, &matrix, &rhs);
	VecDuplicate(rhs, &solution);
	std::vector<double> previous_pressure(static_cast<std::size_t>(graph.nodes), flow.wall.reference_pressure);
	for (int i = 0; i < graph.original_nodes && i < static_cast<int>(state.node_pressure.size()); ++i)
		previous_pressure[static_cast<std::size_t>(i)] = state.node_pressure[static_cast<std::size_t>(i)];
	const auto compliance = OneDNodeCompliance(graph, network, flow, previous_pressure);
	PetscInt first = 0, last = 0;
	MatGetOwnershipRange(matrix, &first, &last);
	for (PetscInt row = first; row < last; ++row) {
		if (row < graph.nodes) {
			const auto* outlet = OneDOutletAtGraphNode(graph, state.outlets, static_cast<int>(row));
			if (outlet && outlet->kind == OneDOutletKind::Pressure) {
				MatSetValue(matrix, row, row, 1.0, ADD_VALUES);
				VecSetValue(rhs, row, outlet->pressure, INSERT_VALUES);
				continue;
			}
			double diagonal = compliance[static_cast<std::size_t>(row)]/dt;
			double value = diagonal*previous_pressure[static_cast<std::size_t>(row)];
			for (int position = graph.node_edge_offsets[static_cast<std::size_t>(row)];
				position < graph.node_edge_offsets[static_cast<std::size_t>(row)+1]; ++position) {
				const int edge_index = graph.node_edges[static_cast<std::size_t>(position)];
				const auto& edge = graph.edges[static_cast<std::size_t>(edge_index)];
				MatSetValue(matrix, row, graph.nodes+edge_index,
					edge.from == row ? 1.0 : -1.0, ADD_VALUES);
			}
			if (row == graph.root) value += inlet_flow;
			if (outlet) {
				const double resistance = OutletEffectiveResistance(*outlet, dt);
				const double pressure = OutletEffectivePressure(*outlet, dt);
				diagonal += 1.0/resistance;
				value += pressure/resistance;
			}
			MatSetValue(matrix, row, row, diagonal, ADD_VALUES);
			VecSetValue(rhs, row, value, INSERT_VALUES);
		} else {
			const int edge_index = static_cast<int>(row)-graph.nodes;
			const auto& edge = graph.edges[static_cast<std::size_t>(edge_index)];
			const auto& segment = network.segments[static_cast<std::size_t>(edge.segment)];
			const double resistance = expand_cells
				? 8.0*OneDPi*flow.dynamic_viscosity*edge.length/(edge.area0*edge.area0)
				: OneDLumpedSegmentResistance(network, flow, segment);
			const double inertance = flow.density*edge.length/edge.area0;
			double previous_flow = 0.0;
			if (expand_cells) {
				if (state.flow.size() == static_cast<std::size_t>(network.cells))
					previous_flow = state.flow[static_cast<std::size_t>(segment.cell_offset+edge.cell)];
			} else if (state.segment_flow.size() == network.segments.size())
				previous_flow = state.segment_flow[static_cast<std::size_t>(edge.segment)];
			MatSetValue(matrix, row, edge.from, -1.0, ADD_VALUES);
			MatSetValue(matrix, row, edge.to, 1.0, ADD_VALUES);
			MatSetValue(matrix, row, row, resistance+inertance/dt, ADD_VALUES);
			VecSetValue(rhs, row, inertance*previous_flow/dt, INSERT_VALUES);
		}
	}
	MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(rhs); VecAssemblyEnd(rhs);
	KSPCreate(PETSC_COMM_WORLD, &solver); KSPSetOperators(solver, matrix, matrix); KSPSetFromOptions(solver);
	KSPSolve(solver, rhs, solution);
	KSPConvergedReason reason; KSPGetConvergedReason(solver, &reason);
	if (reason <= 0) throw std::runtime_error(expand_cells
		? "PETSc 1D linear predictor did not converge"
		: "PETSc 0D transient_rlc did not converge");
	std::vector<double> values; OneDGetVectorAll(solution, values);
	std::vector<double> pressure(values.begin(), values.begin()+graph.nodes);
	std::vector<double> edge_flow(values.begin()+graph.nodes, values.end());
	OneDUpdateStateFromImplicitSolution(network, graph, flow, pressure, edge_flow, state, dt);
	KSPDestroy(&solver); VecDestroy(&solution); VecDestroy(&rhs); MatDestroy(&matrix);
}

struct OneDNonlinearContext {
	const OneDNetwork* network = nullptr;
	const OneDImplicitGraph* graph = nullptr;
	const OneDFlowSystemDefinition* flow = nullptr;
	OneDFlowState* state = nullptr;
	double inlet_flow = 0.0;
	double dt = 0.0;
	bool expand_cells = false;
	std::vector<double> old_pressure;
	std::vector<double> old_flow;
	std::vector<double> compliance;
};

inline PetscErrorCode OneDNonlinearResidual(SNES, Vec input, Vec residual, void* raw)
{
	auto& context = *static_cast<OneDNonlinearContext*>(raw);
	std::vector<double> x;
	OneDGetVectorAll(input, x);
	const auto& graph = *context.graph;
	PetscInt first = 0, last = 0;
	VecGetOwnershipRange(residual, &first, &last);
	for (PetscInt row = first; row < last; ++row) {
		double value = 0.0;
		if (row < graph.nodes) {
			const int node = static_cast<int>(row);
			const auto* outlet = OneDOutletAtGraphNode(graph, context.state->outlets, node);
			if (outlet && outlet->kind == OneDOutletKind::Pressure)
				value = x[static_cast<std::size_t>(node)]-outlet->pressure;
			else {
				value = context.compliance[static_cast<std::size_t>(node)]
					*(x[static_cast<std::size_t>(node)]
						-context.old_pressure[static_cast<std::size_t>(node)])/context.dt;
				for (int position = graph.node_edge_offsets[static_cast<std::size_t>(node)];
					position < graph.node_edge_offsets[static_cast<std::size_t>(node)+1]; ++position) {
					const int edge_index = graph.node_edges[static_cast<std::size_t>(position)];
					const auto& edge = graph.edges[static_cast<std::size_t>(edge_index)];
					value += (edge.from == node ? 1.0 : -1.0)
						*x[static_cast<std::size_t>(graph.nodes+edge_index)];
				}
				if (node == graph.root) value -= context.inlet_flow;
				if (outlet) value += (x[static_cast<std::size_t>(node)]
					-OutletEffectivePressure(*outlet, context.dt))
					/OutletEffectiveResistance(*outlet, context.dt);
			}
		} else {
			const std::size_t edge_index = static_cast<std::size_t>(row-graph.nodes);
			const auto& edge = graph.edges[edge_index];
			const auto& segment = context.network->segments[static_cast<std::size_t>(edge.segment)];
			const double p_from = x[static_cast<std::size_t>(edge.from)];
			const double p_to = x[static_cast<std::size_t>(edge.to)];
			const double pressure = 0.5*(p_from+p_to);
			const double area = OneDAreaFromPressure(pressure,
				edge.area0, edge.radius0, context.flow->wall);
			const int child_id = context.network->nodes[static_cast<std::size_t>(segment.child)].id;
			const auto configured = context.flow->lumped.segment_resistance.find(child_id);
			const bool fixed_resistance = !context.expand_cells
				&& configured != context.flow->lumped.segment_resistance.end();
			double resistance = fixed_resistance ? configured->second
				: 8.0*OneDPi*context.flow->dynamic_viscosity*edge.length/(area*area);
			if (!context.expand_cells && !fixed_resistance)
				resistance *= context.flow->lumped.resistance_scale;
			const double inertance = context.flow->density*edge.length/area;
			const double q = x[static_cast<std::size_t>(graph.nodes)+edge_index];
			value = inertance*(q-context.old_flow[edge_index])/context.dt
				+resistance*q+p_to-p_from;
			if (edge.cell == 0) {
				const int incoming = OneDSegmentIntoNode(*context.network, edge.from);
				if (incoming >= 0 && context.network->nodes[
					static_cast<std::size_t>(edge.from)].children.size() > 1) {
					const auto& parent = context.network->segments[static_cast<std::size_t>(incoming)];
					const auto& child = context.network->segments[static_cast<std::size_t>(edge.segment)];
					const double k = OneDJunctionLossCoefficient(*context.flow,
						*context.network, edge.from, parent, child, 0.5);
					value += 0.5*context.flow->density*k*q*std::abs(q)/(area*area);
				}
			}
		}
		VecSetValue(residual, row, value, INSERT_VALUES);
	}
	VecAssemblyBegin(residual); VecAssemblyEnd(residual);
	return 0;
}

inline PetscErrorCode OneDNonlinearJacobian(SNES, Vec input, Mat jacobian, Mat, void* raw)
{
	auto& context = *static_cast<OneDNonlinearContext*>(raw);
	std::vector<double> x; OneDGetVectorAll(input, x);
	const auto& graph = *context.graph;
	MatZeroEntries(jacobian);
	PetscInt first = 0, last = 0; MatGetOwnershipRange(jacobian, &first, &last);
	for (PetscInt row = first; row < last; ++row) {
		if (row < graph.nodes) {
			const auto* outlet = OneDOutletAtGraphNode(graph, context.state->outlets, static_cast<int>(row));
			if (outlet && outlet->kind == OneDOutletKind::Pressure) {
				MatSetValue(jacobian, row, row, 1.0, INSERT_VALUES);
				continue;
			}
			double diagonal = context.compliance[static_cast<std::size_t>(row)]/context.dt;
			if (outlet) diagonal += 1.0/OutletEffectiveResistance(*outlet, context.dt);
			MatSetValue(jacobian, row, row, diagonal, INSERT_VALUES);
			for (int position = graph.node_edge_offsets[static_cast<std::size_t>(row)];
				position < graph.node_edge_offsets[static_cast<std::size_t>(row)+1]; ++position) {
				const int edge_index = graph.node_edges[static_cast<std::size_t>(position)];
				const auto& edge = graph.edges[static_cast<std::size_t>(edge_index)];
				MatSetValue(jacobian, row, graph.nodes+edge_index,
					edge.from == row ? 1.0 : -1.0, INSERT_VALUES);
			}
		} else {
			const std::size_t edge_index = static_cast<std::size_t>(row-graph.nodes);
			const auto& edge = graph.edges[edge_index];
			const auto& segment = context.network->segments[static_cast<std::size_t>(edge.segment)];
			const double pressure = 0.5*(x[static_cast<std::size_t>(edge.from)]+x[static_cast<std::size_t>(edge.to)]);
			const double area = OneDAreaFromPressure(pressure, edge.area0, edge.radius0, context.flow->wall);
			const double da_dp = OneDWallAreaDerivative(pressure, edge.area0, edge.radius0, context.flow->wall);
			const int child_id = context.network->nodes[static_cast<std::size_t>(segment.child)].id;
			const auto configured = context.flow->lumped.segment_resistance.find(child_id);
			const bool fixed_resistance = !context.expand_cells
				&& configured != context.flow->lumped.segment_resistance.end();
			double resistance = fixed_resistance ? configured->second
				: 8.0*OneDPi*context.flow->dynamic_viscosity*edge.length/(area*area);
			if (!context.expand_cells && !fixed_resistance)
				resistance *= context.flow->lumped.resistance_scale;
			const double inertance = context.flow->density*edge.length/area;
			const double q = x[static_cast<std::size_t>(graph.nodes)+edge_index];
			const double dr_dp_endpoint = fixed_resistance ? 0.0
				: -resistance*da_dp/area;
			const double di_dp_endpoint = -0.5*inertance*da_dp/area;
			const double material = di_dp_endpoint*(q-context.old_flow[edge_index])/context.dt+dr_dp_endpoint*q;
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
			MatSetValue(jacobian, row, edge.from, pressure_derivative-1.0, INSERT_VALUES);
			MatSetValue(jacobian, row, edge.to, pressure_derivative+1.0, INSERT_VALUES);
			MatSetValue(jacobian, row, row, dq, INSERT_VALUES);
		}
	}
	MatAssemblyBegin(jacobian, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(jacobian, MAT_FINAL_ASSEMBLY);
	return 0;
}

inline void SolveOneDNonlinearAQPetsc(const OneDNetwork& network,
	const OneDFlowSystemDefinition& flow, OneDFlowState& state,
	double inlet_flow, double dt, bool expand_cells)
{
	OneDFlowState linear_guess = state;
	SolveOneDLinearizedAQPetsc(network, flow, linear_guess, inlet_flow, dt, expand_cells);
	const auto graph = BuildOneDImplicitGraph(network, expand_cells);
	const PetscInt unknowns = graph.nodes+static_cast<int>(graph.edges.size());
	Vec solution = nullptr, residual = nullptr;
	Mat jacobian = nullptr;
	SNES solver = nullptr;
	OneDCreateMixedSparseSystem(graph, &jacobian, &solution);
	VecDuplicate(solution, &residual);
	OneDNonlinearContext context;
	context.network = &network; context.graph = &graph; context.flow = &flow;
	context.state = &state; context.inlet_flow = inlet_flow; context.dt = dt;
	context.expand_cells = expand_cells;
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
	std::vector<double> initial(static_cast<std::size_t>(unknowns), 0.0);
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
	OneDSetInitialVector(solution, initial);
	SNESCreate(PETSC_COMM_WORLD, &solver);
	SNESSetFunction(solver, residual, OneDNonlinearResidual, &context);
	SNESSetJacobian(solver, jacobian, jacobian, OneDNonlinearJacobian, &context);
	KSP nonlinear_ksp = nullptr;
	PC nonlinear_pc = nullptr;
	SNESGetKSP(solver, &nonlinear_ksp);
	KSPGetPC(nonlinear_ksp, &nonlinear_pc);
	int mpi_size = 1;
	MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size);
	if (mpi_size == 1) {
		KSPSetType(nonlinear_ksp, KSPPREONLY);
		PCSetType(nonlinear_pc, PCLU);
	} else {
		KSPSetType(nonlinear_ksp, KSPPREONLY);
		PCSetType(nonlinear_pc, PCLU);
		PCFactorSetMatSolverType(nonlinear_pc, MATSOLVERMUMPS);
	}
	SNESSetFromOptions(solver);
	SNESSolve(solver, nullptr, solution);
	SNESConvergedReason reason; SNESGetConvergedReason(solver, &reason);
	if (reason <= 0) throw std::runtime_error(expand_cells
		? "PETSc implicit_1d_pde did not converge"
		: "PETSc 0D nonlinear_rlc did not converge");
	std::vector<double> values; OneDGetVectorAll(solution, values);
	std::vector<double> pressure(values.begin(), values.begin()+graph.nodes);
	std::vector<double> edge_flow(values.begin()+graph.nodes, values.end());
	OneDUpdateStateFromImplicitSolution(network, graph, flow, pressure, edge_flow, state, dt);
	SNESDestroy(&solver); MatDestroy(&jacobian); VecDestroy(&residual); VecDestroy(&solution);
}

inline void AdvanceImplicitOneD(const OneDNetwork& network,
	const OneDFlowSystemDefinition& flow, OneDFlowState& state,
	double inlet_flow, double dt, OneDPressureNetworkWorkspace* pressure_workspace = nullptr)
{
	state.inlet_flow = inlet_flow;
	if (flow.formulation == OneDImplicitFormulation::SteadyR
		|| flow.formulation == OneDImplicitFormulation::TransientRc)
		SolveOneDPressureNetworkPetsc(network, flow, state, inlet_flow, dt,
			pressure_workspace);
	else if (flow.formulation == OneDImplicitFormulation::TransientRlc)
		SolveOneDLinearizedAQPetsc(network, flow, state, inlet_flow, dt, false);
	else if (flow.formulation == OneDImplicitFormulation::NonlinearRlc)
		SolveOneDNonlinearAQPetsc(network, flow, state, inlet_flow, dt, false);
	else SolveOneDNonlinearAQPetsc(network, flow, state, inlet_flow, dt, true);
}

} // namespace iga

#endif
