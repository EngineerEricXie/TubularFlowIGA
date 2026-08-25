#ifndef IGA_ONE_D_IMPLICIT_HPP
#define IGA_ONE_D_IMPLICIT_HPP

#include "OneDFlow.hpp"

#include <petscksp.h>
#include <petscsnes.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
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
	return graph;
}

inline std::vector<double> OneDNodeCompliance(const OneDImplicitGraph& graph,
	const OneDFlowSystemDefinition& flow, const std::vector<double>& pressure)
{
	std::vector<double> result(static_cast<std::size_t>(graph.nodes), 0.0);
	for (const auto& edge : graph.edges) {
		const double mean_pressure = 0.5*(pressure[static_cast<std::size_t>(edge.from)]
			+pressure[static_cast<std::size_t>(edge.to)]);
		const double derivative = OneDWallAreaDerivative(mean_pressure,
			edge.area0, edge.radius0, flow.wall)*edge.length;
		result[static_cast<std::size_t>(edge.from)] += 0.5*derivative;
		result[static_cast<std::size_t>(edge.to)] += 0.5*derivative;
	}
	return result;
}

inline const OneDOutletState* OneDOutletAtGraphNode(const OneDImplicitGraph& graph,
	const std::vector<OneDOutletState>& outlets, int node)
{
	for (std::size_t i = 0; i < graph.outlet_graph_nodes.size(); ++i)
		if (graph.outlet_graph_nodes[i] == node)
			return &outlets[static_cast<std::size_t>(graph.outlet_state_indices[i])];
	return nullptr;
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
		AdvanceOutletState(outlet, state.segment_flow[static_cast<std::size_t>(incoming)], dt);
	}
}

inline void SolveOneDPressureNetworkPetsc(const OneDNetwork& network,
	const OneDFlowSystemDefinition& flow, OneDFlowState& state,
	double inlet_flow, double dt)
{
	const auto graph = BuildOneDImplicitGraph(network, false);
	const PetscInt n = graph.nodes;
	Mat matrix = nullptr;
	Vec rhs = nullptr, solution = nullptr;
	KSP solver = nullptr;
	OneDPetscCheck(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, n, n,
		8, nullptr, 8, nullptr, &matrix), "MatCreateAIJ");
	OneDPetscCheck(VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, n, &rhs), "VecCreateMPI");
	OneDPetscCheck(VecDuplicate(rhs, &solution), "VecDuplicate");
	PetscInt first = 0, last = 0;
	MatGetOwnershipRange(matrix, &first, &last);
	std::vector<double> previous = state.node_pressure;
	if (previous.size() != static_cast<std::size_t>(n)) previous.assign(static_cast<std::size_t>(n), flow.wall.reference_pressure);
	const auto compliance = OneDNodeCompliance(graph, flow, previous);
	for (PetscInt row = first; row < last; ++row) {
		const auto* outlet = OneDOutletAtGraphNode(graph, state.outlets, static_cast<int>(row));
		if (outlet && outlet->kind == OneDOutletKind::Pressure) {
			MatSetValue(matrix, row, row, 1.0, ADD_VALUES);
			VecSetValue(rhs, row, outlet->pressure, INSERT_VALUES);
			continue;
		}
		double diagonal = compliance[static_cast<std::size_t>(row)]/dt;
		double value = diagonal*previous[static_cast<std::size_t>(row)];
		for (const auto& edge : graph.edges) {
			if (edge.from != row && edge.to != row) continue;
			const double resistance = 8.0*OneDPi*flow.dynamic_viscosity*edge.length/(edge.area0*edge.area0);
			const int other = edge.from == row ? edge.to : edge.from;
			diagonal += 1.0/resistance;
			MatSetValue(matrix, row, other, -1.0/resistance, ADD_VALUES);
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
	}
	MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(rhs); VecAssemblyEnd(rhs);
	KSPCreate(PETSC_COMM_WORLD, &solver);
	KSPSetOperators(solver, matrix, matrix);
	KSPSetFromOptions(solver);
	KSPSolve(solver, rhs, solution);
	KSPConvergedReason reason;
	KSPGetConvergedReason(solver, &reason);
	if (reason <= 0) throw std::runtime_error("PETSc pressure_network did not converge");
	std::vector<double> pressure;
	OneDGetVectorAll(solution, pressure);
	std::vector<double> edge_flow(graph.edges.size());
	for (std::size_t i = 0; i < graph.edges.size(); ++i) {
		const auto& edge = graph.edges[i];
		const double resistance = 8.0*OneDPi*flow.dynamic_viscosity*edge.length/(edge.area0*edge.area0);
		edge_flow[i] = (pressure[static_cast<std::size_t>(edge.from)]-pressure[static_cast<std::size_t>(edge.to)])/resistance;
	}
	OneDUpdateStateFromImplicitSolution(network, graph, flow, pressure, edge_flow, state, dt);
	KSPDestroy(&solver); VecDestroy(&solution); VecDestroy(&rhs); MatDestroy(&matrix);
}

inline void SolveOneDLinearizedAQPetsc(const OneDNetwork& network,
	const OneDFlowSystemDefinition& flow, OneDFlowState& state,
	double inlet_flow, double dt, bool expand_cells)
{
	const auto graph = BuildOneDImplicitGraph(network, expand_cells);
	const PetscInt unknowns = graph.nodes+static_cast<int>(graph.edges.size());
	Mat matrix = nullptr;
	Vec rhs = nullptr, solution = nullptr;
	KSP solver = nullptr;
	MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, unknowns, unknowns,
		10, nullptr, 10, nullptr, &matrix);
	VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, unknowns, &rhs);
	VecDuplicate(rhs, &solution);
	std::vector<double> previous_pressure(static_cast<std::size_t>(graph.nodes), flow.wall.reference_pressure);
	for (int i = 0; i < graph.original_nodes && i < static_cast<int>(state.node_pressure.size()); ++i)
		previous_pressure[static_cast<std::size_t>(i)] = state.node_pressure[static_cast<std::size_t>(i)];
	const auto compliance = OneDNodeCompliance(graph, flow, previous_pressure);
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
			for (std::size_t edge_index = 0; edge_index < graph.edges.size(); ++edge_index) {
				const auto& edge = graph.edges[edge_index];
				if (edge.from == row) MatSetValue(matrix, row, graph.nodes+edge_index, 1.0, ADD_VALUES);
				if (edge.to == row) MatSetValue(matrix, row, graph.nodes+edge_index, -1.0, ADD_VALUES);
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
			const double resistance = 8.0*OneDPi*flow.dynamic_viscosity*edge.length/(edge.area0*edge.area0);
			const double inertance = flow.density*edge.length/edge.area0;
			double previous_flow = 0.0;
			if (expand_cells) {
				const auto& segment = network.segments[static_cast<std::size_t>(edge.segment)];
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
	if (reason <= 0) throw std::runtime_error("PETSc linearized_aq did not converge");
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
		double value = inertance*(q-context.old_flow[edge_index])/context.dt+resistance*q+p_to-p_from;
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
	VecGetOwnershipRange(residual, &first, &last);
	for (PetscInt row = first; row < last; ++row) VecSetValue(residual, row, values[static_cast<std::size_t>(row)], INSERT_VALUES);
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
			for (std::size_t edge = 0; edge < graph.edges.size(); ++edge) {
				if (graph.edges[edge].from == row) MatSetValue(jacobian, row, graph.nodes+edge, 1.0, INSERT_VALUES);
				if (graph.edges[edge].to == row) MatSetValue(jacobian, row, graph.nodes+edge, -1.0, INSERT_VALUES);
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
	VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, unknowns, &solution);
	VecDuplicate(solution, &residual);
	MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, unknowns, unknowns,
		10, nullptr, 10, nullptr, &jacobian);
	OneDNonlinearContext context;
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
	context.compliance = OneDNodeCompliance(graph, flow, context.old_pressure);
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
	if (reason <= 0) throw std::runtime_error("PETSc nonlinear 1d solve did not converge");
	std::vector<double> values; OneDGetVectorAll(solution, values);
	std::vector<double> pressure(values.begin(), values.begin()+graph.nodes);
	std::vector<double> edge_flow(values.begin()+graph.nodes, values.end());
	OneDUpdateStateFromImplicitSolution(network, graph, flow, pressure, edge_flow, state, dt);
	SNESDestroy(&solver); MatDestroy(&jacobian); VecDestroy(&residual); VecDestroy(&solution);
}

inline void AdvanceImplicitOneD(const OneDNetwork& network,
	const OneDFlowSystemDefinition& flow, OneDFlowState& state,
	double inlet_flow, double dt)
{
	state.inlet_flow = inlet_flow;
	if (flow.formulation == OneDImplicitFormulation::PressureNetwork)
		SolveOneDPressureNetworkPetsc(network, flow, state, inlet_flow, dt);
	else if (flow.formulation == OneDImplicitFormulation::LinearizedAQ)
		SolveOneDLinearizedAQPetsc(network, flow, state, inlet_flow, dt, false);
	else if (flow.formulation == OneDImplicitFormulation::NonlinearAQ)
		SolveOneDNonlinearAQPetsc(network, flow, state, inlet_flow, dt, false);
	else SolveOneDNonlinearAQPetsc(network, flow, state, inlet_flow, dt, true);
}

} // namespace iga

#endif
