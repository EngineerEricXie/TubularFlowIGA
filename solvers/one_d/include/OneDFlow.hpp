#ifndef IGA_ONE_D_FLOW_HPP
#define IGA_ONE_D_FLOW_HPP

#include "OneDNetwork.hpp"
#include "TemporalFunction.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

struct OneDOutletState {
	int node = -1;
	OneDOutletKind kind = OneDOutletKind::Pressure;
	double pressure = 0.0;
	double resistance = 0.0;
	double proximal_resistance = 0.0;
	double distal_resistance = 0.0;
	double capacitance = 0.0;
	double reference_pressure = 0.0;
	double capacitor_pressure = 0.0;
	double flow = 0.0;
};

struct OneDInletState {
	std::string quantity = "flow_rate";
	double value = 0.0;
	std::string waveform;
};

struct OneDFlowState {
	std::vector<double> area;
	std::vector<double> flow;
	std::vector<double> pressure;
	std::vector<double> node_pressure;
	std::vector<double> segment_flow;
	std::vector<OneDOutletState> outlets;
	int completed_step = 0;
	double physical_time = 0.0;
	long long internal_substeps = 0;
	double inlet_flow = 0.0;
};

inline double OneDWallStiffness(double radius0, const OneDWallDefinition& wall)
{
	if (wall.model == OneDWallModel::Olufsen)
		return (4.0/3.0)*(wall.olufsen_k1*std::exp(wall.olufsen_k2*radius0)+wall.olufsen_k3);
	return 4.0*wall.young_modulus*wall.thickness_ratio/3.0;
}

inline double OneDPressureFromArea(double area, double area0, double radius0,
	const OneDWallDefinition& wall)
{
	if (!(area > 0.0) || !(area0 > 0.0)) throw std::runtime_error("1d wall law requires positive area");
	const double stiffness = OneDWallStiffness(radius0, wall);
	if (wall.model == OneDWallModel::Olufsen)
		return wall.reference_pressure + stiffness*(1.0-std::sqrt(area0/area));
	return wall.reference_pressure + stiffness*(std::sqrt(area/area0)-1.0);
}

inline double OneDAreaFromPressure(double pressure, double area0, double radius0,
	const OneDWallDefinition& wall)
{
	const double stiffness = OneDWallStiffness(radius0, wall);
	if (wall.model == OneDWallModel::Olufsen) {
		const double denominator = 1.0-(pressure-wall.reference_pressure)/stiffness;
		if (!(denominator > 0.0)) throw std::runtime_error("Olufsen pressure exceeds invertible wall-law range");
		return area0/(denominator*denominator);
	}
	const double scale = 1.0+(pressure-wall.reference_pressure)/stiffness;
	if (!(scale > 0.0)) throw std::runtime_error("linear wall-law pressure produces non-positive radius");
	return area0*scale*scale;
}

inline double OneDWallAreaDerivative(double pressure, double area0, double radius0,
	const OneDWallDefinition& wall)
{
	const double stiffness = OneDWallStiffness(radius0, wall);
	if (wall.model == OneDWallModel::Olufsen) {
		const double denominator = 1.0-(pressure-wall.reference_pressure)/stiffness;
		if (!(denominator > 0.0)) throw std::runtime_error("Olufsen derivative is outside wall-law range");
		return 2.0*area0/(stiffness*denominator*denominator*denominator);
	}
	const double scale = 1.0+(pressure-wall.reference_pressure)/stiffness;
	return 2.0*area0*scale/stiffness;
}

inline double OneDWaveSpeed(double area, double area0, double radius0,
	const OneDWallDefinition& wall, double density)
{
	const double stiffness = OneDWallStiffness(radius0, wall);
	double dp_da = 0.0;
	if (wall.model == OneDWallModel::Olufsen)
		dp_da = 0.5*stiffness*std::sqrt(area0)*std::pow(area, -1.5);
	else dp_da = stiffness/(2.0*std::sqrt(area0*area));
	return std::sqrt(std::max(0.0, area*dp_da/density));
}

inline double OneDPressurePotential(double area, double area0, double radius0,
	const OneDWallDefinition& wall, double density)
{
	const double stiffness = OneDWallStiffness(radius0, wall);
	if (wall.model == OneDWallModel::Olufsen)
		return stiffness*std::sqrt(area0*area)/density;
	return stiffness*std::pow(area, 1.5)/(3.0*density*std::sqrt(area0));
}

inline std::vector<double> SolveDenseSystem(std::vector<double> matrix,
	std::vector<double> rhs)
{
	const int n = static_cast<int>(rhs.size());
	if (matrix.size() != static_cast<std::size_t>(n*n))
		throw std::runtime_error("dense system dimensions do not agree");
	for (int column = 0; column < n; ++column) {
		int pivot = column;
		for (int row = column+1; row < n; ++row)
			if (std::abs(matrix[static_cast<std::size_t>(row*n+column)])
				> std::abs(matrix[static_cast<std::size_t>(pivot*n+column)])) pivot = row;
		if (std::abs(matrix[static_cast<std::size_t>(pivot*n+column)]) < 1.0e-30)
			throw std::runtime_error("singular 1d network system");
		if (pivot != column) {
			for (int item = column; item < n; ++item)
				std::swap(matrix[static_cast<std::size_t>(column*n+item)],
					matrix[static_cast<std::size_t>(pivot*n+item)]);
			std::swap(rhs[static_cast<std::size_t>(column)], rhs[static_cast<std::size_t>(pivot)]);
		}
		const double diagonal = matrix[static_cast<std::size_t>(column*n+column)];
		for (int item = column; item < n; ++item)
			matrix[static_cast<std::size_t>(column*n+item)] /= diagonal;
		rhs[static_cast<std::size_t>(column)] /= diagonal;
		for (int row = 0; row < n; ++row) {
			if (row == column) continue;
			const double factor = matrix[static_cast<std::size_t>(row*n+column)];
			if (factor == 0.0) continue;
			for (int item = column; item < n; ++item)
				matrix[static_cast<std::size_t>(row*n+item)] -= factor*matrix[static_cast<std::size_t>(column*n+item)];
			rhs[static_cast<std::size_t>(row)] -= factor*rhs[static_cast<std::size_t>(column)];
		}
	}
	return rhs;
}

inline OneDInletState ResolveOneDInlet(const OneDConfiguration& configuration)
{
	for (const auto& boundary : configuration.boundaries) {
		if (boundary.role != "inlet") continue;
		for (const auto& condition : boundary.conditions) {
			if (condition.type != "dirichlet") continue;
			if (condition.field != "flow_rate" && condition.field != "velocity") continue;
			OneDInletState inlet;
			inlet.quantity = condition.quantity.empty()
				? (condition.field == "velocity" ? "centerline_velocity" : "flow_rate")
				: condition.quantity;
			if (inlet.quantity != "flow_rate" && inlet.quantity != "centerline_velocity")
				throw std::runtime_error("1d inlet quantity must be flow_rate or centerline_velocity");
			inlet.value = condition.value;
			inlet.waveform = condition.waveform;
			return inlet;
		}
	}
	throw std::runtime_error("1d inlet requires a flow_rate or velocity Dirichlet condition");
}

inline std::vector<OneDOutletState> ResolveOneDOutlets(
	const OneDConfiguration& configuration, const OneDNetwork& network)
{
	std::map<int, OneDBoundaryCondition> conditions;
	const OneDBoundaryCondition* default_condition = nullptr;
	for (const auto& boundary : configuration.boundaries) {
		if (boundary.role != "outlet") continue;
		const OneDBoundaryCondition* selected = nullptr;
		for (const auto& condition : boundary.conditions)
			if (condition.field == "pressure" && (condition.type == "pressure"
				|| condition.type == "resistance" || condition.type == "windkessel_rcr")) {
				if (selected) throw std::runtime_error("an outlet boundary may define only one pressure closure");
				selected = &condition;
			}
		if (!selected) continue;
		if (boundary.node_ids.empty()) {
			if (default_condition) throw std::runtime_error("only one default 1d outlet closure is allowed");
			default_condition = selected;
		} else for (const int id : boundary.node_ids) {
			if (!conditions.emplace(network.node_index.at(id), *selected).second)
				throw std::runtime_error("a 1d outlet node has multiple pressure closures");
		}
	}
	std::vector<OneDOutletState> result;
	for (const int node : network.outlet_nodes) {
		OneDBoundaryCondition condition;
		const auto found = conditions.find(node);
		if (found != conditions.end()) condition = found->second;
		else if (default_condition) condition = *default_condition;
		else throw std::runtime_error("every 1d outlet leaf requires a pressure, resistance, or RCR closure");
		OneDOutletState outlet;
		outlet.node = node;
		if (condition.type == "pressure") {
			outlet.kind = OneDOutletKind::Pressure;
			outlet.pressure = condition.value;
			outlet.reference_pressure = condition.value;
			outlet.capacitor_pressure = condition.value;
		} else if (condition.type == "resistance") {
			outlet.kind = OneDOutletKind::Resistance;
			outlet.resistance = condition.resistance;
			outlet.reference_pressure = condition.reference_pressure;
			outlet.capacitor_pressure = condition.reference_pressure;
		} else {
			outlet.kind = OneDOutletKind::WindkesselRcr;
			outlet.proximal_resistance = condition.proximal_resistance;
			outlet.distal_resistance = condition.distal_resistance;
			outlet.capacitance = condition.capacitance;
			outlet.reference_pressure = condition.reference_pressure;
			outlet.capacitor_pressure = condition.initial_pressure;
		}
		result.push_back(outlet);
	}
	return result;
}

inline double EvaluateOneDInlet(const OneDConfiguration& configuration,
	const OneDInletState& inlet, const std::filesystem::path& case_directory,
	double time, double inlet_area)
{
	double value = inlet.value;
	if (!inlet.waveform.empty()) {
		const auto& function = FindOneDTemporalFunction(configuration, inlet.waveform);
		std::vector<TemporalSample> samples;
		const std::vector<TemporalSample>* pointer = nullptr;
		if (function.kind == TemporalFunctionKind::PeriodicTable) {
			samples = ReadTemporalCsv((case_directory/function.file).string(), function.period);
			pointer = &samples;
		}
		value = EvaluateTemporalFunction(function, time, pointer);
	}
	if (inlet.quantity == "centerline_velocity") return 0.5*inlet_area*value;
	return value;
}

inline double OutletEffectiveResistance(const OneDOutletState& outlet, double dt)
{
	if (outlet.kind == OneDOutletKind::Pressure) return 0.0;
	if (outlet.kind == OneDOutletKind::Resistance) return outlet.resistance;
	const double denominator = 1.0+dt/(outlet.distal_resistance*outlet.capacitance);
	return outlet.proximal_resistance + dt/(outlet.capacitance*denominator);
}

inline double OutletEffectivePressure(const OneDOutletState& outlet, double dt)
{
	if (outlet.kind == OneDOutletKind::Pressure) return outlet.pressure;
	if (outlet.kind == OneDOutletKind::Resistance) return outlet.reference_pressure;
	const double denominator = 1.0+dt/(outlet.distal_resistance*outlet.capacitance);
	return (outlet.capacitor_pressure
		+dt*outlet.reference_pressure/(outlet.distal_resistance*outlet.capacitance))/denominator;
}

inline void AdvanceOutletState(OneDOutletState& outlet, double flow, double dt)
{
	outlet.flow = flow;
	if (outlet.kind == OneDOutletKind::Pressure) return;
	if (outlet.kind == OneDOutletKind::Resistance) {
		outlet.pressure = outlet.reference_pressure+outlet.resistance*flow;
		return;
	}
	const double denominator = 1.0+dt/(outlet.distal_resistance*outlet.capacitance);
	outlet.capacitor_pressure = (outlet.capacitor_pressure
		+dt*(flow+outlet.reference_pressure/outlet.distal_resistance)/outlet.capacitance)/denominator;
	outlet.pressure = outlet.proximal_resistance*flow+outlet.capacitor_pressure;
}

inline void SolveRigidOneD(const OneDNetwork& network, const OneDFlowSystemDefinition& flow,
	OneDFlowState& state, double inlet_flow, double dt)
{
	state.inlet_flow = inlet_flow;
	const int n = static_cast<int>(network.nodes.size());
	std::vector<double> matrix(static_cast<std::size_t>(n*n), 0.0);
	std::vector<double> rhs(static_cast<std::size_t>(n), 0.0);
	auto add_conductance = [&](const OneDSegment& segment) {
		const double conductance = 1.0/segment.resistance;
		matrix[static_cast<std::size_t>(segment.parent*n+segment.parent)] += conductance;
		matrix[static_cast<std::size_t>(segment.parent*n+segment.child)] -= conductance;
		matrix[static_cast<std::size_t>(segment.child*n+segment.child)] += conductance;
		matrix[static_cast<std::size_t>(segment.child*n+segment.parent)] -= conductance;
	};
	for (const auto& segment : network.segments) add_conductance(segment);
	rhs[static_cast<std::size_t>(network.root)] += inlet_flow;
	for (auto& outlet : state.outlets) {
		const int node = outlet.node;
		if (outlet.kind == OneDOutletKind::Pressure) {
			for (int column = 0; column < n; ++column)
				matrix[static_cast<std::size_t>(node*n+column)] = 0.0;
			matrix[static_cast<std::size_t>(node*n+node)] = 1.0;
			rhs[static_cast<std::size_t>(node)] = outlet.pressure;
		} else {
			const double resistance = OutletEffectiveResistance(outlet, dt);
			const double pressure = OutletEffectivePressure(outlet, dt);
			matrix[static_cast<std::size_t>(node*n+node)] += 1.0/resistance;
			rhs[static_cast<std::size_t>(node)] += pressure/resistance;
		}
	}
	state.node_pressure = SolveDenseSystem(std::move(matrix), std::move(rhs));
	state.segment_flow.resize(network.segments.size());
	for (const auto& segment : network.segments)
		state.segment_flow[static_cast<std::size_t>(segment.index)] =
			(state.node_pressure[static_cast<std::size_t>(segment.parent)]
			-state.node_pressure[static_cast<std::size_t>(segment.child)])/segment.resistance;
	for (auto& outlet : state.outlets) {
		const int incoming = OneDSegmentIntoNode(network, outlet.node);
		AdvanceOutletState(outlet, state.segment_flow[static_cast<std::size_t>(incoming)], dt);
	}
	state.area.assign(static_cast<std::size_t>(network.cells), 0.0);
	state.flow.assign(static_cast<std::size_t>(network.cells), 0.0);
	state.pressure.assign(static_cast<std::size_t>(network.cells), 0.0);
	for (const auto& segment : network.segments)
		for (int cell = 0; cell < segment.cells; ++cell) {
			const std::size_t index = static_cast<std::size_t>(segment.cell_offset+cell);
			const double fraction = (static_cast<double>(cell)+0.5)/segment.cells;
			state.area[index] = segment.area0;
			state.flow[index] = state.segment_flow[static_cast<std::size_t>(segment.index)];
			state.pressure[index] = (1.0-fraction)*state.node_pressure[static_cast<std::size_t>(segment.parent)]
				+fraction*state.node_pressure[static_cast<std::size_t>(segment.child)];
		}
	(void)flow;
}

struct OneDConservativeState { double area = 0.0; double flow = 0.0; };

inline OneDConservativeState OneDFlux(const OneDConservativeState& state,
	const OneDSegment& segment, const OneDFlowSystemDefinition& flow)
{
	return {state.flow,
		flow.discretization.alpha*state.flow*state.flow/state.area
		+OneDPressurePotential(state.area, segment.area0, segment.radius0,
			flow.wall, flow.density)};
}

inline OneDConservativeState OneDRusanovFlux(const OneDConservativeState& left,
	const OneDConservativeState& right, const OneDSegment& segment,
	const OneDFlowSystemDefinition& flow)
{
	const auto flux_left = OneDFlux(left, segment, flow);
	const auto flux_right = OneDFlux(right, segment, flow);
	const double speed_left = std::abs(left.flow/left.area)
		+OneDWaveSpeed(left.area, segment.area0, segment.radius0, flow.wall, flow.density);
	const double speed_right = std::abs(right.flow/right.area)
		+OneDWaveSpeed(right.area, segment.area0, segment.radius0, flow.wall, flow.density);
	const double speed = std::max(speed_left, speed_right);
	return {0.5*(flux_left.area+flux_right.area)-0.5*speed*(right.area-left.area),
		0.5*(flux_left.flow+flux_right.flow)-0.5*speed*(right.flow-left.flow)};
}

inline double OneDJunctionLossCoefficient(const OneDFlowSystemDefinition& flow,
	const OneDNetwork& network, int node, const OneDSegment& parent,
	const OneDSegment& child, double child_fraction)
{
	if (flow.junctions.loss_model == "none") return 0.0;
	const int node_id = network.nodes[static_cast<std::size_t>(node)].id;
	const auto configured = flow.junctions.node_coefficients.find(node_id);
	if (configured != flow.junctions.node_coefficients.end()) return configured->second;
	if (flow.junctions.loss_model == "constant")
		return flow.junctions.coefficient;
	const double angle_degrees = OneDBranchAngleDegrees(parent, child);
	if (flow.junctions.loss_model == "table") {
		const auto& table = flow.junctions.angle_table;
		if (table.empty()) return flow.junctions.coefficient;
		if (angle_degrees <= table.front().first) return table.front().second;
		if (angle_degrees >= table.back().first) return table.back().second;
		const auto upper = std::upper_bound(table.begin(), table.end(), angle_degrees,
			[](double value, const std::pair<double, double>& row) { return value < row.first; });
		const auto lower = upper-1;
		const double fraction = (angle_degrees-lower->first)/(upper->first-lower->first);
		return lower->second+fraction*(upper->second-lower->second);
	}
	const double angle = angle_degrees*OneDPi/180.0;
	if (flow.junctions.loss_model == "angle_sin2")
		return flow.junctions.coefficient*std::sin(angle)*std::sin(angle);
	const double geometry_factor = std::pow(std::max(child.area0/parent.area0, 1.0e-12), -0.5);
	return flow.junctions.coefficient*(1.0-child_fraction)*(1.0-child_fraction)
		*geometry_factor*std::sin(0.5*angle)*std::sin(0.5*angle);
}

inline void InitializeCompliantOneD(const OneDNetwork& network,
	const OneDFlowSystemDefinition& flow, OneDFlowState& state)
{
	state.area.resize(static_cast<std::size_t>(network.cells));
	state.flow.assign(static_cast<std::size_t>(network.cells), 0.0);
	state.pressure.assign(static_cast<std::size_t>(network.cells), flow.wall.reference_pressure);
	state.node_pressure.assign(network.nodes.size(), flow.wall.reference_pressure);
	state.segment_flow.assign(network.segments.size(), 0.0);
	for (const auto& segment : network.segments)
		std::fill(state.area.begin()+segment.cell_offset,
			state.area.begin()+segment.cell_offset+segment.cells, segment.area0);
}

inline void InitializeCompliantOneDFromRigid(const OneDNetwork& network,
	const OneDFlowSystemDefinition& flow, OneDFlowState& state,
	double inlet_flow, double dt)
{
	OneDFlowState rigid;
	rigid.outlets = state.outlets;
	SolveRigidOneD(network, flow, rigid, inlet_flow, dt);
	InitializeCompliantOneD(network, flow, state);
	state.outlets = rigid.outlets;
	state.node_pressure = rigid.node_pressure;
	state.segment_flow = rigid.segment_flow;
	for (const auto& segment : network.segments)
		for (int cell = 0; cell < segment.cells; ++cell) {
			const auto index = static_cast<std::size_t>(segment.cell_offset+cell);
			const double fraction = (static_cast<double>(cell)+0.5)/segment.cells;
			const double pressure = (1.0-fraction)*rigid.node_pressure[static_cast<std::size_t>(segment.parent)]
				+fraction*rigid.node_pressure[static_cast<std::size_t>(segment.child)];
			state.pressure[index] = pressure;
			state.area[index] = OneDAreaFromPressure(pressure, segment.area0, segment.radius0, flow.wall);
			state.flow[index] = rigid.segment_flow[static_cast<std::size_t>(segment.index)];
		}
}

inline double OneDExplicitStableDt(const OneDNetwork& network,
	const OneDFlowSystemDefinition& flow, const OneDFlowState& state)
{
	double result = std::numeric_limits<double>::infinity();
	#ifdef _OPENMP
	#pragma omp parallel for reduction(min:result) schedule(static) if(network.segments.size() >= 64)
	#endif
	for (long long segment_index = 0;
		segment_index < static_cast<long long>(network.segments.size()); ++segment_index) {
		const auto& segment = network.segments[static_cast<std::size_t>(segment_index)];
		const double dx = segment.length/segment.cells;
		for (int cell = 0; cell < segment.cells; ++cell) {
			const std::size_t index = static_cast<std::size_t>(segment.cell_offset+cell);
			const double speed = std::abs(state.flow[index]/state.area[index])
				+OneDWaveSpeed(state.area[index], segment.area0, segment.radius0, flow.wall, flow.density);
			if (speed > 0.0) result = std::min(result, flow.discretization.cfl*dx/speed);
		}
	}
	return result;
}

inline void AdvanceExplicitOneD(const OneDNetwork& network,
	const OneDFlowSystemDefinition& flow, OneDFlowState& state,
	double inlet_flow, double requested_dt)
{
	state.inlet_flow = inlet_flow;
	double remaining = requested_dt;
	std::vector<double> new_area(state.area.size());
	std::vector<double> new_flow(state.flow.size());
	std::vector<OneDConservativeState> left_boundary(network.segments.size());
	std::vector<OneDConservativeState> right_boundary(network.segments.size());
	std::vector<OneDConservativeState> interfaces(
		static_cast<std::size_t>(network.cells)+network.segments.size());
	while (remaining > 0.0) {
		const double stable = OneDExplicitStableDt(network, flow, state);
		if (!(stable > 0.0) || !std::isfinite(stable)) throw std::runtime_error("invalid explicit 1d CFL time step");
		const double dt = std::min(remaining, stable);
		for (const auto& segment : network.segments) {
			const std::size_t first = static_cast<std::size_t>(segment.cell_offset);
			const std::size_t last = static_cast<std::size_t>(segment.cell_offset+segment.cells-1);
			left_boundary[static_cast<std::size_t>(segment.index)] = {state.area[first], state.flow[first]};
			right_boundary[static_cast<std::size_t>(segment.index)] = {state.area[last], state.flow[last]};
		}
		for (const auto& segment : network.segments)
			if (segment.parent == network.root)
				left_boundary[static_cast<std::size_t>(segment.index)].flow = inlet_flow;
		for (auto& outlet : state.outlets) {
			const int incoming = OneDSegmentIntoNode(network, outlet.node);
			auto boundary = right_boundary[static_cast<std::size_t>(incoming)];
			const double pressure = OutletEffectivePressure(outlet, dt)
				+OutletEffectiveResistance(outlet, dt)*boundary.flow;
			const auto& segment = network.segments[static_cast<std::size_t>(incoming)];
			boundary.area = OneDAreaFromPressure(pressure, segment.area0, segment.radius0, flow.wall);
			right_boundary[static_cast<std::size_t>(incoming)] = boundary;
		}
		for (std::size_t node = 0; node < network.nodes.size(); ++node) {
			if (static_cast<int>(node) == network.root || network.nodes[node].children.empty()) continue;
			const int incoming = OneDSegmentIntoNode(network, static_cast<int>(node));
			const auto outgoing = OneDSegmentsOutOfNode(network, static_cast<int>(node));
			if (outgoing.empty()) continue;
			const auto parent_state = right_boundary[static_cast<std::size_t>(incoming)];
			double weight_sum = 0.0;
			for (const int child : outgoing)
				weight_sum += 1.0/network.segments[static_cast<std::size_t>(child)].resistance;
			double total = parent_state.flow;
			for (const int child : outgoing) {
				const auto& child_segment = network.segments[static_cast<std::size_t>(child)];
				const double fraction = (1.0/child_segment.resistance)/weight_sum;
				auto child_state = left_boundary[static_cast<std::size_t>(child)];
				child_state.flow = total*fraction;
				double pressure = OneDPressureFromArea(parent_state.area,
					network.segments[static_cast<std::size_t>(incoming)].area0,
					network.segments[static_cast<std::size_t>(incoming)].radius0, flow.wall);
				const double parent_velocity = parent_state.flow/parent_state.area;
				const double child_velocity = child_state.flow/child_state.area;
				if (flow.junctions.pressure_balance == "total")
					pressure += 0.5*flow.density*(parent_velocity*parent_velocity
						-child_velocity*child_velocity);
				const double reference_velocity = flow.junctions.reference_velocity == "child"
					? child_velocity : parent_velocity;
				const double k = OneDJunctionLossCoefficient(flow, network, static_cast<int>(node),
					network.segments[static_cast<std::size_t>(incoming)], child_segment, fraction);
				pressure -= 0.5*flow.density*k*reference_velocity*reference_velocity;
				child_state.area = OneDAreaFromPressure(pressure, child_segment.area0, child_segment.radius0, flow.wall);
				left_boundary[static_cast<std::size_t>(child)] = child_state;
			}
			right_boundary[static_cast<std::size_t>(incoming)].flow = total;
		}

		#ifdef _OPENMP
		#pragma omp parallel for schedule(static) if(network.segments.size() >= 64)
		#endif
		for (long long segment_index = 0;
			segment_index < static_cast<long long>(network.segments.size()); ++segment_index) {
			const auto& segment = network.segments[static_cast<std::size_t>(segment_index)];
			const double dx = segment.length/segment.cells;
			const std::size_t interface_offset = static_cast<std::size_t>(
				segment.cell_offset+segment.index);
			const std::size_t first = static_cast<std::size_t>(segment.cell_offset);
			interfaces[interface_offset] = OneDRusanovFlux(left_boundary[static_cast<std::size_t>(segment.index)],
				{state.area[first], state.flow[first]}, segment, flow);
			for (int cell = 1; cell < segment.cells; ++cell) {
				const std::size_t left = static_cast<std::size_t>(segment.cell_offset+cell-1);
				const std::size_t right = left+1;
				interfaces[interface_offset+static_cast<std::size_t>(cell)] = OneDRusanovFlux(
					{state.area[left], state.flow[left]}, {state.area[right], state.flow[right]}, segment, flow);
			}
			const std::size_t last = static_cast<std::size_t>(segment.cell_offset+segment.cells-1);
			interfaces[interface_offset+static_cast<std::size_t>(segment.cells)] = OneDRusanovFlux(
				{state.area[last], state.flow[last]}, right_boundary[static_cast<std::size_t>(segment.index)], segment, flow);
			for (int cell = 0; cell < segment.cells; ++cell) {
				const std::size_t index = static_cast<std::size_t>(segment.cell_offset+cell);
				new_area[index] = state.area[index]-dt/dx*(interfaces[interface_offset+static_cast<std::size_t>(cell+1)].area
					-interfaces[interface_offset+static_cast<std::size_t>(cell)].area);
				new_flow[index] = state.flow[index]-dt/dx*(interfaces[interface_offset+static_cast<std::size_t>(cell+1)].flow
					-interfaces[interface_offset+static_cast<std::size_t>(cell)].flow)
					-dt*(8.0*OneDPi*flow.dynamic_viscosity/flow.density)*state.flow[index]/state.area[index];
				if (!std::isfinite(new_area[index]) || !std::isfinite(new_flow[index])
					|| !(new_area[index] >= segment.area0*flow.discretization.min_area_fraction))
					throw std::runtime_error("explicit 1d update produced a non-physical state");
			}
		}
		state.area.swap(new_area);
		state.flow.swap(new_flow);
		remaining -= dt;
		++state.internal_substeps;
	}
	for (const auto& segment : network.segments) {
		double mean_flow = 0.0;
		for (int cell = 0; cell < segment.cells; ++cell) {
			const std::size_t index = static_cast<std::size_t>(segment.cell_offset+cell);
			state.pressure[index] = OneDPressureFromArea(state.area[index], segment.area0,
				segment.radius0, flow.wall);
			mean_flow += state.flow[index];
		}
		state.segment_flow[static_cast<std::size_t>(segment.index)] = mean_flow/segment.cells;
	}
	state.node_pressure.assign(network.nodes.size(), flow.wall.reference_pressure);
	for (const auto& segment : network.segments) {
		state.node_pressure[static_cast<std::size_t>(segment.parent)] =
			state.pressure[static_cast<std::size_t>(segment.cell_offset)];
		state.node_pressure[static_cast<std::size_t>(segment.child)] =
			state.pressure[static_cast<std::size_t>(segment.cell_offset+segment.cells-1)];
	}
	for (auto& outlet : state.outlets) {
		const int segment = OneDSegmentIntoNode(network, outlet.node);
		const auto& definition = network.segments[static_cast<std::size_t>(segment)];
		const double q = state.flow[static_cast<std::size_t>(definition.cell_offset+definition.cells-1)];
		AdvanceOutletState(outlet, q, requested_dt);
	}
}

} // namespace iga

#endif
