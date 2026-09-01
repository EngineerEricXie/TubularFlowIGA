#ifndef IGA_PRESSURE_FLOW_COMPONENT_EXECUTOR_HPP
#define IGA_PRESSURE_FLOW_COMPONENT_EXECUTOR_HPP

#include "DomainRuntimeRegistry.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <functional>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

enum class PressureFlowIterationMethod { Explicit, Fixed, Aitken };

struct PressureFlowExecutionControls {
	PressureFlowIterationMethod method = PressureFlowIterationMethod::Explicit;
	int maximum_iterations = 1;
	double pressure_relative_tolerance = 1.0e-6;
	double pressure_reference_pa = 1.0;
	double flow_relative_tolerance = 1.0e-10;
	double relaxation_factor = 0.5;
	double minimum_relaxation = 0.05;
	double maximum_relaxation = 1.0;
};

inline void ValidatePressureFlowExecutionControls(const PressureFlowExecutionControls& controls)
{
	if (controls.method != PressureFlowIterationMethod::Explicit
		&& controls.method != PressureFlowIterationMethod::Fixed
		&& controls.method != PressureFlowIterationMethod::Aitken)
		throw std::runtime_error("pressure-flow executor has an unknown iteration method");
	if (controls.maximum_iterations < 1
		|| !(controls.pressure_relative_tolerance > 0.0)
		|| !std::isfinite(controls.pressure_relative_tolerance)
		|| !(controls.pressure_reference_pa > 0.0)
		|| !std::isfinite(controls.pressure_reference_pa)
		|| !(controls.flow_relative_tolerance > 0.0)
		|| !std::isfinite(controls.flow_relative_tolerance)
		|| !(controls.relaxation_factor > 0.0) || controls.relaxation_factor > 1.0
		|| !std::isfinite(controls.relaxation_factor)
		|| !(controls.minimum_relaxation > 0.0)
		|| controls.minimum_relaxation > controls.relaxation_factor
		|| controls.relaxation_factor > controls.maximum_relaxation
		|| controls.maximum_relaxation > 1.0
		|| !std::isfinite(controls.minimum_relaxation)
		|| !std::isfinite(controls.maximum_relaxation))
		throw std::runtime_error("pressure-flow executor controls are invalid");
	if (controls.method == PressureFlowIterationMethod::Explicit
		&& controls.maximum_iterations != 1)
		throw std::runtime_error("explicit pressure-flow execution requires exactly one iteration");
}

struct PressureFlowEdgeState {
	std::string edge_id;
	double applied_pressure_pa = 0.0;
	double measured_pressure_pa = 0.0;
	double pressure_residual_pa = 0.0;
	double normalized_pressure_residual = 0.0;
	double first_outward_flow_m3_s = 0.0;
	double second_outward_flow_m3_s = 0.0;
	double flow_residual_m3_s = 0.0;
	double normalized_flow_residual = 0.0;
};

struct PressureFlowIterationState {
	int iteration = 0;
	double applied_relaxation = 0.0;
	bool converged = false;
	std::vector<PressureFlowEdgeState> edges;
};

struct PressureFlowStepResult {
	std::vector<PressureFlowIterationState> iterations;
	std::map<PortRef, PortState> accepted_ports;
};

class PressureFlowConvergenceError : public std::runtime_error {
public:
	PressureFlowConvergenceError(std::string message, PressureFlowStepResult result)
		: std::runtime_error(std::move(message)), result_(std::move(result)) {}

	const PressureFlowStepResult& Result() const noexcept { return result_; }

private:
	PressureFlowStepResult result_;
};

class PressureFlowComponentExecutor {
public:
	PressureFlowComponentExecutor(DomainRuntimeRegistry& registry,
		std::string start_domain_id, PressureFlowExecutionControls controls)
		: registry_(registry), plan_(MakeSequentialPlan(registry.Graph(), start_domain_id)),
		  controls_(controls)
	{
		ValidatePressureFlowExecutionControls(controls_);
		for (std::size_t i = 0; i < plan_.edge_ids.size(); ++i)
			interfaces_.push_back(ResolveInterface(plan_.edge_ids[i],
				plan_.domain_ids[i], plan_.domain_ids[i+1]));
	}

	const SequentialCouplingPlan& Plan() const noexcept { return plan_; }

	PressureFlowStepResult Advance(const DomainStepContext& step,
		const std::map<std::string, double>& initial_pressure_pa,
		const std::function<void(const PressureFlowStepResult&)>& before_commit = {})
	{
		step.Validate();
		std::vector<double> pressure;
		pressure.reserve(interfaces_.size());
		for (const auto& interface : interfaces_) {
			const auto found = initial_pressure_pa.find(interface.edge_id);
			if (found == initial_pressure_pa.end() || !std::isfinite(found->second))
				throw std::runtime_error("pressure-flow executor requires a finite initial pressure for edge '"
					+interface.edge_id+"'");
			pressure.push_back(found->second);
		}
		PressureFlowStepResult result;
		try {
			for (const auto& domain_id : plan_.domain_ids)
				registry_.Runtime(domain_id).BeginStep(step);
			std::vector<double> previous_residual;
			double aitken_relaxation = controls_.relaxation_factor;
			for (int iteration = 1; iteration <= controls_.maximum_iterations; ++iteration) {
				if (iteration > 1)
					for (auto domain = plan_.domain_ids.rbegin(); domain != plan_.domain_ids.rend(); ++domain)
						registry_.Runtime(*domain).RollbackTrial();
				for (std::size_t domain_index = 0; domain_index < plan_.domain_ids.size(); ++domain_index) {
					auto& runtime = registry_.Runtime(plan_.domain_ids[domain_index]);
					if (domain_index < interfaces_.size()) {
						PortBoundaryData input;
						input.time_s = step.EndTime();
						input.mean_pressure_pa = pressure[domain_index];
						runtime.SetPortInput(interfaces_[domain_index].pressure_receiver.port_id, input);
					}
					runtime.SolveTrial();
					if (domain_index < interfaces_.size()) {
						const auto& interface = interfaces_[domain_index];
						const auto state = PortStateFor(interface.flow_provider);
						if (!state.outward_flow_m3_s)
							throw std::runtime_error("pressure-flow provider omitted outward flow on edge '"
								+interface.edge_id+"'");
						PortBoundaryData input;
						input.time_s = step.EndTime();
						input.outward_flow_m3_s = -*state.outward_flow_m3_s;
						registry_.Runtime(interface.flow_receiver.domain_id).SetPortInput(
							interface.flow_receiver.port_id, input);
					}
				}

				PressureFlowIterationState iteration_state;
				iteration_state.iteration = iteration;
				iteration_state.applied_relaxation = aitken_relaxation;
				std::vector<double> residual;
				bool converged = controls_.method != PressureFlowIterationMethod::Explicit;
				for (std::size_t edge_index = 0; edge_index < interfaces_.size(); ++edge_index) {
					const auto& interface = interfaces_[edge_index];
					const auto first = PortStateFor(interface.first);
					const auto second = PortStateFor(interface.second);
					const auto measured = PortStateFor(interface.pressure_provider);
					if (!first.outward_flow_m3_s || !second.outward_flow_m3_s
						|| !measured.mean_pressure_pa)
						throw std::runtime_error("pressure-flow edge state is incomplete for edge '"
							+interface.edge_id+"'");
					PressureFlowEdgeState edge_state;
					edge_state.edge_id = interface.edge_id;
					edge_state.applied_pressure_pa = pressure[edge_index];
					edge_state.measured_pressure_pa = *measured.mean_pressure_pa;
					edge_state.pressure_residual_pa = edge_state.measured_pressure_pa-pressure[edge_index];
					edge_state.normalized_pressure_residual = std::abs(edge_state.pressure_residual_pa)
						/std::max({controls_.pressure_reference_pa,
							std::abs(edge_state.applied_pressure_pa),
							std::abs(edge_state.measured_pressure_pa)});
					edge_state.first_outward_flow_m3_s = *first.outward_flow_m3_s;
					edge_state.second_outward_flow_m3_s = *second.outward_flow_m3_s;
					edge_state.flow_residual_m3_s = edge_state.first_outward_flow_m3_s
						+edge_state.second_outward_flow_m3_s;
					const double flow_scale = std::max(std::abs(edge_state.first_outward_flow_m3_s),
						std::abs(edge_state.second_outward_flow_m3_s));
					edge_state.normalized_flow_residual = flow_scale == 0.0 ? 0.0
						: std::abs(edge_state.flow_residual_m3_s)/flow_scale;
					residual.push_back(edge_state.pressure_residual_pa);
					if (edge_state.normalized_pressure_residual > controls_.pressure_relative_tolerance
						|| edge_state.normalized_flow_residual > controls_.flow_relative_tolerance)
						converged = false;
					iteration_state.edges.push_back(edge_state);
					result.accepted_ports[interface.first] = first;
					result.accepted_ports[interface.second] = second;
				}
				iteration_state.converged = converged;
				result.iterations.push_back(iteration_state);
				if (controls_.method == PressureFlowIterationMethod::Explicit || converged) break;
				if (iteration == controls_.maximum_iterations)
					throw PressureFlowConvergenceError(
						"pressure-flow component failed to converge", result);
				if (controls_.method == PressureFlowIterationMethod::Fixed) {
					for (std::size_t i = 0; i < pressure.size(); ++i)
						pressure[i] += controls_.relaxation_factor*residual[i];
				} else {
					if (!previous_residual.empty()) {
						double numerator = 0.0;
						double denominator = 0.0;
						double scale = controls_.pressure_reference_pa;
						for (std::size_t i = 0; i < residual.size(); ++i)
							scale = std::max({scale, std::abs(residual[i]), std::abs(previous_residual[i])});
						for (std::size_t i = 0; i < residual.size(); ++i) {
							const double difference = (residual[i]-previous_residual[i])/scale;
							numerator += previous_residual[i]/scale*difference;
							denominator += difference*difference;
						}
						if (denominator > 64.0*std::numeric_limits<double>::epsilon()) {
							const double candidate = -aitken_relaxation*numerator/denominator;
							if (std::isfinite(candidate)) aitken_relaxation = std::max(
								controls_.minimum_relaxation,
								std::min(controls_.maximum_relaxation, candidate));
						}
					}
					for (std::size_t i = 0; i < pressure.size(); ++i)
						pressure[i] += aitken_relaxation*residual[i];
					previous_residual = residual;
				}
			}
			for (const auto& domain_id : plan_.domain_ids)
				for (const auto& port : registry_.Graph().Domain(domain_id).ports)
					result.accepted_ports[{domain_id, port.id}]
						= PortStateFor({domain_id, port.id});
			if (before_commit) before_commit(result);
			for (const auto& domain_id : plan_.domain_ids)
				registry_.Runtime(domain_id).PrepareCommitStep();
			for (const auto& domain_id : plan_.domain_ids)
				registry_.Runtime(domain_id).FinalizeCommitStep();
			return result;
		} catch (...) {
			const auto primary = std::current_exception();
			std::string cleanup;
			try { cleanup = AbortAll(); }
			catch (...) { cleanup = "cleanup reporting failed: "+ExceptionText(std::current_exception()); }
			if (cleanup.empty()) std::rethrow_exception(primary);
			throw std::runtime_error(ExceptionText(primary)+"; abort failures: "+cleanup);
		}
	}

private:
	struct Interface {
		std::string edge_id;
		PortRef first;
		PortRef second;
		PortRef flow_provider;
		PortRef flow_receiver;
		PortRef pressure_provider;
		PortRef pressure_receiver;
	};

	Interface ResolveInterface(const std::string& edge_id, const std::string& earlier,
		const std::string& later) const
	{
		const auto& edge = registry_.Graph().Edge(edge_id);
		const auto earlier_ref = edge.first.domain_id == earlier ? edge.first : edge.second;
		const auto later_ref = edge.first.domain_id == later ? edge.first : edge.second;
		if (earlier_ref.domain_id != earlier || later_ref.domain_id != later)
			throw std::runtime_error("sequential edge does not join adjacent plan domains");
		const auto& earlier_port = registry_.Graph().Port(earlier_ref);
		const auto& later_port = registry_.Graph().Port(later_ref);
		if (!earlier_port.requires.count(PortQuantity::MeanPressure)
			|| !later_port.requires.count(PortQuantity::FlowRate))
			throw std::runtime_error("sequential pressure-flow plan requires pressure receivers before flow receivers");
		return {edge_id, edge.first, edge.second, earlier_ref, later_ref, later_ref, earlier_ref};
	}

	PortState PortStateFor(const PortRef& reference) const
	{
		auto state = registry_.Runtime(reference.domain_id).GetPortState(reference.port_id);
		ValidatePortState(state);
		return state;
	}

	std::string AbortAll()
	{
		std::ostringstream errors;
		bool first = true;
		for (auto domain = plan_.domain_ids.rbegin(); domain != plan_.domain_ids.rend(); ++domain) {
			try { registry_.Runtime(*domain).AbortStep(); }
			catch (...) {
				if (!first) errors << "; ";
				first = false;
				errors << *domain << ": " << ExceptionText(std::current_exception());
			}
		}
		return errors.str();
	}

	static std::string ExceptionText(const std::exception_ptr& exception)
	{
		try { if (exception) std::rethrow_exception(exception); }
		catch (const std::exception& error) { return error.what(); }
		catch (...) { return "non-standard exception"; }
		return "unknown exception";
	}

	DomainRuntimeRegistry& registry_;
	SequentialCouplingPlan plan_;
	PressureFlowExecutionControls controls_;
	std::vector<Interface> interfaces_;
};

} // namespace iga

#endif
