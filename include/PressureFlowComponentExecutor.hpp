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

// Optional execution-boundary policy keeps this executor usable without MPI.
// Distributed callers must provide both callbacks on every participant, with
// identical graph/control inputs. Runtime methods still coordinate their own
// internal collectives: synchronizing their returned outcomes cannot rescue a
// rank stuck inside a runtime. before_commit must perform local work only.
struct PressureFlowExecutionSynchronization {
	std::function<void(const char*, std::exception_ptr)> outcome;
	std::function<bool(bool)> all_converged;
};

class PressureFlowComponentExecutor {
public:
	PressureFlowComponentExecutor(DomainRuntimeRegistry& registry,
		std::string start_domain_id, PressureFlowExecutionControls controls,
		PressureFlowExecutionSynchronization synchronization = {})
		: registry_(registry),
		  plan_(MakeAcyclicPressureFlowPlan(registry.Graph(), start_domain_id)),
		  controls_(controls), synchronization_(std::move(synchronization))
	{
		ValidatePressureFlowExecutionControls(controls_);
		if (bool(synchronization_.outcome) != bool(synchronization_.all_converged))
			throw std::runtime_error("pressure-flow synchronization requires both callbacks");
		interfaces_ = plan_.interfaces;
	}

	const PressureFlowComponentPlan& Plan() const noexcept { return plan_; }

	PressureFlowStepResult Advance(const DomainStepContext& step,
		const std::map<std::string, double>& initial_pressure_pa,
		const std::function<void(const PressureFlowStepResult&)>& before_commit = {})
	{
		std::vector<double> pressure;
		std::vector<std::exception_ptr> abort_errors;
		Stage("pressure-flow step input", [&] {
			step.Validate();
			abort_errors.resize(plan_.domain_order.size());
			pressure.reserve(interfaces_.size());
			for (const auto& interface : interfaces_) {
				const auto found = initial_pressure_pa.find(interface.edge_id);
				if (found == initial_pressure_pa.end() || !std::isfinite(found->second))
					throw std::runtime_error("pressure-flow executor requires a finite initial pressure for edge '"
						+interface.edge_id+"'");
				pressure.push_back(found->second);
			}
		});
		PressureFlowStepResult result;
		try {
			for (const auto& domain_id : plan_.domain_order)
				Stage("pressure-flow begin", [&] { registry_.Runtime(domain_id).BeginStep(step); });
			std::vector<double> previous_residual;
			double aitken_relaxation = controls_.relaxation_factor;
			for (int iteration = 1; iteration <= controls_.maximum_iterations; ++iteration) {
				if (iteration > 1)
					for (auto domain = plan_.domain_order.rbegin(); domain != plan_.domain_order.rend(); ++domain)
						Stage("pressure-flow rollback", [&] { registry_.Runtime(*domain).RollbackTrial(); });
				// Pressure conditions are the iteration unknowns and may live on a
				// downstream port of a domain.  Supply every one before the
				// topological flow sweep so a multi-port 3D/1D runtime sees all of
				// its boundary inputs when its upstream flow arrives.
				for (std::size_t edge_index = 0; edge_index < interfaces_.size(); ++edge_index) {
					PortBoundaryData input;
					input.time_s = step.EndTime();
					input.mean_pressure_pa = pressure[edge_index];
					Stage("pressure-flow pressure input", [&] {
						registry_.Runtime(interfaces_[edge_index].pressure_receiver.domain_id)
							.SetPortInput(interfaces_[edge_index].pressure_receiver.port_id, input);
					});
				}
				for (const auto& domain_id : plan_.domain_order) {
					auto& runtime = registry_.Runtime(domain_id);
					Stage("pressure-flow solve", [&] { runtime.SolveTrial(); });
					for (const auto& interface : interfaces_) {
						if (interface.flow_provider.domain_id != domain_id) continue;
						const auto state = PortStateFor(interface.flow_provider);
						Stage("pressure-flow provider validation", [&] {
							if (!state.outward_flow_m3_s)
								throw std::runtime_error("pressure-flow provider omitted outward flow on edge '"
									+interface.edge_id+"'");
						});
						PortBoundaryData input;
						input.time_s = step.EndTime();
						input.outward_flow_m3_s = -*state.outward_flow_m3_s;
						Stage("pressure-flow flow input", [&] {
							registry_.Runtime(interface.flow_receiver.domain_id).SetPortInput(
								interface.flow_receiver.port_id, input);
						});
					}
				}

				PressureFlowIterationState iteration_state;
				iteration_state.iteration = iteration;
				iteration_state.applied_relaxation = aitken_relaxation;
				std::vector<double> residual;
				bool converged = controls_.method != PressureFlowIterationMethod::Explicit;
				for (std::size_t edge_index = 0; edge_index < interfaces_.size(); ++edge_index) {
					const auto& interface = interfaces_[edge_index];
					const auto& edge = registry_.Graph().Edge(interface.edge_id);
					const auto first = PortStateFor(edge.first);
					const auto second = PortStateFor(edge.second);
					const auto measured = PortStateFor(interface.pressure_provider);
					Stage("pressure-flow edge result", [&] {
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
						result.accepted_ports[edge.first] = first;
						result.accepted_ports[edge.second] = second;
					});
				}
				if (synchronization_.all_converged)
					converged = synchronization_.all_converged(converged);
				Stage("pressure-flow iteration result", [&] {
					iteration_state.converged = converged;
					result.iterations.push_back(iteration_state);
				});
				if (controls_.method == PressureFlowIterationMethod::Explicit || converged) break;
				if (iteration == controls_.maximum_iterations)
					throw PressureFlowConvergenceError(
						"pressure-flow component failed to converge", result);
				Stage("pressure-flow relaxation", [&] {
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
				});
			}
			for (const auto& domain_id : plan_.domain_order)
				for (const auto& port : registry_.Graph().Domain(domain_id).ports) {
					PortRef reference;
					Stage("pressure-flow observation reference", [&] { reference = {domain_id, port.id}; });
					auto state = PortStateFor(reference);
					Stage("pressure-flow accepted port", [&] {
						result.accepted_ports[reference] = std::move(state);
					});
				}
			Stage("pressure-flow before commit", [&] { if (before_commit) before_commit(result); });
			for (const auto& domain_id : plan_.domain_order)
				Stage("pressure-flow prepare commit", [&] { registry_.Runtime(domain_id).PrepareCommitStep(); });
			for (const auto& domain_id : plan_.domain_order)
				registry_.Runtime(domain_id).FinalizeCommitStep();
			return result;
		} catch (...) {
			const auto primary = std::current_exception();
			std::string cleanup;
			try { cleanup = AbortAll(abort_errors); }
			catch (...) { cleanup = "cleanup reporting failed: "+ExceptionText(std::current_exception()); }
			if (cleanup.empty()) std::rethrow_exception(primary);
			throw std::runtime_error(ExceptionText(primary)+"; abort failures: "+cleanup);
		}
	}

private:
	PortState PortStateFor(const PortRef& reference) const
	{
		PortState state;
		Stage("pressure-flow port state", [&] {
			state = registry_.Runtime(reference.domain_id).GetPortState(reference.port_id);
		});
		Stage("pressure-flow port validation", [&] { ValidatePortState(state); });
		return state;
	}

	// Invoke at most one runtime method per stage. Runtime implementations own
	// their internal failure protocol; executor-only callbacks are local.
	template<class Work>
	void Stage(const char* name, Work&& work) const
	{
		std::exception_ptr error;
		try { std::forward<Work>(work)(); }
		catch (...) { error = std::current_exception(); }
		if (synchronization_.outcome) synchronization_.outcome(name, error);
		if (error) std::rethrow_exception(error);
	}

	std::string AbortAll(std::vector<std::exception_ptr>& errors)
	{
		// Storage was prepared before BeginStep. Recording a failure does not
		// allocate and cannot skip later abort calls, even under memory pressure.
		std::size_t index = 0;
		for (auto domain = plan_.domain_order.rbegin(); domain != plan_.domain_order.rend(); ++domain, ++index) {
			try { Stage("pressure-flow abort", [&] { registry_.Runtime(*domain).AbortStep(); }); }
			catch (...) { errors[index] = std::current_exception(); }
		}
		std::string result;
		Stage("pressure-flow abort reporting", [&] {
			std::ostringstream report;
			report.exceptions(std::ios::badbit | std::ios::failbit);
			bool first = true;
			for (std::size_t i = 0; i < errors.size(); ++i) {
				if (!errors[i]) continue;
				if (!first) report << "; ";
				first = false;
				report << plan_.domain_order[errors.size()-1-i] << ": " << ExceptionText(errors[i]);
			}
			result = report.str();
		});
		return result;
	}

	static std::string ExceptionText(const std::exception_ptr& exception)
	{
		try { if (exception) std::rethrow_exception(exception); }
		catch (const std::exception& error) { return error.what(); }
		catch (...) { return "non-standard exception"; }
		return "unknown exception";
	}

	DomainRuntimeRegistry& registry_;
	PressureFlowComponentPlan plan_;
	PressureFlowExecutionControls controls_;
	PressureFlowExecutionSynchronization synchronization_;
	std::vector<PressureFlowInterfacePlan> interfaces_;
};

} // namespace iga

#endif
