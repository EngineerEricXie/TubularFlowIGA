#ifndef IGA_SPECIES_PRESSURE_FLOW_COMPONENT_EXECUTOR_HPP
#define IGA_SPECIES_PRESSURE_FLOW_COMPONENT_EXECUTOR_HPP

#include "PressureFlowComponentExecutor.hpp"

#include <algorithm>
#include <exception>
#include <functional>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

struct SpeciesPressureFlowExecutionControls {
	PressureFlowExecutionControls hydraulic;
	SpeciesRoutingControls routing;
	std::map<std::string, SpeciesAmountTolerance> amount_tolerances;

	void Validate() const
	{
		ValidatePressureFlowExecutionControls(hydraulic);
		routing.Validate();
		for (const auto& tolerance : amount_tolerances) {
			if (tolerance.first.empty())
				throw std::runtime_error("species amount tolerance requires a species id");
			tolerance.second.Validate();
		}
	}
};

struct SpeciesEdgeAmountDiagnostic {
	std::string edge_id;
	std::string species_id;
	PortRef first;
	PortRef second;
	SpeciesDonor donor = SpeciesDonor::First;
	double first_outward_amount = 0.0;
	double second_outward_amount = 0.0;
	double residual = 0.0;
	double normalized_residual = 0.0;
};

struct SpeciesDomainBalanceDiagnostic {
	std::string domain_id;
	std::string species_id;
	SpeciesStepAccounting accounting;
	double recomputed_residual = 0.0;
	double normalized_residual = 0.0;
};

struct SpeciesGlobalBalanceDiagnostic {
	std::string species_id;
	double initial_mass = 0.0;
	double final_mass = 0.0;
	double source_amount = 0.0;
	double outward_amount = 0.0;
	double gross_activity = 0.0;
	double residual = 0.0;
	double normalized_residual = 0.0;
};

struct SpeciesPressureFlowStepResult {
	std::vector<PressureFlowIterationState> hydraulic_iterations;
	std::map<PortRef, PortState> accepted_ports;
	std::vector<std::string> transport_domain_order;
	std::map<std::pair<std::string, std::string>, SpeciesDonor> donor_ownership;
	std::vector<SpeciesEdgeAmountDiagnostic> edge_amounts;
	std::vector<SpeciesDomainBalanceDiagnostic> domain_balances;
	std::map<std::string, SpeciesGlobalBalanceDiagnostic> global_balances;
};

class SpeciesPressureFlowConvergenceError : public std::runtime_error {
public:
	SpeciesPressureFlowConvergenceError(std::string message,
		SpeciesPressureFlowStepResult result)
		: std::runtime_error(std::move(message)), result_(std::move(result)) {}

	const SpeciesPressureFlowStepResult& Result() const noexcept { return result_; }

private:
	SpeciesPressureFlowStepResult result_;
};

class SpeciesPressureFlowComponentExecutor {
public:
	SpeciesPressureFlowComponentExecutor(DomainRuntimeRegistry& registry,
		std::string start_domain_id, SpeciesPressureFlowExecutionControls controls)
		: registry_(registry), plan_(MakeAcyclicPressureFlowPlan(registry.Graph(), start_domain_id)),
		  controls_(std::move(controls))
	{
		controls_.Validate();
		interfaces_ = plan_.interfaces;
		for (const auto& domain_id : plan_.domain_order) {
			auto* staged = dynamic_cast<StagedFlowTransportDomainRuntime*>(
				&registry_.Runtime(domain_id));
			if (!staged)
				throw std::runtime_error("species pressure-flow executor requires staged flow/transport runtime for domain '"
					+domain_id+"'");
			staged_.emplace(domain_id, staged);
		}
		for (const auto& edge : registry_.Graph().Edges())
			for (const auto& species : edge.species)
				if (!controls_.amount_tolerances.count(species))
					throw std::runtime_error("species pressure-flow executor is missing an amount tolerance for species '"
						+species+"'");
	}

	const PressureFlowComponentPlan& Plan() const noexcept { return plan_; }
	const std::map<std::pair<std::string, std::string>, SpeciesDonor>&
	CommittedDonorOwnership() const noexcept { return committed_donors_; }

	SpeciesPressureFlowStepResult Advance(const DomainStepContext& step,
		const std::map<std::string, double>& initial_pressure_pa,
		const std::function<void(const SpeciesPressureFlowStepResult&)>& before_commit = {})
	{
		step.Validate();
		std::vector<double> pressure = InitialPressure(initial_pressure_pa);
		SpeciesPressureFlowStepResult result;
		try {
			for (const auto& domain_id : plan_.domain_order)
				registry_.Runtime(domain_id).BeginStep(step);
			SolveHydraulics(step, pressure, result);
			CaptureAcceptedPorts(result);
			auto candidate_donors = ResolveDonors(result);
			result.donor_ownership = candidate_donors;
			const auto routes = MakeRoutes(candidate_donors);
			result.transport_domain_order = TransportOrder(routes);
			SolveTransport(step, routes, result.transport_domain_order);
			VerifyAmounts(result);
			if (before_commit) before_commit(result);
			for (const auto& domain_id : plan_.domain_order)
				registry_.Runtime(domain_id).PrepareCommitStep();
			for (const auto& domain_id : plan_.domain_order)
				registry_.Runtime(domain_id).FinalizeCommitStep();
			committed_donors_ = std::move(candidate_donors);
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
	struct Route {
		std::string edge_id;
		std::string species_id;
		PortRef donor;
		PortRef receiver;
		SpeciesDonor ownership = SpeciesDonor::First;
	};

	std::vector<double> InitialPressure(const std::map<std::string, double>& values) const
	{
		std::vector<double> result;
		for (const auto& interface : interfaces_) {
			const auto found = values.find(interface.edge_id);
			if (found == values.end() || !std::isfinite(found->second))
				throw std::runtime_error("species pressure-flow executor requires a finite initial pressure for edge '"
					+interface.edge_id+"'");
			result.push_back(found->second);
		}
		return result;
	}

	void SolveHydraulics(const DomainStepContext& step, std::vector<double>& pressure,
		SpeciesPressureFlowStepResult& result)
	{
		std::vector<double> previous_residual;
		double aitken_relaxation = controls_.hydraulic.relaxation_factor;
		for (int iteration = 1; iteration <= controls_.hydraulic.maximum_iterations; ++iteration) {
			if (iteration > 1)
				for (auto domain = plan_.domain_order.rbegin(); domain != plan_.domain_order.rend(); ++domain)
					staged_.at(*domain)->RollbackHydraulicTrial();
			for (std::size_t i = 0; i < interfaces_.size(); ++i) {
				PortBoundaryData input;
				input.time_s = step.EndTime();
				input.mean_pressure_pa = pressure[i];
				registry_.Runtime(interfaces_[i].pressure_receiver.domain_id).SetPortInput(
					interfaces_[i].pressure_receiver.port_id, input);
			}
			for (const auto& domain_id : plan_.domain_order) {
				staged_.at(domain_id)->SolveHydraulicTrial();
				for (const auto& interface : interfaces_) if (interface.flow_provider.domain_id == domain_id) {
					const auto state = HydraulicState(interface.flow_provider);
					if (!state.outward_flow_m3_s)
						throw std::runtime_error("species pressure-flow provider omitted outward flow on edge '"
							+interface.edge_id+"'");
					PortBoundaryData input;
					input.time_s = step.EndTime();
					input.outward_flow_m3_s = -*state.outward_flow_m3_s;
					registry_.Runtime(interface.flow_receiver.domain_id).SetPortInput(
						interface.flow_receiver.port_id, input);
				}
			}
			PressureFlowIterationState state;
			state.iteration = iteration;
			state.applied_relaxation = aitken_relaxation;
			std::vector<double> residual;
			bool converged = controls_.hydraulic.method != PressureFlowIterationMethod::Explicit;
			for (std::size_t i = 0; i < interfaces_.size(); ++i) {
				const auto& interface = interfaces_[i];
				const auto& edge = registry_.Graph().Edge(interface.edge_id);
				const auto first = HydraulicState(edge.first);
				const auto second = HydraulicState(edge.second);
				const auto measured = HydraulicState(interface.pressure_provider);
				if (!first.outward_flow_m3_s || !second.outward_flow_m3_s || !measured.mean_pressure_pa)
					throw std::runtime_error("species pressure-flow edge state is incomplete for edge '"+interface.edge_id+"'");
				PressureFlowEdgeState edge_state;
				edge_state.edge_id = interface.edge_id;
				edge_state.applied_pressure_pa = pressure[i];
				edge_state.measured_pressure_pa = *measured.mean_pressure_pa;
				edge_state.pressure_residual_pa = edge_state.measured_pressure_pa-pressure[i];
				edge_state.normalized_pressure_residual = std::abs(edge_state.pressure_residual_pa)
					/std::max({controls_.hydraulic.pressure_reference_pa,
						std::abs(edge_state.applied_pressure_pa), std::abs(edge_state.measured_pressure_pa)});
				edge_state.first_outward_flow_m3_s = *first.outward_flow_m3_s;
				edge_state.second_outward_flow_m3_s = *second.outward_flow_m3_s;
				edge_state.flow_residual_m3_s = edge_state.first_outward_flow_m3_s
					+edge_state.second_outward_flow_m3_s;
				const double flow_scale = std::max(std::abs(edge_state.first_outward_flow_m3_s),
					std::abs(edge_state.second_outward_flow_m3_s));
				edge_state.normalized_flow_residual = flow_scale == 0.0 ? 0.0
					: std::abs(edge_state.flow_residual_m3_s)/flow_scale;
				residual.push_back(edge_state.pressure_residual_pa);
				if (edge_state.normalized_pressure_residual > controls_.hydraulic.pressure_relative_tolerance
					|| edge_state.normalized_flow_residual > controls_.hydraulic.flow_relative_tolerance)
					converged = false;
				state.edges.push_back(edge_state);
			}
			state.converged = converged;
			result.hydraulic_iterations.push_back(state);
			if (controls_.hydraulic.method == PressureFlowIterationMethod::Explicit || converged) return;
			if (iteration == controls_.hydraulic.maximum_iterations)
				throw SpeciesPressureFlowConvergenceError("species pressure-flow component failed to converge", result);
			if (controls_.hydraulic.method == PressureFlowIterationMethod::Fixed) {
				for (std::size_t i = 0; i < pressure.size(); ++i)
					pressure[i] += controls_.hydraulic.relaxation_factor*residual[i];
			} else {
				if (!previous_residual.empty()) {
					double numerator = 0.0;
					double denominator = 0.0;
					double scale = controls_.hydraulic.pressure_reference_pa;
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
							controls_.hydraulic.minimum_relaxation,
							std::min(controls_.hydraulic.maximum_relaxation, candidate));
					}
				}
				for (std::size_t i = 0; i < pressure.size(); ++i)
					pressure[i] += aitken_relaxation*residual[i];
				previous_residual = residual;
			}
		}
	}

	void CaptureAcceptedPorts(SpeciesPressureFlowStepResult& result) const
	{
		for (const auto& domain_id : plan_.domain_order)
			for (const auto& port : registry_.Graph().Domain(domain_id).ports)
				result.accepted_ports[{domain_id, port.id}] = HydraulicState({domain_id, port.id});
	}

	std::map<std::pair<std::string, std::string>, SpeciesDonor>
	ResolveDonors(const SpeciesPressureFlowStepResult& result) const
	{
		std::map<std::pair<std::string, std::string>, SpeciesDonor> result_donors;
		auto edges = registry_.Graph().Edges();
		std::sort(edges.begin(), edges.end(), [](const CouplingEdge& a, const CouplingEdge& b) {
			return a.id < b.id;
		});
		for (const auto& edge : edges) for (const auto& species : edge.species) {
			const auto& first = result.accepted_ports.at(edge.first);
			const auto& second = result.accepted_ports.at(edge.second);
			if (!first.outward_flow_m3_s || !second.outward_flow_m3_s)
				throw std::runtime_error("species routing requires accepted outward flow on edge '"+edge.id+"'");
			const auto key = std::make_pair(edge.id, species);
			const auto prior = committed_donors_.find(key);
			result_donors.emplace(key, ResolveSpeciesDonor(*first.outward_flow_m3_s,
				*second.outward_flow_m3_s, controls_.routing,
				prior == committed_donors_.end() ? std::nullopt : std::optional<SpeciesDonor>(prior->second)));
		}
		return result_donors;
	}

	std::vector<Route> MakeRoutes(const std::map<std::pair<std::string, std::string>, SpeciesDonor>& donors) const
	{
		std::vector<Route> result;
		for (const auto& item : donors) {
			const auto& edge = registry_.Graph().Edge(item.first.first);
			const bool first = item.second == SpeciesDonor::First;
			result.push_back({edge.id, item.first.second, first ? edge.first : edge.second,
				first ? edge.second : edge.first, item.second});
		}
		std::sort(result.begin(), result.end(), [](const Route& a, const Route& b) {
			return std::tie(a.edge_id, a.species_id, a.donor, a.receiver)
				< std::tie(b.edge_id, b.species_id, b.donor, b.receiver);
		});
		return result;
	}

	std::vector<std::string> TransportOrder(const std::vector<Route>& routes) const
	{
		std::map<std::string, int> indegree;
		std::map<std::string, std::set<std::string>> children;
		for (const auto& domain : registry_.Graph().Domains()) {
			indegree.emplace(domain.first, 0);
			children.emplace(domain.first, std::set<std::string>{});
		}
		for (const auto& route : routes)
			if (route.donor.domain_id != route.receiver.domain_id
				&& children.at(route.donor.domain_id).insert(route.receiver.domain_id).second)
				++indegree.at(route.receiver.domain_id);
		std::set<std::string> ready;
		for (const auto& item : indegree) if (item.second == 0) ready.insert(item.first);
		std::vector<std::string> result;
		while (!ready.empty()) {
			const auto domain = *ready.begin();
			ready.erase(ready.begin());
			result.push_back(domain);
			for (const auto& child : children.at(domain)) if (--indegree.at(child) == 0)
				ready.insert(child);
		}
		if (result.size() != indegree.size())
			throw std::runtime_error("species transport dependencies contain a cycle");
		return result;
	}

	void SolveTransport(const DomainStepContext& step, const std::vector<Route>& routes,
		const std::vector<std::string>& order) const
	{
		for (const auto& domain_id : order) {
			std::map<std::string, std::map<std::string, double>> inputs;
			for (const auto& route : routes) if (route.receiver.domain_id == domain_id) {
				const auto state = staged_.at(route.donor.domain_id)->GetTransportPortState(route.donor.port_id);
				ValidatePortState(state);
				const auto& edge = registry_.Graph().Edge(route.edge_id);
				const auto& required = edge.species;
				if (!state.concentration.count(route.species_id)
					|| !std::isfinite(state.concentration.at(route.species_id)))
					throw std::runtime_error("species donor port omitted a finite concentration for '"+route.species_id+"'");
				inputs[route.receiver.port_id].emplace(route.species_id,
					state.concentration.at(route.species_id));
				if (inputs.at(route.receiver.port_id).size() == required.size())
					staged_.at(domain_id)->SetTransportConcentration(route.receiver.port_id,
						step.EndTime(), inputs.at(route.receiver.port_id));
			}
			staged_.at(domain_id)->SolveTransportTrial();
		}
	}

	void VerifyAmounts(SpeciesPressureFlowStepResult& result) const
	{
		std::map<std::string, std::map<std::string, SpeciesStepAccounting>> accounting;
		for (const auto& domain_id : plan_.domain_order)
			accounting.emplace(domain_id, staged_.at(domain_id)->GetSpeciesStepAccounting());
		for (const auto& edge : registry_.Graph().Edges()) for (const auto& species : edge.species) {
			const PortRef* first_port = &edge.first;
			const PortRef* second_port = &edge.second;
			bool swapped = false;
			if (*second_port < *first_port) {
				std::swap(first_port, second_port);
				swapped = true;
			}
			const auto& first = AccountingFor(accounting, first_port->domain_id, species);
			const auto& second = AccountingFor(accounting, second_port->domain_id, species);
			const double first_amount = PortAmount(first, first_port->port_id, species);
			const double second_amount = PortAmount(second, second_port->port_id, species);
			const double residual = first_amount+second_amount;
			const auto& tolerance = controls_.amount_tolerances.at(species);
			const double scale = std::max({tolerance.reference_amount, std::abs(first_amount), std::abs(second_amount)});
			Gate(species, residual, scale, "edge '"+edge.id+"'");
			auto donor = result.donor_ownership.at({edge.id, species});
			if (swapped)
				donor = donor == SpeciesDonor::First
					? SpeciesDonor::Second : SpeciesDonor::First;
			result.edge_amounts.push_back({edge.id, species, *first_port, *second_port,
				donor, first_amount, second_amount,
				residual, Normalized(residual, scale)});
		}
		std::sort(result.edge_amounts.begin(), result.edge_amounts.end(),
			[](const SpeciesEdgeAmountDiagnostic& first,
				const SpeciesEdgeAmountDiagnostic& second) {
				return first.edge_id == second.edge_id
					? first.species_id < second.species_id
					: first.edge_id < second.edge_id;
			});
		for (const auto& item : accounting) for (const auto& species : UsedSpecies(item.first)) {
			const auto& value = AccountingFor(accounting, item.first, species);
			ValidateAccounting(value, item.first, species);
			const double scale = AccountingScale(value, controls_.amount_tolerances.at(species));
			const double recomputed = RecomputeResidual(value);
			Gate(species, recomputed, scale, "domain '"+item.first+"'");
			Gate(species, value.residual-recomputed, scale,
				"native accounting residual disagreement in domain '"+item.first+"'");
			result.domain_balances.push_back({item.first, species, value, recomputed,
				Normalized(recomputed, scale)});
			auto& global = result.global_balances[species];
			global.species_id = species;
			global.initial_mass += value.initial_mass;
			global.final_mass += value.final_mass;
			global.source_amount += value.source_amount;
			global.gross_activity += std::abs(value.final_mass-value.initial_mass)
				+std::abs(value.source_amount);
			for (const auto& amount : value.outward_port_amount) {
				global.outward_amount += amount.second;
				global.gross_activity += std::abs(amount.second);
			}
		}
		for (auto& global : result.global_balances) {
			global.second.residual = global.second.final_mass-global.second.initial_mass
				+global.second.outward_amount-global.second.source_amount;
			const double scale = std::max(controls_.amount_tolerances.at(global.first).reference_amount,
				global.second.gross_activity);
			Gate(global.first, global.second.residual, scale, "global balance");
			global.second.normalized_residual = Normalized(global.second.residual, scale);
		}
	}

	const SpeciesStepAccounting& AccountingFor(const std::map<std::string,
		std::map<std::string, SpeciesStepAccounting>>& values, const std::string& domain,
		const std::string& species) const
	{
		const auto domain_found = values.find(domain);
		if (domain_found == values.end() || !domain_found->second.count(species))
			throw std::runtime_error("species accounting is missing '"+species+"' for domain '"+domain+"'");
		return domain_found->second.at(species);
	}

	double PortAmount(const SpeciesStepAccounting& value, const std::string& port,
		const std::string& species) const
	{
		const auto found = value.outward_port_amount.find(port);
		if (found == value.outward_port_amount.end() || !std::isfinite(found->second))
			throw std::runtime_error("species accounting is missing a finite port amount for '"+species+"'");
		return found->second;
	}

	std::set<std::string> UsedSpecies(const std::string& domain_id) const
	{
		std::set<std::string> result;
		for (const auto& edge : registry_.Graph().Edges())
			if (edge.first.domain_id == domain_id || edge.second.domain_id == domain_id)
				result.insert(edge.species.begin(), edge.species.end());
		return result;
	}

	static void ValidateAccounting(const SpeciesStepAccounting& value,
		const std::string& domain, const std::string& species)
	{
		if (!std::isfinite(value.initial_mass) || !std::isfinite(value.final_mass)
			|| !std::isfinite(value.source_amount) || !std::isfinite(value.residual))
			throw std::runtime_error("species accounting is non-finite for '"+species+"' in domain '"+domain+"'");
		for (const auto& amount : value.outward_port_amount)
			if (!std::isfinite(amount.second))
				throw std::runtime_error("species accounting contains a non-finite port amount");
	}

	double AccountingScale(const SpeciesStepAccounting& value,
		const SpeciesAmountTolerance& tolerance) const
	{
		double scale = std::max({tolerance.reference_amount,
			std::abs(value.final_mass-value.initial_mass), std::abs(value.source_amount)});
		for (const auto& amount : value.outward_port_amount) scale = std::max(scale, std::abs(amount.second));
		return scale;
	}

	static double RecomputeResidual(const SpeciesStepAccounting& value)
	{
		double residual = value.final_mass-value.initial_mass-value.source_amount;
		for (const auto& amount : value.outward_port_amount) residual += amount.second;
		return residual;
	}

	void Gate(const std::string& species, double residual, double scale,
		const std::string& what) const
	{
		const auto& tolerance = controls_.amount_tolerances.at(species);
		if (std::abs(residual) > tolerance.absolute_tolerance+tolerance.relative_tolerance*scale)
			throw std::runtime_error("species amount balance failed for "+what+" and species '"+species+"'");
	}

	static double Normalized(double residual, double scale)
	{
		return scale == 0.0 ? 0.0 : std::abs(residual)/scale;
	}

	PortState HydraulicState(const PortRef& reference) const
	{
		auto state = staged_.at(reference.domain_id)->GetHydraulicPortState(reference.port_id);
		ValidatePortState(state);
		return state;
	}

	std::string AbortAll() const
	{
		std::ostringstream errors;
		bool first = true;
		for (auto domain = plan_.domain_order.rbegin(); domain != plan_.domain_order.rend(); ++domain) {
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
	PressureFlowComponentPlan plan_;
	SpeciesPressureFlowExecutionControls controls_;
	std::vector<PressureFlowInterfacePlan> interfaces_;
	std::map<std::string, StagedFlowTransportDomainRuntime*> staged_;
	std::map<std::pair<std::string, std::string>, SpeciesDonor> committed_donors_;
};

} // namespace iga

#endif
