#ifndef IGA_GRAPH_ACCEPTED_HISTORY_HPP
#define IGA_GRAPH_ACCEPTED_HISTORY_HPP

#include "CheckpointMetadataCodec.hpp"
#include "SpeciesPressureFlowComponentExecutor.hpp"
#include "ZeroDFlowDomain.hpp"
#include <climits>

namespace iga {
struct GraphAcceptedFlowStep {
	int step = 0;
	double time_s = 0.0;
	int iterations = 0;
	double three_d_mass_m3_s = 0.0;
	double external_outward_flow_m3_s = 0.0;
	std::map<std::string, std::pair<double, double>> three_d_balance;
	std::map<std::string, ZeroDFlowState> zero_d_states;
	std::map<std::string, ZeroDFlowStepAccounting> zero_d_accounting;
	PressureFlowStepResult result;
};
struct GraphAcceptedSpeciesStep {
	int step = 0;
	double time_s = 0.0;
	int hydraulic_iterations = 0;
	double three_d_mass_m3_s = 0.0;
	double external_outward_flow_m3_s = 0.0;
	std::map<std::string, std::pair<double, double>> three_d_balance;
	std::map<PortRef, PortState> transport_ports;
	SpeciesPressureFlowStepResult result;
};
namespace graph_history_detail {
using namespace checkpoint_metadata;
inline std::uint64_t Count(Reader& in)
{
	const auto n = in.Unsigned(); Require(n <= maximum_map_entries, "history collection exceeds limit"); return n;
}
inline void Count(Writer& out, std::size_t n)
{
	Require(n <= maximum_map_entries, "history collection exceeds limit"); out.Unsigned(n);
}
inline int Integer(Reader& in)
{
	const auto n = in.Unsigned(); Require(n <= INT_MAX, "history integer overflows"); return static_cast<int>(n);
}
inline void Reference(Writer& out, const PortRef& ref)
{
	out.Text(ref.domain_id); out.Text(ref.port_id);
}
inline PortRef Reference(Reader& in)
{
	PortRef ref; ref.domain_id = in.Text(); ref.port_id = in.Text();
	Require(!ref.domain_id.empty() && !ref.port_id.empty(), "empty history port"); return ref;
}
inline void Ports(Writer& out, const std::map<PortRef, PortState>& values)
{
	Count(out, values.size()); for (const auto& item : values) { Reference(out, item.first); WritePort(out, item.second); }
}
inline std::map<PortRef, PortState> Ports(Reader& in)
{
	std::map<PortRef, PortState> values; PortRef previous;
	for (auto n = Count(in); n; --n) {
		auto ref = Reference(in); Require(previous < ref, "history ports are unordered or duplicate"); previous = ref;
		auto value = ReadPort(in); values.emplace(std::move(ref), std::move(value));
	}
	return values;
}
inline void Iterations(Writer& out, const std::vector<PressureFlowIterationState>& values)
{
	Count(out, values.size());
	for (const auto& iteration : values) {
		Require(iteration.iteration >= 0, "negative history iteration"); out.Unsigned(iteration.iteration);
		out.Real(iteration.applied_relaxation); out.Unsigned(iteration.converged); Count(out, iteration.edges.size());
		for (const auto& edge : iteration.edges) {
			out.Text(edge.edge_id);
			for (double value : {edge.applied_pressure_pa, edge.measured_pressure_pa, edge.pressure_residual_pa,
				edge.normalized_pressure_residual, edge.first_outward_flow_m3_s, edge.second_outward_flow_m3_s,
				edge.flow_residual_m3_s, edge.normalized_flow_residual}) out.Real(value);
		}
	}
}
inline std::vector<PressureFlowIterationState> Iterations(Reader& in)
{
	std::vector<PressureFlowIterationState> values;
	for (auto n = Count(in); n; --n) {
		PressureFlowIterationState iteration; iteration.iteration = Integer(in); iteration.applied_relaxation = in.Real();
		const auto converged = in.Unsigned(); Require(converged <= 1, "invalid history convergence flag"); iteration.converged = converged;
		for (auto edges = Count(in); edges; --edges) {
			PressureFlowEdgeState edge; edge.edge_id = in.Text(); edge.applied_pressure_pa = in.Real();
			edge.measured_pressure_pa = in.Real(); edge.pressure_residual_pa = in.Real(); edge.normalized_pressure_residual = in.Real();
			edge.first_outward_flow_m3_s = in.Real(); edge.second_outward_flow_m3_s = in.Real(); edge.flow_residual_m3_s = in.Real();
			edge.normalized_flow_residual = in.Real(); iteration.edges.push_back(std::move(edge));
		}
		values.push_back(std::move(iteration));
	}
	Require(!values.empty(), "history has no accepted iteration"); return values;
}
template<class Step> void Common(Writer& out, const Step& step, int iterations)
{
	Require(step.step > 0 && iterations > 0, "invalid history counters"); out.Unsigned(step.step); out.Real(step.time_s);
	out.Unsigned(iterations); out.Real(step.three_d_mass_m3_s); out.Real(step.external_outward_flow_m3_s);
	Count(out, step.three_d_balance.size());
	for (const auto& balance : step.three_d_balance) { out.Text(balance.first); out.Real(balance.second.first); out.Real(balance.second.second); }
}
template<class Step> int Common(Reader& in, Step& step)
{
	step.step = Integer(in); step.time_s = in.Real(); const int iterations = Integer(in);
	Require(step.step > 0 && iterations > 0, "invalid history counters");
	step.three_d_mass_m3_s = in.Real(); step.external_outward_flow_m3_s = in.Real(); std::string previous;
	for (auto n = Count(in); n; --n) {
		auto key = in.Text(); Require(previous < key, "invalid history balance catalog"); previous = key;
		const double sum = in.Real(), absolute = in.Real(); step.three_d_balance.emplace(std::move(key), std::make_pair(sum, absolute));
	}
	return iterations;
}
inline void Accounting(Writer& out, const ZeroDFlowStepAccounting& value)
{
	for (double scalar : {value.initial_stored_volume_m3, value.final_stored_volume_m3, value.prescribed_source_amount_m3,
		value.distal_sink_amount_m3, value.outward_graph_port_amount_m3, value.residual_m3}) out.Real(scalar);
}
inline ZeroDFlowStepAccounting Accounting(Reader& in)
{
	ZeroDFlowStepAccounting value; value.initial_stored_volume_m3 = in.Real(); value.final_stored_volume_m3 = in.Real();
	value.prescribed_source_amount_m3 = in.Real(); value.distal_sink_amount_m3 = in.Real();
	value.outward_graph_port_amount_m3 = in.Real(); value.residual_m3 = in.Real(); return value;
}
} // namespace graph_history_detail

inline std::string SerializeGraphHistoryRecord(const GraphAcceptedFlowStep& step)
{
	using namespace graph_history_detail; Writer out; out.Text("IGA_GRAPH_FLOW_HISTORY/1"); Common(out, step, step.iterations);
	Require(step.zero_d_states.size() == step.zero_d_accounting.size(), "0D history catalogs differ"); Count(out, step.zero_d_states.size());
	for (const auto& item : step.zero_d_states) {
		out.Text(item.first); out.Real(item.second.stored_pressure_pa); Accounting(out, step.zero_d_accounting.at(item.first));
	}
	Iterations(out, step.result.iterations); Ports(out, step.result.accepted_ports); return out.Bytes();
}
inline GraphAcceptedFlowStep ParseGraphFlowHistoryRecord(std::string_view bytes)
{
	using namespace graph_history_detail; Reader in(bytes); Require(in.Text() == "IGA_GRAPH_FLOW_HISTORY/1", "unknown flow history schema");
	GraphAcceptedFlowStep step; step.iterations = Common(in, step); std::string previous;
	for (auto n = Count(in); n; --n) {
		auto key = in.Text(); Require(previous < key, "invalid 0D history catalog"); previous = key;
		step.zero_d_states.emplace(key, ZeroDFlowState{in.Real()}); step.zero_d_accounting.emplace(key, Accounting(in));
	}
	step.result.iterations = Iterations(in); step.result.accepted_ports = Ports(in); in.Finish();
	Require(static_cast<std::size_t>(step.iterations) == step.result.iterations.size(), "history iteration count differs"); return step;
}
inline std::string SerializeGraphHistoryRecord(const GraphAcceptedSpeciesStep& step)
{
	using namespace graph_history_detail; Writer out; out.Text("IGA_GRAPH_SPECIES_HISTORY/1"); Common(out, step, step.hydraulic_iterations);
	Ports(out, step.transport_ports); const auto& result = step.result;
	Iterations(out, result.hydraulic_iterations); Ports(out, result.accepted_ports); Ports(out, result.transport_ports);
	Count(out, result.transport_domain_order.size()); for (const auto& domain : result.transport_domain_order) out.Text(domain);
	Count(out, result.donor_ownership.size());
	for (const auto& donor : result.donor_ownership) {
		Require(donor.second == SpeciesDonor::First || donor.second == SpeciesDonor::Second, "unknown history donor");
		out.Text(donor.first.first); out.Text(donor.first.second); out.Unsigned(donor.second == SpeciesDonor::First ? 0 : 1);
	}
	Count(out, result.edge_amounts.size());
	for (const auto& edge : result.edge_amounts) {
		out.Text(edge.edge_id); out.Text(edge.species_id); Reference(out, edge.first); Reference(out, edge.second);
		Require(edge.donor == SpeciesDonor::First || edge.donor == SpeciesDonor::Second, "unknown history edge donor");
		out.Unsigned(edge.donor == SpeciesDonor::First ? 0 : 1);
		for (double value : {edge.first_outward_amount, edge.second_outward_amount, edge.residual, edge.normalized_residual}) out.Real(value);
	}
	Count(out, result.domain_balances.size());
	for (const auto& domain : result.domain_balances) {
		out.Text(domain.domain_id); out.Text(domain.species_id); const auto& a = domain.accounting;
		out.Real(a.initial_mass); out.Real(a.final_mass); out.Real(a.source_amount); out.Reals(a.outward_port_amount); out.Real(a.residual);
		out.Real(domain.recomputed_residual); out.Real(domain.normalized_residual);
	}
	Count(out, result.global_balances.size());
	for (const auto& item : result.global_balances) {
		out.Text(item.first); const auto& a = item.second; out.Text(a.species_id);
		for (double value : {a.initial_mass, a.final_mass, a.source_amount, a.outward_amount, a.gross_activity, a.residual, a.normalized_residual}) out.Real(value);
	}
	return out.Bytes();
}
inline GraphAcceptedSpeciesStep ParseGraphSpeciesHistoryRecord(std::string_view bytes)
{
	using namespace graph_history_detail; Reader in(bytes); Require(in.Text() == "IGA_GRAPH_SPECIES_HISTORY/1", "unknown species history schema");
	GraphAcceptedSpeciesStep step; step.hydraulic_iterations = Common(in, step); step.transport_ports = Ports(in); auto& result = step.result;
	result.hydraulic_iterations = Iterations(in); result.accepted_ports = Ports(in); result.transport_ports = Ports(in);
	for (auto n = Count(in); n; --n) result.transport_domain_order.push_back(in.Text());
	std::pair<std::string, std::string> previous_donor;
	for (auto n = Count(in); n; --n) {
		const auto edge = in.Text(), species = in.Text(); const auto key = std::make_pair(edge, species); const auto donor = in.Unsigned();
		Require(previous_donor < key && donor <= 1, "invalid history donor catalog"); previous_donor = key;
		result.donor_ownership.emplace(key, donor ? SpeciesDonor::Second : SpeciesDonor::First);
	}
	for (auto n = Count(in); n; --n) {
		SpeciesEdgeAmountDiagnostic edge; edge.edge_id = in.Text(); edge.species_id = in.Text(); edge.first = Reference(in); edge.second = Reference(in);
		const auto donor = in.Unsigned(); Require(donor <= 1, "invalid history edge donor"); edge.donor = donor ? SpeciesDonor::Second : SpeciesDonor::First;
		edge.first_outward_amount = in.Real(); edge.second_outward_amount = in.Real(); edge.residual = in.Real(); edge.normalized_residual = in.Real();
		result.edge_amounts.push_back(std::move(edge));
	}
	for (auto n = Count(in); n; --n) {
		SpeciesDomainBalanceDiagnostic domain; domain.domain_id = in.Text(); domain.species_id = in.Text(); auto& a = domain.accounting;
		a.initial_mass = in.Real(); a.final_mass = in.Real(); a.source_amount = in.Real(); a.outward_port_amount = in.Reals(); a.residual = in.Real();
		domain.recomputed_residual = in.Real(); domain.normalized_residual = in.Real(); result.domain_balances.push_back(std::move(domain));
	}
	std::string previous;
	for (auto n = Count(in); n; --n) {
		auto key = in.Text(); Require(previous < key, "invalid history global catalog"); previous = key;
		SpeciesGlobalBalanceDiagnostic a; a.species_id = in.Text(); a.initial_mass = in.Real(); a.final_mass = in.Real(); a.source_amount = in.Real();
		a.outward_amount = in.Real(); a.gross_activity = in.Real(); a.residual = in.Real(); a.normalized_residual = in.Real(); result.global_balances.emplace(key, std::move(a));
	}
	in.Finish(); Require(static_cast<std::size_t>(step.hydraulic_iterations) == result.hydraulic_iterations.size(), "species history iteration count differs"); return step;
}
} // namespace iga
#endif
