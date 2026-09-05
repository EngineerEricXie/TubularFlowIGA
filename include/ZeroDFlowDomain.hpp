#ifndef IGA_ZERO_D_FLOW_DOMAIN_HPP
#define IGA_ZERO_D_FLOW_DOMAIN_HPP

// Dependency-free 0D hydraulic contracts.  Flow is always positive outward
// from this 0D subsystem, and all quantities use SI units.
#include "CouplingPort.hpp"
#include "Sha256.hpp"

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

enum class ZeroDFlowRole : std::uint8_t {
	SourceReservoir,
	TerminalRcr
};

inline const char* ZeroDFlowRoleName(ZeroDFlowRole role)
{
	if (role == ZeroDFlowRole::SourceReservoir) return "source_reservoir";
	if (role == ZeroDFlowRole::TerminalRcr) return "terminal_rcr";
	return "unknown";
}

inline bool IsKnownZeroDFlowRole(ZeroDFlowRole role)
{
	return role == ZeroDFlowRole::SourceReservoir || role == ZeroDFlowRole::TerminalRcr;
}

struct ZeroDSourceReservoirModel {
	double capacitance_m3_pa = 0.0;
	double resistance_pa_s_m3 = 0.0;
	double prescribed_flow_m3_s = 0.0;
};

struct ZeroDTerminalRcrModel {
	double proximal_resistance_pa_s_m3 = 0.0;
	double distal_resistance_pa_s_m3 = 0.0;
	double capacitance_m3_pa = 0.0;
	double distal_pressure_pa = 0.0;
};

struct ZeroDFlowModel {
	ZeroDFlowRole role = ZeroDFlowRole::SourceReservoir;
	ZeroDSourceReservoirModel source;
	ZeroDTerminalRcrModel terminal;
};

struct ZeroDFlowState {
	double stored_pressure_pa = 0.0;
};

struct ZeroDFlowModelConfiguration {
	ZeroDFlowModel model;
	ZeroDFlowState initial_state;
};

struct ZeroDFlowStorageBalance {
	double stored_volume_change_m3 = 0.0;
	double prescribed_source_amount_m3 = 0.0;
	double distal_sink_amount_m3 = 0.0;
	double outward_graph_port_amount_m3 = 0.0;
	double residual_m3 = 0.0;
};

struct ZeroDFlowTrial {
	ZeroDFlowState state;
	PortState port;
	ZeroDFlowStorageBalance storage;
};

namespace zero_d_flow_detail {

inline void AppendString(Sha256& hash, const std::string& value)
{
	hash.AppendLittleEndian64(static_cast<std::uint64_t>(value.size()));
	hash.Append(value.data(), value.size());
}

inline void RequirePositiveFinite(const std::string& name, double value)
{
	if (!(value > 0.0) || !std::isfinite(value))
		throw std::runtime_error(name+" must be positive and finite");
}

inline void RequireNonnegativeFinite(const std::string& name, double value)
{
	if (!(value >= 0.0) || !std::isfinite(value))
		throw std::runtime_error(name+" must be nonnegative and finite");
}

} // namespace zero_d_flow_detail

inline void ValidateZeroDFlowModel(const ZeroDFlowModel& model)
{
	if (!IsKnownZeroDFlowRole(model.role))
		throw std::runtime_error("0D flow model has an unknown role");
	if (model.role == ZeroDFlowRole::SourceReservoir) {
		zero_d_flow_detail::RequirePositiveFinite("0D source capacitance_m3_pa",
			model.source.capacitance_m3_pa);
		zero_d_flow_detail::RequirePositiveFinite("0D source resistance_pa_s_m3",
			model.source.resistance_pa_s_m3);
		RequireFinitePortValue("0D source prescribed_flow_m3_s",
			model.source.prescribed_flow_m3_s);
	} else {
		zero_d_flow_detail::RequireNonnegativeFinite("0D terminal proximal_resistance_pa_s_m3",
			model.terminal.proximal_resistance_pa_s_m3);
		zero_d_flow_detail::RequirePositiveFinite("0D terminal distal_resistance_pa_s_m3",
			model.terminal.distal_resistance_pa_s_m3);
		zero_d_flow_detail::RequirePositiveFinite("0D terminal capacitance_m3_pa",
			model.terminal.capacitance_m3_pa);
		RequireFinitePortValue("0D terminal distal_pressure_pa", model.terminal.distal_pressure_pa);
	}
}

inline void ValidateZeroDFlowState(const ZeroDFlowState& state)
{
	RequireFinitePortValue("0D flow stored_pressure_pa", state.stored_pressure_pa);
}

inline std::string BuildZeroDFlowModelIdentitySha256(const ZeroDFlowModel& model)
{
	ValidateZeroDFlowModel(model);
	Sha256 hash;
	zero_d_flow_detail::AppendString(hash, "ZeroDFlowModel/v1");
	hash.AppendLittleEndian32(static_cast<std::uint32_t>(model.role));
	if (model.role == ZeroDFlowRole::SourceReservoir) {
		hash.AppendNormalizedDouble(model.source.capacitance_m3_pa);
		hash.AppendNormalizedDouble(model.source.resistance_pa_s_m3);
		hash.AppendNormalizedDouble(model.source.prescribed_flow_m3_s);
	} else {
		hash.AppendNormalizedDouble(model.terminal.proximal_resistance_pa_s_m3);
		hash.AppendNormalizedDouble(model.terminal.distal_resistance_pa_s_m3);
		hash.AppendNormalizedDouble(model.terminal.capacitance_m3_pa);
		hash.AppendNormalizedDouble(model.terminal.distal_pressure_pa);
	}
	return hash.Hex();
}

inline std::string BuildZeroDFlowStateIdentitySha256(const ZeroDFlowModel& model,
	const ZeroDFlowState& state)
{
	ValidateZeroDFlowModel(model);
	ValidateZeroDFlowState(state);
	Sha256 hash;
	zero_d_flow_detail::AppendString(hash, "ZeroDFlowState/v1");
	zero_d_flow_detail::AppendString(hash, BuildZeroDFlowModelIdentitySha256(model));
	hash.AppendNormalizedDouble(state.stored_pressure_pa);
	return hash.Hex();
}

inline CouplingPort MakeZeroDFlowPort(const std::string& domain_id, ZeroDFlowRole role)
{
	if (domain_id.empty()) throw std::runtime_error("0D flow domain id must be nonempty");
	if (!IsKnownZeroDFlowRole(role)) throw std::runtime_error("0D flow role is unknown");
	CouplingPort port;
	port.id = "port";
	port.subsystem_id = domain_id;
	port.locator_kind = "zero_d_port";
	port.locator = "port";
	port.provides = {PortQuantity::MeanPressure, PortQuantity::FlowRate};
	port.requires = role == ZeroDFlowRole::SourceReservoir
		? std::set<PortQuantity>{PortQuantity::MeanPressure}
		: std::set<PortQuantity>{PortQuantity::FlowRate};
	return port;
}

inline void ValidateZeroDFlowDomainMetadata(const std::string& domain_id,
	const std::vector<CouplingPort>& ports, ZeroDFlowRole role)
{
	if (ports.size() != 1) throw std::runtime_error("0D flow domain requires exactly one port");
	const auto expected = MakeZeroDFlowPort(domain_id, role);
	const auto& port = ports.front();
	ValidateCouplingPort(port);
	if (port.id != expected.id || port.subsystem_id != expected.subsystem_id
		|| port.locator_kind != expected.locator_kind || port.locator != expected.locator
		|| port.orientation.native_to_outward_sign != 1 || port.provides != expected.provides
		|| port.requires != expected.requires || !port.species.empty())
		throw std::runtime_error("0D flow domain port does not match its model role");
}

inline ZeroDFlowTrial EvaluateZeroDFlowTrial(const ZeroDFlowModel& model,
	const ZeroDFlowState& committed, double port_input, double dt_s, double time_s = 0.0)
{
	ValidateZeroDFlowModel(model);
	ValidateZeroDFlowState(committed);
	zero_d_flow_detail::RequirePositiveFinite("0D flow dt_s", dt_s);
	RequireFinitePortValue("0D flow port input", port_input);
	RequireFinitePortValue("0D flow time_s", time_s);
	ZeroDFlowTrial result;
	result.port.time_s = time_s+dt_s;
	if (model.role == ZeroDFlowRole::SourceReservoir) {
		const auto& source = model.source;
		const double pressure = ((source.capacitance_m3_pa/dt_s)*committed.stored_pressure_pa
			+source.prescribed_flow_m3_s+port_input/source.resistance_pa_s_m3)
			/((source.capacitance_m3_pa/dt_s)+1.0/source.resistance_pa_s_m3);
		const double outward_flow = (pressure-port_input)/source.resistance_pa_s_m3;
		result.state.stored_pressure_pa = pressure;
		// The source consumes the interface pressure; it is reported unchanged
		// because this one-port contract has no fabricated internal boundary area.
		result.port.mean_pressure_pa = port_input;
		result.port.outward_flow_m3_s = outward_flow;
		result.storage.stored_volume_change_m3 = source.capacitance_m3_pa
			*(pressure-committed.stored_pressure_pa);
		result.storage.prescribed_source_amount_m3 = source.prescribed_flow_m3_s*dt_s;
		result.storage.outward_graph_port_amount_m3 = outward_flow*dt_s;
		result.storage.residual_m3 = result.storage.stored_volume_change_m3
			-result.storage.prescribed_source_amount_m3
			+result.storage.outward_graph_port_amount_m3;
	} else {
		const auto& terminal = model.terminal;
		// Input uses the 0D port convention: Q_port=-Qin.
		const double inlet_flow = -port_input;
		const double pressure = ((terminal.capacitance_m3_pa/dt_s)
			*committed.stored_pressure_pa+inlet_flow
			+terminal.distal_pressure_pa/terminal.distal_resistance_pa_s_m3)
			/((terminal.capacitance_m3_pa/dt_s)+1.0/terminal.distal_resistance_pa_s_m3);
		const double distal_flow = (pressure-terminal.distal_pressure_pa)
			/terminal.distal_resistance_pa_s_m3;
		result.state.stored_pressure_pa = pressure;
		result.port.mean_pressure_pa = pressure
			+terminal.proximal_resistance_pa_s_m3*inlet_flow;
		result.port.outward_flow_m3_s = port_input;
		result.storage.stored_volume_change_m3 = terminal.capacitance_m3_pa
			*(pressure-committed.stored_pressure_pa);
		result.storage.distal_sink_amount_m3 = distal_flow*dt_s;
		result.storage.outward_graph_port_amount_m3 = port_input*dt_s;
		result.storage.residual_m3 = result.storage.stored_volume_change_m3
			+result.storage.distal_sink_amount_m3
			+result.storage.outward_graph_port_amount_m3;
	}
	ValidateZeroDFlowState(result.state);
	ValidatePortState(result.port);
	return result;
}

} // namespace iga

#endif
