#ifndef IGA_ZERO_D_TERMINAL_RCR_SPECIES_HPP
#define IGA_ZERO_D_TERMINAL_RCR_SPECIES_HPP

#include "ZeroDFlowDomain.hpp"
#include "ZeroDSpeciesReservoir.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>

namespace iga {

struct ZeroDTerminalRcrSpeciesTrial
{
	ZeroDFlowTrial hydraulic;
	ZeroDSpeciesReservoirTrial transport;
};

// The RCR capacitor's pressure-volume law determines only delta-V. Its
// absolute mixing volume is supplied independently by the transport state.
inline ZeroDTerminalRcrSpeciesTrial EvaluateZeroDTerminalRcrSpeciesTrial(
	const ZeroDFlowModel& hydraulic_model,
	const ZeroDFlowState& hydraulic_committed,
	const ZeroDSpeciesReservoirModel& species_model,
	const ZeroDSpeciesReservoirState& species_committed,
	double graph_outward_flow_m3_s,
	const std::map<std::string,double>& graph_donor_concentration_mol_m3,
	const std::map<std::string,double>& distal_donor_concentration_mol_m3,
	double dt_s,double start_time_s)
{
	if(hydraulic_model.role!=ZeroDFlowRole::TerminalRcr
		||species_model.port_ids!=std::set<std::string>{"graph","distal"})
		throw std::invalid_argument("terminal RCR species requires graph and distal ports");
	ZeroDTerminalRcrSpeciesTrial result;
	result.hydraulic=EvaluateZeroDFlowTrial(hydraulic_model,hydraulic_committed,
		graph_outward_flow_m3_s,dt_s,start_time_s);
	const double distal_flow=result.hydraulic.storage.distal_sink_amount_m3/dt_s;
	const std::map<std::string,ZeroDSpeciesReservoirPortInput> inputs{
		{"graph",{graph_outward_flow_m3_s,graph_donor_concentration_mol_m3}},
		{"distal",{distal_flow,distal_donor_concentration_mol_m3}}};
	result.transport=EvaluateZeroDSpeciesReservoirTrial(species_model,
		species_committed,inputs,dt_s,start_time_s+dt_s);
	const double transport_change=result.transport.state.volume_m3
		-species_committed.volume_m3;
	const double hydraulic_change=result.hydraulic.storage.stored_volume_change_m3;
	const double scale=std::max({std::abs(species_committed.volume_m3),
		std::abs(result.transport.state.volume_m3),std::abs(hydraulic_change)});
	const double tolerance=64.0*std::numeric_limits<double>::epsilon()*scale;
	if(!std::isfinite(transport_change)||!std::isfinite(hydraulic_change)
		||std::abs(transport_change-hydraulic_change)>tolerance)
		throw std::runtime_error("terminal RCR and species stored volume changes differ");
	return result;
}

} // namespace iga

#endif
