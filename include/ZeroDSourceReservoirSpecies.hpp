#ifndef IGA_ZERO_D_SOURCE_RESERVOIR_SPECIES_HPP
#define IGA_ZERO_D_SOURCE_RESERVOIR_SPECIES_HPP

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

struct ZeroDSourceReservoirSpeciesTrial
{
	ZeroDFlowTrial hydraulic;
	ZeroDSpeciesReservoirTrial transport;
};

// The prescribed hydraulic source is an external pump port. Positive pump
// flow enters the compartment; reverse pump flow leaves with mixed species.
inline ZeroDSourceReservoirSpeciesTrial EvaluateZeroDSourceReservoirSpeciesTrial(
	const ZeroDFlowModel& hydraulic_model,
	const ZeroDFlowState& hydraulic_committed,
	const ZeroDSpeciesReservoirModel& species_model,
	const ZeroDSpeciesReservoirState& species_committed,
	double graph_pressure_pa,
	const std::map<std::string,double>& graph_donor_concentration_mol_m3,
	const std::map<std::string,double>& pump_donor_concentration_mol_m3,
	double dt_s,double start_time_s)
{
	if(hydraulic_model.role!=ZeroDFlowRole::SourceReservoir
		||species_model.port_ids!=std::set<std::string>{"graph","pump"})
		throw std::invalid_argument("source reservoir species requires graph and pump ports");
	ZeroDSourceReservoirSpeciesTrial result;
	result.hydraulic=EvaluateZeroDFlowTrial(hydraulic_model,hydraulic_committed,
		graph_pressure_pa,dt_s,start_time_s);
	const double graph_flow=*result.hydraulic.port.outward_flow_m3_s;
	const double pump_flow=-hydraulic_model.source.prescribed_flow_m3_s;
	const std::map<std::string,ZeroDSpeciesReservoirPortInput> inputs{
		{"graph",{graph_flow,graph_donor_concentration_mol_m3}},
		{"pump",{pump_flow,pump_donor_concentration_mol_m3}}};
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
		throw std::runtime_error("source reservoir hydraulic and species volume changes differ");
	return result;
}

} // namespace iga

#endif
