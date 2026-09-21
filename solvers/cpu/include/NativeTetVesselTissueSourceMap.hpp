#ifndef IGA_NATIVE_TET_VESSEL_TISSUE_SOURCE_MAP_HPP
#define IGA_NATIVE_TET_VESSEL_TISSUE_SOURCE_MAP_HPP

#include "NativeTetFem.hpp"
#include "NativeTetBoundaryFlow.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

// Positive vessel outward flow [m^3/s] enters the tissue as a Darcy source.
// Negative flow reverses the direction. Each port explicitly names its target
// tissue tetrahedra and nonnegative weights; no anatomical mapping is inferred.
struct NativeTetVesselTissuePort
{
	std::string name;
	double vessel_outward_flow_m3_s=0.;
	std::vector<std::pair<std::uint64_t,double>> tissue_cell_id_weights;
};

struct NativeTetVesselTissueSourceResult
{
	std::vector<double> tissue_source_s_inv;
	std::map<std::string,double> vessel_outward_flow_m3_s;
	std::map<std::string,double> tissue_inward_flow_m3_s;
	double total_vessel_outward_flow_m3_s=0.;
	double total_tissue_inward_flow_m3_s=0.;
	double maximum_port_balance_defect_m3_s=0.;
};

// The label must identify an explicitly declared exchange surface in the
// vessel mesh. It is never inferred from a name, organ, or proximity.
struct NativeTetVesselTissueBoundaryPort
{
	std::string name;
	int vessel_boundary_label=0;
	std::vector<std::pair<std::uint64_t,double>> tissue_cell_id_weights;
};

inline NativeTetVesselTissueSourceResult MapNativeTetVesselFlowToTissueSource(
	const NativeTetMesh& tissue_mesh,
	const std::vector<NativeTetVesselTissuePort>& ports)
{
	if(tissue_mesh.cells.empty()||ports.empty())
		throw std::invalid_argument("vessel-tissue source map requires cells and ports");
	std::map<std::uint64_t,std::size_t> cell_index;
	std::vector<double> volumes(tissue_mesh.cells.size());
	for(std::size_t index=0;index<tissue_mesh.cells.size();++index){
		if(!cell_index.emplace(tissue_mesh.cells[index].id,index).second)
			throw std::invalid_argument("vessel-tissue source map has duplicate cell id");
		volumes[index]=EvaluateNativeTetGeometry(tissue_mesh,
			tissue_mesh.cells[index]).determinant/6.;
	}
	NativeTetVesselTissueSourceResult result;
	result.tissue_source_s_inv.assign(tissue_mesh.cells.size(),0.);
	for(const auto& port:ports){
		if(port.name.empty()||!std::isfinite(port.vessel_outward_flow_m3_s)
			||port.tissue_cell_id_weights.empty()
			||!result.vessel_outward_flow_m3_s.emplace(port.name,
				port.vessel_outward_flow_m3_s).second)
			throw std::invalid_argument("vessel-tissue source port is invalid");
		std::set<std::uint64_t> used_cells;
		double weight_sum=0.,distributed=0.;
		for(const auto& target:port.tissue_cell_id_weights){
			const auto found=cell_index.find(target.first);
			if(found==cell_index.end()||!used_cells.insert(target.first).second
				||!std::isfinite(target.second)||!(target.second>0.))
				throw std::invalid_argument("vessel-tissue source target is invalid");
			weight_sum+=target.second;
			const double flow=port.vessel_outward_flow_m3_s*target.second;
			result.tissue_source_s_inv[found->second]+=flow/volumes[found->second];
			distributed+=flow;
		}
		if(!std::isfinite(weight_sum)||std::abs(weight_sum-1.)>1e-12)
			throw std::invalid_argument("vessel-tissue source weights must sum to one");
		result.tissue_inward_flow_m3_s.emplace(port.name,distributed);
		result.total_vessel_outward_flow_m3_s+=port.vessel_outward_flow_m3_s;
		result.total_tissue_inward_flow_m3_s+=distributed;
		result.maximum_port_balance_defect_m3_s=std::max(
			result.maximum_port_balance_defect_m3_s,
			std::abs(port.vessel_outward_flow_m3_s-distributed));
	}
	for(const double source:result.tissue_source_s_inv)
		if(!std::isfinite(source))
			throw std::invalid_argument("vessel-tissue source is nonfinite");
	return result;
}

// Extract each signed Q [m^3/s] from the project-owned P2/P1 vessel field,
// then use the same conservative cell-source map as externally supplied Q.
// The caller remains responsible for proving the labelled vessel surface is
// geometrically an exchange interface with the named tissue cells.
inline NativeTetVesselTissueSourceResult MapNativeTetVesselStateToTissueSource(
	const NativeTetMesh& vessel_mesh,const NativeTaylorHoodTopology& topology,
	const std::vector<double>& vessel_state,const NativeTetMesh& tissue_mesh,
	const std::vector<NativeTetVesselTissueBoundaryPort>& boundaries)
{
	if(boundaries.empty())
		throw std::invalid_argument("vessel-tissue boundary map requires ports");
	for(const double value:vessel_state)
		if(!std::isfinite(value))
			throw std::invalid_argument("vessel-tissue boundary state is nonfinite");
	const auto measured=EvaluateNativeTetBoundaryFlows(vessel_mesh,topology,vessel_state);
	std::set<int> used_labels;
	std::vector<NativeTetVesselTissuePort> ports;
	ports.reserve(boundaries.size());
	for(const auto& boundary:boundaries){
		const auto found=measured.find(boundary.vessel_boundary_label);
		if(found==measured.end()||!used_labels.insert(boundary.vessel_boundary_label).second
			||!(found->second.area_m2>0.)
			||!std::isfinite(found->second.outward_flow_m3_s))
			throw std::invalid_argument("vessel-tissue exchange boundary label is invalid");
		ports.push_back({boundary.name,found->second.outward_flow_m3_s,
			boundary.tissue_cell_id_weights});
	}
	return MapNativeTetVesselFlowToTissueSource(tissue_mesh,ports);
}

} // namespace iga

#endif
