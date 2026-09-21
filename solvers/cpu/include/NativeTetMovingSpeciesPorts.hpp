#ifndef IGA_NATIVE_TET_MOVING_SPECIES_PORTS_HPP
#define IGA_NATIVE_TET_MOVING_SPECIES_PORTS_HPP

#include "NativeTetMovingSpeciesTransport.hpp"
#include "CoupledDomainRuntime.hpp"
#include "FlowDomainPortMetadata.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

// Maps one logical graph species to native tetra boundary labels. This is a
// scalar port contract, not a CoupledDomainRuntime implementation.
class NativeTetMovingSpeciesPorts
{
public:
	NativeTetMovingSpeciesPorts(std::string domain_id,std::string species_id,
		const NativeTetMesh& mesh,std::vector<CouplingPort> ports)
		:domain_id_(std::move(domain_id)),species_id_(std::move(species_id)),
		ports_(std::move(ports))
	{
		if(domain_id_.empty()||species_id_.empty())
			throw std::invalid_argument("native moving species port identity is empty");
		ValidateThreeDBodyFittedFlowDomainMetadata(domain_id_,ports_);
		std::set<int> mesh_labels;
		for(const auto& face:mesh.boundary_triangles)
			mesh_labels.insert(face.boundary_label);
		for(const auto& port:ports_){
			if(port.subsystem_id!=domain_id_
				||port.locator_kind!="boundary_label"
				||port.orientation.native_to_outward_sign!=1
				||port.species!=std::set<std::string>{species_id_}
				||!port.provides.count(PortQuantity::FlowRate)
				||!port.provides.count(PortQuantity::MeanPressure)
				||!port.provides.count(PortQuantity::SpeciesConcentration)
				||!port.provides.count(PortQuantity::SpeciesFlux)
				||!port.requires.count(PortQuantity::SpeciesConcentration)
				||!port.requires.count(PortQuantity::SpeciesFlux)
				||(port.requires.count(PortQuantity::FlowRate)
					+port.requires.count(PortQuantity::MeanPressure)!=1))
				throw std::invalid_argument("native moving species port metadata is unsupported");
			const int label=ParseThreeDFlowBoundaryLabel(port);
			if(!mesh_labels.count(label)||!label_to_port_.emplace(label,port.id).second)
				throw std::invalid_argument("native moving species port label is absent or duplicate");
		}
		if(label_to_port_.empty())
			throw std::invalid_argument("native moving species port catalog is empty");
	}

	std::map<int,double> InflowByLabel(
		const std::map<std::string,double>& outward_flow_m3_s,
		const std::map<std::string,PortBoundaryData>& inputs,
		double time_s,double flow_epsilon_m3_s=0.) const
	{
		if(!std::isfinite(time_s)||!(flow_epsilon_m3_s>=0.)
			||!std::isfinite(flow_epsilon_m3_s))
			throw std::invalid_argument("native moving species port control is invalid");
		std::map<int,double> result;
		for(const auto& entry:inputs)
			if(!HasPort(entry.first))
				throw std::invalid_argument("native moving species input port is unknown");
		for(const auto& entry:outward_flow_m3_s)
			if(!HasPort(entry.first))
				throw std::invalid_argument("native moving species flow port is unknown");
		for(const auto& port:ports_){
			const auto flow=outward_flow_m3_s.find(port.id);
			if(flow==outward_flow_m3_s.end()||!std::isfinite(flow->second))
				throw std::invalid_argument("native moving species port flow is missing");
			const auto input=inputs.find(port.id);
			bool has_concentration=false;
			if(input!=inputs.end()){
				ValidatePortBoundaryData(input->second);
				if(input->second.outward_flow_m3_s
					||input->second.mean_pressure_pa
					||input->second.mean_normal_traction_pa
					||input->second.total_pressure_pa
					||!input->second.outward_species_flux.empty())
					throw std::invalid_argument("native moving species scalar input contains unsupported quantities");
				if(std::abs(input->second.time_s-time_s)
					>1e-12*std::max({1.,std::abs(time_s),
						std::abs(input->second.time_s)}))
					throw std::invalid_argument("native moving species port input time differs");
				for(const auto& species:input->second.concentration)
					if(species.first!=species_id_)
						throw std::invalid_argument("native moving species input contains another species");
				has_concentration=input->second.concentration.count(species_id_)!=0;
			}
			if(flow->second < -flow_epsilon_m3_s){
				if(!has_concentration||input->second.concentration.at(species_id_)<0.)
					throw std::invalid_argument("native moving species inward port requires nonnegative concentration");
				result.emplace(ParseThreeDFlowBoundaryLabel(port),
					input->second.concentration.at(species_id_));
			}else if(flow->second>flow_epsilon_m3_s&&has_concentration)
				throw std::invalid_argument("native moving species non-inward port cannot impose concentration");
		}
		return result;
	}

	std::map<std::string,PortState> Observe(const NativeTetMesh& current_mesh,
		const NativeTetMovingSpeciesStep& step,
		double time_s,double flow_epsilon_m3_s=0.) const
	{
		if(!std::isfinite(time_s)||!(flow_epsilon_m3_s>=0.)
			||!std::isfinite(flow_epsilon_m3_s))
			throw std::invalid_argument("native moving species port observation control is invalid");
		if(step.concentration_mol_m3.size()!=current_mesh.points.size())
			throw std::invalid_argument("native moving species port concentration shape is invalid");
		for(const double value:step.concentration_mol_m3)
			if(!std::isfinite(value)||value<0.)
				throw std::runtime_error("native moving species port concentration is invalid");
		std::map<std::string,PortState> result;
		for(const auto& entry:label_to_port_){
			const double flow=step.outward_relative_flow_by_label_m3_s.at(entry.first);
			const double flux=step.outward_advective_flux_by_label_mol_s.at(entry.first);
			const double positive=step.positive_relative_flow_by_label_m3_s.at(entry.first);
			const double negative=step.negative_relative_flow_by_label_m3_s.at(entry.first);
			if(!std::isfinite(flow)||!std::isfinite(flux))
				throw std::runtime_error("native moving species port observation is nonfinite");
			if(!std::isfinite(positive)||!std::isfinite(negative)
				||positive < 0.||negative > 0.
				||(positive>flow_epsilon_m3_s
					&&negative<-flow_epsilon_m3_s))
				throw std::runtime_error("native moving species port has mixed or invalid flow directions");
			PortState state;state.time_s=time_s;state.outward_flow_m3_s=flow;
			state.outward_species_flux.emplace(species_id_,flux);
			if(std::abs(flow)>flow_epsilon_m3_s){
				const double concentration=flux/flow;
				if(!std::isfinite(concentration)||concentration<0.)
					throw std::runtime_error("native moving species port effective concentration is invalid");
				state.concentration.emplace(species_id_,concentration);
			}else{
				if(std::abs(flux)>1e-12)
					throw std::runtime_error("native moving species near-zero flow has nonzero flux");
				state.concentration.emplace(species_id_,
					AreaMeanConcentration(current_mesh,step.concentration_mol_m3,
						entry.first));
			}
			ValidatePortState(state);result.emplace(entry.second,std::move(state));
		}
		return result;
	}

	SpeciesStepAccounting Accounting(const NativeTetMovingSpeciesStep& step,
		double dt_s,double unported_flux_tolerance_mol_s=1e-12) const
	{
		if(!(dt_s>0.)||!std::isfinite(dt_s)
			||!(unported_flux_tolerance_mol_s>=0.)
			||!std::isfinite(unported_flux_tolerance_mol_s))
			throw std::invalid_argument("native moving species accounting control is invalid");
		SpeciesStepAccounting result;
		result.initial_mass=step.previous_inventory_mol;
		result.final_mass=step.current_inventory_mol;
		result.source_amount=step.source_mol_s*dt_s;
		for(const auto& flux:step.outward_advective_flux_by_label_mol_s){
			const auto port=label_to_port_.find(flux.first);
			if(port==label_to_port_.end()){
				if(std::abs(flux.second)>unported_flux_tolerance_mol_s)
					throw std::runtime_error("native moving species unported boundary has flux");
				continue;
			}
			result.outward_port_amount.emplace(port->second,flux.second*dt_s);
		}
		if(result.outward_port_amount.size()!=label_to_port_.size())
			throw std::runtime_error("native moving species port accounting is incomplete");
		result.residual=result.final_mass-result.initial_mass-result.source_amount;
		for(const auto& amount:result.outward_port_amount)
			result.residual+=amount.second;
		return result;
	}

private:
	static double AreaMeanConcentration(const NativeTetMesh& mesh,
		const std::vector<double>& concentration,int label)
	{
		if(concentration.size()!=mesh.points.size())
			throw std::invalid_argument("native moving species port concentration size is invalid");
		double area=0.,integral=0.;
		for(const auto& face:mesh.boundary_triangles){
			if(face.boundary_label!=label)continue;
			for(const auto node:face.nodes)
				if(node>=mesh.points.size()||!std::isfinite(concentration[node])
					||concentration[node]<0.)
					throw std::invalid_argument("native moving species port face concentration is invalid");
			const auto& a=mesh.points[face.nodes[0]];
			const auto& b=mesh.points[face.nodes[1]];
			const auto& c=mesh.points[face.nodes[2]];
			const std::array<double,3> u{{b[0]-a[0],b[1]-a[1],b[2]-a[2]}};
			const std::array<double,3> v{{c[0]-a[0],c[1]-a[1],c[2]-a[2]}};
			const std::array<double,3> cross{{u[1]*v[2]-u[2]*v[1],
				u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]}};
			const double patch=0.5*std::sqrt(cross[0]*cross[0]
				+cross[1]*cross[1]+cross[2]*cross[2]);
			if(!(patch>0.)||!std::isfinite(patch))
				throw std::invalid_argument("native moving species port face area is invalid");
			area+=patch;
			integral+=patch*(concentration[face.nodes[0]]
				+concentration[face.nodes[1]]+concentration[face.nodes[2]])/3.;
		}
		if(!(area>0.)||!std::isfinite(integral))
			throw std::invalid_argument("native moving species port label has no valid area");
		return integral/area;
	}

	bool HasPort(const std::string& id) const
	{
		for(const auto& port:ports_)if(port.id==id)return true;
		return false;
	}
	std::string domain_id_,species_id_;
	std::vector<CouplingPort> ports_;
	std::map<int,std::string> label_to_port_;
};

} // namespace iga

#endif
