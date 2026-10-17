#ifndef IGA_PERFUSED_ORGAN_TRANSPORT_HPP
#define IGA_PERFUSED_ORGAN_TRANSPORT_HPP

#include "NativeTetFem.hpp"
#include "PerfusedOrganKinetics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace iga {

// Cell graph for conservative upwind transport. Every quantity is a volume
// flow [m^3/s]; `upstream[i]` lists (neighbour j, flow j -> i). The exchange
// inflow carries the perfusate inlet concentration; exchange outflow and face
// outflow leave at the cell concentration. The graph is mesh-agnostic, so the
// same solver serves tetrahedra or aggregated compartments.
struct OrganTransportGraph
{
	std::vector<double> volume_m3;
	std::vector<double> face_outflow_m3_s;     // to neighbours and the boundary
	std::vector<double> boundary_outflow_m3_s; // the boundary part of face_outflow
	std::vector<double> exchange_in_m3_s,exchange_out_m3_s;
	std::vector<std::vector<std::pair<std::size_t,double>>> upstream;
};

// Builds the graph from conservative tetra face flows (outward from each cell,
// on the face opposite each local node). Flows below `relative_cutoff` times
// the largest face flow are dropped, which removes round-off cycles.
inline OrganTransportGraph BuildTetTransportGraph(const NativeTetMesh& mesh,
	const std::vector<std::array<double,4>>& face_flow_m3_s,
	const std::vector<double>& exchange_in_m3_s,const std::vector<double>& exchange_out_m3_s,
	double relative_cutoff=1e-12)
{
	const std::size_t cells=mesh.cells.size();
	if(face_flow_m3_s.size()!=cells||exchange_in_m3_s.size()!=cells||exchange_out_m3_s.size()!=cells)
		throw std::invalid_argument("organ transport graph inputs do not match the cells");
	OrganTransportGraph graph;
	graph.volume_m3.resize(cells);
	graph.face_outflow_m3_s.assign(cells,0.);
	graph.boundary_outflow_m3_s.assign(cells,0.);
	graph.exchange_in_m3_s=exchange_in_m3_s;
	graph.exchange_out_m3_s=exchange_out_m3_s;
	graph.upstream.assign(cells,{});
	double largest=0.;
	for(const auto& flows:face_flow_m3_s)
		for(const double flow:flows)largest=std::max(largest,std::abs(flow));
	const double cutoff=relative_cutoff*largest;
	std::map<std::array<std::uint32_t,3>,std::pair<std::size_t,std::size_t>> faces;
	for(std::size_t i=0;i<cells;++i){
		graph.volume_m3[i]=EvaluateNativeTetGeometry(mesh,mesh.cells[i]).determinant/6.;
		for(std::size_t opposite=0;opposite<4;++opposite){
			std::array<std::uint32_t,3> face{};
			std::size_t next=0;
			for(std::size_t local=0;local<4;++local)
				if(local!=opposite)face[next++]=mesh.cells[i].nodes[local];
			std::sort(face.begin(),face.end());
			const auto found=faces.find(face);
			if(found==faces.end()){faces.emplace(face,std::make_pair(i,opposite));continue;}
			const std::size_t j=found->second.first;
			const double flow=face_flow_m3_s[i][opposite]; // i -> j positive
			faces.erase(found);
			if(std::abs(flow)<=cutoff)continue;
			if(flow>0.){graph.face_outflow_m3_s[i]+=flow;graph.upstream[j].push_back({i,flow});}
			else{graph.face_outflow_m3_s[j]-=flow;graph.upstream[i].push_back({j,-flow});}
		}
	}
	// Remaining faces are on the boundary; only outflow is admitted there.
	for(const auto& item:faces){
		const double flow=face_flow_m3_s[item.second.first][item.second.second];
		if(flow>cutoff){
			graph.face_outflow_m3_s[item.second.first]+=flow;
			graph.boundary_outflow_m3_s[item.second.first]+=flow;
		}else if(flow<-cutoff)
			throw std::invalid_argument("organ transport boundary inflow is not supported");
	}
	return graph;
}

// Splits every cell into `zones` well-mixed lobular zones in series (tanks in
// series for the sinusoidal path; zones = 3 mirrors hepatic zones 1-3). Node
// i*zones+z is zone z of cell i with volume V_i/zones. Exchange inflow and
// inflow from neighbouring cells enter zone 0; the cell's throughflow passes
// zone by zone; exchange outflow and outflow to neighbours (into their zone 0)
// leave from the last zone. zones = 1 returns the cell graph unchanged.
inline OrganTransportGraph SplitIntoLobularZones(const OrganTransportGraph& cells,std::size_t zones)
{
	if(zones==0)throw std::invalid_argument("lobular zone count must be positive");
	if(zones==1)return cells;
	const std::size_t n=cells.volume_m3.size();
	OrganTransportGraph graph;
	graph.volume_m3.resize(n*zones);
	graph.face_outflow_m3_s.assign(n*zones,0.);
	graph.boundary_outflow_m3_s.assign(n*zones,0.);
	graph.exchange_in_m3_s.assign(n*zones,0.);
	graph.exchange_out_m3_s.assign(n*zones,0.);
	graph.upstream.assign(n*zones,{});
	for(std::size_t i=0;i<n;++i){
		double through=cells.exchange_in_m3_s[i];
		for(const auto& edge:cells.upstream[i])through+=edge.second;
		const std::size_t first=i*zones,last=first+zones-1;
		for(std::size_t z=0;z<zones;++z){
			graph.volume_m3[first+z]=cells.volume_m3[i]/static_cast<double>(zones);
			if(z>0)graph.upstream[first+z].push_back({first+z-1,through});
			if(z+1<zones)graph.face_outflow_m3_s[first+z]=through;
		}
		graph.exchange_in_m3_s[first]=cells.exchange_in_m3_s[i];
		for(const auto& edge:cells.upstream[i])
			graph.upstream[first].push_back({edge.first*zones+zones-1,edge.second});
		graph.exchange_out_m3_s[last]=cells.exchange_out_m3_s[i];
		graph.face_outflow_m3_s[last]=cells.face_outflow_m3_s[i];
		graph.boundary_outflow_m3_s[last]=cells.boundary_outflow_m3_s[i];
	}
	return graph;
}

struct OrganTransportInput
{
	double dt_s=0.;
	OrganConcentration inlet{};
	OrganConcentration distribution_volume_fraction{};
	std::vector<OrganConcentration> previous; // per cell
	std::vector<HepaticCellState> state;      // per cell
	HepaticKinetics kinetics;
	PerfusateChemistry chemistry;
};

struct OrganTransportResult
{
	std::vector<OrganConcentration> concentration;
	std::vector<HepaticRates> rates;
	// Per-species totals over the organ [mol/s or mU L^-1 m^3 s^-1].
	OrganConcentration inflow{},outflow{},storage_rate{},reaction{},balance_defect{};
	std::size_t components=0,largest_component=0,newton_failures=0;
	std::size_t gauss_seidel_sweeps=0;
};

namespace organ_transport_detail {

// Solves the dense n x n system a x = b in place by partial pivoting.
template<std::size_t N>
inline bool SolveDense(std::array<std::array<double,N>,N> a,std::array<double,N>& b)
{
	for(std::size_t column=0;column<N;++column){
		std::size_t pivot=column;
		for(std::size_t row=column+1;row<N;++row)
			if(std::abs(a[row][column])>std::abs(a[pivot][column]))pivot=row;
		if(!(std::abs(a[pivot][column])>0.))return false;
		std::swap(a[pivot],a[column]);
		std::swap(b[pivot],b[column]);
		for(std::size_t row=column+1;row<N;++row){
			const double factor=a[row][column]/a[column][column];
			for(std::size_t k=column;k<N;++k)a[row][k]-=factor*a[column][k];
			b[row]-=factor*b[column];
		}
	}
	for(std::size_t column=N;column-->0;){
		for(std::size_t k=column+1;k<N;++k)b[column]-=a[column][k]*b[k];
		b[column]/=a[column][column];
	}
	return true;
}

// Strongly connected components of the upstream graph in upstream-first order
// (iterative Tarjan; components are emitted downstream-first and reversed).
inline std::vector<std::vector<std::size_t>> UpstreamFirstComponents(const OrganTransportGraph& graph)
{
	const std::size_t n=graph.volume_m3.size();
	std::vector<std::vector<std::size_t>> downstream(n);
	for(std::size_t i=0;i<n;++i)
		for(const auto& edge:graph.upstream[i])downstream[edge.first].push_back(i);
	const std::size_t unvisited=static_cast<std::size_t>(-1);
	std::vector<std::size_t> index(n,unvisited),low(n,0);
	std::vector<char> on_stack(n,0);
	std::vector<std::size_t> stack;
	std::vector<std::vector<std::size_t>> components;
	std::size_t counter=0;
	std::vector<std::pair<std::size_t,std::size_t>> call;
	for(std::size_t start=0;start<n;++start){
		if(index[start]!=unvisited)continue;
		call.push_back({start,0});
		index[start]=low[start]=counter++;
		stack.push_back(start);on_stack[start]=1;
		while(!call.empty()){
			auto& frame=call.back();
			const std::size_t v=frame.first;
			if(frame.second<downstream[v].size()){
				const std::size_t w=downstream[v][frame.second++];
				if(index[w]==unvisited){
					index[w]=low[w]=counter++;
					stack.push_back(w);on_stack[w]=1;
					call.push_back({w,0});
				}else if(on_stack[w])low[v]=std::min(low[v],index[w]);
				continue;
			}
			if(low[v]==index[v]){
				std::vector<std::size_t> component;
				std::size_t w=0;
				do{w=stack.back();stack.pop_back();on_stack[w]=0;component.push_back(w);}while(w!=v);
				components.push_back(component);
			}
			call.pop_back();
			if(!call.empty())low[call.back().first]=std::min(low[call.back().first],low[v]);
		}
	}
	std::reverse(components.begin(),components.end());
	return components;
}

} // namespace organ_transport_detail

// Backward Euler over dt for every species in every cell:
//   phi_s V (c - c_old)/dt + (F_out + X_out) c - sum_j F_ji c_j - X_in c_inlet = V R(c)
// With perfusion residence times far below dt this is the quasi-steady state;
// unperfused cells relax over dt instead of having no steady state.
inline OrganTransportResult SolveOrganTransport(const OrganTransportGraph& graph,
	const OrganTransportInput& input)
{
	using namespace organ_transport_detail;
	const std::size_t n=graph.volume_m3.size();
	if(!(input.dt_s>0.)||!std::isfinite(input.dt_s))
		throw std::invalid_argument("organ transport time step must be positive");
	if(input.previous.size()!=n||input.state.size()!=n||graph.upstream.size()!=n
		||graph.face_outflow_m3_s.size()!=n||graph.boundary_outflow_m3_s.size()!=n
		||graph.exchange_in_m3_s.size()!=n
		||graph.exchange_out_m3_s.size()!=n)
		throw std::invalid_argument("organ transport inputs do not match the graph");
	for(std::size_t s=0;s<OrganSpeciesCount;++s)
		if(!(input.distribution_volume_fraction[s]>0.)||!(input.inlet[s]>=0.))
			throw std::invalid_argument("organ transport species data are invalid");
	ValidateHepaticKinetics(input.kinetics,input.chemistry);
	OrganTransportResult result;
	result.concentration=input.previous;
	result.rates.assign(n,{});
	OrganConcentration scale{};
	for(std::size_t s=0;s<OrganSpeciesCount;++s){
		scale[s]=std::max(input.inlet[s],1e-6);
		for(const auto& value:input.previous)scale[s]=std::max(scale[s],std::abs(value[s]));
	}

	const auto residual=[&](std::size_t i,const OrganConcentration& c,OrganConcentration& r){
		const double volume=graph.volume_m3[i];
		const double leaving=graph.face_outflow_m3_s[i]+graph.exchange_out_m3_s[i];
		const auto rates=EvaluateHepaticRates(c,input.state[i],input.kinetics,input.chemistry);
		for(std::size_t s=0;s<OrganSpeciesCount;++s){
			const double storage=input.distribution_volume_fraction[s]*volume/input.dt_s;
			double arriving=graph.exchange_in_m3_s[i]*input.inlet[s];
			for(const auto& edge:graph.upstream[i])
				arriving+=edge.second*result.concentration[edge.first][s];
			r[s]=(storage+leaving)*c[s]-arriving-storage*input.previous[i][s]
				-volume*rates.reaction[s];
		}
		return rates;
	};
	// Local Newton with a finite-difference Jacobian and a nonnegativity-
	// preserving backtracking line search. Returns the largest scaled update.
	const auto solve_cell=[&](std::size_t i){
		OrganConcentration c=result.concentration[i],r{};
		residual(i,c,r);
		const auto norm=[&](const OrganConcentration& value){
			double sum=0.;
			for(std::size_t s=0;s<OrganSpeciesCount;++s){
				const double reference=(input.distribution_volume_fraction[s]*graph.volume_m3[i]
					/input.dt_s+graph.face_outflow_m3_s[i]+graph.exchange_out_m3_s[i])*scale[s];
				sum+=std::pow(value[s]/reference,2);
			}
			return std::sqrt(sum);
		};
		double current=norm(r),change=0.;
		bool converged=current<=1e-12;
		for(int iteration=0;iteration<60&&!converged;++iteration){
			std::array<std::array<double,OrganSpeciesCount>,OrganSpeciesCount> jacobian{};
			for(std::size_t s=0;s<OrganSpeciesCount;++s){
				OrganConcentration shifted=c,shifted_r{};
				const double step=1e-7*std::max(std::abs(c[s]),1e-3*scale[s]);
				shifted[s]+=step;
				residual(i,shifted,shifted_r);
				for(std::size_t row=0;row<OrganSpeciesCount;++row)
					jacobian[row][s]=(shifted_r[row]-r[row])/step;
			}
			OrganConcentration delta=r;
			if(!SolveDense(jacobian,delta))break;
			double alpha=1.;
			OrganConcentration trial{},trial_r{};
			double trial_norm=current;
			for(int search=0;search<30;++search){
				for(std::size_t s=0;s<OrganSpeciesCount;++s)trial[s]=std::max(c[s]-alpha*delta[s],0.);
				residual(i,trial,trial_r);
				trial_norm=norm(trial_r);
				if(trial_norm<current||trial_norm<=1e-12)break;
				alpha*=0.5;
			}
			change=0.;
			for(std::size_t s=0;s<OrganSpeciesCount;++s)
				change=std::max(change,std::abs(trial[s]-c[s])/scale[s]);
			c=trial;r=trial_r;current=trial_norm;
			converged=current<=1e-10||change<=1e-13;
		}
		if(!converged)++result.newton_failures;
		const double moved=[&]{
			double value=0.;
			for(std::size_t s=0;s<OrganSpeciesCount;++s)
				value=std::max(value,std::abs(c[s]-result.concentration[i][s])/scale[s]);
			return value;}();
		result.concentration[i]=c;
		return moved;
	};

	const auto components=UpstreamFirstComponents(graph);
	result.components=components.size();
	for(const auto& component:components){
		result.largest_component=std::max(result.largest_component,component.size());
		if(component.size()==1){solve_cell(component.front());continue;}
		for(int sweep=0;sweep<500;++sweep){
			double moved=0.;
			for(const auto cell:component)moved=std::max(moved,solve_cell(cell));
			++result.gauss_seidel_sweeps;
			if(moved<=1e-11)break;
		}
	}

	for(std::size_t i=0;i<n;++i){
		const auto& c=result.concentration[i];
		result.rates[i]=EvaluateHepaticRates(c,input.state[i],input.kinetics,input.chemistry);
		for(std::size_t s=0;s<OrganSpeciesCount;++s){
			result.inflow[s]+=graph.exchange_in_m3_s[i]*input.inlet[s];
			result.outflow[s]+=(graph.exchange_out_m3_s[i]+graph.boundary_outflow_m3_s[i])*c[s];
			result.storage_rate[s]+=input.distribution_volume_fraction[s]*graph.volume_m3[i]
				*(c[s]-input.previous[i][s])/input.dt_s;
			result.reaction[s]+=graph.volume_m3[i]*result.rates[i].reaction[s];
		}
	}
	for(std::size_t s=0;s<OrganSpeciesCount;++s)
		result.balance_defect[s]=result.inflow[s]-result.outflow[s]-result.storage_rate[s]
			+result.reaction[s];
	return result;
}

} // namespace iga

#endif
