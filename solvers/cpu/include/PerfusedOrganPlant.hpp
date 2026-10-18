#ifndef IGA_PERFUSED_ORGAN_PLANT_HPP
#define IGA_PERFUSED_ORGAN_PLANT_HPP

#include "PerfusedOrganHemodynamics.hpp"
#include "PerfusedOrganKinetics.hpp"
#include "PerfusedOrganTransport.hpp"
#include "VascularCoupling.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

// A machine-perfused organ: the shared VCA external circuit (reservoir,
// pressure-controlled pump, dialysis, infusion, boluses) feeds an organ-layer
// gas exchanger and the tree/Darcy/transport organ; slow tissue states evolve
// between quasi-steady organ solves.

// Geometry and discretisation shared by every step.
struct PerfusedOrganModel
{
	NativeTetMesh mesh;
	std::vector<PerfusedOrganTree> trees;
	std::size_t zones=1;
	std::vector<double> cell_volume;
	// Pool of each cell for the lumped resolutions (identity for "cell").
	std::vector<std::size_t> territory_pool;
	std::size_t territory_pools=0;
};

inline void FinalizePerfusedOrganModel(PerfusedOrganModel& model)
{
	const std::size_t cells=model.mesh.cells.size();
	model.cell_volume.resize(cells);
	for(std::size_t i=0;i<cells;++i)model.cell_volume[i]=PerfusedOrganCellVolume(model.mesh,i);
	model.territory_pool.assign(cells,0);
	model.territory_pools=0;
	const auto inflow=std::find_if(model.trees.begin(),model.trees.end(),
		[](const PerfusedOrganTree& tree){return tree.inflow;});
	if(inflow==model.trees.end())throw std::invalid_argument("perfused organ needs an inflow tree");
	std::vector<char> assigned(cells,0);
	for(const auto& territory:inflow->territories){
		for(const auto& entry:territory.second){
			if(assigned[entry.first])throw std::invalid_argument("inflow territories overlap");
			assigned[entry.first]=1;
			model.territory_pool[entry.first]=model.territory_pools;
		}
		++model.territory_pools;
	}
	if(std::find(assigned.begin(),assigned.end(),0)!=assigned.end())
		throw std::invalid_argument("inflow territories do not cover every cell");
}

struct OrganStepResult
{
	PerfusedOrganFlowResult hemodynamics;
	OrganTransportResult species;  // node view: cell-major, zone-minor
	std::vector<double> node_volume;
	std::size_t pools=0;
};

// One quasi-steady organ solve at the requested transport resolution
// ("cell", "portal_territory" or "organ"). Lumped resolutions volume-average
// the node inputs into pools and broadcast pool results back.
inline OrganStepResult StepPerfusedOrgan(const PerfusedOrganModel& model,
	const PerfusedOrganFlowInput& flow,OrganTransportInput transport,const std::string& resolution)
{
	const std::size_t cells=model.mesh.cells.size(),zones=model.zones,nodes=cells*zones;
	OrganStepResult result;
	result.hemodynamics=SolvePerfusedOrganFlow(model.mesh,model.trees,flow);
	std::vector<std::size_t> group(cells,0);
	std::size_t groups=cells;
	if(resolution=="cell")std::iota(group.begin(),group.end(),0);
	else if(resolution=="organ")groups=1;
	else if(resolution=="portal_territory"){group=model.territory_pool;groups=model.territory_pools;}
	else throw std::invalid_argument("resolution must be cell, portal_territory or organ");
	result.pools=groups;
	const auto& h=result.hemodynamics;
	const auto cell_graph=BuildTetTransportGraph(model.mesh,h.darcy.conservative_face_flow_m3_s,
		h.cell_exchange_in_m3_s,h.cell_exchange_out_m3_s);
	const bool lumped=resolution!="cell";
	const auto graph=SplitIntoLobularZones(lumped?AggregateTransportGraph(cell_graph,group,groups):cell_graph,zones);
	if(transport.previous.size()!=nodes||transport.state.size()!=nodes)
		throw std::invalid_argument("organ step node inputs do not match the model");
	if(lumped){
		OrganTransportInput pooled=transport;
		pooled.previous.assign(groups*zones,OrganConcentration{});
		pooled.state.assign(groups*zones,{0.,0.});
		std::vector<double> pool_volume(groups*zones,0.);
		for(std::size_t i=0;i<cells;++i)
			for(std::size_t z=0;z<zones;++z){
				const std::size_t node=i*zones+z,pool=group[i]*zones+z;
				const double v=model.cell_volume[i];
				pool_volume[pool]+=v;
				for(std::size_t s=0;s<OrganSpeciesCount;++s)pooled.previous[pool][s]+=v*transport.previous[node][s];
				pooled.state[pool].viability+=v*transport.state[node].viability;
				pooled.state[pool].glycogen_mol_m3+=v*transport.state[node].glycogen_mol_m3;
			}
		for(std::size_t pool=0;pool<groups*zones;++pool){
			for(auto& value:pooled.previous[pool])value/=pool_volume[pool];
			pooled.state[pool].viability/=pool_volume[pool];
			pooled.state[pool].glycogen_mol_m3/=pool_volume[pool];
		}
		transport=std::move(pooled);
	}
	const auto pooled_species=SolveOrganTransport(graph,transport);
	result.species=pooled_species;
	if(lumped){
		result.species.concentration.resize(nodes);
		result.species.rates.resize(nodes);
		for(std::size_t i=0;i<cells;++i)
			for(std::size_t z=0;z<zones;++z){
				result.species.concentration[i*zones+z]=pooled_species.concentration[group[i]*zones+z];
				result.species.rates[i*zones+z]=pooled_species.rates[group[i]*zones+z];
			}
	}
	result.node_volume.resize(nodes);
	for(std::size_t i=0;i<cells;++i)
		for(std::size_t z=0;z<zones;++z)result.node_volume[i*zones+z]=model.cell_volume[i]/double(zones);
	return result;
}

// Physical parameters of the plant (hidden in the truth, guessed by tools).
struct OrganPlantParameters
{
	double viscosity_pa_s=3e-3;
	std::map<std::string,double> outflow_pressure_pa;          // per outflow tree
	std::map<std::string,double> specific_conductance_inv_pa_s; // per tree
	double mobility_m2_pa_s=5e-11;
	OrganConcentration distribution_volume_fraction{};
	HepaticKinetics kinetics;
	PerfusateChemistry chemistry;     // oxygen.hemoglobin_g_dl is replaced by Hct/3 of the reservoir
	double hemoglobin_per_hematocrit=1./3.; // g/dL per % haematocrit
	// Membrane gas exchanger: gas pO2 = FiO2 (Patm - pH2O); the O2 gas is
	// air blended with carbogen (o2 carbogen_o2_fraction, CO2 carbogen_co2_fraction),
	// so the CO2 partial pressure of the gas rises with FiO2. Arterial O2 and
	// CO2 contents are efficiency-weighted blends of equilibrium and venous values.
	double atmospheric_pressure_mmhg=760.,water_vapour_mmhg=47.;
	double carbogen_o2_fraction=0.95,carbogen_co2_fraction=0.05,air_o2_fraction=0.21;
	double oxygenator_efficiency=0.95;
	// Injury: viable fraction falls while the energy deficit exceeds a threshold.
	double injury_deficit_threshold=0.2,injury_rate_s_inv=2e-4;
	// Sinusoidal oedema (per cell) grows with energy deficit and with tissue
	// pressure above a threshold, resolves first-order, and scales exchange
	// conductance and mobility by (1 - E)^exponent.
	double oedema_deficit_rate_s_inv=5e-5,oedema_pressure_threshold_pa=1000.;
	double oedema_pressure_rate_inv_pa_s=1e-7,oedema_resolution_s_inv=2e-5;
	double oedema_maximum=0.8,oedema_exponent=2.;
	// Vasodilator in the reservoir [mg/L]: first-order elimination and an
	// Emax effect on exchange conductance.
	double vasodilator_elimination_s_inv=2e-4,vasodilator_emax=0.5,vasodilator_ec50_mg_l=5.;
	// Optional costs of aggressive control (all off by default):
	//  - hyperoxic vasoconstriction: exchange conductance divided by
	//    1 + c (pO2_arterial - threshold)_+ / 100 mmHg;
	//  - oxidative injury: viability decays at k (pO2_node - threshold)_+ / 100 mmHg;
	//  - shear injury: viability decays at k (tissue pressure - threshold)_+ per Pa.
	double hyperoxic_constriction_per_100mmhg=0.,hyperoxic_constriction_threshold_mmhg=150.;
	double oxidative_injury_rate_s_inv=0.,oxidative_injury_threshold_mmhg=200.;
	double shear_injury_rate_inv_pa_s=0.,shear_injury_threshold_pa=1000.;
};

struct OrganPlantControls
{
	double inlet_pressure_pa=0.;
	double fio2=0.5;
	double dialysis_in_m3_s=0.,dialysis_out_m3_s=0.;
	std::map<std::string,double> dialysate;      // composition of dialysis inflow
	std::map<std::string,double> infusion_rates; // per species, amount/s
	std::map<std::string,double> bolus;          // per species, amount at the step start
	double bolus_volume_m3=0.;
	double vasodilator_bolus_mg=0.;
	double temperature_c=37.;
	// Optional independent CO2 fraction of the dry sweep gas. NaN keeps the
	// Model41 air/carbogen blend, in which gas CO2 follows FiO2.
	double gas_co2_fraction=std::numeric_limits<double>::quiet_NaN();
};

struct OrganPlantState
{
	double time_s=0.;
	ReservoirState reservoir;            // species by PerfusedOrganSpeciesNames
	double vasodilator_mg_l=0.;
	std::vector<OrganConcentration> tissue; // per node
	std::vector<double> viability,glycogen_mol_m3; // per node
	std::vector<double> oedema;          // per cell
};

struct PerfusateReading
{
	OrganConcentration concentration{};
	double po2_mmhg=0.,saturation=0.,ph=0.,pco2_mmhg=0.,bicarbonate_mol_m3=0.;
};

struct OrganPlantRecord
{
	double time_s=0.,flow_m3_s=0.,inlet_pressure_pa=0.,reservoir_volume_m3=0.,hematocrit_percent=0.;
	double vasodilator_mg_l=0.;
	PerfusateReading arterial,venous;
	double hypoxic_volume_fraction=0.,mean_viability=0.,mean_oedema=0.;
	double last_zone_po2_p05=0.,energy_deficit_mean=0.,energy_deficit_p95=0.,o2_uptake_mol_s=0.;
	OrganConcentration organ_uptake{};
	// Diagnostics of the substep's organ solve; the substep is committed either
	// way, so callers decide whether a nonzero failure count is acceptable.
	std::size_t newton_failures=0,gauss_seidel_sweeps=0,gauss_seidel_unconverged=0;
	double flow_residual_norm=0.,exchange_imbalance_m3_s=0.;
	OrganConcentration balance_defect{};
};

// Rejects parameters that would otherwise be clamped or turn into NaN silently.
inline void ValidatePerfusedOrganPlantParameters(const OrganPlantParameters& p)
{
	const auto nonnegative=[](double value){return value>=0.&&std::isfinite(value);};
	if(!(p.air_o2_fraction>0.&&p.air_o2_fraction<p.carbogen_o2_fraction&&p.carbogen_o2_fraction<=1.))
		throw std::invalid_argument("plant gas fractions need 0 < air O2 < carbogen O2 <= 1");
	if(!nonnegative(p.carbogen_co2_fraction)||p.carbogen_o2_fraction+p.carbogen_co2_fraction>1.+1e-12)
		throw std::invalid_argument("plant carbogen CO2 fraction is invalid");
	if(!nonnegative(p.water_vapour_mmhg)||!(p.atmospheric_pressure_mmhg>p.water_vapour_mmhg)
		||!std::isfinite(p.atmospheric_pressure_mmhg))
		throw std::invalid_argument("plant gas pressures need atmospheric > water vapour >= 0");
	if(!(p.oxygenator_efficiency>=0.&&p.oxygenator_efficiency<=1.))
		throw std::invalid_argument("oxygenator efficiency must lie in [0,1]");
	if(!(p.injury_deficit_threshold>=0.&&p.injury_deficit_threshold<=1.))
		throw std::invalid_argument("injury deficit threshold must lie in [0,1]");
	if(!(p.oedema_maximum>=0.&&p.oedema_maximum<1.))
		throw std::invalid_argument("oedema maximum must lie in [0,1)");
	for(const double value:{p.hemoglobin_per_hematocrit,p.injury_rate_s_inv,p.oedema_deficit_rate_s_inv,
		p.oedema_pressure_rate_inv_pa_s,p.oedema_resolution_s_inv,p.oedema_exponent,
		p.vasodilator_elimination_s_inv,p.vasodilator_emax})
		if(!nonnegative(value))throw std::invalid_argument("plant rates and exponents must be finite and nonnegative");
	if(!std::isfinite(p.oedema_pressure_threshold_pa))
		throw std::invalid_argument("oedema pressure threshold must be finite");
	if(!(p.vasodilator_ec50_mg_l>0.)||!std::isfinite(p.vasodilator_ec50_mg_l))
		throw std::invalid_argument("vasodilator EC50 must be positive");
	if(!nonnegative(p.hyperoxic_constriction_per_100mmhg)||!nonnegative(p.oxidative_injury_rate_s_inv)
		||!nonnegative(p.shear_injury_rate_inv_pa_s))
		throw std::invalid_argument("plant control-cost rates must be finite and nonnegative");
	if(!std::isfinite(p.hyperoxic_constriction_threshold_mmhg)||!std::isfinite(p.oxidative_injury_threshold_mmhg)
		||!std::isfinite(p.shear_injury_threshold_pa))
		throw std::invalid_argument("plant control-cost thresholds must be finite");
}

// Restored or initial states must be ones the advance itself can produce.
inline void ValidatePerfusedOrganPlantState(const OrganPlantState& state)
{
	if(!std::isfinite(state.time_s)||!(state.vasodilator_mg_l>=0.)||!std::isfinite(state.vasodilator_mg_l))
		throw std::invalid_argument("plant state time or vasodilator is invalid");
	for(const auto& c:state.tissue)
		for(const double value:c)
			if(!(value>=0.)||!std::isfinite(value))
				throw std::invalid_argument("plant tissue concentrations must be finite and nonnegative");
	for(const double value:state.viability)
		if(!(value>=0.&&value<=1.))throw std::invalid_argument("plant viability must lie in [0,1]");
	for(const double value:state.glycogen_mol_m3)
		if(!(value>=0.)||!std::isfinite(value))
			throw std::invalid_argument("plant glycogen must be finite and nonnegative");
	for(const double value:state.oedema)
		if(!(value>=0.&&value<1.))throw std::invalid_argument("plant oedema must lie in [0,1)");
}

inline PerfusateReading ReadPerfusate(const OrganConcentration& c,const PerfusateChemistry& chemistry)
{
	PerfusateReading reading;
	reading.concentration=c;
	const auto acid=SolveAcidBase(c[OrganBufferBase],c[OrganCarbonDioxide],chemistry);
	reading.ph=acid.ph;
	reading.pco2_mmhg=acid.pco2_mmhg;
	reading.bicarbonate_mol_m3=acid.bicarbonate_mol_m3;
	reading.po2_mmhg=PartialPressureO2(c[OrganOxygen],acid.ph,chemistry);
	reading.saturation=OxygenHillSaturation(reading.po2_mmhg,BohrShifted(acid.ph,chemistry).p50_mmhg,
		chemistry.oxygen.hill_exponent);
	return reading;
}

// Total CO2 at which the perfusate with this buffer base has the given pCO2.
inline double Co2TotalAtPco2(double buffer_base,double pco2_mmhg,const PerfusateChemistry& chemistry)
{
	const double dissolved=chemistry.co2_solubility_mol_m3_mmhg*std::max(pco2_mmhg,0.);
	if(dissolved<=0.)return 0.;
	// HCO3-(pH) = BB - A-(pH) must equal dissolved 10^(pH - pK); bisection on pH.
	const double beta=NonBicarbonateBuffer(chemistry),anion=BufferAnionAt74(chemistry);
	double low=5.5,high=9.;
	for(int i=0;i<200;++i){
		const double mid=0.5*(low+high);
		const double difference=(buffer_base-anion-beta*(mid-7.4))
			-dissolved*std::pow(10.,mid-chemistry.bicarbonate_pk);
		(difference>0.?low:high)=mid;
	}
	const double ph=0.5*(low+high);
	return std::max(buffer_base-anion-beta*(ph-7.4),0.)+dissolved;
}

inline OrganConcentration Oxygenate(const OrganConcentration& venous,double fio2,
	const OrganPlantParameters& p,const PerfusateChemistry& chemistry,
	double gas_co2_fraction=std::numeric_limits<double>::quiet_NaN())
{
	if(!(fio2>=p.air_o2_fraction-1e-12&&fio2<=p.carbogen_o2_fraction+1e-12))
		throw std::invalid_argument("FiO2 must lie between air and carbogen");
	const double dry=p.atmospheric_pressure_mmhg-p.water_vapour_mmhg;
	const double carbogen=(fio2-p.air_o2_fraction)/(p.carbogen_o2_fraction-p.air_o2_fraction);
	const double gas_po2=fio2*dry;
	if(!std::isnan(gas_co2_fraction)&&!(gas_co2_fraction>=0.&&gas_co2_fraction<=0.2))
		throw std::invalid_argument("gas CO2 fraction must be in [0, 0.2]");
	if(!std::isnan(gas_co2_fraction)&&fio2+gas_co2_fraction>1.+1e-12)
		throw std::invalid_argument("FiO2 plus the gas CO2 fraction must not exceed 1");
	const double gas_pco2=std::isnan(gas_co2_fraction)
		?std::max(carbogen,0.)*p.carbogen_co2_fraction*dry:gas_co2_fraction*dry;
	OrganConcentration arterial=venous;
	const double eta=std::clamp(p.oxygenator_efficiency,0.,1.);
	const double co2_equilibrium=Co2TotalAtPco2(venous[OrganBufferBase],gas_pco2,chemistry);
	arterial[OrganCarbonDioxide]=eta*co2_equilibrium+(1.-eta)*venous[OrganCarbonDioxide];
	const double ph=SolveAcidBase(arterial[OrganBufferBase],arterial[OrganCarbonDioxide],chemistry).ph;
	arterial[OrganOxygen]=eta*OxygenContent(gas_po2,ph,chemistry)+(1.-eta)*venous[OrganOxygen];
	return arterial;
}

inline OrganConcentration ReservoirComposition(const ReservoirState& reservoir)
{
	OrganConcentration c{};
	const auto& names=PerfusedOrganSpeciesNames();
	for(std::size_t s=0;s<OrganSpeciesCount;++s){
		const auto found=reservoir.species.find(names[s]);
		if(found==reservoir.species.end())
			throw std::invalid_argument(std::string("reservoir lacks species ")+names[s]);
		c[s]=found->second;
	}
	return c;
}

inline CouplingDefinition PerfusedOrganCircuitDefinition(const OrganPlantState& state,
	const OrganPlantControls& controls)
{
	CouplingDefinition definition;
	definition.mode=SimulationScopeMode::VcaClosedLoop;
	auto& circuit=definition.external_circuit;
	circuit.reservoir.enabled=true;
	circuit.reservoir.volume_m3=state.reservoir.volume_m3;
	circuit.reservoir.temperature_c=controls.temperature_c;
	circuit.reservoir.species=state.reservoir.species;
	circuit.pump.mode="pressure_control";
	circuit.pump.pressure_pa=controls.inlet_pressure_pa;
	circuit.dialyzer.enabled=controls.dialysis_in_m3_s>0.||controls.dialysis_out_m3_s>0.;
	circuit.dialyzer.inflow_m3_s=controls.dialysis_in_m3_s;
	circuit.dialyzer.outflow_m3_s=controls.dialysis_out_m3_s;
	circuit.dialyzer.dialysate=controls.dialysate;
	circuit.infusion_rates=controls.infusion_rates;
	definition.perfusate.oxygen.hematocrit_percent=state.reservoir.hematocrit_percent;
	return definition;
}

// Advances the plant by `duration_s` in `substeps` equal organ solves.
// `mobility_scale` and `exchange_scale` are the static per-cell baseline
// heterogeneity (empty = uniform). Returns one record per substep.
inline std::vector<OrganPlantRecord> AdvancePerfusedOrganPlant(const PerfusedOrganModel& model,
	const OrganPlantParameters& p,const OrganPlantControls& controls,
	const std::vector<double>& mobility_scale,const std::vector<double>& exchange_scale,
	OrganPlantState& state,double duration_s,std::size_t substeps,const std::string& resolution,
	OrganStepResult* last_step=nullptr)
{
	const std::size_t cells=model.mesh.cells.size(),nodes=cells*model.zones;
	if(!(duration_s>0.)||substeps==0)throw std::invalid_argument("plant advance needs a positive duration and substeps");
	if(state.tissue.size()!=nodes||state.viability.size()!=nodes||state.glycogen_mol_m3.size()!=nodes
		||state.oedema.size()!=cells)
		throw std::invalid_argument("plant state does not match the model");
	ValidatePerfusedOrganPlantParameters(p);
	ValidatePerfusedOrganPlantState(state);
	if(!mobility_scale.empty()&&mobility_scale.size()!=cells)throw std::invalid_argument("mobility scale size");
	if(!exchange_scale.empty()&&exchange_scale.size()!=cells)throw std::invalid_argument("exchange scale size");
	VcaExternalCircuit circuit(PerfusedOrganCircuitDefinition(state,controls));
	circuit.RestoreState(state.reservoir);
	circuit.ApplyBolus(controls.bolus,controls.bolus_volume_m3);
	if(controls.vasodilator_bolus_mg>0.)
		state.vasodilator_mg_l+=controls.vasodilator_bolus_mg/(circuit.State().volume_m3*1000.);
	const double dt=duration_s/double(substeps);
	std::vector<OrganPlantRecord> records;
	for(std::size_t step=0;step<substeps;++step){
		const auto& reservoir=circuit.State();
		PerfusateChemistry chemistry=p.chemistry;
		chemistry.oxygen.hematocrit_percent=reservoir.hematocrit_percent;
		chemistry.oxygen.hemoglobin_g_dl=reservoir.hematocrit_percent*p.hemoglobin_per_hematocrit;
		const auto venous_in=ReservoirComposition(reservoir);
		const auto arterial=Oxygenate(venous_in,controls.fio2,p,chemistry,controls.gas_co2_fraction);
		const double drug=state.vasodilator_mg_l;
		const double dilation=1.+p.vasodilator_emax*drug/(p.vasodilator_ec50_mg_l+drug);
		const double arterial_po2=PartialPressureO2(arterial[OrganOxygen],
			SolveAcidBase(arterial[OrganBufferBase],arterial[OrganCarbonDioxide],chemistry).ph,chemistry);
		const double constriction=1./(1.+p.hyperoxic_constriction_per_100mmhg
			*std::max(arterial_po2-p.hyperoxic_constriction_threshold_mmhg,0.)/100.);

		PerfusedOrganFlowInput flow;
		flow.viscosity_pa_s=p.viscosity_pa_s;
		for(const auto& tree:model.trees){
			const auto g=p.specific_conductance_inv_pa_s.find(tree.id);
			if(g==p.specific_conductance_inv_pa_s.end())throw std::invalid_argument("no conductance for tree "+tree.id);
			flow.specific_conductance_inv_pa_s.push_back(g->second);
			if(tree.inflow)flow.root_pressure_pa.push_back(controls.inlet_pressure_pa);
			else{
				const auto out=p.outflow_pressure_pa.find(tree.id);
				if(out==p.outflow_pressure_pa.end())throw std::invalid_argument("no outflow pressure for tree "+tree.id);
				flow.root_pressure_pa.push_back(out->second);
			}
		}
		flow.mobility_by_cell_m2_pa_s.resize(cells);
		flow.exchange_scale_by_cell.resize(cells);
		for(std::size_t i=0;i<cells;++i){
			const double compression=std::pow(1.-std::clamp(state.oedema[i],0.,p.oedema_maximum),p.oedema_exponent);
			flow.mobility_by_cell_m2_pa_s[i]=p.mobility_m2_pa_s*(mobility_scale.empty()?1.:mobility_scale[i])*compression;
			flow.exchange_scale_by_cell[i]=(exchange_scale.empty()?1.:exchange_scale[i])*compression*dilation
				*constriction;
		}
		OrganTransportInput transport;
		transport.dt_s=dt;
		transport.inlet=arterial;
		transport.distribution_volume_fraction=p.distribution_volume_fraction;
		transport.previous=state.tissue;
		transport.state.resize(nodes);
		for(std::size_t n=0;n<nodes;++n)transport.state[n]={state.viability[n],state.glycogen_mol_m3[n]};
		transport.kinetics=p.kinetics;
		transport.kinetics.temperature_c=controls.temperature_c;
		transport.chemistry=chemistry;
		auto organ=StepPerfusedOrgan(model,flow,transport,resolution);
		const auto& species=organ.species;
		const double q=std::accumulate(organ.hemodynamics.cell_exchange_in_m3_s.begin(),
			organ.hemodynamics.cell_exchange_in_m3_s.end(),0.);

		// Reservoir: organ uptake at the arterial composition; O2 and CO2 are
		// reset by the gas exchanger on every pass.
		std::map<std::string,double> uptake,reset;
		const auto& names=PerfusedOrganSpeciesNames();
		OrganConcentration organ_uptake{};
		for(std::size_t s=0;s<OrganSpeciesCount;++s){
			organ_uptake[s]=species.inflow[s]-species.outflow[s];
			uptake[names[s]]=organ_uptake[s];
		}
		reset[names[OrganOxygen]]=arterial[OrganOxygen];
		reset[names[OrganCarbonDioxide]]=arterial[OrganCarbonDioxide];
		// Reservoir mass balance: dM/dt = q c_arterial - U - q c_reservoir, which
		// AdvanceWithUptake applies for the reset species; for the others the
		// arterial value is the reservoir value and dM/dt = -U.
		circuit.AdvanceWithUptake(uptake,reset,q,dt,state.time_s);

		// Slow tissue states.
		std::vector<double> cell_deficit(cells,0.);
		for(std::size_t n=0;n<nodes;++n){
			const auto& rates=species.rates[n];
			const double excess=std::max(rates.energy_deficit-p.injury_deficit_threshold,0.);
			const double hyperoxia=std::max(rates.po2_mmhg-p.oxidative_injury_threshold_mmhg,0.)/100.;
			const double shear=std::max(organ.hemodynamics.cell_mean_pressure_pa[n/model.zones]
				-p.shear_injury_threshold_pa,0.);
			state.viability[n]*=std::exp(-(p.injury_rate_s_inv*excess+p.oxidative_injury_rate_s_inv*hyperoxia
				+p.shear_injury_rate_inv_pa_s*shear)*dt);
			state.glycogen_mol_m3[n]=std::max(state.glycogen_mol_m3[n]+dt*rates.glycogen_rate,0.);
			cell_deficit[n/model.zones]+=rates.energy_deficit/double(model.zones);
		}
		for(std::size_t i=0;i<cells;++i){
			const double pressure_excess=std::max(organ.hemodynamics.cell_mean_pressure_pa[i]
				-p.oedema_pressure_threshold_pa,0.);
			const double drive=p.oedema_deficit_rate_s_inv*cell_deficit[i]
				+p.oedema_pressure_rate_inv_pa_s*pressure_excess;
			const double rate=p.oedema_resolution_s_inv+drive;
			const double target=rate>0.?drive/rate:0.;
			state.oedema[i]=std::clamp(target+(state.oedema[i]-target)*std::exp(-rate*dt),0.,p.oedema_maximum);
		}
		state.tissue=species.concentration;
		state.vasodilator_mg_l*=std::exp(-p.vasodilator_elimination_s_inv*dt);
		state.time_s+=dt;

		OrganPlantRecord record;
		record.time_s=state.time_s;
		record.flow_m3_s=q;
		record.inlet_pressure_pa=controls.inlet_pressure_pa;
		record.reservoir_volume_m3=circuit.State().volume_m3;
		record.hematocrit_percent=circuit.State().hematocrit_percent;
		record.vasodilator_mg_l=state.vasodilator_mg_l;
		record.arterial=ReadPerfusate(arterial,chemistry);
		// Venous effluent: the organ outflow mixed (before the reservoir).
		OrganConcentration venous{};
		double out=0.;
		for(const auto& exchange:organ.hemodynamics.exchanges){
			if(exchange.flow_m3_s>=0.)continue;
			out-=exchange.flow_m3_s;
			for(std::size_t s=0;s<OrganSpeciesCount;++s)
				venous[s]-=exchange.flow_m3_s*species.concentration[exchange.cell*model.zones+model.zones-1][s];
		}
		for(auto& value:venous)value=out>0.?value/out:0.;
		record.venous=ReadPerfusate(venous,chemistry);
		double hypoxic=0.,total=0.,viable=0.,deficit=0.,uptake_o2=0.;
		std::vector<std::pair<double,double>> deficits,last_zone;
		for(std::size_t n=0;n<nodes;++n){
			const double v=organ.node_volume[n];
			total+=v;
			if(species.rates[n].po2_mmhg<5.)hypoxic+=v;
			viable+=v*state.viability[n];
			deficit+=v*species.rates[n].energy_deficit;
			uptake_o2+=v*species.rates[n].o2_uptake;
			deficits.push_back({species.rates[n].energy_deficit,v});
			if(n%model.zones==model.zones-1)last_zone.push_back({species.rates[n].po2_mmhg,v});
		}
		const auto quantile=[](std::vector<std::pair<double,double>> values,double fraction){
			std::sort(values.begin(),values.end());
			double sum=0.,running=0.;
			for(const auto& item:values)sum+=item.second;
			for(const auto& item:values){running+=item.second;if(running>=fraction*sum)return item.first;}
			return values.back().first;
		};
		record.hypoxic_volume_fraction=hypoxic/total;
		record.mean_viability=viable/total;
		record.energy_deficit_mean=deficit/total;
		record.energy_deficit_p95=quantile(deficits,0.95);
		record.last_zone_po2_p05=quantile(last_zone,0.05);
		record.o2_uptake_mol_s=uptake_o2;
		double oedema=0.,volume=0.;
		for(std::size_t i=0;i<cells;++i){oedema+=model.cell_volume[i]*state.oedema[i];volume+=model.cell_volume[i];}
		record.mean_oedema=oedema/volume;
		record.organ_uptake=organ_uptake;
		record.newton_failures=species.newton_failures;
		record.gauss_seidel_sweeps=species.gauss_seidel_sweeps;
		record.gauss_seidel_unconverged=species.gauss_seidel_unconverged;
		record.flow_residual_norm=organ.hemodynamics.residual_norm;
		record.exchange_imbalance_m3_s=organ.hemodynamics.exchange_imbalance_m3_s;
		record.balance_defect=species.balance_defect;
		records.push_back(record);
		if(last_step&&step+1==substeps)*last_step=std::move(organ);
	}
	state.reservoir=circuit.State();
	return records;
}

} // namespace iga

#endif
