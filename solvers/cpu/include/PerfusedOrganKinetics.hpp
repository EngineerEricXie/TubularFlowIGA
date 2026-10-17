#ifndef IGA_PERFUSED_ORGAN_KINETICS_HPP
#define IGA_PERFUSED_ORGAN_KINETICS_HPP

#include "OxygenCapacity.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

namespace iga {

// Species carried by the perfusate. Concentrations are per volume of
// perfusate; reaction rates are per volume of tissue.
enum PerfusedOrganSpecies : std::size_t {
	OrganOxygen=0,      // total O2 content, dissolved + haemoglobin-bound [mol/m^3]
	OrganGlucose=1,     // [mol/m^3] (= mmol/L)
	OrganLactate=2,     // [mol/m^3]
	OrganBufferBase=3,  // HCO3- plus non-bicarbonate buffer anions [mol/m^3]
	OrganCarbonDioxide=4, // total CO2, HCO3- plus dissolved CO2 [mol/m^3]
	OrganInsulin=5,     // [mU/L]
	OrganSpeciesCount=6
};

using OrganConcentration=std::array<double,OrganSpeciesCount>;

inline const std::array<const char*,OrganSpeciesCount>& PerfusedOrganSpeciesNames()
{
	static const std::array<const char*,OrganSpeciesCount> names{{
		"o2_total_mol_m3","glucose_mol_m3","lactate_mol_m3",
		"buffer_base_mol_m3","co2_total_mol_m3","insulin_mU_L"}};
	return names;
}

// Perfusate chemistry that does not react in the tissue. Oxygen binding reuses
// the shared OxygenCapacity model (Hill saturation, 1.34 mL O2/g Hb, plasma
// solubility); its P50 is the value at pH 7.4 and shifts with the Bohr effect.
struct PerfusateChemistry
{
	OxygenCapacityParameters oxygen;
	double bohr_coefficient=-0.48;         // d log10(P50) / d pH
	double co2_solubility_mol_m3_mmhg=0.0307;
	double bicarbonate_pk=6.1;
	// Non-bicarbonate buffer A-(pH) = A0 + beta (pH - 7.4). NaN selects the
	// whole-blood estimates from haemoglobin: Van Slyke beta = 1.43 Hb + 7.7
	// mmol/(L pH) and A0 = 17.7 + 0.42 Hb mmol/L (normal buffer base
	// 41.7 + 0.42 Hb with HCO3- 24 mmol/L), Hb in g/dL.
	double non_bicarbonate_buffer_mol_m3_ph=std::numeric_limits<double>::quiet_NaN();
	double buffer_anion_at_7_4_mol_m3=std::numeric_limits<double>::quiet_NaN();
};

inline double NonBicarbonateBuffer(const PerfusateChemistry& chemistry)
{
	return std::isnan(chemistry.non_bicarbonate_buffer_mol_m3_ph)
		?1.43*chemistry.oxygen.hemoglobin_g_dl+7.7:chemistry.non_bicarbonate_buffer_mol_m3_ph;
}

inline double BufferAnionAt74(const PerfusateChemistry& chemistry)
{
	return std::isnan(chemistry.buffer_anion_at_7_4_mol_m3)
		?17.7+0.42*chemistry.oxygen.hemoglobin_g_dl:chemistry.buffer_anion_at_7_4_mol_m3;
}

struct PerfusateAcidBase
{
	double ph=7.4,bicarbonate_mol_m3=0.,co2_dissolved_mol_m3=0.,pco2_mmhg=0.;
};

// Solves HCO3-(pH) (1 + 10^(pK - pH)) = TCO2 with HCO3- = BB - A-(pH), which
// is strictly decreasing in pH, by safeguarded Newton on [5.5, 9].
inline PerfusateAcidBase SolveAcidBase(double buffer_base,double co2_total,
	const PerfusateChemistry& chemistry)
{
	const double beta=NonBicarbonateBuffer(chemistry),anion=BufferAnionAt74(chemistry);
	const double pk=chemistry.bicarbonate_pk;
	const auto bicarbonate=[&](double ph){return buffer_base-anion-beta*(ph-7.4);};
	const auto residual=[&](double ph){
		return bicarbonate(ph)*(1.+std::pow(10.,pk-ph))-std::max(co2_total,0.);};
	double low=5.5,high=9.;
	if(beta>0.)high=std::clamp(7.4+(buffer_base-anion)/beta,low,9.);
	double ph=std::clamp(7.4,low,high);
	if(residual(low)<=0.)ph=low;
	else if(residual(high)>=0.)ph=high;
	else for(int iteration=0;iteration<100;++iteration){
		const double r=residual(ph);
		if(r>0.)low=ph;else high=ph;
		const double factor=std::pow(10.,pk-ph);
		const double slope=-beta*(1.+factor)-bicarbonate(ph)*factor*std::log(10.);
		double next=slope<0.?ph-r/slope:0.5*(low+high);
		if(!(next>low&&next<high))next=0.5*(low+high);
		if(std::abs(next-ph)<1e-12){ph=next;break;}
		ph=next;
	}
	PerfusateAcidBase result;
	result.ph=ph;
	result.bicarbonate_mol_m3=std::max(bicarbonate(ph),0.);
	result.co2_dissolved_mol_m3=std::max(co2_total,0.)-result.bicarbonate_mol_m3;
	if(result.co2_dissolved_mol_m3<0.)result.co2_dissolved_mol_m3=0.;
	result.pco2_mmhg=result.co2_dissolved_mol_m3/chemistry.co2_solubility_mol_m3_mmhg;
	return result;
}

// Hepatic metabolism, model "hepatic_v1". Rates scale with the local viable
// fraction v in [0,1] and with Q10^((T-37)/10).
//  - O2 uptake  VO2 = v Vmax pO2/(Km + pO2) (supply-limited, Michaelis-Menten).
//  - ATP demand D = v Vmax atp_per_o2; the aerobic shortfall is covered by
//    anaerobic glycolysis (Pasteur effect) up to v glycolysis_capacity, one ATP
//    per lactate, fuelled by glycogen first and then perfusate glucose. The
//    remaining shortfall fraction is the energy deficit reported to the
//    caller's injury model.
//  - Lactate clearance v Vmax_L f_O2 L/(Km_L + L) needs oxygen; a fraction
//    becomes glucose (gluconeogenesis, two lactate per glucose).
//  - Glycogenolysis releases glucose, suppressed by insulin; insulin promotes
//    glucose uptake into glycogen. Insulin is cleared first-order.
//  - Lactate production releases one H+ per lactate, lowering buffer base;
//    lactate clearance consumes it. Aerobic metabolism adds total CO2 at the
//    respiratory quotient. pH, HCO3- and pCO2 follow from buffer base and total
//    CO2 with a linear non-bicarbonate buffer (SolveAcidBase).
struct HepaticKinetics
{
	double vmax_o2_mol_m3_s=0.;
	double km_o2_mmhg=2.;
	double atp_per_o2=5.;
	double glycolysis_capacity_mol_m3_s=0.;
	double glycogen_half_saturation_mol_m3=20.;
	double km_glucose_mol_m3=1.;
	double lactate_vmax_mol_m3_s=0.;
	double lactate_km_mol_m3=2.;
	double gluconeogenic_fraction=0.5;
	double glycogenolysis_vmax_mol_m3_s=0.;
	double glucose_uptake_vmax_mol_m3_s=0.;
	double insulin_half_effect_mU_L=50.;
	double insulin_clearance_s_inv=0.;
	double respiratory_quotient=0.8;
	double q10=2.;
	double temperature_c=37.;
};

struct HepaticCellState
{
	double viability=1.;     // viable fraction in [0,1]
	double glycogen_mol_m3=0.; // glucosyl units per tissue volume
};

struct HepaticRates
{
	OrganConcentration reaction{}; // mol/(m^3 tissue s), insulin mU/(L s)
	double po2_mmhg=0.,ph=0.,pco2_mmhg=0.,bicarbonate=0.;
	double o2_uptake=0.,lactate_production=0.,lactate_clearance=0.;
	double glycogen_rate=0.;       // d(glycogen)/dt
	double energy_deficit=0.;      // unmet fraction of ATP demand in [0,1]
};

inline OxygenCapacityParameters BohrShifted(double ph,const PerfusateChemistry& chemistry)
{
	auto parameters=chemistry.oxygen;
	parameters.p50_mmhg*=std::pow(10.,chemistry.bohr_coefficient*(std::clamp(ph,6.5,8.0)-7.4));
	return parameters;
}

inline double OxygenContent(double po2_mmhg,double ph,const PerfusateChemistry& chemistry)
{
	return OxygenFromPo2(po2_mmhg,BohrShifted(ph,chemistry)).total_oxygen_mol_m3;
}

// Inverts the monotone content curve by safeguarded bisection/secant.
inline double PartialPressureO2(double content,double ph,const PerfusateChemistry& chemistry)
{
	if(content<=0.)return 0.;
	const double solubility=chemistry.oxygen.oxygen_solubility_ml_dl_mmhg*10.
		/chemistry.oxygen.gas_molar_volume_ml_mmol; // mol/(m^3 mmHg)
	double low=0.,high=content/solubility;
	double value=std::min(high,chemistry.oxygen.p50_mmhg);
	for(int iteration=0;iteration<200;++iteration){
		const double residual=OxygenContent(value,ph,chemistry)-content;
		if(std::abs(residual)<=1e-13*std::max(content,1e-12))break;
		if(residual>0.)high=value;else low=value;
		const double step=std::max(1e-9,1e-7*value);
		const double slope=(OxygenContent(value+step,ph,chemistry)-OxygenContent(value,ph,chemistry))/step;
		double next=slope>0.?value-residual/slope:0.5*(low+high);
		if(!(next>low&&next<high))next=0.5*(low+high);
		if(std::abs(next-value)<=1e-12*std::max(1.,value)){value=next;break;}
		value=next;
	}
	return value;
}

inline void ValidateHepaticKinetics(const HepaticKinetics& k,const PerfusateChemistry& chemistry)
{
	if(!(NonBicarbonateBuffer(chemistry)>=0.)||!std::isfinite(BufferAnionAt74(chemistry)))
		throw std::invalid_argument("perfusate buffer constants are invalid");
	ValidateOxygenCapacity(chemistry.oxygen);
	const double values[]={k.vmax_o2_mol_m3_s,k.atp_per_o2,k.glycolysis_capacity_mol_m3_s,
		k.lactate_vmax_mol_m3_s,k.glycogenolysis_vmax_mol_m3_s,
		k.glucose_uptake_vmax_mol_m3_s,k.insulin_clearance_s_inv,k.respiratory_quotient};
	for(const double value:values)
		if(!(value>=0.)||!std::isfinite(value))
			throw std::invalid_argument("hepatic kinetics rates must be finite and nonnegative");
	const double positive[]={k.km_o2_mmhg,k.glycogen_half_saturation_mol_m3,k.km_glucose_mol_m3,
		k.lactate_km_mol_m3,k.insulin_half_effect_mU_L,k.q10,chemistry.co2_solubility_mol_m3_mmhg};
	for(const double value:positive)
		if(!(value>0.)||!std::isfinite(value))
			throw std::invalid_argument("hepatic kinetics constants must be positive");
	if(!(k.gluconeogenic_fraction>=0.&&k.gluconeogenic_fraction<=1.))
		throw std::invalid_argument("gluconeogenic fraction must be in [0,1]");
	if(!std::isfinite(k.temperature_c))
		throw std::invalid_argument("hepatic kinetics temperature is invalid");
}

inline HepaticRates EvaluateHepaticRates(const OrganConcentration& concentration,
	const HepaticCellState& state,const HepaticKinetics& k,const PerfusateChemistry& chemistry)
{
	HepaticRates rates;
	OrganConcentration c=concentration;
	for(auto& value:c)value=std::max(value,0.);
	const double v=std::clamp(state.viability,0.,1.);
	const double glycogen=std::max(state.glycogen_mol_m3,0.);
	const double thermal=std::pow(k.q10,(k.temperature_c-37.)/10.);
	const auto acid_base=SolveAcidBase(c[OrganBufferBase],c[OrganCarbonDioxide],chemistry);
	rates.ph=acid_base.ph;
	rates.pco2_mmhg=acid_base.pco2_mmhg;
	rates.bicarbonate=acid_base.bicarbonate_mol_m3;
	rates.po2_mmhg=PartialPressureO2(c[OrganOxygen],rates.ph,chemistry);
	const double oxygen_fraction=rates.po2_mmhg/(k.km_o2_mmhg+rates.po2_mmhg);
	const double vmax=v*k.vmax_o2_mol_m3_s*thermal;
	rates.o2_uptake=vmax*oxygen_fraction;
	const double demand=vmax*k.atp_per_o2;
	const double shortfall=demand-rates.o2_uptake*k.atp_per_o2;
	const double glycogen_share=glycogen/(glycogen+k.glycogen_half_saturation_mol_m3);
	const double glucose_share=(1.-glycogen_share)
		*c[OrganGlucose]/(k.km_glucose_mol_m3+c[OrganGlucose]);
	const double fuel=glycogen_share+glucose_share;
	rates.lactate_production=std::min(shortfall,v*k.glycolysis_capacity_mol_m3_s*thermal)*fuel;
	rates.energy_deficit=demand>0.?std::clamp((shortfall-rates.lactate_production)/demand,0.,1.):0.;
	const double glycolysis_glucose=fuel>0.?0.5*rates.lactate_production*glucose_share/fuel:0.;
	const double glycolysis_glycogen=fuel>0.?0.5*rates.lactate_production*glycogen_share/fuel:0.;
	rates.lactate_clearance=v*k.lactate_vmax_mol_m3_s*thermal*oxygen_fraction
		*c[OrganLactate]/(k.lactate_km_mol_m3+c[OrganLactate]);
	const double insulin_effect=c[OrganInsulin]/(c[OrganInsulin]+k.insulin_half_effect_mU_L);
	const double release=v*k.glycogenolysis_vmax_mol_m3_s*thermal*glycogen
		/(glycogen+k.glycogen_half_saturation_mol_m3)*(1.-insulin_effect);
	const double uptake=v*k.glucose_uptake_vmax_mol_m3_s*thermal*insulin_effect
		*c[OrganGlucose]/(k.km_glucose_mol_m3+c[OrganGlucose]);
	const double gluconeogenesis=0.5*k.gluconeogenic_fraction*rates.lactate_clearance;
	rates.glycogen_rate=uptake-release-glycolysis_glycogen;
	rates.reaction[OrganOxygen]=-rates.o2_uptake;
	rates.reaction[OrganGlucose]=release+gluconeogenesis-uptake-glycolysis_glucose;
	rates.reaction[OrganLactate]=rates.lactate_production-rates.lactate_clearance;
	rates.reaction[OrganBufferBase]=rates.lactate_clearance-rates.lactate_production;
	rates.reaction[OrganCarbonDioxide]=k.respiratory_quotient*rates.o2_uptake;
	rates.reaction[OrganInsulin]=-v*k.insulin_clearance_s_inv*thermal*c[OrganInsulin];
	return rates;
}

} // namespace iga

#endif
