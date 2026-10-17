#ifndef IGA_PERFUSED_ORGAN_KINETICS_HPP
#define IGA_PERFUSED_ORGAN_KINETICS_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace iga {

// Species carried by the perfusate. Concentrations are per volume of
// perfusate; reaction rates are per volume of tissue.
enum PerfusedOrganSpecies : std::size_t {
	OrganOxygen=0,      // total O2 content, dissolved + haemoglobin-bound [mol/m^3]
	OrganGlucose=1,     // [mol/m^3] (= mmol/L)
	OrganLactate=2,     // [mol/m^3]
	OrganBicarbonate=3, // [mol/m^3]
	OrganCarbonDioxide=4, // dissolved CO2 [mol/m^3]
	OrganInsulin=5,     // [mU/L]
	OrganSpeciesCount=6
};

using OrganConcentration=std::array<double,OrganSpeciesCount>;

inline const std::array<const char*,OrganSpeciesCount>& PerfusedOrganSpeciesNames()
{
	static const std::array<const char*,OrganSpeciesCount> names{{
		"o2_total_mol_m3","glucose_mol_m3","lactate_mol_m3",
		"bicarbonate_mol_m3","co2_dissolved_mol_m3","insulin_mU_L"}};
	return names;
}

// Perfusate chemistry that does not react in the tissue.
struct PerfusateChemistry
{
	double haemoglobin_g_dl=0.;            // Hb in the perfusate
	double o2_solubility_mol_m3_mmhg=1.34e-3; // plasma at 37 C
	double p50_mmhg=26.8;                  // at pH 7.4, 37 C
	double hill_coefficient=2.7;
	double bohr_coefficient=-0.48;         // d log10(P50) / d pH
	double co2_solubility_mol_m3_mmhg=0.0307;
	double bicarbonate_pk=6.1;
};

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
//  - Lactate production releases one H+ per lactate (HCO3- -> CO2), lactate
//    clearance consumes it; aerobic metabolism makes CO2 at the respiratory
//    quotient. pH follows Henderson-Hasselbalch; no other buffer is modelled.
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
	double po2_mmhg=0.,ph=0.,o2_uptake=0.,lactate_production=0.,lactate_clearance=0.;
	double glycogen_rate=0.;       // d(glycogen)/dt
	double energy_deficit=0.;      // unmet fraction of ATP demand in [0,1]
};

inline double OxygenSaturation(double po2_mmhg,double p50_mmhg,double hill)
{
	if(po2_mmhg<=0.)return 0.;
	const double x=std::pow(po2_mmhg/p50_mmhg,hill);
	return x/(1.+x);
}

inline double PerfusatePh(const OrganConcentration& c,const PerfusateChemistry& chemistry)
{
	const double bicarbonate=std::max(c[OrganBicarbonate],1e-9);
	const double co2=std::max(c[OrganCarbonDioxide],1e-9);
	return chemistry.bicarbonate_pk+std::log10(bicarbonate/co2);
}

inline double HaemoglobinSites(const PerfusateChemistry& chemistry)
{
	// g/dL -> g/L, divided by 16114.5 g per mol of haem, -> mol/m^3.
	return chemistry.haemoglobin_g_dl*10./16114.5*1000.;
}

inline double OxygenContent(double po2_mmhg,double ph,const PerfusateChemistry& chemistry)
{
	const double p50=chemistry.p50_mmhg*std::pow(10.,
		chemistry.bohr_coefficient*(std::clamp(ph,6.5,8.0)-7.4));
	return chemistry.o2_solubility_mol_m3_mmhg*po2_mmhg
		+HaemoglobinSites(chemistry)*OxygenSaturation(po2_mmhg,p50,chemistry.hill_coefficient);
}

// Inverts the monotone content curve by safeguarded bisection/secant.
inline double PartialPressureO2(double content,double ph,const PerfusateChemistry& chemistry)
{
	if(content<=0.)return 0.;
	double low=0.,high=content/chemistry.o2_solubility_mol_m3_mmhg;
	double value=std::min(high,chemistry.p50_mmhg);
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
	const double values[]={k.vmax_o2_mol_m3_s,k.atp_per_o2,k.glycolysis_capacity_mol_m3_s,
		k.lactate_vmax_mol_m3_s,k.glycogenolysis_vmax_mol_m3_s,
		k.glucose_uptake_vmax_mol_m3_s,k.insulin_clearance_s_inv,k.respiratory_quotient,
		chemistry.haemoglobin_g_dl};
	for(const double value:values)
		if(!(value>=0.)||!std::isfinite(value))
			throw std::invalid_argument("hepatic kinetics rates must be finite and nonnegative");
	const double positive[]={k.km_o2_mmhg,k.glycogen_half_saturation_mol_m3,k.km_glucose_mol_m3,
		k.lactate_km_mol_m3,k.insulin_half_effect_mU_L,k.q10,
		chemistry.o2_solubility_mol_m3_mmhg,chemistry.p50_mmhg,chemistry.hill_coefficient,
		chemistry.co2_solubility_mol_m3_mmhg};
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
	rates.ph=PerfusatePh(c,chemistry);
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
	rates.reaction[OrganBicarbonate]=rates.lactate_clearance-rates.lactate_production;
	rates.reaction[OrganCarbonDioxide]=k.respiratory_quotient*rates.o2_uptake
		+rates.lactate_production-rates.lactate_clearance;
	rates.reaction[OrganInsulin]=-v*k.insulin_clearance_s_inv*thermal*c[OrganInsulin];
	return rates;
}

} // namespace iga

#endif
