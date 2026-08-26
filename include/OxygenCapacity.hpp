#ifndef IGA_OXYGEN_CAPACITY_HPP
#define IGA_OXYGEN_CAPACITY_HPP

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace iga {

enum class OxygenTransportState { Dissolved, Total };

inline OxygenTransportState ParseOxygenTransportState(const std::string& value)
{
	if (value == "dissolved_oxygen") return OxygenTransportState::Dissolved;
	if (value == "total_oxygen") return OxygenTransportState::Total;
	throw std::runtime_error(
		"oxygen_state must be 'dissolved_oxygen' or 'total_oxygen'");
}

inline const char* OxygenTransportStateName(OxygenTransportState value)
{
	return value == OxygenTransportState::Total ? "total_oxygen" : "dissolved_oxygen";
}

struct OxygenCapacityParameters {
	double hematocrit_percent = 0.0;
	double oxygen_solubility_ml_dl_mmhg = 0.0031;
	double hemoglobin_g_dl = 15.0;
	double p50_mmhg = 26.8;
	double hill_exponent = 2.7;
	double gas_molar_volume_ml_mmol = 22.4;
};

struct OxygenEquilibriumState {
	double total_oxygen_mol_m3 = 0.0;
	double dissolved_oxygen_mol_m3 = 0.0;
	double bound_oxygen_mol_m3 = 0.0;
	double po2_mmhg = 0.0;
	double saturation = 0.0;
};

inline void ValidateOxygenCapacity(const OxygenCapacityParameters& value)
{
	if (!std::isfinite(value.hematocrit_percent)
		|| value.hematocrit_percent < 0.0 || value.hematocrit_percent > 100.0
		|| !std::isfinite(value.oxygen_solubility_ml_dl_mmhg)
		|| !(value.oxygen_solubility_ml_dl_mmhg > 0.0)
		|| !std::isfinite(value.hemoglobin_g_dl) || value.hemoglobin_g_dl < 0.0
		|| !std::isfinite(value.p50_mmhg) || !(value.p50_mmhg > 0.0)
		|| !std::isfinite(value.hill_exponent) || !(value.hill_exponent > 0.0)
		|| !std::isfinite(value.gas_molar_volume_ml_mmol)
		|| !(value.gas_molar_volume_ml_mmol > 0.0))
		throw std::runtime_error("oxygen-capacity parameters are invalid");
}

inline double OxygenHillSaturation(double po2_mmhg,
	double p50_mmhg, double exponent)
{
	const double pressure = std::max(po2_mmhg, 0.0);
	const double numerator = std::pow(pressure, exponent);
	return numerator > 0.0
		? numerator/(numerator+std::pow(p50_mmhg, exponent)) : 0.0;
}

inline double DissolvedOxygenMolM3(double po2_mmhg,
	const OxygenCapacityParameters& parameters)
{
	ValidateOxygenCapacity(parameters);
	return parameters.oxygen_solubility_ml_dl_mmhg
		*std::max(po2_mmhg, 0.0)*10.0/parameters.gas_molar_volume_ml_mmol;
}

inline double BoundOxygenMolM3(double po2_mmhg,
	const OxygenCapacityParameters& parameters)
{
	ValidateOxygenCapacity(parameters);
	const double saturation = OxygenHillSaturation(po2_mmhg,
		parameters.p50_mmhg, parameters.hill_exponent);
	// 1.34 mL O2/g Hb * g Hb/dL * 10 dL/L / (mL/mmol) = mmol/L = mol/m^3.
	return 1.34*parameters.hemoglobin_g_dl*saturation*10.0
		/parameters.gas_molar_volume_ml_mmol;
}

inline OxygenEquilibriumState OxygenFromPo2(double po2_mmhg,
	const OxygenCapacityParameters& parameters)
{
	OxygenEquilibriumState result;
	result.po2_mmhg = std::max(po2_mmhg, 0.0);
	result.saturation = OxygenHillSaturation(result.po2_mmhg,
		parameters.p50_mmhg, parameters.hill_exponent);
	result.dissolved_oxygen_mol_m3 = DissolvedOxygenMolM3(result.po2_mmhg, parameters);
	result.bound_oxygen_mol_m3 = BoundOxygenMolM3(result.po2_mmhg, parameters);
	result.total_oxygen_mol_m3 = result.dissolved_oxygen_mol_m3
		+result.bound_oxygen_mol_m3;
	return result;
}

inline double Po2FromDissolvedOxygen(double dissolved_mol_m3,
	const OxygenCapacityParameters& parameters)
{
	ValidateOxygenCapacity(parameters);
	if (!std::isfinite(dissolved_mol_m3) || dissolved_mol_m3 < 0.0)
		throw std::runtime_error("dissolved oxygen must be finite and nonnegative");
	return dissolved_mol_m3*parameters.gas_molar_volume_ml_mmol
		/(10.0*parameters.oxygen_solubility_ml_dl_mmhg);
}

inline double Po2FromTotalOxygen(double total_mol_m3,
	const OxygenCapacityParameters& parameters)
{
	ValidateOxygenCapacity(parameters);
	if (!std::isfinite(total_mol_m3) || total_mol_m3 < 0.0)
		throw std::runtime_error("total oxygen must be finite and nonnegative");
	if (total_mol_m3 == 0.0) return 0.0;
	double lower = 0.0;
	double upper = 5000.0;
	while (OxygenFromPo2(upper, parameters).total_oxygen_mol_m3 < total_mol_m3
		&& upper < 1.0e7) upper *= 2.0;
	if (OxygenFromPo2(upper, parameters).total_oxygen_mol_m3 < total_mol_m3)
		throw std::runtime_error("total oxygen exceeds the equilibrium inversion bound");
	for (int iteration = 0; iteration < 100; ++iteration) {
		const double middle = 0.5*(lower+upper);
		if (OxygenFromPo2(middle, parameters).total_oxygen_mol_m3 < total_mol_m3)
			lower = middle;
		else upper = middle;
		if (upper-lower <= 1.0e-10*std::max(1.0, upper)) break;
	}
	return 0.5*(lower+upper);
}

inline OxygenEquilibriumState OxygenFromTransported(double concentration,
	OxygenTransportState state, const OxygenCapacityParameters& parameters)
{
	const double po2 = state == OxygenTransportState::Total
		? Po2FromTotalOxygen(concentration, parameters)
		: Po2FromDissolvedOxygen(concentration, parameters);
	return OxygenFromPo2(po2, parameters);
}

} // namespace iga

#endif
