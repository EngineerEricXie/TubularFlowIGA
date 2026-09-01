#ifndef IGA_EXPLICIT_ONE_D_THREE_D_COUPLING_HPP
#define IGA_EXPLICIT_ONE_D_THREE_D_COUPLING_HPP

#include "CouplingPort.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>

namespace iga {

struct OneDSubcyclingPlan {
	int configured_substeps_per_macro_step = 1;
	double configured_dt_s = 0.0;
	double macro_dt_s = 0.0;
	int macro_steps = 0;
	int configured_steps = 0;
};

inline OneDSubcyclingPlan MakeOneDSubcyclingPlan(double configured_dt_s, int configured_steps,
	double macro_dt_s, int macro_steps)
{
	if (!(configured_dt_s > 0.0) || !std::isfinite(configured_dt_s)
		|| configured_steps < 1 || !(macro_dt_s > 0.0) || !std::isfinite(macro_dt_s) || macro_steps < 1)
		throw std::runtime_error("1D subcycling requires finite positive dt and step counts");
	const long double ratio = static_cast<long double>(macro_dt_s)/static_cast<long double>(configured_dt_s);
	if (!std::isfinite(ratio) || ratio < 1.0L
		|| ratio > static_cast<long double>(std::numeric_limits<int>::max()))
		throw std::runtime_error("1D subcycling requires configured dt no larger than macro dt and an int substep count");
	const long long rounded = std::llround(ratio);
	if (rounded < 1 || rounded > std::numeric_limits<int>::max())
		throw std::runtime_error("1D subcycling ratio is out of range");
	const long double reconstructed = static_cast<long double>(rounded)*static_cast<long double>(configured_dt_s);
	const long double scale = std::max({1.0L, std::abs(static_cast<long double>(macro_dt_s)),
		std::abs(reconstructed)});
	if (std::abs(ratio-static_cast<long double>(rounded)) > 1.0e-12L*std::max(1.0L, std::abs(ratio))
		|| std::abs(reconstructed-static_cast<long double>(macro_dt_s)) > 1.0e-12L*scale)
		throw std::runtime_error("1D subcycling requires an integer macro/configured dt ratio");
	if (static_cast<long long>(macro_steps) > std::numeric_limits<int>::max()/rounded
		|| configured_steps != macro_steps*rounded)
		throw std::runtime_error("1D subcycling configured steps must equal macro steps times substeps");
	const long double configured_horizon = static_cast<long double>(configured_steps)*configured_dt_s;
	const long double macro_horizon = static_cast<long double>(macro_steps)*macro_dt_s;
	if (std::abs(configured_horizon-macro_horizon) > 1.0e-12L*
		std::max({1.0L, std::abs(configured_horizon), std::abs(macro_horizon)}))
		throw std::runtime_error("1D subcycling configured and macro horizons differ");
	return {static_cast<int>(rounded), configured_dt_s, macro_dt_s, macro_steps, configured_steps};
}

// Per-macro-step work is accumulated from each actual 1D trial attempt.  In
// particular, strong iterations must not infer it from a nominal substep plan:
// a failed attempt can complete only part of that plan.
struct OneDTrialWorkAccumulator {
	long long configured_substeps = 0;
	long long explicit_cfl_substeps = 0;

	void Add(int attempted_configured_substeps, long long attempted_explicit_cfl_substeps)
	{
		if (attempted_configured_substeps < 0 || attempted_explicit_cfl_substeps < 0)
			throw std::runtime_error("1D trial work requires nonnegative attempt counters");
		if (configured_substeps > std::numeric_limits<long long>::max()-attempted_configured_substeps
			|| explicit_cfl_substeps > std::numeric_limits<long long>::max()-attempted_explicit_cfl_substeps)
			throw std::runtime_error("1D trial work counter overflow");
		configured_substeps += attempted_configured_substeps;
		explicit_cfl_substeps += attempted_explicit_cfl_substeps;
	}

	long long RejectedConfiguredSubsteps(int accepted_configured_substeps) const
	{
		if (accepted_configured_substeps < 0 || accepted_configured_substeps > configured_substeps)
			throw std::runtime_error("1D accepted configured work exceeds trial work");
		return configured_substeps-accepted_configured_substeps;
	}

	long long RejectedExplicitCflSubsteps(long long accepted_explicit_cfl_substeps) const
	{
		if (accepted_explicit_cfl_substeps < 0 || accepted_explicit_cfl_substeps > explicit_cfl_substeps)
			throw std::runtime_error("1D accepted explicit CFL work exceeds trial work");
		return explicit_cfl_substeps-accepted_explicit_cfl_substeps;
	}
};

struct ExplicitCouplingScalarPreflight {
	double dt_s = 0.0;
	int steps = 0;
	double density_kg_m3 = 0.0;
	double dynamic_viscosity_pa_s = 0.0;
	double normalized_length_m = 0.0;
	double reference_inlet_outward_flow_m3_s = 0.0;
};

// Application options use a double dash.  A single-dash option is passed
// through to PETSc, whose conventional separated-value form is
// "-option value".
inline bool IsExplicitCouplingPetscOption(const std::string& argument)
{
	return argument.size() > 1 && argument[0] == '-' && argument[1] != '-';
}

inline bool ExplicitCouplingPetscOptionConsumesNextValue(const std::string& option,
	const std::string& next)
{
	return IsExplicitCouplingPetscOption(option) && option.find('=') == std::string::npos
		&& !next.empty() && next[0] != '-';
}

inline void ValidateExplicitCouplingScalarPreflight(
	const ExplicitCouplingScalarPreflight& values)
{
	if (!(values.dt_s > 0.0) || !std::isfinite(values.dt_s) || values.steps < 1)
		throw std::runtime_error("explicit 1D--3D coupling requires finite positive dt and steps");
	if (!(values.density_kg_m3 > 0.0) || !std::isfinite(values.density_kg_m3)
		|| !(values.dynamic_viscosity_pa_s > 0.0)
		|| !std::isfinite(values.dynamic_viscosity_pa_s))
		throw std::runtime_error("explicit 1D--3D coupling requires finite positive density and viscosity");
	if (!std::isfinite(values.normalized_length_m)
		|| std::abs(values.normalized_length_m-1.0) > 1.0e-12)
		throw std::runtime_error("explicit 1D--3D coupling requires .ntiga normalized length exactly 1 m");
	if (!std::isfinite(values.reference_inlet_outward_flow_m3_s)
		|| values.reference_inlet_outward_flow_m3_s == 0.0)
		throw std::runtime_error("explicit 1D--3D coupling requires finite nonzero reference inlet flow");
}

inline double ExplicitCouplingNormalizedResidual(double first_outward_flow_m3_s,
	double second_outward_flow_m3_s)
{
	RequireFinitePortValue("explicit coupling first outward flow", first_outward_flow_m3_s);
	RequireFinitePortValue("explicit coupling second outward flow", second_outward_flow_m3_s);
	const double scale = std::max({std::abs(first_outward_flow_m3_s),
		std::abs(second_outward_flow_m3_s), 1.0e-30});
	return (first_outward_flow_m3_s+second_outward_flow_m3_s)/scale;
}

struct ExplicitCouplingHistoryRow {
	double time_s = 0.0;
	double upstream_root_pressure_pa = 0.0;
	double upstream_root_outward_flow_m3_s = 0.0;
	double upstream_root_area_m2 = 0.0;
	double upstream_terminal_pressure_pa = 0.0;
	double upstream_terminal_outward_flow_m3_s = 0.0;
	double upstream_terminal_area_m2 = 0.0;
	double three_d_inlet_pressure_pa = 0.0;
	double three_d_inlet_outward_flow_m3_s = 0.0;
	double three_d_inlet_area_m2 = 0.0;
	double three_d_outlet_pressure_pa = 0.0;
	double three_d_outlet_outward_flow_m3_s = 0.0;
	double three_d_outlet_area_m2 = 0.0;
	double downstream_root_pressure_pa = 0.0;
	double downstream_root_outward_flow_m3_s = 0.0;
	double downstream_root_area_m2 = 0.0;
	double downstream_terminal_pressure_pa = 0.0;
	double downstream_terminal_outward_flow_m3_s = 0.0;
	double downstream_terminal_area_m2 = 0.0;
	double upstream_three_d_flow_residual_m3_s = 0.0;
	double three_d_downstream_flow_residual_m3_s = 0.0;
	double upstream_three_d_normalized_residual = 0.0;
	double three_d_downstream_normalized_residual = 0.0;
	double three_d_mass_imbalance_m3_s = 0.0;
	double three_d_wall_outward_flow_m3_s = 0.0;
	double net_external_outward_flow_m3_s = 0.0;
	double external_pressure_drop_pa = 0.0;
	double upstream_three_d_pressure_jump_pa = 0.0;
	double three_d_downstream_pressure_jump_pa = 0.0;
	int iteration_count = 1;
	double relaxation_factor = 1.0;
	long long three_d_trial_linear_iterations = 0;
	int upstream_configured_substeps_attempted = 0;
	long long upstream_explicit_cfl_substeps = 0;
	int downstream_configured_substeps_attempted = 0;
	long long downstream_explicit_cfl_substeps = 0;
	long long upstream_accepted_configured_substeps = 0;
	long long upstream_all_configured_substeps = 0;
	long long upstream_rejected_configured_substeps = 0;
	long long upstream_accepted_explicit_cfl_substeps = 0;
	long long upstream_all_explicit_cfl_substeps = 0;
	long long upstream_rejected_explicit_cfl_substeps = 0;
	long long downstream_accepted_configured_substeps = 0;
	long long downstream_all_configured_substeps = 0;
	long long downstream_rejected_configured_substeps = 0;
	long long downstream_accepted_explicit_cfl_substeps = 0;
	long long downstream_all_explicit_cfl_substeps = 0;
	long long downstream_rejected_explicit_cfl_substeps = 0;
};

inline void ValidateExplicitCouplingHistoryRow(const ExplicitCouplingHistoryRow& row)
{
	const double* values[] = {&row.time_s, &row.upstream_root_pressure_pa,
		&row.upstream_root_outward_flow_m3_s, &row.upstream_root_area_m2,
		&row.upstream_terminal_pressure_pa,
		&row.upstream_terminal_outward_flow_m3_s, &row.upstream_terminal_area_m2,
		&row.three_d_inlet_pressure_pa, &row.three_d_inlet_outward_flow_m3_s,
		&row.three_d_inlet_area_m2, &row.three_d_outlet_pressure_pa,
		&row.three_d_outlet_outward_flow_m3_s, &row.three_d_outlet_area_m2,
		&row.downstream_root_pressure_pa, &row.downstream_root_outward_flow_m3_s,
		&row.downstream_root_area_m2, &row.downstream_terminal_pressure_pa,
		&row.downstream_terminal_outward_flow_m3_s, &row.downstream_terminal_area_m2,
		&row.upstream_three_d_flow_residual_m3_s,
		&row.three_d_downstream_flow_residual_m3_s,
		&row.upstream_three_d_normalized_residual,
		&row.three_d_downstream_normalized_residual, &row.three_d_mass_imbalance_m3_s,
		&row.three_d_wall_outward_flow_m3_s,
		&row.net_external_outward_flow_m3_s, &row.external_pressure_drop_pa,
		&row.upstream_three_d_pressure_jump_pa,
		&row.three_d_downstream_pressure_jump_pa, &row.relaxation_factor};
	for (const auto* value : values) RequireFinitePortValue("explicit coupling history value", *value);
	if (!(row.upstream_root_area_m2 > 0.0) || !(row.upstream_terminal_area_m2 > 0.0)
		|| !(row.three_d_inlet_area_m2 > 0.0)
		|| !(row.three_d_outlet_area_m2 > 0.0) || !(row.downstream_root_area_m2 > 0.0))
		throw std::runtime_error("explicit coupling history requires positive port areas");
	if (!(row.downstream_terminal_area_m2 > 0.0))
		throw std::runtime_error("explicit coupling history requires positive external terminal area");
	if (row.iteration_count < 1 || !(row.relaxation_factor > 0.0) || row.relaxation_factor > 1.0
		|| row.three_d_trial_linear_iterations < 0)
		throw std::runtime_error("coupling history requires a positive iteration count and relaxation in (0,1]");
	if (row.upstream_configured_substeps_attempted < 1 || row.downstream_configured_substeps_attempted < 1
		|| row.upstream_explicit_cfl_substeps < 0 || row.downstream_explicit_cfl_substeps < 0)
		throw std::runtime_error("coupling history requires nonnegative 1D subcycling diagnostics");
	if (row.upstream_accepted_configured_substeps < 0 || row.upstream_all_configured_substeps < row.upstream_accepted_configured_substeps
		|| row.upstream_rejected_configured_substeps != row.upstream_all_configured_substeps-row.upstream_accepted_configured_substeps
		|| row.downstream_accepted_configured_substeps < 0 || row.downstream_all_configured_substeps < row.downstream_accepted_configured_substeps
		|| row.downstream_rejected_configured_substeps != row.downstream_all_configured_substeps-row.downstream_accepted_configured_substeps
		|| row.upstream_accepted_explicit_cfl_substeps < 0 || row.upstream_all_explicit_cfl_substeps < row.upstream_accepted_explicit_cfl_substeps
		|| row.upstream_rejected_explicit_cfl_substeps != row.upstream_all_explicit_cfl_substeps-row.upstream_accepted_explicit_cfl_substeps
		|| row.downstream_accepted_explicit_cfl_substeps < 0 || row.downstream_all_explicit_cfl_substeps < row.downstream_accepted_explicit_cfl_substeps
		|| row.downstream_rejected_explicit_cfl_substeps != row.downstream_all_explicit_cfl_substeps-row.downstream_accepted_explicit_cfl_substeps)
		throw std::runtime_error("coupling history has inconsistent 1D work accounting");
}

inline void WriteExplicitCouplingHistoryHeader(std::ostream& output)
{
	output << "time_s,upstream_root_pressure_pa,upstream_root_outward_flow_m3_s,upstream_root_area_m2,"
		"upstream_terminal_pressure_pa,upstream_terminal_outward_flow_m3_s,upstream_terminal_area_m2,"
		"three_d_inlet_pressure_pa,three_d_inlet_outward_flow_m3_s,three_d_inlet_area_m2,"
		"three_d_outlet_pressure_pa,three_d_outlet_outward_flow_m3_s,three_d_outlet_area_m2,"
		"downstream_root_pressure_pa,downstream_root_outward_flow_m3_s,downstream_root_area_m2,"
		"downstream_terminal_pressure_pa,downstream_terminal_outward_flow_m3_s,downstream_terminal_area_m2,"
		"upstream_three_d_flow_residual_m3_s,three_d_downstream_flow_residual_m3_s,"
		"upstream_three_d_normalized_residual,three_d_downstream_normalized_residual,"
		"three_d_mass_imbalance_m3_s,three_d_wall_outward_flow_m3_s,net_external_outward_flow_m3_s,external_pressure_drop_pa,"
		"upstream_three_d_pressure_jump_pa,three_d_downstream_pressure_jump_pa,"
		"iteration_count,relaxation_factor,three_d_trial_linear_iterations,"
		"upstream_configured_substeps_attempted,upstream_explicit_cfl_substeps,"
		"downstream_configured_substeps_attempted,downstream_explicit_cfl_substeps,"
		"upstream_accepted_configured_substeps,upstream_all_configured_substeps,upstream_rejected_configured_substeps,"
		"upstream_accepted_explicit_cfl_substeps,upstream_all_explicit_cfl_substeps,upstream_rejected_explicit_cfl_substeps,"
		"downstream_accepted_configured_substeps,downstream_all_configured_substeps,downstream_rejected_configured_substeps,"
		"downstream_accepted_explicit_cfl_substeps,downstream_all_explicit_cfl_substeps,downstream_rejected_explicit_cfl_substeps\n";
}

inline void WriteExplicitCouplingHistoryRow(std::ostream& output,
	const ExplicitCouplingHistoryRow& row)
{
	ValidateExplicitCouplingHistoryRow(row);
	output << std::setprecision(17) << row.time_s << ',' << row.upstream_root_pressure_pa << ','
		<< row.upstream_root_outward_flow_m3_s << ',' << row.upstream_root_area_m2 << ','
		<< row.upstream_terminal_pressure_pa << ','
		<< row.upstream_terminal_outward_flow_m3_s << ',' << row.upstream_terminal_area_m2 << ','
		<< row.three_d_inlet_pressure_pa << ',' << row.three_d_inlet_outward_flow_m3_s << ','
		<< row.three_d_inlet_area_m2 << ',' << row.three_d_outlet_pressure_pa << ','
		<< row.three_d_outlet_outward_flow_m3_s << ',' << row.three_d_outlet_area_m2 << ','
		<< row.downstream_root_pressure_pa << ',' << row.downstream_root_outward_flow_m3_s << ','
		<< row.downstream_root_area_m2 << ',' << row.downstream_terminal_pressure_pa << ','
		<< row.downstream_terminal_outward_flow_m3_s << ',' << row.downstream_terminal_area_m2 << ','
		<< row.upstream_three_d_flow_residual_m3_s << ','
		<< row.three_d_downstream_flow_residual_m3_s << ','
		<< row.upstream_three_d_normalized_residual << ','
		<< row.three_d_downstream_normalized_residual << ',' << row.three_d_mass_imbalance_m3_s << ','
		<< row.three_d_wall_outward_flow_m3_s << ',' << row.net_external_outward_flow_m3_s << ','
		<< row.external_pressure_drop_pa << ','
		<< row.upstream_three_d_pressure_jump_pa << ','
		<< row.three_d_downstream_pressure_jump_pa << ',' << row.iteration_count << ','
		<< row.relaxation_factor << ',' << row.three_d_trial_linear_iterations << ','
		<< row.upstream_configured_substeps_attempted << ',' << row.upstream_explicit_cfl_substeps << ','
		<< row.downstream_configured_substeps_attempted << ',' << row.downstream_explicit_cfl_substeps << ','
		<< row.upstream_accepted_configured_substeps << ',' << row.upstream_all_configured_substeps << ',' << row.upstream_rejected_configured_substeps << ','
		<< row.upstream_accepted_explicit_cfl_substeps << ',' << row.upstream_all_explicit_cfl_substeps << ',' << row.upstream_rejected_explicit_cfl_substeps << ','
		<< row.downstream_accepted_configured_substeps << ',' << row.downstream_all_configured_substeps << ',' << row.downstream_rejected_configured_substeps << ','
		<< row.downstream_accepted_explicit_cfl_substeps << ',' << row.downstream_all_explicit_cfl_substeps << ',' << row.downstream_rejected_explicit_cfl_substeps << '\n';
}

} // namespace iga

#endif
