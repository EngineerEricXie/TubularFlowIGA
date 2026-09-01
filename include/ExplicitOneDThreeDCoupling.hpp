#ifndef IGA_EXPLICIT_ONE_D_THREE_D_COUPLING_HPP
#define IGA_EXPLICIT_ONE_D_THREE_D_COUPLING_HPP

#include "CouplingPort.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>

namespace iga {

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
	if (row.iteration_count != 1 || row.relaxation_factor != 1.0
		|| row.three_d_trial_linear_iterations < 0)
		throw std::runtime_error("explicit coupling history requires one unrelaxed staggered iteration");
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
		"iteration_count,relaxation_factor,three_d_trial_linear_iterations\n";
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
		<< row.relaxation_factor << ',' << row.three_d_trial_linear_iterations << '\n';
}

} // namespace iga

#endif
