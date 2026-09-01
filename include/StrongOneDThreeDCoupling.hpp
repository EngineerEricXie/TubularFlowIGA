#ifndef IGA_STRONG_ONE_D_THREE_D_COUPLING_HPP
#define IGA_STRONG_ONE_D_THREE_D_COUPLING_HPP

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <ostream>
#include <stdexcept>

namespace iga {

struct StrongCouplingControls {
	int maximum_iterations = 50;
	double pressure_relative_tolerance = 1.0e-6;
	double pressure_reference_pa = 0.0;
	double flow_relative_tolerance = 1.0e-10;
	double relaxation_factor = 0.5;
};

inline void ValidateStrongCouplingControls(const StrongCouplingControls& controls)
{
	if (controls.maximum_iterations < 1 || !(controls.pressure_relative_tolerance > 0.0)
		|| !std::isfinite(controls.pressure_relative_tolerance)
		|| !(controls.pressure_reference_pa > 0.0) || !std::isfinite(controls.pressure_reference_pa)
		|| !(controls.flow_relative_tolerance > 0.0) || !std::isfinite(controls.flow_relative_tolerance)
		|| !(controls.relaxation_factor > 0.0) || !(controls.relaxation_factor <= 1.0)
		|| !std::isfinite(controls.relaxation_factor))
		throw std::runtime_error("strong coupling controls require positive finite tolerances, pressure reference, and relaxation in (0,1]");
}

inline double StrongCouplingSignedNormalizedPressureResidual(double applied_pa, double measured_pa,
	double reference_pa)
{
	if (!std::isfinite(applied_pa) || !std::isfinite(measured_pa)
		|| !(reference_pa > 0.0) || !std::isfinite(reference_pa))
		throw std::runtime_error("strong coupling pressure residual requires finite values and positive reference");
	return (measured_pa-applied_pa)/std::max({reference_pa, std::abs(measured_pa), std::abs(applied_pa)});
}

inline double StrongCouplingPressureResidual(double applied_pa, double measured_pa,
	double reference_pa)
{
	return std::abs(StrongCouplingSignedNormalizedPressureResidual(applied_pa, measured_pa, reference_pa));
}

inline double StrongCouplingFixedUpdate(double applied_pa, double measured_pa, double relaxation)
{
	if (!std::isfinite(applied_pa) || !std::isfinite(measured_pa)
		|| !(relaxation > 0.0) || !(relaxation <= 1.0) || !std::isfinite(relaxation))
		throw std::runtime_error("strong coupling fixed update requires finite pressures and relaxation in (0,1]");
	return applied_pa+relaxation*(measured_pa-applied_pa);
}

struct StrongCouplingIterationRow {
	int physical_step = 0;
	double time_s = 0.0;
	int iteration = 0;
	double applied_upstream_terminal_pressure_pa = 0.0;
	double applied_three_d_outlet_traction_pressure_pa = 0.0;
	double measured_three_d_inlet_pressure_pa = 0.0;
	double measured_downstream_root_pressure_pa = 0.0;
	double signed_upstream_pressure_residual_pa = 0.0;
	double signed_downstream_pressure_residual_pa = 0.0;
	double normalized_upstream_pressure_residual = 0.0;
	double normalized_downstream_pressure_residual = 0.0;
	double next_upstream_terminal_pressure_pa = 0.0;
	double next_three_d_outlet_traction_pressure_pa = 0.0;
	double upstream_terminal_outward_flow_m3_s = 0.0;
	double three_d_inlet_outward_flow_m3_s = 0.0;
	double three_d_outlet_outward_flow_m3_s = 0.0;
	double downstream_root_outward_flow_m3_s = 0.0;
	double upstream_three_d_flow_residual_m3_s = 0.0;
	double three_d_downstream_flow_residual_m3_s = 0.0;
	double normalized_upstream_three_d_flow_residual = 0.0;
	double normalized_three_d_downstream_flow_residual = 0.0;
	double three_d_wall_outward_flow_m3_s = 0.0;
	double three_d_mass_imbalance_m3_s = 0.0;
	long long three_d_attempt_linear_iterations = 0;
	long long three_d_cumulative_step_linear_iterations = 0;
	int upstream_1d_attempt_configured_substeps = 0;
	long long upstream_1d_attempt_explicit_cfl_substeps = 0;
	long long upstream_1d_cumulative_configured_substeps = 0;
	long long upstream_1d_cumulative_explicit_cfl_substeps = 0;
	int downstream_1d_attempt_configured_substeps = 0;
	long long downstream_1d_attempt_explicit_cfl_substeps = 0;
	long long downstream_1d_cumulative_configured_substeps = 0;
	long long downstream_1d_cumulative_explicit_cfl_substeps = 0;
	double relaxation_factor_for_next_guess = 1.0;
	double unclamped_relaxation_factor = 1.0;
	double aitken_scaled_numerator = 0.0;
	double aitken_scaled_denominator = 0.0;
	int aitken_status_code = -1;
	bool relaxation_update_applied = false;
	bool aitken_has_previous_residual = false;
	bool converged = false;
};

inline void ValidateStrongCouplingIterationRow(const StrongCouplingIterationRow& row)
{
	const double values[] = {row.time_s, row.applied_upstream_terminal_pressure_pa,
		row.applied_three_d_outlet_traction_pressure_pa, row.measured_three_d_inlet_pressure_pa,
		row.measured_downstream_root_pressure_pa, row.signed_upstream_pressure_residual_pa,
		row.signed_downstream_pressure_residual_pa, row.normalized_upstream_pressure_residual,
		row.normalized_downstream_pressure_residual, row.next_upstream_terminal_pressure_pa,
		row.next_three_d_outlet_traction_pressure_pa, row.upstream_three_d_flow_residual_m3_s,
		row.three_d_inlet_outward_flow_m3_s, row.upstream_terminal_outward_flow_m3_s,
		row.three_d_outlet_outward_flow_m3_s, row.downstream_root_outward_flow_m3_s,
		row.three_d_downstream_flow_residual_m3_s,
		row.normalized_upstream_three_d_flow_residual,
		row.normalized_three_d_downstream_flow_residual,
		row.three_d_wall_outward_flow_m3_s, row.three_d_mass_imbalance_m3_s,
		row.relaxation_factor_for_next_guess, row.unclamped_relaxation_factor,
		row.aitken_scaled_numerator, row.aitken_scaled_denominator};
	for (const auto value : values) if (!std::isfinite(value))
		throw std::runtime_error("strong coupling iteration row requires finite values");
	if (row.physical_step < 1 || row.iteration < 1
		|| row.normalized_upstream_pressure_residual < 0.0
		|| row.normalized_downstream_pressure_residual < 0.0
		|| row.three_d_attempt_linear_iterations < 0
		|| row.three_d_cumulative_step_linear_iterations < row.three_d_attempt_linear_iterations
		|| row.upstream_1d_attempt_configured_substeps < 1 || row.downstream_1d_attempt_configured_substeps < 1
		|| row.upstream_1d_attempt_explicit_cfl_substeps < 0 || row.downstream_1d_attempt_explicit_cfl_substeps < 0
		|| row.upstream_1d_cumulative_configured_substeps < row.upstream_1d_attempt_configured_substeps
		|| row.downstream_1d_cumulative_configured_substeps < row.downstream_1d_attempt_configured_substeps
		|| row.upstream_1d_cumulative_explicit_cfl_substeps < row.upstream_1d_attempt_explicit_cfl_substeps
		|| row.downstream_1d_cumulative_explicit_cfl_substeps < row.downstream_1d_attempt_explicit_cfl_substeps
		|| !(row.relaxation_factor_for_next_guess > 0.0) || row.relaxation_factor_for_next_guess > 1.0
		|| row.aitken_scaled_denominator < 0.0 || row.aitken_status_code < -1 || row.aitken_status_code > 5
		|| (row.converged && row.relaxation_update_applied))
		throw std::runtime_error("strong coupling iteration row has invalid iteration or residual");
	const bool no_previous_status = row.aitken_status_code == -1 || row.aitken_status_code == 0;
	if ((no_previous_status && (row.aitken_has_previous_residual || row.aitken_scaled_numerator != 0.0
		|| row.aitken_scaled_denominator != 0.0
		|| row.relaxation_factor_for_next_guess != row.unclamped_relaxation_factor))
		|| (!no_previous_status && !row.aitken_has_previous_residual)
		|| ((row.aitken_status_code == 4 || row.aitken_status_code == 5)
			&& row.relaxation_factor_for_next_guess != row.unclamped_relaxation_factor))
		throw std::runtime_error("strong coupling iteration row has inconsistent Aitken status");
}

inline void WriteStrongCouplingIterationHeader(std::ostream& output)
{
	output << "physical_step,time_s,iteration,applied_upstream_terminal_pressure_pa,applied_three_d_outlet_traction_pressure_pa,"
		"measured_three_d_inlet_pressure_pa,measured_downstream_root_pressure_pa,"
		"signed_upstream_pressure_residual_pa,signed_downstream_pressure_residual_pa,"
		"normalized_upstream_pressure_residual,normalized_downstream_pressure_residual,"
		"next_upstream_terminal_pressure_pa,next_three_d_outlet_traction_pressure_pa,"
		"upstream_terminal_outward_flow_m3_s,three_d_inlet_outward_flow_m3_s,"
		"three_d_outlet_outward_flow_m3_s,downstream_root_outward_flow_m3_s,"
		"upstream_three_d_flow_residual_m3_s,three_d_downstream_flow_residual_m3_s,"
		"normalized_upstream_three_d_flow_residual,normalized_three_d_downstream_flow_residual,"
		"three_d_wall_outward_flow_m3_s,three_d_mass_imbalance_m3_s,"
		"three_d_attempt_linear_iterations,three_d_cumulative_step_linear_iterations,"
		"upstream_1d_attempt_configured_substeps,upstream_1d_attempt_explicit_cfl_substeps,upstream_1d_cumulative_configured_substeps,upstream_1d_cumulative_explicit_cfl_substeps,"
		"downstream_1d_attempt_configured_substeps,downstream_1d_attempt_explicit_cfl_substeps,downstream_1d_cumulative_configured_substeps,downstream_1d_cumulative_explicit_cfl_substeps,"
		"relaxation_factor_for_next_guess,unclamped_relaxation_factor,aitken_scaled_numerator,aitken_scaled_denominator,"
		"aitken_status_code,relaxation_update_applied,aitken_has_previous_residual,converged\n";
}

inline void WriteStrongCouplingIterationRow(std::ostream& output, const StrongCouplingIterationRow& row)
{
	ValidateStrongCouplingIterationRow(row);
	output << std::setprecision(17) << row.physical_step << ',' << row.time_s << ',' << row.iteration << ','
		<< row.applied_upstream_terminal_pressure_pa << ',' << row.applied_three_d_outlet_traction_pressure_pa << ','
		<< row.measured_three_d_inlet_pressure_pa << ',' << row.measured_downstream_root_pressure_pa << ','
		<< row.signed_upstream_pressure_residual_pa << ',' << row.signed_downstream_pressure_residual_pa << ','
		<< row.normalized_upstream_pressure_residual << ',' << row.normalized_downstream_pressure_residual << ','
		<< row.next_upstream_terminal_pressure_pa << ',' << row.next_three_d_outlet_traction_pressure_pa << ','
		<< row.upstream_terminal_outward_flow_m3_s << ',' << row.three_d_inlet_outward_flow_m3_s << ','
		<< row.three_d_outlet_outward_flow_m3_s << ',' << row.downstream_root_outward_flow_m3_s << ','
		<< row.upstream_three_d_flow_residual_m3_s << ',' << row.three_d_downstream_flow_residual_m3_s << ','
		<< row.normalized_upstream_three_d_flow_residual << ','
		<< row.normalized_three_d_downstream_flow_residual << ','
		<< row.three_d_wall_outward_flow_m3_s << ',' << row.three_d_mass_imbalance_m3_s << ','
		<< row.three_d_attempt_linear_iterations << ',' << row.three_d_cumulative_step_linear_iterations << ','
		<< row.upstream_1d_attempt_configured_substeps << ',' << row.upstream_1d_attempt_explicit_cfl_substeps << ',' << row.upstream_1d_cumulative_configured_substeps << ',' << row.upstream_1d_cumulative_explicit_cfl_substeps << ','
		<< row.downstream_1d_attempt_configured_substeps << ',' << row.downstream_1d_attempt_explicit_cfl_substeps << ',' << row.downstream_1d_cumulative_configured_substeps << ',' << row.downstream_1d_cumulative_explicit_cfl_substeps << ','
		<< row.relaxation_factor_for_next_guess << ',' << row.unclamped_relaxation_factor << ','
		<< row.aitken_scaled_numerator << ',' << row.aitken_scaled_denominator << ',' << row.aitken_status_code << ','
		<< (row.relaxation_update_applied ? 1 : 0) << ',' << (row.aitken_has_previous_residual ? 1 : 0) << ','
		<< (row.converged ? 1 : 0) << '\n';
}

} // namespace iga

#endif
