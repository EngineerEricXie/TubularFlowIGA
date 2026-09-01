#include "ExplicitOneDThreeDCoupling.hpp"
#include "AitkenRelaxation.hpp"
#include "StrongOneDThreeDCoupling.hpp"

#include <cassert>
#include <cmath>
#include <sstream>
#include <stdexcept>

template <class Function>
void RequireRejected(Function&& function)
{
	bool rejected = false;
	try { function(); }
	catch (const std::runtime_error&) { rejected = true; }
	assert(rejected);
}

int main()
{
	iga::ExplicitCouplingScalarPreflight preflight;
	preflight.dt_s = 0.01;
	preflight.steps = 2;
	preflight.density_kg_m3 = 1060.0;
	preflight.dynamic_viscosity_pa_s = 0.0035;
	preflight.normalized_length_m = 1.0;
	preflight.reference_inlet_outward_flow_m3_s = -1.0e-8;
	iga::ValidateExplicitCouplingScalarPreflight(preflight);
	preflight.normalized_length_m = 0.999;
	RequireRejected([&preflight] { iga::ValidateExplicitCouplingScalarPreflight(preflight); });
	preflight.normalized_length_m = 1.0;
	preflight.reference_inlet_outward_flow_m3_s = 0.0;
	RequireRejected([&preflight] { iga::ValidateExplicitCouplingScalarPreflight(preflight); });

	assert(std::abs(iga::ExplicitCouplingNormalizedResidual(2.0, -2.0)) < 1.0e-15);
	assert(std::abs(iga::ExplicitCouplingNormalizedResidual(2.0, -1.0)-0.5) < 1.0e-15);
	assert(std::abs(iga::ExplicitCouplingNormalizedResidual(-3.0, 1.0)+2.0/3.0) < 1.0e-15);
	assert(iga::IsExplicitCouplingPetscOption("-ksp_type"));
	assert(iga::IsExplicitCouplingPetscOption("-ksp_type=gmres"));
	assert(!iga::IsExplicitCouplingPetscOption("--output-dir"));
	assert(iga::ExplicitCouplingPetscOptionConsumesNextValue("-ksp_type", "gmres"));
	assert(!iga::ExplicitCouplingPetscOptionConsumesNextValue("-ksp_type=gmres", "ignored"));
	assert(!iga::ExplicitCouplingPetscOptionConsumesNextValue("-ksp_monitor", "-snes_monitor"));

	iga::ExplicitCouplingHistoryRow row;
	row.upstream_root_area_m2 = 1.0;
	row.upstream_terminal_area_m2 = 1.0;
	row.three_d_inlet_area_m2 = 1.0;
	row.three_d_outlet_area_m2 = 1.0;
	row.downstream_root_area_m2 = 1.0;
	row.downstream_terminal_area_m2 = 1.0;
	iga::ValidateExplicitCouplingHistoryRow(row);
	std::ostringstream output;
	iga::WriteExplicitCouplingHistoryHeader(output);
	iga::WriteExplicitCouplingHistoryRow(output, row);
	assert(output.str().find("three_d_mass_imbalance_m3_s") != std::string::npos);
	assert(output.str().find("net_external_outward_flow_m3_s") != std::string::npos);
	assert(output.str().find("external_pressure_drop_pa") != std::string::npos);
	row.iteration_count = 2;
	row.relaxation_factor = 0.5;
	iga::ValidateExplicitCouplingHistoryRow(row);
	row.iteration_count = 0;
	RequireRejected([&row] { iga::ValidateExplicitCouplingHistoryRow(row); });
	row.iteration_count = 1;
	row.relaxation_factor = 1.0;
	row.downstream_root_area_m2 = 0.0;
	RequireRejected([&row] { iga::ValidateExplicitCouplingHistoryRow(row); });

	iga::StrongCouplingControls controls;
	controls.pressure_reference_pa = 100.0;
	iga::ValidateStrongCouplingControls(controls);
	controls.relaxation_factor = 0.0;
	RequireRejected([&controls] { iga::ValidateStrongCouplingControls(controls); });
	controls.relaxation_factor = 0.5;
	assert(std::abs(iga::StrongCouplingSignedNormalizedPressureResidual(100.0, 125.0, 100.0)-0.2) < 1.0e-15);
	assert(std::abs(iga::StrongCouplingPressureResidual(100.0, 75.0, 100.0)-0.25) < 1.0e-15);
	assert(std::abs(iga::StrongCouplingFixedUpdate(100.0, 140.0, 0.5)-120.0) < 1.0e-15);
	RequireRejected([] { iga::StrongCouplingFixedUpdate(1.0, 2.0, 1.1); });
	iga::StrongCouplingIterationRow iteration;
	iteration.physical_step = 1;
	iteration.time_s = 0.01;
	iteration.iteration = 2;
	iteration.applied_upstream_terminal_pressure_pa = 100.0;
	iteration.applied_three_d_outlet_traction_pressure_pa = 90.0;
	iteration.measured_three_d_inlet_pressure_pa = 110.0;
	iteration.measured_downstream_root_pressure_pa = 95.0;
	iteration.signed_upstream_pressure_residual_pa = 10.0;
	iteration.signed_downstream_pressure_residual_pa = 5.0;
	iteration.normalized_upstream_pressure_residual = 1.0/11.0;
	iteration.normalized_downstream_pressure_residual = 0.05;
	iteration.next_upstream_terminal_pressure_pa = 105.0;
	iteration.next_three_d_outlet_traction_pressure_pa = 92.5;
	iteration.upstream_terminal_outward_flow_m3_s = 1.0e-3;
	iteration.three_d_inlet_outward_flow_m3_s = -1.0e-3;
	iteration.three_d_outlet_outward_flow_m3_s = 1.0e-3;
	iteration.downstream_root_outward_flow_m3_s = -1.0e-3;
	iteration.upstream_three_d_flow_residual_m3_s = 1.0e-12;
	iteration.three_d_downstream_flow_residual_m3_s = -1.0e-12;
	iteration.normalized_upstream_three_d_flow_residual = 1.0e-12;
	iteration.normalized_three_d_downstream_flow_residual = -1.0e-12;
	iteration.three_d_wall_outward_flow_m3_s = 0.0;
	iteration.three_d_mass_imbalance_m3_s = 1.0e-12;
	iteration.three_d_attempt_linear_iterations = 3;
	iteration.three_d_cumulative_step_linear_iterations = 6;
	iga::ValidateStrongCouplingIterationRow(iteration);
	std::ostringstream iteration_output;
	iga::WriteStrongCouplingIterationHeader(iteration_output);
	iga::WriteStrongCouplingIterationRow(iteration_output, iteration);
	assert(iteration_output.str().find("signed_upstream_pressure_residual_pa") != std::string::npos);
	iteration.iteration = 0;
	RequireRejected([&iteration] { iga::ValidateStrongCouplingIterationRow(iteration); });
	iteration.iteration = 2;
	iteration.aitken_status_code = 1;
	iteration.aitken_has_previous_residual = false;
	RequireRejected([&iteration] { iga::ValidateStrongCouplingIterationRow(iteration); });
	iteration.aitken_status_code = -1;
	iteration.aitken_scaled_numerator = 1.0;
	RequireRejected([&iteration] { iga::ValidateStrongCouplingIterationRow(iteration); });
	iga::AitkenRelaxationControls aitken_controls;
	iga::ValidateAitkenRelaxationControls(aitken_controls);
	iga::AitkenRelaxation<2> aitken(aitken_controls);
	const std::array<double, 2> x{{0.0, 0.0}};
	const std::array<double, 2> first_residual{{2.0, 0.0}};
	const auto first = aitken.Propose(x, first_residual, 1.0);
	assert(first.status == iga::AitkenRelaxationStatus::Initial && std::abs(first.relaxation_factor-0.5) < 1.0e-15);
	aitken.AcceptApplied(first_residual, first.relaxation_factor);
	const auto dynamic = aitken.Propose(x, std::array<double, 2>{{1.0, 0.0}}, 1.0);
	assert(dynamic.status == iga::AitkenRelaxationStatus::Dynamic && std::abs(dynamic.relaxation_factor-1.0) < 1.0e-15);
	aitken.Reset();
	assert(!aitken.Propose(x, first_residual, 1.0).has_previous_residual);
	aitken.AcceptApplied(first_residual, 0.5);
	const auto equal = aitken.Propose(x, first_residual, 1.0);
	assert(equal.status == iga::AitkenRelaxationStatus::TinyDifferenceFallback);
	iga::AitkenRelaxationControls minimum_controls;
	minimum_controls.minimum_relaxation = 0.2;
	minimum_controls.initial_relaxation = 0.5;
	iga::AitkenRelaxation<2> minimum_aitken(minimum_controls);
	minimum_aitken.AcceptApplied(std::array<double, 2>{{1.0, 0.0}}, 0.5);
	const auto minimum = minimum_aitken.Propose(x, std::array<double, 2>{{2.0, 0.0}}, 1.0);
	assert(minimum.status == iga::AitkenRelaxationStatus::ClampedMinimum && minimum.relaxation_factor == 0.2);
	iga::AitkenRelaxationControls maximum_controls;
	maximum_controls.maximum_relaxation = 0.75;
	iga::AitkenRelaxation<2> maximum_aitken(maximum_controls);
	maximum_aitken.AcceptApplied(std::array<double, 2>{{2.0, 0.0}}, 0.5);
	const auto maximum = maximum_aitken.Propose(x, std::array<double, 2>{{1.0, 0.0}}, 1.0);
	assert(maximum.status == iga::AitkenRelaxationStatus::ClampedMaximum && maximum.relaxation_factor == 0.75);
	iga::AitkenRelaxation<2> default_maximum(iga::AitkenRelaxationControls{});
	default_maximum.AcceptApplied(std::array<double, 2>{{2.0, 0.0}}, 0.5);
	const auto default_clamp = default_maximum.Propose(x, std::array<double, 2>{{1.5, 0.0}}, 1.0);
	assert(default_clamp.status == iga::AitkenRelaxationStatus::ClampedMaximum && default_clamp.relaxation_factor == 1.0);
	const auto huge = aitken.Propose(std::array<double, 2>{{1.0e200, -1.0e200}},
		std::array<double, 2>{{1.0e100, -1.0e100}}, 1.0);
	assert(std::isfinite(huge.next[0]) && std::isfinite(huge.next[1]));
	RequireRejected([&aitken] { aitken.Propose(std::array<double, 2>{{NAN, 0.0}}, std::array<double, 2>{{1.0, 0.0}}, 1.0); });
	iga::AitkenRelaxationControls invalid_controls;
	invalid_controls.minimum_relaxation = 0.8;
	invalid_controls.initial_relaxation = 0.5;
	RequireRejected([&invalid_controls] { iga::ValidateAitkenRelaxationControls(invalid_controls); });
}
