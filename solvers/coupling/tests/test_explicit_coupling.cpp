#include "ExplicitOneDThreeDCoupling.hpp"

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
	RequireRejected([&row] { iga::ValidateExplicitCouplingHistoryRow(row); });
	row.iteration_count = 1;
	row.downstream_root_area_m2 = 0.0;
	RequireRejected([&row] { iga::ValidateExplicitCouplingHistoryRow(row); });
}
