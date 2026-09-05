#ifndef IGA_FLUID_CAUCHY_STRESS_HPP
#define IGA_FLUID_CAUCHY_STRESS_HPP

// Shared incompressible Newtonian Cauchy-stress kernel.  Surface catalogs use
// the outward normal of the fluid domain, so SigmaN is the traction exerted by
// the exterior on the fluid; fluid-on-structure is its negative.
#include <array>
#include <cmath>
#include <stdexcept>

namespace iga {

inline std::array<double, 3> IncompressibleNewtonianCauchyTraction(
	double pressure_pa, const std::array<std::array<double, 3>, 3>& velocity_gradient_per_s,
	double dynamic_viscosity_pa_s, const std::array<double, 3>& outward_fluid_normal)
{
	if (!std::isfinite(pressure_pa) || !std::isfinite(dynamic_viscosity_pa_s)
		|| !(dynamic_viscosity_pa_s > 0.0))
		throw std::invalid_argument("fluid Cauchy stress pressure or viscosity is invalid");
	std::array<double, 3> result{};
	for (int row = 0; row < 3; ++row) {
		if (!std::isfinite(outward_fluid_normal[row]))
			throw std::invalid_argument("fluid Cauchy stress normal is nonfinite");
		result[row] = -pressure_pa*outward_fluid_normal[row];
		for (int column = 0; column < 3; ++column) {
			if (!std::isfinite(velocity_gradient_per_s[row][column]))
				throw std::invalid_argument("fluid Cauchy stress velocity gradient is nonfinite");
			result[row] += dynamic_viscosity_pa_s
				*(velocity_gradient_per_s[row][column]+velocity_gradient_per_s[column][row])
				*outward_fluid_normal[column];
		}
		if (!std::isfinite(result[row])) throw std::overflow_error("fluid Cauchy traction is nonfinite");
	}
	return result;
}

inline std::array<double, 3> FluidOnStructureCauchyTraction(
	double pressure_pa, const std::array<std::array<double, 3>, 3>& velocity_gradient_per_s,
	double dynamic_viscosity_pa_s, const std::array<double, 3>& outward_fluid_normal)
{
	auto result = IncompressibleNewtonianCauchyTraction(pressure_pa, velocity_gradient_per_s,
		dynamic_viscosity_pa_s, outward_fluid_normal);
	for (double& value : result) value = -value;
	return result;
}

} // namespace iga

#endif
