#include "WomersleyReference.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>

int main()
{
	using Complex = std::complex<long double>;
	assert(std::abs(iga::ComplexBesselJ0(Complex(0.0L, 0.0L))-Complex(1.0L, 0.0L))
		< 1e-18L);
	for (const long double value : {1.0L, 5.0L, 10.0L})
		assert(std::abs(std::real(iga::ComplexBesselJ0(Complex(value, 0.0L)))
			-std::cyl_bessel_j(0.0L, value)) < 1e-14L);
	const auto moderate = iga::ComplexBesselJ0(Complex(10.0L, -10.0L));
	const auto moderate_reference = Complex(-2314.9753144452134L, -411.5628570253805L);
	assert(std::abs(moderate-moderate_reference)/std::abs(moderate_reference) < 1e-14L);
	const auto high = iga::ComplexBesselJ0(Complex(30.0L, -30.0L));
	const auto high_reference = Complex(-155873574616.3195L, -637098819459.4559L);
	assert(std::abs(high-high_reference)/std::abs(high_reference) < 1e-12L);

	const auto parsed = iga::ParseWomersleyReferenceConfiguration(R"({
		"schema_version": 1,
		"axis_origin": [0, 0, 0],
		"axis_direction": [0, 0, 2],
		"radius": 1,
		"density": 1,
		"dynamic_viscosity": 1,
		"period": 2,
		"mean_pressure_gradient": 4,
		"cosine_pressure_gradient": [2],
		"sine_pressure_gradient": [],
		"sample_times": [0, 0.5, 1, 1.5]
	})");
	assert(std::abs(parsed.axis_direction[2]-1.0) < 1e-15);
	assert(parsed.cosine_pressure_gradient.size() == 1);
	assert(parsed.sine_pressure_gradient.size() == 1);
	assert(std::abs(iga::WomersleyAxialVelocity(parsed, 1.0, 0.25)) < 1e-14);
	assert(std::abs(iga::WomersleyAxialVelocity(parsed, 0.4, 0.25)
		-iga::WomersleyAxialVelocity(parsed, 0.4, 2.25)) < 1e-13);
	const auto vector_velocity = iga::WomersleyVelocity(parsed, {{0.5, 0.0, 0.2}}, 0.0);
	assert(vector_velocity[0] == 0.0 && vector_velocity[1] == 0.0);

	auto steady = parsed;
	steady.cosine_pressure_gradient.clear();
	steady.sine_pressure_gradient.clear();
	assert(std::abs(iga::WomersleyAxialVelocity(steady, 0.5, 17.0)-0.75) < 1e-15);

	auto quasi_steady = parsed;
	quasi_steady.mean_pressure_gradient = 0.0;
	quasi_steady.period = 1e7;
	const auto center_velocity = iga::WomersleyAxialVelocity(quasi_steady, 0.0, 0.0);
	assert(std::abs(center_velocity-0.5) < 1e-7);

	std::cout << "Womersley reference tests passed\n";
	return 0;
}
