#include "ImmersedNitscheWall.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

iga::RawSurfaceTriangle Face(std::int64_t a, std::int64_t b, std::int64_t c, std::uint32_t label)
{
	iga::RawSurfaceTriangle face; face.indices = {{a, b, c}}; face.boundary_id = label; return face;
}
iga::RawSurfaceSoup Cube(double lower, double upper, bool mixed_labels = false)
{
	iga::RawSurfaceSoup result;
	result.vertices = {{{{lower,lower,lower}}, {{upper,lower,lower}}, {{upper,upper,lower}}, {{lower,upper,lower}},
		{{lower,lower,upper}}, {{upper,lower,upper}}, {{upper,upper,upper}}, {{lower,upper,upper}}}};
	const auto other = mixed_labels ? 8u : 7u;
	result.triangles = {Face(0,2,1,7), Face(0,3,2,7), Face(4,5,6,other), Face(4,6,7,other),
		Face(0,1,5,7), Face(0,5,4,7), Face(1,2,6,other), Face(1,6,5,other),
		Face(2,3,7,7), Face(2,7,6,7), Face(3,0,4,other), Face(3,4,7,other)};
	return result;
}
iga::CartesianDomainClassification Domain(iga::RawSurfaceSoup soup, iga::CubicCartesianGridSpec spec)
{
	return iga::CartesianDomainClassification(iga::CubicCartesianBackground(spec),
		iga::SurfaceSpatialIndex(iga::ClosedTriangulatedSurface::Build(std::move(soup))));
}
template <class Function> void Reject(Function&& function)
{
	bool rejected = false; try { function(); } catch (const std::exception&) { rejected = true; } assert(rejected);
}
double MaxAbs(const std::vector<PetscScalar>& values)
{
	double result = 0.0; for (const auto value : values) result = std::max(result, std::abs(PetscRealPart(value))); return result;
}
std::vector<PetscScalar> Difference(const std::vector<PetscScalar>& left, const std::vector<PetscScalar>& right)
{
	assert(left.size() == right.size()); std::vector<PetscScalar> result(left.size());
	for (std::size_t i = 0; i < left.size(); ++i) result[i] = left[i]-right[i];
	return result;
}
double WallVelocityFrobeniusNorm(const std::vector<PetscScalar>& matrix, std::size_t nen)
{
	const auto ndof = 4*nen;
	double squared = 0.0;
	for (std::size_t a = 0; a < nen; ++a)
		for (std::size_t b = 0; b < nen; ++b)
			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j) {
					const double value = PetscRealPart(matrix[(4*a+i)*ndof+4*b+j]);
					squared += value*value;
				}
	return std::sqrt(squared);
}
double RelativeError(double actual, double expected)
{
	return std::abs(actual-expected)/std::max({1.0, std::abs(actual), std::abs(expected)});
}
std::vector<std::array<double, 4>> State(const iga::Element& element, double shear = 0.0, double pressure = 0.0)
{
	std::vector<std::array<double, 4>> state(element.connectivity.size());
	for (std::size_t a = 0; a < state.size(); ++a) {
		state[a][0] = shear*element.bezier_points[a][1];
		state[a][3] = pressure;
	}
	return state;
}

} // namespace

int main()
{
	const iga::CubicCartesianGridSpec one{{{0,0,0}}, {{1,1,1}}, {{1,1,1}}};
	const auto domain = Domain(Cube(.25, .75, true), one);
	const iga::CutCellVolumeQuadratureCatalog volume(domain, {5,500000,500000,3000000});
	const iga::ImmersedSurfaceQuadratureCatalog surface(domain);
	const auto element = domain.Background().MaterializeElement(0);
	const iga::NavierStokesParameters parameters{1.0, 0.25, 0.0};
	const auto zero = State(element);

	// The immersed rule is genuinely interior to the Cartesian reference cell;
	// assembly must not route it through a body-fitted face assumption.
	for (const auto& point : surface.UsableRule(domain, 0).Points())
		for (double coordinate : point.parametric) assert(coordinate > 0.0 && coordinate < 1.0);

	const auto without_wall = iga::BuildImmersedNitscheWallElement(domain, volume, surface, 0,
		zero, {}, parameters, {}, 2.0);
	const auto direct_volume = iga::BuildNavierStokesElement(element, zero, {}, parameters,
		volume.UsableRule(domain, 0));
	assert(without_wall.system.jacobian == direct_volume.jacobian);
	assert(without_wall.system.negative_residual == direct_volume.negative_residual);
	// Materializing this exact compact emitted stream is an independent check of
	// the visitor path used by global assembly; it intentionally does not compare
	// against a separately generated expanded catalog.
	const iga::CutCellVolumeQuadratureCatalog compact_volume(domain, {5,500000,500000,3000000},
		iga::CutCellVolumeQuadratureStorageMode::Compact);
	std::vector<iga::VolumeQuadraturePoint> compact_points;
	const auto& compact_rule = compact_volume.UsableCompactRule(domain, 0);
	iga::ForEachVolumePoint(compact_rule, [&](const iga::VolumeQuadraturePoint& point) { compact_points.push_back(point); });
	const auto compact_visitor = iga::BuildNavierStokesElementFromPoints(element, zero, {}, parameters,
		[&compact_rule](const auto& consume) { iga::ForEachVolumePoint(compact_rule, consume); },
		[](const std::array<double, 3>&) { return std::array<double, 3>{{0.0,0.0,0.0}}; });
	const auto compact_materialized = iga::BuildNavierStokesElement(element, zero, {}, parameters,
		iga::VolumeQuadratureRule(compact_points));
	assert(compact_visitor.jacobian == compact_materialized.jacobian);
	assert(compact_visitor.negative_residual == compact_materialized.negative_residual);
	const auto streamed_wall = iga::BuildImmersedNitscheWallElementFromVolumeSystem(domain, compact_volume,
		surface, 0, zero, {}, parameters, {7,8}, compact_visitor, 2.0);
	assert(streamed_wall.system.jacobian.size() == compact_visitor.jacobian.size());
	assert(streamed_wall.system.negative_residual.size() == compact_visitor.negative_residual.size());
	iga::NavierStokesSystem malformed = compact_visitor;
	malformed.negative_residual.pop_back();
	Reject([&] { iga::BuildImmersedNitscheWallElementFromVolumeSystem(domain, compact_volume, surface, 0,
		zero, {}, parameters, {7}, malformed); });
	malformed = compact_visitor; malformed.jacobian.pop_back();
	Reject([&] { iga::BuildImmersedNitscheWallElementFromVolumeSystem(domain, compact_volume, surface, 0,
		zero, {}, parameters, {7}, malformed); });
	malformed = compact_visitor; malformed.negative_residual.push_back(0.0);
	Reject([&] { iga::BuildImmersedNitscheWallElementFromVolumeSystem(domain, compact_volume, surface, 0,
		zero, {}, parameters, {7}, malformed); });
	malformed = compact_visitor; malformed.jacobian[0] = std::numeric_limits<double>::infinity();
	Reject([&] { iga::BuildImmersedNitscheWallElementFromVolumeSystem(domain, compact_volume, surface, 0,
		zero, {}, parameters, {7}, malformed); });
	for (const auto bad : {iga::NavierStokesParameters{0.0, .25, 0.0},
		iga::NavierStokesParameters{1.0, std::numeric_limits<double>::infinity(), 0.0},
		iga::NavierStokesParameters{1.0, .25, -1.0},
		iga::NavierStokesParameters{1.0, .25, std::numeric_limits<double>::quiet_NaN()}})
		Reject([&] { iga::BuildImmersedNitscheWallElementFromVolumeSystem(domain, compact_volume, surface, 0,
			zero, {}, bad, {7}, compact_visitor); });
	const auto rigid = iga::BuildImmersedNitscheWallElement(domain, volume, surface, 0,
		zero, {}, parameters, {7,8}, 2.0);
	assert(MaxAbs(Difference(rigid.system.negative_residual, without_wall.system.negative_residual)) < 2e-13);
	const auto wall_jacobian = Difference(rigid.system.jacobian, without_wall.system.jacobian);
	const auto ndof = 4*element.connectivity.size();
	for (std::size_t a = 0; a < element.connectivity.size(); ++a)
		for (std::size_t b = 0; b < element.connectivity.size(); ++b)
			for (int i = 0; i < 3; ++i) {
				for (int j = 0; j < 3; ++j)
					assert(std::abs(PetscRealPart(wall_jacobian[(4*a+i)*ndof+4*b+j]
						-wall_jacobian[(4*b+j)*ndof+4*a+i])) < 2e-12);
				assert(std::abs(PetscRealPart(wall_jacobian[(4*a+i)*ndof+4*b+3]
					+wall_jacobian[(4*b+3)*ndof+4*a+i])) < 2e-12);
			}
	const auto rigid_high_gamma = iga::BuildImmersedNitscheWallElement(domain, volume, surface, 0,
		zero, {}, parameters, {7,8}, 4.0);
	const auto high_wall_jacobian = Difference(rigid_high_gamma.system.jacobian, without_wall.system.jacobian);
	auto Energy = [ndof](const std::vector<PetscScalar>& matrix, const std::vector<double>& velocity_trial) {
		double value = 0.0;
		for (std::size_t row = 0; row < ndof; ++row)
			for (std::size_t column = 0; column < ndof; ++column)
				value += velocity_trial[row]*PetscRealPart(matrix[row*ndof+column])*velocity_trial[column];
		return value;
	};
	std::vector<std::vector<double>> velocity_trials(3, std::vector<double>(ndof, 0.0));
	for (std::size_t a = 0; a < element.connectivity.size(); ++a) {
		velocity_trials[0][4*a] = 1.0+.01*a;
		velocity_trials[1][4*a+1] = .7-.03*a;
		velocity_trials[2][4*a] = .2+.07*a;
		velocity_trials[2][4*a+1] = -.4+.02*a;
		velocity_trials[2][4*a+2] = .3-.01*a;
	}
	for (const auto& trial : velocity_trials) {
		const double total_energy = Energy(rigid.system.jacobian, trial);
		const double high_total_energy = Energy(rigid_high_gamma.system.jacobian, trial);
		assert(std::isfinite(total_energy) && total_energy >= -3e-11);
		assert(high_total_energy > total_energy+1e-10);
	}
	assert(Energy(wall_jacobian, velocity_trials[0]) >= -2e-12);
	assert(Energy(high_wall_jacobian, velocity_trials[0]) > Energy(wall_jacobian, velocity_trials[0]));
	assert(std::abs(rigid_high_gamma.diagnostics.maximum_eta-2.0*rigid.diagnostics.maximum_eta) < 2e-12);

	// An affine shear with its exact moving wall has zero gap; the volume and
	// consistent traction terms cancel for the closed dyadic cube.
	const auto shear = State(element, 1.5);
	const auto shear_wall = iga::BuildImmersedNitscheWallElement(domain, volume, surface, 0, shear,
		{}, parameters, {7,8}, 2.0, [](const std::array<double, 3>& x, int) {
			return std::array<double, 3>{{1.5*x[1], 0.0, 0.0}};
		});
	assert(MaxAbs(shear_wall.system.negative_residual) < 3e-11);

	// Finite differences apply only to the additive wall contribution: the
	// established full Navier--Stokes tangent is intentionally approximate.
	auto nonzero = State(element, .3, .7);
	for (std::size_t a = 0; a < nonzero.size(); ++a) nonzero[a][1] = .01*(a+1);
	const auto moving = [](const std::array<double, 3>& x, int label) {
		return std::array<double, 3>{{.2+x[0], -.1*x[1], label == 7 ? .3 : -.2}};
	};
	const auto base = iga::BuildImmersedNitscheWallElement(domain, volume, surface, 0, nonzero,
		{}, parameters, {7,8}, 2.0, moving);
	const auto base_volume = iga::BuildImmersedNitscheWallElement(domain, volume, surface, 0, nonzero,
		{}, parameters, {}, 2.0);
	auto direction = State(element);
	for (std::size_t a = 0; a < direction.size(); ++a)
		for (int i = 0; i < 4; ++i) direction[a][i] = .001*(1+a+3*i);
	const double epsilon = 1e-6;
	auto perturbed = nonzero;
	for (std::size_t a = 0; a < perturbed.size(); ++a)
		for (int i = 0; i < 4; ++i) perturbed[a][i] += epsilon*direction[a][i];
	const auto plus = iga::BuildImmersedNitscheWallElement(domain, volume, surface, 0, perturbed,
		{}, parameters, {7,8}, 2.0, moving);
	const auto plus_volume = iga::BuildImmersedNitscheWallElement(domain, volume, surface, 0, perturbed,
		{}, parameters, {}, 2.0);
	const auto base_wall_residual = Difference(base.system.negative_residual, base_volume.system.negative_residual);
	const auto plus_wall_residual = Difference(plus.system.negative_residual, plus_volume.system.negative_residual);
	const auto base_wall_jacobian = Difference(base.system.jacobian, base_volume.system.jacobian);
	for (std::size_t row = 0; row < ndof; ++row) {
		double expected = 0.0;
		for (std::size_t column = 0; column < ndof; ++column)
			expected -= PetscRealPart(base_wall_jacobian[row*ndof+column])*direction[column/4][column%4];
		const double actual = PetscRealPart((plus_wall_residual[row]-base_wall_residual[row])/epsilon);
		assert(std::abs(actual-expected) < 4e-8);
	}
	const auto alternate_moving = iga::BuildImmersedNitscheWallElement(domain, volume, surface, 0, nonzero,
		{}, parameters, {7,8}, 2.0, [](const std::array<double, 3>&, int) {
			return std::array<double, 3>{{3.0, -2.0, 1.0}};
		});
	assert(MaxAbs(Difference(alternate_moving.system.jacobian, base.system.jacobian)) < 2e-13);
	assert(MaxAbs(Difference(alternate_moving.system.negative_residual, base.system.negative_residual)) > 1e-5);

	const auto filtered = iga::BuildImmersedNitscheWallElement(domain, volume, surface, 0, zero,
		{}, parameters, {7}, 2.0);
	assert(filtered.diagnostics.by_boundary_id.at(7).selected_points > 0);
	assert(filtered.diagnostics.by_boundary_id.at(8).skipped_points > 0);
	assert(filtered.diagnostics.fraction_lower > 0.0 && filtered.diagnostics.fraction_estimate > 0.0);
	assert(filtered.diagnostics.minimum_h_n_m > 0.0 && filtered.diagnostics.minimum_eta > 0.0);
	assert(std::abs(filtered.diagnostics.maximum_eta_h_n_over_mu-32.0/filtered.diagnostics.fraction_estimate) < 2e-12);

	const auto pressure = iga::BuildImmersedNitscheWallElement(domain, volume, surface, 0,
		State(element, 0.0, 2.0), {}, parameters, {7,8}, 2.0);
	const auto pressure_volume = iga::BuildImmersedNitscheWallElement(domain, volume, surface, 0,
		State(element, 0.0, 2.0), {}, parameters, {}, 2.0);
	const auto pressure_wall = Difference(pressure.system.negative_residual, pressure_volume.system.negative_residual);
	double manual_pressure = 0.0;
	for (const auto& point : surface.UsableRule(domain, 0).Points()) {
		const auto basis = iga::EvaluateBasis(element, point.parametric[0], point.parametric[1], point.parametric[2]);
		manual_pressure -= 2.0*basis.value[0]*point.normal[0]*point.weight;
	}
	assert(std::abs(PetscRealPart(pressure_wall[0])-manual_pressure) < 2e-12);

	// The pressure test equation receives the signed normal gap directly.
	const auto gap_pressure = iga::BuildImmersedNitscheWallElement(domain, volume, surface, 0,
		zero, {}, parameters, {7,8}, 2.0, [](const std::array<double, 3>& x, int) {
			return std::array<double, 3>{{x[0], 0.0, 0.0}};
		});
	const auto gap_pressure_wall = Difference(gap_pressure.system.negative_residual,
		without_wall.system.negative_residual);
	double manual_pressure_gap = 0.0;
	for (const auto& point : surface.UsableRule(domain, 0).Points()) {
		const auto basis = iga::EvaluateBasis(element, point.parametric[0], point.parametric[1], point.parametric[2]);
		manual_pressure_gap += basis.value[0]*(-point.normal[0]*point.physical[0])*point.weight;
	}
	assert(std::abs(manual_pressure_gap) > 1e-5);
	assert(std::abs(PetscRealPart(gap_pressure_wall[3])-manual_pressure_gap) < 2e-12);

	// Similar immersed cubes in scaled Cartesian elements preserve the cut
	// fraction while h_n, eta, and the isolated velocity wall stiffness have
	// their physical s, 1/s, and s scalings, respectively.
	const auto ScaledWall = [&](double scale) {
		const iga::CubicCartesianGridSpec grid{{{0,0,0}}, {{scale,scale,scale}}, {{1,1,1}}};
		const auto scaled_domain = Domain(Cube(.25*scale, .75*scale, true), grid);
		const iga::CutCellVolumeQuadratureCatalog scaled_volume(scaled_domain, {5,500000,500000,3000000});
		const iga::ImmersedSurfaceQuadratureCatalog scaled_surface(scaled_domain);
		const auto scaled_element = scaled_domain.Background().MaterializeElement(0);
		const auto scaled_zero = State(scaled_element);
		const auto scaled_wall = iga::BuildImmersedNitscheWallElement(scaled_domain, scaled_volume, scaled_surface,
			0, scaled_zero, {}, parameters, {7,8}, 2.0);
		const auto scaled_volume_only = iga::BuildImmersedNitscheWallElement(scaled_domain, scaled_volume, scaled_surface,
			0, scaled_zero, {}, parameters, {}, 2.0);
		return std::array<double, 6>{{scaled_wall.diagnostics.fraction_estimate,
			scaled_wall.diagnostics.minimum_h_n_m, scaled_wall.diagnostics.maximum_h_n_m,
			scaled_wall.diagnostics.minimum_eta, scaled_wall.diagnostics.maximum_eta,
			WallVelocityFrobeniusNorm(Difference(scaled_wall.system.jacobian,
				scaled_volume_only.system.jacobian), scaled_element.connectivity.size())}};
	};
	const auto scale_one = ScaledWall(1.0), scale_two = ScaledWall(2.0), scale_half = ScaledWall(.5);
	for (const auto& scaled : {scale_two, scale_half}) assert(RelativeError(scaled[0], scale_one[0]) < 5e-13);
	for (std::size_t i = 1; i < 3; ++i) {
		assert(RelativeError(scale_two[i], 2.0*scale_one[i]) < 5e-12);
		assert(RelativeError(scale_half[i], .5*scale_one[i]) < 5e-12);
	}
	for (std::size_t i = 3; i < 5; ++i) {
		assert(RelativeError(scale_two[i], .5*scale_one[i]) < 5e-12);
		assert(RelativeError(scale_half[i], 2.0*scale_one[i]) < 5e-12);
	}
	assert(RelativeError(scale_two[5], 2.0*scale_one[5]) < 5e-12);
	assert(RelativeError(scale_half[5], .5*scale_one[5]) < 5e-12);

	// Axis-aligned immersed fragments in an anisotropic element exercise the
	// physical h_n=1/||J^{-1}n|| metric rather than a reference-cell length.
	const iga::CubicCartesianGridSpec anisotropic_grid{{{0,0,0}}, {{2,3,4}}, {{1,1,1}}};
	const auto anisotropic_domain = Domain(Cube(.5, 1.5, true), anisotropic_grid);
	const iga::CutCellVolumeQuadratureCatalog anisotropic_volume(anisotropic_domain, {5,500000,500000,3000000});
	const iga::ImmersedSurfaceQuadratureCatalog anisotropic_surface(anisotropic_domain);
	const auto anisotropic_element = anisotropic_domain.Background().MaterializeElement(0);
	const auto anisotropic_wall = iga::BuildImmersedNitscheWallElement(anisotropic_domain, anisotropic_volume,
		anisotropic_surface, 0, State(anisotropic_element), {}, parameters, {7,8}, 2.0);
	assert(std::abs(anisotropic_wall.diagnostics.minimum_h_n_m-2.0) < 2e-12);
	assert(std::abs(anisotropic_wall.diagnostics.maximum_h_n_m-4.0) < 2e-12);
	assert(RelativeError(anisotropic_wall.diagnostics.maximum_eta,
		8.0/(anisotropic_wall.diagnostics.fraction_estimate*2.0)) < 2e-12);
	assert(RelativeError(anisotropic_wall.diagnostics.minimum_eta,
		8.0/(anisotropic_wall.diagnostics.fraction_estimate*4.0)) < 2e-12);

	Reject([&] { iga::BuildImmersedNitscheWallElement(domain, volume, surface, 0, zero, {}, parameters, {7}, 0.0); });
	Reject([&] { iga::BuildImmersedNitscheWallElement(domain, volume, surface, 0, zero, {}, parameters, {8,7}, 2.0); });
	Reject([&] { iga::BuildImmersedNitscheWallElement(domain, volume, surface, 0, zero, {}, parameters, {9}, 2.0); });
	Reject([&] { iga::BuildImmersedNitscheWallElement(domain, volume, surface, 0, zero, {}, parameters, {7}, 2.0,
		[](const std::array<double, 3>&, int) { return std::array<double, 3>{{std::numeric_limits<double>::infinity(),0,0}}; }); });
	const double largest_finite = std::numeric_limits<double>::max();
	Reject([&] { iga::BuildImmersedNitscheWallElement(domain, volume, surface, 0, zero, {}, parameters, {7}, 2.0,
		[largest_finite](const std::array<double, 3>&, int) { return std::array<double, 3>{{largest_finite,0,0}}; }); });
	auto opposite_largest = zero;
	for (auto& node : opposite_largest) node[0] = -largest_finite;
	Reject([&] { iga::BuildImmersedNitscheWallElement(domain, volume, surface, 0, opposite_largest,
		{}, parameters, {7}, 2.0, [largest_finite](const std::array<double, 3>&, int) {
			return std::array<double, 3>{{largest_finite,0,0}};
		}); });
	Reject([&] { iga::BuildImmersedNitscheWallElement(domain, volume, surface, 0, zero, {}, parameters,
		{7}, largest_finite); });
	const auto wrong_surface_domain = Domain(Cube(.20, .70, true), one);
	const iga::ImmersedSurfaceQuadratureCatalog wrong_surface(wrong_surface_domain);
	Reject([&] { iga::BuildImmersedNitscheWallElement(domain, volume, wrong_surface, 0, zero, {}, parameters, {7}, 2.0); });
	const auto wrong_domain = Domain(Cube(.25, .75), {{{0,0,0}}, {{2,1,1}}, {{1,1,1}}});
	Reject([&] { iga::BuildImmersedNitscheWallElement(wrong_domain, volume, surface, 0, zero, {}, parameters, {7}, 2.0); });
	return 0;
}
