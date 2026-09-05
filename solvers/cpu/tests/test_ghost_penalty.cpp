#include "CutCellGhostPenalty.hpp"
#include "ImmersedNitscheWall.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace {

iga::RawSurfaceTriangle Face(std::int64_t a, std::int64_t b, std::int64_t c)
{
	iga::RawSurfaceTriangle result;
	result.indices = {{a, b, c}};
	return result;
}

iga::RawSurfaceSoup Box(const std::array<double, 3>& lower,
	const std::array<double, 3>& upper)
{
	iga::RawSurfaceSoup result;
	result.vertices = {{{{lower[0],lower[1],lower[2]}}, {{upper[0],lower[1],lower[2]}},
		{{upper[0],upper[1],lower[2]}}, {{lower[0],upper[1],lower[2]}},
		{{lower[0],lower[1],upper[2]}}, {{upper[0],lower[1],upper[2]}},
		{{upper[0],upper[1],upper[2]}}, {{lower[0],upper[1],upper[2]}}}};
	result.triangles = {Face(0,2,1), Face(0,3,2), Face(4,5,6), Face(4,6,7),
		Face(0,1,5), Face(0,5,4), Face(1,2,6), Face(1,6,5),
		Face(2,3,7), Face(2,7,6), Face(3,0,4), Face(3,4,7)};
	return result;
}

iga::CartesianDomainClassification Domain(const iga::CubicCartesianGridSpec& spec,
	const std::array<double, 3>& lower, const std::array<double, 3>& upper)
{
	return iga::CartesianDomainClassification(iga::CubicCartesianBackground(spec),
		iga::SurfaceSpatialIndex(iga::ClosedTriangulatedSurface::Build(Box(lower, upper))));
}

template <class Function> void Reject(Function&& function)
{
	bool rejected = false;
	try {
		function();
	} catch (const std::exception&) {
		rejected = true;
	}
	assert(rejected);
}

double RelativeError(double actual, double expected)
{
	return std::abs(actual-expected)/std::max({1.0, std::abs(actual), std::abs(expected)});
}

double Norm(const std::vector<PetscScalar>& values)
{
	double squared = 0.0;
	for (const auto value : values) squared += PetscRealPart(value)*PetscRealPart(value);
	return std::sqrt(squared);
}

double BlockNorm(const iga::CutCellGhostPenaltyAssembly& assembly, int field)
{
	const auto n = assembly.connectivity.size();
	const auto ndof = 4*n;
	double squared = 0.0;
	for (std::size_t a = 0; a < n; ++a)
		for (std::size_t b = 0; b < n; ++b) {
			const double value = PetscRealPart(assembly.jacobian[(4*a+field)*ndof+4*b+field]);
			squared += value*value;
		}
	return std::sqrt(squared);
}

void AssertCoefficientDiagnostics(const iga::CutCellGhostPenaltyCatalog& catalog)
{
	const auto& diagnostics = catalog.Diagnostics();
	if (catalog.Faces().empty()) {
		assert(diagnostics.minimum_h_normal_m == 0.0 && diagnostics.maximum_h_normal_m == 0.0);
		assert(diagnostics.minimum_face_area_m2 == 0.0 && diagnostics.maximum_face_area_m2 == 0.0);
		assert(diagnostics.maximum_abs_trace_jump == 0.0);
		assert(diagnostics.maximum_abs_coefficient == 0.0);
		return;
	}
	assert(diagnostics.selected_faces == catalog.Faces().size());
	assert(std::isfinite(diagnostics.minimum_h_normal_m) && diagnostics.minimum_h_normal_m > 0.0);
	assert(std::isfinite(diagnostics.maximum_h_normal_m) && diagnostics.maximum_h_normal_m > 0.0);
	assert(std::isfinite(diagnostics.minimum_face_area_m2) && diagnostics.minimum_face_area_m2 > 0.0);
	assert(std::isfinite(diagnostics.maximum_face_area_m2) && diagnostics.maximum_face_area_m2 > 0.0);
	assert(std::isfinite(diagnostics.maximum_abs_trace_jump) && diagnostics.maximum_abs_trace_jump > 0.0);
	assert(std::isfinite(diagnostics.maximum_abs_coefficient));
	assert(diagnostics.maximum_abs_coefficient > 0.0);
}

double Knot(std::uint32_t index, std::uint32_t cells)
{
	return index < 4 ? 0.0 : (index <= cells+2 ? static_cast<double>(index-3) : static_cast<double>(cells));
}

double CubicControl(std::uint32_t index, std::uint32_t cells, int degree)
{
	const std::array<double, 3> knots{{Knot(index+1, cells), Knot(index+2, cells), Knot(index+3, cells)}};
	if (degree == 0) return 1.0;
	if (degree == 1) return (knots[0]+knots[1]+knots[2])/3.0;
	if (degree == 2) return (knots[0]*knots[1]+knots[0]*knots[2]+knots[1]*knots[2])/3.0;
	return knots[0]*knots[1]*knots[2];
}

std::vector<double> Eigenvalues(std::vector<double> matrix, std::size_t n)
{
	for (std::size_t sweep = 0; sweep < 80*n*n; ++sweep) {
		std::size_t p = 0, q = 0;
		double largest = 0.0;
		for (std::size_t row = 0; row < n; ++row)
			for (std::size_t column = row+1; column < n; ++column)
				if (std::abs(matrix[row*n+column]) > largest) {
					largest = std::abs(matrix[row*n+column]);
					p = row;
					q = column;
				}
		if (largest <= 2e-14) break;
		const double app = matrix[p*n+p], aqq = matrix[q*n+q], apq = matrix[p*n+q];
		const double angle = .5*std::atan2(2.0*apq, aqq-app);
		const double c = std::cos(angle), s = std::sin(angle);
		for (std::size_t k = 0; k < n; ++k) if (k != p && k != q) {
			const double mkp = matrix[k*n+p], mkq = matrix[k*n+q];
			matrix[k*n+p] = matrix[p*n+k] = c*mkp-s*mkq;
			matrix[k*n+q] = matrix[q*n+k] = s*mkp+c*mkq;
		}
		matrix[p*n+p] = c*c*app-2.0*s*c*apq+s*s*aqq;
		matrix[q*n+q] = s*s*app+2.0*s*c*apq+c*c*aqq;
		matrix[p*n+q] = matrix[q*n+p] = 0.0;
	}
	std::vector<double> result(n);
	for (std::size_t i = 0; i < n; ++i) result[i] = matrix[i*n+i];
	std::sort(result.begin(), result.end());
	return result;
}

void AssertAlgebra(const iga::CutCellGhostPenaltyCatalog& catalog, std::size_t face,
	const iga::CartesianDomainClassification& domain, const iga::CutCellVolumeQuadratureCatalog& volume)
{
	std::vector<std::array<double, 4>> state(domain.Background().NodeCount());
	std::vector<std::array<double, 4>> direction(state.size());
	for (std::size_t i = 0; i < state.size(); ++i) {
		state[i] = {{.03*(i+1), -.02*(i+2), .01*(i+3), .04*(i+4)}};
		direction[i] = {{.07*(i+1), -.05*(i+2), .03*(i+3), -.02*(i+4)}};
	}
	const auto assembly = catalog.AssembleFace(face, domain, volume, state, 1.7);
	const auto n = 4*assembly.connectivity.size();
	assert(BlockNorm(assembly, 0) > 0.0 && BlockNorm(assembly, 3) > 0.0);
	std::vector<double> symmetric(n*n);
	for (std::size_t row = 0; row < n; ++row)
		for (std::size_t column = 0; column < n; ++column) {
			const double left = PetscRealPart(assembly.jacobian[row*n+column]);
			const double right = PetscRealPart(assembly.jacobian[column*n+row]);
			assert(std::abs(left-right) <= 5e-12);
			symmetric[row*n+column] = .5*(left+right);
		}
	const auto eigenvalues = Eigenvalues(symmetric, n);
	assert(eigenvalues.back() > 0.0);
	assert(eigenvalues.front() >= -1e-11*eigenvalues.back());

	std::vector<std::array<double, 4>> plus = state, minus = state;
	const double epsilon = 1e-6;
	for (std::size_t i = 0; i < state.size(); ++i)
		for (int component = 0; component < 4; ++component) {
			plus[i][component] += epsilon*direction[i][component];
			minus[i][component] -= epsilon*direction[i][component];
		}
	const auto plus_assembly = catalog.AssembleFace(face, domain, volume, plus, 1.7);
	const auto minus_assembly = catalog.AssembleFace(face, domain, volume, minus, 1.7);
	for (std::size_t row = 0; row < n; ++row) {
		double expected = 0.0;
		for (std::size_t column = 0; column < n; ++column)
			expected -= PetscRealPart(assembly.jacobian[row*n+column])
				*direction[assembly.connectivity[column/4]][column%4];
		const double actual = PetscRealPart((plus_assembly.negative_residual[row]
			-minus_assembly.negative_residual[row])/(2.0*epsilon));
		assert(RelativeError(actual, expected) <= 1e-8);
	}
	std::vector<std::array<double, 4>> constant_pressure(state.size());
	for (auto& value : constant_pressure) value[3] = 1.0;
	assert(Norm(catalog.AssembleFace(face, domain, volume, constant_pressure, 1.7).negative_residual) < 2e-11);
}

std::uint64_t CellId(const iga::CubicCartesianGridSpec& spec, std::uint32_t x,
	std::uint32_t y, std::uint32_t z)
{
	return x+spec.cells[0]*(y+spec.cells[1]*z);
}

std::vector<std::array<double, 4>> ElementState(const iga::Element& element)
{
	return std::vector<std::array<double, 4>>(element.connectivity.size());
}

} // namespace

int main()
{
	const iga::CubicCartesianGridSpec spec{{{0,0,0}}, {{3,3,3}}, {{3,3,3}}};
	const auto domain = Domain(spec, {{.2,.2,.2}}, {{2.25,2.8,2.8}});
	const iga::CutCellVolumeQuadratureCatalog volume(domain, {4,300000,300000,3000000});
	const iga::CutCellGhostPenaltyCatalog catalog(domain, volume);
	assert(catalog.Usable() && !catalog.Faces().empty());
	AssertCoefficientDiagnostics(catalog);
	assert(catalog.Diagnostics().quadrature_points == 16*catalog.Faces().size());
	for (std::size_t i = 1; i < catalog.Faces().size(); ++i) {
		const auto& previous = catalog.Faces()[i-1];
		const auto& current = catalog.Faces()[i];
		assert(std::tie(previous.minus_cell, previous.plus_cell, previous.axis)
			< std::tie(current.minus_cell, current.plus_cell, current.axis));
	}
	for (const auto& face : catalog.Faces()) assert(face.minus_cell < face.plus_cell);
	AssertAlgebra(catalog, 0, domain, volume);

	for (std::size_t face = 0; face < catalog.Faces().size(); ++face) {
		for (std::uint32_t order = 0; order < 3; ++order) {
			const auto trace = catalog.EvaluateFaceJump(face, order, .37, .61, domain, volume);
			for (double value : trace.coefficients)
				assert(std::abs(value) <= 1e-11*std::pow(catalog.Faces()[face].h_normal_m, -static_cast<int>(order)));
		}
		const auto trace = catalog.EvaluateFaceJump(face, 3, .37, .61, domain, volume);
		assert(std::any_of(trace.coefficients.begin(), trace.coefficients.end(),
			[](double value) { return std::abs(value) > 1e-6; }));
	}
	for (int degree = 0; degree <= 3; ++degree) {
		std::vector<std::array<double, 4>> polynomial(domain.Background().NodeCount());
		for (std::uint32_t k = 0; k < spec.cells[2]+3; ++k)
			for (std::uint32_t j = 0; j < spec.cells[1]+3; ++j)
				for (std::uint32_t i = 0; i < spec.cells[0]+3; ++i) {
					const auto node = i+(spec.cells[0]+3)*(j+(spec.cells[1]+3)*k);
					const auto value = CubicControl(i, spec.cells[0], degree);
					polynomial[node] = {{value,value,value,value}};
				}
		for (std::size_t face = 0; face < catalog.Faces().size(); ++face)
			assert(Norm(catalog.AssembleFace(face, domain, volume, polynomial, 1.7).negative_residual) < 2e-10);
	}

	// Independent gamma and viscosity factors in the two physical blocks.
	const auto base = catalog.AssembleFace(0, domain, volume,
		std::vector<std::array<double, 4>>(domain.Background().NodeCount()), 1.0);
	iga::CutCellGhostPenaltyOptions tuned;
	tuned.gamma_u = .04;
	tuned.gamma_p = .03;
	const iga::CutCellGhostPenaltyCatalog tuned_catalog(domain, volume, tuned);
	AssertCoefficientDiagnostics(tuned_catalog);
	const auto tuned_assembly = tuned_catalog.AssembleFace(0, domain, volume,
		std::vector<std::array<double, 4>>(domain.Background().NodeCount()), 3.0);
	assert(RelativeError(BlockNorm(tuned_assembly, 0)/BlockNorm(base, 0), 12.0) < 2e-11);
	assert(RelativeError(BlockNorm(tuned_assembly, 3)/BlockNorm(base, 3), 1.0) < 2e-11);
	const auto ScaledNorms = [](double scale) {
		const iga::CubicCartesianGridSpec grid{{{0,0,0}}, {{3*scale,3*scale,3*scale}}, {{3,3,3}}};
		const auto scaled_domain = Domain(grid, {{.2*scale,.2*scale,.2*scale}},
			{{2.25*scale,2.8*scale,2.8*scale}});
		const iga::CutCellVolumeQuadratureCatalog scaled_volume(scaled_domain, {4,300000,300000,3000000});
		const iga::CutCellGhostPenaltyCatalog scaled_catalog(scaled_domain, scaled_volume);
		AssertCoefficientDiagnostics(scaled_catalog);
		const auto scaled = scaled_catalog.AssembleFace(0, scaled_domain, scaled_volume,
			std::vector<std::array<double, 4>>(scaled_domain.Background().NodeCount()), 1.0);
		return std::array<double, 2>{{BlockNorm(scaled, 0), BlockNorm(scaled, 3)}};
	};
	const auto half = ScaledNorms(.5), one = ScaledNorms(1.0), twice = ScaledNorms(2.0);
	assert(RelativeError(half[0], .5*one[0]) < 3e-11 && RelativeError(twice[0], 2.0*one[0]) < 3e-11);
	assert(RelativeError(half[1], .125*one[1]) < 3e-11 && RelativeError(twice[1], 8.0*one[1]) < 3e-11);
	const iga::CubicCartesianGridSpec anisotropic_spec{{{0,0,0}}, {{2,3,5}}, {{3,3,3}}};
	const auto anisotropic_domain = Domain(anisotropic_spec, {{.2,.3,.5}}, {{1.5,2.8,4.5}});
	const iga::CutCellVolumeQuadratureCatalog anisotropic_volume(anisotropic_domain, {4,300000,300000,3000000});
	const iga::CutCellGhostPenaltyCatalog anisotropic(anisotropic_domain, anisotropic_volume);
	AssertCoefficientDiagnostics(anisotropic);
	const auto anisotropic_assembly = anisotropic.AssembleFace(0, anisotropic_domain, anisotropic_volume,
		std::vector<std::array<double, 4>>(anisotropic_domain.Background().NodeCount()), 1.0);
	const auto& base_face = catalog.Faces()[0];
	const auto& anisotropic_face = anisotropic.Faces()[0];
	assert(RelativeError(BlockNorm(anisotropic_assembly, 0)/BlockNorm(base, 0),
		(anisotropic_face.area_m2/anisotropic_face.h_normal_m)/(base_face.area_m2/base_face.h_normal_m)) < 3e-11);
	assert(RelativeError(BlockNorm(anisotropic_assembly, 3)/BlockNorm(base, 3),
		(anisotropic_face.area_m2*anisotropic_face.h_normal_m)/(base_face.area_m2*base_face.h_normal_m)) < 3e-11);

	const auto permuted_domain = Domain(spec, {{.2,.2,.2}}, {{2.8,2.25,2.8}});
	const iga::CutCellVolumeQuadratureCatalog permuted_volume(permuted_domain, {4,300000,300000,3000000});
	const iga::CutCellGhostPenaltyCatalog permuted(permuted_domain, permuted_volume);
	AssertCoefficientDiagnostics(permuted);
	assert(permuted.Faces().size() == catalog.Faces().size());
	assert(permuted.Diagnostics().selected_by_axis[0] == catalog.Diagnostics().selected_by_axis[1]);
	assert(permuted.Diagnostics().selected_by_axis[1] == catalog.Diagnostics().selected_by_axis[0]);
	assert(permuted.Diagnostics().selected_by_axis[2] == catalog.Diagnostics().selected_by_axis[2]);
	assert(iga::CutCellGhostPenaltyCatalog::SelectFace(true, false));
	assert(iga::CutCellGhostPenaltyCatalog::SelectFace(false, true));
	assert(!iga::CutCellGhostPenaltyCatalog::SelectFace(false, false));

	// A near-full terminal Cut cell must not exceed its certified unit upper
	// bound, or ghost catalog construction would fail closed on a valid rule.
	const auto near_full_domain = Domain(spec, {{0.0,0.0,0.0}}, {{2.999,2.9,2.9}});
	const iga::CutCellVolumeQuadratureCatalog near_full_volume(near_full_domain, {4,300000,300000,3000000});
	const auto near_full_cut = CellId(spec, 2, 1, 1);
	const auto& near_full_diagnostics = near_full_volume.Cell(near_full_cut).diagnostics;
	assert(near_full_domain.Cells()[near_full_cut].classification == iga::CellClassification::Cut);
	assert(near_full_diagnostics.lower_reference_volume <= near_full_diagnostics.estimated_reference_volume);
	assert(near_full_diagnostics.estimated_reference_volume <= near_full_diagnostics.upper_reference_volume);
	assert(near_full_diagnostics.upper_reference_volume == 1.0);
	const iga::CutCellGhostPenaltyCatalog near_full_catalog(near_full_domain, near_full_volume);
	assert(near_full_catalog.Covered(near_full_cut));

	// A one-cell closed surface is a positive but intentionally uncovered Cut cell.
	const auto isolated = Domain(spec, {{.2,.2,.2}}, {{.8,.8,.8}});
	const iga::CutCellVolumeQuadratureCatalog isolated_volume(isolated, {4,300000,300000,3000000});
	const iga::CutCellGhostPenaltyCatalog isolated_catalog(isolated, isolated_volume);
	AssertCoefficientDiagnostics(isolated_catalog);
	const iga::ImmersedSurfaceQuadratureCatalog isolated_surface(isolated);
	assert(isolated_catalog.Faces().empty() && !isolated_catalog.UncoveredCells().empty());
	const iga::NavierStokesParameters parameters{1.0, .25, 0.0};
	const auto isolated_cut = isolated_catalog.UncoveredCells().front();
	assert(!isolated_catalog.Covered(isolated_cut));

	Reject([&] { iga::CutCellGhostPenaltyOptions bad; bad.gamma_u = 0.0;
		iga::CutCellGhostPenaltyCatalog rejected(domain, volume, bad); });
	Reject([&] { iga::CutCellGhostPenaltyOptions bad; bad.gamma_p = std::numeric_limits<double>::infinity();
		iga::CutCellGhostPenaltyCatalog rejected(domain, volume, bad); });
	const iga::CubicCartesianGridSpec subunit_spec{{{0,0,0}}, {{1.5,1.5,1.5}}, {{3,3,3}}};
	const auto subunit_domain = Domain(subunit_spec, {{.1,.1,.1}}, {{1.125,1.4,1.4}});
	const iga::CutCellVolumeQuadratureCatalog subunit_volume(subunit_domain, {4,300000,300000,3000000});
	Reject([&] { iga::CutCellGhostPenaltyOptions bad; bad.gamma_u = std::numeric_limits<double>::denorm_min();
		iga::CutCellGhostPenaltyCatalog rejected(subunit_domain, subunit_volume, bad); });
	Reject([&] { iga::CutCellGhostPenaltyOptions bad; bad.gamma_p = std::numeric_limits<double>::denorm_min();
		iga::CutCellGhostPenaltyCatalog rejected(subunit_domain, subunit_volume, bad); });
	const double large_scale = 10.0;
	const iga::CubicCartesianGridSpec near_overflow_spec{{{0,0,0}}, {{3*large_scale,3*large_scale,3*large_scale}}, {{3,3,3}}};
	const auto near_overflow_domain = Domain(near_overflow_spec, {{.2*large_scale,.2*large_scale,.2*large_scale}},
		{{2.25*large_scale,2.8*large_scale,2.8*large_scale}});
	const iga::CutCellVolumeQuadratureCatalog near_overflow_volume(near_overflow_domain, {4,300000,300000,3000000});
	iga::CutCellGhostPenaltyOptions near_overflow_options;
	near_overflow_options.gamma_u = .25*std::numeric_limits<double>::max()/std::pow(large_scale, 5);
	near_overflow_options.gamma_p = .25*std::numeric_limits<double>::max()/std::pow(large_scale, 7);
	const iga::CutCellGhostPenaltyCatalog near_overflow_catalog(near_overflow_domain, near_overflow_volume, near_overflow_options);
	AssertCoefficientDiagnostics(near_overflow_catalog);
	assert(near_overflow_catalog.Diagnostics().maximum_abs_coefficient > std::numeric_limits<double>::max()/8.0);
	Reject([&] { iga::CutCellGhostPenaltyOptions bad; bad.max_faces = catalog.Faces().size()-1;
		iga::CutCellGhostPenaltyCatalog rejected(domain, volume, bad); });
	Reject([&] { iga::CutCellGhostPenaltyOptions bad; bad.max_quadrature_points = 1;
		iga::CutCellGhostPenaltyCatalog rejected(domain, volume, bad); });
	Reject([&] { iga::CutCellGhostPenaltyOptions bad; bad.max_quadrature_points = 15;
		iga::CutCellGhostPenaltyCatalog rejected(domain, volume, bad); });
	Reject([&] { iga::CutCellGhostPenaltyOptions bad; bad.max_trace_entries = 1;
		iga::CutCellGhostPenaltyCatalog rejected(domain, volume, bad); });
	Reject([&] { catalog.EvaluateFaceJump(catalog.Faces().size(), 3, .2, .3, domain, volume); });
	Reject([&] { catalog.EvaluateFaceJump(0, 4, .2, .3, domain, volume); });
	Reject([&] { catalog.AssembleFace(0, domain, volume,
		std::vector<std::array<double, 4>>(domain.Background().NodeCount()), 0.0); });
	auto nonfinite = std::vector<std::array<double, 4>>(domain.Background().NodeCount());
	nonfinite[0][0] = std::numeric_limits<double>::infinity();
	Reject([&] { catalog.AssembleFace(0, domain, volume, nonfinite, 1.0); });
	const auto same_domain = Domain(spec, {{.2,.2,.2}}, {{2.25,2.8,2.8}});
	const iga::CutCellVolumeQuadratureCatalog other_volume(domain, {4,300000,300000,3000000});
	Reject([&] { catalog.EvaluateFaceJump(0, 3, .2, .3, same_domain, volume); });
	Reject([&] { catalog.EvaluateFaceJump(0, 3, .2, .3, domain, other_volume); });

	// The deepest bounded real sliver catalog that is feasible with this
	// retained-rule octree is m=1..4 at depth four; each plane is dyadic.
	for (int m = 1; m <= 4; ++m) {
		const double fraction = std::ldexp(1.0, -m);
		const auto sliver_domain = Domain(spec, {{.1,.1,.1}}, {{2.0+fraction,2.9,2.9}});
		const iga::CutCellVolumeQuadratureCatalog sliver_volume(sliver_domain, {4,800000,800000,8000000});
		const iga::CutCellGhostPenaltyCatalog sliver_catalog(sliver_domain, sliver_volume);
		const iga::ImmersedSurfaceQuadratureCatalog sliver_surface(sliver_domain);
		const auto cut = CellId(spec, 2, 1, 1);
		const auto anchor = CellId(spec, 1, 1, 1);
		assert(sliver_domain.Cells()[cut].classification == iga::CellClassification::Cut);
		assert(RelativeError(sliver_volume.Cell(cut).diagnostics.estimated_reference_volume, fraction) < 2e-12);
		assert(sliver_catalog.Covered(cut));
		auto selected = std::find_if(sliver_catalog.Faces().begin(), sliver_catalog.Faces().end(),
			[anchor, cut](const iga::CutCellGhostPenaltyFace& face) {
				return face.minus_cell == anchor && face.plus_cell == cut && face.axis == 0;
			});
		assert(selected != sliver_catalog.Faces().end());
		const auto face = static_cast<std::size_t>(selected-sliver_catalog.Faces().begin());
		const auto assembly = sliver_catalog.AssembleFace(face, sliver_domain, sliver_volume,
			std::vector<std::array<double, 4>>(sliver_domain.Background().NodeCount()), 1.0);
		assert(BlockNorm(assembly, 0) > 0.0 && BlockNorm(assembly, 3) > 0.0);
		const auto element = sliver_domain.Background().MaterializeElement(cut);
		const auto covered = iga::BuildImmersedNitscheWallElement(sliver_domain, sliver_volume,
			sliver_surface, cut, ElementState(element), {}, parameters, {0}, sliver_catalog, 2.0);
		const auto legacy = iga::BuildImmersedNitscheWallElement(sliver_domain, sliver_volume,
			sliver_surface, cut, ElementState(element), {}, parameters, {0}, 2.0);
		assert(covered.diagnostics.ghost_covered_policy);
		assert(RelativeError(covered.diagnostics.maximum_eta_h_n_over_mu, 32.0) < 2e-12);
		assert(RelativeError(legacy.diagnostics.maximum_eta_h_n_over_mu, 32.0/fraction) < 2e-12);
	}

	std::cout << "ghost_penalty_tests=passed faces=" << catalog.Faces().size() << '\n';
}
