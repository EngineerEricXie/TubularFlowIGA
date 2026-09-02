#include "ImmersedStaticFlowRuntime.hpp"

#include <petscblaslapack.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void Require(bool condition, const char* expression)
{
	if (!condition) throw std::runtime_error(std::string("closure assertion failed: ")+expression);
}

#undef assert
#define assert(expression) Require((expression), #expression)

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kBoxLower = 0.375;
constexpr double kBoxUpper = 2.625;
constexpr double kLength = kBoxUpper-kBoxLower;
constexpr double kWave = kPi/kLength;
constexpr double kVelocityAmplitude = 1.0e-2;
constexpr double kDynamicViscosity = 1.0;
constexpr double kPressureAmplitude = kDynamicViscosity*kWave*kVelocityAmplitude;
constexpr double kPressureMean = 0.0;

iga::RawSurfaceTriangle Face(std::int64_t a, std::int64_t b, std::int64_t c)
{
	iga::RawSurfaceTriangle result; result.indices = {{a,b,c}}; result.boundary_id = 7; return result;
}

iga::RawSurfaceSoup Box(const std::array<double, 3>& lower, const std::array<double, 3>& upper)
{
	iga::RawSurfaceSoup result;
	result.vertices = {{{{lower[0],lower[1],lower[2]}},{{upper[0],lower[1],lower[2]}},{{upper[0],upper[1],lower[2]}},{{lower[0],upper[1],lower[2]}},
		{{lower[0],lower[1],upper[2]}},{{upper[0],lower[1],upper[2]}},{{upper[0],upper[1],upper[2]}},{{lower[0],upper[1],upper[2]}}}};
	result.triangles = {Face(0,2,1),Face(0,3,2),Face(4,5,6),Face(4,6,7),Face(0,1,5),Face(0,5,4),
		Face(1,2,6),Face(1,6,5),Face(2,3,7),Face(2,7,6),Face(3,0,4),Face(3,4,7)};
	return result;
}

iga::CartesianDomainClassification Domain(const iga::CubicCartesianGridSpec& spec,
	const std::array<double, 3>& lower, const std::array<double, 3>& upper)
{
	return iga::CartesianDomainClassification(iga::CubicCartesianBackground(spec),
		iga::SurfaceSpatialIndex(iga::ClosedTriangulatedSurface::Build(Box(lower, upper))));
}

std::array<double, 3> ExactVelocity(const std::array<double, 3>& x)
{
	const double y = kWave*(x[1]-kBoxLower), z = kWave*(x[2]-kBoxLower);
	return {{kVelocityAmplitude*std::sin(y)*std::sin(z), 0.0, 0.0}};
}

std::array<std::array<double, 3>, 3> ExactVelocityGradient(const std::array<double, 3>& x)
{
	std::array<std::array<double, 3>, 3> gradient{};
	const double y = kWave*(x[1]-kBoxLower), z = kWave*(x[2]-kBoxLower);
	gradient[0] = {{0.0, kVelocityAmplitude*kWave*std::cos(y)*std::sin(z),
		kVelocityAmplitude*kWave*std::sin(y)*std::cos(z)}};
	return gradient;
}

double ExactPressure(const std::array<double, 3>& x)
{
	const double X = kWave*(x[0]-kBoxLower), Y = kWave*(x[1]-kBoxLower), Z = kWave*(x[2]-kBoxLower);
	return kPressureAmplitude*std::sin(X)*std::cos(Y)*std::sin(Z)-kPressureMean;
}

std::array<double, 3> ExactForce(const std::array<double, 3>& x)
{
	const double sx = std::sin(kWave*(x[0]-kBoxLower)), cx = std::cos(kWave*(x[0]-kBoxLower));
	const double sy = std::sin(kWave*(x[1]-kBoxLower)), cy = std::cos(kWave*(x[1]-kBoxLower));
	const double sz = std::sin(kWave*(x[2]-kBoxLower)), cz = std::cos(kWave*(x[2]-kBoxLower));
	const double ux = kVelocityAmplitude*sy*sz;
	return {{2.0*kDynamicViscosity*kWave*kWave*ux + kPressureAmplitude*kWave*cx*cy*sz,
		-kPressureAmplitude*kWave*sx*sy*sz, kPressureAmplitude*kWave*sx*cy*cz}};
}

iga::ImmersedStaticFlowOptions ManufacturedOptions()
{
	iga::ImmersedStaticFlowOptions result;
	result.parameters = {1.0, kDynamicViscosity, 0.0};
	result.wall_labels = {7};
	result.wall_gamma0 = 4.0;
	result.wall_velocity = [](const std::array<double, 3>& x, int label) { assert(label == 7); return ExactVelocity(x); };
	result.body_force = [](const std::array<double, 3>& x) { return ExactForce(x); };
	result.nonlinear_maximum_iterations = 96;
	result.ksp_maximum_iterations = 4000;
	result.ksp_relative_tolerance = 1e-10;
	result.nonlinear_relative_tolerance = 1e-9;
	result.nonlinear_absolute_tolerance = 1e-14;
	result.nonlinear_block_reduction = 1e9;
	// Closure-only preconditioner safeguard; the assembled operator is unshifted.
	result.lu_pivot_shift = 1e-12;
	return result;
}

iga::OctreeCutQuadratureOptions CompactOptions(std::uint32_t depth)
{
	iga::OctreeCutQuadratureOptions result;
	result.max_depth = depth;
	result.max_nodes = 3000000;
	result.max_leaves = 2600000;
	result.max_points = 4000000;
	result.max_records = 3000000;
	result.max_retained_bytes = 512u*1024u*1024u;
	result.max_logical_points = 128u*1024u*1024u;
	return result;
}

struct ErrorNorms {
	double velocity_l2 = 0.0, velocity_h1 = 0.0, pressure_l2 = 0.0, divergence_l2 = 0.0;
	double pressure_mean = 0.0, pressure_error_shift = 0.0, measure = 0.0;
	double exact_velocity_l2 = 0.0, exact_velocity_h1 = 0.0, exact_pressure_l2 = 0.0;
	double numerical_velocity_l2 = 0.0, velocity_dot_exact = 0.0;
};

template <class Callback> void VisitPositive(const iga::CartesianDomainClassification& domain,
	const iga::CutCellVolumeQuadratureCatalog& volume, Callback&& callback)
{
	assert(volume.StorageMode() == iga::CutCellVolumeQuadratureStorageMode::Compact);
	for (std::uint64_t id = 0; id < domain.Cells().size(); ++id) {
		const auto classification = domain.Cells()[static_cast<std::size_t>(id)].classification;
		const auto& cell = volume.Cell(id);
		if ((classification != iga::CellClassification::Inside && classification != iga::CellClassification::Cut)
			|| !cell.usable || iga::CompactCutCellVolumeLogicalPointCount(cell.compact_rule) == 0) continue;
		const auto element = domain.Background().MaterializeElement(id);
		iga::ForEachVolumePoint(volume.UsableCompactRule(domain, id), [&](const iga::VolumeQuadraturePoint& point) {
			callback(element, point);
		});
	}
}

// This deliberately does not reuse the catalog's 4-point integration over
// certified blocks.  It streams a tensor 8x8x8 Gauss rule per compact block.
// Unresolved sample leaves are likewise re-integrated and reclassified point
// by point against the exact box surface; no global logical-point expansion is
// materialized.  Thus the rule is independent in both quadrature order and
// (for unresolved leaves) occupancy treatment.
template <class Callback> void VisitIndependentEight(const iga::CartesianDomainClassification& domain,
	const iga::CutCellVolumeQuadratureCatalog& volume, Callback&& callback)
{
	constexpr std::array<double, 8> points{{0.019855071751231884, 0.10166676129318664, 0.2372337950418355, 0.4082826787521751,
		0.5917173212478249, 0.7627662049581645, 0.8983332387068134, 0.9801449282487681}};
	constexpr std::array<double, 8> weights{{0.05061426814518813, 0.11119051722668724, 0.15685332293894365, 0.18134189168918099,
		0.18134189168918099, 0.15685332293894365, 0.11119051722668724, 0.05061426814518813}};
	assert(volume.StorageMode() == iga::CutCellVolumeQuadratureStorageMode::Compact);
	for (std::uint64_t id = 0; id < domain.Cells().size(); ++id) {
		const auto classification = domain.Cells()[static_cast<std::size_t>(id)].classification;
		const auto& cell = volume.Cell(id);
		if ((classification != iga::CellClassification::Inside && classification != iga::CellClassification::Cut)
			|| !cell.usable) continue;
		// Do not reject a compact cut leaf because its original 4^3 mask was
		// empty: the independent 8^3 rule classifies every unresolved point
		// against the exact surface predicate below.
		const auto element = domain.Background().MaterializeElement(id); const auto& rule = cell.compact_rule;
		const double lattice = static_cast<double>(std::uint64_t(1) << rule.max_depth);
		const auto integrate = [&](const std::array<std::uint32_t,3>& lower, const std::array<std::uint32_t,3>& upper, bool reclassify) {
			const double reference_volume = (static_cast<double>(upper[0]-lower[0])/lattice)*(static_cast<double>(upper[1]-lower[1])/lattice)*(static_cast<double>(upper[2]-lower[2])/lattice);
			for (std::size_t qz = 0; qz < 8; ++qz) for (std::size_t qy = 0; qy < 8; ++qy) for (std::size_t qx = 0; qx < 8; ++qx) {
				std::array<double,3> parametric{{static_cast<double>(lower[0])/lattice+static_cast<double>(upper[0]-lower[0])/lattice*points[qx],
					static_cast<double>(lower[1])/lattice+static_cast<double>(upper[1]-lower[1])/lattice*points[qy],
					static_cast<double>(lower[2])/lattice+static_cast<double>(upper[2]-lower[2])/lattice*points[qz]}};
				if (reclassify && domain.SurfaceIndex().LocatePoint(iga::EvaluateElementGeometry(element, parametric).physical) != iga::PointLocation::Inside) continue;
				callback(element, iga::VolumeQuadraturePoint{parametric, reference_volume*weights[qx]*weights[qy]*weights[qz]});
			}
		};
		for (const auto& block : rule.certified_blocks) integrate(block.lower, block.upper, false);
		for (const auto& leaf : rule.sample_leaves) {
			const std::uint32_t width = std::uint32_t(1) << (rule.max_depth-leaf.depth); std::array<std::uint32_t,3> upper{};
			for (int axis = 0; axis < 3; ++axis) upper[axis] = leaf.key[axis]*width+width;
			integrate({{leaf.key[0]*width,leaf.key[1]*width,leaf.key[2]*width}}, upper, true);
		}
	}
}

ErrorNorms MeasureErrors(const iga::CartesianDomainClassification& domain,
	const iga::CutCellVolumeQuadratureCatalog& independent_volume, const iga::ImmersedStaticFlowRuntime& runtime, bool eight_point = false)
{
	const auto state = runtime.CommittedState();
	ErrorNorms result;
	const auto visit = [&](const auto& consume) { if (eight_point) VisitIndependentEight(domain, independent_volume, consume); else VisitPositive(domain, independent_volume, consume); };
	visit([&](const iga::Element& element, const iga::VolumeQuadraturePoint& point) {
		const auto basis = iga::EvaluateBasis(element, point.parametric[0], point.parametric[1], point.parametric[2], false);
		const auto physical = iga::EvaluateElementGeometry(element, point.parametric).physical;
		const double weight = point.weight*basis.raw_determinant;
		double ph = 0.0, pe = ExactPressure(physical);
		for (std::size_t a = 0; a < element.connectivity.size(); ++a)
			ph += basis.value[a]*PetscRealPart(state[static_cast<std::size_t>(runtime.Dof(element.connectivity[a], 3))]);
		result.measure += weight; result.pressure_mean += weight*ph; result.pressure_error_shift += weight*(ph-pe);
	});
	assert(result.measure > 0.0);
	result.pressure_mean /= result.measure; result.pressure_error_shift /= result.measure;
	visit([&](const iga::Element& element, const iga::VolumeQuadraturePoint& point) {
		const auto basis = iga::EvaluateBasis(element, point.parametric[0], point.parametric[1], point.parametric[2], false);
		const auto physical = iga::EvaluateElementGeometry(element, point.parametric).physical;
		const double weight = point.weight*basis.raw_determinant;
		std::array<double, 3> uh{}; std::array<std::array<double, 3>, 3> gh{}; double ph = 0.0, divergence = 0.0;
		for (std::size_t a = 0; a < element.connectivity.size(); ++a) {
			for (int c = 0; c < 3; ++c) {
				const double value = PetscRealPart(state[static_cast<std::size_t>(runtime.Dof(element.connectivity[a], c))]);
				uh[c] += basis.value[a]*value;
				for (int d = 0; d < 3; ++d) gh[c][d] += basis.gradient[a][d]*value;
			}
			ph += basis.value[a]*PetscRealPart(state[static_cast<std::size_t>(runtime.Dof(element.connectivity[a], 3))]);
		}
		for (int c = 0; c < 3; ++c) divergence += gh[c][c];
		const auto ue = ExactVelocity(physical); const auto ge = ExactVelocityGradient(physical);
		for (int c = 0; c < 3; ++c) {
			result.velocity_l2 += weight*(uh[c]-ue[c])*(uh[c]-ue[c]);
			result.numerical_velocity_l2 += weight*uh[c]*uh[c]; result.exact_velocity_l2 += weight*ue[c]*ue[c]; result.velocity_dot_exact += weight*uh[c]*ue[c];
			for (int d = 0; d < 3; ++d) result.velocity_h1 += weight*(gh[c][d]-ge[c][d])*(gh[c][d]-ge[c][d]);
			for (int d = 0; d < 3; ++d) result.exact_velocity_h1 += weight*ge[c][d]*ge[c][d];
		}
		const double pd = ph-ExactPressure(physical)-result.pressure_error_shift;
		result.pressure_l2 += weight*pd*pd; result.divergence_l2 += weight*divergence*divergence; result.exact_pressure_l2 += weight*ExactPressure(physical)*ExactPressure(physical);
	});
	result.velocity_l2 = std::sqrt(result.velocity_l2); result.velocity_h1 = std::sqrt(result.velocity_h1);
	result.pressure_l2 = std::sqrt(result.pressure_l2); result.divergence_l2 = std::sqrt(result.divergence_l2); result.exact_velocity_l2 = std::sqrt(result.exact_velocity_l2); result.exact_velocity_h1 = std::sqrt(result.exact_velocity_h1); result.exact_pressure_l2 = std::sqrt(result.exact_pressure_l2); result.numerical_velocity_l2 = std::sqrt(result.numerical_velocity_l2);
	for (double value : {result.velocity_l2,result.velocity_h1,result.pressure_l2,result.divergence_l2,result.pressure_mean,result.pressure_error_shift,result.exact_velocity_l2,result.exact_velocity_h1,result.exact_pressure_l2,result.numerical_velocity_l2,result.velocity_dot_exact}) assert(std::isfinite(value));
	return result;
}

double MaxDifference(const std::vector<PetscScalar>& left, const std::vector<PetscScalar>& right)
{
	assert(left.size() == right.size()); double result = 0.0;
	for (std::size_t i = 0; i < left.size(); ++i) result = std::max(result, std::abs(PetscRealPart(left[i]-right[i])));
	return result;
}

double VectorNorm(const std::vector<PetscScalar>& values)
{
	double squared = 0.0; for (const auto value : values) squared += PetscRealPart(value)*PetscRealPart(value);
	return std::sqrt(squared);
}

struct ResidualBlocks { double momentum = 0.0, continuity = 0.0, gauge = 0.0; };
ResidualBlocks SplitResidualBlocks(const iga::ImmersedStaticFlowRuntime& runtime, const std::vector<PetscScalar>& residual)
{
	assert(residual.size() == runtime.Diagnostics().total_dofs && residual.size() == runtime.Diagnostics().physical_dofs+1);
	ResidualBlocks result;
	for (std::size_t row = 0; row < residual.size(); ++row) {
		const double value = PetscRealPart(residual[row]);
		if (row == static_cast<std::size_t>(runtime.GaugeDof())) result.gauge += value*value;
		else if ((row%4) == 3) result.continuity += value*value;
		else result.momentum += value*value;
	}
	result.momentum = std::sqrt(result.momentum); result.continuity = std::sqrt(result.continuity); result.gauge = std::sqrt(result.gauge); return result;
}

void PrintResidualBlocks(unsigned int cells, const char* stage, const ResidualBlocks& blocks)
{
	std::cout << "manufactured_residual_blocks cells=" << cells << " stage=" << stage
		<< " ordering=active_node_x4_u0_u1_u2_p_then_gauge negative_residual"
		<< " momentum=" << blocks.momentum << " continuity=" << blocks.continuity << " gauge=" << blocks.gauge << '\n';
}

struct FiniteDifferenceDefect {
	double mismatch_norm = 0.0, action_norm = 0.0, relative_defect = 0.0;
};

FiniteDifferenceDefect ManufacturedJacobianConsistency(iga::ImmersedStaticFlowRuntime& runtime, double epsilon)
{
	std::vector<PetscScalar> zero(runtime.Diagnostics().total_dofs, 0.0);
	runtime.SetCommittedState(zero); runtime.Assemble(); const auto base = runtime.AssembledNegativeResidual();
	std::vector<PetscScalar> direction(base.size(), 0.0), plus(base.size(), 0.0), minus(base.size(), 0.0);
	for (std::size_t i = 0; i < direction.size(); ++i) direction[i] = std::sin(static_cast<double>(17*i+3));
	const double direction_norm = VectorNorm(direction); assert(direction_norm > 0.0);
	for (auto& value : direction) value /= direction_norm;
	// J is dR/dx, while the stored vector is b=-R.  Save J(0)d before
	// perturbing the state; no perturbed Jacobian enters this diagnostic.
	const auto action = runtime.AssembledJacobianAction(direction);
	for (std::size_t i = 0; i < plus.size(); ++i) { plus[i] = epsilon*direction[i]; minus[i] = -epsilon*direction[i]; }
	runtime.SetCommittedState(plus); runtime.Assemble(); plus = runtime.AssembledNegativeResidual();
	runtime.SetCommittedState(minus); runtime.Assemble(); minus = runtime.AssembledNegativeResidual();
	double error = 0.0, scale = 0.0;
	for (std::size_t i = 0; i < base.size(); ++i) {
		const double centered_residual_action = -PetscRealPart((plus[i]-minus[i])/(2.0*epsilon));
		const double expected = PetscRealPart(action[i]);
		error += (centered_residual_action-expected)*(centered_residual_action-expected); scale += expected*expected;
	}
	runtime.SetCommittedState(zero); runtime.Assemble();
	const double mismatch = std::sqrt(error), denominator = std::sqrt(scale);
	return {mismatch, denominator, mismatch/std::max(denominator, std::numeric_limits<double>::min())};
}

struct AssembledContribution { std::vector<PetscScalar> matrix, residual; };

AssembledContribution AssembleContribution(const iga::CartesianDomainClassification& domain,
	const iga::CutCellVolumeQuadratureCatalog& volume, const iga::ImmersedSurfaceQuadratureCatalog& surface,
	const iga::CutCellGhostPenaltyCatalog& ghost, iga::ImmersedStaticFlowOptions options,
	const std::vector<PetscScalar>& state)
{
	iga::ImmersedStaticFlowRuntime runtime(domain, volume, surface, ghost, std::move(options));
	runtime.SetCommittedState(state); runtime.Assemble();
	return {runtime.AssembledJacobianDense(), runtime.AssembledNegativeResidual()};
}

std::vector<PetscScalar> SumContributions(const std::vector<AssembledContribution>& parts, bool matrix)
{
	assert(!parts.empty()); const std::size_t size = matrix ? parts.front().matrix.size() : parts.front().residual.size();
	std::vector<PetscScalar> result(size, 0.0);
	for (const auto& part : parts) {
		const auto& values = matrix ? part.matrix : part.residual; assert(values.size() == size);
		for (std::size_t i = 0; i < size; ++i) result[i] += values[i];
	}
	return result;
}

void AssertContributionDecomposition(const iga::CartesianDomainClassification& domain,
	const iga::CutCellVolumeQuadratureCatalog& volume, const iga::ImmersedSurfaceQuadratureCatalog& surface,
	const iga::CutCellGhostPenaltyCatalog& ghost, const iga::ImmersedStaticFlowOptions& base)
{
	std::vector<PetscScalar> state(4*base.wall_labels.size(), 0.0); // resized below from a reference assembly.
	iga::ImmersedStaticFlowRuntime reference(domain, volume, surface, ghost, base);
	state.assign(reference.Diagnostics().total_dofs, 0.0);
	for (std::size_t i = 0; i < state.size(); ++i) state[i] = .031*std::sin(static_cast<double>(11*i+5));
	auto options = base; options.assemble_gauge = false;
	auto volume_only = options; volume_only.assemble_volume = true; volume_only.assemble_wall = volume_only.assemble_ghost = false;
	auto wall_only = options; wall_only.assemble_volume = false; wall_only.assemble_wall = true; wall_only.assemble_ghost = false;
	auto ghost_only = options; ghost_only.assemble_volume = ghost_only.assemble_wall = false; ghost_only.assemble_ghost = true;
	auto gauge_only = options; gauge_only.assemble_volume = gauge_only.assemble_wall = gauge_only.assemble_ghost = false; gauge_only.assemble_gauge = true;
	auto full_ungauged = options;
	auto full_gauged = base;
	const auto v = AssembleContribution(domain, volume, surface, ghost, volume_only, state);
	const auto w = AssembleContribution(domain, volume, surface, ghost, wall_only, state);
	const auto g = AssembleContribution(domain, volume, surface, ghost, ghost_only, state);
	const auto q = AssembleContribution(domain, volume, surface, ghost, gauge_only, state);
	const auto u = AssembleContribution(domain, volume, surface, ghost, full_ungauged, state);
	const auto f = AssembleContribution(domain, volume, surface, ghost, full_gauged, state);
	const double tolerance = 2e-11;
	assert(MaxDifference(u.matrix, SumContributions({v,w,g}, true)) < tolerance);
	assert(MaxDifference(u.residual, SumContributions({v,w,g}, false)) < tolerance);
	assert(MaxDifference(f.matrix, SumContributions({u,q}, true)) < tolerance);
	assert(MaxDifference(f.residual, SumContributions({u,q}, false)) < tolerance);
	std::cout << "manufactured_contribution_decomposition state=nonzero matrix_residual_tolerance=" << tolerance << '\n';
}

void PrintNewtonDiagnostics(unsigned int cells, const iga::ImmersedStaticFlowRuntime& runtime,
	const iga::ImmersedStaticFlowOptions& options, double wall_seconds)
{
	const auto& diagnostics = runtime.Diagnostics();
	const double initial = diagnostics.newton_steps.empty() ? diagnostics.residual_norm : diagnostics.newton_steps.front().residual_norm;
	const double threshold = std::max(options.nonlinear_absolute_tolerance, options.nonlinear_relative_tolerance*initial);
	std::cout << "manufactured_solver cells=" << cells << " ksp_rtol=" << options.ksp_relative_tolerance
		<< " lu_shift=" << options.lu_pivot_shift << " initial=" << initial << " target=" << threshold
		<< " final=" << diagnostics.residual_norm << " newton=" << diagnostics.nonlinear_iterations
		<< " ksp_total=" << diagnostics.ksp_iterations << " reason=" << static_cast<int>(diagnostics.ksp_reason)
		<< " wall_s=" << wall_seconds << '\n';
	for (const auto& step : diagnostics.newton_steps)
		std::cout << "manufactured_newton cells=" << cells << " i=" << step.iteration << " r=" << step.residual_norm
			<< " update=" << step.update_norm << " ksp_r=" << step.ksp_residual_norm << " linear_r=" << step.linear_residual_norm << " linear_rel=" << step.linear_relative_residual << " ksp=" << step.ksp_iterations
			<< " reason=" << static_cast<int>(step.ksp_reason) << " damping=" << step.damping << " candidate=" << step.candidate_residual_norm << '\n';
}

struct Spectrum { double min_abs = 0.0, max_abs = 0.0, condition = 0.0, symmetry_defect = 0.0; };
Spectrum MixedGaugedSpectrum(const iga::ImmersedStaticFlowRuntime& runtime)
{
	const std::size_t n = runtime.Diagnostics().total_dofs;
	auto matrix = runtime.AssembledJacobianDense();
	if (matrix.size() != n*n || n > static_cast<std::size_t>(std::numeric_limits<PetscBLASInt>::max())) throw std::runtime_error("dense mixed spectrum dimension is invalid");
	// At zero state the unscaled mixed operator becomes symmetric after the
	// orthogonal left sign flip on pressure and gauge rows.  This preserves all
	// singular values; its real eigenvalues are the requested raw spectrum.
	for (std::size_t row = 0; row < n; ++row) if ((row%4) == 3 || row+1 == n)
		for (std::size_t column = 0; column < n; ++column) matrix[row*n+column] = -matrix[row*n+column];
	double scale = 0.0, defect = 0.0;
	for (std::size_t row = 0; row < n; ++row) for (std::size_t column = 0; column < n; ++column) {
		scale = std::max(scale, std::abs(PetscRealPart(matrix[row*n+column])));
		defect = std::max(defect, std::abs(PetscRealPart(matrix[row*n+column]-matrix[column*n+row])));
	}
	const double relative_defect = defect/std::max(1.0, scale);
	if (!std::isfinite(relative_defect) || relative_defect > 2e-10) throw std::runtime_error("sign-flipped mixed operator is not symmetric");
	PetscBLASInt size = static_cast<PetscBLASInt>(n), leading = size, info = 0, lwork = std::max<PetscBLASInt>(1, 3*size-1);
	std::vector<PetscReal> eigenvalues(n); std::vector<PetscScalar> workspace(static_cast<std::size_t>(lwork));
	char job = 'N', triangle = 'U';
	LAPACKsyev_(&job, &triangle, &size, matrix.data(), &leading, eigenvalues.data(), workspace.data(), &lwork, &info);
	if (info != 0) throw std::runtime_error("LAPACK SYEV failed to converge for mixed operator");
	Spectrum result; result.symmetry_defect = relative_defect; result.min_abs = std::numeric_limits<double>::infinity();
	for (const auto eigenvalue : eigenvalues) {
		const double magnitude = std::abs(static_cast<double>(eigenvalue));
		if (!std::isfinite(magnitude)) throw std::runtime_error("mixed spectrum has a non-finite eigenvalue");
		result.min_abs = std::min(result.min_abs, magnitude); result.max_abs = std::max(result.max_abs, magnitude);
	}
	result.condition = result.max_abs/result.min_abs;
	if (!(result.min_abs > 0.0) || !(result.max_abs >= result.min_abs) || !std::isfinite(result.condition)) throw std::runtime_error("mixed spectrum is singular or invalid");
	return result;
}

struct SliverResult { int m = 0; double fraction = 0.0, defect = 0.0, gauge = 0.0, box_measure = 0.0, normalized_defect = 0.0; std::size_t dofs = 0, faces = 0, records = 0, logical_points = 0; Spectrum spectrum; };
SliverResult RunSliver(int m)
{
	const iga::CubicCartesianGridSpec spec{{{0.0,0.0,0.0}},{{3.0,3.0,3.0}},{{3,3,3}}};
	const double tip = 2.0+std::ldexp(1.0, -m);
	const auto domain = Domain(spec, {{0.125,0.125,0.125}}, {{tip,2.875,2.875}});
	const iga::CutCellVolumeQuadratureCatalog volume(domain, CompactOptions(8), iga::CutCellVolumeQuadratureStorageMode::Compact);
	const iga::ImmersedSurfaceQuadratureCatalog surface(domain); const iga::CutCellGhostPenaltyCatalog ghost(domain, volume);
	iga::ImmersedStaticFlowOptions options; options.parameters = {1.0,1.0,0.0}; options.wall_labels = {7}; options.wall_gamma0 = 2.0;
	iga::ImmersedStaticFlowRuntime runtime(domain, volume, surface, ghost, options); runtime.Assemble();
	std::size_t records = 0, logical = 0, expected_volume_cells = 0, expected_surface_cells = 0; double minimum = 1.0;
	for (std::uint64_t id = 0; id < domain.Cells().size(); ++id) {
		const auto& cell = volume.Cell(id); assert(cell.rule.Points().empty()); records += cell.diagnostics.certified_blocks+cell.diagnostics.sample_leaves; logical += cell.diagnostics.logical_output_points;
		const auto classification = domain.Cells()[static_cast<std::size_t>(id)].classification;
		const bool positive = cell.usable && (classification == iga::CellClassification::Inside || classification == iga::CellClassification::Cut) && iga::CompactCutCellVolumeLogicalPointCount(cell.compact_rule) != 0;
		if (positive) ++expected_volume_cells;
		if (classification == iga::CellClassification::Cut && positive) {
			++expected_surface_cells;
			assert(ghost.Covered(id)); minimum = std::min(minimum, cell.diagnostics.estimated_reference_volume);
		}
	}
	const double expected_fraction = (49.0/64.0)*std::ldexp(1.0, -m);
	assert(std::abs(minimum-expected_fraction) <= 5e-5*expected_fraction);
	const double box_measure = (tip-.125)*2.75*2.75;
	const auto& initial = runtime.Diagnostics();
	assert(initial.volume_cells == expected_volume_cells && initial.surface_cells == expected_surface_cells && initial.ghost_faces == ghost.Faces().size() && !ghost.Faces().empty());
	assert(std::abs(initial.pressure_measure-box_measure) <= 5e-8*box_measure);
	const auto residual = runtime.AssembledNegativeResidual(); std::vector<PetscScalar> probe(runtime.Diagnostics().total_dofs, 0.0);
	for (std::size_t i = 0; i < probe.size(); ++i) probe[i] = static_cast<double>((13*i)%37)/37.0;
	const auto action = runtime.AssembledJacobianAction(probe);
	runtime.Assemble(); // deterministic repeated global assembly without a second runtime.
	assert(MaxDifference(residual, runtime.AssembledNegativeResidual()) == 0.0 && MaxDifference(action, runtime.AssembledJacobianAction(probe)) == 0.0);
	const auto& d = runtime.Diagnostics(); const auto spectrum = MixedGaugedSpectrum(runtime);
	const double normalized_defect = d.constant_pressure_defect/std::max(1.0, spectrum.max_abs);
	assert(d.total_dofs == 4*d.active_nodes+1 && normalized_defect < 2e-10);
	return {m,minimum,d.constant_pressure_defect,d.pressure_measure,box_measure,normalized_defect,d.total_dofs,ghost.Faces().size(),records,logical,spectrum};
}

} // namespace

int main(int argc, char** argv)
{
	if (PetscInitialize(&argc, &argv, nullptr, nullptr) != 0) return 1;
	try {
		std::cout << std::setprecision(10) << std::scientific;
		int manufactured_level = 0, sliver_m = 0; bool aggregate = argc == 1, manufactured_triplet = false;
		double manufactured_ksp_rtol = 1e-10, manufactured_lu_shift = 1e-12;
		double manufactured_outer_rtol = 1e-9, manufactured_outer_atol = 1e-14;
		for (int argument = 1; argument < argc; ++argument) {
			const std::string option(argv[argument]);
			if (option == "--all") aggregate = true;
			else if (option == "--manufactured-triplet") manufactured_triplet = true;
			else if (option == "--manufactured-level" && argument+1 < argc) manufactured_level = std::stoi(argv[++argument]);
			else if (option == "--manufactured-ksp-rtol" && argument+1 < argc) manufactured_ksp_rtol = std::stod(argv[++argument]);
			else if (option == "--manufactured-outer-rtol" && argument+1 < argc) manufactured_outer_rtol = std::stod(argv[++argument]);
			else if (option == "--manufactured-outer-atol" && argument+1 < argc) manufactured_outer_atol = std::stod(argv[++argument]);
			else if (option == "--manufactured-lu-shift" && argument+1 < argc) manufactured_lu_shift = std::stod(argv[++argument]);
			else if (option == "--sliver-m" && argument+1 < argc) sliver_m = std::stoi(argv[++argument]);
			else throw std::invalid_argument("usage: --all | --manufactured-triplet | --manufactured-level {3,6,12} [--manufactured-ksp-rtol value] [--manufactured-outer-rtol value] [--manufactured-outer-atol value] [--manufactured-lu-shift value] | --sliver-m {1..8}");
		}
		if ((manufactured_level != 0 && manufactured_level != 3 && manufactured_level != 6 && manufactured_level != 12) || (sliver_m < 0 || sliver_m > 8)) throw std::invalid_argument("requested closure case is outside the supported nested family");
		if ((manufactured_level != 0 && sliver_m != 0) || (aggregate && (manufactured_level != 0 || sliver_m != 0 || manufactured_triplet)) || (manufactured_triplet && (manufactured_level != 0 || sliver_m != 0))) throw std::invalid_argument("select exactly one closure selector");
		if (!std::isfinite(manufactured_ksp_rtol) || !(manufactured_ksp_rtol > 0.0) || !std::isfinite(manufactured_outer_rtol) || !(manufactured_outer_rtol > 0.0) || !std::isfinite(manufactured_outer_atol) || !(manufactured_outer_atol > 0.0) || !std::isfinite(manufactured_lu_shift) || manufactured_lu_shift < 0.0) throw std::invalid_argument("manufactured solver overrides are invalid");
		const auto run_manufactured = [&](unsigned int cells) {
			const iga::CubicCartesianGridSpec spec{{{0.0,0.0,0.0}},{{3.0,3.0,3.0}},{{cells,cells,cells}}};
			// These dyadic, non-grid-aligned bounds retain cut cells and a fully
			// interior anchor at every level; finite octree depth can resolve them.
			const auto domain = Domain(spec, {{kBoxLower,kBoxLower,kBoxLower}}, {{kBoxUpper,kBoxUpper,kBoxUpper}});
			const iga::CutCellVolumeQuadratureCatalog volume(domain, CompactOptions(6), iga::CutCellVolumeQuadratureStorageMode::Compact);
			const iga::ImmersedSurfaceQuadratureCatalog surface(domain); const iga::CutCellGhostPenaltyCatalog ghost(domain, volume);
			auto options = ManufacturedOptions(); options.ksp_relative_tolerance = manufactured_ksp_rtol; options.nonlinear_relative_tolerance = manufactured_outer_rtol; options.nonlinear_absolute_tolerance = manufactured_outer_atol; options.lu_pivot_shift = manufactured_lu_shift;
			iga::ImmersedStaticFlowRuntime runtime(domain, volume, surface, ghost, options);
			std::size_t expected_volume_cells = 0, expected_surface_cells = 0; std::vector<std::int32_t> expected_active;
			for (std::uint64_t id = 0; id < domain.Cells().size(); ++id) {
				const auto classification = domain.Cells()[static_cast<std::size_t>(id)].classification; const auto& cell = volume.Cell(id);
				const bool positive = cell.usable && (classification == iga::CellClassification::Inside || classification == iga::CellClassification::Cut) && iga::CompactCutCellVolumeLogicalPointCount(cell.compact_rule) != 0;
				if (!positive) continue;
				++expected_volume_cells; if (classification == iga::CellClassification::Cut) ++expected_surface_cells;
				const auto element = domain.Background().MaterializeElement(id); expected_active.insert(expected_active.end(), element.connectivity.begin(), element.connectivity.end());
			}
			std::sort(expected_active.begin(), expected_active.end()); expected_active.erase(std::unique(expected_active.begin(), expected_active.end()), expected_active.end());
			if (cells == 3) {
				AssertContributionDecomposition(domain, volume, surface, ghost, options);
				const std::array<double,4> epsilons{{1e-4,1e-5,1e-6,1e-7}};
				for (const auto& contribution : std::array<std::pair<const char*, std::array<bool,3>>,3>{{
					{"full", {{true,true,true}}}, {"volume_only", {{true,false,false}}}, {"wall_ghost_only", {{false,true,true}}}}}) {
					for (const double epsilon : epsilons) {
						auto diagnostic_options = options; diagnostic_options.assemble_volume = contribution.second[0]; diagnostic_options.assemble_wall = contribution.second[1]; diagnostic_options.assemble_ghost = contribution.second[2];
						iga::ImmersedStaticFlowRuntime diagnostic_runtime(domain, volume, surface, ghost, diagnostic_options);
						const auto defect = ManufacturedJacobianConsistency(diagnostic_runtime, epsilon);
						std::cout << "manufactured_fd cells=" << cells << " contribution=" << contribution.first << " epsilon=" << epsilon
							<< " mismatch_norm=" << defect.mismatch_norm << " jacobian_action_norm=" << defect.action_norm
							<< " relative_defect=" << defect.relative_defect << " convention=minus_centered_negative_residual_equals_Jd" << '\n';
					}
				}
			}
			const auto jacobian_consistency = ManufacturedJacobianConsistency(runtime, 1e-6);
			std::cout << "manufactured_jacobian cells=" << cells << " relative_fd_error=" << jacobian_consistency.relative_defect << '\n';
			runtime.Assemble(); const auto initial_residual = runtime.AssembledNegativeResidual(); const auto initial_blocks = SplitResidualBlocks(runtime, initial_residual); const double initial_unshifted = VectorNorm(initial_residual);
			PrintResidualBlocks(cells, "initial", initial_blocks);
			const auto solve_start = std::chrono::steady_clock::now();
			try { assert(runtime.SolveTrial()); }
			catch (...) { PrintNewtonDiagnostics(cells, runtime, options, std::chrono::duration<double>(std::chrono::steady_clock::now()-solve_start).count()); throw; }
			PrintNewtonDiagnostics(cells, runtime, options, std::chrono::duration<double>(std::chrono::steady_clock::now()-solve_start).count()); runtime.Commit();
			const iga::CutCellVolumeQuadratureCatalog error_volume(domain, CompactOptions(7), iga::CutCellVolumeQuadratureStorageMode::Compact);
			const auto baseline_norms = MeasureErrors(domain, error_volume, runtime);
			const auto norms = MeasureErrors(domain, error_volume, runtime, true);
			const double true_unshifted_residual = VectorNorm(runtime.AssembledNegativeResidual());
			const auto final_blocks = SplitResidualBlocks(runtime, runtime.AssembledNegativeResidual()); PrintResidualBlocks(cells, "final", final_blocks);
			const double relative_target = options.nonlinear_relative_tolerance*initial_unshifted;
			const double gauge_scale_aware_tolerance = options.nonlinear_absolute_tolerance+relative_target;
			std::cout << "manufactured_convergence_control cells=" << cells << " initial=" << initial_unshifted << " relative_target=" << relative_target << " absolute_target=" << options.nonlinear_absolute_tolerance << " final=" << true_unshifted_residual << " gauge_scale_aware_target=" << gauge_scale_aware_tolerance << '\n';
			assert(relative_target > options.nonlinear_absolute_tolerance && std::isfinite(true_unshifted_residual) && true_unshifted_residual <= relative_target);
			for (const auto& pair : std::array<std::pair<double,double>,3>{{{initial_blocks.momentum,final_blocks.momentum},{initial_blocks.continuity,final_blocks.continuity},{initial_blocks.gauge,final_blocks.gauge}}}) if (pair.first > 0.0) assert(pair.second <= pair.first/1e9); else assert(pair.second <= gauge_scale_aware_tolerance);
			const auto integration_change = [](double baseline, double independent) { return std::abs(independent-baseline)/std::max({std::abs(baseline),std::abs(independent),std::numeric_limits<double>::min()}); };
			std::cout << "manufactured_error_integration cells=" << cells << " baseline_uL2=" << baseline_norms.velocity_l2 << " high_uL2=" << norms.velocity_l2 << " baseline_uH1=" << baseline_norms.velocity_h1 << " high_uH1=" << norms.velocity_h1 << " baseline_pL2=" << baseline_norms.pressure_l2 << " high_pL2=" << norms.pressure_l2 << " baseline_div=" << baseline_norms.divergence_l2 << " high_div=" << norms.divergence_l2 << '\n';
			assert(integration_change(baseline_norms.velocity_l2,norms.velocity_l2) < 2e-3 && integration_change(baseline_norms.velocity_h1,norms.velocity_h1) < 2e-3 && integration_change(baseline_norms.pressure_l2,norms.pressure_l2) < 2e-3 && integration_change(baseline_norms.divergence_l2,norms.divergence_l2) < 2e-3);
			const auto& diagnostics = runtime.Diagnostics();
			assert(diagnostics.converged && diagnostics.volume_cells == expected_volume_cells && diagnostics.surface_cells == expected_surface_cells && diagnostics.ghost_faces == ghost.Faces().size());
			assert(diagnostics.active_nodes == expected_active.size() && diagnostics.total_dofs == 4*diagnostics.active_nodes+1);
			assert(std::abs(diagnostics.pressure_measure-kLength*kLength*kLength) <= 5e-8*kLength*kLength*kLength && std::abs(norms.pressure_mean) < 1e-10);
			return std::make_pair(norms, runtime.Diagnostics().total_dofs);
		};
		if (aggregate || manufactured_triplet || manufactured_level != 0) {
			const std::array<unsigned int,3> levels{{3,6,12}}; std::vector<ErrorNorms> errors; std::vector<std::size_t> dofs;
			std::cout << "manufactured h dofs uL2 uH1 pL2 divL2 pmean pshift\n";
			for (const auto cells : levels) if (aggregate || manufactured_triplet || static_cast<int>(cells) == manufactured_level) {
				const auto result = run_manufactured(cells); errors.push_back(result.first); dofs.push_back(result.second); const auto& norms = result.first;
				std::cout << (3.0/cells) << ' ' << result.second << ' ' << norms.velocity_l2 << ' ' << norms.velocity_h1 << ' ' << norms.pressure_l2 << ' ' << norms.divergence_l2 << ' ' << norms.pressure_mean << ' ' << norms.pressure_error_shift << '\n';
				std::cout << "manufactured_error_norms cells=" << cells << " u_exact_l2=" << norms.exact_velocity_l2 << " u_exact_h1=" << norms.exact_velocity_h1 << " p_exact_l2=" << norms.exact_pressure_l2 << " uh_l2=" << norms.numerical_velocity_l2 << " uh_dot_ue=" << norms.velocity_dot_exact
					<< " normalized_uL2=" << norms.velocity_l2/norms.exact_velocity_l2 << " normalized_uH1=" << norms.velocity_h1/norms.exact_velocity_h1 << " normalized_pL2=" << norms.pressure_l2/norms.exact_pressure_l2 << " normalized_div=" << norms.divergence_l2/norms.exact_velocity_h1 << '\n';
			}
			if (aggregate || manufactured_triplet) {
				std::cout << "manufactured pair_orders uL2 uH1 pL2 divL2\n";
				for (std::size_t i = 1; i < errors.size(); ++i) {
					const auto order = [&](double coarse, double fine) { return std::log(coarse/fine)/std::log(2.0); };
					const std::array<double,4> values{{order(errors[i-1].velocity_l2,errors[i].velocity_l2),order(errors[i-1].velocity_h1,errors[i].velocity_h1),order(errors[i-1].pressure_l2,errors[i].pressure_l2),order(errors[i-1].divergence_l2,errors[i].divergence_l2)}};
					std::cout << values[0] << ' ' << values[1] << ' ' << values[2] << ' ' << values[3] << '\n';
					assert(errors[i].velocity_l2 < errors[i-1].velocity_l2 && errors[i].velocity_h1 < errors[i-1].velocity_h1 && errors[i].pressure_l2 < errors[i-1].pressure_l2 && errors[i].divergence_l2 < errors[i-1].divergence_l2);
				}
				assert(dofs.size() == 3 && dofs[0] < dofs[1] && dofs[1] < dofs[2]);
				const auto order = [&](double coarse, double fine) { return std::log(coarse/fine)/std::log(2.0); };
				assert(order(errors[1].velocity_l2,errors[2].velocity_l2) >= 2.5 && order(errors[1].velocity_h1,errors[2].velocity_h1) >= 2.0 && order(errors[1].pressure_l2,errors[2].pressure_l2) >= 2.0 && order(errors[1].divergence_l2,errors[2].divergence_l2) >= 1.5);
				assert(errors.back().velocity_l2/errors.back().exact_velocity_l2 < .02 && errors.back().velocity_h1/errors.back().exact_velocity_h1 < .05 && errors.back().pressure_l2/errors.back().exact_pressure_l2 < .05 && errors.back().divergence_l2/errors.back().exact_velocity_h1 < .05);
			} else if (manufactured_level == 12) {
				const auto& finest = errors.back();
				assert(finest.velocity_l2/finest.exact_velocity_l2 < .02 && finest.velocity_h1/finest.exact_velocity_h1 < .05 && finest.pressure_l2/finest.exact_pressure_l2 < .05 && finest.divergence_l2/finest.exact_velocity_h1 < .05);
				std::cout << "manufactured_level12_local_caps_only not_triplet_convergence_evidence\n";
			}
		}
		if (aggregate || sliver_m != 0) {
			std::cout << "sliver m alpha_min defect gauge box_measure normalized_defect dofs faces records logical min_abs max_abs condition symmetry\n";
			std::vector<SliverResult> slivers;
			for (int m = 1; m <= 8; ++m) if (aggregate || m == sliver_m) {
				slivers.push_back(RunSliver(m)); const auto& r = slivers.back();
				std::cout << r.m << ' ' << r.fraction << ' ' << r.defect << ' ' << r.gauge << ' ' << r.box_measure << ' ' << r.normalized_defect << ' ' << r.dofs << ' ' << r.faces << ' ' << r.records << ' ' << r.logical_points << ' ' << r.spectrum.min_abs << ' ' << r.spectrum.max_abs << ' ' << r.spectrum.condition << ' ' << r.spectrum.symmetry_defect << '\n';
			}
			if (aggregate) { const double growth = slivers.back().spectrum.condition/slivers.front().spectrum.condition; assert(std::isfinite(growth) && growth < 100.0); std::cout << "sliver_condition_growth " << growth << '\n'; }
		}
		return PetscFinalize() == 0 ? 0 : 1;
	} catch (const std::exception& error) {
		std::cerr << "phase5 closure failure: " << error.what() << '\n';
		PetscFinalize(); return 1;
	} catch (...) {
		std::cerr << "phase5 closure failure: non-standard exception\n";
		PetscFinalize(); return 1;
	}
}
