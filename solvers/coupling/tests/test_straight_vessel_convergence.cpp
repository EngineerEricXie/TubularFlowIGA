#include "CouplingFixture.hpp"
#include "BoundarySupport.hpp"
#include "OneDRuntime.hpp"

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <map>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unistd.h>

namespace {

struct TemporaryDirectory {
	std::filesystem::path path;
	~TemporaryDirectory()
	{
		if (std::getenv("TUBULARFLOWIGA_KEEP_TEST_OUTPUT")) return;
		std::error_code error;
		std::filesystem::remove_all(path, error);
	}
};

} // namespace

int Run()
{
	const auto path = std::filesystem::temp_directory_path()/"tubularflowiga-c0-square-duct-test.ntiga";
	for (int axial_elements : {1, 2})
		for (std::uint32_t ranks : {std::uint32_t{1}, std::uint32_t{2}}) {
			iga::test::WriteC0SquareDuctDatabase(path, 2, ranks, axial_elements);
			iga::Database database(path.string());
			assert(database.header().elements == static_cast<std::uint64_t>(4*axial_elements));
			assert(database.header().nodes == static_cast<std::uint64_t>((3*axial_elements+1)*7*7));
			for (std::uint32_t rank = 0; rank < ranks; ++rank) {
				assert(database.RequiredElementIndices(static_cast<int>(rank)).size() == static_cast<std::size_t>(4*axial_elements));
				assert(database.LoadOwned(static_cast<int>(rank)).size() == static_cast<std::size_t>(ranks == 1 ? 4*axial_elements : 2*axial_elements));
			}
			auto element_index = [](int ex, int ez, int ey) { return (ex*2+ez)*2+ey; };
			for (int ex = 0; ex < axial_elements; ++ex)
				for (int ez = 0; ez < 2; ++ez)
					for (int ey = 0; ey < 2; ++ey) {
						const auto element = database.Load(static_cast<std::uint64_t>(element_index(ex, ez, ey)));
						const auto& labels = element.boundary_labels;
						assert(labels[0] == (ez == 0 ? 0 : -1));
						assert(labels[1] == (ey == 0 ? 0 : -1));
						assert(labels[2] == (ex+1 == axial_elements ? 2 : -1));
						assert(labels[3] == (ey == 1 ? 0 : -1));
						assert(labels[4] == (ex == 0 ? 1 : -1));
						assert(labels[5] == (ez == 1 ? 0 : -1));
						assert(element.connectivity.size() == 64);
						for (int k = 0; k < 4; ++k)
							for (int j = 0; j < 4; ++j)
								for (int i = 0; i < 4; ++i) {
									const auto& point = element.bezier_points[static_cast<std::size_t>(i+4*j+16*k)];
									assert(std::abs(point[0]-(ex+i/3.0)/axial_elements) < 1.0e-15);
									assert(std::abs(point[1]-(ey+j/3.0)/2.0) < 1.0e-15);
									assert(std::abs(point[2]-(ez+k/3.0)/2.0) < 1.0e-15);
								}
						const auto& origin = element.bezier_points.front();
						const auto& x_point = element.bezier_points[1];
						const auto& y_point = element.bezier_points[4];
						const auto& z_point = element.bezier_points[16];
						const double determinant = (x_point[0]-origin[0])*(y_point[1]-origin[1])*(z_point[2]-origin[2]);
						assert(determinant > 0.0);
					}
			for (int ex = 0; ex < axial_elements; ++ex)
				for (int ez = 0; ez < 2; ++ez)
					for (int ey = 0; ey < 2; ++ey) {
						const auto current = database.Load(static_cast<std::uint64_t>(element_index(ex, ez, ey)));
						if (ex+1 < axial_elements) {
							const auto next = database.Load(static_cast<std::uint64_t>(element_index(ex+1, ez, ey)));
							for (int k = 0; k < 4; ++k) for (int j = 0; j < 4; ++j)
								assert(current.connectivity[static_cast<std::size_t>(3+4*j+16*k)] == next.connectivity[static_cast<std::size_t>(4*j+16*k)]);
						}
						if (ey+1 < 2) {
							const auto next = database.Load(static_cast<std::uint64_t>(element_index(ex, ez, ey+1)));
							for (int k = 0; k < 4; ++k) for (int i = 0; i < 4; ++i)
								assert(current.connectivity[static_cast<std::size_t>(i+12+16*k)] == next.connectivity[static_cast<std::size_t>(i+16*k)]);
						}
						if (ez+1 < 2) {
							const auto next = database.Load(static_cast<std::uint64_t>(element_index(ex, ez+1, ey)));
							for (int j = 0; j < 4; ++j) for (int i = 0; i < 4; ++i)
								assert(current.connectivity[static_cast<std::size_t>(i+4*j+48)] == next.connectivity[static_cast<std::size_t>(i+4*j)]);
						}
					}
		}
	assert(std::abs(iga::test::kSquareDuctResistanceCoefficient-28.45415376956191) < 1.0e-14);
	assert(std::abs(iga::test::kEquivalentSquareDuctHydraulicLengthM-1.1321548059806663) < 1.0e-14);
	const double mean = iga::test::SquareDuctSeriesMean();
	assert(mean > 0.0 && std::isfinite(mean));
	assert(std::abs(iga::test::SquareDuctNormalizedReferenceProfile(0.0, 0.5)) < 1.0e-13);
	assert(std::abs(iga::test::SquareDuctNormalizedReferenceProfile(0.5, 0.0)) < 1.0e-13);
	assert(iga::test::SquareDuctNormalizedReferenceProfile(0.5, 0.5) > 1.0);
	assert(std::abs(iga::test::RefinementObservedOrder(0.08, 0.02)-2.0) < 1.0e-14);
	assert(std::abs(iga::test::WrapPhaseRadians(3.0*3.14159265358979323846)+3.14159265358979323846) < 1.0e-14);
	assert(std::abs(iga::test::HarmonicPhaseDifferenceRadians(-3.13, 3.13)-0.02318530717958645) < 1.0e-12);
	for (int degree = 0; degree <= 3; ++degree) {
		std::array<double, 4> coefficients{};
		for (int row = 0; row < 4; ++row) {
			const auto inverse = iga::test::CubicBernsteinCollocationInverseRow(row);
			for (int sample = 0; sample < 4; ++sample) coefficients[static_cast<std::size_t>(row)] += inverse[static_cast<std::size_t>(sample)]*std::pow(sample/3.0, degree);
		}
		for (double t : {0.0, 0.17, 0.5, 0.91, 1.0}) assert(std::abs(iga::test::CubicBernsteinValue(coefficients, t)-std::pow(t, degree)) < 1.0e-13);
	}
	for (int n : {2, 4})
		for (int j = 0; j <= 3*n; ++j)
			for (int k = 0; k <= 3*n; ++k) {
				const double coefficient = iga::test::SquareDuctBernsteinCoefficient(n, j, k);
				if (j == 0 || k == 0 || j == 3*n || k == 3*n) assert(coefficient == 0.0);
				assert(std::abs(coefficient-iga::test::SquareDuctBernsteinCoefficient(n, k, j)) < 64.0*std::numeric_limits<double>::epsilon()*std::max(1.0, std::abs(coefficient)));
			}
	// The tensor-product collocation inverse must reproduce every polynomial in
	// Q_3, independently of the square-duct profile used by the fixture.
	for (int exponent_y = 0; exponent_y <= 3; ++exponent_y)
		for (int exponent_z = 0; exponent_z <= 3; ++exponent_z)
			for (int patch_y = 0; patch_y < 2; ++patch_y)
				for (int patch_z = 0; patch_z < 2; ++patch_z) {
					std::array<std::array<double, 4>, 4> coefficients{};
					for (int j = 0; j < 4; ++j)
						for (int k = 0; k < 4; ++k)
							coefficients[static_cast<std::size_t>(j)][static_cast<std::size_t>(k)] =
								iga::test::CubicTensorBernsteinCoefficient([&](int sample_y, int sample_z) {
									return std::pow((patch_y+sample_y/3.0)/2.0, exponent_y)
										*std::pow((patch_z+sample_z/3.0)/2.0, exponent_z);
								}, j, k);
					for (double local_y : {0.0, 0.13, 0.5, 0.87, 1.0})
						for (double local_z : {0.0, 0.21, 0.5, 0.79, 1.0}) {
							const double expected = std::pow((patch_y+local_y)/2.0, exponent_y)
								*std::pow((patch_z+local_z)/2.0, exponent_z);
							assert(std::abs(iga::test::CubicTensorBernsteinValue(coefficients, local_y, local_z)-expected) < 2.0e-13);
						}
				}
	for (int n : {2, 4}) {
		auto patch_choices = [n](int global) {
			std::vector<std::pair<int, int>> choices;
			if (global > 0 && global < 3*n && global%3 == 0) {
				choices.emplace_back(global/3-1, 3);
				choices.emplace_back(global/3, 0);
			} else {
				const int patch = std::min(global/3, n-1);
				choices.emplace_back(patch, global-3*patch);
			}
			return choices;
		};
		for (int global_y = 0; global_y <= 3*n; ++global_y)
			for (int global_z = 0; global_z <= 3*n; ++global_z) {
				std::vector<double> proposals;
				for (const auto& y : patch_choices(global_y))
					for (const auto& z : patch_choices(global_z))
						proposals.push_back(iga::test::SquareDuctBernsteinPatchCoefficient(n, y.first, z.first,
							y.second, z.second));
				for (const double proposal : proposals)
					assert(std::abs(proposal-proposals.front()) <= 64.0*std::numeric_limits<double>::epsilon()
						*std::max({1.0, std::abs(proposal), std::abs(proposals.front())}));
				const bool wall = global_y == 0 || global_z == 0 || global_y == 3*n || global_z == 3*n;
				if (wall) for (const double proposal : proposals) assert(proposal == 0.0);
			}
	}
	std::array<double, 4> corrected_l2{}, sampled_l2{};
	for (std::size_t index = 0; index < 4; ++index) {
		const int n = 1 << static_cast<int>(index);
		corrected_l2[index] = iga::test::SquareDuctProfileL2ErrorGauss(n, iga::test::SquareDuctProfileMode::CorrectedBernstein);
		sampled_l2[index] = iga::test::SquareDuctProfileL2ErrorGauss(n, iga::test::SquareDuctProfileMode::SampledBaseline);
		const double coefficient_integral = [&] {
			double sum = 0.0;
			for (int patch_y = 0; patch_y < n; ++patch_y)
				for (int patch_z = 0; patch_z < n; ++patch_z)
					for (int j = 0; j < 4; ++j)
						for (int k = 0; k < 4; ++k)
							sum += iga::test::SquareDuctBernsteinPatchCoefficient(n, patch_y, patch_z, j, k);
			return sum/(16.0*n*n);
		}();
		const double gauss_integral = iga::test::SquareDuctProfileIntegralGauss(n, iga::test::SquareDuctProfileMode::CorrectedBernstein);
		assert(coefficient_integral > 0.0 && gauss_integral > 0.0);
		assert(std::abs(coefficient_integral-gauss_integral) < 5.0e-13);
		assert(corrected_l2[index] < sampled_l2[index]);
		if (index > 0) assert(corrected_l2[index] < corrected_l2[index-1]);
	}
	// The profile is smooth away from exact walls; a conservative L2 order-two
	// threshold confirms the intended cubic interpolation without overfitting a
	// finite series truncation or quadrature roundoff plateau.
	assert(iga::test::RefinementObservedOrder(corrected_l2[1], corrected_l2[2]) >= 2.0);
	assert(iga::test::RefinementObservedOrder(corrected_l2[2], corrected_l2[3]) >= 2.0);
	const double mean_64 = iga::test::SquareDuctSeriesMeanGauss(64, 64);
	const double mean_128 = iga::test::SquareDuctSeriesMeanGauss(128, 64);
	const double high_accuracy_mean = iga::test::SquareDuctSeriesMeanGauss(256, 64);
	const double resistance_64 = iga::test::kPi*iga::test::kPi*iga::test::kPi/(4.0*mean_64);
	const double resistance_128 = iga::test::kPi*iga::test::kPi*iga::test::kPi/(4.0*mean_128);
	const double resistance_from_mean = iga::test::kPi*iga::test::kPi*iga::test::kPi/(4.0*high_accuracy_mean);
	assert(std::abs(resistance_128-iga::test::kSquareDuctResistanceCoefficient)
		< std::abs(resistance_64-iga::test::kSquareDuctResistanceCoefficient));
	assert(std::abs(resistance_from_mean-iga::test::kSquareDuctResistanceCoefficient)
		< std::abs(resistance_128-iga::test::kSquareDuctResistanceCoefficient));
	assert(std::abs(resistance_from_mean-iga::test::kSquareDuctResistanceCoefficient)
		< 2.0e-9*iga::test::kSquareDuctResistanceCoefficient);
	const auto fixture_mesh = std::filesystem::temp_directory_path()/"tubularflowiga-c0-square-duct-wall-trace.vtk";
	iga::test::WriteC0SquareDuctDatabase(path, 2, 1, 2);
	iga::test::WriteC0SquareDuctControlMesh(fixture_mesh, 2, 2);
	iga::Database wall_database(path.string());
	const auto wall_mesh = iga::ReadLabeledHexMesh(fixture_mesh.string(), wall_database.header().nodes, wall_database.header().elements);
	const auto wall_trace = iga::WallTraceBasis(wall_database, wall_mesh, 0);
	std::set<std::int32_t> expected_wall_trace;
	for (int k = 0; k < 7; ++k)
		for (int j = 0; j < 7; ++j)
			for (int i = 0; i < 7; ++i)
				if (j == 0 || j == 6 || k == 0 || k == 6) expected_wall_trace.insert(i+7*j+49*k);
	assert(wall_trace.size() == 168 && wall_trace == expected_wall_trace);
	assert(iga::ExternalFaces(wall_mesh).size() == 2*2*2+4*2*2);
	std::filesystem::remove(fixture_mesh);
	std::cout << std::setprecision(17) << "fixture_profile_l2 corrected="
		<< corrected_l2[0] << ',' << corrected_l2[1] << ',' << corrected_l2[2] << ',' << corrected_l2[3]
		<< " sampled=" << sampled_l2[0] << ',' << sampled_l2[1] << ',' << sampled_l2[2] << ',' << sampled_l2[3]
		<< " orders=" << iga::test::RefinementObservedOrder(corrected_l2[1], corrected_l2[2])
		<< ',' << iga::test::RefinementObservedOrder(corrected_l2[2], corrected_l2[3])
		<< " resistance_64_128_256=" << resistance_64 << ',' << resistance_128 << ',' << resistance_from_mean
		<< " wall_trace_basis=" << wall_trace.size() << '\n';
	for (int n : {1, 2, 4}) {
		const int bases = iga::test::OpenUniformCubicBasisCount(n);
		for (int element = 0; element < n; ++element)
			for (double local : {0.0, 0.17, 0.5, 0.83, 1.0}) {
				const double coordinate = local == 1.0 && element+1 < n
					? std::nextafter((element+1.0)/n, -std::numeric_limits<double>::infinity())
					: (element+local)/n;
				double partition = 0.0;
				for (int basis = 0; basis < bases; ++basis) {
					const double extracted = iga::test::OpenUniformCubicExtractedBasisValue(n, element, basis, local);
					const double direct = iga::test::OpenUniformCubicBasisValue(n, basis, coordinate);
					assert(extracted >= -1.0e-13);
					assert(std::abs(extracted-direct) < 2.0e-12);
					partition += extracted;
				}
				assert(std::abs(partition-1.0) < 2.0e-12);
			}
		for (int interface = 1; interface < n; ++interface)
			for (int basis = 0; basis < bases; ++basis)
				for (int derivative = 0; derivative <= 2; ++derivative)
					assert(std::abs(iga::test::OpenUniformCubicExtractedBasisDerivative(n, interface-1, basis, 1.0, derivative)
						-iga::test::OpenUniformCubicExtractedBasisDerivative(n, interface, basis, 0.0, derivative)) < 2.0e-10);
		for (double coordinate : {0.0, 0.07, 0.31, 0.69, 1.0}) {
			double partition = 0.0, linear = 0.0;
			for (int basis = 0; basis < bases; ++basis) {
				const double value = iga::test::OpenUniformCubicBasisValue(n, basis, coordinate);
				partition += value;
				linear += iga::test::OpenUniformCubicGreville(n, basis)*value;
			}
			assert(std::abs(partition-1.0) < 2.0e-12 && std::abs(linear-coordinate) < 2.0e-12);
		}
	}
	{
		const int n = 2, bases = iga::test::OpenUniformCubicBasisCount(n);
		const auto extraction = iga::test::OpenUniformCubicExtraction(n, 1);
		auto nonsymmetric_control = [bases](int x, int y, int z) { return 100.0*x+10.0*y+z; };
		for (int k = 0; k < 4; ++k)
			for (int j = 0; j < 4; ++j)
				for (int i = 0; i < 4; ++i) {
					double extracted = 0.0;
					for (int c = 0; c < 4; ++c)
						for (int b = 0; b < 4; ++b)
							for (int a = 0; a < 4; ++a)
								extracted += extraction[static_cast<std::size_t>(a)][static_cast<std::size_t>(i)]
									*extraction[static_cast<std::size_t>(b)][static_cast<std::size_t>(j)]
									*extraction[static_cast<std::size_t>(c)][static_cast<std::size_t>(k)]
									*nonsymmetric_control(1+a, 1+b, 1+c);
					double x = 0.0, y = 0.0, z = 0.0;
					for (int a = 0; a < 4; ++a) x += extraction[static_cast<std::size_t>(a)][static_cast<std::size_t>(i)]*(1+a);
					for (int b = 0; b < 4; ++b) y += extraction[static_cast<std::size_t>(b)][static_cast<std::size_t>(j)]*(1+b);
					for (int c = 0; c < 4; ++c) z += extraction[static_cast<std::size_t>(c)][static_cast<std::size_t>(k)]*(1+c);
					const double expected = 100.0*x+10.0*y+z;
					assert(std::abs(extracted-expected) < 2.0e-12);
				}
	}
	std::array<double, 4> c2_l2{}, c2_energy_resistance{}, c2_ritz_resistance{}, c2_weak_residual{};
	for (std::size_t index = 0; index < c2_l2.size(); ++index) {
		const int n = 1 << static_cast<int>(index);
		const auto profile = iga::test::BuildC2SquareDuctProfile(n);
		c2_l2[index] = iga::test::C2SquareDuctProfileL2ErrorGauss(profile);
		c2_energy_resistance[index] = iga::test::C2SquareDuctProfileEnergyResistance(profile);
		const auto ritz = iga::test::AuditC2SquareDuctRitzProfile(profile);
		c2_ritz_resistance[index] = ritz.ritz_resistance;
		c2_weak_residual[index] = ritz.normalized_weak_residual;
		assert(std::abs(ritz.profile_energy_resistance-c2_energy_resistance[index])
			< 1.0e-11*c2_energy_resistance[index]);
		assert(c2_energy_resistance[index]+1.0e-11 >= c2_ritz_resistance[index]);
		assert(c2_ritz_resistance[index]+1.0e-11 >= iga::test::kSquareDuctResistanceCoefficient);
		assert(c2_l2[index] > 0.0);
		const double integral = iga::test::GaussLegendre8UnitSquare([&](double y, double z) {
			const double value = iga::test::EvaluateC2SquareDuctProfile(profile, y, z);
			assert(value >= -1.0e-12);
			return value;
		}, n);
		assert(integral > 0.0);
		if (index > 0) assert(c2_l2[index] < c2_l2[index-1]);
	}
	assert(iga::test::RefinementObservedOrder(c2_l2[1], c2_l2[2]) >= 2.0);
	assert(iga::test::RefinementObservedOrder(c2_l2[2], c2_l2[3]) >= 2.0);
	const auto c2_database_path = std::filesystem::temp_directory_path()/"tubularflowiga-c2-square-duct-test.ntiga";
	const auto c2_mesh_path = std::filesystem::temp_directory_path()/"tubularflowiga-c2-square-duct-test.vtk";
	for (int transverse_elements : {1, 2})
		for (std::uint32_t ranks : {std::uint32_t{1}, std::uint32_t{2}})
			for (double length_m : {1.0, 1.5}) {
				const int axial_elements = 2*transverse_elements;
				const int bases_x = iga::test::OpenUniformCubicBasisCount(axial_elements);
				const int bases_yz = iga::test::OpenUniformCubicBasisCount(transverse_elements);
				const std::uint64_t elements = static_cast<std::uint64_t>(axial_elements)
					*transverse_elements*transverse_elements;
				iga::test::WriteC2SquareDuctDatabase(c2_database_path, transverse_elements,
					ranks, axial_elements, length_m);
				iga::Database audit(c2_database_path.string());
				assert(audit.header().elements == elements);
				assert(audit.header().nodes == static_cast<std::uint64_t>(bases_x)*bases_yz*bases_yz);
				for (std::uint32_t rank = 0; rank < ranks; ++rank) {
					std::size_t expected_owned = 0;
					for (std::uint64_t element = 0; element < elements; ++element)
						if (element%ranks == rank) ++expected_owned;
					assert(audit.LoadOwned(static_cast<int>(rank)).size() == expected_owned);
					assert(audit.RequiredElementIndices(static_cast<int>(rank)).size()
						== static_cast<std::size_t>(elements));
				}
				auto element_index = [transverse_elements](int ex, int ez, int ey) {
					return static_cast<std::uint64_t>((ex*transverse_elements+ez)*transverse_elements+ey);
				};
				for (int ex = 0; ex < axial_elements; ++ex)
					for (int ez = 0; ez < transverse_elements; ++ez)
						for (int ey = 0; ey < transverse_elements; ++ey) {
							const auto element = audit.Load(element_index(ex, ez, ey));
							assert(element.boundary_labels[0] == (ez == 0 ? 0 : -1));
							assert(element.boundary_labels[1] == (ey == 0 ? 0 : -1));
							assert(element.boundary_labels[2] == (ex+1 == axial_elements ? 2 : -1));
							assert(element.boundary_labels[3] == (ey+1 == transverse_elements ? 0 : -1));
							assert(element.boundary_labels[4] == (ex == 0 ? 1 : -1));
							assert(element.boundary_labels[5] == (ez+1 == transverse_elements ? 0 : -1));
							for (int k = 0; k < 4; ++k)
								for (int j = 0; j < 4; ++j)
									for (int i = 0; i < 4; ++i) {
										const auto& point = element.bezier_points[static_cast<std::size_t>(i+4*j+16*k)];
										assert(std::abs(point[0]-length_m*(ex+i/3.0)/axial_elements) < 3.0e-13);
										assert(std::abs(point[1]-(ey+j/3.0)/transverse_elements) < 3.0e-13);
										assert(std::abs(point[2]-(ez+k/3.0)/transverse_elements) < 3.0e-13);
									}
							const auto& origin = element.bezier_points[0];
							const double determinant = (element.bezier_points[1][0]-origin[0])
								*(element.bezier_points[4][1]-origin[1])
								*(element.bezier_points[16][2]-origin[2]);
							assert(determinant > 0.0);
							if (ex+1 < axial_elements) {
								const auto next = audit.Load(element_index(ex+1, ez, ey));
								for (int k = 0; k < 4; ++k) for (int j = 0; j < 4; ++j) for (int i = 1; i < 4; ++i)
									assert(element.connectivity[static_cast<std::size_t>(i+4*j+16*k)]
										== next.connectivity[static_cast<std::size_t>(i-1+4*j+16*k)]);
							}
							if (ey+1 < transverse_elements) {
								const auto next = audit.Load(element_index(ex, ez, ey+1));
								for (int k = 0; k < 4; ++k) for (int j = 1; j < 4; ++j) for (int i = 0; i < 4; ++i)
									assert(element.connectivity[static_cast<std::size_t>(i+4*j+16*k)]
										== next.connectivity[static_cast<std::size_t>(i+4*(j-1)+16*k)]);
							}
							if (ez+1 < transverse_elements) {
								const auto next = audit.Load(element_index(ex, ez+1, ey));
								for (int k = 1; k < 4; ++k) for (int j = 0; j < 4; ++j) for (int i = 0; i < 4; ++i)
									assert(element.connectivity[static_cast<std::size_t>(i+4*j+16*k)]
										== next.connectivity[static_cast<std::size_t>(i+4*j+16*(k-1))]);
							}
						}
			}
	iga::test::WriteC2SquareDuctDatabase(c2_database_path, 2, 1, 2);
	iga::test::WriteC2SquareDuctControlMesh(c2_mesh_path, 2, 2);
	iga::Database c2_database(c2_database_path.string());
	assert(c2_database.header().elements == 8 && c2_database.header().nodes == 125);
	const auto c2_mesh = iga::ReadLabeledHexMesh(c2_mesh_path.string(), c2_database.header().nodes, c2_database.header().elements);
	assert(iga::ExternalFaces(c2_mesh).size() == 24);
	const auto c2_wall_trace = iga::WallTraceBasis(c2_database, c2_mesh, 0);
	std::set<std::int32_t> expected_c2_wall_trace;
	int c2_inlet_nonwall = 0;
	for (int k = 0; k < 5; ++k)
		for (int j = 0; j < 5; ++j)
			for (int i = 0; i < 5; ++i) {
				const std::int32_t id = i+5*(j+5*k);
				if (j == 0 || j == 4 || k == 0 || k == 4) expected_c2_wall_trace.insert(id);
				if (i == 0 && j > 0 && j < 4 && k > 0 && k < 4) ++c2_inlet_nonwall;
			}
	assert(c2_wall_trace == expected_c2_wall_trace && c2_wall_trace.size() == 80 && c2_inlet_nonwall == 9);
	bool found_sparse_row = false;
	for (std::uint64_t element = 0; element < c2_database.header().elements; ++element) {
		const auto c2_element = c2_database.Load(element);
		assert(c2_element.connectivity.size() == 64 && c2_element.extraction.size() == 64);
		for (const auto& row : c2_element.extraction) {
			int nonzeros = 0;
			for (double value : row) if (value != 0.0) ++nonzeros;
			assert(nonzeros > 0 && nonzeros <= 64);
			found_sparse_row = found_sparse_row || nonzeros < 64;
		}
	}
	assert(found_sparse_row);
	std::filesystem::remove(c2_database_path);
	std::filesystem::remove(c2_mesh_path);
	std::cout << std::setprecision(17) << "fixture_c2_profile_l2=" << c2_l2[0] << ',' << c2_l2[1] << ','
		<< c2_l2[2] << ',' << c2_l2[3] << " orders="
		<< iga::test::RefinementObservedOrder(c2_l2[1], c2_l2[2]) << ','
		<< iga::test::RefinementObservedOrder(c2_l2[2], c2_l2[3]) << " energy_resistance="
		<< c2_energy_resistance[0] << ',' << c2_energy_resistance[1] << ','
		<< c2_energy_resistance[2] << ',' << c2_energy_resistance[3] << " ritz_resistance="
		<< c2_ritz_resistance[0] << ',' << c2_ritz_resistance[1] << ','
		<< c2_ritz_resistance[2] << ',' << c2_ritz_resistance[3] << " weak_residual="
		<< c2_weak_residual[0] << ',' << c2_weak_residual[1] << ','
		<< c2_weak_residual[2] << ',' << c2_weak_residual[3]
		<< " wall_trace_basis=" << c2_wall_trace.size() << '\n';
	if (std::getenv("TUBULARFLOWIGA_FIXTURE_ONLY")) return 0;
	const auto mesh = std::filesystem::temp_directory_path()/"tubularflowiga-c0-square-duct-test.vtk";
	const auto velocity = std::filesystem::temp_directory_path()/"tubularflowiga-c0-square-duct-test-velocity.txt";
	iga::test::WriteC0SquareDuctControlMesh(mesh, 2);
	iga::test::WriteSquareDuctReferenceVelocity(velocity, 2, 1.0e-3);
	assert(std::filesystem::file_size(mesh) > 0 && std::filesystem::file_size(velocity) > 0);
	const auto case_root = std::filesystem::temp_directory_path()/"tubularflowiga-c0-square-duct-case";
	iga::test::WriteSquareDuctThreeDCase(case_root/"three_d", 2, 0.01, 3, 1.0e-6, 1.0, 1.0e-3);
	iga::test::WriteRigidOneDStraightCase(case_root/"one_d", 1.0, std::sqrt(1.0/3.14159265358979323846), 0.01, 3, 1.0e-6, 1.0, 1.0e-3);
	iga::test::WriteRigidOneDStraightCase(case_root/"one_d_sinusoid", 1.0, std::sqrt(1.0/3.14159265358979323846), 0.01, 3, 1.0e-6, 1.0, 1.0e-3, "sinusoid", 1.0);
	assert(std::filesystem::exists(case_root/"three_d/simulation_config.json"));
	assert(std::filesystem::exists(case_root/"one_d/simulation_config.json"));
	std::ifstream sinusoid_input(case_root/"one_d_sinusoid/simulation_config.json");
	const std::string sinusoid_text((std::istreambuf_iterator<char>(sinusoid_input)), std::istreambuf_iterator<char>());
	const auto sinusoid = sinusoid_text.find("\"kind\":\"sinusoid\"");
	const auto equation_systems = sinusoid_text.find("\"equation_systems\"");
	assert(sinusoid != std::string::npos && equation_systems != std::string::npos);
	assert(sinusoid_text.substr(sinusoid, equation_systems-sinusoid).find("\"value\"") == std::string::npos);
	// A single serial runtime is the legitimate all-1D reference: all three
	// portions share the area-matched circular radius, and the middle segment
	// has the square duct's hydraulic-equivalent, not geometric, length.
	const double pi = 3.14159265358979323846;
	const double radius = std::sqrt(1.0/pi);
	const double middle = iga::test::kEquivalentSquareDuctHydraulicLengthM;
	iga::test::WriteRigidOneDStraightCase(case_root/"all_one_d", 1.0, radius, 0.01, 1, 1.0e-6, 1.0, 1.0e-3, "constant", 1.0, 4);
	{
		std::ofstream tree(case_root/"all_one_d/tree.swc", std::ios::trunc);
		tree << std::setprecision(17) << "1 2 0 0 0 " << radius << " -1\n"
			<< "2 2 1 0 0 " << radius << " 1\n"
			<< "3 2 " << 1.0+middle << " 0 0 " << radius << " 2\n"
			<< "4 2 " << 2.0+middle << " 0 0 " << radius << " 3\n";
	}
	std::ifstream reference_config_input(case_root/"all_one_d/simulation_config.json");
	const std::string reference_config((std::istreambuf_iterator<char>(reference_config_input)), std::istreambuf_iterator<char>());
	const auto configuration = iga::ParseOneDConfiguration(reference_config);
	const auto flow = configuration.flow_systems.front();
	const auto network = iga::ReadOneDNetwork(case_root/"all_one_d/tree.swc", 1.0, 1, flow.dynamic_viscosity);
	iga::OneDFlowRuntime reference(configuration, flow, network, iga::ResolveOneDInlet(configuration), case_root/"all_one_d");
	reference.InitializeOpenLoop(1.0e-3);
	reference.BeginStep(0.0, 0.01);
	reference.SetOpenLoopInlet(reference.OpenLoopInlet(0.01, 1.0e-3));
	reference.SolveTrial();
	const auto root = reference.GetPortState("root");
	const auto terminal = reference.GetPortState("outlet:4");
	const double expected_drop = (2.0*8.0*pi+iga::test::kSquareDuctResistanceCoefficient)*1.0e-3;
	auto number = [](double value) { std::ostringstream output; output << std::setprecision(17) << value; return output.str(); };
	auto quote = [](const std::filesystem::path& value) { return "'"+value.string()+"'"; };
	const double observed_drop = *root.mean_pressure_pa-*terminal.mean_pressure_pa;
	assert(std::abs(observed_drop-expected_drop) <= 1.0e-13*std::max(1.0, expected_drop));
	assert(*root.outward_flow_m3_s < 0.0 && *terminal.outward_flow_m3_s > 0.0);
	assert(std::abs(*root.outward_flow_m3_s+*terminal.outward_flow_m3_s) < 1.0e-15);
	std::cout << std::setprecision(17) << "all_1d_hydraulic_reference_drop_pa=" << observed_drop
		<< " expected_pa=" << expected_drop << std::endl;
	reference.CommitStep();
	const bool axial_sweep = std::getenv("TUBULARFLOWIGA_AXIAL_ONLY") != nullptr;
	const bool isotropic = std::getenv("TUBULARFLOWIGA_ISOTROPIC_ONLY") != nullptr;
	const bool length_sweep = std::getenv("TUBULARFLOWIGA_LENGTH_ONLY") != nullptr;
	const bool bulk_sweep = std::getenv("TUBULARFLOWIGA_BULK_ONLY") != nullptr
		|| (!axial_sweep && !isotropic && !length_sweep);
	if (axial_sweep || isotropic || length_sweep || bulk_sweep) {
		struct AxialMetric {
			int transverse_elements = 0;
			int axial_elements = 0;
			double length_m = 0.0;
			double three_d_resistance_pa_s_m3 = 0.0;
			double external_resistance_pa_s_m3 = 0.0;
			double three_d_relative_error = 0.0;
			double external_relative_error = 0.0;
			double final_step_resistance_change = 0.0;
			double upstream_pressure_jump_pa = 0.0;
			double downstream_pressure_jump_pa = 0.0;
			double outlet_static_minus_traction_pa = 0.0;
			double inlet_fixed_point_mismatch_pa = 0.0;
			double max_interface = 0.0;
			double max_mass = 0.0;
			double max_wall = 0.0;
			long long total_sweeps = 0;
			long long total_ksp = 0;
		};
		auto parse_csv = [](const std::filesystem::path& csv) {
			std::ifstream input(csv);
			if (!input) throw std::runtime_error("missing axial diagnostic CSV "+csv.string());
			std::string header;
			if (!std::getline(input, header)) throw std::runtime_error("empty axial diagnostic CSV");
			std::vector<std::string> names;
			std::istringstream headings(header);
			for (std::string name; std::getline(headings, name, ',');) names.push_back(name);
			std::vector<std::map<std::string, double>> rows;
			for (std::string line; std::getline(input, line);) {
				if (line.empty()) continue;
				std::istringstream values(line);
				std::map<std::string, double> row;
				for (const auto& name : names) {
					std::string value;
					if (!std::getline(values, value, ',')) throw std::runtime_error("truncated axial diagnostic CSV");
					row.emplace(name, std::stod(value));
				}
				rows.push_back(std::move(row));
			}
			return rows;
		};
		const auto axial_root = std::filesystem::temp_directory_path()/(
			"tubularflowiga-straight-vessel-axial-"+std::to_string(static_cast<long long>(getpid())));
		TemporaryDirectory axial_cleanup{axial_root};
		if (std::getenv("TUBULARFLOWIGA_KEEP_TEST_OUTPUT"))
			std::cout << "retained_test_output=" << axial_root << std::endl;
		std::error_code axial_error;
		std::filesystem::remove_all(axial_root, axial_error);
		std::vector<AxialMetric> metrics;
		struct CaseSpecification { int transverse_elements; int axial_elements; double length_m; };
		std::vector<CaseSpecification> cases;
		if (bulk_sweep) {
			for (int transverse_elements : {2, 4, 8})
				for (double length_m : {1.0, 1.5, 2.0})
					cases.push_back({transverse_elements,
						static_cast<int>(std::lround(length_m*transverse_elements)), length_m});
		} else {
			for (int level_value : {1, 2, 4}) {
				const int axial_elements = length_sweep ? 4*level_value : level_value;
				cases.push_back({isotropic ? axial_elements : 4, axial_elements,
					length_sweep ? static_cast<double>(level_value) : 1.0});
			}
		}
		for (const auto& specification : cases) {
			const int axial_elements = specification.axial_elements;
			const int transverse_elements = specification.transverse_elements;
			const double length_m = specification.length_m;
			const double level_expected_drop = (2.0*8.0*pi
				+iga::test::kSquareDuctResistanceCoefficient*length_m)*1.0e-3;
			const auto level = axial_root/("n"+std::to_string(transverse_elements)
				+"-nx"+std::to_string(axial_elements)+"-L"+number(length_m));
			std::filesystem::create_directories(level/"three_d");
			std::filesystem::create_directories(level/"upstream");
			std::filesystem::create_directories(level/"downstream");
			iga::test::WriteC2SquareDuctDatabase(level/"case.ntiga", transverse_elements, 1, axial_elements, length_m);
			iga::test::WriteC2SquareDuctThreeDCase(level/"three_d", transverse_elements, 0.01, 3, 1.0e-6, 1.0, 1.0e-3, axial_elements, length_m);
			iga::test::WriteRigidOneDStraightCase(level/"upstream", 1.0, radius, 0.01, 3, 1.0e-6, 1.0, 1.0e-3);
			iga::test::WriteRigidOneDStraightCase(level/"downstream", 1.0, radius, 0.01, 3, 1.0e-6, 1.0, 1.0e-3);
			const auto output = level/"output";
			const std::string command = "mpiexec -np 1 ./iga_1d_3d_explicit "+quote(level/"case.ntiga")+" "+quote(level/"three_d")+" "+quote(level/"upstream")+" "+quote(level/"downstream")
				+" --upstream-terminal-node 2 --output-dir "+quote(output)
				+" --coupling-mode strong-aitken --strong-pressure-reference-pa "+number(level_expected_drop)
				+" --strong-pressure-relative-tol 1e-9 --strong-flow-relative-tol 1e-10 --strong-max-iterations 30 --three-d-max-newton 60 -ksp_type preonly -pc_type lu > /dev/null";
			if (std::system(command.c_str()) != 0)
				throw std::runtime_error("axial diagnostic production run failed nx="+std::to_string(axial_elements));
			const auto history = parse_csv(output/"strong_coupling_history.csv");
			const auto iterations = parse_csv(output/"strong_coupling_iterations.csv");
			if (history.size() != 3 || iterations.empty()) throw std::runtime_error("axial diagnostic artifacts are incomplete");
			AxialMetric metric;
			metric.transverse_elements = transverse_elements;
			metric.axial_elements = axial_elements;
			metric.length_m = length_m;
			for (const auto& row : history) {
				for (const auto& value : row) if (!std::isfinite(value.second)) throw std::runtime_error("nonfinite axial diagnostic history");
				metric.max_interface = std::max(metric.max_interface, std::max(
					std::abs(row.at("upstream_three_d_normalized_residual")), std::abs(row.at("three_d_downstream_normalized_residual"))));
				metric.max_mass = std::max(metric.max_mass, std::abs(row.at("three_d_mass_imbalance_m3_s"))/1.0e-3);
				metric.max_wall = std::max(metric.max_wall, std::abs(row.at("three_d_wall_outward_flow_m3_s"))/1.0e-3);
				metric.total_sweeps += static_cast<long long>(row.at("iteration_count"));
				metric.total_ksp += static_cast<long long>(row.at("all_three_d_ksp_iterations"));
			}
			const auto& final_history = history.back();
			metric.three_d_resistance_pa_s_m3 = (final_history.at("three_d_inlet_pressure_pa")-final_history.at("three_d_outlet_pressure_pa"))/1.0e-3;
			metric.external_resistance_pa_s_m3 = final_history.at("external_pressure_drop_pa")/1.0e-3;
			metric.three_d_relative_error = std::abs(metric.three_d_resistance_pa_s_m3
				-iga::test::kSquareDuctResistanceCoefficient*length_m)
				/(iga::test::kSquareDuctResistanceCoefficient*length_m);
			metric.external_relative_error = std::abs(final_history.at("external_pressure_drop_pa")
				-level_expected_drop)/level_expected_drop;
			const auto& previous_history = history[history.size()-2];
			const double previous_resistance = (previous_history.at("three_d_inlet_pressure_pa")
				-previous_history.at("three_d_outlet_pressure_pa"))/1.0e-3;
			metric.final_step_resistance_change = std::abs(metric.three_d_resistance_pa_s_m3-previous_resistance)
				/std::max(std::abs(metric.three_d_resistance_pa_s_m3), 1.0e-30);
			metric.upstream_pressure_jump_pa = final_history.at("upstream_three_d_pressure_jump_pa");
			metric.downstream_pressure_jump_pa = final_history.at("three_d_downstream_pressure_jump_pa");
			const auto final_iteration = std::find_if(iterations.rbegin(), iterations.rend(), [](const auto& row) {
				auto found = row.find("converged"); return found != row.end() && found->second == 1.0;
			});
			if (final_iteration == iterations.rend()) throw std::runtime_error("axial diagnostic has no converged iteration");
			if (final_iteration->at("physical_step") != static_cast<double>(history.size())
				|| final_iteration->at("time_s") != final_history.at("time_s"))
				throw std::runtime_error("axial diagnostic final iteration does not match the final physical step");
			metric.outlet_static_minus_traction_pa = final_history.at("three_d_outlet_pressure_pa")
				-final_iteration->at("applied_three_d_outlet_traction_pressure_pa");
			metric.inlet_fixed_point_mismatch_pa = final_iteration->at("measured_three_d_inlet_pressure_pa")
				-final_iteration->at("applied_upstream_terminal_pressure_pa");
			if (metric.max_interface > 1.0e-10 || metric.max_mass > 1.0e-6 || metric.max_wall > 1.0e-10)
				throw std::runtime_error("axial diagnostic conservation gate failed nx="+std::to_string(axial_elements));
			metrics.push_back(metric);
			std::cout << std::setprecision(17) << (bulk_sweep ? "bulk" : (length_sweep ? "length" : (isotropic ? "isotropic" : "axial")))
				<< " L=" << length_m
				<< " n=" << transverse_elements << " nx=" << axial_elements
				<< " elements=" << axial_elements*transverse_elements*transverse_elements
				<< " nodes=" << (axial_elements+3)*(transverse_elements+3)*(transverse_elements+3)
				<< " three_d_resistance_pa_s_m3=" << metric.three_d_resistance_pa_s_m3
				<< " three_d_relative_error=" << metric.three_d_relative_error
				<< " external_relative_error=" << metric.external_relative_error
				<< " final_step_resistance_change=" << metric.final_step_resistance_change
				<< " upstream_pressure_jump_pa=" << metric.upstream_pressure_jump_pa
				<< " downstream_pressure_jump_pa=" << metric.downstream_pressure_jump_pa
				<< " outlet_static_minus_applied_traction_pa=" << metric.outlet_static_minus_traction_pa
				<< " inlet_fixed_point_mismatch_pa=" << metric.inlet_fixed_point_mismatch_pa
				<< " max_interface=" << metric.max_interface << " max_mass=" << metric.max_mass
				<< " max_wall=" << metric.max_wall << " total_sweeps=" << metric.total_sweeps
				<< " total_three_d_ksp=" << metric.total_ksp << std::endl;
		}
		if (bulk_sweep) {
			if (metrics.size() != 9) throw std::runtime_error("bulk convergence ladder is incomplete");
			std::array<double, 3> three_d_bulk{}, external_bulk{}, three_d_error{}, external_error{};
			std::array<double, 3> slope_disagreement{}, traction_gap_range{};
			for (std::size_t level = 0; level < three_d_bulk.size(); ++level) {
				const auto& short_case = metrics[3*level];
				const auto& middle_case = metrics[3*level+1];
				const auto& long_case = metrics[3*level+2];
				if (short_case.length_m != 1.0 || middle_case.length_m != 1.5 || long_case.length_m != 2.0
					|| short_case.transverse_elements != middle_case.transverse_elements
					|| short_case.transverse_elements != long_case.transverse_elements)
					throw std::runtime_error("bulk convergence case pairing is invalid");
				three_d_bulk[level] = long_case.three_d_resistance_pa_s_m3-short_case.three_d_resistance_pa_s_m3;
				external_bulk[level] = long_case.external_resistance_pa_s_m3-short_case.external_resistance_pa_s_m3;
				three_d_error[level] = std::abs(three_d_bulk[level]-iga::test::kSquareDuctResistanceCoefficient)
					/iga::test::kSquareDuctResistanceCoefficient;
				external_error[level] = std::abs(external_bulk[level]-iga::test::kSquareDuctResistanceCoefficient)
					/iga::test::kSquareDuctResistanceCoefficient;
				const double first_three_d_slope = 2.0*(middle_case.three_d_resistance_pa_s_m3
					-short_case.three_d_resistance_pa_s_m3);
				const double second_three_d_slope = 2.0*(long_case.three_d_resistance_pa_s_m3
					-middle_case.three_d_resistance_pa_s_m3);
				const double first_external_slope = 2.0*(middle_case.external_resistance_pa_s_m3
					-short_case.external_resistance_pa_s_m3);
				const double second_external_slope = 2.0*(long_case.external_resistance_pa_s_m3
					-middle_case.external_resistance_pa_s_m3);
				slope_disagreement[level] = std::max(std::abs(first_three_d_slope-second_three_d_slope),
					std::abs(first_external_slope-second_external_slope))
					/iga::test::kSquareDuctResistanceCoefficient;
				const double minimum_traction_gap = std::min({short_case.outlet_static_minus_traction_pa,
					middle_case.outlet_static_minus_traction_pa, long_case.outlet_static_minus_traction_pa});
				const double maximum_traction_gap = std::max({short_case.outlet_static_minus_traction_pa,
					middle_case.outlet_static_minus_traction_pa, long_case.outlet_static_minus_traction_pa});
				traction_gap_range[level] = maximum_traction_gap-minimum_traction_gap;
				if (short_case.final_step_resistance_change > 1.0e-8
					|| middle_case.final_step_resistance_change > 1.0e-8
					|| long_case.final_step_resistance_change > 1.0e-8)
					throw std::runtime_error("bulk convergence final resistance is not steady n="
						+std::to_string(short_case.transverse_elements));
				std::cout << std::setprecision(17) << "bulk_pair n=" << short_case.transverse_elements
					<< " three_d_coefficient=" << three_d_bulk[level]
					<< " three_d_relative_error=" << three_d_error[level]
					<< " external_coefficient=" << external_bulk[level]
					<< " external_relative_error=" << external_error[level]
					<< " slope_disagreement=" << slope_disagreement[level]
					<< " outlet_traction_gap_range_pa=" << traction_gap_range[level] << std::endl;
			}
			for (std::size_t level = 1; level < three_d_error.size(); ++level)
				if (!(three_d_error[level] < three_d_error[level-1]
					&& external_error[level] < external_error[level-1]
					&& slope_disagreement[level] < slope_disagreement[level-1]
					&& traction_gap_range[level] < traction_gap_range[level-1]))
					throw std::runtime_error("bulk convergence errors are not strictly decreasing");
			if (slope_disagreement.back() > 0.005
				|| traction_gap_range.back() > 0.005*iga::test::kSquareDuctResistanceCoefficient*1.0e-3)
				throw std::runtime_error("bulk convergence fine length-linearity gate failed");
			const double three_d_order_2_4 = iga::test::RefinementObservedOrder(three_d_error[0], three_d_error[1]);
			const double three_d_order_4_8 = iga::test::RefinementObservedOrder(three_d_error[1], three_d_error[2]);
			const double external_order_2_4 = iga::test::RefinementObservedOrder(external_error[0], external_error[1]);
			const double external_order_4_8 = iga::test::RefinementObservedOrder(external_error[1], external_error[2]);
			if (std::min({three_d_order_2_4, three_d_order_4_8, external_order_2_4, external_order_4_8}) < 1.5)
				throw std::runtime_error("bulk convergence observed order is below 1.5");
			if (three_d_error.back() > 0.005 || external_error.back() > 0.005)
				throw std::runtime_error("bulk convergence fine accuracy gate failed");
			for (std::size_t index = 6; index < metrics.size(); ++index) {
				const auto& fine = metrics[index];
				const double expected_direct_drop = (2.0*8.0*pi
					+iga::test::kSquareDuctResistanceCoefficient*fine.length_m)*1.0e-3;
				if (fine.external_relative_error > 0.02)
					throw std::runtime_error("bulk convergence fine direct all-1D pressure-drop gate failed");
				if (std::abs(fine.upstream_pressure_jump_pa) > 1.0e-8*expected_direct_drop
					|| std::abs(fine.inlet_fixed_point_mismatch_pa) > 1.0e-8*expected_direct_drop)
					throw std::runtime_error("bulk convergence fine upstream pressure-continuity gate failed");
				if (std::abs(fine.downstream_pressure_jump_pa) > 0.005*expected_direct_drop
					|| std::abs(fine.outlet_static_minus_traction_pa) > 0.005*expected_direct_drop)
					throw std::runtime_error("bulk convergence fine outlet pressure-definition gate failed");
			}
			const double observed_order = std::log(std::abs(three_d_bulk[1]-three_d_bulk[0])
				/std::abs(three_d_bulk[2]-three_d_bulk[1]))/std::log(2.0);
			if (observed_order < 1.5 || !std::isfinite(observed_order))
				throw std::runtime_error("bulk convergence Richardson order is invalid");
			const double richardson_denominator = std::pow(2.0, observed_order)-1.0;
			const double extrapolated = three_d_bulk[2]+(three_d_bulk[2]-three_d_bulk[1])/richardson_denominator;
			const double fine_gci = 1.25*std::abs(three_d_bulk[2]-three_d_bulk[1])
				/(std::abs(three_d_bulk[2])*richardson_denominator);
			if (fine_gci > 0.01) throw std::runtime_error("bulk convergence fine GCI exceeds one percent");
			std::cout << std::setprecision(17) << "bulk_convergence orders_three_d=" << three_d_order_2_4
				<< ',' << three_d_order_4_8 << " orders_external=" << external_order_2_4 << ',' << external_order_4_8
				<< " richardson_order=" << observed_order << " extrapolated_coefficient=" << extrapolated
				<< " fine_gci=" << fine_gci << std::endl;
			std::filesystem::remove_all(case_root);
			std::filesystem::remove(path);
			return 0;
		}
		if (length_sweep) {
			const double slope_12 = metrics[1].three_d_resistance_pa_s_m3-metrics[0].three_d_resistance_pa_s_m3;
			const double slope_24 = 0.5*(metrics[2].three_d_resistance_pa_s_m3-metrics[1].three_d_resistance_pa_s_m3);
			const double end_correction = metrics[0].three_d_resistance_pa_s_m3-slope_24;
			std::cout << std::setprecision(17) << "length fixed_h=0.25 slope_L1_to_2=" << slope_12
				<< " slope_L2_to_4_per_m=" << slope_24 << " end_correction=" << end_correction
				<< " slope_relative_difference=" << std::abs(slope_12-slope_24)/iga::test::kSquareDuctResistanceCoefficient
				<< " fine_slope_relative_error=" << std::abs(slope_24-iga::test::kSquareDuctResistanceCoefficient)
					/iga::test::kSquareDuctResistanceCoefficient << std::endl;
			std::filesystem::remove_all(case_root);
			std::filesystem::remove(path);
			return 0;
		}
		const double nx1_to_2 = std::abs(metrics[1].three_d_resistance_pa_s_m3-metrics[0].three_d_resistance_pa_s_m3)
			/std::max(std::abs(metrics[1].three_d_resistance_pa_s_m3), 1.0e-30);
		const double nx2_to_4 = std::abs(metrics[2].three_d_resistance_pa_s_m3-metrics[1].three_d_resistance_pa_s_m3)
			/std::max(std::abs(metrics[2].three_d_resistance_pa_s_m3), 1.0e-30);
		std::cout << std::setprecision(17) << (isotropic ? "isotropic" : "axial n=4")
			<< " resistance_relative_change_nx1_to_2=" << nx1_to_2
			<< " nx2_to_4=" << nx2_to_4 << " diagnostic_threshold=0.0025"
			<< " nx2_to_4_below_threshold=" << (nx2_to_4 <= 0.0025 ? 1 : 0)
			<< " (diagnostic only; transverse spatial gates remain unchanged)" << std::endl;
		std::filesystem::remove_all(case_root);
		std::filesystem::remove(path);
		return 0;
	}
	std::filesystem::remove_all(case_root);
	std::filesystem::remove(mesh);
	std::filesystem::remove(velocity);
	std::filesystem::remove(path);
	return 0;
}

int main()
{
	try { return Run(); }
	catch (const std::exception& error) { std::cerr << "straight-vessel convergence test failed: " << error.what() << '\n'; return 1; }
}
