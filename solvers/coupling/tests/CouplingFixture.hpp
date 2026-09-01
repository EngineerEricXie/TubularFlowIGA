#ifndef IGA_COUPLING_TEST_FIXTURE_HPP
#define IGA_COUPLING_TEST_FIXTURE_HPP

#include "IgaDatabase.hpp"

#include <array>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga::test {

constexpr double kSquareDuctResistanceCoefficient = 28.45415376956191;
constexpr double kEquivalentSquareDuctHydraulicLengthM = 1.1321548059806663;
constexpr double kPi = 3.141592653589793238462643383279502884;

enum class SquareDuctProfileMode { SampledBaseline, CorrectedBernstein };

inline double SquareDuctSeriesRaw(double y, double z, int odd_terms)
{
	double value = 0.0;
	for (int index = 0; index < odd_terms; ++index) {
		const int mode = 2*index+1;
		const double a = mode*kPi;
		// This is cosh(a*(z-1/2))/cosh(a/2), rearranged to avoid
		// overflowing cosh for the high-accuracy independent quadrature.
		const double wall_ratio = (std::exp(-a*z)+std::exp(-a*(1.0-z)))/(1.0+std::exp(-a));
		value += std::sin(a*y)*(1.0-wall_ratio)/(mode*mode*mode);
	}
	return value;
}

inline double SquareDuctSeriesUnnormalized(double y, double z, int odd_terms = 64)
{
	if (!(y >= 0.0 && y <= 1.0 && z >= 0.0 && z <= 1.0) || odd_terms < 1)
		throw std::runtime_error("square-duct series requires unit-square coordinates and positive terms");
	if (y == 0.0 || y == 1.0 || z == 0.0 || z == 1.0) return 0.0;
	return 0.5*(SquareDuctSeriesRaw(y, z, odd_terms)+SquareDuctSeriesRaw(z, y, odd_terms));
}

inline double SquareDuctSeriesMean(int odd_terms = 64, int quadrature = 192)
{
	if (quadrature < 1) throw std::runtime_error("square-duct series quadrature must be positive");
	double sum = 0.0;
	for (int j = 0; j < quadrature; ++j)
		for (int k = 0; k < quadrature; ++k)
			sum += SquareDuctSeriesUnnormalized((j+0.5)/quadrature, (k+0.5)/quadrature, odd_terms);
	return sum/(quadrature*quadrature);
}

inline double SquareDuctNormalizedReferenceProfile(double y, double z, int odd_terms = 64)
{
	static const double default_mean = SquareDuctSeriesMean();
	const double mean = odd_terms == 64 ? default_mean : SquareDuctSeriesMean(odd_terms);
	if (!(mean > 0.0) || !std::isfinite(mean)) throw std::runtime_error("invalid square-duct series mean");
	return SquareDuctSeriesUnnormalized(y, z, odd_terms)/mean;
}

inline double RefinementObservedOrder(double coarse_error, double fine_error)
{
	if (!(coarse_error > 0.0) || !(fine_error > 0.0) || !std::isfinite(coarse_error) || !std::isfinite(fine_error))
		throw std::runtime_error("refinement order requires finite positive errors");
	return std::log(coarse_error/fine_error)/std::log(2.0);
}

inline double WrapPhaseRadians(double phase)
{
	if (!std::isfinite(phase)) throw std::runtime_error("phase must be finite");
	const double pi = 3.141592653589793238462643383279502884;
	phase = std::fmod(phase+pi, 2.0*pi);
	if (phase < 0.0) phase += 2.0*pi;
	return phase-pi;
}

inline double HarmonicPhaseDifferenceRadians(double first, double second)
{
	return WrapPhaseRadians(first-second);
}

inline std::array<double, 4> CubicBernsteinCollocationInverseRow(int row)
{
	static constexpr std::array<std::array<double, 4>, 4> inverse{{
		{{1.0, 0.0, 0.0, 0.0}}, {{-5.0/6.0, 3.0, -3.0/2.0, 1.0/3.0}},
		{{1.0/3.0, -3.0/2.0, 3.0, -5.0/6.0}}, {{0.0, 0.0, 0.0, 1.0}}}};
	if (row < 0 || row > 3) throw std::runtime_error("invalid cubic Bernstein collocation row");
	return inverse[static_cast<std::size_t>(row)];
}

template <class Sample>
inline double CubicTensorBernsteinCoefficient(Sample&& sample, int local_y, int local_z)
{
	const auto ay = CubicBernsteinCollocationInverseRow(local_y);
	const auto az = CubicBernsteinCollocationInverseRow(local_z);
	double coefficient = 0.0;
	for (int y = 0; y < 4; ++y)
		for (int z = 0; z < 4; ++z)
			coefficient += ay[static_cast<std::size_t>(y)]*sample(y, z)*az[static_cast<std::size_t>(z)];
	return coefficient;
}

inline double SquareDuctBernsteinPatchCoefficient(int n, int patch_y, int patch_z,
	int local_y, int local_z, int odd_terms = 64)
{
	if (n < 1 || patch_y < 0 || patch_y >= n || patch_z < 0 || patch_z >= n
		|| local_y < 0 || local_y > 3 || local_z < 0 || local_z > 3)
		throw std::runtime_error("invalid square-duct Bernstein patch coefficient index");
	const int global_y = 3*patch_y+local_y;
	const int global_z = 3*patch_z+local_z;
	if (global_y == 0 || global_z == 0 || global_y == 3*n || global_z == 3*n) return 0.0;
	return CubicTensorBernsteinCoefficient([&](int y, int z) {
		return SquareDuctSeriesUnnormalized((patch_y+y/3.0)/n, (patch_z+z/3.0)/n, odd_terms);
	}, local_y, local_z);
}

inline double SquareDuctBernsteinCoefficient(int n, int global_y, int global_z, int odd_terms = 64)
{
	if (n < 1 || global_y < 0 || global_z < 0 || global_y > 3*n || global_z > 3*n)
		throw std::runtime_error("invalid square-duct Bernstein coefficient index");
	if (global_y == 0 || global_z == 0 || global_y == 3*n || global_z == 3*n) return 0.0;
	const int patch_y = std::min(global_y/3, n-1), patch_z = std::min(global_z/3, n-1);
	return SquareDuctBernsteinPatchCoefficient(n, patch_y, patch_z,
		global_y-3*patch_y, global_z-3*patch_z, odd_terms);
}

inline double CubicBernsteinValue(const std::array<double, 4>& coefficients, double t)
{
	const double s = 1.0-t;
	return coefficients[0]*s*s*s+3.0*coefficients[1]*t*s*s+3.0*coefficients[2]*t*t*s+coefficients[3]*t*t*t;
}

inline double CubicTensorBernsteinValue(const std::array<std::array<double, 4>, 4>& coefficients,
	double y, double z)
{
	std::array<double, 4> along_y{};
	for (int k = 0; k < 4; ++k) {
		std::array<double, 4> row{};
		for (int j = 0; j < 4; ++j) row[static_cast<std::size_t>(j)] = coefficients[static_cast<std::size_t>(j)][static_cast<std::size_t>(k)];
		along_y[static_cast<std::size_t>(k)] = CubicBernsteinValue(row, y);
	}
	return CubicBernsteinValue(along_y, z);
}

inline double SquareDuctPiecewiseBernsteinProfile(int n, double y, double z,
	SquareDuctProfileMode mode, int odd_terms = 64)
{
	if (n < 1 || !(y >= 0.0 && y <= 1.0) || !(z >= 0.0 && z <= 1.0))
		throw std::runtime_error("square-duct Bernstein profile requires unit-square coordinates and positive refinement");
	if (y == 0.0 || z == 0.0 || y == 1.0 || z == 1.0) return 0.0;
	const int patch_y = std::min(static_cast<int>(std::floor(n*y)), n-1);
	const int patch_z = std::min(static_cast<int>(std::floor(n*z)), n-1);
	const double local_y = n*y-patch_y, local_z = n*z-patch_z;
	std::array<std::array<double, 4>, 4> coefficients{};
	for (int j = 0; j < 4; ++j)
		for (int k = 0; k < 4; ++k)
			coefficients[static_cast<std::size_t>(j)][static_cast<std::size_t>(k)] =
				mode == SquareDuctProfileMode::CorrectedBernstein
					? SquareDuctBernsteinPatchCoefficient(n, patch_y, patch_z, j, k, odd_terms)
					: SquareDuctSeriesUnnormalized((patch_y+j/3.0)/n, (patch_z+k/3.0)/n, odd_terms);
	return CubicTensorBernsteinValue(coefficients, local_y, local_z);
}

inline const std::array<double, 8>& GaussLegendre8Nodes()
{
	static const std::array<double, 8> nodes{{-0.96028985649753623168, -0.79666647741362673959,
		-0.52553240991632898582, -0.18343464249564980494, 0.18343464249564980494,
		0.52553240991632898582, 0.79666647741362673959, 0.96028985649753623168}};
	return nodes;
}

inline const std::array<double, 8>& GaussLegendre8Weights()
{
	static const std::array<double, 8> weights{{0.10122853629037625915, 0.22238103445337447054,
		0.31370664587788728734, 0.36268378337836198297, 0.36268378337836198297,
		0.31370664587788728734, 0.22238103445337447054, 0.10122853629037625915}};
	return weights;
}

template <class Function>
inline double GaussLegendre8UnitSquare(Function&& function, int panels_per_axis = 1)
{
	if (panels_per_axis < 1) throw std::runtime_error("Gauss quadrature panels must be positive");
	const auto& nodes = GaussLegendre8Nodes();
	const auto& weights = GaussLegendre8Weights();
	const double scale = 0.5/panels_per_axis;
	double value = 0.0;
	for (int panel_y = 0; panel_y < panels_per_axis; ++panel_y)
		for (int panel_z = 0; panel_z < panels_per_axis; ++panel_z)
			for (int y = 0; y < 8; ++y)
				for (int z = 0; z < 8; ++z) {
					const double physical_y = (panel_y+0.5+0.5*nodes[static_cast<std::size_t>(y)])/panels_per_axis;
					const double physical_z = (panel_z+0.5+0.5*nodes[static_cast<std::size_t>(z)])/panels_per_axis;
					value += scale*scale*weights[static_cast<std::size_t>(y)]*weights[static_cast<std::size_t>(z)]
						*function(physical_y, physical_z);
				}
	return value;
}

inline double SquareDuctProfileIntegralGauss(int n, SquareDuctProfileMode mode, int odd_terms = 64)
{
	return GaussLegendre8UnitSquare([&](double y, double z) {
		return SquareDuctPiecewiseBernsteinProfile(n, y, z, mode, odd_terms);
	}, n);
}

inline double SquareDuctProfileL2ErrorGauss(int n, SquareDuctProfileMode mode, int odd_terms = 64)
{
	const double squared = GaussLegendre8UnitSquare([&](double y, double z) {
		const double difference = SquareDuctPiecewiseBernsteinProfile(n, y, z, mode, odd_terms)
			-SquareDuctSeriesUnnormalized(y, z, odd_terms);
		return difference*difference;
	}, n);
	return std::sqrt(squared);
}

inline double SquareDuctSeriesMeanGauss(int odd_terms, int panels_per_axis)
{
	return GaussLegendre8UnitSquare([&](double y, double z) {
		return SquareDuctSeriesUnnormalized(y, z, odd_terms);
	}, panels_per_axis);
}

inline int OpenUniformCubicBasisCount(int elements)
{
	if (elements < 1) throw std::runtime_error("open-uniform cubic fixture requires positive elements");
	return elements+3;
}

inline std::vector<double> OpenUniformCubicKnots(int elements)
{
	const int bases = OpenUniformCubicBasisCount(elements);
	std::vector<double> knots(static_cast<std::size_t>(bases+4), 1.0);
	for (int index = 0; index < 4; ++index) knots[static_cast<std::size_t>(index)] = 0.0;
	for (int index = 1; index < elements; ++index) knots[static_cast<std::size_t>(index+3)] = index/static_cast<double>(elements);
	return knots;
}

inline double OpenUniformCubicBasisValue(int elements, int basis, double coordinate)
{
	const int bases = OpenUniformCubicBasisCount(elements);
	if (basis < 0 || basis >= bases || !std::isfinite(coordinate) || coordinate < 0.0 || coordinate > 1.0)
		throw std::runtime_error("invalid open-uniform cubic basis evaluation");
	if (coordinate == 1.0) return basis+1 == bases ? 1.0 : 0.0;
	const auto knots = OpenUniformCubicKnots(elements);
	std::vector<double> values(knots.size()-1, 0.0), next(knots.size()-1, 0.0);
	for (std::size_t index = 0; index+1 < knots.size(); ++index)
		if ((coordinate >= knots[index] && coordinate < knots[index+1])
			|| (coordinate == 1.0 && index+2 == knots.size())) values[index] = 1.0;
	for (int degree = 1; degree <= 3; ++degree) {
		std::fill(next.begin(), next.end(), 0.0);
		for (int index = 0; index+degree+1 < static_cast<int>(knots.size()); ++index) {
			const double left_denominator = knots[static_cast<std::size_t>(index+degree)]-knots[static_cast<std::size_t>(index)];
			const double right_denominator = knots[static_cast<std::size_t>(index+degree+1)]-knots[static_cast<std::size_t>(index+1)];
			if (left_denominator > 0.0) next[static_cast<std::size_t>(index)] +=
				(coordinate-knots[static_cast<std::size_t>(index)])/left_denominator*values[static_cast<std::size_t>(index)];
			if (right_denominator > 0.0) next[static_cast<std::size_t>(index)] +=
				(knots[static_cast<std::size_t>(index+degree+1)]-coordinate)/right_denominator*values[static_cast<std::size_t>(index+1)];
		}
		values.swap(next);
	}
	return values[static_cast<std::size_t>(basis)];
}

inline double OpenUniformCubicGreville(int elements, int basis)
{
	const int bases = OpenUniformCubicBasisCount(elements);
	if (basis < 0 || basis >= bases) throw std::runtime_error("invalid open-uniform cubic Greville index");
	const auto knots = OpenUniformCubicKnots(elements);
	return (knots[static_cast<std::size_t>(basis+1)]+knots[static_cast<std::size_t>(basis+2)]
		+knots[static_cast<std::size_t>(basis+3)])/3.0;
}

inline int OpenUniformCubicCellCornerBasis(int endpoint, int elements)
{
	if (endpoint < 0 || endpoint > elements) throw std::runtime_error("invalid open-uniform cubic cell endpoint");
	if (endpoint == 0) return 0;
	if (endpoint == elements) return elements+2;
	return endpoint+1;
}

inline std::array<std::array<double, 4>, 4> OpenUniformCubicExtraction(int elements, int element)
{
	if (element < 0 || element >= elements) throw std::runtime_error("invalid open-uniform cubic element");
	std::array<std::array<double, 4>, 4> extraction{};
	for (int local_basis = 0; local_basis < 4; ++local_basis)
		for (int bernstein = 0; bernstein < 4; ++bernstein) {
			double coefficient = 0.0;
			for (int sample = 0; sample < 4; ++sample) {
				double coordinate = (element+sample/3.0)/elements;
				// The last collocation point is evaluated from the current
				// element's left limit, including the final clamped span.
				if (sample == 3) coordinate = std::nextafter(coordinate, -std::numeric_limits<double>::infinity());
				coefficient += CubicBernsteinCollocationInverseRow(bernstein)[static_cast<std::size_t>(sample)]
					*OpenUniformCubicBasisValue(elements, element+local_basis, coordinate);
			}
			if (!std::isfinite(coefficient) || coefficient < -1.0e-11)
				throw std::runtime_error("open-uniform cubic extraction has an invalid coefficient");
			extraction[static_cast<std::size_t>(local_basis)][static_cast<std::size_t>(bernstein)] =
				std::abs(coefficient) < 1.0e-13 ? 0.0 : coefficient;
		}
	return extraction;
}

inline double OpenUniformCubicExtractedBasisValue(int elements, int element, int global_basis, double local_coordinate)
{
	if (!(local_coordinate >= 0.0 && local_coordinate <= 1.0)) throw std::runtime_error("invalid local cubic coordinate");
	if (global_basis < element || global_basis > element+3) return 0.0;
	const auto extraction = OpenUniformCubicExtraction(elements, element);
	return CubicBernsteinValue(extraction[static_cast<std::size_t>(global_basis-element)], local_coordinate);
}

inline double OpenUniformCubicExtractedBasisDerivative(int elements, int element, int global_basis,
	double local_coordinate, int derivative)
{
	if (derivative < 0 || derivative > 2) throw std::runtime_error("invalid cubic derivative order");
	if (global_basis < element || global_basis > element+3) return 0.0;
	const auto extraction = OpenUniformCubicExtraction(elements, element);
	const auto& c = extraction[static_cast<std::size_t>(global_basis-element)];
	if (derivative == 0) return CubicBernsteinValue(c, local_coordinate);
	std::array<double, 3> first{{3.0*(c[1]-c[0]), 3.0*(c[2]-c[1]), 3.0*(c[3]-c[2])}};
	if (derivative == 1) {
		const double s = 1.0-local_coordinate;
		return elements*(first[0]*s*s+2.0*first[1]*local_coordinate*s+first[2]*local_coordinate*local_coordinate);
	}
	return elements*elements*(2.0*(first[1]-first[0])*(1.0-local_coordinate)
		+2.0*(first[2]-first[1])*local_coordinate);
}

inline std::vector<std::vector<double>> InvertDeterministicMatrix(std::vector<std::vector<double>> matrix)
{
	const std::size_t dimension = matrix.size();
	if (dimension == 0) throw std::runtime_error("cannot invert empty matrix");
	for (const auto& row : matrix) if (row.size() != dimension) throw std::runtime_error("matrix must be square");
	std::vector<std::vector<double>> inverse(dimension, std::vector<double>(dimension, 0.0));
	for (std::size_t row = 0; row < dimension; ++row) inverse[row][row] = 1.0;
	for (std::size_t column = 0; column < dimension; ++column) {
		std::size_t pivot = column;
		for (std::size_t row = column+1; row < dimension; ++row)
			if (std::abs(matrix[row][column]) > std::abs(matrix[pivot][column])) pivot = row;
		if (!std::isfinite(matrix[pivot][column]) || std::abs(matrix[pivot][column]) < 1.0e-12)
			throw std::runtime_error("open-uniform cubic collocation matrix has a small pivot");
		if (pivot != column) { std::swap(matrix[pivot], matrix[column]); std::swap(inverse[pivot], inverse[column]); }
		const double scale = matrix[column][column];
		for (std::size_t entry = 0; entry < dimension; ++entry) { matrix[column][entry] /= scale; inverse[column][entry] /= scale; }
		for (std::size_t row = 0; row < dimension; ++row) if (row != column) {
			const double factor = matrix[row][column];
			for (std::size_t entry = 0; entry < dimension; ++entry) {
				matrix[row][entry] -= factor*matrix[column][entry];
				inverse[row][entry] -= factor*inverse[column][entry];
			}
		}
	}
	for (const auto& row : inverse) for (double value : row)
		if (!std::isfinite(value)) throw std::runtime_error("open-uniform cubic inverse is nonfinite");
	return inverse;
}

struct C2SquareDuctProfile {
	int elements = 0;
	std::vector<double> coefficients;

	double At(int y, int z) const
	{
		const int bases = OpenUniformCubicBasisCount(elements);
		if (y < 0 || y >= bases || z < 0 || z >= bases) throw std::runtime_error("invalid C2 square-duct coefficient index");
		return coefficients[static_cast<std::size_t>(y+bases*z)];
	}
};

inline C2SquareDuctProfile BuildC2SquareDuctProfile(int elements, int odd_terms = 64)
{
	const int bases = OpenUniformCubicBasisCount(elements);
	std::vector<std::vector<double>> collocation(static_cast<std::size_t>(bases), std::vector<double>(static_cast<std::size_t>(bases), 0.0));
	for (int point = 0; point < bases; ++point)
		for (int basis = 0; basis < bases; ++basis)
			collocation[static_cast<std::size_t>(point)][static_cast<std::size_t>(basis)] =
				OpenUniformCubicBasisValue(elements, basis, OpenUniformCubicGreville(elements, point));
	const auto inverse = InvertDeterministicMatrix(collocation);
	std::vector<std::vector<double>> samples(static_cast<std::size_t>(bases), std::vector<double>(static_cast<std::size_t>(bases), 0.0));
	for (int y = 0; y < bases; ++y)
		for (int z = 0; z < bases; ++z)
			samples[static_cast<std::size_t>(y)][static_cast<std::size_t>(z)] =
				SquareDuctSeriesUnnormalized(OpenUniformCubicGreville(elements, y), OpenUniformCubicGreville(elements, z), odd_terms);
	C2SquareDuctProfile profile;
	profile.elements = elements;
	profile.coefficients.assign(static_cast<std::size_t>(bases*bases), 0.0);
	for (int y = 0; y < bases; ++y)
		for (int z = 0; z < bases; ++z) {
			double value = 0.0;
			for (int first = 0; first < bases; ++first)
				for (int second = 0; second < bases; ++second)
					value += inverse[static_cast<std::size_t>(y)][static_cast<std::size_t>(first)]
						*samples[static_cast<std::size_t>(first)][static_cast<std::size_t>(second)]
						*inverse[static_cast<std::size_t>(z)][static_cast<std::size_t>(second)];
			if (y == 0 || z == 0 || y+1 == bases || z+1 == bases) value = 0.0;
			if (!std::isfinite(value) || value < -1.0e-10)
				throw std::runtime_error("C2 square-duct profile has a meaningful negative coefficient");
			profile.coefficients[static_cast<std::size_t>(y+bases*z)] = std::abs(value) < 1.0e-13 ? 0.0 : value;
		}
	for (int y = 0; y < bases; ++y)
		for (int z = 0; z < bases; ++z) {
			double reconstructed = 0.0;
			for (int first = 0; first < bases; ++first)
				for (int second = 0; second < bases; ++second)
					reconstructed += collocation[static_cast<std::size_t>(y)][static_cast<std::size_t>(first)]
						*profile.At(first, second)*collocation[static_cast<std::size_t>(z)][static_cast<std::size_t>(second)];
			if (std::abs(reconstructed-samples[static_cast<std::size_t>(y)][static_cast<std::size_t>(z)]) > 2.0e-10)
				throw std::runtime_error("C2 square-duct profile collocation residual is too large");
		}
	return profile;
}

inline double EvaluateC2SquareDuctProfile(const C2SquareDuctProfile& profile, double y, double z)
{
	const int bases = OpenUniformCubicBasisCount(profile.elements);
	double value = 0.0;
	for (int first = 0; first < bases; ++first)
		for (int second = 0; second < bases; ++second)
			value += profile.At(first, second)*OpenUniformCubicBasisValue(profile.elements, first, y)
				*OpenUniformCubicBasisValue(profile.elements, second, z);
	return value;
}

inline std::array<double, 3> EvaluateC2SquareDuctProfileWithGradient(
	const C2SquareDuctProfile& profile, double y, double z)
{
	if (!(y >= 0.0 && y <= 1.0) || !(z >= 0.0 && z <= 1.0))
		throw std::runtime_error("C2 square-duct profile gradient requires unit-square coordinates");
	const int n = profile.elements;
	const int element_y = std::min(static_cast<int>(n*y), n-1);
	const int element_z = std::min(static_cast<int>(n*z), n-1);
	const double local_y = n*y-element_y, local_z = n*z-element_z;
	std::array<double, 3> result{};
	for (int local_basis_y = 0; local_basis_y < 4; ++local_basis_y)
		for (int local_basis_z = 0; local_basis_z < 4; ++local_basis_z) {
			const int basis_y = element_y+local_basis_y, basis_z = element_z+local_basis_z;
			const double coefficient = profile.At(basis_y, basis_z);
			const double value_y = OpenUniformCubicExtractedBasisDerivative(
				n, element_y, basis_y, local_y, 0);
			const double value_z = OpenUniformCubicExtractedBasisDerivative(
				n, element_z, basis_z, local_z, 0);
			result[0] += coefficient*value_y*value_z;
			result[1] += coefficient*OpenUniformCubicExtractedBasisDerivative(
				n, element_y, basis_y, local_y, 1)*value_z;
			result[2] += coefficient*value_y*OpenUniformCubicExtractedBasisDerivative(
				n, element_z, basis_z, local_z, 1);
		}
	return result;
}

inline double C2SquareDuctProfileEnergyResistance(const C2SquareDuctProfile& profile)
{
	const double flow = GaussLegendre8UnitSquare([&](double y, double z) {
		return EvaluateC2SquareDuctProfileWithGradient(profile, y, z)[0];
	}, profile.elements);
	const double dissipation = GaussLegendre8UnitSquare([&](double y, double z) {
		const auto value = EvaluateC2SquareDuctProfileWithGradient(profile, y, z);
		return value[1]*value[1]+value[2]*value[2];
	}, profile.elements);
	if (!(flow > 0.0) || !(dissipation > 0.0) || !std::isfinite(flow) || !std::isfinite(dissipation))
		throw std::runtime_error("C2 square-duct energy resistance is invalid");
	return dissipation/(flow*flow);
}

struct C2SquareDuctRitzMetrics {
	double profile_energy_resistance = 0.0;
	double ritz_resistance = 0.0;
	double normalized_weak_residual = 0.0;
};

inline C2SquareDuctRitzMetrics AuditC2SquareDuctRitzProfile(
	const C2SquareDuctProfile& profile)
{
	const int n = profile.elements;
	const int bases = OpenUniformCubicBasisCount(n);
	const int interior_axis = bases-2;
	const int dimension = interior_axis*interior_axis;
	std::vector<std::vector<double>> stiffness(static_cast<std::size_t>(dimension),
		std::vector<double>(static_cast<std::size_t>(dimension), 0.0));
	std::vector<double> load(static_cast<std::size_t>(dimension), 0.0);
	auto interior_index = [interior_axis](int y, int z) {
		return (y-1)+interior_axis*(z-1);
	};
	const auto& nodes = GaussLegendre8Nodes();
	const auto& weights = GaussLegendre8Weights();
	const double scale = 0.5/n;
	for (int element_y = 0; element_y < n; ++element_y)
		for (int element_z = 0; element_z < n; ++element_z)
			for (int quadrature_y = 0; quadrature_y < 8; ++quadrature_y)
				for (int quadrature_z = 0; quadrature_z < 8; ++quadrature_z) {
					const double local_y = 0.5+0.5*nodes[static_cast<std::size_t>(quadrature_y)];
					const double local_z = 0.5+0.5*nodes[static_cast<std::size_t>(quadrature_z)];
					const double weight = scale*scale*weights[static_cast<std::size_t>(quadrature_y)]
						*weights[static_cast<std::size_t>(quadrature_z)];
					struct ActiveBasis { int index; double value; double dy; double dz; };
					std::vector<ActiveBasis> active;
					for (int local_basis_y = 0; local_basis_y < 4; ++local_basis_y)
						for (int local_basis_z = 0; local_basis_z < 4; ++local_basis_z) {
							const int basis_y = element_y+local_basis_y;
							const int basis_z = element_z+local_basis_z;
							if (basis_y == 0 || basis_z == 0 || basis_y+1 == bases || basis_z+1 == bases) continue;
							const double value_y = OpenUniformCubicExtractedBasisDerivative(
								n, element_y, basis_y, local_y, 0);
							const double value_z = OpenUniformCubicExtractedBasisDerivative(
								n, element_z, basis_z, local_z, 0);
							active.push_back({interior_index(basis_y, basis_z), value_y*value_z,
								OpenUniformCubicExtractedBasisDerivative(n, element_y, basis_y, local_y, 1)*value_z,
								value_y*OpenUniformCubicExtractedBasisDerivative(n, element_z, basis_z, local_z, 1)});
						}
					for (const auto& first : active) {
						load[static_cast<std::size_t>(first.index)] += weight*first.value;
						for (const auto& second : active)
							stiffness[static_cast<std::size_t>(first.index)][static_cast<std::size_t>(second.index)]
								+= weight*(first.dy*second.dy+first.dz*second.dz);
					}
				}
	std::vector<double> coefficients(static_cast<std::size_t>(dimension), 0.0);
	for (int z = 1; z+1 < bases; ++z)
		for (int y = 1; y+1 < bases; ++y)
			coefficients[static_cast<std::size_t>(interior_index(y, z))] = profile.At(y, z);
	std::vector<double> stiffness_profile(static_cast<std::size_t>(dimension), 0.0);
	double flow = 0.0, energy = 0.0;
	for (int row = 0; row < dimension; ++row) {
		flow += load[static_cast<std::size_t>(row)]*coefficients[static_cast<std::size_t>(row)];
		for (int column = 0; column < dimension; ++column)
			stiffness_profile[static_cast<std::size_t>(row)]
				+= stiffness[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)]
					*coefficients[static_cast<std::size_t>(column)];
		energy += coefficients[static_cast<std::size_t>(row)]*stiffness_profile[static_cast<std::size_t>(row)];
	}
	const auto inverse = InvertDeterministicMatrix(stiffness);
	double compliance = 0.0;
	for (int row = 0; row < dimension; ++row)
		for (int column = 0; column < dimension; ++column)
			compliance += load[static_cast<std::size_t>(row)]
				*inverse[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)]
				*load[static_cast<std::size_t>(column)];
	if (!(flow > 0.0) || !(energy > 0.0) || !(compliance > 0.0))
		throw std::runtime_error("C2 square-duct Ritz audit is non-positive");
	const double lambda = energy/flow;
	double residual_squared = 0.0, reference_squared = 0.0;
	for (int row = 0; row < dimension; ++row) {
		const double reference = lambda*load[static_cast<std::size_t>(row)];
		const double residual = stiffness_profile[static_cast<std::size_t>(row)]-reference;
		residual_squared += residual*residual;
		reference_squared += reference*reference;
	}
	return {energy/(flow*flow), 1.0/compliance,
		std::sqrt(residual_squared)/std::max(std::sqrt(reference_squared), 1.0e-30)};
}

inline double C2SquareDuctProfileL2ErrorGauss(const C2SquareDuctProfile& profile, int odd_terms = 64)
{
	const double squared = GaussLegendre8UnitSquare([&](double y, double z) {
		const double difference = EvaluateC2SquareDuctProfile(profile, y, z)-SquareDuctSeriesUnnormalized(y, z, odd_terms);
		return difference*difference;
	}, profile.elements);
	return std::sqrt(squared);
}

// Test-only C0 cubic tensor-product unit square duct.  Each transverse cell
// owns a local cubic patch, while adjacent patches share their face basis IDs.
// Keeping the extraction identity makes its h-level geometry exactly affine.
inline void WriteC0SquareDuctDatabase(const std::filesystem::path& path, int transverse_elements,
	std::uint32_t ranks, int axial_elements = 1)
{
	if (transverse_elements < 1 || ranks < 1 || axial_elements < 1)
		throw std::runtime_error("C0 square-duct fixture requires positive element and rank counts");
	const int axial_points = 3*axial_elements+1;
	const int transverse_points = 3*transverse_elements+1;
	const std::uint64_t nodes = static_cast<std::uint64_t>(axial_points)*transverse_points*transverse_points;
	const std::uint64_t elements = static_cast<std::uint64_t>(axial_elements)*transverse_elements*transverse_elements;
	std::vector<std::ostringstream> records(static_cast<std::size_t>(elements));
	std::vector<std::int32_t> owners(static_cast<std::size_t>(elements));
	for (int ex = 0; ex < axial_elements; ++ex)
		for (int ez = 0; ez < transverse_elements; ++ez)
			for (int ey = 0; ey < transverse_elements; ++ey) {
			const std::uint64_t element = static_cast<std::uint64_t>((ex*transverse_elements+ez)*transverse_elements+ey);
			auto& output = records[static_cast<std::size_t>(element)];
			const std::int32_t owner = static_cast<std::int32_t>(element%ranks);
			owners[static_cast<std::size_t>(element)] = owner;
			Write(output, element);
			Write(output, std::int32_t{0});
			Write(output, owner);
			Write(output, std::uint32_t{64});
			// FaceBezierColumns order is z=0, y=0, x=1, y=1, x=0, z=1.
			const std::array<std::int32_t, 6> labels{{ez == 0 ? 0 : -1, ey == 0 ? 0 : -1,
				ex+1 == axial_elements ? 2 : -1, ey+1 == transverse_elements ? 0 : -1,
				ex == 0 ? 1 : -1, ez+1 == transverse_elements ? 0 : -1}};
			output.write(reinterpret_cast<const char*>(labels.data()), static_cast<std::streamsize>(sizeof(labels)));
			for (int k = 0; k < 4; ++k)
				for (int j = 0; j < 4; ++j)
					for (int i = 0; i < 4; ++i) {
						const std::int32_t id = ex*3+i+axial_points*(ey*3+j)
							+axial_points*transverse_points*(ez*3+k);
						Write(output, id);
					}
			for (std::uint8_t row = 0; row < 64; ++row) {
				Write(output, std::uint8_t{1});
				Write(output, row);
				Write(output, 1.0);
			}
			for (int k = 0; k < 4; ++k)
				for (int j = 0; j < 4; ++j)
					for (int i = 0; i < 4; ++i) {
						const std::array<double, 3> point{{(ex+i/3.0)/axial_elements,
							(ey+j/3.0)/transverse_elements, (ez+k/3.0)/transverse_elements}};
						output.write(reinterpret_cast<const char*>(point.data()), static_cast<std::streamsize>(sizeof(point)));
					}
		}

	std::ofstream output(path, std::ios::binary | std::ios::trunc);
	if (!output) throw std::runtime_error("cannot create C0 square-duct fixture database");
	constexpr std::uint64_t header_size = 88;
	const std::uint64_t index_size = (elements+1)*sizeof(std::uint64_t)+elements*sizeof(std::int32_t);
	const std::uint64_t records_offset = header_size+index_size;
	std::vector<std::uint64_t> offsets(static_cast<std::size_t>(elements+1), records_offset);
	for (std::uint64_t e = 0; e < elements; ++e)
		offsets[static_cast<std::size_t>(e+1)] = offsets[static_cast<std::size_t>(e)]
			+static_cast<std::uint64_t>(records[static_cast<std::size_t>(e)].str().size());
	const std::uint64_t rank_index_offset = offsets.back();
	output.write(kMagic.data(), static_cast<std::streamsize>(kMagic.size()));
	Write(output, kVersion); Write(output, ranks); Write(output, elements); Write(output, nodes);
	Write(output, kBezierPointCount); Write(output, std::uint32_t{0}); Write(output, rank_index_offset);
	for (int axis = 0; axis < 3; ++axis) Write(output, 0.0);
	Write(output, 1.0); Write(output, 1.0);
	output.write(reinterpret_cast<const char*>(offsets.data()), static_cast<std::streamsize>(offsets.size()*sizeof(std::uint64_t)));
	output.write(reinterpret_cast<const char*>(owners.data()), static_cast<std::streamsize>(owners.size()*sizeof(std::int32_t)));
	for (const auto& record : records) output << record.str();
	std::vector<std::uint64_t> rank_offsets(static_cast<std::size_t>(ranks)+1, 0);
	for (std::uint32_t rank = 0; rank < ranks; ++rank)
		rank_offsets[static_cast<std::size_t>(rank+1)] = rank_offsets[static_cast<std::size_t>(rank)]+elements;
	output.write(reinterpret_cast<const char*>(rank_offsets.data()), static_cast<std::streamsize>(rank_offsets.size()*sizeof(std::uint64_t)));
	for (std::uint32_t rank = 0; rank < ranks; ++rank)
		for (std::uint64_t element = 0; element < elements; ++element) Write(output, element);
	if (!output) throw std::runtime_error("cannot finalize C0 square-duct fixture database");
}

inline void WriteC0SquareDuctControlMesh(const std::filesystem::path& path, int transverse_elements,
	int axial_elements = 1)
{
	if (transverse_elements < 1 || axial_elements < 1)
		throw std::runtime_error("C0 square-duct control mesh requires positive element count");
	const int points_x = 3*axial_elements+1;
	const int points_yz = 3*transverse_elements+1;
	const int points = points_x*points_yz*points_yz;
	std::ofstream output(path);
	if (!output) throw std::runtime_error("cannot create C0 square-duct control mesh");
	output << "# vtk DataFile Version 3.0\nC0 square duct\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS " << points << " double\n";
	for (int k = 0; k < points_yz; ++k)
		for (int j = 0; j < points_yz; ++j)
			for (int i = 0; i < points_x; ++i) output << i/(3.0*axial_elements) << ' ' << j/(3.0*transverse_elements) << ' ' << k/(3.0*transverse_elements) << '\n';
	const int cells = axial_elements*transverse_elements*transverse_elements;
	output << "CELLS " << cells << ' ' << 9*cells << "\n";
	for (int ex = 0; ex < axial_elements; ++ex)
		for (int ez = 0; ez < transverse_elements; ++ez)
			for (int ey = 0; ey < transverse_elements; ++ey) {
			auto id = [points_x, points_yz](int i, int j, int k) { return i+points_x*j+points_x*points_yz*k; };
			output << "8 " << id(3*ex, 3*ey, 3*ez) << ' ' << id(3*(ex+1), 3*ey, 3*ez) << ' '
				<< id(3*(ex+1), 3*(ey+1), 3*ez) << ' ' << id(3*ex, 3*(ey+1), 3*ez) << ' '
				<< id(3*ex, 3*ey, 3*(ez+1)) << ' ' << id(3*(ex+1), 3*ey, 3*(ez+1)) << ' '
				<< id(3*(ex+1), 3*(ey+1), 3*(ez+1)) << ' ' << id(3*ex, 3*(ey+1), 3*(ez+1)) << '\n';
		}
	output << "CELL_TYPES " << cells << "\n";
	for (int cell = 0; cell < cells; ++cell) output << "12\n";
	output << "POINT_DATA " << points << "\nSCALARS boundary_label int 1\nLOOKUP_TABLE default\n";
	for (int k = 0; k < points_yz; ++k)
		for (int j = 0; j < points_yz; ++j)
			for (int i = 0; i < points_x; ++i) {
				const bool wall = j == 0 || j+1 == points_yz || k == 0 || k+1 == points_yz;
				output << (wall ? 0 : (i == 0 ? 1 : (i+1 == points_x ? 2 : -1))) << '\n';
			}
}

inline void WriteC2SquareDuctDatabase(const std::filesystem::path& path, int transverse_elements,
	std::uint32_t ranks, int axial_elements = 1, double length_m = 1.0)
{
	if (transverse_elements < 1 || axial_elements < 1 || ranks < 1
		|| !(length_m > 0.0) || !std::isfinite(length_m))
		throw std::runtime_error("C2 square-duct fixture requires positive element and rank counts");
	const int bases_x = OpenUniformCubicBasisCount(axial_elements);
	const int bases_yz = OpenUniformCubicBasisCount(transverse_elements);
	const std::uint64_t nodes = static_cast<std::uint64_t>(bases_x)*bases_yz*bases_yz;
	const std::uint64_t elements = static_cast<std::uint64_t>(axial_elements)*transverse_elements*transverse_elements;
	std::vector<std::ostringstream> records(static_cast<std::size_t>(elements));
	std::vector<std::int32_t> owners(static_cast<std::size_t>(elements));
	for (int ex = 0; ex < axial_elements; ++ex)
		for (int ez = 0; ez < transverse_elements; ++ez)
			for (int ey = 0; ey < transverse_elements; ++ey) {
				const std::uint64_t element = static_cast<std::uint64_t>((ex*transverse_elements+ez)*transverse_elements+ey);
				auto& output = records[static_cast<std::size_t>(element)];
				const std::int32_t owner = static_cast<std::int32_t>(element%ranks);
				owners[static_cast<std::size_t>(element)] = owner;
				const auto extraction_x = OpenUniformCubicExtraction(axial_elements, ex);
				const auto extraction_y = OpenUniformCubicExtraction(transverse_elements, ey);
				const auto extraction_z = OpenUniformCubicExtraction(transverse_elements, ez);
				Write(output, element); Write(output, std::int32_t{0}); Write(output, owner); Write(output, std::uint32_t{64});
				const std::array<std::int32_t, 6> labels{{ez == 0 ? 0 : -1, ey == 0 ? 0 : -1,
					ex+1 == axial_elements ? 2 : -1, ey+1 == transverse_elements ? 0 : -1,
					ex == 0 ? 1 : -1, ez+1 == transverse_elements ? 0 : -1}};
				output.write(reinterpret_cast<const char*>(labels.data()), static_cast<std::streamsize>(sizeof(labels)));
				for (int c = 0; c < 4; ++c)
					for (int b = 0; b < 4; ++b)
						for (int a = 0; a < 4; ++a) {
							const std::int32_t id = ex+a+bases_x*(ey+b+bases_yz*(ez+c));
							Write(output, id);
						}
				for (int c = 0; c < 4; ++c)
					for (int b = 0; b < 4; ++b)
						for (int a = 0; a < 4; ++a) {
							std::vector<std::pair<std::uint8_t, double>> entries;
							for (int k = 0; k < 4; ++k)
								for (int j = 0; j < 4; ++j)
									for (int i = 0; i < 4; ++i) {
										const double coefficient = extraction_x[static_cast<std::size_t>(a)][static_cast<std::size_t>(i)]
											*extraction_y[static_cast<std::size_t>(b)][static_cast<std::size_t>(j)]
											*extraction_z[static_cast<std::size_t>(c)][static_cast<std::size_t>(k)];
										if (!std::isfinite(coefficient) || coefficient < -1.0e-11)
											throw std::runtime_error("C2 tensor extraction has an invalid coefficient");
										if (std::abs(coefficient) >= 1.0e-13)
											entries.emplace_back(static_cast<std::uint8_t>(i+4*j+16*k), coefficient);
									}
							if (entries.empty()) throw std::runtime_error("C2 tensor extraction has an empty row");
							Write(output, static_cast<std::uint8_t>(entries.size()));
							for (const auto& entry : entries) { Write(output, entry.first); Write(output, entry.second); }
						}
				for (int k = 0; k < 4; ++k)
					for (int j = 0; j < 4; ++j)
						for (int i = 0; i < 4; ++i) {
							std::array<double, 3> point{{0.0, 0.0, 0.0}};
							for (int c = 0; c < 4; ++c)
								for (int b = 0; b < 4; ++b)
									for (int a = 0; a < 4; ++a) {
										const double coefficient = extraction_x[static_cast<std::size_t>(a)][static_cast<std::size_t>(i)]
											*extraction_y[static_cast<std::size_t>(b)][static_cast<std::size_t>(j)]
											*extraction_z[static_cast<std::size_t>(c)][static_cast<std::size_t>(k)];
										point[0] += coefficient*length_m*OpenUniformCubicGreville(axial_elements, ex+a);
										point[1] += coefficient*OpenUniformCubicGreville(transverse_elements, ey+b);
										point[2] += coefficient*OpenUniformCubicGreville(transverse_elements, ez+c);
									}
							output.write(reinterpret_cast<const char*>(point.data()), static_cast<std::streamsize>(sizeof(point)));
						}
			}
	std::ofstream output(path, std::ios::binary | std::ios::trunc);
	if (!output) throw std::runtime_error("cannot create C2 square-duct fixture database");
	constexpr std::uint64_t header_size = 88;
	const std::uint64_t index_size = (elements+1)*sizeof(std::uint64_t)+elements*sizeof(std::int32_t);
	const std::uint64_t records_offset = header_size+index_size;
	std::vector<std::uint64_t> offsets(static_cast<std::size_t>(elements+1), records_offset);
	for (std::uint64_t e = 0; e < elements; ++e) offsets[static_cast<std::size_t>(e+1)] = offsets[static_cast<std::size_t>(e)]
		+static_cast<std::uint64_t>(records[static_cast<std::size_t>(e)].str().size());
	const std::uint64_t rank_index_offset = offsets.back();
	output.write(kMagic.data(), static_cast<std::streamsize>(kMagic.size()));
	Write(output, kVersion); Write(output, ranks); Write(output, elements); Write(output, nodes);
	Write(output, kBezierPointCount); Write(output, std::uint32_t{0}); Write(output, rank_index_offset);
	for (int axis = 0; axis < 3; ++axis) Write(output, 0.0);
	Write(output, 1.0); Write(output, 1.0);
	output.write(reinterpret_cast<const char*>(offsets.data()), static_cast<std::streamsize>(offsets.size()*sizeof(std::uint64_t)));
	output.write(reinterpret_cast<const char*>(owners.data()), static_cast<std::streamsize>(owners.size()*sizeof(std::int32_t)));
	for (const auto& record : records) output << record.str();
	std::vector<std::uint64_t> rank_offsets(static_cast<std::size_t>(ranks)+1, 0);
	for (std::uint32_t rank = 0; rank < ranks; ++rank) rank_offsets[static_cast<std::size_t>(rank+1)] = rank_offsets[static_cast<std::size_t>(rank)]+elements;
	output.write(reinterpret_cast<const char*>(rank_offsets.data()), static_cast<std::streamsize>(rank_offsets.size()*sizeof(std::uint64_t)));
	for (std::uint32_t rank = 0; rank < ranks; ++rank)
		for (std::uint64_t element = 0; element < elements; ++element) Write(output, element);
	if (!output) throw std::runtime_error("cannot finalize C2 square-duct fixture database");
}

inline void WriteC2SquareDuctControlMesh(const std::filesystem::path& path, int transverse_elements,
	int axial_elements = 1, double length_m = 1.0)
{
	if (transverse_elements < 1 || axial_elements < 1
		|| !(length_m > 0.0) || !std::isfinite(length_m))
		throw std::runtime_error("C2 square-duct control mesh requires positive elements and length");
	const int bases_x = OpenUniformCubicBasisCount(axial_elements);
	const int bases_yz = OpenUniformCubicBasisCount(transverse_elements);
	const int points = bases_x*bases_yz*bases_yz;
	std::ofstream output(path);
	if (!output) throw std::runtime_error("cannot create C2 square-duct control mesh");
	output << "# vtk DataFile Version 3.0\nC2 square duct\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS " << points << " double\n";
	for (int k = 0; k < bases_yz; ++k)
		for (int j = 0; j < bases_yz; ++j)
			for (int i = 0; i < bases_x; ++i)
				output << length_m*OpenUniformCubicGreville(axial_elements, i) << ' '
					<< OpenUniformCubicGreville(transverse_elements, j) << ' '
					<< OpenUniformCubicGreville(transverse_elements, k) << '\n';
	const int cells = axial_elements*transverse_elements*transverse_elements;
	output << "CELLS " << cells << ' ' << 9*cells << "\n";
	auto id = [bases_x, bases_yz](int i, int j, int k) { return i+bases_x*(j+bases_yz*k); };
	for (int ex = 0; ex < axial_elements; ++ex)
		for (int ez = 0; ez < transverse_elements; ++ez)
			for (int ey = 0; ey < transverse_elements; ++ey) {
				const int x0 = OpenUniformCubicCellCornerBasis(ex, axial_elements), x1 = OpenUniformCubicCellCornerBasis(ex+1, axial_elements);
				const int y0 = OpenUniformCubicCellCornerBasis(ey, transverse_elements), y1 = OpenUniformCubicCellCornerBasis(ey+1, transverse_elements);
				const int z0 = OpenUniformCubicCellCornerBasis(ez, transverse_elements), z1 = OpenUniformCubicCellCornerBasis(ez+1, transverse_elements);
				output << "8 " << id(x0, y0, z0) << ' ' << id(x1, y0, z0) << ' ' << id(x1, y1, z0) << ' ' << id(x0, y1, z0) << ' '
					<< id(x0, y0, z1) << ' ' << id(x1, y0, z1) << ' ' << id(x1, y1, z1) << ' ' << id(x0, y1, z1) << '\n';
			}
	output << "CELL_TYPES " << cells << "\n";
	for (int cell = 0; cell < cells; ++cell) output << "12\n";
	output << "POINT_DATA " << points << "\nSCALARS boundary_label int 1\nLOOKUP_TABLE default\n";
	for (int k = 0; k < bases_yz; ++k)
		for (int j = 0; j < bases_yz; ++j)
			for (int i = 0; i < bases_x; ++i) {
				const bool wall = j == 0 || j+1 == bases_yz || k == 0 || k+1 == bases_yz;
				output << (wall ? 0 : (i == 0 ? 1 : (i+1 == bases_x ? 2 : -1))) << '\n';
			}
}

inline void WriteC2SquareDuctReferenceVelocity(const std::filesystem::path& path, int transverse_elements,
	double reference_profile_scale_m_s, int axial_elements = 1, int odd_terms = 64)
{
	if (!std::isfinite(reference_profile_scale_m_s) || transverse_elements < 1 || axial_elements < 1)
		throw std::runtime_error("C2 square-duct reference velocity requires finite scale and positive elements");
	const auto profile = BuildC2SquareDuctProfile(transverse_elements, odd_terms);
	const int bases_x = OpenUniformCubicBasisCount(axial_elements);
	const int bases_yz = OpenUniformCubicBasisCount(transverse_elements);
	std::ofstream output(path);
	if (!output) throw std::runtime_error("cannot create C2 square-duct reference velocity");
	for (int k = 0; k < bases_yz; ++k)
		for (int j = 0; j < bases_yz; ++j)
			for (int i = 0; i < bases_x; ++i)
				output << std::setprecision(17) << (i == 0 ? reference_profile_scale_m_s*profile.At(j, k) : 0.0) << " 0 0\n";
}

inline void WriteC2SquareDuctThreeDCase(const std::filesystem::path& directory, int transverse_elements,
	double dt_s, int steps, double density_kg_m3, double dynamic_viscosity_pa_s,
	double reference_profile_scale_m_s, int axial_elements = 1, double length_m = 1.0)
{
	if (!(dt_s > 0.0) || steps < 1 || !(density_kg_m3 > 0.0) || !(dynamic_viscosity_pa_s > 0.0)
		|| !(length_m > 0.0) || !std::isfinite(length_m))
		throw std::runtime_error("C2 square-duct 3D case requires positive time and material controls");
	std::filesystem::create_directories(directory);
	WriteC2SquareDuctControlMesh(directory/"controlmesh.vtk", transverse_elements, axial_elements, length_m);
	WriteC2SquareDuctReferenceVelocity(directory/"initial_velocityfield.txt", transverse_elements,
		reference_profile_scale_m_s, axial_elements);
	std::ofstream output(directory/"simulation_config.json");
	if (!output) throw std::runtime_error("cannot create C2 square-duct 3D configuration");
	output << std::setprecision(17)
		<< "{\n  \"schema_version\": 3, \"dimension\": \"3d\",\n"
		<< "  \"simulation_scope\": {\"mode\": \"flow_only\"},\n"
		<< "  \"coupling\": {\"scheme\": \"explicit_staggered\", \"flow_epsilon_m3_s\": 1e-14, \"three_d_ports\": {\"inlet_label\": 1, \"outlet_labels\": [2]}},\n"
		<< "  \"fields\": [{\"name\": \"velocity\", \"kind\": \"vector3\"}, {\"name\": \"pressure\", \"kind\": \"pressure\"}],\n"
		<< "  \"time\": {\"dt\": " << dt_s << ", \"steps\": " << steps << "},\n"
		<< "  \"equation_systems\": [{\"name\": \"flow\", \"kind\": \"navier_stokes\", \"unknowns\": [\"velocity\", \"pressure\"], \"viscosity\": " << dynamic_viscosity_pa_s << ", \"density\": " << density_kg_m3 << ", \"time_integration\": \"backward_euler\"}],\n"
		<< "  \"boundaries\": [{\"label\": 0, \"name\": \"no_slip_wall\", \"conditions\": [{\"field\": \"velocity\", \"type\": \"dirichlet\", \"value\": [0,0,0]}]},\n"
		<< "    {\"label\": 1, \"name\": \"inlet\", \"conditions\": [{\"field\": \"velocity\", \"type\": \"dirichlet\", \"profile\": \"initial_velocityfield.txt\", \"scale\": 1}]},\n"
		<< "    {\"label\": 2, \"name\": \"outlet\", \"conditions\": [{\"field\": \"pressure\", \"type\": \"pressure_traction\", \"value\": 0}]}]\n}\n";
}

inline void WriteSquareDuctReferenceVelocity(const std::filesystem::path& path, int transverse_elements,
	double reference_profile_scale_m_s, int odd_terms = 64,
	SquareDuctProfileMode mode = SquareDuctProfileMode::CorrectedBernstein, int axial_elements = 1)
{
	if (!std::isfinite(reference_profile_scale_m_s) || axial_elements < 1)
		throw std::runtime_error("square-duct reference velocity must be finite with positive axial elements");
	const int points_x = 3*axial_elements+1;
	const int points_yz = 3*transverse_elements+1;
	std::ofstream output(path);
	if (!output) throw std::runtime_error("cannot create square-duct reference velocity");
	for (int k = 0; k < points_yz; ++k)
		for (int j = 0; j < points_yz; ++j)
			for (int i = 0; i < points_x; ++i) {
				const double y = j/(3.0*transverse_elements);
				const double z = k/(3.0*transverse_elements);
				const double coefficient = mode == SquareDuctProfileMode::CorrectedBernstein
					? SquareDuctBernsteinCoefficient(transverse_elements, j, k, odd_terms)
					: SquareDuctSeriesUnnormalized(y, z, odd_terms);
				const double velocity = i == 0 ? reference_profile_scale_m_s*coefficient : 0.0;
				output << std::setprecision(17) << velocity << " 0 0\n";
			}
}

inline void WriteSquareDuctThreeDCase(const std::filesystem::path& directory, int transverse_elements,
	double dt_s, int steps, double density_kg_m3, double dynamic_viscosity_pa_s,
	double reference_profile_scale_m_s, int axial_elements = 1)
{
	if (!(dt_s > 0.0) || steps < 1 || !(density_kg_m3 > 0.0) || !(dynamic_viscosity_pa_s > 0.0))
		throw std::runtime_error("square-duct 3D case requires positive time and material controls");
	std::filesystem::create_directories(directory);
	WriteC0SquareDuctControlMesh(directory/"controlmesh.vtk", transverse_elements, axial_elements);
	WriteSquareDuctReferenceVelocity(directory/"initial_velocityfield.txt", transverse_elements,
		reference_profile_scale_m_s, 64, SquareDuctProfileMode::CorrectedBernstein, axial_elements);
	std::ofstream output(directory/"simulation_config.json");
	if (!output) throw std::runtime_error("cannot create square-duct 3D configuration");
	output << std::setprecision(17)
		<< "{\n  \"schema_version\": 3, \"dimension\": \"3d\",\n"
		<< "  \"simulation_scope\": {\"mode\": \"flow_only\"},\n"
		<< "  \"coupling\": {\"scheme\": \"explicit_staggered\", \"flow_epsilon_m3_s\": 1e-14, \"three_d_ports\": {\"inlet_label\": 1, \"outlet_labels\": [2]}},\n"
		<< "  \"fields\": [{\"name\": \"velocity\", \"kind\": \"vector3\"}, {\"name\": \"pressure\", \"kind\": \"pressure\"}],\n"
		<< "  \"time\": {\"dt\": " << dt_s << ", \"steps\": " << steps << "},\n"
		<< "  \"equation_systems\": [{\"name\": \"flow\", \"kind\": \"navier_stokes\", \"unknowns\": [\"velocity\", \"pressure\"], \"viscosity\": " << dynamic_viscosity_pa_s << ", \"density\": " << density_kg_m3 << ", \"time_integration\": \"backward_euler\"}],\n"
		<< "  \"boundaries\": [{\"label\": 0, \"name\": \"no_slip_wall\", \"conditions\": [{\"field\": \"velocity\", \"type\": \"dirichlet\", \"value\": [0,0,0]}]},\n"
		<< "    {\"label\": 1, \"name\": \"inlet\", \"conditions\": [{\"field\": \"velocity\", \"type\": \"dirichlet\", \"profile\": \"initial_velocityfield.txt\", \"scale\": 1}]},\n"
		<< "    {\"label\": 2, \"name\": \"outlet\", \"conditions\": [{\"field\": \"pressure\", \"type\": \"pressure_traction\", \"value\": 0}]}]\n}\n";
}

inline void WriteRigidOneDStraightCase(const std::filesystem::path& directory, double length_m,
	double radius_m, double dt_s, int steps, double density_kg_m3, double dynamic_viscosity_pa_s,
	double inlet_flow_m3_s, const std::string& temporal_kind = "constant", double period_s = 1.0,
	int outlet_node_id = 2, const std::string& scheme = "steady_poiseuille")
{
	if (!(length_m > 0.0) || !(radius_m > 0.0) || !(dt_s > 0.0) || steps < 1 || !(period_s > 0.0)
		|| outlet_node_id < 2 || (temporal_kind != "constant" && temporal_kind != "sinusoid")
		|| (scheme != "steady_poiseuille" && scheme != "rigid_inertance"))
		throw std::runtime_error("rigid 1D fixture requires positive geometry and time controls");
	std::filesystem::create_directories(directory);
	std::ofstream network(directory/"tree.swc");
	if (!network) throw std::runtime_error("cannot create rigid 1D fixture network");
	network << std::setprecision(17) << "1 2 0 0 0 " << radius_m << " -1\n2 2 " << length_m << " 0 0 " << radius_m << " 1\n";
	std::ofstream output(directory/"simulation_config.json");
	if (!output) throw std::runtime_error("cannot create rigid 1D fixture configuration");
	output << std::setprecision(17)
		<< "{\"schema_version\":3,\"dimension\":\"1d\",\"simulation_scope\":{\"mode\":\"flow_only\"},\n"
		<< "\"geometry\":{\"kind\":\"swc_network\",\"file\":\"tree.swc\",\"length_scale_to_m\":1},\n"
		<< "\"fields\":[{\"name\":\"area\",\"kind\":\"scalar\"},{\"name\":\"flow_rate\",\"kind\":\"scalar\"},{\"name\":\"pressure\",\"kind\":\"pressure\"}],\n"
		<< "\"time\":{\"dt\":" << dt_s << ",\"steps\":" << steps << ",\"output_every\":1},\n"
		<< "\"temporal_functions\":[{\"name\":\"inlet_flow\",\"kind\":\"" << temporal_kind << "\",\"units\":\"m3/s\"";
	if (temporal_kind == "constant") output << ",\"value\":" << inlet_flow_m3_s;
	else output << ",\"mean\":" << inlet_flow_m3_s << ",\"amplitude\":" << 0.2*inlet_flow_m3_s << ",\"period\":" << period_s << ",\"phase\":0";
	output << "}],\n\"equation_systems\":[{\"name\":\"flow\",\"kind\":\"network_flow_1d\",\"unknowns\":[\"area\",\"flow_rate\",\"pressure\"],\"model\":\"rigid\",\"scheme\":\"" << scheme << "\",\"dynamic_viscosity\":" << dynamic_viscosity_pa_s << ",\"density\":" << density_kg_m3 << ",\"discretization\":{\"cells_per_segment\":1}}],\n"
		<< "\"boundaries\":[{\"name\":\"inlet\",\"role\":\"inlet\",\"node_ids\":[1],\"conditions\":[{\"field\":\"flow_rate\",\"type\":\"dirichlet\",\"quantity\":\"flow_rate\",\"waveform\":\"inlet_flow\"}]},{\"name\":\"outlet\",\"role\":\"outlet\",\"node_ids\":[" << outlet_node_id << "],\"conditions\":[{\"field\":\"pressure\",\"type\":\"pressure\",\"value\":0}]}]}\n";
}

} // namespace iga::test

#endif
