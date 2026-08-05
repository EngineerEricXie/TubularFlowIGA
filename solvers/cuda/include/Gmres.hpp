#ifndef IGA_CUDA_GMRES_HPP
#define IGA_CUDA_GMRES_HPP

#include "IgaCudaKernels.cuh"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace iga::cuda {

struct GmresResult {
	int iterations = 0;
	double residual = std::numeric_limits<double>::infinity();
	unsigned int singular_diagonal_blocks = 0;
	double device_used_gib = 0.0;
	bool converged = false;
};

template <int Fields>
class GmresWorkspace {
public:
	GmresWorkspace(int nodes, int restart_size)
		: size(nodes*Fields), restart(restart_size),
		  inverse(static_cast<std::size_t>(nodes)*Fields*Fields), singular(1),
		  residual(size), raw(size), work(size), preconditioner_residual(size), correction(size),
		  vectors(static_cast<std::size_t>(restart_size+1)*size),
		  h(static_cast<std::size_t>(restart_size+1)*restart_size),
		  cosine(restart_size), sine(restart_size), g(restart_size+1), y(restart_size)
	{
		if (nodes <= 0 || restart_size <= 0)
			throw std::invalid_argument("GMRES workspace dimensions must be positive");
	}

	int size;
	int restart;
	BlasHandle blas;
	DeviceBuffer<double> inverse;
	DeviceBuffer<unsigned int> singular;
	DeviceBuffer<double> residual, raw, work, preconditioner_residual, correction, vectors;
	std::vector<double> h, cosine, sine, g, y;
};

template <int Fields>
inline void ApplyPolynomialBlockJacobi(BlasHandle& blas, const BlockMatrix<Fields>& matrix,
	const double* inverse, const double* input, double* output, double* residual,
	double* correction, int sweeps)
{
	const auto pattern = matrix.pattern();
	const int size = pattern.nodes*Fields;
	const int node_blocks = (pattern.nodes+255)/256;
	ApplyBlockInverseKernel<Fields><<<node_blocks,256>>>(pattern.nodes, inverse, input, output);
	CheckKernel("ApplyBlockInverseKernel");
	const double minus_one = -1.0, one = 1.0;
	for (int sweep = 0; sweep < sweeps; ++sweep) {
		BlockSpmvKernel<Fields><<<node_blocks,256>>>(pattern, matrix.values(), output, residual);
		CheckKernel("BlockSpmvKernel preconditioner");
		Check(cublasDscal(blas, size, &minus_one, residual, 1), "cublasDscal preconditioner");
		Check(cublasDaxpy(blas, size, &one, input, 1, residual, 1), "cublasDaxpy preconditioner residual");
		ApplyBlockInverseKernel<Fields><<<node_blocks,256>>>(pattern.nodes, inverse, residual, correction);
		CheckKernel("ApplyBlockInverseKernel correction");
		Check(cublasDaxpy(blas, size, &one, correction, 1, output, 1), "cublasDaxpy preconditioner correction");
	}
}

template <int Fields>
inline GmresResult SolveGmres(const BlockMatrix<Fields>& matrix, GmresWorkspace<Fields>& workspace,
	const double* rhs, double* solution, int maximum_iterations, double relative_tolerance,
	bool verbose = false, int preconditioner_sweeps = 0, bool rebuild_preconditioner = true)
{
	const auto pattern = matrix.pattern();
	const int size = pattern.nodes*Fields;
	if (workspace.size != size) throw std::invalid_argument("GMRES workspace size does not match matrix");
	const int restart = workspace.restart;
	const int node_blocks = (pattern.nodes+255)/256;
	auto& blas = workspace.blas;
	auto& inverse = workspace.inverse;
	auto& singular = workspace.singular;
	auto& residual = workspace.residual;
	auto& raw = workspace.raw;
	auto& work = workspace.work;
	auto& preconditioner_residual = workspace.preconditioner_residual;
	auto& correction = workspace.correction;
	auto& vectors = workspace.vectors;
	auto& h = workspace.h;
	auto& cosine = workspace.cosine;
	auto& sine = workspace.sine;
	auto& g = workspace.g;
	auto& y = workspace.y;
	if (rebuild_preconditioner) {
		singular.Clear();
		BuildBlockInverseKernel<Fields><<<node_blocks,256>>>(
			pattern, matrix.values(), inverse.data(), singular.data());
		CheckKernel("BuildBlockInverseKernel");
	}
	BlockSpmvKernel<Fields><<<node_blocks,256>>>(pattern, matrix.values(), solution, raw.data());
	CheckKernel("BlockSpmvKernel initial");
	Check(cublasDcopy(blas, size, rhs, 1, residual.data(), 1), "cublasDcopy");
	const double minus_one = -1.0;
	Check(cublasDaxpy(blas, size, &minus_one, raw.data(), 1, residual.data(), 1), "cublasDaxpy residual");
	ApplyPolynomialBlockJacobi(blas, matrix, inverse.data(), residual.data(), work.data(),
		preconditioner_residual.data(), correction.data(), preconditioner_sweeps);
	double initial_norm = 0.0;
	Check(cublasDnrm2(blas, size, work.data(), 1, &initial_norm), "cublasDnrm2 initial");
	GmresResult result;
	std::size_t free_device = 0, total_device = 0;
	Check(cudaMemGetInfo(&free_device, &total_device), "cudaMemGetInfo GMRES");
	result.device_used_gib = Gibibytes(total_device-free_device);
	singular.CopyToHost(&result.singular_diagonal_blocks, 1);
	if (initial_norm == 0.0) {
		result.residual = 0.0;
		result.converged = true;
		return result;
	}
	const double tolerance = std::max(1e-30, relative_tolerance*initial_norm);
	double beta = initial_norm;
	while (result.iterations < maximum_iterations) {
		Check(cublasDcopy(blas, size, work.data(), 1, vectors.data(), 1), "cublasDcopy V0");
		const double inverse_beta = 1.0/beta;
		Check(cublasDscal(blas, size, &inverse_beta, vectors.data(), 1), "cublasDscal V0");
		std::fill(h.begin(), h.end(), 0.0);
		std::fill(g.begin(), g.end(), 0.0);
		g[0] = beta;
		int used = 0;
		bool estimated_converged = false;
		for (int j = 0; j < restart && result.iterations < maximum_iterations; ++j) {
			const double* vj = vectors.data()+static_cast<std::size_t>(j)*size;
			BlockSpmvKernel<Fields><<<node_blocks,256>>>(pattern, matrix.values(), vj, raw.data());
			CheckKernel("BlockSpmvKernel Arnoldi");
			ApplyPolynomialBlockJacobi(blas, matrix, inverse.data(), raw.data(), work.data(),
				preconditioner_residual.data(), correction.data(), preconditioner_sweeps);
			for (int i = 0; i <= j; ++i) {
				const double* vi = vectors.data()+static_cast<std::size_t>(i)*size;
				double value = 0.0;
				Check(cublasDdot(blas, size, vi, 1, work.data(), 1, &value), "cublasDdot");
				h[static_cast<std::size_t>(i)*restart+j] = value;
				const double negative = -value;
				Check(cublasDaxpy(blas, size, &negative, vi, 1, work.data(), 1), "cublasDaxpy Arnoldi");
			}
			double next_norm = 0.0;
			Check(cublasDnrm2(blas, size, work.data(), 1, &next_norm), "cublasDnrm2 Arnoldi");
			h[static_cast<std::size_t>(j+1)*restart+j] = next_norm;
			if (next_norm > 1e-30) {
				double* next = vectors.data()+static_cast<std::size_t>(j+1)*size;
				Check(cublasDcopy(blas, size, work.data(), 1, next, 1), "cublasDcopy Vnext");
				const double inverse_norm = 1.0/next_norm;
				Check(cublasDscal(blas, size, &inverse_norm, next, 1), "cublasDscal Vnext");
			}
			for (int i = 0; i < j; ++i) {
				const double first = h[static_cast<std::size_t>(i)*restart+j];
				const double second = h[static_cast<std::size_t>(i+1)*restart+j];
				h[static_cast<std::size_t>(i)*restart+j] = cosine[i]*first+sine[i]*second;
				h[static_cast<std::size_t>(i+1)*restart+j] = -sine[i]*first+cosine[i]*second;
			}
			const double first = h[static_cast<std::size_t>(j)*restart+j];
			const double second = h[static_cast<std::size_t>(j+1)*restart+j];
			const double magnitude = std::hypot(first, second);
			cosine[j] = magnitude > 0.0 ? first/magnitude : 1.0;
			sine[j] = magnitude > 0.0 ? second/magnitude : 0.0;
			h[static_cast<std::size_t>(j)*restart+j] = magnitude;
			h[static_cast<std::size_t>(j+1)*restart+j] = 0.0;
			const double old_g = g[j];
			g[j] = cosine[j]*old_g;
			g[j+1] = -sine[j]*old_g;
			++result.iterations;
			used = j+1;
			result.residual = std::abs(g[j+1]);
			if (verbose && (result.iterations == 1 || result.iterations%25 == 0))
				std::cout << "  gmres iteration=" << result.iterations
					<< " preconditioned_residual=" << result.residual << '\n';
			if (result.residual <= tolerance || next_norm <= 1e-30) {
				estimated_converged = true;
				break;
			}
		}
		for (int row = used-1; row >= 0; --row) {
			double value = g[row];
			for (int column = row+1; column < used; ++column)
				value -= h[static_cast<std::size_t>(row)*restart+column]*y[column];
			const double diagonal = h[static_cast<std::size_t>(row)*restart+row];
			if (std::abs(diagonal) < 1e-30) throw std::runtime_error("GMRES Hessenberg matrix is singular");
			y[row] = value/diagonal;
		}
		for (int i = 0; i < used; ++i)
			Check(cublasDaxpy(blas, size, &y[i],
				vectors.data()+static_cast<std::size_t>(i)*size, 1, solution, 1), "cublasDaxpy solution");

		BlockSpmvKernel<Fields><<<node_blocks,256>>>(pattern, matrix.values(), solution, raw.data());
		CheckKernel("BlockSpmvKernel restart");
		Check(cublasDcopy(blas, size, rhs, 1, residual.data(), 1), "cublasDcopy restart residual");
		Check(cublasDaxpy(blas, size, &minus_one, raw.data(), 1, residual.data(), 1), "cublasDaxpy restart residual");
		ApplyPolynomialBlockJacobi(blas, matrix, inverse.data(), residual.data(), work.data(),
			preconditioner_residual.data(), correction.data(), preconditioner_sweeps);
		Check(cublasDnrm2(blas, size, work.data(), 1, &beta), "cublasDnrm2 restart");
		result.residual = beta;
		if (beta <= tolerance) {
			result.converged = true;
			break;
		}
		if (estimated_converged && beta > tolerance && verbose)
			std::cout << "  gmres residual replacement=" << beta << '\n';
	}
	return result;
}

} // namespace iga::cuda

#endif
