#ifndef IGA_CUDA_SCALAR_KRYLOV_HPP
#define IGA_CUDA_SCALAR_KRYLOV_HPP

#include "SparseKernels.cuh"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

// Jacobi-preconditioned CG and BiCGStab for scalar (Fields=1) CSR systems.
// Every scalar stays on the device, so an iteration launches kernels without
// waiting for the GPU; the host reads the residual norm only every
// kCheckInterval iterations. Both solvers stop on the Jacobi-preconditioned
// residual norm ||D^-1 r|| <= max(rtol*||D^-1 b||, atol), PETSc's default
// test for left preconditioning from a zero initial guess (the solution
// array may hold any initial guess).

namespace iga::cuda {

constexpr int kCheckInterval = 8;

struct KrylovResult {
	int iterations = 0;
	double residual = std::numeric_limits<double>::infinity();
	double tolerance = 0.0;
	bool converged = false;
};

// Scalars of one solve, kept in device memory.
struct KrylovScalars {
	double rho, rho_previous, alpha, omega, beta, curvature, ts, tt, norm, load_norm;
	int breakdown;
};

class ScalarKrylovWorkspace {
public:
	explicit ScalarKrylovWorkspace(int size)
		: size(size), inverse(size), r(size), z(size), p(size), v(size), s(size), t(size), shadow(size),
		  scalars(1), singular(1)
	{
		if (size <= 0) throw std::invalid_argument("Krylov workspace size must be positive");
		Check(cublasSetPointerMode(blas, CUBLAS_POINTER_MODE_DEVICE), "cublasSetPointerMode device");
	}

	int size;
	BlasHandle blas;
	DeviceBuffer<double> inverse, r, z, p, v, s, t, shadow;
	DeviceBuffer<KrylovScalars> scalars;
	DeviceBuffer<unsigned int> singular;
};

__global__ void ScaledSpmvKernel(DevicePatternView pattern, const double* values, const double* inverse,
	const double* x, double* y)
{
	const int row = blockIdx.x*blockDim.x+threadIdx.x;
	if (row >= pattern.nodes) return;
	double sum = 0.0;
	for (int entry = pattern.row_offsets[row]; entry < pattern.row_offsets[row+1]; ++entry)
		sum += values[entry]*x[pattern.columns[entry]];
	y[row] = inverse[row]*sum;
}

// r = D^-1 (b - A x); with left preconditioning the solver works on D^-1 A.
__global__ void ScaledResidualKernel(DevicePatternView pattern, const double* values, const double* inverse,
	const double* b, const double* x, double* r)
{
	const int row = blockIdx.x*blockDim.x+threadIdx.x;
	if (row >= pattern.nodes) return;
	double sum = 0.0;
	for (int entry = pattern.row_offsets[row]; entry < pattern.row_offsets[row+1]; ++entry)
		sum += values[entry]*x[pattern.columns[entry]];
	r[row] = inverse[row]*(b[row]-sum);
}

__global__ void ScaleKernel(int size, const double* inverse, const double* x, double* y)
{
	const int row = blockIdx.x*blockDim.x+threadIdx.x;
	if (row < size) y[row] = inverse[row]*x[row];
}

// CG on the SPD matrix A with z = D^-1 r. The curvature is p.Ap.
__global__ void ConjugateGradientStepKernel(int size, KrylovScalars* scalars, const double* inverse,
	const double* p, const double* q, double* x, double* r, double* z)
{
	const int row = blockIdx.x*blockDim.x+threadIdx.x;
	if (row >= size) return;
	const double curvature = scalars->curvature;
	if (!(curvature > 0.0) || !isfinite(curvature)) {
		if (row == 0) scalars->breakdown = 1;
		return;
	}
	const double alpha = scalars->rho/curvature;
	x[row] += alpha*p[row];
	r[row] -= alpha*q[row];
	z[row] = inverse[row]*r[row];
}

__global__ void ConjugateGradientDirectionScalarKernel(KrylovScalars* scalars)
{
	scalars->beta = scalars->rho_previous == 0.0 ? 0.0 : scalars->rho/scalars->rho_previous;
}

__global__ void ConjugateGradientDirectionKernel(int size, const KrylovScalars* scalars, const double* z, double* p)
{
	const int row = blockIdx.x*blockDim.x+threadIdx.x;
	if (row < size) p[row] = z[row]+scalars->beta*p[row];
}

__global__ void ShiftScalarKernel(KrylovScalars* scalars)
{
	scalars->rho_previous = scalars->rho;
}

inline void SetupJacobi(const BlockMatrix<1>& matrix, ScalarKrylovWorkspace& workspace)
{
	const auto pattern = matrix.pattern();
	workspace.singular.Clear();
	BuildScalarJacobiKernel<<<(pattern.nodes+255)/256,256>>>(pattern, matrix.values(), workspace.inverse.data(),
		workspace.singular.data());
	CheckKernel("BuildScalarJacobiKernel");
	unsigned int singular = 0;
	workspace.singular.CopyToHost(&singular, 1);
	if (singular) throw std::runtime_error("Jacobi preconditioner met zero or non-finite diagonal entries");
}

inline KrylovScalars ReadScalars(const ScalarKrylovWorkspace& workspace)
{
	KrylovScalars scalars{};
	workspace.scalars.CopyToHost(&scalars, 1);
	return scalars;
}

// Jacobi-preconditioned CG for a symmetric positive definite matrix.
inline KrylovResult SolveConjugateGradient(const BlockMatrix<1>& matrix, ScalarKrylovWorkspace& workspace,
	const double* rhs, double* solution, int maximum_iterations, double relative_tolerance, double absolute_tolerance)
{
	const auto pattern = matrix.pattern();
	const int size = pattern.nodes, blocks = (size+255)/256;
	if (workspace.size != size) throw std::invalid_argument("Krylov workspace size does not match matrix");
	SetupJacobi(matrix, workspace);
	auto& blas = workspace.blas;
	auto* scalars = workspace.scalars.data();
	workspace.scalars.Clear();
	double *r = workspace.r.data(), *z = workspace.z.data(), *p = workspace.p.data(), *q = workspace.v.data();
	// Reference norm ||D^-1 b|| and the initial residual r = b - A x, z = D^-1 r.
	ScaleKernel<<<blocks,256>>>(size, workspace.inverse.data(), rhs, z);
	Check(cublasDnrm2(blas, size, z, 1, &scalars->load_norm), "cublasDnrm2 CG load");
	BlockSpmvKernel<1><<<blocks,256>>>(pattern, matrix.values(), solution, q);
	Check(cudaMemcpy(r, rhs, static_cast<std::size_t>(size)*sizeof(double), cudaMemcpyDeviceToDevice), "copy CG rhs");
	const double minus_one = -1.0;
	Check(cublasSetPointerMode(blas, CUBLAS_POINTER_MODE_HOST), "cublasSetPointerMode host");
	Check(cublasDaxpy(blas, size, &minus_one, q, 1, r, 1), "cublasDaxpy CG residual");
	Check(cublasSetPointerMode(blas, CUBLAS_POINTER_MODE_DEVICE), "cublasSetPointerMode device");
	ScaleKernel<<<blocks,256>>>(size, workspace.inverse.data(), r, z);
	Check(cublasDnrm2(blas, size, z, 1, &scalars->norm), "cublasDnrm2 CG initial");
	Check(cublasDdot(blas, size, r, 1, z, 1, &scalars->rho), "cublasDdot CG rho");
	Check(cudaMemcpy(p, z, static_cast<std::size_t>(size)*sizeof(double), cudaMemcpyDeviceToDevice), "copy CG direction");
	CheckKernel("CG setup");
	auto host = ReadScalars(workspace);
	KrylovResult result;
	result.tolerance = std::max(relative_tolerance*host.load_norm, absolute_tolerance);
	result.residual = host.norm;
	if (host.norm <= result.tolerance) {
		result.converged = true;
		return result;
	}
	while (result.iterations < maximum_iterations) {
		ShiftScalarKernel<<<1,1>>>(scalars);
		BlockSpmvKernel<1><<<blocks,256>>>(pattern, matrix.values(), p, q);
		Check(cublasDdot(blas, size, p, 1, q, 1, &scalars->curvature), "cublasDdot CG curvature");
		ConjugateGradientStepKernel<<<blocks,256>>>(size, scalars, workspace.inverse.data(), p, q, solution, r, z);
		Check(cublasDdot(blas, size, r, 1, z, 1, &scalars->rho), "cublasDdot CG rho");
		Check(cublasDnrm2(blas, size, z, 1, &scalars->norm), "cublasDnrm2 CG");
		ConjugateGradientDirectionScalarKernel<<<1,1>>>(scalars);
		ConjugateGradientDirectionKernel<<<blocks,256>>>(size, scalars, z, p);
		++result.iterations;
		if (result.iterations % kCheckInterval && result.iterations < maximum_iterations) continue;
		CheckKernel("CG iteration");
		host = ReadScalars(workspace);
		result.residual = host.norm;
		if (host.breakdown) throw std::runtime_error("CG met a non-positive curvature; the matrix is not SPD");
		if (!std::isfinite(host.norm)) break;
		if (host.norm <= result.tolerance) {
			result.converged = true;
			break;
		}
	}
	return result;
}

__global__ void BiCgStabDirectionKernel(int size, KrylovScalars* scalars, const double* r, const double* v, double* p)
{
	const int row = blockIdx.x*blockDim.x+threadIdx.x;
	if (row >= size) return;
	const double beta = (scalars->rho/scalars->rho_previous)*(scalars->alpha/scalars->omega);
	p[row] = r[row]+beta*(p[row]-scalars->omega*v[row]);
}

__global__ void BiCgStabAlphaKernel(KrylovScalars* scalars)
{
	// curvature holds shadow.v here.
	if (scalars->curvature == 0.0 || !isfinite(scalars->curvature)) scalars->breakdown = 1;
	else scalars->alpha = scalars->rho/scalars->curvature;
}

__global__ void BiCgStabHalfStepKernel(int size, const KrylovScalars* scalars, const double* r, const double* v, double* s)
{
	const int row = blockIdx.x*blockDim.x+threadIdx.x;
	if (row < size) s[row] = r[row]-scalars->alpha*v[row];
}

__global__ void BiCgStabOmegaKernel(KrylovScalars* scalars)
{
	// t = D^-1 A s vanishes only when s does (A is nonsingular): the half step
	// already converged, so the update keeps x + alpha p and r = s.
	if (!isfinite(scalars->tt)) scalars->breakdown = 1;
	else scalars->omega = scalars->tt == 0.0 ? 0.0 : scalars->ts/scalars->tt;
}

__global__ void BiCgStabUpdateKernel(int size, const KrylovScalars* scalars, const double* p, const double* s,
	const double* t, double* x, double* r)
{
	const int row = blockIdx.x*blockDim.x+threadIdx.x;
	if (row >= size) return;
	x[row] += scalars->alpha*p[row]+scalars->omega*s[row];
	r[row] = s[row]-scalars->omega*t[row];
}

__global__ void BiCgStabStartKernel(KrylovScalars* scalars)
{
	scalars->rho_previous = 1.0;
	scalars->alpha = 1.0;
	scalars->omega = 1.0;
}

// Left Jacobi-preconditioned BiCGStab for a general matrix: the iteration runs
// on D^-1 A, so its residual is the preconditioned residual.
inline KrylovResult SolveBiCgStab(const BlockMatrix<1>& matrix, ScalarKrylovWorkspace& workspace,
	const double* rhs, double* solution, int maximum_iterations, double relative_tolerance, double absolute_tolerance)
{
	const auto pattern = matrix.pattern();
	const int size = pattern.nodes, blocks = (size+255)/256;
	if (workspace.size != size) throw std::invalid_argument("Krylov workspace size does not match matrix");
	SetupJacobi(matrix, workspace);
	auto& blas = workspace.blas;
	auto* scalars = workspace.scalars.data();
	const double* inverse = workspace.inverse.data();
	double *r = workspace.r.data(), *p = workspace.p.data(), *v = workspace.v.data(), *s = workspace.s.data();
	double *t = workspace.t.data(), *shadow = workspace.shadow.data();
	const auto bytes = static_cast<std::size_t>(size)*sizeof(double);
	workspace.scalars.Clear();
	ScaleKernel<<<blocks,256>>>(size, inverse, rhs, s);
	Check(cublasDnrm2(blas, size, s, 1, &scalars->load_norm), "cublasDnrm2 BiCGStab load");
	ScaledResidualKernel<<<blocks,256>>>(pattern, matrix.values(), inverse, rhs, solution, r);
	Check(cublasDnrm2(blas, size, r, 1, &scalars->norm), "cublasDnrm2 BiCGStab initial");
	Check(cudaMemcpy(shadow, r, bytes, cudaMemcpyDeviceToDevice), "copy BiCGStab shadow residual");
	Check(cudaMemset(p, 0, bytes), "clear BiCGStab direction");
	Check(cudaMemset(v, 0, bytes), "clear BiCGStab v");
	BiCgStabStartKernel<<<1,1>>>(scalars);
	CheckKernel("BiCGStab setup");
	auto host = ReadScalars(workspace);
	KrylovResult result;
	result.tolerance = std::max(relative_tolerance*host.load_norm, absolute_tolerance);
	result.residual = host.norm;
	if (host.norm <= result.tolerance) {
		result.converged = true;
		return result;
	}
	while (result.iterations < maximum_iterations) {
		Check(cublasDdot(blas, size, shadow, 1, r, 1, &scalars->rho), "cublasDdot BiCGStab rho");
		BiCgStabDirectionKernel<<<blocks,256>>>(size, scalars, r, v, p);
		ScaledSpmvKernel<<<blocks,256>>>(pattern, matrix.values(), inverse, p, v);
		Check(cublasDdot(blas, size, shadow, 1, v, 1, &scalars->curvature), "cublasDdot BiCGStab shadow.v");
		BiCgStabAlphaKernel<<<1,1>>>(scalars);
		BiCgStabHalfStepKernel<<<blocks,256>>>(size, scalars, r, v, s);
		ScaledSpmvKernel<<<blocks,256>>>(pattern, matrix.values(), inverse, s, t);
		Check(cublasDdot(blas, size, t, 1, s, 1, &scalars->ts), "cublasDdot BiCGStab t.s");
		Check(cublasDdot(blas, size, t, 1, t, 1, &scalars->tt), "cublasDdot BiCGStab t.t");
		BiCgStabOmegaKernel<<<1,1>>>(scalars);
		BiCgStabUpdateKernel<<<blocks,256>>>(size, scalars, p, s, t, solution, r);
		Check(cublasDnrm2(blas, size, r, 1, &scalars->norm), "cublasDnrm2 BiCGStab");
		ShiftScalarKernel<<<1,1>>>(scalars);
		++result.iterations;
		if (result.iterations % kCheckInterval && result.iterations < maximum_iterations) continue;
		CheckKernel("BiCGStab iteration");
		host = ReadScalars(workspace);
		result.residual = host.norm;
		if (host.breakdown) throw std::runtime_error("BiCGStab broke down");
		if (!std::isfinite(host.norm)) break;
		if (host.norm <= result.tolerance) {
			result.converged = true;
			break;
		}
	}
	return result;
}

} // namespace iga::cuda

#endif
