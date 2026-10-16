#ifndef IGA_CUDA_SADDLE_POINT_SOLVER_HPP
#define IGA_CUDA_SADDLE_POINT_SOLVER_HPP

#include "AggregationAmg.hpp"
#include "ScalarKrylov.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// Iterative solver for the P2/P1 Newton system [A B'; B C] with the velocity
// unknowns first. FGMRES uses a block upper-triangular right preconditioner
// with PETSc's "selfp" Schur approximation S = C - B diag(A)^-1 B':
//   y_p = S^-1 r_p,  y_u = A^-1 (r_u - B' y_p),
// where both inner solves are left-Jacobi BiCGStab to a loose tolerance. This
// is the GPU counterpart of the CPU fieldsplit configuration (full Schur,
// selfp) without its algebraic multigrid, so it suits mass-dominated
// transient steps.

namespace iga::cuda {

struct SaddlePointOptions {
	double relative_tolerance = 1e-10;
	int maximum_iterations = 1000;
	int restart = 50;
	double velocity_relative_tolerance = 1e-3;
	double pressure_relative_tolerance = 1e-3;
	int inner_maximum_iterations = 400;
	// Inner solves: "jacobi" (left-Jacobi BiCGStab), "amg" (BiCGStab
	// preconditioned by aggregation multigrid on the velocity block and on S),
	// or "amg-velocity" / "amg-pressure" for multigrid on one block only.
	std::string preconditioner = "jacobi";
	int amg_cycles = 1;
	AggregationAmgOptions amg;
};

struct SaddlePointResult {
	int iterations = 0;
	double relative_residual = 0.0;
	bool converged = false;
};

// y[row-row_begin] = sum over columns in [column_begin, column_end) of
// J(row, column) x[column-column_begin], for rows in [row_begin, row_end).
__global__ void RangeSpmvKernel(DevicePatternView pattern, const double* values, int row_begin, int row_end,
	int column_begin, int column_end, const double* x, double* y)
{
	const int row = row_begin+blockIdx.x*blockDim.x+threadIdx.x;
	if (row >= row_end) return;
	double sum = 0.0;
	for (int entry = pattern.row_offsets[row]; entry < pattern.row_offsets[row+1]; ++entry) {
		const int column = pattern.columns[entry];
		if (column >= column_begin && column < column_end) sum += values[entry]*x[column-column_begin];
	}
	y[row-row_begin] = sum;
}

__global__ void GatherEntriesKernel(int count, const int* source, const double* values, double* destination)
{
	const int index = blockIdx.x*blockDim.x+threadIdx.x;
	if (index < count) destination[index] = values[source[index]];
}

// S(i,:) = C(i,:) - sum_k B(i,k) / A(k,k) * B'(k,:) for pressure row i; J holds
// all four blocks, `velocity` is the number of velocity unknowns.
__global__ void SelfpSchurKernel(DevicePatternView jacobian, const double* values, int velocity,
	const double* velocity_diagonal_inverse, DevicePatternView schur, double* schur_values)
{
	const int row = blockIdx.x*blockDim.x+threadIdx.x;
	if (row >= schur.nodes) return;
	for (int entry = schur.row_offsets[row]; entry < schur.row_offsets[row+1]; ++entry) schur_values[entry] = 0.0;
	const int jacobian_row = velocity+row;
	for (int entry = jacobian.row_offsets[jacobian_row]; entry < jacobian.row_offsets[jacobian_row+1]; ++entry) {
		const int column = jacobian.columns[entry];
		if (column >= velocity) {
			schur_values[FindBlock(schur, row, column-velocity)] += values[entry];
			continue;
		}
		const double scale = values[entry]*velocity_diagonal_inverse[column];
		for (int inner = jacobian.row_offsets[column]; inner < jacobian.row_offsets[column+1]; ++inner) {
			const int pressure = jacobian.columns[inner];
			if (pressure >= velocity)
				schur_values[FindBlock(schur, row, pressure-velocity)] -= scale*values[inner];
		}
	}
}

__global__ void SubtractKernel(int size, const double* x, double* y)
{
	const int index = blockIdx.x*blockDim.x+threadIdx.x;
	if (index < size) y[index] -= x[index];
}

class SaddlePointSolver {
public:
	// `host` is the pattern of `matrix` (the full Newton Jacobian).
	SaddlePointSolver(const BlockPattern& host, BlockMatrix<1>& matrix, int velocity, SaddlePointOptions options = {})
		: matrix_(matrix), size_(static_cast<int>(host.diagonal.size())), velocity_(velocity),
		  pressure_(size_-velocity), options_(options),
		  velocity_host_(BuildVelocityPattern(host, velocity)), velocity_pattern_(velocity_host_),
		  velocity_block_(velocity_pattern_), schur_host_(BuildSchurPattern(host, velocity)),
		  schur_pattern_(schur_host_), schur_(schur_pattern_),
		  velocity_source_(velocity_entries_.size()), velocity_workspace_(velocity), schur_workspace_(size_-velocity),
		  basis_(static_cast<std::size_t>(options.restart+1)*size_), preconditioned_(static_cast<std::size_t>(options.restart)*size_),
		  residual_(size_), work_(size_), solution_(size_), velocity_rhs_(velocity), product_(size_),
		  projection_(options.restart+1)
	{
		if (velocity <= 0 || pressure_ <= 0 || options.restart < 1)
			throw std::invalid_argument("saddle-point solver needs velocity and pressure unknowns");
		velocity_source_.CopyFromHost(velocity_entries_.data(), velocity_entries_.size());
		Check(cublasSetPointerMode(blas_, CUBLAS_POINTER_MODE_HOST), "cublasSetPointerMode host");
		const bool amg_velocity = options.preconditioner == "amg" || options.preconditioner == "amg-velocity";
		const bool amg_pressure = options.preconditioner == "amg" || options.preconditioner == "amg-pressure";
		if (amg_velocity) velocity_amg_ = std::make_unique<AggregationAmg>(velocity_host_, velocity_block_, 3, options.amg);
		if (amg_pressure) pressure_amg_ = std::make_unique<AggregationAmg>(schur_host_, schur_, 1, options.amg);
		if (!amg_velocity && !amg_pressure && options.preconditioner != "jacobi")
			throw std::invalid_argument("saddle-point preconditioner must be jacobi, amg, amg-velocity, or amg-pressure");
	}

	int VelocityLevels() const { return velocity_amg_ ? velocity_amg_->Levels() : 1; }
	int AmgFallbacks() const
	{
		return (velocity_amg_ ? velocity_amg_->FallbackCount() : 0)+(pressure_amg_ ? pressure_amg_->FallbackCount() : 0);
	}
	int PressureLevels() const { return pressure_amg_ ? pressure_amg_->Levels() : 1; }
	int AmgGraphs() const
	{
		return (velocity_amg_ ? velocity_amg_->GraphCount() : 0)+(pressure_amg_ ? pressure_amg_->GraphCount() : 0);
	}

	// Solves J x = b for the current Jacobian values; b is overwritten with x.
	SaddlePointResult Solve(double* right_hand_side)
	{
		UpdateBlocks();
		const auto bytes = static_cast<std::size_t>(size_)*sizeof(double);
		SaddlePointResult result;
		double load = 0.0;
		Check(cublasDnrm2(blas_, size_, right_hand_side, 1, &load), "cublasDnrm2 FGMRES load");
		Check(cudaMemset(solution_.data(), 0, bytes), "clear FGMRES solution");
		if (load == 0.0) {
			result.converged = true;
			return result;
		}
		const double tolerance = options_.relative_tolerance*load;
		const int blocks = (size_+255)/256, restart = options_.restart;
		std::vector<double> h(static_cast<std::size_t>(restart+1)*restart), cosine(restart), sine(restart), g(restart+1);
		while (result.iterations < options_.maximum_iterations) {
			// r = b - J x
			BlockSpmvKernel<1><<<blocks,256>>>(matrix_.pattern(), matrix_.values(), solution_.data(), product_.data());
			Check(cudaMemcpy(residual_.data(), right_hand_side, bytes, cudaMemcpyDeviceToDevice), "copy FGMRES rhs");
			SubtractKernel<<<blocks,256>>>(size_, product_.data(), residual_.data());
			double beta = 0.0;
			Check(cublasDnrm2(blas_, size_, residual_.data(), 1, &beta), "cublasDnrm2 FGMRES residual");
			result.relative_residual = beta/load;
			if (beta <= tolerance) {
				result.converged = true;
				break;
			}
			std::fill(h.begin(), h.end(), 0.0);
			std::fill(g.begin(), g.end(), 0.0);
			g[0] = beta;
			const double inverse_beta = 1.0/beta;
			Check(cudaMemcpy(basis_.data(), residual_.data(), bytes, cudaMemcpyDeviceToDevice), "copy FGMRES v0");
			Check(cublasDscal(blas_, size_, &inverse_beta, basis_.data(), 1), "cublasDscal FGMRES v0");
			int used = 0;
			for (int j = 0; j < restart && result.iterations < options_.maximum_iterations; ++j) {
				double* v = basis_.data()+static_cast<std::size_t>(j)*size_;
				double* z = preconditioned_.data()+static_cast<std::size_t>(j)*size_;
				ApplyPreconditioner(v, z);
				double* w = basis_.data()+static_cast<std::size_t>(j+1)*size_;
				BlockSpmvKernel<1><<<blocks,256>>>(matrix_.pattern(), matrix_.values(), z, w);
				// Classical Gram-Schmidt, applied twice for stability.
				const double one = 1.0, zero = 0.0, minus_one = -1.0;
				for (int pass = 0; pass < 2; ++pass) {
					Check(cublasDgemv(blas_, CUBLAS_OP_T, size_, j+1, &one, basis_.data(), size_, w, 1, &zero,
						projection_.data(), 1), "cublasDgemv FGMRES projection");
					Check(cublasDgemv(blas_, CUBLAS_OP_N, size_, j+1, &minus_one, basis_.data(), size_,
						projection_.data(), 1, &one, w, 1), "cublasDgemv FGMRES orthogonalize");
					std::vector<double> projection(static_cast<std::size_t>(j+1));
					projection_.CopyToHost(projection.data(), projection.size());
					for (int i = 0; i <= j; ++i) h[static_cast<std::size_t>(j)*(restart+1)+i] += projection[i];
				}
				double norm = 0.0;
				Check(cublasDnrm2(blas_, size_, w, 1, &norm), "cublasDnrm2 FGMRES w");
				h[static_cast<std::size_t>(j)*(restart+1)+j+1] = norm;
				if (norm > 0.0) {
					const double inverse = 1.0/norm;
					Check(cublasDscal(blas_, size_, &inverse, w, 1), "cublasDscal FGMRES w");
				}
				double* hj = h.data()+static_cast<std::size_t>(j)*(restart+1);
				for (int i = 0; i < j; ++i) {
					const double first = cosine[i]*hj[i]+sine[i]*hj[i+1];
					hj[i+1] = -sine[i]*hj[i]+cosine[i]*hj[i+1];
					hj[i] = first;
				}
				const double radius = std::hypot(hj[j], hj[j+1]);
				cosine[j] = radius == 0.0 ? 1.0 : hj[j]/radius;
				sine[j] = radius == 0.0 ? 0.0 : hj[j+1]/radius;
				hj[j] = radius;
				hj[j+1] = 0.0;
				g[j+1] = -sine[j]*g[j];
				g[j] = cosine[j]*g[j];
				++used;
				++result.iterations;
				result.relative_residual = std::abs(g[j+1])/load;
				if (std::abs(g[j+1]) <= tolerance || norm == 0.0) break;
			}
			// x += Z y with H y = g (upper triangular).
			std::vector<double> y(static_cast<std::size_t>(used));
			for (int i = used-1; i >= 0; --i) {
				double value = g[i];
				for (int k = i+1; k < used; ++k) value -= h[static_cast<std::size_t>(k)*(restart+1)+i]*y[k];
				y[i] = value/h[static_cast<std::size_t>(i)*(restart+1)+i];
			}
			projection_.CopyFromHost(y.data(), y.size());
			const double one = 1.0;
			Check(cublasDgemv(blas_, CUBLAS_OP_N, size_, used, &one, preconditioned_.data(), size_,
				projection_.data(), 1, &one, solution_.data(), 1), "cublasDgemv FGMRES update");
		}
		if (!result.converged) {
			// Confirm the final residual of the last restart cycle.
			BlockSpmvKernel<1><<<blocks,256>>>(matrix_.pattern(), matrix_.values(), solution_.data(), product_.data());
			Check(cudaMemcpy(residual_.data(), right_hand_side, bytes, cudaMemcpyDeviceToDevice), "copy FGMRES rhs");
			SubtractKernel<<<blocks,256>>>(size_, product_.data(), residual_.data());
			double beta = 0.0;
			Check(cublasDnrm2(blas_, size_, residual_.data(), 1, &beta), "cublasDnrm2 FGMRES final");
			result.relative_residual = beta/load;
			result.converged = beta <= tolerance;
		}
		CheckKernel("FGMRES");
		Check(cudaMemcpy(right_hand_side, solution_.data(), bytes, cudaMemcpyDeviceToDevice), "copy FGMRES solution");
		if (!result.converged)
			throw std::runtime_error("saddle-point FGMRES did not converge: iterations="
				+std::to_string(result.iterations)+" relative_residual="+std::to_string(result.relative_residual));
		return result;
	}

	int VelocityIterations() const { return velocity_iterations_; }
	int PressureIterations() const { return pressure_iterations_; }

private:
	BlockPattern BuildVelocityPattern(const BlockPattern& host, int velocity)
	{
		BlockPattern pattern;
		pattern.row_offsets.push_back(0);
		pattern.diagonal.assign(static_cast<std::size_t>(velocity), -1);
		for (int row = 0; row < velocity; ++row) {
			for (int entry = host.row_offsets[row]; entry < host.row_offsets[row+1]; ++entry) {
				const int column = host.columns[entry];
				if (column >= velocity) continue;
				if (column == row) pattern.diagonal[row] = static_cast<int>(pattern.columns.size());
				pattern.columns.push_back(column);
				velocity_entries_.push_back(entry);
			}
			if (pattern.diagonal[row] < 0) throw std::runtime_error("velocity block lacks a diagonal entry");
			pattern.row_offsets.push_back(static_cast<int>(pattern.columns.size()));
		}
		return pattern;
	}

	static BlockPattern BuildSchurPattern(const BlockPattern& host, int velocity)
	{
		const int size = static_cast<int>(host.diagonal.size()), pressure = size-velocity;
		BlockPattern pattern;
		pattern.row_offsets.push_back(0);
		pattern.diagonal.assign(static_cast<std::size_t>(pressure), -1);
		std::vector<int> row;
		for (int p = 0; p < pressure; ++p) {
			row.assign(1, p);
			const int jacobian_row = velocity+p;
			for (int entry = host.row_offsets[jacobian_row]; entry < host.row_offsets[jacobian_row+1]; ++entry) {
				const int column = host.columns[entry];
				if (column >= velocity) {
					row.push_back(column-velocity);
					continue;
				}
				for (int inner = host.row_offsets[column]; inner < host.row_offsets[column+1]; ++inner)
					if (host.columns[inner] >= velocity) row.push_back(host.columns[inner]-velocity);
			}
			std::sort(row.begin(), row.end());
			row.erase(std::unique(row.begin(), row.end()), row.end());
			pattern.diagonal[p] = static_cast<int>(pattern.columns.size()
				+(std::lower_bound(row.begin(), row.end(), p)-row.begin()));
			pattern.columns.insert(pattern.columns.end(), row.begin(), row.end());
			pattern.row_offsets.push_back(static_cast<int>(pattern.columns.size()));
		}
		return pattern;
	}

	void UpdateBlocks()
	{
		const int count = static_cast<int>(velocity_entries_.size());
		GatherEntriesKernel<<<(count+255)/256,256>>>(count, velocity_source_.data(), matrix_.values(),
			velocity_block_.values());
		SetupJacobi(velocity_block_, velocity_workspace_);
		SelfpSchurKernel<<<(pressure_+255)/256,256>>>(matrix_.pattern(), matrix_.values(), velocity_,
			velocity_workspace_.inverse.data(), schur_pattern_.view(), schur_.values());
		CheckKernel("saddle-point blocks");
		if (velocity_amg_) velocity_amg_->Update();
		if (pressure_amg_) pressure_amg_->Update();
	}

	void ApplyPreconditioner(const double* v, double* z)
	{
		double* z_velocity = z;
		double* z_pressure = z+velocity_;
		Check(cudaMemset(z, 0, static_cast<std::size_t>(size_)*sizeof(double)), "clear preconditioned vector");
		if (pressure_amg_) {
			const auto pressure = SolvePreconditionedBiCgStab(schur_, schur_workspace_, v+velocity_, z_pressure,
				options_.inner_maximum_iterations, options_.pressure_relative_tolerance, 0.0,
				[&](const double* in, double* out) { pressure_amg_->Cycle(in, out, options_.amg_cycles); });
			pressure_iterations_ += pressure.iterations;
		} else {
			const auto pressure = SolveBiCgStab(schur_, schur_workspace_, v+velocity_, z_pressure,
				options_.inner_maximum_iterations, options_.pressure_relative_tolerance, 0.0);
			pressure_iterations_ += pressure.iterations;
		}
		// velocity rhs = r_u - B' y_p
		RangeSpmvKernel<<<(velocity_+255)/256,256>>>(matrix_.pattern(), matrix_.values(), 0, velocity_,
			velocity_, size_, z_pressure, velocity_rhs_.data());
		Check(cudaMemcpy(work_.data(), v, static_cast<std::size_t>(velocity_)*sizeof(double),
			cudaMemcpyDeviceToDevice), "copy velocity residual");
		SubtractKernel<<<(velocity_+255)/256,256>>>(velocity_, velocity_rhs_.data(), work_.data());
		if (velocity_amg_) {
			const auto velocity = SolvePreconditionedBiCgStab(velocity_block_, velocity_workspace_, work_.data(),
				z_velocity, options_.inner_maximum_iterations, options_.velocity_relative_tolerance, 0.0,
				[&](const double* in, double* out) { velocity_amg_->Cycle(in, out, options_.amg_cycles); });
			velocity_iterations_ += velocity.iterations;
		} else {
			const auto velocity = SolveBiCgStab(velocity_block_, velocity_workspace_, work_.data(), z_velocity,
				options_.inner_maximum_iterations, options_.velocity_relative_tolerance, 0.0);
			velocity_iterations_ += velocity.iterations;
		}
	}

	BlockMatrix<1>& matrix_;
	int size_, velocity_, pressure_;
	SaddlePointOptions options_;
	std::vector<int> velocity_entries_;
	BlockPattern velocity_host_;
	DevicePattern velocity_pattern_;
	BlockMatrix<1> velocity_block_;
	BlockPattern schur_host_;
	DevicePattern schur_pattern_;
	BlockMatrix<1> schur_;
	DeviceBuffer<int> velocity_source_;
	ScalarKrylovWorkspace velocity_workspace_, schur_workspace_;
	BlasHandle blas_;
	DeviceBuffer<double> basis_, preconditioned_, residual_, work_, solution_, velocity_rhs_, product_;
	DeviceBuffer<double> projection_;
	int velocity_iterations_ = 0, pressure_iterations_ = 0;
	std::unique_ptr<AggregationAmg> velocity_amg_, pressure_amg_;
};

} // namespace iga::cuda

#endif
