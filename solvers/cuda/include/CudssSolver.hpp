#ifndef IGA_CUDA_CUDSS_SOLVER_HPP
#define IGA_CUDA_CUDSS_SOLVER_HPP

#include "BlockCsr.hpp"
#include "CudaRuntime.hpp"

#include <cudss.h>

#include <cstddef>
#include <cstdio>
#include <stdexcept>
#include <string>

namespace iga::cuda {

inline void Check(cudssStatus_t status, const char* operation)
{
	if (status != CUDSS_STATUS_SUCCESS)
		throw std::runtime_error(std::string(operation)+" failed with cuDSS status "
			+std::to_string(static_cast<int>(status)));
}

// Route cuDSS workspace through the DeviceBuffer allocation counter so the
// reported CUDA peak includes the sparse factors.
inline int CountedCudssAllocate(void*, void** pointer, std::size_t bytes, cudaStream_t)
{
	if (cudaMalloc(pointer, bytes) != cudaSuccess) return 1;
	DeviceAllocationCounter::Add(bytes);
	return 0;
}

inline int CountedCudssFree(void*, void* pointer, std::size_t bytes, cudaStream_t)
{
	DeviceAllocationCounter::Remove(bytes);
	return cudaFree(pointer) == cudaSuccess ? 0 : 1;
}

// Sparse direct LU (cuDSS) of a scalar CSR matrix whose pattern stays fixed.
// The first Solve reorders and factorizes symbolically; later solves only
// refactorize the current values, which the caller updates in place.
class CudssSolver {
public:
	explicit CudssSolver(BlockMatrix<1>& matrix)
		: size_(matrix.pattern().nodes), right_hand_side_(static_cast<std::size_t>(size_)),
		  solution_(static_cast<std::size_t>(size_))
	{
		Check(cudssCreate(&handle_), "cudssCreate");
		cudssDeviceMemHandler_t handler{};
		handler.device_alloc = CountedCudssAllocate;
		handler.device_free = CountedCudssFree;
		std::snprintf(handler.name, CUDSS_ALLOCATOR_NAME_LEN, "TubularFlowIGA counted");
		Check(cudssSetDeviceMemHandler(handle_, &handler), "cudssSetDeviceMemHandler");
		Check(cudssConfigCreate(&config_), "cudssConfigCreate");
		Check(cudssDataCreate(handle_, &data_), "cudssDataCreate");
		const auto pattern = matrix.pattern();
		Check(cudssMatrixCreateCsr(&matrix_, size_, size_, pattern.blocks,
			const_cast<int*>(pattern.row_offsets), nullptr, const_cast<int*>(pattern.columns),
			matrix.values(), CUDSS_R_32I, CUDSS_R_32I, CUDSS_R_64F, CUDSS_MTYPE_GENERAL,
			CUDSS_MVIEW_FULL, CUDSS_BASE_ZERO), "cudssMatrixCreateCsr");
		Check(cudssMatrixCreateDn(&right_hand_side_matrix_, size_, 1, size_, right_hand_side_.data(),
			CUDSS_R_64F, CUDSS_LAYOUT_COL_MAJOR), "cudssMatrixCreateDn rhs");
		Check(cudssMatrixCreateDn(&solution_matrix_, size_, 1, size_, solution_.data(),
			CUDSS_R_64F, CUDSS_LAYOUT_COL_MAJOR), "cudssMatrixCreateDn solution");
	}

	~CudssSolver()
	{
		if (solution_matrix_) cudssMatrixDestroy(solution_matrix_);
		if (right_hand_side_matrix_) cudssMatrixDestroy(right_hand_side_matrix_);
		if (matrix_) cudssMatrixDestroy(matrix_);
		if (data_) cudssDataDestroy(handle_, data_);
		if (config_) cudssConfigDestroy(config_);
		if (handle_) cudssDestroy(handle_);
	}

	CudssSolver(const CudssSolver&) = delete;
	CudssSolver& operator=(const CudssSolver&) = delete;

	// Solve A x = b for the current matrix values; b is overwritten with x.
	void Solve(double* right_hand_side)
	{
		const auto bytes = static_cast<std::size_t>(size_)*sizeof(double);
		Check(cudaMemcpy(right_hand_side_.data(), right_hand_side, bytes, cudaMemcpyDeviceToDevice),
			"copy cuDSS right-hand side");
		if (!analysed_) {
			Execute(CUDSS_PHASE_ANALYSIS, "cuDSS analysis");
			Execute(CUDSS_PHASE_FACTORIZATION, "cuDSS factorization");
			analysed_ = true;
		} else Execute(CUDSS_PHASE_REFACTORIZATION, "cuDSS refactorization");
		int info = 0;
		std::size_t written = 0;
		Check(cudssDataGet(handle_, data_, CUDSS_DATA_INFO, &info, sizeof(info), &written),
			"cudssDataGet info");
		if (info != 0) throw std::runtime_error("cuDSS factorization failed with info "+std::to_string(info));
		Execute(CUDSS_PHASE_SOLVE, "cuDSS solve");
		Check(cudaMemcpy(right_hand_side, solution_.data(), bytes, cudaMemcpyDeviceToDevice),
			"copy cuDSS solution");
	}

private:
	void Execute(cudssPhase_t phase, const char* operation)
	{
		Check(cudssExecute(handle_, phase, config_, data_, matrix_, solution_matrix_, right_hand_side_matrix_),
			operation);
	}

	int size_;
	bool analysed_ = false;
	cudssHandle_t handle_ = nullptr;
	cudssConfig_t config_ = nullptr;
	cudssData_t data_ = nullptr;
	cudssMatrix_t matrix_ = nullptr, right_hand_side_matrix_ = nullptr, solution_matrix_ = nullptr;
	DeviceBuffer<double> right_hand_side_, solution_;
};

} // namespace iga::cuda

#endif
