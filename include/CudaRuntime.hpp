#ifndef IGA_CUDA_RUNTIME_HPP
#define IGA_CUDA_RUNTIME_HPP

#include <cublas_v2.h>
#include <cuda_runtime.h>

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>

namespace iga::cuda {

inline void Check(cudaError_t code, const char* operation)
{
	if (code != cudaSuccess)
		throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(code));
}

inline void Check(cublasStatus_t code, const char* operation)
{
	if (code != CUBLAS_STATUS_SUCCESS)
		throw std::runtime_error(std::string(operation) + ": cuBLAS status " + std::to_string(code));
}

inline void CheckKernel(const char* operation)
{
	Check(cudaGetLastError(), operation);
}

template <class T>
class DeviceBuffer {
public:
	DeviceBuffer() = default;
	explicit DeviceBuffer(std::size_t count) { Allocate(count); }
	~DeviceBuffer() { if (data_) cudaFree(data_); }

	DeviceBuffer(const DeviceBuffer&) = delete;
	DeviceBuffer& operator=(const DeviceBuffer&) = delete;

	DeviceBuffer(DeviceBuffer&& other) noexcept
		: data_(std::exchange(other.data_, nullptr)), size_(std::exchange(other.size_, 0)) {}

	DeviceBuffer& operator=(DeviceBuffer&& other) noexcept
	{
		if (this != &other) {
			if (data_) cudaFree(data_);
			data_ = std::exchange(other.data_, nullptr);
			size_ = std::exchange(other.size_, 0);
		}
		return *this;
	}

	void Allocate(std::size_t count)
	{
		if (data_) throw std::logic_error("DeviceBuffer is already allocated");
		size_ = count;
		if (count) Check(cudaMalloc(reinterpret_cast<void**>(&data_), count * sizeof(T)), "cudaMalloc");
	}

	void Clear()
	{
		if (data_) Check(cudaMemset(data_, 0, size_ * sizeof(T)), "cudaMemset");
	}

	void CopyFromHost(const T* source, std::size_t count)
	{
		if (count > size_) throw std::out_of_range("host-to-device copy exceeds buffer");
		if (count) Check(cudaMemcpy(data_, source, count * sizeof(T), cudaMemcpyHostToDevice), "cudaMemcpy H2D");
	}

	void CopyToHost(T* destination, std::size_t count) const
	{
		if (count > size_) throw std::out_of_range("device-to-host copy exceeds buffer");
		if (count) Check(cudaMemcpy(destination, data_, count * sizeof(T), cudaMemcpyDeviceToHost), "cudaMemcpy D2H");
	}

	T* data() { return data_; }
	const T* data() const { return data_; }
	std::size_t size() const { return size_; }
	std::size_t bytes() const { return size_ * sizeof(T); }

private:
	T* data_ = nullptr;
	std::size_t size_ = 0;
};

class BlasHandle {
public:
	BlasHandle() { Check(cublasCreate(&handle_), "cublasCreate"); }
	~BlasHandle() { if (handle_) cublasDestroy(handle_); }
	BlasHandle(const BlasHandle&) = delete;
	BlasHandle& operator=(const BlasHandle&) = delete;
	operator cublasHandle_t() const { return handle_; }

private:
	cublasHandle_t handle_ = nullptr;
};

inline double Gibibytes(std::size_t bytes)
{
	return static_cast<double>(bytes) / (1024.0 * 1024.0 * 1024.0);
}

} // namespace iga::cuda

#endif
