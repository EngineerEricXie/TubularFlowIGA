#ifndef IGA_DEVICE_ALLOCATION_HPP
#define IGA_DEVICE_ALLOCATION_HPP

#include <atomic>
#include <cstddef>

namespace iga::cuda {

// Requested bytes owned by DeviceBuffer, excluding driver-wide memory usage
// and allocations internal to cuBLAS. A move transfers existing ownership.
class DeviceAllocationCounter {
public:
	static void Add(std::size_t bytes) noexcept
	{
		const auto current = current_.fetch_add(bytes) + bytes;
		auto peak = peak_.load();
		while (peak < current && !peak_.compare_exchange_weak(peak, current)) {}
	}

	static void Remove(std::size_t bytes) noexcept { current_.fetch_sub(bytes); }
	static std::size_t Current() noexcept { return current_.load(); }
	static std::size_t Peak() noexcept { return peak_.load(); }

private:
	inline static std::atomic<std::size_t> current_{0};
	inline static std::atomic<std::size_t> peak_{0};
};

} // namespace iga::cuda

#endif
