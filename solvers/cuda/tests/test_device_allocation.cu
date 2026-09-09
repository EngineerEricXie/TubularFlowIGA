#include "CudaRuntime.hpp"

#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>

namespace {

void Require(bool condition)
{
	if (!condition) throw std::runtime_error("device allocation ownership accounting failed");
}

} // namespace

int main()
{
	try {
		using iga::cuda::DeviceAllocationCounter;
		using iga::cuda::DeviceBuffer;
		Require(DeviceAllocationCounter::Current() == 0);
		{
			DeviceBuffer<double> first(10);
			DeviceBuffer<int> second(5);
			Require(DeviceAllocationCounter::Current() == 100);
			DeviceBuffer<double> moved(std::move(first));
			Require(first.size() == 0 && moved.size() == 10);
			second = DeviceBuffer<int>(2);
			Require(DeviceAllocationCounter::Current() == 88);
			Require(DeviceAllocationCounter::Peak() == 108);
			bool rejected = false;
			try { first.Allocate(std::numeric_limits<std::size_t>::max()); }
			catch (const std::overflow_error&) { rejected = true; }
			Require(rejected && first.size() == 0 && DeviceAllocationCounter::Current() == 88);
			DeviceBuffer<double> empty(0);
			Require(empty.data() == nullptr && DeviceAllocationCounter::Current() == 88);
		}
		Require(DeviceAllocationCounter::Current() == 0 && DeviceAllocationCounter::Peak() == 108);
		std::cout << "device_allocation_test passed peak_requested_bytes=108 live_bytes=0\n";
		return 0;
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n';
		return 1;
	}
}
