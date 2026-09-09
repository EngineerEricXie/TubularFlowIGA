#ifndef IGA_TEST_ABORT_REPORT_ALLOCATION_HPP
#define IGA_TEST_ABORT_REPORT_ALLOCATION_HPP

#include <cstdlib>
#include <cstring>
#include <new>

namespace abort_report_allocation {
inline thread_local bool active = false, injected = false;
inline thread_local std::size_t allocations = 0, fail_at = 0;

inline void Arm(std::size_t ordinal)
{
	allocations = 0;
	fail_at = ordinal;
	injected = false;
	active = true;
}
inline void Disarm() { active = false; }
}

// This header is used by single-translation-unit test executables only.
// Arm after the final abort outcome, disarm before report outcome agreement.
[[gnu::noinline]] void* operator new(std::size_t size)
{
	if (abort_report_allocation::active
		&& ++abort_report_allocation::allocations == abort_report_allocation::fail_at) {
		abort_report_allocation::active = false;
		abort_report_allocation::injected = true;
		throw std::bad_alloc();
	}
	if (auto* value = std::malloc(size ? size : 1)) return value;
	throw std::bad_alloc();
}

[[gnu::noinline]] void operator delete(void* value) noexcept { std::free(value); }
[[gnu::noinline]] void operator delete(void* value, std::size_t) noexcept { std::free(value); }

#endif
