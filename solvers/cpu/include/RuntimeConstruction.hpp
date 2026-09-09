#ifndef IGA_RUNTIME_CONSTRUCTION_HPP
#define IGA_RUNTIME_CONSTRUCTION_HPP

#include "CollectiveFailure.hpp"
#include <cstddef>
#include <memory>
#include <new>
#include <optional>
#include <type_traits>
#include <utility>

#ifdef IGA_RUNTIME_CONSTRUCTION_TESTING
#include <functional>
#endif

namespace iga {

#ifdef IGA_RUNTIME_CONSTRUCTION_TESTING
struct RuntimeConstructionHooks {
	std::function<void(const char*)> after_local_stage;
	std::function<void(PetscObject)> object_created;
};

inline RuntimeConstructionHooks& RuntimeConstructionHooksForTesting()
{
	static thread_local RuntimeConstructionHooks hooks;
	return hooks;
}
#endif

inline void ProbeRuntimeConstruction(const char* stage)
{
#ifdef IGA_RUNTIME_CONSTRUCTION_TESTING
	const auto& probe = RuntimeConstructionHooksForTesting().after_local_stage;
	if (probe) probe(stage);
#else
	(void)stage;
#endif
}

// Like CollectiveLocalStage, work must contain only local operations. The
// test hook runs after the work and before the common failure decision.
template<class Work>
void RuntimeConstructionStage(MPI_Comm communicator, const char* stage, Work&& work)
{
	CollectiveLocalStage(communicator, stage, [&] {
		std::forward<Work>(work)();
		ProbeRuntimeConstruction(stage);
	});
}

// Used before a member initializer that enters MPI. Copy allocation belongs
// inside the local stage; returning the prepared value cannot throw.
template<class Value, class... Arguments>
Value PrepareRuntimeConstructionInput(MPI_Comm communicator, const char* stage, Arguments&&... arguments)
{
	static_assert(std::is_nothrow_move_constructible<Value>::value,
		"prepared construction input must move without throwing");
	std::optional<Value> result;
	RuntimeConstructionStage(communicator, stage, [&] {
		result.emplace(std::forward<Arguments>(arguments)...);
	});
	return std::move(*result);
}

// Allocate host storage collectively before entering a collective constructor.
// Arguments are forwarded by reference: callers must prepare any allocating
// expressions first, and the constructor must coordinate its own input copies.
// This factory is for ordinary host runtime types without class-specific
// operator new/delete. Their normal delete expression owns successful objects.
template<class Runtime, class... Arguments>
std::unique_ptr<Runtime> AllocateCollectiveRuntime(MPI_Comm communicator, Arguments&&... arguments)
{
	static_assert(alignof(Runtime) <= alignof(std::max_align_t),
		"collective runtime storage requires an ordinary host alignment");
	auto release_storage = [](void* pointer) noexcept { ::operator delete(pointer); };
	std::unique_ptr<void, decltype(release_storage)> storage(nullptr, release_storage);
	RuntimeConstructionStage(communicator, "runtime object storage ready", [&] {
		ProbeRuntimeConstruction("runtime object allocation");
		storage.reset(::operator new(sizeof(Runtime)));
	});
	// Deliberately outside a local-stage callback: the constructor enters MPI.
	auto* runtime = ::new (storage.get()) Runtime(std::forward<Arguments>(arguments)...);
	storage.release();
	return std::unique_ptr<Runtime>(runtime);
}

// Observe only handles created successfully on all ranks. Production makes
// no extra PETSc call here. Tests may retain references to prove that failed
// construction releases the runtime's ownership before a retry.
inline void ObserveConstructedObject(MPI_Comm communicator, const char* stage, PetscObject object)
{
#ifdef IGA_RUNTIME_CONSTRUCTION_TESTING
	RuntimeConstructionStage(communicator, stage, [&] {
		const auto& observer = RuntimeConstructionHooksForTesting().object_created;
		if (observer) observer(object);
	});
#else
	(void)communicator;
	(void)stage;
	(void)object;
#endif
}

} // namespace iga

#endif
