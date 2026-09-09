#include "CollectiveDomainRuntimeRegistry.hpp"
#include <array>
#include <petscvec.h>

#define main RegistryFixtureMain
#include "test_pressure_flow_executor.cpp"
#undef main

namespace {

void Check(bool value, const char* message)
{
	if (!value) throw std::runtime_error(message);
}

class TrackedRuntime : public FakeRuntime {
public:
	TrackedRuntime(const iga::DomainNode& node, int& destroyed)
		: FakeRuntime(node.id, node.kind, node.ports), id(node.id), kind(node.kind),
		  ports(node.ports), destroyed_(destroyed) {}
	~TrackedRuntime() override { if (vector) VecDestroy(&vector); ++destroyed_; }
	const std::string& DomainId() const noexcept override { return id; }
	iga::DomainKind Kind() const noexcept override { return kind; }
	const std::vector<iga::CouplingPort>& Ports() const noexcept override { return ports; }
	std::string id;
	iga::DomainKind kind;
	std::vector<iga::CouplingPort> ports;
	Vec vector = nullptr;
private:
	int& destroyed_;
};

void RunCases(MPI_Comm comm, int seed)
{
	int rank = 0, ranks = 0;
	MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &ranks);
	const auto graph = Chain();
	const std::array<std::string, 3> ids{{"up", "mid", "down"}};
	const std::vector<std::string> faults{"kind", "domain", "ports", "count", "null",
		"registry synchronization ready", "runtime object allocation",
		"runtime object storage ready", "registry index ready"};
	for (const auto& fault : faults) {
		int destroyed = 0;
		std::vector<std::unique_ptr<iga::CoupledDomainRuntime>> owners;
		std::array<TrackedRuntime*, 3> raw{};
		std::array<Vec, 3> retained{};
		iga::CollectiveLocalStage(comm, "registry test owners", [&] {
			for (std::size_t i = 0; i < ids.size(); ++i) {
				auto runtime = std::make_unique<TrackedRuntime>(graph.Domain(ids[i]), destroyed);
				raw[i] = runtime.get(); owners.push_back(std::move(runtime));
			}
		});
		for (std::size_t i = 0; i < ids.size(); ++i) {
			iga::RequireCollectivePetscSuccess(comm, "registry test vector", VecCreateMPI(comm, 1, ranks, &raw[i]->vector));
			retained[i] = raw[i]->vector;
			iga::RequireCollectivePetscSuccess(comm, "registry test retained vector",
				PetscObjectReference(reinterpret_cast<PetscObject>(retained[i])));
		}
		std::unique_ptr<iga::CoupledDomainRuntime> spare;
		iga::CollectiveLocalStage(comm, "registry test mutation", [&] {
			if (rank != ranks-1) return;
			if (fault == "kind") raw[0]->kind = iga::DomainKind::ZeroDFlow;
			if (fault == "domain") raw[0]->id = "missing";
			if (fault == "ports") raw[0]->ports.front().locator = "invalid";
			if (fault == "count" || fault == "null") spare = std::move(owners.back());
			if (fault == "count") owners.pop_back();
		});
		const auto count = owners.size();
		std::array<iga::CoupledDomainRuntime*, 3> prior{};
		for (std::size_t i = 0; i < count; ++i) prior[i] = owners[i].get();
		bool injected = false;
		iga::RuntimeConstructionHooksForTesting().after_local_stage = [&](const char* stage) {
			if (rank == ranks-1 && !injected && fault == stage) {
				injected = true; throw std::bad_alloc();
			}
		};
		std::string diagnostic;
		try { auto registry = iga::CreateCollectiveDomainRuntimeRegistry(comm, graph, owners); }
		catch (const std::exception& error) { diagnostic = error.what(); }
		iga::RuntimeConstructionHooksForTesting() = {};
		iga::CollectiveLocalStage(comm, "registry retained ownership", [&] {
			Check(!diagnostic.empty(), "invalid registry was accepted");
			Check(destroyed == 0 && owners.size() == count, "failed registry consumed runtime ownership");
			for (std::size_t i = 0; i < count; ++i)
				Check(owners[i].get() == prior[i], "failed registry changed a runtime pointer");
		});
		iga::RequireCollectiveSameText(comm, "registry common diagnostic", diagnostic);
		iga::CollectiveLocalStage(comm, "registry retry input", [&] {
			for (std::size_t i = 0; i < ids.size(); ++i) {
				raw[i]->id = ids[i]; raw[i]->kind = graph.Domain(ids[i]).kind;
				raw[i]->ports = graph.Domain(ids[i]).ports;
			}
			if (spare) {
				if (fault == "count") owners.push_back(std::move(spare));
				else owners.back() = std::move(spare);
			}
		});
		auto registry = iga::CreateCollectiveDomainRuntimeRegistry(comm, graph, owners);
		iga::CollectiveLocalStage(comm, "registry transfer verification", [&] {
			Check(owners.empty() && destroyed == 0, "successful registry did not adopt owners");
			for (std::size_t i = 0; i < ids.size(); ++i)
				Check(&registry->Runtime(ids[i]) == raw[i], "registry index changed runtime identity");
		});
		registry.reset();
		iga::CollectiveLocalStage(comm, "registry destruction verification", [&] {
			Check(destroyed == 3, "registry did not destroy each runtime exactly once");
			for (Vec vector : retained) {
				PetscInt references = 0;
				Check(PetscObjectGetReference(reinterpret_cast<PetscObject>(vector), &references) == 0
					&& references == 1, "registry leaked PETSc ownership");
			}
		});
		for (auto& vector : retained)
			iga::RequireCollectivePetscSuccess(comm, "registry retained release", VecDestroy(&vector));
	}
	if (rank == 0) std::cout << "registry ranks=" << ranks << " seed=" << seed
		<< " faults=9 retries=9 vectors_released=27\n";
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	try {
		int rank = 0, ranks = 0;
		MPI_Comm_rank(PETSC_COMM_WORLD, &rank); MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
		Check(ranks == 3, "run with three ranks");
		RunCases(PETSC_COMM_WORLD, 0);
		MPI_Comm subgroup = MPI_COMM_NULL;
		MPI_Comm_split(PETSC_COMM_WORLD, rank == 0 ? 0 : 1, rank, &subgroup);
		RunCases(subgroup, rank == 0 ? 1 : 2);
		MPI_Comm_free(&subgroup);
		PetscFinalize(); return 0;
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n'; MPI_Abort(PETSC_COMM_WORLD, 2); return 2;
	}
}
