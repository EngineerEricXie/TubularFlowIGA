#ifndef IGA_DOMAIN_GROUP_RUNTIME_PROXY_HPP
#define IGA_DOMAIN_GROUP_RUNTIME_PROXY_HPP

#include "CollectiveFailure.hpp"
#include "CoupledDomainRuntime.hpp"

#include <mpi.h>

#include <cstdint>
#include <cstring>
#include <exception>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace iga {

enum class DomainGroupSolveKind { Trial, Hydraulic, Transport };

namespace domain_group_detail {

template<class Value>
void Append(std::string& bytes, const Value& value)
{
	static_assert(std::is_trivially_copyable<Value>::value, "binary proxy values must be trivial");
	bytes.append(reinterpret_cast<const char*>(&value), sizeof(Value));
}

inline void Append(std::string& bytes, const std::string& value)
{
	Append(bytes, static_cast<std::uint64_t>(value.size()));
	bytes.append(value);
}

template<class Value>
Value Read(std::string_view bytes, std::size_t& offset)
{
	static_assert(std::is_trivially_copyable<Value>::value, "binary proxy values must be trivial");
	if (offset > bytes.size() || bytes.size()-offset < sizeof(Value))
		throw std::runtime_error("truncated domain-group message");
	Value value;
	std::memcpy(&value, bytes.data()+offset, sizeof(Value));
	offset += sizeof(Value);
	return value;
}

inline std::string ReadString(std::string_view bytes, std::size_t& offset)
{
	const auto size = Read<std::uint64_t>(bytes, offset);
	if (size > bytes.size()-offset) throw std::runtime_error("truncated domain-group string");
	std::string result(bytes.substr(offset, size));
	offset += static_cast<std::size_t>(size);
	return result;
}

inline void AppendOptional(std::string& bytes, const std::optional<double>& value)
{
	Append(bytes, static_cast<std::uint8_t>(bool(value)));
	if (value) Append(bytes, *value);
}

inline std::optional<double> ReadOptional(std::string_view bytes, std::size_t& offset)
{
	const auto present = Read<std::uint8_t>(bytes, offset);
	if (present > 1) throw std::runtime_error("invalid domain-group optional flag");
	return present ? std::optional<double>(Read<double>(bytes, offset)) : std::nullopt;
}

inline void AppendMap(std::string& bytes, const std::map<std::string, double>& values)
{
	Append(bytes, static_cast<std::uint64_t>(values.size()));
	for (const auto& value : values) {
		Append(bytes, value.first);
		Append(bytes, value.second);
	}
}

inline std::map<std::string, double> ReadMap(std::string_view bytes, std::size_t& offset)
{
	const auto size = Read<std::uint64_t>(bytes, offset);
	if (size > bytes.size()) throw std::runtime_error("invalid domain-group map size");
	std::map<std::string, double> result;
	for (std::uint64_t i = 0; i < size; ++i) {
		const auto key = ReadString(bytes, offset);
		const auto value = Read<double>(bytes, offset);
		if (!result.emplace(key, value).second)
			throw std::runtime_error("duplicate domain-group map key");
	}
	return result;
}

inline std::string Serialize(const PortState& value)
{
	std::string bytes;
	Append(bytes, value.time_s);
	AppendOptional(bytes, value.area_m2);
	AppendOptional(bytes, value.outward_flow_m3_s);
	AppendOptional(bytes, value.mean_pressure_pa);
	AppendOptional(bytes, value.mean_normal_traction_pa);
	AppendOptional(bytes, value.total_pressure_pa);
	AppendMap(bytes, value.concentration);
	AppendMap(bytes, value.outward_species_flux);
	return bytes;
}

inline PortState ParsePortState(std::string_view bytes)
{
	std::size_t offset = 0;
	PortState result;
	result.time_s = Read<double>(bytes, offset);
	result.area_m2 = ReadOptional(bytes, offset);
	result.outward_flow_m3_s = ReadOptional(bytes, offset);
	result.mean_pressure_pa = ReadOptional(bytes, offset);
	result.mean_normal_traction_pa = ReadOptional(bytes, offset);
	result.total_pressure_pa = ReadOptional(bytes, offset);
	result.concentration = ReadMap(bytes, offset);
	result.outward_species_flux = ReadMap(bytes, offset);
	if (offset != bytes.size()) throw std::runtime_error("domain-group port message has trailing bytes");
	ValidatePortState(result);
	return result;
}

inline std::string Serialize(const std::map<std::string, SpeciesStepAccounting>& values)
{
	std::string bytes;
	Append(bytes, static_cast<std::uint64_t>(values.size()));
	for (const auto& value : values) {
		Append(bytes, value.first);
		Append(bytes, value.second.initial_mass);
		Append(bytes, value.second.final_mass);
		Append(bytes, value.second.source_amount);
		AppendMap(bytes, value.second.outward_port_amount);
		Append(bytes, value.second.residual);
	}
	return bytes;
}

inline std::map<std::string, SpeciesStepAccounting> ParseAccounting(std::string_view bytes)
{
	std::size_t offset = 0;
	const auto size = Read<std::uint64_t>(bytes, offset);
	if (size > bytes.size()) throw std::runtime_error("invalid domain-group accounting size");
	std::map<std::string, SpeciesStepAccounting> result;
	for (std::uint64_t i = 0; i < size; ++i) {
		const auto key = ReadString(bytes, offset);
		SpeciesStepAccounting value;
		value.initial_mass = Read<double>(bytes, offset);
		value.final_mass = Read<double>(bytes, offset);
		value.source_amount = Read<double>(bytes, offset);
		value.outward_port_amount = ReadMap(bytes, offset);
		value.residual = Read<double>(bytes, offset);
		if (!result.emplace(key, std::move(value)).second)
			throw std::runtime_error("duplicate domain-group accounting key");
	}
	if (offset != bytes.size()) throw std::runtime_error("domain-group accounting message has trailing bytes");
	return result;
}

} // namespace domain_group_detail

class DomainGroupRuntimeProxy : public CoupledDomainRuntime {
public:
	DomainGroupRuntimeProxy(MPI_Comm world, int owner_world_rank, std::string domain_id,
		DomainKind kind, std::vector<CouplingPort> ports,
		std::unique_ptr<CoupledDomainRuntime> local_runtime)
		: world_(world), owner_world_rank_(owner_world_rank), domain_id_(std::move(domain_id)),
		  kind_(kind), ports_(std::move(ports)), local_(std::move(local_runtime))
	{
		int rank = 0, ranks = 1;
		MPI_Comm_rank(world_, &rank);
		MPI_Comm_size(world_, &ranks);
		CollectiveLocalStage(world_, "domain-group proxy construction", [&] {
			if (owner_world_rank_ < 0 || owner_world_rank_ >= ranks)
				throw std::runtime_error("domain-group owner rank is outside the allocation");
			if (rank == owner_world_rank_ && !local_)
				throw std::runtime_error("domain-group owner rank requires a local runtime");
			if (local_ && (local_->DomainId() != domain_id_ || local_->Kind() != kind_))
				throw std::runtime_error("domain-group proxy runtime metadata differs");
		});
	}

	const std::string& DomainId() const noexcept override { return domain_id_; }
	DomainKind Kind() const noexcept override { return kind_; }
	const std::vector<CouplingPort>& Ports() const noexcept override { return ports_; }

	void BeginStep(const DomainStepContext& step) override
	{
		Execute("domain-group begin", [&](auto& runtime) { runtime.BeginStep(step); });
	}
	void SetPortInput(const std::string& port_id, const PortBoundaryData& input) override
	{
		Execute("domain-group port input", [&](auto& runtime) { runtime.SetPortInput(port_id, input); });
	}
	void SolveTrial() override
	{
		Execute("domain-group solve", [&](auto& runtime) { runtime.SolveTrial(); });
	}
	PortState GetPortState(const std::string& port_id) const override
	{
		PortState result;
		Execute("domain-group port state", [&](auto& runtime) { result = runtime.GetPortState(port_id); });
		return Broadcast(result, "domain-group port state exchange");
	}
	void RollbackTrial() override
	{
		Execute("domain-group rollback", [&](auto& runtime) { runtime.RollbackTrial(); });
	}
	void AbortStep() override
	{
		Execute("domain-group abort", [&](auto& runtime) { runtime.AbortStep(); });
	}
	void PrepareCommitStep() override
	{
		Execute("domain-group prepare commit", [&](auto& runtime) { runtime.PrepareCommitStep(); });
	}
	void FinalizeCommitStep() noexcept override
	{
		if (local_) local_->FinalizeCommitStep();
	}

protected:
	template<class Work>
	void Execute(const char* stage, Work&& work) const
	{
		std::exception_ptr error;
		if (local_) {
			try { std::forward<Work>(work)(*local_); }
			catch (...) { error = std::current_exception(); }
		}
		CollectiveLocalStage(world_, stage, [&] { if (error) std::rethrow_exception(error); });
	}

	PortState Broadcast(const PortState& value, const char* stage) const
	{
		std::string bytes;
		CollectiveLocalStage(world_, stage, [&] {
			int rank = 0;
			MPI_Comm_rank(world_, &rank);
			if (rank == owner_world_rank_) bytes = domain_group_detail::Serialize(value);
		});
		BroadcastBytes(bytes, stage);
		PortState result;
		CollectiveLocalStage(world_, stage, [&] { result = domain_group_detail::ParsePortState(bytes); });
		return result;
	}

	std::map<std::string, SpeciesStepAccounting> Broadcast(
		const std::map<std::string, SpeciesStepAccounting>& value, const char* stage) const
	{
		std::string bytes;
		CollectiveLocalStage(world_, stage, [&] {
			int rank = 0;
			MPI_Comm_rank(world_, &rank);
			if (rank == owner_world_rank_) bytes = domain_group_detail::Serialize(value);
		});
		BroadcastBytes(bytes, stage);
		std::map<std::string, SpeciesStepAccounting> result;
		CollectiveLocalStage(world_, stage, [&] { result = domain_group_detail::ParseAccounting(bytes); });
		return result;
	}

	CoupledDomainRuntime* Local() const noexcept { return local_.get(); }
	friend void SolveDomainGroupBatch(MPI_Comm, const std::vector<CoupledDomainRuntime*>&,
		DomainGroupSolveKind);

private:
	void BroadcastBytes(std::string& bytes, const char* stage) const
	{
		int rank = 0;
		MPI_Comm_rank(world_, &rank);
		std::uint64_t size = rank == owner_world_rank_ ? bytes.size() : 0;
		MPI_Bcast(&size, 1, MPI_UINT64_T, owner_world_rank_, world_);
		CollectiveLocalStage(world_, stage, [&] {
			if (size > 64ULL*1024ULL*1024ULL)
				throw std::runtime_error("domain-group message exceeds 64 MiB");
			bytes.resize(static_cast<std::size_t>(size));
		});
		MPI_Bcast(bytes.data(), static_cast<int>(size), MPI_CHAR, owner_world_rank_, world_);
	}

	MPI_Comm world_;
	int owner_world_rank_ = 0;
	std::string domain_id_;
	DomainKind kind_ = DomainKind::OneDFlow;
	std::vector<CouplingPort> ports_;
	std::unique_ptr<CoupledDomainRuntime> local_;
};

inline void SolveDomainGroupBatch(MPI_Comm world,
	const std::vector<CoupledDomainRuntime*>& runtimes, DomainGroupSolveKind kind)
{
	std::exception_ptr error;
	try {
		CoupledDomainRuntime* local = nullptr;
		for (auto* runtime : runtimes) {
			auto* proxy = dynamic_cast<DomainGroupRuntimeProxy*>(runtime);
			if (!proxy) throw std::runtime_error("domain-group batch contains a non-proxy runtime");
			if (!proxy->Local()) continue;
			if (local) throw std::runtime_error("domain-group batch schedules two domains on one rank group");
			local = proxy->Local();
		}
		if (local) {
			if (kind == DomainGroupSolveKind::Trial) local->SolveTrial();
			else {
				auto* staged = dynamic_cast<StagedFlowTransportDomainRuntime*>(local);
				if (!staged) throw std::runtime_error("domain-group staged batch contains a non-staged runtime");
				if (kind == DomainGroupSolveKind::Hydraulic) staged->SolveHydraulicTrial();
				else staged->SolveTransportTrial();
			}
		}
	} catch (...) { error = std::current_exception(); }
	CollectiveLocalStage(world, "domain-group solve batch", [&] {
		if (error) std::rethrow_exception(error);
	});
}

class StagedDomainGroupRuntimeProxy final : public DomainGroupRuntimeProxy,
	public StagedFlowTransportDomainRuntime {
public:
	using DomainGroupRuntimeProxy::DomainGroupRuntimeProxy;

	void SolveHydraulicTrial() override
	{
		Execute("domain-group hydraulic solve", [&](auto& runtime) { Staged(runtime).SolveHydraulicTrial(); });
	}
	PortState GetHydraulicPortState(const std::string& port_id) const override
	{
		PortState result;
		Execute("domain-group hydraulic state", [&](auto& runtime) {
			result = Staged(runtime).GetHydraulicPortState(port_id);
		});
		return Broadcast(result, "domain-group hydraulic state exchange");
	}
	void RollbackHydraulicTrial() override
	{
		Execute("domain-group hydraulic rollback", [&](auto& runtime) { Staged(runtime).RollbackHydraulicTrial(); });
	}
	void SetTransportConcentration(const std::string& port_id, double time_s,
		const std::map<std::string, double>& concentration) override
	{
		Execute("domain-group transport input", [&](auto& runtime) {
			Staged(runtime).SetTransportConcentration(port_id, time_s, concentration);
		});
	}
	void SolveTransportTrial() override
	{
		Execute("domain-group transport solve", [&](auto& runtime) { Staged(runtime).SolveTransportTrial(); });
	}
	PortState GetTransportPortState(const std::string& port_id) const override
	{
		PortState result;
		Execute("domain-group transport state", [&](auto& runtime) {
			result = Staged(runtime).GetTransportPortState(port_id);
		});
		return Broadcast(result, "domain-group transport state exchange");
	}
	void RollbackTransportTrial() override
	{
		Execute("domain-group transport rollback", [&](auto& runtime) { Staged(runtime).RollbackTransportTrial(); });
	}
	std::map<std::string, SpeciesStepAccounting> GetSpeciesStepAccounting() const override
	{
		std::map<std::string, SpeciesStepAccounting> result;
		Execute("domain-group species accounting", [&](auto& runtime) {
			result = Staged(runtime).GetSpeciesStepAccounting();
		});
		return Broadcast(result, "domain-group species accounting exchange");
	}

private:
	static StagedFlowTransportDomainRuntime& Staged(CoupledDomainRuntime& runtime)
	{
		auto* staged = dynamic_cast<StagedFlowTransportDomainRuntime*>(&runtime);
		if (!staged) throw std::runtime_error("domain-group runtime does not support staged transport");
		return *staged;
	}
};

} // namespace iga

#endif
