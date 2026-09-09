#ifndef IGA_IMMERSED_DISTRIBUTED_TRANSIENT_VOLUME_HPP
#define IGA_IMMERSED_DISTRIBUTED_TRANSIENT_VOLUME_HPP

#include "ImmersedDistributedVelocityHistory.hpp"
#include "RuntimeConstruction.hpp"
#include <map>

namespace iga {

// Owned fixed-geometry backward-Euler volume integration. This component
// freezes history and body force together; it does not own Newton, wall/port
// terms, or the accepted clock. Geometry, layout and communicator are borrowed
// and immutable. Construction, Freeze and Close are collective. BuildVolume
// is local and can be called from an ImmersedDistributedAssembly callback.
class ImmersedDistributedTransientVolume {
public:
	ImmersedDistributedTransientVolume(MPI_Comm communicator,
		const CartesianDomainClassification& domain, const CutCellVolumeQuadratureCatalog& volume,
		const ImmersedActiveLayout& layout, Vec state_template, const std::vector<std::uint64_t>& owned_cells)
		: communicator_(communicator), domain_(domain), volume_(volume), layout_(layout)
	{
		std::vector<int> local_owners,global_owners,active;
		std::vector<std::int32_t> required;
		std::string signature;
		CollectiveLocalStage(communicator_,"transient volume topology",[&] {
			if (domain.Cells().size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
				throw std::overflow_error("transient volume ownership exceeds MPI count capacity");
			const auto expected = ImmersedActiveLayout::Build(domain,volume,layout.GeometryIdentity(),layout.PortIds(),layout.HasGaugeRow());
			if (expected.HashSha256() != layout.HashSha256()) throw std::invalid_argument("transient volume active layout differs");
			local_owners.assign(domain.Cells().size(),0); global_owners.resize(local_owners.size()); active.resize(local_owners.size());
			for (std::size_t c = 0; c < active.size(); ++c) {
				const auto classification = domain.Cells()[c].classification;
				active[c] = classification == CellClassification::Inside || (classification == CellClassification::Cut
					&& volume.Cell(c).diagnostics.estimated_reference_volume > 0);
			}
			cells_ = owned_cells; std::sort(cells_.begin(),cells_.end());
			for (auto c : cells_) {
				if (c >= active.size() || !active[c] || local_owners[c]++) throw std::invalid_argument("invalid or duplicate transient volume owner cell");
				if (volume.StorageMode() == CutCellVolumeQuadratureStorageMode::Expanded) volume.ValidateUsableRule(domain,c);
				else volume.ValidateUsableCompactRule(domain,c);
				if (!PointCount(c)) throw std::invalid_argument("active transient volume has no quadrature points");
				const auto element = domain.Background().MaterializeElement(c);
				required.insert(required.end(),element.connectivity.begin(),element.connectivity.end());
			}
			signature = Signature();
		});
		RequireCollectiveSameText(communicator_,"transient volume configuration agreement",signature);
		MPI_Allreduce(local_owners.data(),global_owners.data(),static_cast<int>(local_owners.size()),MPI_INT,MPI_SUM,communicator_);
		CollectiveLocalStage(communicator_,"transient volume unique ownership",[&] {
			if (global_owners != active) throw std::invalid_argument("transient volume cells require exactly one owner");
		});
		history_ = AllocateCollectiveRuntime<ImmersedDistributedVelocityHistory>(communicator_,communicator_,layout,state_template,required);
	}
	ImmersedDistributedTransientVolume(const ImmersedDistributedTransientVolume&) = delete;
	ImmersedDistributedTransientVolume& operator=(const ImmersedDistributedTransientVolume&) = delete;

	// Callbacks must be local and represent the same physical field on all
	// members. Only the cell owner evaluates them, once per canonical point.
	// A failed callback discards both force candidate and newly frozen history.
	void Freeze(Vec committed, const ImmersedActiveLayout& source_layout,
		double source_time, std::uint64_t source_index, double target_time, std::uint64_t target_index,
		const NavierStokesParameters& parameters, const NavierStokesBodyForceEvaluator& body_force)
	{
		std::string signature;
		CollectiveLocalStage(communicator_,"transient volume freeze preflight",[&] {
			RequireOpen();
			if (history_->Active()) throw std::logic_error("transient volume inputs are already frozen");
			ValidateTransientNavierStokesPreflight(parameters);
			if (!body_force || !std::isfinite(parameters.dynamic_viscosity) || !(parameters.dynamic_viscosity > 0))
				throw std::invalid_argument("invalid transient volume viscosity or force evaluator");
			std::ostringstream text; text.exceptions(std::ios::badbit | std::ios::failbit);
			text << std::setprecision(std::numeric_limits<double>::max_digits10)
				<< parameters.density << ':' << parameters.dynamic_viscosity << ':' << parameters.dt;
			signature = text.str();
		});
		RequireCollectiveSameText(communicator_,"transient volume parameter agreement",signature);
		history_->Freeze(committed,source_layout,source_time,source_index,target_time,target_index,parameters.dt);
		std::map<std::uint64_t,std::vector<std::array<double,3>>> candidate;
		try {
			CollectiveLocalStage(communicator_,"transient volume force freeze",[&] {
				for (auto c : cells_) {
					const auto element = domain_.Background().MaterializeElement(c);
					auto& forces = candidate[c]; forces.reserve(PointCount(c));
					VolumePoints(c,[&](const VolumeQuadraturePoint& point) {
						const auto force = body_force(EvaluateElementGeometry(element,point.parametric).physical);
						for (double value : force) if (!std::isfinite(value)) throw std::invalid_argument("nonfinite transient volume force");
						forces.push_back(force);
					});
					if (forces.size() != PointCount(c)) throw std::logic_error("transient volume frozen force count differs");
				}
			});
		} catch (...) { history_->ReleaseTrial(); throw; }
		forces_.swap(candidate); parameters_ = parameters;
	}

	NavierStokesSystem BuildVolume(std::uint64_t cell, const std::vector<std::array<double,4>>& nodal_state) const
	{
		RequireOpen();
		if (!history_->Active()) throw std::logic_error("transient volume inputs are not frozen");
		const auto found = forces_.find(cell);
		if (found == forces_.end()) throw std::out_of_range("transient volume cell is not locally owned");
		const auto element = domain_.Background().MaterializeElement(cell);
		const auto velocity = history_->Localize(element,history_->TargetTimeS());
		std::vector<std::array<double,4>> old(velocity.size());
		for (std::size_t i = 0; i < old.size(); ++i)
			for (int c = 0; c < 3; ++c) old[i][c] = velocity[i][c];
		std::size_t index = 0;
		const auto force = [&](const std::array<double,3>&) {
			if (index >= found->second.size()) throw std::logic_error("transient volume force replay exceeded frozen points");
			return found->second[index++];
		};
		auto system = BuildNavierStokesElementFromPoints(element,nodal_state,old,parameters_,
			[&](const auto& consume) { VolumePoints(cell,consume); },force,NavierStokesResolvedMixedForm::Conservative);
		if (index != found->second.size()) throw std::logic_error("transient volume force replay missed frozen points");
		return system;
	}
	const ImmersedDistributedVelocityHistory& History() const noexcept { return *history_; }
	const std::vector<std::uint64_t>& OwnedCells() const noexcept { return cells_; }
	std::size_t FrozenPointCount() const noexcept
	{
		std::size_t count = 0;
		for (const auto& cell : forces_) count += cell.second.size();
		return count;
	}
	void ReleaseTrial() noexcept { history_->ReleaseTrial(); forces_.clear(); }
	void Close() { closed_ = true; ReleaseTrial(); history_->Close(); }
private:
	std::size_t PointCount(std::uint64_t cell) const
	{
		return volume_.StorageMode() == CutCellVolumeQuadratureStorageMode::Expanded ? volume_.Cell(cell).rule.Points().size()
			: CompactCutCellVolumeLogicalPointCount(volume_.Cell(cell).compact_rule);
	}
	template<class Function> void VolumePoints(std::uint64_t cell, Function&& function) const
	{
		if (volume_.StorageMode() == CutCellVolumeQuadratureStorageMode::Expanded)
			for (const auto& point : volume_.UsableRule(domain_,cell).Points()) function(point);
		else ForEachVolumePoint(volume_.UsableCompactRule(domain_,cell),function);
	}
	std::string Signature() const
	{
		std::ostringstream text; text.exceptions(std::ios::badbit | std::ios::failbit);
		text << std::setprecision(std::numeric_limits<double>::max_digits10)
			<< layout_.HashSha256() << ':' << domain_.SurfaceCanonicalHash() << ':' << static_cast<int>(volume_.StorageMode());
		const auto& grid = domain_.Background().Spec();
		for (int d = 0; d < 3; ++d) text << ':' << grid.lower_m[d] << ':' << grid.upper_m[d] << ':' << grid.cells[d];
		const auto& q = volume_.Options();
		text << ':' << q.max_depth << ':' << q.max_nodes << ':' << q.max_leaves << ':' << q.max_points
			<< ':' << q.max_records << ':' << q.max_retained_bytes << ':' << q.max_logical_points << ':' << q.empty_rule_rescue_max_depth;
		return text.str();
	}
	void RequireOpen() const { if (closed_) throw std::logic_error("transient volume is closed"); }
	MPI_Comm communicator_;
	const CartesianDomainClassification& domain_;
	const CutCellVolumeQuadratureCatalog& volume_;
	const ImmersedActiveLayout& layout_;
	std::vector<std::uint64_t> cells_;
	std::unique_ptr<ImmersedDistributedVelocityHistory> history_;
	std::map<std::uint64_t,std::vector<std::array<double,3>>> forces_;
	NavierStokesParameters parameters_;
	bool closed_ = false;
};

} // namespace iga
#endif
