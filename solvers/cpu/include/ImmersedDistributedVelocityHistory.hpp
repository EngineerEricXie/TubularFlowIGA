#ifndef IGA_IMMERSED_DISTRIBUTED_VELOCITY_HISTORY_HPP
#define IGA_IMMERSED_DISTRIBUTED_VELOCITY_HISTORY_HPP

#include "CollectiveFailure.hpp"
#include "ExecutionResources.hpp"
#include "ImmersedTransientState.hpp"
#include "NavierStokesElement.hpp"
#include "PetscReadArray.hpp"
#include "RuntimeCleanup.hpp"
#include <petscvec.h>
#include <iomanip>
#include <sstream>

namespace iga {

// Fixed-layout identity history. Only three velocity fields per owned node
// and the requested local halo are stored. Pressure and scalar controller
// rows never enter the history. The immutable layout and communicator must
// outlive this object. Construction, Freeze and Close are collective.
// ReleaseTrial is a local, noexcept publication operation: the runtime calls
// it on every member after a coordinated commit/abort. It does not commit a
// state or advance the runtime clock itself.
class ImmersedDistributedVelocityHistory {
public:
	ImmersedDistributedVelocityHistory(MPI_Comm communicator, const ImmersedActiveLayout& layout,
		Vec state_template, const std::vector<std::int32_t>& required_node_ids)
		: communicator_(communicator), layout_(layout)
	{
		std::vector<PetscInt> extraction_rows, halo_rows;
		std::string signature;
		CollectiveLocalStage(communicator_, "immersed history topology", [&] {
			RequirePetscRealDouble();
			if (!layout.Valid() || layout.NodeIds().empty()) throw std::invalid_argument("invalid immersed history layout");
			InspectSource(state_template, begin_, end_);
			const auto physical = static_cast<PetscInt>(layout.NodeFieldRows());
			const auto first = std::min(begin_,physical), last = std::min(end_,physical);
			if (first%4 || last%4) throw std::invalid_argument("immersed history ownership splits a node");
			owned_nodes_ = (last-first)/4;
			for (PetscInt row = first; row < last; row += 4)
				for (PetscInt field = 0; field < 3; ++field) extraction_rows.push_back(row+field);
			required_ids_ = required_node_ids;
			std::sort(required_ids_.begin(),required_ids_.end());
			required_ids_.erase(std::unique(required_ids_.begin(),required_ids_.end()),required_ids_.end());
			for (const auto id : required_ids_)
				for (PetscInt field = 0; field < 3; ++field)
					halo_rows.push_back(static_cast<PetscInt>(3*layout.LocalNode(id))+field);
			frozen_.resize(halo_rows.size()); candidate_.resize(halo_rows.size());
			signature = layout.HashSha256();
		});
		RequireCollectiveSameText(communicator_, "immersed history layout agreement", signature);
		try {
			Check("history owned velocity create",VecCreateMPI(communicator_,3*owned_nodes_,
				static_cast<PetscInt>(3*layout.NodeIds().size()),&owned_));
			Check("history extraction rows",ISCreateGeneral(PETSC_COMM_SELF,static_cast<PetscInt>(extraction_rows.size()),extraction_rows.data(),PETSC_COPY_VALUES,&source_rows_));
			Check("history extraction offsets",ISCreateStride(PETSC_COMM_SELF,3*owned_nodes_,3*(std::min(begin_,static_cast<PetscInt>(layout.NodeFieldRows()))/4),1,&owned_rows_));
			Check("history extraction create",VecScatterCreate(state_template,source_rows_,owned_,owned_rows_,&extract_));
			Check("history halo create",VecCreateSeq(PETSC_COMM_SELF,static_cast<PetscInt>(halo_rows.size()),&halo_));
			Check("history halo rows",ISCreateGeneral(PETSC_COMM_SELF,static_cast<PetscInt>(halo_rows.size()),halo_rows.data(),PETSC_COPY_VALUES,&halo_rows_));
			Check("history halo offsets",ISCreateStride(PETSC_COMM_SELF,static_cast<PetscInt>(halo_rows.size()),0,1,&halo_offsets_));
			Check("history halo scatter",VecScatterCreate(owned_,halo_rows_,halo_,halo_offsets_,&scatter_));
		} catch (...) { Release(); throw; }
	}
	~ImmersedDistributedVelocityHistory() { Release(); }
	ImmersedDistributedVelocityHistory(const ImmersedDistributedVelocityHistory&) = delete;
	ImmersedDistributedVelocityHistory& operator=(const ImmersedDistributedVelocityHistory&) = delete;

	// All validation and exchanges precede publication. A failed Freeze leaves
	// the object idle and retryable. Later Newton mutations of the source Vec
	// cannot change this snapshot. Source provenance is always Committed;
	// remapping to a different geometry/layout is deliberately rejected here.
	void Freeze(Vec committed, const ImmersedActiveLayout& source_layout, double source_time,
		std::uint64_t source_index, double target_time, std::uint64_t target_index, double dt)
	{
		std::string signature;
		CollectiveLocalStage(communicator_, "immersed history freeze preflight", [&] {
			RequireOpen();
			if (active_) throw std::logic_error("immersed velocity history is already frozen");
			if (source_layout.HashSha256() != layout_.HashSha256()) throw std::invalid_argument("immersed history source layout differs");
			PetscInt first = 0,last = 0; InspectSource(committed,first,last);
			if (first != begin_ || last != end_) throw std::invalid_argument("immersed history source ownership differs");
			if (!std::isfinite(dt) || !(dt > 0) || target_time != CheckedTransientTargetTime(source_time,dt))
				throw std::invalid_argument("immersed history target time differs");
			if (source_index == std::numeric_limits<std::uint64_t>::max() || target_index != source_index+1)
				throw std::invalid_argument("immersed history target index differs");
			std::ostringstream text; text.exceptions(std::ios::badbit | std::ios::failbit);
			text << std::setprecision(std::numeric_limits<double>::max_digits10)
				<< source_time << ':' << source_index << ':' << target_time << ':' << target_index << ':' << dt;
			signature = text.str();
		});
		RequireCollectiveSameText(communicator_, "immersed history time agreement", signature);
		Check("history extraction begin",VecScatterBegin(extract_,committed,owned_,INSERT_VALUES,SCATTER_FORWARD));
		Check("history extraction end",VecScatterEnd(extract_,committed,owned_,INSERT_VALUES,SCATTER_FORWARD));
		CollectiveLocalStage(communicator_, "immersed history owned velocity validation", [&] {
			PetscReadArray view; view.Acquire(owned_);
			for (PetscInt i = 0; i < 3*owned_nodes_; ++i)
				if (!std::isfinite(PetscRealPart(view.Data()[i]))) throw std::invalid_argument("nonfinite immersed committed velocity");
			view.Restore();
		});
		Check("history exchange begin",VecScatterBegin(scatter_,owned_,halo_,INSERT_VALUES,SCATTER_FORWARD));
		Check("history exchange end",VecScatterEnd(scatter_,owned_,halo_,INSERT_VALUES,SCATTER_FORWARD));
		CollectiveLocalStage(communicator_, "immersed history halo validation", [&] {
			PetscReadArray view; view.Acquire(halo_);
			for (std::size_t i = 0; i < candidate_.size(); ++i) {
				candidate_[i] = PetscRealPart(view.Data()[i]);
				if (!std::isfinite(candidate_[i])) throw std::invalid_argument("nonfinite immersed history halo");
			}
			view.Restore();
		});
		frozen_.swap(candidate_); source_time_ = source_time; target_time_ = target_time;
		source_index_ = source_index; target_index_ = target_index; active_ = true;
	}

	std::vector<std::array<double,3>> Localize(const Element& element, double target_time) const
	{
		RequireActive();
		if (target_time != target_time_) throw std::invalid_argument("immersed history assembly time differs");
		if (element.connectivity.empty()) throw std::invalid_argument("immersed history element is empty");
		std::vector<std::array<double,3>> values; values.reserve(element.connectivity.size());
		for (const auto id : element.connectivity) {
			const auto found = std::lower_bound(required_ids_.begin(),required_ids_.end(),id);
			if (found == required_ids_.end() || *found != id) throw std::out_of_range("immersed history node is outside required halo");
			const auto offset = 3*static_cast<std::size_t>(found-required_ids_.begin());
			values.push_back({{frozen_[offset],frozen_[offset+1],frozen_[offset+2]}});
		}
		return values;
	}
	bool Active() const noexcept { return active_; }
	PetscInt OwnedNodes() const noexcept { return owned_nodes_; }
	const std::vector<std::int32_t>& RequiredNodeIds() const noexcept { return required_ids_; }
	double SourceTimeS() const { RequireActive(); return source_time_; }
	double TargetTimeS() const { RequireActive(); return target_time_; }
	std::uint64_t SourceIndex() const { RequireActive(); return source_index_; }
	std::uint64_t TargetIndex() const { RequireActive(); return target_index_; }
	const ImmersedActiveLayout& Layout() const noexcept { return layout_; }
	void ReleaseTrial() noexcept { active_ = false; }
	void Close() { Release(); cleanup_.Check(communicator_,"immersed history close"); }
private:
	void InspectSource(Vec source, PetscInt& first, PetscInt& last) const
	{
		if (!source) throw std::invalid_argument("immersed history source vector is null");
		MPI_Comm source_comm = MPI_COMM_NULL; int relation = MPI_UNEQUAL; PetscInt rows = 0;
		if (PetscObjectGetComm(reinterpret_cast<PetscObject>(source),&source_comm)
			|| MPI_Comm_compare(communicator_,source_comm,&relation) != MPI_SUCCESS
			|| (relation != MPI_IDENT && relation != MPI_CONGRUENT))
			throw std::invalid_argument("immersed history source communicator differs");
		if (VecGetSize(source,&rows) || VecGetOwnershipRange(source,&first,&last)
			|| rows != static_cast<PetscInt>(layout_.Rows())) throw std::invalid_argument("immersed history source rows differ");
	}
	void RequireOpen() const { if (closed_) throw std::logic_error("immersed history is closed"); }
	void RequireActive() const { RequireOpen(); if (!active_) throw std::logic_error("immersed history is not frozen"); }
	void Check(const char* stage, PetscErrorCode code) const { RequireCollectivePetscSuccess(communicator_,stage,code); }
	void Release() noexcept
	{
		if (closed_) return;
		closed_ = true; active_ = false;
		if (scatter_) cleanup_.Observe("history VecScatterDestroy halo",VecScatterDestroy(&scatter_));
		if (extract_) cleanup_.Observe("history VecScatterDestroy extraction",VecScatterDestroy(&extract_));
		if (source_rows_) cleanup_.Observe("history ISDestroy source",ISDestroy(&source_rows_));
		if (owned_rows_) cleanup_.Observe("history ISDestroy owned",ISDestroy(&owned_rows_));
		if (halo_rows_) cleanup_.Observe("history ISDestroy halo",ISDestroy(&halo_rows_));
		if (halo_offsets_) cleanup_.Observe("history ISDestroy offsets",ISDestroy(&halo_offsets_));
		if (halo_) cleanup_.Observe("history VecDestroy halo",VecDestroy(&halo_));
		if (owned_) cleanup_.Observe("history VecDestroy owned",VecDestroy(&owned_));
	}
	MPI_Comm communicator_;
	const ImmersedActiveLayout& layout_;
	PetscInt begin_ = 0,end_ = 0,owned_nodes_ = 0;
	Vec owned_ = nullptr,halo_ = nullptr;
	IS source_rows_ = nullptr,owned_rows_ = nullptr,halo_rows_ = nullptr,halo_offsets_ = nullptr;
	VecScatter extract_ = nullptr,scatter_ = nullptr;
	std::vector<std::int32_t> required_ids_;
	std::vector<double> frozen_,candidate_;
	double source_time_ = 0,target_time_ = 0;
	std::uint64_t source_index_ = 0,target_index_ = 0;
	bool active_ = false,closed_ = false;
	RuntimeCleanupResult cleanup_;
};

inline NavierStokesSystem BuildTransientNavierStokesElement(const Element& element,
	const std::vector<std::array<double,4>>& nodal_state,
	const ImmersedDistributedVelocityHistory& history, double target_time,
	const NavierStokesParameters& parameters, const VolumeQuadratureRule& quadrature,
	const NavierStokesBodyForceEvaluator& body_force,
	NavierStokesResolvedMixedForm form = NavierStokesResolvedMixedForm::Conservative)
{
	ValidateTransientNavierStokesPreflight(parameters);
	if (target_time != CheckedTransientTargetTime(history.SourceTimeS(),parameters.dt))
		throw std::invalid_argument("immersed distributed history time step differs");
	return BuildNavierStokesElementFromLocalVelocityHistory(element,nodal_state,
		history.Localize(element,target_time),parameters,quadrature,body_force,form);
}

} // namespace iga
#endif
