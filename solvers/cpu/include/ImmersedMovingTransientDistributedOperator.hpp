#ifndef IGA_IMMERSED_MOVING_TRANSIENT_DISTRIBUTED_OPERATOR_HPP
#define IGA_IMMERSED_MOVING_TRANSIENT_DISTRIBUTED_OPERATOR_HPP

#include "ImmersedTransientDistributedOperator.hpp"

namespace iga {
// One immutable moving-geometry endpoint. Shared sparse topology/port/ghost/
// gauge machinery remains identical to the fixed operator. This class owns
// neither an accepted clock nor an old epoch; its parent runtime owns both.
class ImmersedMovingTransientDistributedOperator
{
public:
	ImmersedMovingTransientDistributedOperator(MPI_Comm comm,const MovingCutGeometry& geometry,
		const ImmersedTransientFlowOptions& options,ImmersedWorkPartition partition=ImmersedWorkPartition::CellCount)
		: comm_(comm),geometry_(geometry),base_(comm,geometry.Domain(),geometry.Volume(),geometry.Surface(),geometry.Ghost(),
			geometry.GeometryIdentitySha256(),options,partition)
	{
		std::string identity;
		CollectiveLocalStage(comm_,"moving operator material binding",[&] {
			geometry_.Evaluation().Validate();identity=geometry_.Evaluation().ContentIdentitySha256();
			for(auto cell:base_.inputs_->OwnedCells()) {
				const auto& rule=geometry_.Surface().UsableRule(geometry_.Domain(),cell);
				const auto& provenance=geometry_.Surface().UsableProvenance(geometry_.Domain(),cell);
				if(rule.Points().size()!=provenance.size())throw std::invalid_argument("moving operator provenance count differs");
				for(const auto& point:provenance)for(double value:MaterialVelocity(point))
					if(!std::isfinite(value))throw std::invalid_argument("moving operator material velocity is nonfinite");
			}
		});
		RequireCollectiveSameText(comm_,"moving operator material agreement",identity);
	}
	ImmersedMovingTransientDistributedOperator(const ImmersedMovingTransientDistributedOperator&)=delete;
	ImmersedMovingTransientDistributedOperator& operator=(const ImmersedMovingTransientDistributedOperator&)=delete;
	const ImmersedActiveLayout& Layout() const noexcept
	{
		return base_.Layout();
	}
	const ImmersedTransientFlowOptions& Options() const noexcept
	{
		return base_.Options();
	}
	const ImmersedStaticFlowDiagnostics& Diagnostics() const noexcept
	{
		return base_.Diagnostics();
	}
	ImmersedDistributedAssembly& Assembly() noexcept
	{
		return base_.Assembly();
	}
	const ImmersedDistributedAssembly& Assembly() const noexcept
	{
		return base_.Assembly();
	}
	const ImmersedDistributedTransientVolume& Inputs() const noexcept
	{
		return base_.Inputs();
	}
	void SetPortControlValue(const std::string& id,double value)
	{
		base_.SetPortControlValue(id,value);
	}
	void Freeze(DistributedImmersedVelocityExtension& extension,Vec committed,const ImmersedActiveLayout& source_layout,
		double source_time,std::uint64_t source_index,double target_time,std::uint64_t target_index,double dt)
	{
		CollectiveLocalStage(comm_,"moving operator epoch preflight",[&] {
			base_.RequireOpen();base_.core_->Topology().ValidateAllFlowCompatibility();
			const auto& material=geometry_.Evaluation();
			if(source_time!=material.StepStartS()||target_time!=material.StepEndS()
				||target_time!=material.EvaluatedTimeS()||dt!=material.DtS())
				throw std::invalid_argument("moving operator trial differs from material epoch");
		});
		auto parameters=base_.options_.parameters;parameters.dt=dt;
		base_.inputs_->Freeze(committed,source_layout,source_time,source_index,target_time,target_index,parameters,base_.options_.body_force,&extension);
		base_.options_.parameters.dt=dt;
	}
	void Assemble()
	{
		CollectiveLocalStage(comm_,"moving operator assembly preflight",[&] {
			base_.RequireOpen();if(!base_.inputs_->History().Active())throw std::logic_error("moving operator inputs are not frozen");
		});
		const ImmersedMaterialWallVelocityEvaluator velocity=[&](const SurfaceQuadraturePoint&,const ImmersedSurfaceQuadraturePointProvenance& point) {
			return MaterialVelocity(point);
		};
		base_.core_->AssembleWithVolume([&](std::uint64_t cell,const auto& nodal) { return base_.BuildCellWithVelocity(cell,nodal,velocity); });
	}
	// Eulerian field flux/divergence diagnostics only. The parent moving
	// runtime additionally needs material flux and geometric conservation.
	ImmersedStaticFlowConservationDiagnostics ConservationDiagnostics() const
	{
		return base_.ConservationDiagnostics();
	}
	void ReleaseTrial() noexcept
	{
		base_.ReleaseTrial();
	}
	void Close()
	{
		base_.Close();
	}
private:
	std::array<double,3> MaterialVelocity(const ImmersedSurfaceQuadraturePointProvenance& point) const
	{
		return geometry_.Evaluation().WallVelocity(point.canonical_triangle,point.canonical_barycentric);
	}
	MPI_Comm comm_;const MovingCutGeometry& geometry_;ImmersedTransientDistributedOperator base_;
};
} // namespace iga
#endif
