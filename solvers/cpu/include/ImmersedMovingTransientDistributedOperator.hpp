#ifndef IGA_IMMERSED_MOVING_TRANSIENT_DISTRIBUTED_OPERATOR_HPP
#define IGA_IMMERSED_MOVING_TRANSIENT_DISTRIBUTED_OPERATOR_HPP

#include "ImmersedTransientDistributedOperator.hpp"
#include "FluidSurfaceTraction.hpp"

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
	// Collective extraction of retained owned-cell coefficients. The existing
	// required-state halo supplies shared control points; no full field gather.
	std::vector<FluidSurfaceElementState> CaptureOwnedPatchState(const MaterialSurfacePatchMap& map)
	{
		std::vector<std::size_t> membership;
		CollectiveLocalStage(comm_, "moving patch capture binding", [&] {
			membership = map.LayoutTrianglesByCanonical(geometry_.Evaluation());
		});
		RequireCollectiveSameText(comm_, "moving patch capture reference", map.ReferenceIdentitySha256());
		std::vector<FluidSurfaceElementState> result;
		Assembly().WithRequiredState([&] {
			for (auto cell : Inputs().OwnedCells()) {
				if (geometry_.Domain().Cells()[cell].classification != CellClassification::Cut) continue;
				const auto& provenance = geometry_.Surface().UsableProvenance(geometry_.Domain(), cell);
				if (!std::any_of(provenance.begin(), provenance.end(), [&](const auto& point) {
					return membership.at(point.canonical_triangle) != std::numeric_limits<std::size_t>::max();
				})) continue;
				const auto element = geometry_.Domain().Background().MaterializeElement(cell);
				FluidSurfaceElementState state; state.cell_id = cell;
				state.nodal_state.resize(element.connectivity.size());
				for (std::size_t node = 0; node < element.connectivity.size(); ++node)
					for (int field = 0; field < 4; ++field) {
						const double value = Assembly().StateAt(4*Layout().LocalNode(element.connectivity[node])+field);
						if (!std::isfinite(value)) throw std::runtime_error("nonfinite moving patch coefficient");
						state.nodal_state[node][field] = value;
					}
				result.push_back(std::move(state));
			}
		});
		return result;
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
	// Collective endpoint diagnostics. Geometry is replicated, but each surface
	// quadrature point contributes only on its volume-cell owner; numerical
	// coefficients are read through the existing required-state halo.
	ImmersedTransientFlowConservationDiagnostics MaterialConservationDiagnostics() const
	{
		const auto fluid=base_.ConservationDiagnostics();
		ImmersedTransientFlowConservationDiagnostics result;
		std::map<int,std::size_t> labels;
		std::vector<double> local,global;
		std::uint64_t local_count=0,global_count=0;
		CollectiveLocalStage(comm_,"moving material conservation storage",[&] {
			result.surface_flow_by_boundary_label_m3_s=fluid.surface_flow_by_boundary_label_m3_s;
			result.endpoint_volume_divergence_m3_s=fluid.volume_divergence_integral_m3_s;
			result.total_surface_outward_flow_m3_s=fluid.total_surface_outward_flow_m3_s;
			result.open_port_outward_flow_m3_s=fluid.open_port_outward_flow_m3_s;
			result.wall_outward_flow_m3_s=fluid.wall_outward_flow_m3_s;
			for(const auto& item:result.surface_flow_by_boundary_label_m3_s) labels.emplace(item.first,labels.size());
			if(labels.size()>=static_cast<std::size_t>(std::numeric_limits<int>::max()))
				throw std::overflow_error("moving material label reduction exceeds MPI count");
			local.resize(labels.size()+1);global.resize(local.size());
		});
		base_.core_->Assembly().WithRequiredState([&] {
			const auto add_raw=[&](double value) {
				if(local_count==std::numeric_limits<std::uint64_t>::max())
					throw std::overflow_error("moving material flux count overflows");
				AddImmersedFlowPortFinite(local.back(),std::abs(value),"moving absolute surface flux");++local_count;
			};
			for(auto cell:Inputs().OwnedCells()) {
				if(geometry_.Domain().Cells()[cell].classification!=CellClassification::Cut) continue;
				const auto element=geometry_.Domain().Background().MaterializeElement(cell);
				const auto& points=geometry_.Surface().UsableRule(geometry_.Domain(),cell).Points();
				const auto& provenance=geometry_.Surface().UsableProvenance(geometry_.Domain(),cell);
				if(points.size()!=provenance.size()) throw std::logic_error("moving material provenance count differs");
				for(std::size_t i=0;i<points.size();++i) {
					const auto& point=points[i];const auto velocity=MaterialVelocity(provenance[i]);
					const auto basis=EvaluateBasis(element,point.parametric[0],point.parametric[1],point.parametric[2],false);
					double material=0,flow=0;
					for(int field=0;field<3;++field)
						AddImmersedFlowPortFinite(material,velocity[field]*point.normal[field]*point.weight,"moving material flux");
					for(std::size_t a=0;a<element.connectivity.size();++a) for(int field=0;field<3;++field)
						AddImmersedFlowPortFinite(flow,base_.core_->Assembly().StateAt(4*Layout().LocalNode(element.connectivity[a])+field)
							*basis.value[a]*point.normal[field]*point.weight,"moving raw fluid flux");
					AddImmersedFlowPortFinite(local.at(labels.at(point.boundary_id)),material,"moving material label flux");
					add_raw(flow);
					if(std::binary_search(Options().wall_labels.begin(),Options().wall_labels.end(),point.boundary_id)) add_raw(material);
				}
			}
		});
		// Each count is bounded before SUM, including ranks with no owned cells.
		int size=0;MPI_Comm_size(comm_,&size);
		CollectiveLocalStage(comm_,"moving material count bound",[&] {
			if(local_count>std::numeric_limits<std::uint64_t>::max()/static_cast<std::uint64_t>(size))
				throw std::overflow_error("moving global flux count may overflow");
		});
		MPI_Allreduce(local.data(),global.data(),static_cast<int>(global.size()),MPI_DOUBLE,MPI_SUM,comm_);
		MPI_Allreduce(&local_count,&global_count,1,MPI_UINT64_T,MPI_SUM,comm_);
		CollectiveLocalStage(comm_,"moving material conservation reduction",[&] {
			for(double value:global) if(!std::isfinite(value)) throw std::runtime_error("moving global material flux is nonfinite");
			result.absolute_surface_flux_sum_m3_s=global.back();result.surface_flux_term_count=global_count;
			for(const auto& item:labels) {
				const double value=global[item.second];result.material_surface_outward_flow_by_boundary_label_m3_s[item.first]=value;
				AddImmersedFlowPortFinite(result.total_material_surface_outward_flow_m3_s,value,"moving total material flux");
				if(std::binary_search(Options().wall_labels.begin(),Options().wall_labels.end(),item.first)) {
					result.material_wall_outward_flow_by_boundary_label_m3_s[item.first]=value;
					AddImmersedFlowPortFinite(result.total_material_wall_outward_flow_m3_s,value,"moving total material wall flux");
				}
			}
			result.wall_relative_leakage_m3_s=result.wall_outward_flow_m3_s-result.total_material_wall_outward_flow_m3_s;
			result.divergence_theorem_defect_m3_s=result.endpoint_volume_divergence_m3_s-result.total_surface_outward_flow_m3_s;
			result.discrete_moving_wall_continuity_defect_m3_s=result.open_port_outward_flow_m3_s+result.total_material_wall_outward_flow_m3_s;
			result.discrete_moving_wall_continuity_normalization_scale_m3_s=std::max({Options().flow_controller_reference_flow_m3_s,
				std::abs(result.open_port_outward_flow_m3_s),std::abs(result.total_material_wall_outward_flow_m3_s)});
			if(!immersed_transient_detail::MovingWallContinuityIdentityReconciles(result.open_port_outward_flow_m3_s,
				result.wall_outward_flow_m3_s,result.total_material_wall_outward_flow_m3_s,result.discrete_moving_wall_continuity_defect_m3_s,
				result.wall_relative_leakage_m3_s,result.total_surface_outward_flow_m3_s,result.absolute_surface_flux_sum_m3_s,result.surface_flux_term_count))
				throw std::logic_error("moving distributed continuity roundoff identity does not reconcile");
			const double scale=std::max(Options().flow_controller_reference_flow_m3_s,std::abs(result.open_port_outward_flow_m3_s));
			result.normalized_wall_leakage=std::abs(result.wall_relative_leakage_m3_s)/scale;
			result.normalized_open_balance=std::abs(result.divergence_theorem_defect_m3_s)/scale;
			result.normalized_discrete_moving_wall_continuity_defect=std::abs(result.discrete_moving_wall_continuity_defect_m3_s)
				/result.discrete_moving_wall_continuity_normalization_scale_m3_s;
			for(double value:{result.divergence_theorem_defect_m3_s,result.discrete_moving_wall_continuity_defect_m3_s,
				result.normalized_wall_leakage,result.normalized_open_balance,result.normalized_discrete_moving_wall_continuity_defect})
				if(!std::isfinite(value)) throw std::runtime_error("moving distributed conservation quotient is nonfinite");
		});
		return result;
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
