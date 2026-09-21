#ifndef IGA_IMMERSED_TRANSIENT_DISTRIBUTED_OPERATOR_HPP
#define IGA_IMMERSED_TRANSIENT_DISTRIBUTED_OPERATOR_HPP

#include "ImmersedStaticDistributedOperator.hpp"
#include "ImmersedDistributedTransientVolume.hpp"
#include "ImmersedTransientFlowRuntime.hpp"

namespace iga {

// Stationary, immutable catalogs define this fixed-geometry operator for all
// trial times. Material wall velocity is identically zero by this interface;
// moving/evaluated-epoch geometry must use a separate moving runtime. Source
// history is owned, while shared port/gauge/ghost assembly retains its exact
// distributed topology. No serial numerical runtime is constructed here.
class ImmersedTransientDistributedOperator {
public:
	ImmersedTransientDistributedOperator(MPI_Comm communicator,
		const CartesianDomainClassification& domain,const CutCellVolumeQuadratureCatalog& volume,
		const ImmersedSurfaceQuadratureCatalog& surface,const CutCellGhostPenaltyCatalog& ghost,
		const std::string& geometry_identity,const ImmersedTransientFlowOptions& options,
		ImmersedWorkPartition partition = ImmersedWorkPartition::CellCount)
		: communicator_(communicator),domain_(domain),volume_(volume),surface_(surface),ghost_(ghost)
	{
		ImmersedStaticFlowOptions topology_options;
		std::string signature;
		CollectiveLocalStage(communicator_,"transient operator setup",[&] {
			options_ = options; options_.parameters.dt = 0;
			if (!options.pressure_port_backflow_beta.empty())
				throw std::invalid_argument("distributed immersed pressure-port backflow is not implemented");
			if (!std::isfinite(options.wall_inertial_gamma0) || options.wall_inertial_gamma0 < 0
				|| !(options.flow_controller_absolute_tolerance_m3_s > 0))
				throw std::invalid_argument("invalid transient wall impedance or flow tolerance");
			const bool pressure = std::any_of(options.ports.begin(),options.ports.end(),[](const auto& p) { return IsImmersedFlowPressureLike(p.control_mode); });
			if (!options.include_pressure_gauge && !pressure) throw std::invalid_argument("transient pressure gauge is mandatory without a pressure-like port");
			std::vector<std::uint64_t> ids;
			for (const auto& p : options.ports) if (p.control_mode == ImmersedFlowPortControlMode::FlowRate) {
				ValidateImmersedFlowPortDefinition(p); ids.push_back(static_cast<std::uint64_t>(p.boundary_label));
			}
			std::sort(ids.begin(),ids.end());
			layout_ = ImmersedActiveLayout::Build(domain,volume,geometry_identity,ids,!pressure);
			topology_options.parameters = options_.parameters; topology_options.wall_labels = options.wall_labels;
			topology_options.ports = options.ports; topology_options.wall_gamma0 = options.wall_gamma0;
			topology_options.body_force = options.body_force;
			topology_options.solver_options_prefix = options.solver_options_prefix;
			topology_options.nonlinear_maximum_iterations = options.nonlinear_maximum_iterations;
			topology_options.ksp_maximum_iterations = options.ksp_maximum_iterations;
			topology_options.ksp_relative_tolerance = options.ksp_relative_tolerance;
			topology_options.nonlinear_relative_tolerance = options.nonlinear_relative_tolerance;
			topology_options.nonlinear_absolute_tolerance = options.nonlinear_absolute_tolerance;
			topology_options.flow_controller_relative_tolerance = options.flow_controller_relative_tolerance;
			topology_options.flow_controller_absolute_tolerance_m3_s = options.flow_controller_absolute_tolerance_m3_s;
			topology_options.flow_controller_reference_flow_m3_s = options.flow_controller_reference_flow_m3_s;
			topology_options.minimum_damping = options.minimum_damping; topology_options.lu_pivot_shift = options.lu_pivot_shift;
			topology_options.nonlinear_block_reduction = options.nonlinear_block_reduction;
			std::ostringstream text; text.exceptions(std::ios::badbit | std::ios::failbit);
			text << layout_.HashSha256() << ':' << std::setprecision(std::numeric_limits<double>::max_digits10)
				<< options.wall_inertial_gamma0 << ':' << options.include_pressure_gauge;
			signature = text.str();
		});
		RequireCollectiveSameText(communicator_,"transient operator configuration agreement",signature);
		core_ = AllocateCollectiveRuntime<ImmersedStaticDistributedOperator>(communicator_,communicator_,domain,volume,surface,ghost,
			topology_options,partition,ImmersedFlowTopologyMode::FixedTransient);
		std::vector<std::uint64_t> cells;
		CollectiveLocalStage(communicator_,"transient operator row binding",[&] {
			if (core_->Diagnostics().total_dofs != layout_.Rows()) throw std::logic_error("transient operator row count differs");
			for (auto id : layout_.NodeIds()) for (int c = 0; c < 4; ++c)
				if (core_->Topology().Dof(id,c) != static_cast<PetscInt>(4*layout_.LocalNode(id)+c)) throw std::logic_error("transient operator node row differs");
			for (const auto& p : core_->Diagnostics().ports) if (p.multiplier_row >= 0)
				if (p.multiplier_row != static_cast<PetscInt>(layout_.ControllerRow(p.boundary_label))) throw std::logic_error("transient operator controller row differs");
			cells = core_->OwnedVolumeCells();
		});
		inputs_ = AllocateCollectiveRuntime<ImmersedDistributedTransientVolume>(communicator_,communicator_,domain,volume,layout_,core_->Assembly().State(),cells);
	}
	ImmersedTransientDistributedOperator(const ImmersedTransientDistributedOperator&) = delete;
	ImmersedTransientDistributedOperator& operator=(const ImmersedTransientDistributedOperator&) = delete;
	const ImmersedActiveLayout& Layout() const noexcept { return layout_; }
	const ImmersedTransientFlowOptions& Options() const noexcept { return options_; }
	const ImmersedStaticFlowDiagnostics& Diagnostics() const noexcept { return core_->Diagnostics(); }
	ImmersedDistributedAssembly& Assembly() noexcept { return core_->Assembly(); }
	const ImmersedDistributedAssembly& Assembly() const noexcept { return core_->Assembly(); }
	const ImmersedDistributedTransientVolume& Inputs() const noexcept { return *inputs_; }
	void Freeze(Vec committed,const ImmersedActiveLayout& source_layout,double source_time,std::uint64_t source_index,
		double target_time,std::uint64_t target_index,double dt)
	{
		CollectiveLocalStage(communicator_,"transient operator freeze guard",[&] { RequireOpen(); core_->Topology().ValidateAllFlowCompatibility(); });
		auto parameters = options_.parameters; parameters.dt = dt;
		inputs_->Freeze(committed,source_layout,source_time,source_index,target_time,target_index,parameters,options_.body_force);
		options_.parameters.dt = dt;
	}
	void SetPortControlValue(const std::string& id,double value)
	{
		CollectiveLocalStage(communicator_,"transient operator port guard",[&] {
			RequireOpen(); if (inputs_->History().Active()) throw std::logic_error("cannot change frozen transient port controls");
		});
		core_->SetPortControlValue(id,value);
		for (auto& p : options_.ports) if (p.id == id) p.value = value;
	}
	void Assemble()
	{
		CollectiveLocalStage(communicator_,"transient operator assembly guard",[&] {
			RequireOpen(); if (!inputs_->History().Active()) throw std::logic_error("transient operator inputs are not frozen");
		});
		core_->AssembleWithVolume([&](std::uint64_t cell,const auto& nodal) { return BuildCell(cell,nodal); });
	}
	ImmersedStaticFlowConservationDiagnostics ConservationDiagnostics() const
	{
		CollectiveLocalStage(communicator_,"transient conservation guard",[&] { RequireOpen(); });
		return core_->ConservationDiagnostics();
	}
	void ReleaseTrial() noexcept { inputs_->ReleaseTrial(); options_.parameters.dt = 0; }
	void Close()
	{
		closed_ = true; options_.parameters.dt = 0; std::exception_ptr error;
		try { inputs_->Close(); } catch (...) { error = std::current_exception(); }
		try { core_->Close(); } catch (...) { if (!error) error = std::current_exception(); }
		if (error) std::rethrow_exception(error);
	}
private:
	friend class ImmersedMovingTransientDistributedOperator;
	ImmersedOwnedVolumeContribution BuildCell(std::uint64_t cell,const std::vector<std::array<double,4>>& nodal) const
	{
		return BuildCellWithVelocity(cell,nodal,
			[](const SurfaceQuadraturePoint&,const ImmersedSurfaceQuadraturePointProvenance&) { return std::array<double,3>{{0,0,0}}; });
	}
	ImmersedOwnedVolumeContribution BuildCellWithVelocity(std::uint64_t cell,const std::vector<std::array<double,4>>& nodal,
		const ImmersedMaterialWallVelocityEvaluator& wall_velocity) const
	{
		ImmersedOwnedVolumeContribution result;
		const auto volume = inputs_->BuildVolume(cell,nodal); result.system = volume;
		if (domain_.Cells()[cell].classification != CellClassification::Cut) return result;
		const auto element = domain_.Background().MaterializeElement(cell);
		const auto& rule = surface_.UsableRule(domain_,cell);
		const auto trace = BuildImmersedConservativeMixedTraceElement(element,rule,nodal);
		for (std::size_t i = 0; i < result.system.jacobian.size(); ++i) result.system.jacobian[i] += trace.jacobian[i];
		for (std::size_t i = 0; i < result.system.negative_residual.size(); ++i) result.system.negative_residual[i] += trace.negative_residual[i];
		if (core_->Topology().HasSelectedWallPoint(rule)) {
			const auto velocity = inputs_->History().Localize(element,inputs_->History().TargetTimeS());
			std::vector<std::array<double,4>> old(velocity.size());
			for (std::size_t i = 0; i < old.size(); ++i) for (int c = 0; c < 3; ++c) old[i][c] = velocity[i][c];
			const auto wall = BuildImmersedNitscheWallElementFromVolumeSystemMaterialAware(domain_,volume_,surface_,cell,nodal,old,
				options_.parameters,options_.wall_labels,volume,ghost_,options_.wall_gamma0,options_.wall_inertial_gamma0,wall_velocity);
			for (std::size_t i = 0; i < result.system.jacobian.size(); ++i) result.system.jacobian[i] += wall.system.jacobian[i]-volume.jacobian[i];
			for (std::size_t i = 0; i < result.system.negative_residual.size(); ++i) result.system.negative_residual[i] += wall.system.negative_residual[i]-volume.negative_residual[i];
			for (const auto& item : wall.diagnostics.by_boundary_id) result.wall_selected_points[item.first] = item.second.selected_points;
		}
		return result;
	}
	void RequireOpen() const { if (closed_) throw std::logic_error("transient operator is closed"); }
	MPI_Comm communicator_;
	const CartesianDomainClassification& domain_;
	const CutCellVolumeQuadratureCatalog& volume_;
	const ImmersedSurfaceQuadratureCatalog& surface_;
	const CutCellGhostPenaltyCatalog& ghost_;
	ImmersedTransientFlowOptions options_;
	ImmersedActiveLayout layout_;
	std::unique_ptr<ImmersedStaticDistributedOperator> core_;
	std::unique_ptr<ImmersedDistributedTransientVolume> inputs_;
	bool closed_ = false;
};

} // namespace iga
#endif
