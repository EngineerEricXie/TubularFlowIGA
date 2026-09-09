#ifndef IGA_IMMERSED_STATIC_DISTRIBUTED_OPERATOR_HPP
#define IGA_IMMERSED_STATIC_DISTRIBUTED_OPERATOR_HPP

#include "ImmersedStaticFlowSetup.hpp"
#include "ImmersedDistributedAssembly.hpp"
#include "RuntimeConstruction.hpp"
#include <chrono>
#include <iomanip>
#include <optional>
#include <memory>

namespace iga {

// Immutable geometry and controls for a steady four-field operator. Local
// evaluators must represent the same physical functions on every MPI member.
// Geometry/controls are validated before distributed resources are created.
// The communicator and catalogs must outlive this object. Construction,
// Assemble and Close are collective over that communicator; local callbacks
// must not call MPI. Assembly().State() is a borrowed distributed vector.
class ImmersedStaticDistributedOperator {
public:
	ImmersedStaticDistributedOperator(MPI_Comm communicator,
		const CartesianDomainClassification& domain, const CutCellVolumeQuadratureCatalog& volume,
		const ImmersedSurfaceQuadratureCatalog& surface, const CutCellGhostPenaltyCatalog& ghost,
		const ImmersedStaticFlowOptions& options)
		: communicator_(communicator), domain_(domain), volume_(volume), surface_(surface), ghost_(ghost)
	{
		MPI_Comm_rank(communicator_, &rank_); MPI_Comm_size(communicator_, &size_);
		std::vector<PetscInt> offsets;
		std::vector<ImmersedAssemblyStencil> stencils;
		std::string signature;
		CollectiveLocalStage(communicator_, "immersed static setup", [&] {
			setup_ = std::make_unique<ImmersedStaticFlowSetup>(domain,volume,surface,ghost,options);
			diagnostics_ = setup_->Diagnostics();
			local_ports_.resize(options.ports.size()); global_ports_.resize(options.ports.size());
			if (options.ports.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())/7)
				throw std::overflow_error("immersed port reduction exceeds MPI count capacity");
			for (int label : options.wall_labels) wall_labels_.emplace(label,0);
			for (const auto& entry : surface.Diagnostics().source_area_by_boundary_id) wall_labels_.emplace(entry.first,0);
			std::size_t label_index = 0;
			for (auto& entry : wall_labels_) entry.second = label_index++;
			if (wall_labels_.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())/2)
				throw std::overflow_error("immersed wall reduction exceeds MPI count capacity");
			local_walls_.resize(2*wall_labels_.size()); global_walls_.resize(local_walls_.size());
			signature = Signature();
			BuildTopology(offsets,stencils);
		});
		RequireCollectiveSameText(communicator_, "immersed static configuration agreement", signature);
		assembly_ = AllocateCollectiveRuntime<ImmersedDistributedAssembly>(communicator_,communicator_,offsets,stencils);
		try {
			double measure = 0.0;
			CollectiveLocalStage(communicator_, "immersed owned gauge weights", [&] {
				for (const auto& stencil : assembly_->OwnedStencils()) {
					const auto& task = work_[stencil.id]; if (task.kind != Kind::Volume) continue;
					const auto element = domain_.Background().MaterializeElement(task.cell);
					auto& weights = gauge_weights_[task.cell]; weights.assign(element.connectivity.size(),0.0);
					VolumePoints(task.cell,[&](const VolumeQuadraturePoint& point) {
						const auto basis = EvaluateBasis(element,point.parametric[0],point.parametric[1],point.parametric[2],false);
						const double weight = point.weight*basis.raw_determinant;
						if (!std::isfinite(weight) || !(weight > 0.0)) throw std::runtime_error("invalid immersed gauge measure");
						measure += weight;
						for (std::size_t a = 0; a < weights.size(); ++a) weights[a] += basis.value[a]*weight;
					});
				}
				if (!std::isfinite(measure)) throw std::runtime_error("nonfinite immersed gauge measure");
			});
			MPI_Allreduce(&measure,&diagnostics_.pressure_measure,1,MPI_DOUBLE,MPI_SUM,communicator_);
			CollectiveLocalStage(communicator_, "immersed global gauge measure", [&] {
				if (!std::isfinite(diagnostics_.pressure_measure) || !(diagnostics_.pressure_measure > 0.0)) throw std::runtime_error("invalid immersed global gauge measure");
			});
			Check("immersed pressure constant create",VecDuplicate(assembly_->State(),&constant_pressure_));
			Check("immersed pressure defect create",VecDuplicate(assembly_->State(),&pressure_defect_));
			Check("immersed pressure constant reset",VecSet(constant_pressure_,0.0));
			CollectiveLocalStage(communicator_, "immersed pressure constant entries", [&] {
				for (PetscInt row = assembly_->RowBegin(); row < assembly_->RowEnd(); ++row)
					if (static_cast<std::size_t>(row) < diagnostics_.physical_dofs && row%4 == 3)
						if (VecSetValue(constant_pressure_,row,1.0,INSERT_VALUES)) throw std::runtime_error("cannot set immersed pressure constant");
			});
			Check("immersed pressure constant begin",VecAssemblyBegin(constant_pressure_));
			Check("immersed pressure constant end",VecAssemblyEnd(constant_pressure_));
		} catch (...) { ReleaseVectors(); throw; }
	}
	~ImmersedStaticDistributedOperator() { ReleaseVectors(); }
	ImmersedStaticDistributedOperator(const ImmersedStaticDistributedOperator&) = delete;
	ImmersedStaticDistributedOperator& operator=(const ImmersedStaticDistributedOperator&) = delete;
	const ImmersedStaticFlowSetup& Topology() const noexcept { return *setup_; }
	const ImmersedStaticFlowOptions& Options() const noexcept { return setup_->Options(); }
	const ImmersedStaticFlowDiagnostics& Diagnostics() const noexcept { return diagnostics_; }
	ImmersedDistributedAssembly& Assembly() noexcept { return *assembly_; }
	const ImmersedDistributedAssembly& Assembly() const noexcept { return *assembly_; }

	void Assemble()
	{
		PhaseScope phase(ProfilePhase::Assembly);
		const auto start = std::chrono::steady_clock::now();
		std::fill(local_ports_.begin(),local_ports_.end(),std::array<double,7>{}); local_counts_.fill(0);
		std::fill(local_walls_.begin(),local_walls_.end(),0);
		assembly_->Assemble([&](const auto& stencil,auto& matrix,auto& residual) { Integrate(stencil,matrix,residual); });
		MPI_Allreduce(local_counts_.data(),global_counts_.data(),3,MPI_UINT64_T,MPI_SUM,communicator_);
		if (!local_ports_.empty()) MPI_Allreduce(local_ports_.data(),global_ports_.data(),static_cast<int>(local_ports_.size()*7),MPI_DOUBLE,MPI_SUM,communicator_);
		if (!local_walls_.empty()) MPI_Allreduce(local_walls_.data(),global_walls_.data(),static_cast<int>(local_walls_.size()),MPI_UINT64_T,MPI_SUM,communicator_);
		Check("immersed constant pressure action",MatMult(assembly_->Matrix(),constant_pressure_,pressure_defect_));
		double local_defect = 0.0, defect = 0.0;
		CollectiveLocalStage(communicator_, "immersed pressure defect measurement", [&] {
			PetscReadArray view; view.Acquire(pressure_defect_);
			for (PetscInt row = assembly_->RowBegin(); row < assembly_->RowEnd(); ++row)
				if (static_cast<std::size_t>(row) < diagnostics_.physical_dofs) {
					const double value = PetscRealPart(view.Data()[row-assembly_->RowBegin()]);
					if (!std::isfinite(value)) throw std::runtime_error("nonfinite immersed pressure action");
					local_defect = std::max(local_defect,std::abs(value));
				}
			view.Restore();
		});
		MPI_Allreduce(&local_defect,&defect,1,MPI_DOUBLE,MPI_MAX,communicator_);
		std::optional<ImmersedStaticFlowDiagnostics> candidate;
		CollectiveLocalStage(communicator_, "immersed static global diagnostics", [&] {
			candidate = diagnostics_;
			candidate->volume_cells = global_counts_[0]; candidate->surface_cells = global_counts_[1]; candidate->ghost_faces = global_counts_[2];
			candidate->constant_pressure_defect = defect;
			candidate->wall_selected_points = 0; candidate->wall_selected_points_by_label.clear();
			for (const auto& entry : wall_labels_) if (global_walls_[2*entry.second+1]) {
				const auto count = global_walls_[2*entry.second];
				candidate->wall_selected_points_by_label[entry.first] = count; candidate->wall_selected_points += count;
			}
			for (PetscInt row = std::max(assembly_->RowBegin(),static_cast<PetscInt>(diagnostics_.physical_dofs)); row < assembly_->RowEnd(); ++row) {
				PetscInt count = 0; const PetscInt* columns = nullptr; const PetscScalar* values = nullptr;
				if (MatGetRow(assembly_->Matrix(),row,&count,&columns,&values)) throw std::runtime_error("cannot inspect immersed scalar row");
				bool diagonal = false;
				for (PetscInt j = 0; j < count; ++j) if (columns[j] == row) diagonal = values[j] == 0.0;
				if (MatRestoreRow(assembly_->Matrix(),row,&count,&columns,&values)) throw std::runtime_error("cannot restore immersed scalar row");
				if (!diagonal) throw std::runtime_error("immersed scalar structural zero diagonal missing");
			}
			candidate->scalar_diagonal_structure_verified = true;
			if (!std::isfinite(defect)) throw std::runtime_error("nonfinite immersed pressure defect");
			for (std::size_t port = 0; port < global_ports_.size(); ++port) {
				const auto& values = global_ports_[port];
				for (double value : values) if (!std::isfinite(value)) throw std::runtime_error("nonfinite immersed global port measurement");
				if (!(values[0] > 0.0)) throw std::runtime_error("immersed port area is not positive");
				auto& diagnostic = candidate->ports[port];
				diagnostic.measurement = {values[0],values[1],values[2]/values[0],values[3]/values[0],values[4]/values[0]};
				if (!std::isfinite(diagnostic.measurement.mean_pressure_pa) || !std::isfinite(diagnostic.measurement.mean_normal_traction_pa)
					|| !std::isfinite(diagnostic.measurement.mean_velocity_squared_m2_s2)) throw std::runtime_error("nonfinite immersed port mean");
				diagnostic.assembled_surface_points = static_cast<std::size_t>(values[5]);
				if (diagnostic.multiplier_row >= 0) {
					diagnostic.multiplier = values[6]; diagnostic.constraint_residual = diagnostic.target-values[1];
					diagnostic.absolute_flow_residual_m3_s = std::abs(diagnostic.constraint_residual);
					diagnostic.flow_tolerance_m3_s = Options().flow_controller_absolute_tolerance_m3_s
						+Options().flow_controller_relative_tolerance*std::max(std::abs(diagnostic.target),Options().flow_controller_reference_flow_m3_s);
					if (!std::isfinite(diagnostic.flow_tolerance_m3_s) || !(diagnostic.flow_tolerance_m3_s > 0.0)) throw std::runtime_error("invalid immersed flow tolerance");
					diagnostic.normalized_flow_residual = diagnostic.absolute_flow_residual_m3_s/diagnostic.flow_tolerance_m3_s;
					if (!std::isfinite(diagnostic.normalized_flow_residual)) throw std::runtime_error("nonfinite immersed controller residual");
				}
			}
			candidate->last_assembly_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
			candidate->aggregate_assembly_seconds += candidate->last_assembly_seconds;
		});
		diagnostics_ = std::move(*candidate);
	}
	ImmersedStaticFlowConservationDiagnostics ConservationDiagnostics() const
	{
		PhaseScope phase(ProfilePhase::Diagnostics);
		std::vector<double> local,global;
		std::vector<std::uint64_t> seen,global_seen;
		CollectiveLocalStage(communicator_,"immersed conservation storage",[&] {
			local.resize(wall_labels_.size()+1); global.resize(local.size());
			seen.resize(wall_labels_.size()); global_seen.resize(seen.size());
		});
		assembly_->WithRequiredState([&] {
			for (const auto& stencil : assembly_->OwnedStencils()) {
				const auto& task = work_[stencil.id]; if (task.kind != Kind::Volume) continue;
				const auto element = domain_.Background().MaterializeElement(task.cell); const auto nodal = Gather(element);
				VolumePoints(task.cell,[&](const VolumeQuadraturePoint& point) {
					const auto basis = EvaluateBasis(element,point.parametric[0],point.parametric[1],point.parametric[2],false);
					double divergence = 0.0;
					for (std::size_t a = 0; a < nodal.size(); ++a) for (int field = 0; field < 3; ++field)
						AddImmersedFlowPortFinite(divergence,nodal[a][field]*basis.gradient[a][field],"distributed conservation divergence");
					AddImmersedFlowPortFinite(local.back(),point.weight*basis.raw_determinant*divergence,"distributed volume divergence");
				});
				if (domain_.Cells()[task.cell].classification != CellClassification::Cut) continue;
				for (const auto& point : surface_.UsableRule(domain_,task.cell).Points()) {
					const auto basis = EvaluateBasis(element,point.parametric[0],point.parametric[1],point.parametric[2],false);
					double flow = 0.0;
					for (std::size_t a = 0; a < nodal.size(); ++a) for (int field = 0; field < 3; ++field)
						AddImmersedFlowPortFinite(flow,nodal[a][field]*basis.value[a]*point.normal[field]*point.weight,"distributed surface flow");
					const auto index = wall_labels_.at(point.boundary_id);
					AddImmersedFlowPortFinite(local[index],flow,"distributed boundary flow"); ++seen[index];
				}
			}
		});
		MPI_Allreduce(local.data(),global.data(),static_cast<int>(local.size()),MPI_DOUBLE,MPI_SUM,communicator_);
		MPI_Allreduce(seen.data(),global_seen.data(),static_cast<int>(seen.size()),MPI_UINT64_T,MPI_SUM,communicator_);
		ImmersedStaticFlowConservationDiagnostics result;
		CollectiveLocalStage(communicator_,"immersed global conservation",[&] {
			for (double value : global) if (!std::isfinite(value)) throw std::runtime_error("nonfinite immersed global conservation");
			result.volume_divergence_integral_m3_s = global.back();
			for (const auto& entry : wall_labels_) if (global_seen[entry.second]) {
				const auto flow = global[entry.second]; result.surface_flow_by_boundary_label_m3_s[entry.first] = flow;
				AddImmersedFlowPortFinite(result.total_surface_outward_flow_m3_s,flow,"global conservation total surface flow");
				if (std::binary_search(Options().wall_labels.begin(),Options().wall_labels.end(),entry.first))
					AddImmersedFlowPortFinite(result.wall_outward_flow_m3_s,flow,"global conservation wall flow");
				else if (std::any_of(Options().ports.begin(),Options().ports.end(),[&](const auto& port) { return port.boundary_label == entry.first; }))
					AddImmersedFlowPortFinite(result.open_port_outward_flow_m3_s,flow,"global conservation open-port flow");
				else throw std::runtime_error("immersed conservation diagnostics found an unconfigured surface label");
			}
		});
		return result;
	}
	void Close()
	{
		ReleaseVectors();
		assembly_->Close();
		cleanup_.Check(communicator_,"immersed static operator close");
	}
private:
	enum class Kind { Volume, Ghost, Gauge, Port, Measurement, Target };
	struct Work { Kind kind; std::uint64_t cell; std::size_t port; };
	template <class Function>
	void VolumePoints(std::uint64_t cell,Function&& function) const
	{
		if (volume_.StorageMode() == CutCellVolumeQuadratureStorageMode::Expanded)
			for (const auto& point : volume_.UsableRule(domain_,cell).Points()) function(point);
		else ForEachVolumePoint(volume_.UsableCompactRule(domain_,cell),function);
	}
	std::string Signature() const
	{
		std::ostringstream text; text.exceptions(std::ios::badbit | std::ios::failbit);
		text << std::setprecision(std::numeric_limits<double>::max_digits10) << domain_.SurfaceCanonicalHash();
		const auto& grid = domain_.Background().Spec();
		for (int d = 0; d < 3; ++d) text << ' ' << grid.lower_m[d] << ' ' << grid.upper_m[d] << ' ' << grid.cells[d];
		const auto& q = volume_.Options(); const auto& s = surface_.Options(); const auto& g = ghost_.Options(); const auto& o = Options();
		text << ' ' << static_cast<int>(volume_.StorageMode()) << ' ' << q.max_depth << ' ' << q.max_nodes << ' ' << q.max_leaves << ' ' << q.max_points
			<< ' ' << q.max_records << ' ' << q.max_retained_bytes << ' ' << q.max_logical_points << ' ' << q.empty_rule_rescue_max_depth
			<< ' ' << s.max_candidates << ' ' << s.max_fragments << ' ' << s.max_points << ' ' << s.max_exact_limbs
			<< ' ' << g.gamma_u << ' ' << g.gamma_p << ' ' << g.max_faces << ' ' << g.max_quadrature_points << ' ' << g.max_trace_entries
			<< ' ' << o.parameters.density << ' ' << o.parameters.dynamic_viscosity << ' ' << o.parameters.dt << ' ' << o.wall_gamma0
			<< ' ' << o.assemble_volume << ' ' << o.assemble_wall << ' ' << o.assemble_ghost << ' ' << o.assemble_gauge
			<< ' ' << o.nonlinear_maximum_iterations << ' ' << o.ksp_maximum_iterations << ' ' << o.ksp_relative_tolerance
			<< ' ' << o.nonlinear_relative_tolerance << ' ' << o.nonlinear_absolute_tolerance << ' ' << o.nonlinear_block_reduction << ' ' << o.minimum_damping << ' ' << o.lu_pivot_shift
			<< ' ' << o.flow_controller_relative_tolerance << ' ' << o.flow_controller_absolute_tolerance_m3_s << ' ' << o.flow_controller_reference_flow_m3_s;
		for (int label : o.wall_labels) text << " w" << label;
		for (const auto& port : o.ports) text << " p" << port.id.size() << ':' << port.id << ':' << port.boundary_label << ':' << static_cast<int>(port.control_mode) << ':' << port.value;
		return text.str();
	}
	void BuildTopology(std::vector<PetscInt>& offsets,std::vector<ImmersedAssemblyStencil>& stencils)
	{
		const auto nodes = setup_->ActiveNodes().size();
		for (int r = 0; r <= size_; ++r) offsets.push_back(r == size_ ? static_cast<PetscInt>(diagnostics_.total_dofs) : static_cast<PetscInt>(4*(nodes*r/size_)));
		std::vector<std::uint64_t> active;
		std::vector<int> owners(domain_.Cells().size(),-1);
		for (std::uint64_t cell = 0; cell < domain_.Cells().size(); ++cell) if (setup_->UsablePositive(cell)) active.push_back(cell);
		const auto add = [&](int owner,ImmersedStencilPattern pattern,std::vector<PetscInt> rows,Work work) {
			stencils.push_back({stencils.size(),owner,pattern,std::move(rows)}); work_.push_back(work);
		};
		for (std::size_t i = 0; i < active.size(); ++i) {
			const auto cell = active[i]; const int owner = static_cast<int>(i*size_/active.size()); owners[cell] = owner;
			const auto element = domain_.Background().MaterializeElement(cell);
			std::vector<PetscInt> rows,pressure,velocity;
			for (auto node : element.connectivity) {
				for (int f = 0; f < 4; ++f) rows.push_back(setup_->Dof(node,f));
				pressure.push_back(setup_->Dof(node,3));
				for (int f = 0; f < 3; ++f) velocity.push_back(setup_->Dof(node,f));
			}
			add(owner,ImmersedStencilPattern::Dense,rows,{Kind::Volume,cell,0});
			if (setup_->HasGauge()) { pressure.push_back(setup_->GaugeDof()); add(owner,ImmersedStencilPattern::Scalar,pressure,{Kind::Gauge,cell,0}); }
			const auto& rule = surface_.UsableRule(domain_,cell);
			for (std::size_t port = 0; port < Options().ports.size(); ++port) {
				if (!setup_->RuleHasLabel(rule,Options().ports[port].boundary_label)) continue;
				const auto scalar = diagnostics_.ports[port].multiplier_row;
				auto port_rows = scalar < 0 ? rows : velocity;
				if (scalar >= 0) port_rows.push_back(scalar);
				add(owner,scalar < 0 ? ImmersedStencilPattern::Dense : ImmersedStencilPattern::Scalar,std::move(port_rows),{Kind::Port,cell,port});
			}
		}
		// The serial port audit visits every surface rule, including a certified-
		// empty cut cell. Keep that measurement when its nodes are represented
		// by neighboring active cells; Dof rejects an unrepresented node just
		// as the serial Gather does. This stencil adds no physical coefficients.
		for (std::uint64_t cell = 0; cell < domain_.Cells().size(); ++cell) if (owners[cell] < 0) {
			const auto& rule = surface_.UsableRule(domain_,cell);
			for (std::size_t port = 0; port < Options().ports.size(); ++port) {
				if (!setup_->RuleHasLabel(rule,Options().ports[port].boundary_label)) continue;
				std::vector<PetscInt> rows;
				for (auto node : domain_.Background().MaterializeElement(cell).connectivity)
					for (int f = 0; f < 4; ++f) rows.push_back(setup_->Dof(node,f));
				add(static_cast<int>(cell%static_cast<std::uint64_t>(size_)),ImmersedStencilPattern::Dense,std::move(rows),{Kind::Measurement,cell,port});
			}
		}
		for (std::size_t face = 0; face < ghost_.Faces().size(); ++face) {
			const auto& pair = ghost_.Faces()[face];
			std::vector<PetscInt> rows;
			for (auto node : CubicCartesianSplineFaceConnectivity(domain_,pair.minus_cell,pair.plus_cell))
				for (int f = 0; f < 4; ++f) rows.push_back(setup_->Dof(node,f));
			add(owners.at(pair.minus_cell),ImmersedStencilPattern::SameField,std::move(rows),{Kind::Ghost,face,0});
		}
		for (std::size_t port = 0; port < Options().ports.size(); ++port)
			if (diagnostics_.ports[port].multiplier_row >= 0)
				add(size_-1,ImmersedStencilPattern::Scalar,{diagnostics_.ports[port].multiplier_row},{Kind::Target,0,port});
	}
	std::vector<std::array<double,4>> Gather(const Element& element) const
	{
		std::vector<std::array<double,4>> nodal(element.connectivity.size());
		for (std::size_t a = 0; a < nodal.size(); ++a)
			for (int f = 0; f < 4; ++f) nodal[a][f] = assembly_->StateAt(setup_->Dof(element.connectivity[a],f));
		return nodal;
	}
	void Integrate(const ImmersedAssemblyStencil& stencil,std::vector<PetscScalar>& matrix,std::vector<PetscScalar>& residual)
	{
		const auto& task = work_[stencil.id]; const auto n = stencil.rows.size(); const auto& options = Options();
		matrix.assign(n*n,0.0); residual.assign(n,0.0);
		if (task.kind == Kind::Ghost) {
			if (!options.assemble_ghost) return;
			auto block = ghost_.AssembleFaceLocal(task.cell,domain_,volume_,[&](std::int32_t node,int f) { return assembly_->StateAt(setup_->Dof(node,f)); },options.parameters.dynamic_viscosity);
			matrix = std::move(block.jacobian); residual = std::move(block.negative_residual); ++local_counts_[2]; return;
		}
		if (task.kind == Kind::Target) {
			residual[0] = options.ports[task.port].value;
			local_ports_[task.port][6] = assembly_->StateAt(stencil.rows[0]); return;
		}
		if (task.kind == Kind::Gauge) {
			if (!options.assemble_gauge) return;
			const auto& coefficients = gauge_weights_.at(task.cell); const double lambda = assembly_->StateAt(stencil.rows.back());
			for (std::size_t a = 0; a < coefficients.size(); ++a) {
				matrix[a*n+n-1] = matrix[(n-1)*n+a] = coefficients[a]; residual[a] = -coefficients[a]*lambda;
				residual.back() -= coefficients[a]*assembly_->StateAt(stencil.rows[a]);
			}
			return;
		}
		const auto element = domain_.Background().MaterializeElement(task.cell); const auto nodal = Gather(element);
		const auto& rule = surface_.UsableRule(domain_,task.cell);
		if (task.kind == Kind::Volume) {
			auto volume = BuildNavierStokesElementFromPoints(element,nodal,{},options.parameters,[&](const auto& consume) { VolumePoints(task.cell,consume); },options.body_force,NavierStokesResolvedMixedForm::Conservative);
			if (options.assemble_volume) { matrix = volume.jacobian; residual = volume.negative_residual; ++local_counts_[0]; }
			if (domain_.Cells()[task.cell].classification == CellClassification::Cut) {
				if (options.assemble_volume) {
					const auto trace = BuildImmersedConservativeMixedTraceElement(element,rule,nodal);
					for (std::size_t i = 0; i < matrix.size(); ++i) matrix[i] += trace.jacobian[i];
					for (std::size_t i = 0; i < residual.size(); ++i) residual[i] += trace.negative_residual[i];
				}
				if (options.ports.empty() || setup_->HasSelectedWallPoint(rule)) {
					const auto wall = BuildImmersedNitscheWallElementFromVolumeSystem(domain_,volume_,surface_,task.cell,nodal,{},options.parameters,options.wall_labels,volume,ghost_,options.wall_gamma0,options.wall_velocity);
					std::size_t wall_points = 0;
					for (const auto& entry : wall.diagnostics.by_boundary_id) wall_points += entry.second.selected_points;
					if (!wall_points) throw std::runtime_error("selected immersed wall labels produced no surface contribution");
					if (options.assemble_wall) {
						for (std::size_t i = 0; i < matrix.size(); ++i) matrix[i] += wall.system.jacobian[i]-volume.jacobian[i];
						for (std::size_t i = 0; i < residual.size(); ++i) residual[i] += wall.system.negative_residual[i]-volume.negative_residual[i];
						++local_counts_[1];
						for (const auto& entry : wall.diagnostics.by_boundary_id) {
							const auto index = 2*wall_labels_.at(entry.first);
							local_walls_[index] += entry.second.selected_points; ++local_walls_[index+1];
						}
					}
				}
			}
			return;
		}
		const auto& port = options.ports[task.port];
		const auto measure = MeasureImmersedFlowPortElement(element,rule,port.boundary_label,nodal,options.parameters.dynamic_viscosity);
		auto& sums = local_ports_[task.port]; sums[0] += measure.area_m2; sums[1] += measure.outward_flow_m3_s;
		sums[2] += measure.mean_pressure_pa*measure.area_m2; sums[3] += measure.mean_normal_traction_pa*measure.area_m2;
		sums[4] += measure.mean_velocity_squared_m2_s2*measure.area_m2;
		if (task.kind == Kind::Measurement) return;
		const auto block = BuildImmersedFlowPortElement(element,rule,port.boundary_label,port.control_mode,port.value,nodal);
		for (const auto& point : rule.Points()) if (point.boundary_id == port.boundary_label) ++sums[5];
		if (IsImmersedFlowPressureLike(port.control_mode)) { residual = block.negative_residual; return; }
		const double lambda = assembly_->StateAt(stencil.rows.back());
		for (std::size_t a = 0; a < element.connectivity.size(); ++a) for (int field = 0; field < 3; ++field) {
			const auto i = 3*a+field; const double coefficient = PetscRealPart(block.flow_coefficient[4*a+field]);
			matrix[i*n+n-1] = matrix[(n-1)*n+i] = coefficient; residual[i] = -coefficient*lambda;
		}
		residual.back() = -measure.outward_flow_m3_s;
	}
	void Check(const char* stage,PetscErrorCode code) const { RequireCollectivePetscSuccess(communicator_,stage,code); }
	void ReleaseVectors() noexcept
	{
		if (pressure_defect_) { cleanup_.Observe("pressure defect destroy",VecDestroy(&pressure_defect_)); pressure_defect_ = nullptr; }
		if (constant_pressure_) { cleanup_.Observe("pressure constant destroy",VecDestroy(&constant_pressure_)); constant_pressure_ = nullptr; }
	}
	MPI_Comm communicator_;
	int rank_ = 0,size_ = 1;
	const CartesianDomainClassification& domain_;
	const CutCellVolumeQuadratureCatalog& volume_;
	const ImmersedSurfaceQuadratureCatalog& surface_;
	const CutCellGhostPenaltyCatalog& ghost_;
	std::unique_ptr<ImmersedStaticFlowSetup> setup_;
	std::unique_ptr<ImmersedDistributedAssembly> assembly_;
	std::vector<Work> work_;
	std::map<std::uint64_t,std::vector<double>> gauge_weights_;
	std::map<int,std::size_t> wall_labels_;
	std::vector<std::uint64_t> local_walls_,global_walls_;
	std::vector<std::array<double,7>> local_ports_,global_ports_;
	std::array<std::uint64_t,3> local_counts_{},global_counts_{};
	ImmersedStaticFlowDiagnostics diagnostics_{};
	Vec constant_pressure_ = nullptr,pressure_defect_ = nullptr;
	RuntimeCleanupResult cleanup_;
};

} // namespace iga
#endif
