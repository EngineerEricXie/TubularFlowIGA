#ifndef IGA_IMMERSED_FLOW_CASE_HPP
#define IGA_IMMERSED_FLOW_CASE_HPP

// Production owner for steady or fixed-geometry transient immersed flow. Keep
// the catalog dependencies in declaration order: the runtime borrows ghost,
// surface and volume, ghost borrows classification and volume, and all of
// them ultimately retain the classified surface.
#include "CheckedText.hpp"
#include "CaseConfig.hpp"
#include "FlowDomainPortMetadata.hpp"
#include "SimulationConfig.hpp"
#include "ThreeDImmersedFlowDomain.hpp"
#include "ThreeDImmersedDistributedFlowDomain.hpp"
#include "ThreeDImmersedTransientDistributedFlowDomain.hpp"
#include "Sha256.hpp"
#include <type_traits>
#include "../solvers/cpu/include/CartesianDomainClassification.hpp"
#include "../solvers/cpu/include/CutCellGhostPenalty.hpp"
#include "../solvers/cpu/include/ImmersedStaticFlowRuntime.hpp"
#include "../solvers/cpu/include/ImmersedSurfaceQuadrature.hpp"
#include "../solvers/cpu/include/SurfaceReaders.hpp"

#include <filesystem>
#include <fstream>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

struct ImmersedCaseDistribution {
	std::size_t global_rows = 0,owned_rows = 0,owned_stencils = 0,required_rows = 0;
};

class ImmersedFlowCase final : public CoupledDomainRuntime {
public:
	static std::unique_ptr<ImmersedFlowCase> Load(const std::filesystem::path& case_directory,
		const std::string& domain_id,const std::vector<CouplingPort>& ports,int mpi_size)
	{
		if (mpi_size != 1) throw std::runtime_error("immersed flow cases require MPI size 1 through the serial Load interface");
		auto result = Preflight(case_directory,domain_id,ports);
		if (result->IsTransient()) {
			result->InitializeDistributed(PETSC_COMM_SELF);
			return result;
		}
		result->runtime_ = std::make_unique<ImmersedStaticFlowRuntime>(*result->classification_,
			*result->volume_,*result->surface_,*result->ghost_,result->runtime_options_);
		result->adapter_ = std::make_unique<ThreeDImmersedFlowDomain>(domain_id,*result->runtime_,ports);
		return result;
	}
	// Local geometry/configuration preparation: safe inside a caller's local
	// failure stage. No distributed resource is created until initialization.
	static std::unique_ptr<ImmersedFlowCase> Preflight(const std::filesystem::path& case_directory,
		const std::string& domain_id,const std::vector<CouplingPort>& ports)
	{
		PhaseScope input_phase(ProfilePhase::Input);
		ValidateThreeDImmersedFlowDomainMetadata(domain_id, ports);
		auto result = std::unique_ptr<ImmersedFlowCase>(new ImmersedFlowCase);
		result->case_directory_ = std::filesystem::canonical(case_directory);
		if (!std::filesystem::is_directory(result->case_directory_))
			throw std::runtime_error("immersed flow case directory is not a directory");
		result->configuration_ = ReadSimulationConfiguration(
			Contained(result->case_directory_, "simulation_config.json").string());
		result->transient_ = ValidateSimulation(result->configuration_);
		const auto root = ParseJson(Contained(result->case_directory_, "immersed_geometry.json"));
		const auto& object = config_detail::RequireObject(root, "immersed_geometry.json");
		config_detail::RequireKnownKeys(object, {"surface", "boundary_array", "grid",
			"volume_quadrature", "surface_quadrature", "ghost_penalty", "wall_labels",
			"runtime"}, "immersed_geometry.json");
		const std::string surface_file = config_detail::RequireString(
			Required(object, "surface"), "immersed_geometry.json.surface");
		const std::string boundary_array = OptionalString(object, "boundary_array", "boundary_id");
		const auto surface_path = Contained(result->case_directory_, surface_file);
		input_phase.Stop();
		PhaseScope geometry_phase(ProfilePhase::Geometry);
		result->classification_ = std::make_unique<CartesianDomainClassification>(
			CubicCartesianBackground(ParseGrid(Required(object, "grid"))),
			SurfaceSpatialIndex(SurfaceReaders::ReadVtpPath(surface_path.string(), {}, boundary_array)));
		const auto volume_options = ParseVolume(Required(object, "volume_quadrature"));
		const auto storage = ParseStorage(Required(object, "volume_quadrature"));
		result->volume_ = std::make_unique<CutCellVolumeQuadratureCatalog>(*result->classification_,
			volume_options, storage);
		result->surface_ = std::make_unique<ImmersedSurfaceQuadratureCatalog>(*result->classification_,
			ParseSurface(Required(object, "surface_quadrature")));
		result->ghost_ = std::make_unique<CutCellGhostPenaltyCatalog>(*result->classification_,
			*result->volume_, ParseGhost(Required(object, "ghost_penalty")));
		const auto wall_labels = ParseLabels(Required(object, "wall_labels"));
		if (result->transient_) {
			result->transient_options_ = ParseRuntime<ImmersedTransientFlowOptions>(Required(object,"runtime"),result->configuration_,wall_labels);
			result->runtime_parameters_ = result->transient_options_.parameters;
			result->runtime_parameters_.dt = result->configuration_.time.dt;
		} else {
			result->runtime_options_ = ParseRuntime<ImmersedStaticFlowOptions>(Required(object,"runtime"),result->configuration_,wall_labels);
			result->runtime_parameters_ = result->runtime_options_.parameters;
		}
		ValidateLabelPartition(result->classification_->SurfaceIndex().Surface(),wall_labels,ports,
			result->transient_ ? result->transient_options_.ports : result->runtime_options_.ports);
		result->domain_id_ = domain_id; result->graph_ports_ = ports;
		result->geometry_identity_ = result->FixedGeometryIdentity();
		return result;
	}
	void InitializeDistributed(MPI_Comm communicator)
	{
		std::string signature;
		CollectiveLocalStage(communicator,"immersed case distributed initialization",[&] {
			if (runtime_ || distributed_runtime_ || transient_runtime_ || adapter_ || !classification_ || !volume_ || !surface_ || !ghost_)
				throw std::logic_error("immersed case is not an uninitialized preflight result");
			std::ostringstream text; text.exceptions(std::ios::badbit | std::ios::failbit);
			text << std::setprecision(std::numeric_limits<double>::max_digits10) << transient_;
			if (transient_) text << ':' << configuration_.time.dt << ':' << configuration_.time.steps;
			signature = text.str();
		});
		RequireCollectiveSameText(communicator,"immersed case time integration agreement",signature);
		if (transient_) {
			auto runtime = AllocateCollectiveRuntime<ImmersedTransientDistributedRuntime>(communicator,communicator,
				*classification_,*volume_,*surface_,*ghost_,geometry_identity_,transient_options_);
			auto adapter = AllocateCollectiveRuntime<ThreeDImmersedTransientDistributedFlowDomain>(communicator,domain_id_,*runtime,graph_ports_);
			transient_runtime_ = std::move(runtime); adapter_ = std::move(adapter); return;
		}
		auto runtime = AllocateCollectiveRuntime<ImmersedStaticDistributedRuntime>(communicator,communicator,
			*classification_,*volume_,*surface_,*ghost_,runtime_options_);
		auto adapter = AllocateCollectiveRuntime<ThreeDImmersedDistributedFlowDomain>(communicator,domain_id_,*runtime,graph_ports_);
		distributed_runtime_ = std::move(runtime); adapter_ = std::move(adapter);
	}
	bool IsTransient() const noexcept { return transient_; }
	bool IsDistributed() const noexcept { return distributed_runtime_ || transient_runtime_; }
	ImmersedCaseDistribution Distribution() const
	{
		if (transient_runtime_) return DistributionOf(*transient_runtime_);
		if (distributed_runtime_) return DistributionOf(*distributed_runtime_);
		throw std::logic_error("immersed case has no distributed ownership");
	}
	// Collective for a distributed case; geometry and audit metadata stay valid.
	void CloseDistributed() const
	{
		if (transient_runtime_) transient_runtime_->Close();
		if (distributed_runtime_) distributed_runtime_->Close();
	}
	ImmersedTransientDistributedRuntime& TransientRuntime()
	{
		if (!transient_runtime_) throw std::logic_error("immersed case has no transient runtime");
		return *transient_runtime_;
	}
	const ImmersedTransientDistributedRuntime& TransientRuntime() const
	{
		if (!transient_runtime_) throw std::logic_error("immersed case has no transient runtime");
		return *transient_runtime_;
	}
	ImmersedStaticDistributedRuntime& DistributedRuntime()
	{
		if (!distributed_runtime_) throw std::logic_error("immersed case has no distributed runtime");
		return *distributed_runtime_;
	}

	ImmersedStaticFlowRuntime& Runtime() { if (!runtime_) throw std::logic_error("immersed case has no serial runtime"); return *runtime_; }
	const ImmersedStaticFlowRuntime& Runtime() const { if (!runtime_) throw std::logic_error("immersed case has no serial runtime"); return *runtime_; }
	const SimulationConfiguration& Configuration() const noexcept { return configuration_; }
	const std::string& SurfaceHash() const noexcept { return classification_->SurfaceCanonicalHash(); }
	const CubicCartesianGridSpec& Grid() const noexcept { return classification_->Background().Spec(); }
	const CartesianDomainDiagnostics& ClassificationDiagnostics() const noexcept
	{ return classification_->Diagnostics(); }
	const CutCellVolumeQuadratureDiagnostics& VolumeDiagnostics() const noexcept
	{ return volume_->Diagnostics(); }
	const ImmersedSurfaceQuadratureDiagnostics& SurfaceDiagnostics() const noexcept
	{ return surface_->Diagnostics(); }
	const CutCellGhostPenaltyDiagnostics& GhostDiagnostics() const noexcept
	{ return ghost_->Diagnostics(); }
	const NavierStokesParameters& RuntimeParameters() const noexcept
	{ return runtime_parameters_; }

	const std::string& DomainId() const noexcept override { return domain_id_; }
	DomainKind Kind() const noexcept override { return DomainKind::ThreeDImmersedFlow; }
	const std::vector<CouplingPort>& Ports() const noexcept override { return graph_ports_; }
	void BeginStep(const DomainStepContext& step) override
	{
		if (transient_runtime_) CollectiveLocalStage(transient_runtime_->Communicator(),"immersed case timestep",[&] {
			if (step.dt_s != configuration_.time.dt) throw std::invalid_argument("transient immersed timestep differs from case configuration");
		});
		Adapter().BeginStep(step);
	}
	void SetPortInput(const std::string& id, const PortBoundaryData& input) override
	{ Adapter().SetPortInput(id, input); }
	void SolveTrial() override { Adapter().SolveTrial(); }
	PortState GetPortState(const std::string& id) const override { return Adapter().GetPortState(id); }
	void RollbackTrial() override { Adapter().RollbackTrial(); }
	void AbortStep() override { Adapter().AbortStep(); }
	void PrepareCommitStep() override { Adapter().PrepareCommitStep(); }
	void FinalizeCommitStep() noexcept override { if (adapter_) adapter_->FinalizeCommitStep(); }

private:
	template<class Runtime> static ImmersedCaseDistribution DistributionOf(const Runtime& runtime)
	{
		return {runtime.Diagnostics().total_dofs,static_cast<std::size_t>(runtime.RowEnd()-runtime.RowBegin()),
			runtime.OwnedStencilCount(),runtime.RequiredStateRows()};
	}
	std::string FixedGeometryIdentity() const
	{
		Sha256 hash;
		immersed_transient_detail::AppendString(hash,"FixedImmersedCaseGeometry/v1");
		immersed_transient_detail::AppendString(hash,SurfaceHash());
		for (int d = 0; d < 3; ++d) {
			hash.AppendNormalizedDouble(Grid().lower_m[d]); hash.AppendNormalizedDouble(Grid().upper_m[d]);
			hash.AppendLittleEndian64(Grid().cells[d]);
		}
		hash.AppendLittleEndian64(volume_->Options().max_depth);
		hash.AppendLittleEndian64(volume_->Options().empty_rule_rescue_max_depth);
		hash.AppendLittleEndian64(static_cast<std::uint64_t>(volume_->StorageMode()));
		return hash.Hex();
	}
	CoupledDomainRuntime& Adapter() const
	{
		if (!adapter_) throw std::logic_error("immersed case runtime has not been initialized");
		return *adapter_;
	}
	static const config_detail::JsonValue& Required(const std::map<std::string, config_detail::JsonValue>& object,
		const std::string& key)
	{
		const auto value = config_detail::Find(object, key);
		if (!value) throw std::runtime_error("immersed_geometry.json: missing required key '"+key+"'");
		return *value;
	}
	static std::filesystem::path Contained(const std::filesystem::path& root, const std::filesystem::path& relative)
	{
		if (relative.empty() || relative.is_absolute()) throw std::runtime_error("immersed case asset path must be relative");
		std::error_code error; const auto path = std::filesystem::canonical(root/relative, error);
		const auto mismatch = std::mismatch(root.begin(), root.end(), path.begin(), path.end());
		if (error || mismatch.first != root.end() || !std::filesystem::is_regular_file(path))
			throw std::runtime_error("immersed case asset must be a contained regular file");
		return path;
	}
	static config_detail::JsonValue ParseJson(const std::filesystem::path& path)
	{
		std::ifstream input(path); if (!input) throw std::runtime_error("cannot open immersed geometry configuration");
		return config_detail::JsonParser(ReadCheckedText(input)).Parse();
	}
	static std::string OptionalString(const std::map<std::string, config_detail::JsonValue>& object,
		const std::string& key, const std::string& fallback)
	{ const auto value = config_detail::Find(object, key); return value ? config_detail::RequireString(*value, "immersed_geometry.json."+key) : fallback; }
	static std::array<double, 3> Vector3(const config_detail::JsonValue& value, const std::string& context)
	{
		const auto& a = config_detail::RequireArray(value, context); if (a.size() != 3) throw std::runtime_error("immersed_geometry.json: "+context+" must have three values");
		return {{config_detail::RequireNumber(a[0], context), config_detail::RequireNumber(a[1], context), config_detail::RequireNumber(a[2], context)}};
	}
	static CubicCartesianGridSpec ParseGrid(const config_detail::JsonValue& value)
	{
		const auto& o = config_detail::RequireObject(value, "immersed_geometry.json.grid");
		config_detail::RequireKnownKeys(o, {"lower_m", "upper_m", "cells"}, "immersed_geometry.json.grid");
		CubicCartesianGridSpec r; r.lower_m = Vector3(Required(o, "lower_m"), "grid.lower_m"); r.upper_m = Vector3(Required(o, "upper_m"), "grid.upper_m");
		const auto& a = config_detail::RequireArray(Required(o, "cells"), "immersed_geometry.json.grid.cells");
		if (a.size() != 3) throw std::runtime_error("immersed_geometry.json: grid.cells must have three values");
		for (std::size_t i = 0; i < 3; ++i)
			r.cells[i] = PositiveInteger<std::uint32_t>(a[i], "grid.cells");
		return r;
	}
	template <class Integer> static Integer PositiveInteger(const config_detail::JsonValue& value,
		const std::string& context)
	{
		static_assert(std::numeric_limits<Integer>::is_integer, "integer target required");
		const double number = config_detail::RequireNumber(value, context);
		const double maximum = static_cast<double>(std::numeric_limits<Integer>::max());
		const bool maximum_is_exact = std::numeric_limits<Integer>::digits
			<= std::numeric_limits<double>::digits;
		if (!(number >= 1.0) || std::floor(number) != number
			|| (maximum_is_exact ? number > maximum : !(number < maximum)))
			throw std::runtime_error("immersed_geometry.json: "+context+" must be a positive integer in range");
		return static_cast<Integer>(number);
	}
	template <class Integer> static Integer NonnegativeInteger(const config_detail::JsonValue& value,
		const std::string& context)
	{
		const double number = config_detail::RequireNumber(value, context);
		const double maximum = static_cast<double>(std::numeric_limits<Integer>::max());
		const bool maximum_is_exact = std::numeric_limits<Integer>::digits
			<= std::numeric_limits<double>::digits;
		if (!(number >= 0.0) || std::floor(number) != number
			|| (maximum_is_exact ? number > maximum : !(number < maximum)))
			throw std::runtime_error("immersed_geometry.json: "+context+" must be a nonnegative integer in range");
		return static_cast<Integer>(number);
	}
	template <class Integer> static Integer Limit(const std::map<std::string, config_detail::JsonValue>& o,
		const std::string& key, Integer fallback, const std::string& context)
	{ const auto v = config_detail::Find(o, key); return v ? PositiveInteger<Integer>(*v, context+"."+key) : fallback; }
	static OctreeCutQuadratureOptions ParseVolume(const config_detail::JsonValue& value)
	{
		const auto& o = config_detail::RequireObject(value, "immersed_geometry.json.volume_quadrature"); config_detail::RequireKnownKeys(o, {"storage", "max_depth", "max_nodes", "max_leaves", "max_points", "max_records", "max_retained_bytes", "max_logical_points"}, "immersed_geometry.json.volume_quadrature");
		OctreeCutQuadratureOptions r; r.max_depth = Limit(o, "max_depth", r.max_depth, "volume_quadrature"); r.max_nodes = Limit(o,"max_nodes",r.max_nodes,"volume_quadrature"); r.max_leaves=Limit(o,"max_leaves",r.max_leaves,"volume_quadrature"); r.max_points=Limit(o,"max_points",r.max_points,"volume_quadrature"); r.max_records=Limit(o,"max_records",r.max_records,"volume_quadrature"); r.max_retained_bytes=Limit(o,"max_retained_bytes",r.max_retained_bytes,"volume_quadrature"); r.max_logical_points=Limit(o,"max_logical_points",r.max_logical_points,"volume_quadrature"); return r;
	}
	static CutCellVolumeQuadratureStorageMode ParseStorage(const config_detail::JsonValue& value)
	{ const auto& o=config_detail::RequireObject(value,"immersed_geometry.json.volume_quadrature"); const auto s=OptionalString(o,"storage","compact"); if(s=="compact")return CutCellVolumeQuadratureStorageMode::Compact; if(s=="expanded")return CutCellVolumeQuadratureStorageMode::Expanded; throw std::runtime_error("immersed_geometry.json: volume_quadrature.storage must be compact or expanded"); }
	static ImmersedSurfaceQuadratureOptions ParseSurface(const config_detail::JsonValue& value)
	{
		const auto& o=config_detail::RequireObject(value,"immersed_geometry.json.surface_quadrature"); config_detail::RequireKnownKeys(o,{"max_candidates","max_fragments","max_points","max_exact_limbs"},"immersed_geometry.json.surface_quadrature"); ImmersedSurfaceQuadratureOptions r; r.max_candidates=Limit(o,"max_candidates",r.max_candidates,"surface_quadrature"); r.max_fragments=Limit(o,"max_fragments",r.max_fragments,"surface_quadrature"); r.max_points=Limit(o,"max_points",r.max_points,"surface_quadrature"); r.max_exact_limbs=Limit(o,"max_exact_limbs",r.max_exact_limbs,"surface_quadrature"); return r;
	}
	static CutCellGhostPenaltyOptions ParseGhost(const config_detail::JsonValue& value)
	{
		const auto& o=config_detail::RequireObject(value,"immersed_geometry.json.ghost_penalty"); config_detail::RequireKnownKeys(o,{"gamma_u","gamma_p","max_faces","max_quadrature_points","max_trace_entries"},"immersed_geometry.json.ghost_penalty"); CutCellGhostPenaltyOptions r; if(const auto v=config_detail::Find(o,"gamma_u"))r.gamma_u=config_detail::RequireNumber(*v,"ghost_penalty.gamma_u"); if(const auto v=config_detail::Find(o,"gamma_p"))r.gamma_p=config_detail::RequireNumber(*v,"ghost_penalty.gamma_p"); r.max_faces=Limit(o,"max_faces",r.max_faces,"ghost_penalty");r.max_quadrature_points=Limit(o,"max_quadrature_points",r.max_quadrature_points,"ghost_penalty");r.max_trace_entries=Limit(o,"max_trace_entries",r.max_trace_entries,"ghost_penalty");return r;
	}
	static std::vector<int> ParseLabels(const config_detail::JsonValue& value)
	{
		const auto& a=config_detail::RequireArray(value,"immersed_geometry.json.wall_labels"); std::vector<int> r; std::set<int> seen; for(const auto& v:a){const int x=NonnegativeInteger<int>(v,"wall_labels");if(!seen.insert(x).second)throw std::runtime_error("immersed_geometry.json: wall_labels must be unique nonnegative integers");r.push_back(x);} return r;
	}
	static double Number(const std::map<std::string, config_detail::JsonValue>& o, const std::string& key, double fallback)
	{ const auto v=config_detail::Find(o,key); return v ? config_detail::RequireNumber(*v,"immersed_geometry.json.runtime."+key) : fallback; }
	template<class Options> static Options ParseRuntime(const config_detail::JsonValue& value,
		const SimulationConfiguration& simulation, const std::vector<int>& wall_labels)
	{
		const auto& o=config_detail::RequireObject(value,"immersed_geometry.json.runtime"); const std::set<std::string> common_keys{"wall_gamma0","nonlinear_maximum_iterations","ksp_maximum_iterations","ksp_relative_tolerance","nonlinear_relative_tolerance","nonlinear_absolute_tolerance","flow_controller_relative_tolerance","flow_controller_absolute_tolerance_m3_s","flow_controller_reference_flow_m3_s","minimum_damping","lu_pivot_shift","ports"};
		auto keys = common_keys;
		if constexpr (std::is_same<Options,ImmersedTransientFlowOptions>::value) keys.insert("wall_inertial_gamma0");
		config_detail::RequireKnownKeys(o,keys,"immersed_geometry.json.runtime");
		const EquationSystemDefinition* flow = nullptr;
		for (const auto& system : simulation.equation_systems)
			if (system.kind == EquationKind::NavierStokes) flow = &system;
		if (!flow) throw std::runtime_error("immersed flow case has no Navier-Stokes system");
		Options r;
		if constexpr (std::is_same<Options,ImmersedTransientFlowOptions>::value)
			r.wall_inertial_gamma0 = Number(o,"wall_inertial_gamma0",r.wall_inertial_gamma0);
		r.parameters = NavierStokesParameters{flow->density, flow->viscosity, 0.0};
		r.wall_labels=wall_labels; r.wall_gamma0=Number(o,"wall_gamma0",r.wall_gamma0); r.nonlinear_maximum_iterations=Limit(o,"nonlinear_maximum_iterations",r.nonlinear_maximum_iterations,"runtime");r.ksp_maximum_iterations=Limit(o,"ksp_maximum_iterations",r.ksp_maximum_iterations,"runtime");r.ksp_relative_tolerance=Number(o,"ksp_relative_tolerance",r.ksp_relative_tolerance);r.nonlinear_relative_tolerance=Number(o,"nonlinear_relative_tolerance",r.nonlinear_relative_tolerance);r.nonlinear_absolute_tolerance=Number(o,"nonlinear_absolute_tolerance",r.nonlinear_absolute_tolerance);r.flow_controller_relative_tolerance=Number(o,"flow_controller_relative_tolerance",r.flow_controller_relative_tolerance);r.flow_controller_absolute_tolerance_m3_s=Number(o,"flow_controller_absolute_tolerance_m3_s",r.flow_controller_absolute_tolerance_m3_s);r.flow_controller_reference_flow_m3_s=Number(o,"flow_controller_reference_flow_m3_s",r.flow_controller_reference_flow_m3_s);r.minimum_damping=Number(o,"minimum_damping",r.minimum_damping);r.lu_pivot_shift=Number(o,"lu_pivot_shift",r.lu_pivot_shift);
		const auto& a=config_detail::RequireArray(Required(o,"ports"),"immersed_geometry.json.runtime.ports"); for(const auto& v:a){const auto& p=config_detail::RequireObject(v,"runtime.ports[]");config_detail::RequireKnownKeys(p,{"id","boundary_label","control_mode","value"},"runtime.ports[]");ImmersedFlowPortDefinition d;d.id=config_detail::RequireString(Required(p,"id"),"runtime.ports[].id");d.boundary_label=PositiveInteger<int>(Required(p,"boundary_label"),"runtime.ports[].boundary_label");const auto mode=config_detail::RequireString(Required(p,"control_mode"),"runtime.ports[].control_mode");if(mode=="pressure")d.control_mode=ImmersedFlowPortControlMode::Pressure;else if(mode=="mean_normal_traction")d.control_mode=ImmersedFlowPortControlMode::MeanNormalTraction;else if(mode=="flow_rate")d.control_mode=ImmersedFlowPortControlMode::FlowRate;else throw std::runtime_error("immersed_geometry.json: total pressure and unknown port controls are unsupported");d.value=config_detail::RequireNumber(Required(p,"value"),"runtime.ports[].value");r.ports.push_back(std::move(d));} return r;
	}
	static bool ValidateSimulation(const SimulationConfiguration& c)
	{
		int flows=0, transports=0; const EquationSystemDefinition* flow=nullptr; for(const auto& s:c.equation_systems){if(s.kind==EquationKind::NavierStokes){++flows;flow=&s;}if(s.kind==EquationKind::LinearTransport)++transports;} if(c.dimension!="3d"||c.has_mesh||flows!=1||transports||!flow||flow->unknowns.size()!=2||c.physiology.enabled||c.coupling.mode!=SimulationScopeMode::FlowOnly||(flow->time_integration!="steady"&&flow->time_integration!="backward_euler"))throw std::runtime_error("immersed flow case requires one steady or backward-Euler flow-only 3D Navier-Stokes system without a body-fitted mesh or transport"); for(const auto& f:c.fields)if(f.kind==FieldKind::Scalar)throw std::runtime_error("immersed flow case rejects transport/species fields"); if(!c.velocity_sources.empty())throw std::runtime_error("immersed flow case rejects body-fitted reference velocity profiles"); for(const auto& b:c.boundaries)for(const auto& condition:b.conditions)if(condition.pressure_gauge||condition.kind==FieldBoundaryKind::Resistance||condition.kind==FieldBoundaryKind::WindkesselRC||condition.kind==FieldBoundaryKind::WindkesselRCR||condition.kind==FieldBoundaryKind::PressureTraction||!condition.profile.empty())throw std::runtime_error("immersed flow case rejects body-fitted boundary, outlet-model, gauge, and reference-profile assumptions");
		return flow->time_integration == "backward_euler";
	}
	static void ValidateLabelPartition(const ClosedTriangulatedSurface& surface, const std::vector<int>& walls,
		const std::vector<CouplingPort>& ports, const std::vector<ImmersedFlowPortDefinition>& runtime_ports)
	{
		std::set<int> labels; for(const auto& item:surface.Diagnostics().boundary_id_histogram) labels.insert(static_cast<int>(item.first)); std::set<int> expected_walls(walls.begin(),walls.end()), expected_ports; std::map<int,const CouplingPort*> graph; for(const auto& port:ports){const int label=ParseThreeDImmersedFlowBoundaryLabel(port);if(label==0)throw std::runtime_error("immersed graph ports must use positive labels");if(!graph.emplace(label,&port).second)throw std::runtime_error("immersed graph ports have duplicate boundary labels");expected_ports.insert(label);} for(const auto& p:runtime_ports){const auto it=graph.find(p.boundary_label);if(it==graph.end()||p.id!=it->second->id)throw std::runtime_error("immersed runtime port id/label does not match graph metadata"); const auto required=it->second->requires;const auto want=p.control_mode==ImmersedFlowPortControlMode::FlowRate?PortQuantity::FlowRate:p.control_mode==ImmersedFlowPortControlMode::Pressure?PortQuantity::MeanPressure:PortQuantity::MeanNormalTraction;if(required!=std::set<PortQuantity>{want})throw std::runtime_error("immersed runtime port control mode does not match graph metadata");} if(runtime_ports.size()!=ports.size())throw std::runtime_error("immersed runtime port catalog does not exactly match graph metadata"); for(int l:expected_walls)if(expected_ports.count(l))throw std::runtime_error("immersed wall and graph port labels overlap"); std::set<int> all=expected_walls;all.insert(expected_ports.begin(),expected_ports.end());if(all!=labels)throw std::runtime_error("immersed wall and graph port labels must exactly partition VTP boundary labels");
	}

	ImmersedFlowCase() = default;
	std::filesystem::path case_directory_; SimulationConfiguration configuration_;
	std::unique_ptr<CartesianDomainClassification> classification_;
	std::unique_ptr<CutCellVolumeQuadratureCatalog> volume_;
	std::unique_ptr<ImmersedSurfaceQuadratureCatalog> surface_;
	std::unique_ptr<CutCellGhostPenaltyCatalog> ghost_;
	ImmersedStaticFlowOptions runtime_options_;
	ImmersedTransientFlowOptions transient_options_;
	bool transient_ = false;
	std::string geometry_identity_;
	std::string domain_id_;
	std::vector<CouplingPort> graph_ports_;
	std::unique_ptr<ImmersedStaticFlowRuntime> runtime_;
	std::unique_ptr<ImmersedStaticDistributedRuntime> distributed_runtime_;
	std::unique_ptr<ImmersedTransientDistributedRuntime> transient_runtime_;
	NavierStokesParameters runtime_parameters_;
	std::unique_ptr<CoupledDomainRuntime> adapter_;
};

} // namespace iga

#endif
