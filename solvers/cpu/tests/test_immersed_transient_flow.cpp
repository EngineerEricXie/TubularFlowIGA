#include "ImmersedTransientFlowRuntime.hpp"
#include "ImmersedVelocityExtension.hpp"
#include "PrescribedSurfaceMotion.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

iga::RawSurfaceTriangle Face(std::int64_t a, std::int64_t b, std::int64_t c, std::uint32_t label)
{
	iga::RawSurfaceTriangle r; r.indices={{a,b,c}}; r.boundary_id=label; return r;
}

iga::RawSurfaceSoup Cube(double lo=.1, double hi=.9)
{
	iga::RawSurfaceSoup r;
	r.vertices={{{{lo,lo,lo}},{{hi,lo,lo}},{{hi,hi,lo}},{{lo,hi,lo}},{{lo,lo,hi}},{{hi,lo,hi}},{{hi,hi,hi}},{{lo,hi,hi}}}};
	// Two flow-controlled patches, walls everywhere else: no pressure-like port
	// means the production gauge row must be appended.
	r.triangles={Face(0,2,1,7),Face(0,3,2,7),Face(4,5,6,7),Face(4,6,7,7),
		Face(0,1,5,7),Face(0,5,4,7),Face(1,2,6,2),Face(1,6,5,2),
		Face(2,3,7,7),Face(2,7,6,7),Face(3,0,4,1),Face(3,4,7,1)};
	return r;
}

template<class F> void Reject(F&& f) { bool rejected=false; try { f(); } catch(const std::exception&) { rejected=true; } assert(rejected); }

template<class F> void RejectWithMessage(F&& f, const char* message)
{
	bool rejected=false;
	try { f(); }
	catch(const std::exception& error) { rejected=std::string(error.what())==message; }
	assert(rejected);
}

bool SamePorts(const std::vector<iga::ImmersedTransientFlowDiagnostics::Port>& a, const std::vector<iga::ImmersedTransientFlowDiagnostics::Port>& b)
{
	if(a.size()!=b.size()) return false;
	for(std::size_t i=0;i<a.size();++i) {
		if(a[i].id!=b[i].id || a[i].boundary_label!=b[i].boundary_label || a[i].control_mode!=b[i].control_mode || a[i].target!=b[i].target || a[i].multiplier!=b[i].multiplier || a[i].controller_error!=b[i].controller_error || a[i].multiplier_row!=b[i].multiplier_row) return false;
		const auto& x=a[i].measurement; const auto& y=b[i].measurement;
		if(x.area_m2!=y.area_m2 || x.outward_flow_m3_s!=y.outward_flow_m3_s || x.mean_pressure_pa!=y.mean_pressure_pa || x.mean_normal_traction_pa!=y.mean_normal_traction_pa || x.mean_velocity_squared_m2_s2!=y.mean_velocity_squared_m2_s2) return false;
	}
	return true;
}

template<class T> bool SameBits(const T& a, const T& b)
{
	static_assert(std::is_trivially_copyable<T>::value,"bit comparison requires a trivial value");
	return std::memcmp(&a,&b,sizeof(T))==0;
}

bool SamePortsBitwise(const std::vector<iga::ImmersedTransientFlowDiagnostics::Port>& a, const std::vector<iga::ImmersedTransientFlowDiagnostics::Port>& b)
{
	if(a.size()!=b.size()) return false;
	for(std::size_t i=0;i<a.size();++i) {
		const auto& x=a[i]; const auto& y=b[i];
		if(x.id!=y.id || x.boundary_label!=y.boundary_label || x.control_mode!=y.control_mode || !SameBits(x.target,y.target) || !SameBits(x.multiplier,y.multiplier) || !SameBits(x.controller_error,y.controller_error) || x.multiplier_row!=y.multiplier_row || x.measurement_valid!=y.measurement_valid) return false;
		if(!SameBits(x.measurement.area_m2,y.measurement.area_m2) || !SameBits(x.measurement.outward_flow_m3_s,y.measurement.outward_flow_m3_s) || !SameBits(x.measurement.mean_pressure_pa,y.measurement.mean_pressure_pa) || !SameBits(x.measurement.mean_normal_traction_pa,y.measurement.mean_normal_traction_pa) || !SameBits(x.measurement.mean_velocity_squared_m2_s2,y.measurement.mean_velocity_squared_m2_s2)) return false;
	}
	return true;
}

bool SameNewtonRecords(const std::vector<iga::ImmersedTransientFlowNewtonStep>& a, const std::vector<iga::ImmersedTransientFlowNewtonStep>& b)
{
	if(a.size()!=b.size()) return false;
	for(std::size_t i=0;i<a.size();++i) {
		const auto& x=a[i]; const auto& y=b[i];
		if(x.iteration!=y.iteration || x.ksp_iterations!=y.ksp_iterations || x.ksp_reason!=y.ksp_reason || !SameBits(x.residual_norm,y.residual_norm) || !SameBits(x.update_norm,y.update_norm) || !SameBits(x.linear_relative_residual,y.linear_relative_residual) || !SameBits(x.damping,y.damping)) return false;
	}
	return true;
}

bool SameDiagnosticsPublic(const iga::ImmersedTransientFlowDiagnostics& a, const iga::ImmersedTransientFlowDiagnostics& b)
{
	return a.active_nodes==b.active_nodes && a.physical_dofs==b.physical_dofs && a.total_dofs==b.total_dofs && a.volume_cells==b.volume_cells && a.surface_cells==b.surface_cells && a.ghost_faces==b.ghost_faces
		&& SameBits(a.target_time_s,b.target_time_s) && SameBits(a.dt_s,b.dt_s) && SameBits(a.pressure_measure,b.pressure_measure) && SameBits(a.pressure_gauge_defect,b.pressure_gauge_defect) && SameBits(a.residual_norm,b.residual_norm) && SameBits(a.true_linear_relative_residual,b.true_linear_relative_residual)
		&& a.nonlinear_iterations==b.nonlinear_iterations && a.ksp_iterations==b.ksp_iterations && a.ksp_reason==b.ksp_reason && a.identity_history_nodes==b.identity_history_nodes && a.committed_history_nodes==b.committed_history_nodes && a.extended_history_nodes==b.extended_history_nodes && a.missing_history_nodes==b.missing_history_nodes
		&& a.idle==b.idle && a.trial_active==b.trial_active && a.converged==b.converged && a.prepared==b.prepared && a.committed==b.committed && a.scalar_diagonal_structure_verified==b.scalar_diagonal_structure_verified
		&& a.attempt_count==b.attempt_count && a.abort_count==b.abort_count && a.rollback_count==b.rollback_count && a.prepare_count==b.prepare_count && a.finalize_count==b.finalize_count && a.commit_count==b.commit_count && SameBits(a.last_assembly_seconds,b.last_assembly_seconds) && SameBits(a.last_linear_solve_seconds,b.last_linear_solve_seconds)
		&& a.geometry_identity_sha256==b.geometry_identity_sha256 && a.layout_hash_sha256==b.layout_hash_sha256 && a.committed_state_hash_sha256==b.committed_state_hash_sha256 && a.trial_state_hash_sha256==b.trial_state_hash_sha256 && a.history_hash_sha256==b.history_hash_sha256 && a.moving_map_identity_sha256==b.moving_map_identity_sha256 && a.input_hash_sha256==b.input_hash_sha256 && a.solved_state_hash_sha256==b.solved_state_hash_sha256 && a.prepared_hash_sha256==b.prepared_hash_sha256 && a.attempt_hash_sha256==b.attempt_hash_sha256
		&& a.attempt_assembly_count==b.attempt_assembly_count && SamePortsBitwise(a.ports,b.ports) && SameNewtonRecords(a.newton_steps,b.newton_steps);
}

iga::ImmersedTransientFlowOptions Options()
{
	iga::ImmersedTransientFlowOptions options;
	options.wall_labels={7}; options.parameters={1.0,1.0,1.0};
	options.ports={{"left",1,iga::ImmersedFlowPortControlMode::FlowRate,0.0},{"right",2,iga::ImmersedFlowPortControlMode::FlowRate,0.0}};
	options.nonlinear_maximum_iterations=16; options.nonlinear_absolute_tolerance=1e-11; options.nonlinear_relative_tolerance=1e-9; options.ksp_relative_tolerance=1e-16; options.lu_pivot_shift=1e-20; options.flow_controller_reference_flow_m3_s=1.0;
	return options;
}

iga::ImmersedTransientFlowOptions PressureLikeOptions(iga::ImmersedFlowPortControlMode mode, double value)
{
	auto options=Options(); options.wall_labels={2,7};
	options.ports={{"port",1,mode,value}};
	return options;
}

std::vector<std::array<double,4>> NonconstantFields(const iga::ImmersedActiveLayout& layout);
std::vector<PetscScalar> NonconstantTrial(const iga::ImmersedTransientFlowRuntime& runtime);

bool WallPenaltyCleared(const iga::ImmersedNitscheWallDiagnostics& value)
{
	return value.by_boundary_id.empty() && value.fraction_lower==0.0 && value.fraction_estimate==0.0
		&& value.fraction_upper==0.0 && value.minimum_h_n_m==0.0 && value.maximum_h_n_m==0.0
		&& value.minimum_eta==0.0 && value.maximum_eta==0.0 && value.minimum_eta_mu==0.0
		&& value.maximum_eta_mu==0.0 && value.minimum_eta_t==0.0 && value.maximum_eta_t==0.0
		&& value.maximum_gap_norm==0.0 && !value.ghost_covered_policy;
}

bool SameWallPenaltyBitwise(const iga::ImmersedNitscheWallDiagnostics& a,
	const iga::ImmersedNitscheWallDiagnostics& b)
{
	if(a.by_boundary_id.size()!=b.by_boundary_id.size() || a.ghost_covered_policy!=b.ghost_covered_policy) return false;
	for(const auto& item:a.by_boundary_id) {
		const auto found=b.by_boundary_id.find(item.first);
		if(found==b.by_boundary_id.end()) return false;
		const auto& x=item.second; const auto& y=found->second;
		if(x.selected_points!=y.selected_points || x.skipped_points!=y.skipped_points
			|| !SameBits(x.selected_area_m2,y.selected_area_m2) || !SameBits(x.skipped_area_m2,y.skipped_area_m2)) return false;
	}
	return SameBits(a.fraction_lower,b.fraction_lower) && SameBits(a.fraction_estimate,b.fraction_estimate)
		&& SameBits(a.fraction_upper,b.fraction_upper) && SameBits(a.minimum_h_n_m,b.minimum_h_n_m)
		&& SameBits(a.maximum_h_n_m,b.maximum_h_n_m) && SameBits(a.minimum_eta,b.minimum_eta)
		&& SameBits(a.maximum_eta,b.maximum_eta) && SameBits(a.maximum_eta_h_n_over_mu,b.maximum_eta_h_n_over_mu)
		&& SameBits(a.minimum_eta_mu,b.minimum_eta_mu) && SameBits(a.maximum_eta_mu,b.maximum_eta_mu)
		&& SameBits(a.minimum_eta_t,b.minimum_eta_t) && SameBits(a.maximum_eta_t,b.maximum_eta_t)
		&& SameBits(a.minimum_eta_mu_h_n_over_mu,b.minimum_eta_mu_h_n_over_mu)
		&& SameBits(a.maximum_eta_mu_h_n_over_mu,b.maximum_eta_mu_h_n_over_mu)
		&& SameBits(a.minimum_eta_t_dt_over_rho_h_n,b.minimum_eta_t_dt_over_rho_h_n)
		&& SameBits(a.maximum_eta_t_dt_over_rho_h_n,b.maximum_eta_t_dt_over_rho_h_n)
		&& SameBits(a.minimum_eta_mu_fraction,b.minimum_eta_mu_fraction)
		&& SameBits(a.maximum_eta_mu_fraction,b.maximum_eta_mu_fraction)
		&& SameBits(a.minimum_eta_t_fraction,b.minimum_eta_t_fraction)
		&& SameBits(a.maximum_eta_t_fraction,b.maximum_eta_t_fraction)
		&& SameBits(a.maximum_gap_norm,b.maximum_gap_norm);
}

std::array<double,3> ExpectedWallReferenceFractionSums(const iga::MovingCutGeometry& geometry,
	int wall_label)
{
	std::array<double,3> result{{0.0,0.0,0.0}};
	for(std::uint64_t cell=0;cell<geometry.Domain().Cells().size();++cell) {
		if(geometry.Domain().Cells()[cell].classification!=iga::CellClassification::Cut) continue;
		const auto& rule=geometry.Surface().UsableRule(geometry.Domain(),cell);
		if(!std::any_of(rule.Points().begin(),rule.Points().end(),[wall_label](const auto& point) { return point.boundary_id==wall_label; })) continue;
		const auto& diagnostics=geometry.Volume().Cell(cell).diagnostics;
		result[0]+=diagnostics.lower_reference_volume;
		result[1]+=diagnostics.estimated_reference_volume;
		result[2]+=diagnostics.upper_reference_volume;
	}
	return result;
}

void CheckTransientNitscheDiagnostics(const iga::MovingCutGeometry& geometry)
{
	auto options=Options(); options.wall_inertial_gamma0=1.0;
	iga::ImmersedTransientFlowRuntime runtime(geometry,options);
	runtime.SetCommittedGlobalState(iga::ImmersedGlobalFlowState(0.0,0,runtime.Layout(),
		NonconstantFields(runtime.Layout()),{0.0,0.0},true,0.0));
	runtime.BeginTrial(1.0,1,1.0); const auto trial=NonconstantTrial(runtime); runtime.SetTrialState(trial); runtime.Assemble();
	const auto wall=runtime.Diagnostics().wall_penalty;
	assert(runtime.Diagnostics().surface_cells>0 && runtime.Diagnostics().ghost_faces>0);
	assert(wall.ghost_covered_policy && wall.by_boundary_id.at(7).selected_points>0);
	assert(wall.minimum_h_n_m>0.0 && wall.maximum_h_n_m>=wall.minimum_h_n_m);
	const auto expected_fractions=ExpectedWallReferenceFractionSums(geometry,7);
	assert(std::isfinite(wall.fraction_lower) && std::isfinite(wall.fraction_estimate)
		&& std::isfinite(wall.fraction_upper) && wall.fraction_lower>=0.0
		&& wall.fraction_lower<=wall.fraction_estimate && wall.fraction_estimate<=wall.fraction_upper);
	assert(SameBits(wall.fraction_lower,expected_fractions[0])
		&& SameBits(wall.fraction_estimate,expected_fractions[1])
		&& SameBits(wall.fraction_upper,expected_fractions[2]));
	assert(wall.minimum_eta>0.0 && wall.maximum_eta>=wall.minimum_eta
		&& wall.minimum_eta_mu>0.0 && wall.maximum_eta_mu>=wall.minimum_eta_mu
		&& wall.minimum_eta_t>0.0 && wall.maximum_eta_t>=wall.minimum_eta_t);
	assert(std::abs(wall.minimum_eta_t_dt_over_rho_h_n-16.0)<2e-12
		&& std::abs(wall.maximum_eta_t_dt_over_rho_h_n-16.0)<2e-12);
	runtime.Assemble(); assert(SameWallPenaltyBitwise(wall,runtime.Diagnostics().wall_penalty));
	runtime.Rollback(); assert(WallPenaltyCleared(runtime.Diagnostics().wall_penalty));
	runtime.SetTrialState(trial);
	runtime.Assemble(); assert(SameWallPenaltyBitwise(wall,runtime.Diagnostics().wall_penalty));
	runtime.AbortTrial(); assert(runtime.Diagnostics().idle && WallPenaltyCleared(runtime.Diagnostics().wall_penalty));
}

std::vector<std::array<double,4>> NonconstantFields(const iga::ImmersedActiveLayout& layout)
{
	std::vector<std::array<double,4>> fields(layout.NodeIds().size());
	for(std::size_t i=0;i<fields.size();++i) for(int q=0;q<4;++q) fields[i][q]=.05*static_cast<double>(static_cast<int>((13*i+5*q)%17)-8);
	return fields;
}

std::vector<PetscScalar> NonconstantTrial(const iga::ImmersedTransientFlowRuntime& runtime)
{
	std::vector<PetscScalar> result(runtime.Layout().Rows());
	for(std::size_t i=0;i<result.size();++i) result[i]=.003*static_cast<double>(static_cast<int>((7*i+3)%19)-9);
	result[static_cast<std::size_t>(runtime.PortMultiplierDof("left"))]=.017; result[static_cast<std::size_t>(runtime.PortMultiplierDof("right"))]=-.011; result[static_cast<std::size_t>(runtime.GaugeDof())]=.013;
	return result;
}

std::vector<PetscScalar> Direction(const iga::ImmersedTransientFlowRuntime& runtime)
{
	std::vector<PetscScalar> result(runtime.Layout().Rows());
	for(std::size_t i=0;i<result.size();++i) result[i]=.02*static_cast<double>(static_cast<int>((11*i+1)%23)-11);
	result[static_cast<std::size_t>(runtime.PortMultiplierDof("left"))]=.031; result[static_cast<std::size_t>(runtime.PortMultiplierDof("right"))]=-.029; result[static_cast<std::size_t>(runtime.GaugeDof())]=.037;
	return result;
}

std::vector<PetscScalar> GenericTrial(const iga::ImmersedTransientFlowRuntime& runtime)
{
	std::vector<PetscScalar> result(runtime.Layout().Rows());
	for(std::size_t i=0;i<result.size();++i) result[i]=.004*static_cast<double>(static_cast<int>((17*i+9)%29)-14);
	return result;
}

std::vector<PetscScalar> ExpectedPressureLikePortLoad(const iga::MovingCutGeometry& geometry,
	const iga::ImmersedTransientFlowRuntime& runtime, int label, iga::ImmersedFlowPortControlMode mode,
	double value, const std::vector<PetscScalar>& state)
{
	std::vector<PetscScalar> result(runtime.Layout().Rows());
	for(std::uint64_t cell=0;cell<geometry.Surface().Cells().size();++cell) {
		const auto& rule=geometry.Surface().UsableRule(geometry.Domain(),cell);
		if(!std::any_of(rule.Points().begin(),rule.Points().end(),[&](const auto& point) { return point.boundary_id==label; })) continue;
		const auto element=geometry.Domain().Background().MaterializeElement(cell);
		std::vector<std::array<double,4>> nodal(element.connectivity.size());
		for(std::size_t a=0;a<nodal.size();++a) for(int q=0;q<4;++q) nodal[a][q]=PetscRealPart(state[static_cast<std::size_t>(runtime.Dof(element.connectivity[a],q))]);
		const auto local=iga::BuildImmersedFlowPortElement(element,rule,label,mode,value,nodal);
		for(std::size_t a=0;a<element.connectivity.size();++a) for(int q=0;q<4;++q) result[static_cast<std::size_t>(runtime.Dof(element.connectivity[a],q))]+=local.negative_residual[4*a+q];
	}
	return result;
}

iga::ImmersedFlowPortMeasurement ExpectedPortMeasurement(const iga::MovingCutGeometry& geometry,
	const iga::ImmersedTransientFlowRuntime& runtime, int label, const std::vector<PetscScalar>& state)
{
	iga::ImmersedFlowPortMeasurement result;
	for(std::uint64_t cell=0;cell<geometry.Surface().Cells().size();++cell) {
		const auto& rule=geometry.Surface().UsableRule(geometry.Domain(),cell);
		if(!std::any_of(rule.Points().begin(),rule.Points().end(),[&](const auto& point) { return point.boundary_id==label; })) continue;
		const auto element=geometry.Domain().Background().MaterializeElement(cell);
		std::vector<std::array<double,4>> nodal(element.connectivity.size());
		for(std::size_t a=0;a<nodal.size();++a) for(int q=0;q<4;++q) nodal[a][q]=PetscRealPart(state[static_cast<std::size_t>(runtime.Dof(element.connectivity[a],q))]);
		const auto local=iga::MeasureImmersedFlowPortElement(element,rule,label,nodal,1.0);
		result.area_m2+=local.area_m2; result.outward_flow_m3_s+=local.outward_flow_m3_s;
		result.mean_pressure_pa+=local.mean_pressure_pa*local.area_m2;
		result.mean_normal_traction_pa+=local.mean_normal_traction_pa*local.area_m2;
		result.mean_velocity_squared_m2_s2+=local.mean_velocity_squared_m2_s2*local.area_m2;
	}
	result.mean_pressure_pa/=result.area_m2; result.mean_normal_traction_pa/=result.area_m2; result.mean_velocity_squared_m2_s2/=result.area_m2;
	return result;
}

void CheckPressureLikePortBranch(const iga::MovingCutGeometry& geometry, iga::ImmersedFlowPortControlMode mode, double value)
{
	auto loaded_options=PressureLikeOptions(mode,value);
	iga::ImmersedTransientFlowRuntime loaded(geometry,loaded_options);
	assert(!loaded.Layout().HasGaugeRow() && loaded.Layout().PortIds().empty());
	loaded.SetCommittedGlobalState(iga::ImmersedGlobalFlowState(0.0,0,loaded.Layout(),NonconstantFields(loaded.Layout()),{},false,0.0));
	loaded.BeginTrial(1.0,1,1.0);
	const auto state=GenericTrial(loaded), direction=GenericTrial(loaded);
	loaded.SetTrialState(state); loaded.Assemble();
	const auto load=loaded.AssembledNegativeResidual(), jl=loaded.AssembledJacobianAction(direction);
	const auto expected=ExpectedPressureLikePortLoad(geometry,loaded,1,mode,value,state);
	const auto measurement=ExpectedPortMeasurement(geometry,loaded,1,state);
	assert(loaded.Diagnostics().ports.front().multiplier_row<0 && measurement.area_m2>0 && std::isfinite(measurement.outward_flow_m3_s) && std::isfinite(measurement.mean_pressure_pa) && std::isfinite(measurement.mean_normal_traction_pa));
	loaded.AbortTrial(); loaded.SetPortControlValue("port",0.0); loaded.BeginTrial(1.0,1,1.0); loaded.SetTrialState(state); loaded.Assemble();
	const auto baseline=loaded.AssembledNegativeResidual(), jz=loaded.AssembledJacobianAction(direction);
	double expected_norm=0.0;
	for(std::size_t i=0;i<load.size();++i) { assert(std::isfinite(PetscRealPart(load[i]))); assert(std::abs(PetscRealPart(load[i]-baseline[i]-expected[i]))<=2e-13); expected_norm+=PetscRealPart(expected[i])*PetscRealPart(expected[i]); }
	assert(expected_norm>0.0);
	for(std::size_t i=0;i<jl.size();++i) assert(std::abs(PetscRealPart(jl[i]-jz[i]))<=2e-13);
	double x_load=0.0; for(std::size_t i=0;i<expected.size();i+=4) x_load+=PetscRealPart(expected[i]);
	// Label 1 is the cube's outward -x face.  The production element defines a
	// pressure load as -p*n and a mean-normal-traction load as +t*n.
	assert(mode==iga::ImmersedFlowPortControlMode::Pressure ? x_load>0.0 : x_load<0.0);
	loaded.AbortTrial();
}

void CheckAbortPublication(const iga::MovingCutGeometry& geometry)
{
	iga::ImmersedTransientFlowRuntime runtime(geometry,Options());
	std::vector<std::array<double,4>> zero(runtime.Layout().NodeIds().size());
	runtime.SetCommittedGlobalState(iga::ImmersedGlobalFlowState(0.0,0,runtime.Layout(),zero,{0.0,0.0},true,0.0));
	const auto committed=runtime.CommittedState(); const auto hash=runtime.CommittedGlobalState().HashSha256(); const auto ports=runtime.Diagnostics().ports;
	runtime.BeginTrial(1.0,1,1.0); runtime.SetTrialState(NonconstantTrial(runtime)); runtime.Assemble(); assert(SamePorts(ports,runtime.Diagnostics().ports));
	assert(!runtime.Diagnostics().input_hash_sha256.empty() && !runtime.Diagnostics().history_hash_sha256.empty());
	runtime.AbortTrial(); assert(runtime.TrialState()==committed && runtime.Diagnostics().idle && runtime.CommittedState()==committed && runtime.CommittedGlobalState().HashSha256()==hash && SamePortsBitwise(ports,runtime.Diagnostics().ports) && runtime.Diagnostics().input_hash_sha256.empty() && runtime.Diagnostics().trial_state_hash_sha256.empty() && runtime.Diagnostics().history_hash_sha256.empty() && runtime.Diagnostics().prepared_hash_sha256.empty());
}

void CheckMovingInnerTrial(const iga::CubicCartesianGridSpec& grid, iga::MovingCutGeometryOptions options)
{
	const auto old_soup=Cube(.18,.45), new_soup=Cube(.18,.95);
	iga::PrescribedSurfaceMotion motion({{0.0,old_soup},{1.0,new_soup}});
	auto old_geometry=iga::MovingCutGeometry::Build(grid,motion.Evaluate(0.0,0.0,1.0),options);
	auto target_geometry=iga::MovingCutGeometry::Build(grid,motion.Evaluate(1.0,0.0,1.0),options,old_geometry.get());
	const auto old_layout=iga::ImmersedActiveLayout::Build(old_geometry->Domain(),old_geometry->Volume(),old_geometry->GeometryIdentitySha256(),{1,2},true);
	iga::ImmersedGlobalFlowState old_state(0.0,0,old_layout,NonconstantFields(old_layout),{0.0,0.0},true,0.0);
	iga::ImmersedTransientFlowRuntime runtime(*target_geometry,Options());
	const auto extension=iga::ImmersedVelocityExtension::Build(*old_geometry,old_layout,old_state,*target_geometry,runtime.Layout(),3);
	std::vector<double> old_pressure(old_state.Coefficients().size());
	for(std::size_t i=0;i<old_pressure.size();++i) old_pressure[i]=old_state.Coefficients()[i][3];
	const auto scalar=extension.ExtendScalar(old_pressure);
	std::vector<std::array<double,4>> seed_fields(runtime.Layout().NodeIds().size());
	for(std::size_t i=0;i<seed_fields.size();++i) { for(int c=0;c<3;++c) seed_fields[i][c]=extension.TargetHistory().Velocities()[i][c]; seed_fields[i][3]=scalar.target_values[i]; }
	const std::vector<double> seed_multipliers{{.021,-.034}};
	const iga::ImmersedGlobalFlowState seed(0.0,0,runtime.Layout(),seed_fields,seed_multipliers,true,0.0);
	const auto map=iga::ImmersedMovingTrialMapIdentity::Create(old_state.HashSha256(),old_geometry->GeometryIdentitySha256(),target_geometry->GeometryIdentitySha256(),target_geometry->PublicationIdentitySha256(),old_layout.HashSha256(),runtime.Layout().HashSha256(),extension.HashSha256(),extension.OperatorHashSha256(),extension.TargetHistory().HashSha256(),scalar.hash_sha256,runtime.Layout().PortIds(),seed_multipliers,true,0.0,0,1.0,1,1.0);
	std::vector<PetscScalar> expected_seed(runtime.Layout().Rows(),0.0);
	for(std::size_t i=0;i<seed.Coefficients().size();++i) for(int field=0;field<4;++field) expected_seed[static_cast<std::size_t>(runtime.Dof(runtime.Layout().NodeIds()[i],field))]=seed.Coefficients()[i][field];
	for(std::size_t i=0;i<seed.PortIds().size();++i) expected_seed[static_cast<std::size_t>(runtime.Layout().ControllerRow(seed.PortIds()[i]))]=seed.PortMultipliers()[i];
	expected_seed[static_cast<std::size_t>(runtime.GaugeDof())]=seed.GaugeMultiplier();
	assert(SameBits(seed.GaugeMultiplier(),0.0) && SameBits(PetscRealPart(expected_seed[static_cast<std::size_t>(runtime.GaugeDof())]),0.0));
	const auto before=runtime.Diagnostics(); const auto before_committed=runtime.CommittedState(); const auto before_trial=runtime.TrialState();
	auto bad_history=iga::ImmersedVelocityHistory(0.0,1.0,old_geometry->GeometryIdentitySha256(),"wrong",extension.TargetHistory().NodeIds(),extension.TargetHistory().Velocities(),extension.TargetHistory().Provenance());
	Reject([&] { runtime.BeginMovingTrial(1.0,1,1.0,bad_history,seed,map); });
	auto short_ids=extension.TargetHistory().NodeIds(); auto short_velocities=extension.TargetHistory().Velocities(); auto short_provenance=extension.TargetHistory().Provenance(); short_ids.pop_back(); short_velocities.pop_back(); short_provenance.pop_back();
	const auto short_history=iga::ImmersedVelocityHistory(0.0,1.0,old_geometry->GeometryIdentitySha256(),target_geometry->GeometryIdentitySha256(),short_ids,short_velocities,short_provenance);
	Reject([&] { runtime.BeginMovingTrial(1.0,1,1.0,short_history,seed,map); });
	const auto late_history=iga::ImmersedVelocityHistory(0.0,2.0,old_geometry->GeometryIdentitySha256(),target_geometry->GeometryIdentitySha256(),extension.TargetHistory().NodeIds(),extension.TargetHistory().Velocities(),extension.TargetHistory().Provenance());
	Reject([&] { runtime.BeginMovingTrial(1.0,1,1.0,late_history,seed,map); });
	const iga::ImmersedGlobalFlowState bad_time(1.0,0,runtime.Layout(),seed_fields,seed_multipliers,true,0.0);
	Reject([&] { runtime.BeginMovingTrial(1.0,1,1.0,extension.TargetHistory(),bad_time,map); });
	const iga::ImmersedGlobalFlowState bad_controller(0.0,0,runtime.Layout(),seed_fields,{.1,seed_multipliers[1]},true,0.0);
	Reject([&] { runtime.BeginMovingTrial(1.0,1,1.0,extension.TargetHistory(),bad_controller,map); });
	Reject([&] { runtime.BeginMovingTrial(1.0,1,1.0,extension.TargetHistory(),old_state,map); });
	const iga::ImmersedGlobalFlowState bad_gauge(0.0,0,runtime.Layout(),seed_fields,seed_multipliers,true,-0.0);
	Reject([&] { runtime.BeginMovingTrial(1.0,1,1.0,extension.TargetHistory(),bad_gauge,map); });
	const auto stale=iga::ImmersedMovingTrialMapIdentity::Create(old_state.HashSha256(),old_geometry->GeometryIdentitySha256(),"stale",target_geometry->PublicationIdentitySha256(),old_layout.HashSha256(),runtime.Layout().HashSha256(),extension.HashSha256(),extension.OperatorHashSha256(),extension.TargetHistory().HashSha256(),scalar.hash_sha256,runtime.Layout().PortIds(),seed_multipliers,true,0.0,0,1.0,1,1.0);
	Reject([&] { runtime.BeginMovingTrial(1.0,1,1.0,extension.TargetHistory(),seed,stale); });
	assert(runtime.CommittedState()==before_committed && runtime.TrialState()==before_trial && SameDiagnosticsPublic(before,runtime.Diagnostics()));
	runtime.BeginMovingTrial(1.0,1,1.0,extension.TargetHistory(),seed,map);
	const auto trial_seed=runtime.TrialState(); assert(trial_seed.size()==expected_seed.size());
	for(std::size_t i=0;i<trial_seed.size();++i) assert(SameBits(trial_seed[i],expected_seed[i]));
	assert(SameBits(PetscRealPart(trial_seed[static_cast<std::size_t>(runtime.GaugeDof())]),0.0));
	assert(runtime.Diagnostics().history_hash_sha256==extension.TargetHistory().HashSha256() && runtime.Diagnostics().moving_map_identity_sha256==map.HashSha256());
	assert(runtime.Diagnostics().identity_history_nodes==runtime.Layout().NodeIds().size() && runtime.Diagnostics().extended_history_nodes>0 && runtime.Diagnostics().missing_history_nodes==0);
	const auto frozen=runtime.TrialState(); const auto input=runtime.Diagnostics().input_hash_sha256;
	runtime.Assemble(); assert(runtime.Diagnostics().surface_cells>0);
	runtime.SetTrialState(GenericTrial(runtime)); runtime.Rollback();
	assert(runtime.TrialState()==frozen && runtime.Diagnostics().input_hash_sha256==input && runtime.Diagnostics().history_hash_sha256==extension.TargetHistory().HashSha256() && runtime.Diagnostics().moving_map_identity_sha256==map.HashSha256());
	runtime.AbortTrial();
	auto throwing_options=Options(); throwing_options.body_force=[](const std::array<double,3>&) -> std::array<double,3> { throw std::runtime_error("moving body force"); };
	iga::ImmersedTransientFlowRuntime throwing(*target_geometry,throwing_options); const auto throwing_before=throwing.Diagnostics(); const auto throwing_state=throwing.CommittedState();
	Reject([&] { throwing.BeginMovingTrial(1.0,1,1.0,extension.TargetHistory(),seed,map); });
	assert(throwing.CommittedState()==throwing_state && throwing.TrialState()==throwing_state && SameDiagnosticsPublic(throwing_before,throwing.Diagnostics()));
}

void CheckFrozenBodyForce(const iga::MovingCutGeometry& geometry)
{
	double scale=1.0; auto options=Options();
	options.body_force=[&scale](const std::array<double,3>& x) { return std::array<double,3>{{scale*(1.0+x[0]),-2.0*scale,0.5*scale}}; };
	iga::ImmersedTransientFlowRuntime runtime(geometry,options);
	runtime.SetCommittedGlobalState(iga::ImmersedGlobalFlowState(0.0,0,runtime.Layout(),NonconstantFields(runtime.Layout()),{0.0,0.0},true,0.0));
	runtime.BeginTrial(1.0,1,1.0); const auto frozen_input=runtime.Diagnostics().input_hash_sha256;
	runtime.SetTrialState(GenericTrial(runtime)); runtime.Assemble(); const auto frozen_residual=runtime.AssembledNegativeResidual();
	scale=17.0; runtime.Assemble(); assert(runtime.AssembledNegativeResidual()==frozen_residual && runtime.Diagnostics().input_hash_sha256==frozen_input);
	runtime.Rollback(); runtime.SetTrialState(GenericTrial(runtime)); runtime.Assemble(); assert(runtime.AssembledNegativeResidual()==frozen_residual && runtime.Diagnostics().input_hash_sha256==frozen_input);
	runtime.AbortTrial(); runtime.BeginTrial(1.0,1,1.0); runtime.SetTrialState(GenericTrial(runtime)); runtime.Assemble();
	assert(runtime.Diagnostics().input_hash_sha256!=frozen_input && runtime.AssembledNegativeResidual()!=frozen_residual); runtime.AbortTrial();
}

void CheckBlockReduction(const iga::MovingCutGeometry& geometry)
{
	auto permissive=Options(); permissive.nonlinear_maximum_iterations=1; permissive.nonlinear_absolute_tolerance=1e9; permissive.nonlinear_block_reduction=0.0;
	iga::ImmersedTransientFlowRuntime accepts(geometry,permissive);
	accepts.SetCommittedGlobalState(iga::ImmersedGlobalFlowState(0.0,0,accepts.Layout(),NonconstantFields(accepts.Layout()),{0.0,0.0},true,0.0)); accepts.BeginTrial(1.0,1,1.0); assert(accepts.SolveTrial()); accepts.AbortTrial();
	auto gated=permissive; gated.nonlinear_block_reduction=1e100;
	iga::ImmersedTransientFlowRuntime rejects(geometry,gated);
	rejects.SetCommittedGlobalState(iga::ImmersedGlobalFlowState(0.0,0,rejects.Layout(),NonconstantFields(rejects.Layout()),{0.0,0.0},true,0.0)); rejects.BeginTrial(1.0,1,1.0); Reject([&] { rejects.SolveTrial(); }); assert(rejects.Diagnostics().trial_active && !rejects.Diagnostics().converged); rejects.AbortTrial();
}

void CheckIdleMutationPublication(const iga::MovingCutGeometry& geometry)
{
	auto options=PressureLikeOptions(iga::ImmersedFlowPortControlMode::Pressure,.125);
	iga::ImmersedTransientFlowRuntime runtime(geometry,options);
	runtime.SetCommittedGlobalState(iga::ImmersedGlobalFlowState(0.0,0,runtime.Layout(),NonconstantFields(runtime.Layout()),{},false,0.0));
	runtime.BeginTrial(1.0,1,1.0); assert(runtime.SolveTrial()); runtime.Commit();
	const auto solved_input=runtime.Diagnostics().input_hash_sha256;
	std::vector<std::array<double,4>> zero(runtime.Layout().NodeIds().size());
	runtime.SetCommittedGlobalState(iga::ImmersedGlobalFlowState(0.0,0,runtime.Layout(),zero, {},false,0.0));
	assert(runtime.Diagnostics().idle && runtime.Diagnostics().input_hash_sha256.empty() && runtime.Diagnostics().solved_state_hash_sha256.empty() && runtime.Diagnostics().prepared_hash_sha256.empty() && runtime.Diagnostics().attempt_hash_sha256.empty() && runtime.Diagnostics().newton_steps.empty() && runtime.Diagnostics().attempt_assembly_count==0);
	for(const auto& port:runtime.Diagnostics().ports) assert(!port.measurement_valid && port.measurement.area_m2==0.0 && port.controller_error==0.0);
	runtime.SetPortControlValue("port",.25); assert(runtime.Diagnostics().ports.front().target==.25 && !runtime.Diagnostics().ports.front().measurement_valid);
	runtime.BeginTrial(1.0,1,1.0); assert(!runtime.Diagnostics().input_hash_sha256.empty() && runtime.Diagnostics().input_hash_sha256!=solved_input); runtime.AbortTrial();
}

void CheckSolverConfigurationIdentity(const iga::MovingCutGeometry& geometry)
{
	auto options=Options(); iga::ImmersedTransientFlowRuntime first(geometry,options), second(geometry,options);
	for(auto* runtime:{&first,&second}) runtime->SetCommittedGlobalState(iga::ImmersedGlobalFlowState(0.0,0,runtime->Layout(),NonconstantFields(runtime->Layout()),{0.0,0.0},true,0.0));
	first.BeginTrial(1.0,1,1.0); second.BeginTrial(1.0,1,1.0); assert(first.Diagnostics().input_hash_sha256==second.Diagnostics().input_hash_sha256); first.AbortTrial(); second.AbortTrial();
	options.ksp_relative_tolerance=.5*options.ksp_relative_tolerance; iga::ImmersedTransientFlowRuntime changed(geometry,options);
	changed.SetCommittedGlobalState(iga::ImmersedGlobalFlowState(0.0,0,changed.Layout(),NonconstantFields(changed.Layout()),{0.0,0.0},true,0.0)); changed.BeginTrial(1.0,1,1.0);
	first.BeginTrial(1.0,1,1.0); assert(changed.Diagnostics().input_hash_sha256!=first.Diagnostics().input_hash_sha256); changed.AbortTrial(); first.AbortTrial();
}

void CheckBeginTrialExceptionTransaction(const iga::MovingCutGeometry& geometry)
{
	bool throw_callback=true; auto options=Options();
	options.body_force=[&throw_callback](const std::array<double,3>&) { if(throw_callback) throw std::runtime_error("injected body-force failure"); return std::array<double,3>{{0.0,0.0,0.0}}; };
	iga::ImmersedTransientFlowRuntime runtime(geometry,options);
	runtime.SetCommittedGlobalState(iga::ImmersedGlobalFlowState(0.0,0,runtime.Layout(),NonconstantFields(runtime.Layout()),{0.0,0.0},true,0.0));
	const auto before_state=runtime.CommittedState(); const auto before_trial=runtime.TrialState(); const auto before_global=runtime.CommittedGlobalState().HashSha256(); const auto before=runtime.Diagnostics();
	Reject([&] { runtime.BeginTrial(1.0,1,1.0); });
	assert(runtime.CommittedState()==before_state && runtime.TrialState()==before_trial && runtime.CommittedGlobalState().HashSha256()==before_global && SameDiagnosticsPublic(before,runtime.Diagnostics()));
	throw_callback=false; runtime.BeginTrial(1.0,1,1.0); runtime.AbortTrial();

	auto nonfinite=Options(); nonfinite.body_force=[](const std::array<double,3>&) { return std::array<double,3>{{0.0,std::numeric_limits<double>::infinity(),0.0}}; };
	iga::ImmersedTransientFlowRuntime invalid(geometry,nonfinite);
	invalid.SetCommittedGlobalState(iga::ImmersedGlobalFlowState(0.0,0,invalid.Layout(),NonconstantFields(invalid.Layout()),{0.0,0.0},true,0.0));
	const auto invalid_state=invalid.CommittedState(); const auto invalid_trial=invalid.TrialState(); const auto invalid_global=invalid.CommittedGlobalState().HashSha256(); const auto invalid_before=invalid.Diagnostics();
	Reject([&] { invalid.BeginTrial(1.0,1,1.0); });
	assert(invalid.CommittedState()==invalid_state && invalid.TrialState()==invalid_trial && invalid.CommittedGlobalState().HashSha256()==invalid_global && SameDiagnosticsPublic(invalid_before,invalid.Diagnostics()));
}

std::size_t LogicalVolumePoints(const iga::MovingCutGeometry& geometry)
{
	std::size_t result=0;
	for(std::uint64_t cell=0;cell<geometry.Domain().Cells().size();++cell) {
		const auto classification=geometry.Domain().Cells()[cell].classification;
		if(classification!=iga::CellClassification::Inside&&classification!=iga::CellClassification::Cut) continue;
		if(geometry.Volume().StorageMode()==iga::CutCellVolumeQuadratureStorageMode::Expanded) result+=geometry.Volume().UsableRule(geometry.Domain(),cell).Points().size();
		else result+=iga::CompactCutCellVolumeLogicalPointCount(geometry.Volume().UsableCompactRule(geometry.Domain(),cell));
	}
	return result;
}

bool CompactRuntimeNear(double left, double right, double tolerance)
{
	return std::abs(left-right)<=tolerance*std::max({1.0,std::abs(left),std::abs(right)});
}

void CheckVolumeDiagnosticParity(const iga::CutCellVolumeQuadratureDiagnostics& expanded,
	const iga::CutCellVolumeQuadratureDiagnostics& compact)
{
	for(const auto values : {std::array<double,2>{{expanded.certified_reference_volume,compact.certified_reference_volume}},
		std::array<double,2>{{expanded.estimated_reference_volume,compact.estimated_reference_volume}},
		std::array<double,2>{{expanded.unresolved_reference_volume,compact.unresolved_reference_volume}},
		std::array<double,2>{{expanded.lower_reference_volume,compact.lower_reference_volume}},
		std::array<double,2>{{expanded.upper_reference_volume,compact.upper_reference_volume}},
		std::array<double,2>{{expanded.certified_physical_volume,compact.certified_physical_volume}},
		std::array<double,2>{{expanded.estimated_physical_volume,compact.estimated_physical_volume}},
		std::array<double,2>{{expanded.unresolved_physical_volume,compact.unresolved_physical_volume}},
		std::array<double,2>{{expanded.lower_physical_volume,compact.lower_physical_volume}},
		std::array<double,2>{{expanded.upper_physical_volume,compact.upper_physical_volume}}})
		assert(CompactRuntimeNear(values[0],values[1],1e-12));
}

void CheckCompactRuntime(const iga::MovingCutGeometry& expanded_geometry, const iga::MovingCutGeometry& compact_geometry)
{
	std::size_t expanded_calls=0, compact_calls=0; auto expanded_options=Options(), compact_options=Options();
	expanded_options.body_force=[&expanded_calls](const std::array<double,3>& x) { ++expanded_calls; return std::array<double,3>{{1.0+x[0],-2.0+x[1],.5+x[2]}}; };
	compact_options.body_force=[&compact_calls](const std::array<double,3>& x) { ++compact_calls; return std::array<double,3>{{1.0+x[0],-2.0+x[1],.5+x[2]}}; };
	iga::ImmersedTransientFlowRuntime expanded(expanded_geometry,expanded_options), compact(compact_geometry,compact_options);
	assert(expanded.Layout().NodeIds()==compact.Layout().NodeIds());
	expanded.SetCommittedGlobalState(iga::ImmersedGlobalFlowState(0.0,0,expanded.Layout(),NonconstantFields(expanded.Layout()),{0.0,0.0},true,0.0));
	compact.SetCommittedGlobalState(iga::ImmersedGlobalFlowState(0.0,0,compact.Layout(),NonconstantFields(compact.Layout()),{0.0,0.0},true,0.0));
	const auto expanded_points=LogicalVolumePoints(expanded_geometry), compact_points=LogicalVolumePoints(compact_geometry);
	assert(expanded_points>0 && compact_points>0 && compact_points<expanded_points);
	CheckVolumeDiagnosticParity(expanded_geometry.Volume().Diagnostics(),compact_geometry.Volume().Diagnostics());
	for(std::uint64_t cell=0;cell<expanded_geometry.Domain().Cells().size();++cell)
		CheckVolumeDiagnosticParity(expanded_geometry.Volume().Cell(cell).diagnostics,compact_geometry.Volume().Cell(cell).diagnostics);
	expanded.BeginTrial(1.0,1,1.0); compact.BeginTrial(1.0,1,1.0);
	const auto compact_input=compact.Diagnostics().input_hash_sha256;
	assert(expanded_calls==expanded_points && compact_calls==compact_points && expanded.Diagnostics().input_hash_sha256!=compact_input);
	const auto x=NonconstantTrial(expanded), direction=Direction(expanded); expanded.SetTrialState(x); compact.SetTrialState(x); expanded.Assemble(); compact.Assemble();
	const auto er=expanded.AssembledNegativeResidual(), cr=compact.AssembledNegativeResidual(), ej=expanded.AssembledJacobianAction(direction), cj=compact.AssembledJacobianAction(direction);
	// Coalesced compact certified blocks preserve degree-seven moments, but the
	// nonlinear convection and stabilization terms are higher degree.
	for(std::size_t i=0;i<er.size();++i) {
		assert(CompactRuntimeNear(PetscRealPart(er[i]),PetscRealPart(cr[i]),5e-8));
		assert(CompactRuntimeNear(PetscRealPart(ej[i]),PetscRealPart(cj[i]),3e-7));
	}
	const auto ec=expanded.ConservationDiagnostics(), cc=compact.ConservationDiagnostics();
	assert(std::abs(ec.endpoint_volume_divergence_m3_s-cc.endpoint_volume_divergence_m3_s)<=2e-11 && std::abs(ec.total_surface_outward_flow_m3_s-cc.total_surface_outward_flow_m3_s)<=2e-11);
	assert(ec.surface_flow_by_boundary_label_m3_s==cc.surface_flow_by_boundary_label_m3_s && ec.material_surface_outward_flow_by_boundary_label_m3_s==cc.material_surface_outward_flow_by_boundary_label_m3_s && ec.material_wall_outward_flow_by_boundary_label_m3_s==cc.material_wall_outward_flow_by_boundary_label_m3_s);
	assert(std::abs(ec.total_material_surface_outward_flow_m3_s-cc.total_material_surface_outward_flow_m3_s)<=2e-11 && std::abs(ec.total_material_wall_outward_flow_m3_s-cc.total_material_wall_outward_flow_m3_s)<=2e-11 && std::abs(ec.wall_relative_leakage_m3_s-cc.wall_relative_leakage_m3_s)<=2e-11);
	assert(std::abs(ec.divergence_theorem_defect_m3_s-cc.divergence_theorem_defect_m3_s)<=2e-11 && std::abs(ec.normalized_open_balance-cc.normalized_open_balance)<=2e-11 && std::abs(ec.normalized_wall_leakage-cc.normalized_wall_leakage)<=2e-11);
	assert(std::abs(ec.discrete_moving_wall_continuity_defect_m3_s-cc.discrete_moving_wall_continuity_defect_m3_s)<=2e-11 && std::abs(ec.discrete_moving_wall_continuity_normalization_scale_m3_s-cc.discrete_moving_wall_continuity_normalization_scale_m3_s)<=2e-11 && std::abs(ec.normalized_discrete_moving_wall_continuity_defect-cc.normalized_discrete_moving_wall_continuity_defect)<=2e-11);
	expanded.Assemble(); compact.Assemble(); assert(expanded_calls==expanded_points && compact_calls==compact_points);
	expanded.AbortTrial(); compact.AbortTrial(); compact.BeginTrial(1.0,1,1.0);
	assert(compact.Diagnostics().input_hash_sha256==compact_input && compact_calls==2*compact_points); compact.AbortTrial();
}

struct FirstSolveFormulationGate {
	double initial_residual=std::numeric_limits<double>::quiet_NaN();
	double convergence_threshold=std::numeric_limits<double>::quiet_NaN();
	bool nonempty_record=false, update_count_within_limit=false, converged=false;
	bool residual_converged=false, linear_residual_satisfied=false, accounting_consistent=false;

	bool Passed() const
	{
		return nonempty_record && update_count_within_limit && converged && residual_converged
			&& linear_residual_satisfied && accounting_consistent;
	}
};

FirstSolveFormulationGate EvaluateFirstSolveFormulationGate(
	const iga::ImmersedTransientFlowDiagnostics& diagnostics,
	const iga::ImmersedTransientFlowOptions& options)
{
	FirstSolveFormulationGate gate;
	const auto& steps=diagnostics.newton_steps;
	gate.nonempty_record=!steps.empty();
	gate.update_count_within_limit=steps.size()<=static_cast<std::size_t>(options.nonlinear_maximum_iterations);
	gate.converged=diagnostics.converged;
	if(gate.nonempty_record) {
		gate.initial_residual=steps.front().residual_norm;
		gate.convergence_threshold=std::max(options.nonlinear_absolute_tolerance,
			options.nonlinear_relative_tolerance*gate.initial_residual);
	}
	gate.residual_converged=gate.nonempty_record && std::isfinite(gate.initial_residual)
		&& std::isfinite(gate.convergence_threshold) && std::isfinite(diagnostics.residual_norm)
		&& diagnostics.residual_norm<=gate.convergence_threshold;
	gate.linear_residual_satisfied=std::isfinite(diagnostics.true_linear_relative_residual)
		&& diagnostics.true_linear_relative_residual<=1e-10;
	PetscInt summed_ksp_iterations=0;
	bool records_consistent=true;
	for(std::size_t i=0;i<steps.size();++i) {
		const auto& step=steps[i];
		summed_ksp_iterations+=step.ksp_iterations;
		records_consistent=records_consistent && step.iteration==static_cast<PetscInt>(i)
			&& step.ksp_reason>0 && std::isfinite(step.residual_norm) && std::isfinite(step.update_norm)
			&& std::isfinite(step.linear_relative_residual) && std::isfinite(step.damping)
			&& step.damping>=options.minimum_damping && step.damping<=1.0;
	}
	gate.accounting_consistent=records_consistent
		&& diagnostics.nonlinear_iterations==static_cast<PetscInt>(steps.size())
		&& diagnostics.ksp_iterations==summed_ksp_iterations
		&& diagnostics.attempt_assembly_count>=2*steps.size();
	return gate;
}

void PrintFirstSolveDiagnostic(const iga::ImmersedTransientFlowDiagnostics& diagnostics,
	const iga::ImmersedTransientFlowOptions& options, const FirstSolveFormulationGate& gate)
{
	std::cout << std::setprecision(17)
		<< "newton_diagnostic updates=" << diagnostics.newton_steps.size()
		<< " configured_max=" << options.nonlinear_maximum_iterations
		<< " initial_residual=" << gate.initial_residual
		<< " final_residual=" << diagnostics.residual_norm
		<< " convergence_threshold=" << gate.convergence_threshold
		<< " true_linear_relative_residual=" << diagnostics.true_linear_relative_residual << '\n';
	for(const auto& step:diagnostics.newton_steps)
		std::cout << "newton_diagnostic step=" << step.iteration << " ksp_iterations=" << step.ksp_iterations
			<< " ksp_reason=" << static_cast<int>(step.ksp_reason) << " residual=" << step.residual_norm
			<< " update=" << step.update_norm << " linear_relative_residual=" << step.linear_relative_residual
			<< " damping=" << step.damping << '\n';
	std::cout << "newton_diagnostic gates nonempty_record=" << (gate.nonempty_record?"true":"false")
		<< " update_count_within_limit=" << (gate.update_count_within_limit?"true":"false")
		<< " converged=" << (gate.converged?"true":"false")
		<< " residual_converged=" << (gate.residual_converged?"true":"false")
		<< " linear_residual_satisfied=" << (gate.linear_residual_satisfied?"true":"false")
		<< " accounting_consistent=" << (gate.accounting_consistent?"true":"false") << '\n';
}

int RunNewtonDiagnostic()
{
	if(iga::CurrentPhaseProfile().Enabled()) std::cout << "fixed_newton_stage=geometry\n";
	const auto soup=Cube(); iga::PrescribedSurfaceMotion motion({{0.0,soup},{1.0,soup}});
	iga::MovingCutGeometryOptions geometry_options;
	geometry_options.volume.max_depth=4; geometry_options.volume.max_nodes=500000;
	geometry_options.volume.max_leaves=500000; geometry_options.volume.max_points=3000000;
	const iga::CubicCartesianGridSpec grid{{{0,0,0}},{{1,1,1}},{{3,3,3}}};
	auto geometry=iga::MovingCutGeometry::Build(grid,motion.Evaluate(1.0,0.0,1.0),geometry_options);
	if(iga::CurrentPhaseProfile().Enabled()) std::cout << "fixed_newton_stage=runtime\n";
	const auto options=Options(); iga::ImmersedTransientFlowRuntime runtime(*geometry,options);
	runtime.SetCommittedGlobalState(iga::ImmersedGlobalFlowState(0.0,0,runtime.Layout(),
		NonconstantFields(runtime.Layout()),{.019,-.023},true,.031));
	runtime.BeginTrial(1.0,1,1.0);
	if(iga::CurrentPhaseProfile().Enabled()) std::cout << "fixed_newton_stage=solve\n";
	const bool solved=runtime.SolveTrial(); const auto& diagnostics=runtime.Diagnostics();
	const auto gate=EvaluateFirstSolveFormulationGate(diagnostics,options);
	PrintFirstSolveDiagnostic(diagnostics,options,gate);
	assert(solved && gate.Passed());
	return 0;
}

}

int main(int argc,char** argv)
{
	bool newton_diagnostic=false;
	for(int i=1;i<argc;++i) if(std::string(argv[i])=="--newton-diagnostic") newton_diagnostic=true;
	PetscInitialize(&argc,&argv,nullptr,nullptr); int status=0;
	iga::CurrentPhaseProfile().EnableFromEnvironment();
	if(newton_diagnostic) std::cout << std::unitbuf;
	try {
		if(newton_diagnostic) status=RunNewtonDiagnostic();
		else {
		const auto soup=Cube(); iga::PrescribedSurfaceMotion motion({{0.0,soup},{1.0,soup}});
		iga::MovingCutGeometryOptions go; go.volume.max_depth=4; go.volume.max_nodes=500000; go.volume.max_leaves=500000; go.volume.max_points=3000000;
		const iga::CubicCartesianGridSpec grid{{{0,0,0}},{{1,1,1}},{{3,3,3}}};
		auto geometry=iga::MovingCutGeometry::Build(grid,motion.Evaluate(1.0,0.0,1.0),go);
		iga::MovingCutGeometryOptions branch_go=go; branch_go.volume.max_depth=3;
		const iga::CubicCartesianGridSpec branch_grid{{{0,0,0}},{{1,1,1}},{{3,3,3}}};
		auto branch_geometry=iga::MovingCutGeometry::Build(branch_grid,motion.Evaluate(1.0,0.0,1.0),branch_go);
		iga::MovingCutGeometryOptions compact_branch_go=branch_go; compact_branch_go.volume_storage=iga::CutCellVolumeQuadratureStorageMode::Compact;
		auto compact_branch_geometry=iga::MovingCutGeometry::Build(branch_grid,motion.Evaluate(1.0,0.0,1.0),compact_branch_go);
		CheckPressureLikePortBranch(*branch_geometry,iga::ImmersedFlowPortControlMode::Pressure,.125);
		CheckPressureLikePortBranch(*branch_geometry,iga::ImmersedFlowPortControlMode::MeanNormalTraction,.125);
		CheckFrozenBodyForce(*branch_geometry);
		CheckBlockReduction(*branch_geometry);
		CheckIdleMutationPublication(*branch_geometry);
		CheckSolverConfigurationIdentity(*branch_geometry);
		CheckBeginTrialExceptionTransaction(*branch_geometry);
		CheckCompactRuntime(*branch_geometry,*compact_branch_geometry);
		CheckTransientNitscheDiagnostics(*branch_geometry);
		CheckMovingInnerTrial(branch_grid,branch_go);
		Reject([&] { auto bad=PressureLikeOptions(iga::ImmersedFlowPortControlMode::Pressure,.1); bad.ports[0].boundary_label=0; iga::ImmersedTransientFlowRuntime rejected(*branch_geometry,bad); });
		Reject([&] { auto bad=PressureLikeOptions(iga::ImmersedFlowPortControlMode::Pressure,.1); bad.ports[0].boundary_label=99; iga::ImmersedTransientFlowRuntime rejected(*branch_geometry,bad); });
		for(const double invalid:{std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()}) {
			Reject([&] { auto bad=Options(); bad.wall_inertial_gamma0=invalid; iga::ImmersedTransientFlowRuntime rejected(*branch_geometry,bad); });
			Reject([&] { auto bad=Options(); bad.flow_controller_relative_tolerance=invalid; iga::ImmersedTransientFlowRuntime rejected(*branch_geometry,bad); });
			Reject([&] { auto bad=Options(); bad.flow_controller_absolute_tolerance_m3_s=invalid; iga::ImmersedTransientFlowRuntime rejected(*branch_geometry,bad); });
			Reject([&] { auto bad=Options(); bad.flow_controller_reference_flow_m3_s=invalid; iga::ImmersedTransientFlowRuntime rejected(*branch_geometry,bad); });
			Reject([&] { auto bad=Options(); bad.lu_pivot_shift=invalid; iga::ImmersedTransientFlowRuntime rejected(*branch_geometry,bad); });
			Reject([&] { auto bad=Options(); bad.minimum_damping=invalid; iga::ImmersedTransientFlowRuntime rejected(*branch_geometry,bad); });
		}
		Reject([&] { auto bad=Options(); bad.ports[0].value=.25; iga::ImmersedTransientFlowRuntime rejected(*branch_geometry,bad); });
		Reject([&] { auto bad=Options(); bad.include_pressure_gauge=false; iga::ImmersedTransientFlowRuntime rejected(*branch_geometry,bad); });
		Reject([&] { auto bad=Options(); bad.ports[1].boundary_label=1; iga::ImmersedTransientFlowRuntime rejected(*branch_geometry,bad); });
		Reject([&] { auto bad=Options(); bad.wall_labels={1,7}; iga::ImmersedTransientFlowRuntime rejected(*branch_geometry,bad); });
		{ auto bad=Options(); bad.body_force=[](const std::array<double,3>&) { return std::array<double,3>{{std::numeric_limits<double>::infinity(),0,0}}; };
			iga::ImmersedTransientFlowRuntime rejected(*branch_geometry,bad);
			Reject([&] { rejected.BeginTrial(1.0,1,1.0); });
			assert(rejected.Diagnostics().idle && !rejected.Diagnostics().trial_active && rejected.Diagnostics().committed && !rejected.Diagnostics().prepared && rejected.Diagnostics().input_hash_sha256.empty() && rejected.Diagnostics().history_hash_sha256.empty()); }

		auto translated=Cube(); for(auto& vertex:translated.vertices) vertex[0]+=.01;
		iga::PrescribedSurfaceMotion moving_motion({{0.0,soup},{1.0,translated}});
		auto moving_geometry=iga::MovingCutGeometry::Build(branch_grid,moving_motion.Evaluate(1.0,0.0,1.0),branch_go);
		const double expected_velocity_m_per_s=.01/(1.0-0.0);
		const double velocity_m_per_s=moving_geometry->Evaluation().SourceVertexVelocitiesMPerS().front()[0];
		const double velocity_tolerance=32.0*std::numeric_limits<double>::epsilon()*std::max(std::abs(expected_velocity_m_per_s),std::numeric_limits<double>::denorm_min());
		assert(std::isfinite(velocity_m_per_s) && std::abs(velocity_m_per_s)>0.5*std::abs(expected_velocity_m_per_s));
		assert(std::abs(velocity_m_per_s-expected_velocity_m_per_s)<=velocity_tolerance);
		iga::ImmersedTransientFlowRuntime moving_runtime(*moving_geometry,Options());
		moving_runtime.SetCommittedGlobalState(iga::ImmersedGlobalFlowState(0.0,0,moving_runtime.Layout(),NonconstantFields(moving_runtime.Layout()),{0.0,0.0},true,0.0));
		const auto moving_before=moving_runtime.Diagnostics(); const auto moving_state=moving_runtime.CommittedState(); const auto moving_trial=moving_runtime.TrialState(); const auto moving_hash=moving_runtime.CommittedGlobalState().HashSha256();
		RejectWithMessage([&] { moving_runtime.BeginTrial(1.0,1,1.0); },"immersed transient fixed geometry requires exactly zero material wall velocity");
		assert(moving_runtime.CommittedState()==moving_state && moving_runtime.TrialState()==moving_trial && moving_runtime.CommittedGlobalState().HashSha256()==moving_hash && SameDiagnosticsPublic(moving_before,moving_runtime.Diagnostics()));
		iga::ImmersedTransientFlowRuntime runtime(*geometry,Options());
		assert(runtime.Layout().PortIds().size()==2 && runtime.Layout().HasGaugeRow());
		assert(runtime.Diagnostics().scalar_diagonal_structure_verified);
		assert(runtime.PortMultiplierDof("left")==static_cast<PetscInt>(runtime.Layout().NodeFieldRows()));
		assert(runtime.PortMultiplierDof("right")==static_cast<PetscInt>(runtime.Layout().NodeFieldRows()+1));
		assert(runtime.GaugeDof()==static_cast<PetscInt>(runtime.Layout().NodeFieldRows()+2));

		const iga::ImmersedGlobalFlowState seeded(0.0,0,runtime.Layout(),NonconstantFields(runtime.Layout()),{.019,-.023},true,.031);
		runtime.SetCommittedGlobalState(seeded); const auto committed_before=runtime.CommittedGlobalState(); const auto committed_vector=runtime.CommittedState();
		Reject([&] { runtime.SetTrialState({}); }); Reject([&] { runtime.BeginTrial(1.0,1,0.0); }); Reject([&] { runtime.BeginTrial(1.0,2,1.0); });

		iga::PrescribedSurfaceMotion candidate_motion({{0.0,Cube(.2,.8)},{1.0,Cube(.2,.8)}});
		auto candidate=iga::MovingCutGeometry::Build(grid,candidate_motion.Evaluate(1.0,0.0,1.0),go);
		iga::ImmersedTransientFlowRuntime candidate_runtime(*candidate,Options());
		assert(candidate_runtime.Layout().GeometryIdentity()!=runtime.Layout().GeometryIdentity());
		Reject([&] { runtime.SetCommittedGlobalState(candidate_runtime.CommittedGlobalState()); });
		iga::ImmersedTransientFlowRuntime overflow_runtime(*geometry,Options());
		overflow_runtime.SetCommittedGlobalState(iga::ImmersedGlobalFlowState(0.0,std::numeric_limits<std::uint64_t>::max(),overflow_runtime.Layout(),NonconstantFields(overflow_runtime.Layout()),{0.0,0.0},true,0.0));
		Reject([&] { overflow_runtime.BeginTrial(1.0,0,1.0); });

		runtime.BeginTrial(1.0,1,1.0);
		assert(runtime.Diagnostics().identity_history_nodes==runtime.Layout().NodeIds().size() && runtime.Diagnostics().committed_history_nodes==runtime.Layout().NodeIds().size());
		assert(runtime.Diagnostics().extended_history_nodes==0 && runtime.Diagnostics().missing_history_nodes==0 && !runtime.Diagnostics().history_hash_sha256.empty());
		const auto frozen=runtime.TrialState(); const auto frozen_hash=runtime.Diagnostics().trial_state_hash_sha256;
		Reject([&] { auto bad=frozen; bad.push_back(0.0); runtime.SetTrialState(bad); }); Reject([&] { auto bad=frozen; bad[0]=std::numeric_limits<double>::quiet_NaN(); runtime.SetTrialState(bad); });

		const auto x=NonconstantTrial(runtime), direction=Direction(runtime); runtime.SetTrialState(x); const auto committed_ports_before_assemble=runtime.Diagnostics().ports; const auto committed_hash_before_assemble=runtime.CommittedGlobalState().HashSha256(); const auto committed_state_before_assemble=runtime.CommittedState(); runtime.Assemble();
		assert(runtime.CommittedGlobalState().HashSha256()==committed_hash_before_assemble && runtime.CommittedState()==committed_state_before_assemble && SamePortsBitwise(committed_ports_before_assemble,runtime.Diagnostics().ports));
		assert(runtime.Diagnostics().volume_cells>0 && runtime.Diagnostics().surface_cells>0 && runtime.Diagnostics().ghost_faces>0);
		const auto jd=runtime.AssembledJacobianAction(direction);
		assert(std::abs(PetscRealPart(jd[static_cast<std::size_t>(runtime.PortMultiplierDof("left"))]))>1e-12 && std::abs(PetscRealPart(jd[static_cast<std::size_t>(runtime.GaugeDof())]))>1e-12);
		double fd_max=0.0, zero_block_max=0.0;
		for(const double epsilon:{1e-4,1e-5,1e-6}) {
			auto plus=x,minus=x; for(std::size_t i=0;i<x.size();++i) { plus[i]+=epsilon*direction[i]; minus[i]-=epsilon*direction[i]; }
			runtime.SetTrialState(plus); runtime.Assemble(); const auto bp=runtime.AssembledNegativeResidual(); runtime.SetTrialState(minus); runtime.Assemble(); const auto bm=runtime.AssembledNegativeResidual();
			for(std::size_t i=0;i<x.size();++i) { const double derivative=PetscRealPart((bp[i]-bm[i])/(2.0*epsilon)); const double error=PetscRealPart(jd[i])+derivative; fd_max=std::max(fd_max,std::abs(error)); if(std::abs(PetscRealPart(jd[i]))+std::abs(derivative)<1e-13) zero_block_max=std::max(zero_block_max,std::abs(error)); }
		}
		runtime.SetTrialState(x); assert(runtime.TrialState()==x && fd_max<=1e-8 && zero_block_max<=1e-11);

		runtime.Rollback(); assert(runtime.TrialState()==frozen && runtime.Diagnostics().trial_state_hash_sha256==frozen_hash);
		const auto committed_ports_before_solve=runtime.Diagnostics().ports; const auto committed_hash_before_solve=runtime.CommittedGlobalState().HashSha256(); const auto committed_state_before_solve=runtime.CommittedState();
		assert(runtime.SolveTrial()); const auto solved=runtime.TrialState(); const auto ports=runtime.Diagnostics().ports; const auto trial_hash=runtime.Diagnostics().trial_state_hash_sha256; const auto history_hash=runtime.Diagnostics().history_hash_sha256;
		const auto input_hash=runtime.Diagnostics().input_hash_sha256, solved_hash=runtime.Diagnostics().solved_state_hash_sha256, attempt_hash=runtime.Diagnostics().attempt_hash_sha256, layout_hash=runtime.Diagnostics().layout_hash_sha256;
		const auto nonlinear_steps=runtime.Diagnostics().newton_steps; const auto first_ksp=runtime.Diagnostics().ksp_iterations; const auto first_assemblies=runtime.Diagnostics().attempt_assembly_count; const auto first_volume=runtime.Diagnostics().volume_cells, first_surface=runtime.Diagnostics().surface_cells, first_ghost=runtime.Diagnostics().ghost_faces;
		const auto first_solve_gate=EvaluateFirstSolveFormulationGate(runtime.Diagnostics(),Options());
		assert(first_solve_gate.Passed());
		assert(runtime.CommittedGlobalState().HashSha256()==committed_hash_before_solve && runtime.CommittedState()==committed_state_before_solve && SamePortsBitwise(committed_ports_before_solve,runtime.Diagnostics().ports));
		assert(SamePorts(ports,runtime.Diagnostics().ports));
		auto stale=solved; stale[0]+=1e-12; runtime.SetTrialState(stale); assert(!runtime.Diagnostics().converged && runtime.Diagnostics().solved_state_hash_sha256.empty() && runtime.Diagnostics().attempt_hash_sha256.empty() && runtime.Diagnostics().newton_steps.empty()); Reject([&] { runtime.PrepareCommit(); });
		runtime.Rollback(); assert(runtime.TrialState()==frozen && runtime.SolveTrial());
		assert(runtime.Diagnostics().trial_state_hash_sha256==trial_hash && runtime.Diagnostics().attempt_hash_sha256.size()==64 && runtime.Diagnostics().attempt_assembly_count>0);
		const auto before_prepare_ports=runtime.Diagnostics().ports; runtime.FailNextPrepareForTesting(); Reject([&] { runtime.PrepareCommit(); });
		assert(runtime.CommittedGlobalState().HashSha256()==committed_before.HashSha256() && runtime.CommittedState()==committed_vector && SamePortsBitwise(before_prepare_ports,runtime.Diagnostics().ports) && !runtime.Diagnostics().prepared && runtime.Diagnostics().prepared_hash_sha256.empty());
		runtime.Rollback(); assert(runtime.TrialState()==frozen && runtime.SolveTrial());
		assert(runtime.TrialState()==solved && SamePortsBitwise(ports,runtime.Diagnostics().ports) && runtime.Diagnostics().trial_state_hash_sha256==trial_hash && runtime.Diagnostics().history_hash_sha256==history_hash);
		assert(runtime.Diagnostics().input_hash_sha256==input_hash && runtime.Diagnostics().layout_hash_sha256==layout_hash && runtime.Diagnostics().solved_state_hash_sha256==solved_hash && runtime.Diagnostics().attempt_hash_sha256==attempt_hash && runtime.Diagnostics().ksp_iterations==first_ksp && runtime.Diagnostics().attempt_assembly_count==first_assemblies && SameNewtonRecords(runtime.Diagnostics().newton_steps,nonlinear_steps) && runtime.Diagnostics().volume_cells==first_volume && runtime.Diagnostics().surface_cells==first_surface && runtime.Diagnostics().ghost_faces==first_ghost);
		runtime.PrepareCommit(); const auto prepared_hash=runtime.Diagnostics().prepared_hash_sha256; assert(!prepared_hash.empty() && SamePortsBitwise(before_prepare_ports,runtime.Diagnostics().ports)); runtime.FinalizeCommit(); const auto finalized=runtime.CommittedState(); const auto final_hash=runtime.CommittedGlobalState().HashSha256(); const auto finalized_ports=runtime.Diagnostics().ports; const auto finalize_count=runtime.Diagnostics().finalize_count;
		assert(finalized==solved && runtime.CommittedGlobalState().TimeS()==1.0 && runtime.CommittedGlobalState().Index()==1 && runtime.Diagnostics().prepared_hash_sha256.empty());
		for(const auto& port:runtime.Diagnostics().ports) { assert(port.measurement.area_m2>0); if(port.multiplier_row>=0) assert(port.multiplier==PetscRealPart(finalized[static_cast<std::size_t>(port.multiplier_row)])); }
		runtime.FinalizeCommit(); assert(runtime.CommittedState()==finalized && runtime.CommittedGlobalState().HashSha256()==final_hash && SamePortsBitwise(finalized_ports,runtime.Diagnostics().ports) && runtime.Diagnostics().finalize_count==finalize_count);
		Reject([&] { runtime.SetTrialState(finalized); }); Reject([&] { runtime.Assemble(); });
		const auto conservation=runtime.ConservationDiagnostics(); assert(std::abs(conservation.normalized_open_balance)<=1e-3 && std::abs(conservation.normalized_wall_leakage)<=1e-3);
		iga::ImmersedTransientFlowRuntime conservation_runtime(*geometry,Options()); conservation_runtime.SetCommittedGlobalState(iga::ImmersedGlobalFlowState(0.0,0,conservation_runtime.Layout(),NonconstantFields(conservation_runtime.Layout()),{0.0,0.0},true,0.0)); conservation_runtime.BeginTrial(1.0,1,1.0); const auto nondivergent=GenericTrial(conservation_runtime); conservation_runtime.SetTrialState(nondivergent); conservation_runtime.Assemble(); const auto nondivergence=conservation_runtime.ConservationDiagnostics();
		double signed_sum=0.0; for(const auto& item:nondivergence.surface_flow_by_boundary_label_m3_s) signed_sum+=item.second;
		assert(SameBits(signed_sum,nondivergence.total_surface_outward_flow_m3_s));
		double material_sum=0.0; for(const auto& item:nondivergence.material_surface_outward_flow_by_boundary_label_m3_s) material_sum+=item.second;
		assert(std::abs(material_sum-nondivergence.total_material_surface_outward_flow_m3_s)<=2e-11);
		double material_wall_sum=0.0; for(const auto& item:nondivergence.material_wall_outward_flow_by_boundary_label_m3_s) material_wall_sum+=item.second;
		assert(std::abs(material_wall_sum-nondivergence.total_material_wall_outward_flow_m3_s)<=2e-11);
		assert(std::abs(nondivergence.total_material_wall_outward_flow_m3_s)<=2e-11);
		assert(SameBits(nondivergence.wall_relative_leakage_m3_s,nondivergence.wall_outward_flow_m3_s));
		assert(SameBits(nondivergence.divergence_theorem_defect_m3_s,nondivergence.endpoint_volume_divergence_m3_s-nondivergence.total_surface_outward_flow_m3_s));
		assert(SameBits(nondivergence.discrete_moving_wall_continuity_defect_m3_s,nondivergence.open_port_outward_flow_m3_s+nondivergence.total_material_wall_outward_flow_m3_s));
		assert(SameBits(nondivergence.discrete_moving_wall_continuity_normalization_scale_m3_s,std::max({Options().flow_controller_reference_flow_m3_s,std::abs(nondivergence.open_port_outward_flow_m3_s),std::abs(nondivergence.total_material_wall_outward_flow_m3_s)})));
		assert(SameBits(nondivergence.normalized_discrete_moving_wall_continuity_defect,std::abs(nondivergence.discrete_moving_wall_continuity_defect_m3_s)/nondivergence.discrete_moving_wall_continuity_normalization_scale_m3_s));
		assert(std::abs((nondivergence.discrete_moving_wall_continuity_defect_m3_s+nondivergence.wall_relative_leakage_m3_s)-nondivergence.total_surface_outward_flow_m3_s)<=2e-11);
		const double corrected=std::abs(nondivergence.endpoint_volume_divergence_m3_s-nondivergence.total_surface_outward_flow_m3_s), old_plus=std::abs(nondivergence.endpoint_volume_divergence_m3_s+nondivergence.total_surface_outward_flow_m3_s);
		assert(corrected<old_plus);
		conservation_runtime.AbortTrial();
		CheckAbortPublication(*geometry);
		std::cout << "immersed_transient_flow_tests=passed active_nodes=" << runtime.Diagnostics().active_nodes << " fd_max=" << fd_max << " zero_block_max=" << zero_block_max << " newton_updates=" << nonlinear_steps.size() << " linear_relative=" << runtime.Diagnostics().true_linear_relative_residual << " open_balance=" << conservation.normalized_open_balance << " wall_leakage=" << conservation.normalized_wall_leakage << '\n';
		}
	} catch(const std::exception& error) { std::cerr << "immersed_transient_flow_test: " << error.what() << '\n'; status=1; }
	int rank=0,ranks=1; MPI_Comm_rank(PETSC_COMM_WORLD,&rank); MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	iga::CurrentPhaseProfile().Write(std::cout,rank,ranks,status);
	PetscFinalize(); return status;
}
