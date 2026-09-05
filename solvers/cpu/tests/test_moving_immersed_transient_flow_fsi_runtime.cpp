#include "MovingImmersedTransientFlowFsiRuntime.hpp"
#include "PrescribedSurfaceMotion.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace {

template <class Function> void Reject(Function&& function)
{ bool rejected=false; try { function(); } catch (const std::exception&) { rejected=true; } assert(rejected); }
void Require(bool value, const char* message)
{ if (!value) throw std::runtime_error(message); }
double Norm(const std::array<double,3>& value)
{ return std::sqrt(value[0]*value[0]+value[1]*value[1]+value[2]*value[2]); }
bool Bits(double left, double right)
{ std::uint64_t a=0,b=0; std::memcpy(&a,&left,sizeof(a)); std::memcpy(&b,&right,sizeof(b)); return a==b; }
bool SameMap(const std::map<int,double>& left, const std::map<int,double>& right)
{
	if(left.size()!=right.size()) return false;
	auto a=left.begin(),b=right.begin();
	for(;a!=left.end();++a,++b) if(a->first!=b->first || !Bits(a->second,b->second)) return false;
	return true;
}

bool SameSurfaceInterface(const iga::SurfaceInterfaceRef& left, const iga::SurfaceInterfaceRef& right)
{ return left.domain_id==right.domain_id && left.subsystem_id==right.subsystem_id && left.interface_id==right.interface_id; }
bool SameStamp(const iga::SurfaceFieldStamp& left, const iga::SurfaceFieldStamp& right)
{
	return Bits(left.time_s,right.time_s) && left.step==right.step && left.coupling_iteration==right.coupling_iteration
		&& left.reference_mesh_identity_sha256==right.reference_mesh_identity_sha256
		&& left.layout_identity_sha256==right.layout_identity_sha256
		&& left.partition_identity_sha256==right.partition_identity_sha256
		&& left.producer_state_identity_sha256==right.producer_state_identity_sha256;
}
bool SameVector3(const std::array<double,3>& left, const std::array<double,3>& right)
{ return Bits(left[0],right[0]) && Bits(left[1],right[1]) && Bits(left[2],right[2]); }
bool SameVector3s(const std::vector<std::array<double,3>>& left, const std::vector<std::array<double,3>>& right)
{
	if(left.size()!=right.size()) return false;
	for(std::size_t i=0;i<left.size();++i) if(!SameVector3(left[i],right[i])) return false;
	return true;
}
bool SameTraction(const iga::SurfaceTraction& left, const iga::SurfaceTraction& right)
{
	return SameSurfaceInterface(left.interface,right.interface) && SameStamp(left.stamp,right.stamp)
		&& left.projection_identity_sha256==right.projection_identity_sha256
		&& SameVector3s(left.traction_on_structure_pa,right.traction_on_structure_pa)
		&& SameVector3s(left.consistent_nodal_force_n,right.consistent_nodal_force_n);
}
bool SameTractionDiagnostics(const iga::FluidSurfaceTractionDiagnostics& left,
	const iga::FluidSurfaceTractionDiagnostics& right)
{
	return SameVector3(left.quadrature_resultant_n,right.quadrature_resultant_n)
		&& SameVector3(left.nodal_resultant_n,right.nodal_resultant_n)
		&& SameVector3(left.quadrature_moment_n_m,right.quadrature_moment_n_m)
		&& SameVector3(left.nodal_moment_n_m,right.nodal_moment_n_m)
		&& left.retained_quadrature_points==right.retained_quadrature_points;
}
bool SameWallLabel(const iga::ImmersedNitscheWallLabelDiagnostics& left,
	const iga::ImmersedNitscheWallLabelDiagnostics& right)
{
	return left.selected_points==right.selected_points && left.skipped_points==right.skipped_points
		&& Bits(left.selected_area_m2,right.selected_area_m2) && Bits(left.skipped_area_m2,right.skipped_area_m2);
}
bool SameWallDiagnostics(const iga::ImmersedNitscheWallDiagnostics& left,
	const iga::ImmersedNitscheWallDiagnostics& right)
{
	if(left.by_boundary_id.size()!=right.by_boundary_id.size()) return false;
	auto a=left.by_boundary_id.begin(),b=right.by_boundary_id.begin();
	for(;a!=left.by_boundary_id.end();++a,++b) if(a->first!=b->first || !SameWallLabel(a->second,b->second)) return false;
	return Bits(left.fraction_lower,right.fraction_lower) && Bits(left.fraction_estimate,right.fraction_estimate)
		&& Bits(left.fraction_upper,right.fraction_upper) && Bits(left.minimum_h_n_m,right.minimum_h_n_m)
		&& Bits(left.maximum_h_n_m,right.maximum_h_n_m) && Bits(left.minimum_eta,right.minimum_eta)
		&& Bits(left.maximum_eta,right.maximum_eta) && Bits(left.maximum_eta_h_n_over_mu,right.maximum_eta_h_n_over_mu)
		&& Bits(left.minimum_eta_mu,right.minimum_eta_mu) && Bits(left.maximum_eta_mu,right.maximum_eta_mu)
		&& Bits(left.minimum_eta_t,right.minimum_eta_t) && Bits(left.maximum_eta_t,right.maximum_eta_t)
		&& Bits(left.minimum_eta_mu_h_n_over_mu,right.minimum_eta_mu_h_n_over_mu)
		&& Bits(left.maximum_eta_mu_h_n_over_mu,right.maximum_eta_mu_h_n_over_mu)
		&& Bits(left.minimum_eta_t_dt_over_rho_h_n,right.minimum_eta_t_dt_over_rho_h_n)
		&& Bits(left.maximum_eta_t_dt_over_rho_h_n,right.maximum_eta_t_dt_over_rho_h_n)
		&& Bits(left.minimum_eta_mu_fraction,right.minimum_eta_mu_fraction)
		&& Bits(left.maximum_eta_mu_fraction,right.maximum_eta_mu_fraction)
		&& Bits(left.minimum_eta_t_fraction,right.minimum_eta_t_fraction)
		&& Bits(left.maximum_eta_t_fraction,right.maximum_eta_t_fraction)
		&& Bits(left.maximum_gap_norm,right.maximum_gap_norm) && left.ghost_covered_policy==right.ghost_covered_policy;
}
bool SamePort(const iga::ImmersedTransientFlowDiagnostics::Port& left,
	const iga::ImmersedTransientFlowDiagnostics::Port& right)
{
	return left.id==right.id && left.boundary_label==right.boundary_label && left.control_mode==right.control_mode
		&& Bits(left.target,right.target) && Bits(left.multiplier,right.multiplier)
		&& Bits(left.controller_error,right.controller_error) && left.multiplier_row==right.multiplier_row
		&& Bits(left.measurement.area_m2,right.measurement.area_m2)
		&& Bits(left.measurement.outward_flow_m3_s,right.measurement.outward_flow_m3_s)
		&& Bits(left.measurement.mean_pressure_pa,right.measurement.mean_pressure_pa)
		&& Bits(left.measurement.mean_normal_traction_pa,right.measurement.mean_normal_traction_pa)
		&& Bits(left.measurement.mean_velocity_squared_m2_s2,right.measurement.mean_velocity_squared_m2_s2)
		&& left.measurement_valid==right.measurement_valid;
}
bool SamePorts(const std::vector<iga::ImmersedTransientFlowDiagnostics::Port>& left,
	const std::vector<iga::ImmersedTransientFlowDiagnostics::Port>& right)
{
	if(left.size()!=right.size()) return false;
	for(std::size_t i=0;i<left.size();++i) if(!SamePort(left[i],right[i])) return false;
	return true;
}
bool SameNewton(const std::vector<iga::ImmersedTransientFlowNewtonStep>& left,
	const std::vector<iga::ImmersedTransientFlowNewtonStep>& right)
{
	if(left.size()!=right.size()) return false;
	for(std::size_t i=0;i<left.size();++i) if(left[i].iteration!=right[i].iteration
		|| left[i].ksp_iterations!=right[i].ksp_iterations || left[i].ksp_reason!=right[i].ksp_reason
		|| !Bits(left[i].residual_norm,right[i].residual_norm) || !Bits(left[i].update_norm,right[i].update_norm)
		|| !Bits(left[i].linear_relative_residual,right[i].linear_relative_residual) || !Bits(left[i].damping,right[i].damping)) return false;
	return true;
}
bool SameFlowTransfer(const iga::ImmersedTransientFlowDiagnostics& left,
	const iga::ImmersedTransientFlowDiagnostics& right)
{
	return left.active_nodes==right.active_nodes && left.physical_dofs==right.physical_dofs && left.total_dofs==right.total_dofs
		&& left.volume_cells==right.volume_cells && left.surface_cells==right.surface_cells && left.ghost_faces==right.ghost_faces
		&& Bits(left.target_time_s,right.target_time_s) && Bits(left.dt_s,right.dt_s)
		&& Bits(left.pressure_measure,right.pressure_measure) && Bits(left.pressure_gauge_defect,right.pressure_gauge_defect)
		&& Bits(left.residual_norm,right.residual_norm) && Bits(left.true_linear_relative_residual,right.true_linear_relative_residual)
		&& left.nonlinear_iterations==right.nonlinear_iterations && left.ksp_iterations==right.ksp_iterations
		&& left.identity_history_nodes==right.identity_history_nodes && left.committed_history_nodes==right.committed_history_nodes
		&& left.extended_history_nodes==right.extended_history_nodes && left.missing_history_nodes==right.missing_history_nodes
		&& left.ksp_reason==right.ksp_reason && left.idle==right.idle && left.trial_active==right.trial_active
		&& left.converged==right.converged && left.prepared==right.prepared && left.committed==right.committed
		&& left.scalar_diagonal_structure_verified==right.scalar_diagonal_structure_verified
		&& SameWallDiagnostics(left.wall_penalty,right.wall_penalty)
		&& left.geometry_identity_sha256==right.geometry_identity_sha256 && left.layout_hash_sha256==right.layout_hash_sha256
		&& left.committed_state_hash_sha256==right.committed_state_hash_sha256 && left.trial_state_hash_sha256==right.trial_state_hash_sha256
		&& left.history_hash_sha256==right.history_hash_sha256 && left.input_hash_sha256==right.input_hash_sha256
		&& left.solved_state_hash_sha256==right.solved_state_hash_sha256 && left.prepared_hash_sha256==right.prepared_hash_sha256
		&& left.attempt_hash_sha256==right.attempt_hash_sha256 && left.moving_map_identity_sha256==right.moving_map_identity_sha256
		&& left.attempt_assembly_count==right.attempt_assembly_count && SamePorts(left.ports,right.ports)
		&& SameNewton(left.newton_steps,right.newton_steps);
}
bool SameConservation(const iga::MovingImmersedTransientFlowConservationDiagnostics& a,
	const iga::MovingImmersedTransientFlowConservationDiagnostics& b)
{
	return a.source_geometry_identity_sha256==b.source_geometry_identity_sha256
		&& a.source_publication_identity_sha256==b.source_publication_identity_sha256
		&& a.target_geometry_identity_sha256==b.target_geometry_identity_sha256
		&& a.target_publication_identity_sha256==b.target_publication_identity_sha256
		&& a.transition_identity_sha256==b.transition_identity_sha256
		&& Bits(a.source_time_s,b.source_time_s) && Bits(a.target_time_s,b.target_time_s) && Bits(a.dt_s,b.dt_s)
		&& a.source_index==b.source_index && a.target_index==b.target_index
		&& Bits(a.source_audited_volume_m3,b.source_audited_volume_m3) && Bits(a.target_audited_volume_m3,b.target_audited_volume_m3)
		&& Bits(a.backward_euler_volume_rate_m3_s,b.backward_euler_volume_rate_m3_s)
		&& SameMap(a.fluid_surface_outward_flow_by_boundary_label_m3_s,b.fluid_surface_outward_flow_by_boundary_label_m3_s)
		&& SameMap(a.material_surface_outward_flow_by_boundary_label_m3_s,b.material_surface_outward_flow_by_boundary_label_m3_s)
		&& SameMap(a.material_wall_outward_flow_by_boundary_label_m3_s,b.material_wall_outward_flow_by_boundary_label_m3_s)
		&& Bits(a.endpoint_volume_divergence_m3_s,b.endpoint_volume_divergence_m3_s)
		&& Bits(a.total_fluid_surface_outward_flow_m3_s,b.total_fluid_surface_outward_flow_m3_s)
		&& Bits(a.open_port_outward_flow_m3_s,b.open_port_outward_flow_m3_s)
		&& Bits(a.wall_outward_flow_m3_s,b.wall_outward_flow_m3_s)
		&& Bits(a.total_material_surface_outward_flow_m3_s,b.total_material_surface_outward_flow_m3_s)
		&& Bits(a.total_material_wall_outward_flow_m3_s,b.total_material_wall_outward_flow_m3_s)
		&& Bits(a.divergence_theorem_defect_m3_s,b.divergence_theorem_defect_m3_s)
		&& Bits(a.reynolds_defect_m3_s,b.reynolds_defect_m3_s) && Bits(a.moving_mass_defect_m3_s,b.moving_mass_defect_m3_s)
		&& Bits(a.normalized_moving_mass_defect,b.normalized_moving_mass_defect)
		&& Bits(a.normalized_wall_relative_leakage,b.normalized_wall_relative_leakage)
		&& Bits(a.normalized_discrete_moving_wall_continuity_defect,b.normalized_discrete_moving_wall_continuity_defect)
		&& Bits(a.wall_relative_leakage_m3_s,b.wall_relative_leakage_m3_s)
		&& Bits(a.discrete_moving_wall_continuity_defect_m3_s,b.discrete_moving_wall_continuity_defect_m3_s)
		&& Bits(a.discrete_moving_wall_continuity_normalization_scale_m3_s,b.discrete_moving_wall_continuity_normalization_scale_m3_s)
		&& Bits(a.legacy_normalization_scale_m3_s,b.legacy_normalization_scale_m3_s)
		&& Bits(a.normalization_scale_m3_s,b.normalization_scale_m3_s)
		&& Bits(a.normalized_divergence_theorem_defect,b.normalized_divergence_theorem_defect)
		&& Bits(a.normalized_reynolds_defect,b.normalized_reynolds_defect)
		&& Bits(a.normalized_open_balance,b.normalized_open_balance) && Bits(a.normalized_wall_leakage,b.normalized_wall_leakage);
}
bool SameMovingTransfer(const iga::MovingImmersedTransientFlowDiagnostics& left,
	const iga::MovingImmersedTransientFlowDiagnostics& right)
{
	return left.committed_geometry_identity_sha256==right.committed_geometry_identity_sha256
		&& left.committed_publication_identity_sha256==right.committed_publication_identity_sha256
		&& left.trial_geometry_identity_sha256==right.trial_geometry_identity_sha256
		&& left.trial_publication_identity_sha256==right.trial_publication_identity_sha256
		&& left.source_geometry_identity_sha256==right.source_geometry_identity_sha256
		&& left.source_publication_identity_sha256==right.source_publication_identity_sha256
		&& left.target_geometry_identity_sha256==right.target_geometry_identity_sha256
		&& left.target_publication_identity_sha256==right.target_publication_identity_sha256
		&& left.extension_hash_sha256==right.extension_hash_sha256
		&& left.scalar_extension_hash_sha256==right.scalar_extension_hash_sha256
		&& left.map_identity_sha256==right.map_identity_sha256
		&& left.transition_identity_sha256==right.transition_identity_sha256
		&& Bits(left.source_time_s,right.source_time_s) && Bits(left.target_time_s,right.target_time_s)
		&& Bits(left.dt_s,right.dt_s) && left.source_index==right.source_index && left.target_index==right.target_index
		&& Bits(left.source_audited_volume_m3,right.source_audited_volume_m3)
		&& Bits(left.target_audited_volume_m3,right.target_audited_volume_m3);
}
struct TrialSnapshot {
	iga::SurfaceTraction traction;
	iga::FluidSurfaceTractionDiagnostics traction_diagnostics;
	iga::MaterialSurfaceKinematics target;
	iga::ImmersedTransientFlowDiagnostics flow;
	iga::MovingImmersedTransientFlowConservationDiagnostics conservation;
	iga::MovingImmersedTransientFlowDiagnostics moving;
	std::string composition;
};
TrialSnapshot TakeTrial(const iga::MovingImmersedTransientFlowFsiRuntime& runtime)
{
	return {runtime.GetSurfaceTraction("patch"),runtime.TrialSurfaceTractionDiagnostics(),
		runtime.TrialMaterialKinematicsSnapshot(),runtime.TrialFlowDiagnostics(),runtime.ConservationDiagnostics(),
		runtime.MovingDiagnostics(),runtime.TrialCompositionIdentitySha256()};
}
bool SameTrial(const TrialSnapshot& left, const TrialSnapshot& right)
{
	return SameTraction(left.traction,right.traction)
		&& SameTractionDiagnostics(left.traction_diagnostics,right.traction_diagnostics)
		&& left.target.ContentIdentitySha256()==right.target.ContentIdentitySha256()
		&& SameFlowTransfer(left.flow,right.flow) && SameConservation(left.conservation,right.conservation)
		&& SameMovingTransfer(left.moving,right.moving) && left.composition==right.composition;
}
template <class Function> std::optional<std::string> TryIdentity(Function&& function)
{
	try { return function(); } catch (const std::logic_error&) { return std::nullopt; }
}
template <class Function> std::optional<iga::FluidSurfaceTractionDiagnostics> TryTractionDiagnostics(Function&& function)
{
	try { return function(); } catch (const std::logic_error&) { return std::nullopt; }
}
struct CommittedSnapshot {
	std::string material, geometry, publication, state;
	std::optional<iga::SurfaceTraction> traction;
	std::optional<iga::FluidSurfaceTractionDiagnostics> traction_diagnostics;
	std::optional<std::string> composition;
};
CommittedSnapshot TakeCommitted(const iga::MovingImmersedTransientFlowFsiRuntime& runtime)
{
	return {runtime.CommittedMaterialKinematicsSnapshot().ContentIdentitySha256(),
		runtime.CommittedGeometry().GeometryIdentitySha256(),runtime.CommittedGeometry().PublicationIdentitySha256(),
		runtime.CommittedGlobalState().HashSha256(),runtime.CommittedSurfaceTractionSnapshot(),
		TryTractionDiagnostics([&]{ return runtime.CommittedSurfaceTractionDiagnostics(); }),
		TryIdentity([&]{ return runtime.CommittedCompositionIdentitySha256(); })};
}
bool SameOptionalTraction(const std::optional<iga::SurfaceTraction>& left,
	const std::optional<iga::SurfaceTraction>& right)
{ return left.has_value()==right.has_value() && (!left.has_value() || SameTraction(*left,*right)); }
bool SameOptionalTractionDiagnostics(const std::optional<iga::FluidSurfaceTractionDiagnostics>& left,
	const std::optional<iga::FluidSurfaceTractionDiagnostics>& right)
{ return left.has_value()==right.has_value() && (!left.has_value() || SameTractionDiagnostics(*left,*right)); }
void CheckCommitted(const iga::MovingImmersedTransientFlowFsiRuntime& runtime, const CommittedSnapshot& expected)
{
	const auto actual=TakeCommitted(runtime);
	Require(actual.material==expected.material && actual.geometry==expected.geometry && actual.publication==expected.publication
		&& actual.state==expected.state && SameOptionalTraction(actual.traction,expected.traction)
		&& SameOptionalTractionDiagnostics(actual.traction_diagnostics,expected.traction_diagnostics)
		&& actual.composition==expected.composition,"fault changed a committed fluid FSI snapshot");
}
void CheckNoTrialLeakage(const iga::MovingImmersedTransientFlowFsiRuntime& runtime)
{
	Reject([&]{ (void)runtime.GetSurfaceTraction("patch"); });
	Reject([&]{ (void)runtime.TrialSurfaceTractionDiagnostics(); });
	Reject([&]{ (void)runtime.TrialFlowDiagnostics(); });
	Reject([&]{ (void)runtime.TrialMaterialKinematicsSnapshot(); });
	Reject([&]{ (void)runtime.TrialCompositionIdentitySha256(); });
	Reject([&]{ (void)runtime.ConservationDiagnostics(); });
}

iga::RawSurfaceTriangle Face(int a, int b, int c, std::uint32_t label)
{ iga::RawSurfaceTriangle value; value.indices={{a,b,c}}; value.boundary_id=label; return value; }

iga::MaterialSurfaceKinematics Initial()
{
	iga::RawSurfaceSoup soup;
	// A closed rectangular channel.  Every face is a 3x3 grid; the top label-7
	// patch has eight triangles, a clamped perimeter, and one free center node.
	std::vector<std::array<int,3>> lattice;
	auto node=[&](int x,int y,int z) {
		const std::array<int,3> key{{x,y,z}};
		for(std::size_t i=0;i<lattice.size();++i) if(lattice[i]==key) return static_cast<int>(i);
		lattice.push_back(key); soup.vertices.push_back({{.2+.3*x,.2+.3*y,.2+.3*z}});
		return static_cast<int>(lattice.size()-1);
	};
	auto add_face=[&](std::uint32_t label, auto coordinate) {
		int grid[3][3];
		for(int i=0;i<3;++i) for(int j=0;j<3;++j) { const auto p=coordinate(i,j); grid[i][j]=node(p[0],p[1],p[2]); }
		for(int i=0;i<2;++i) for(int j=0;j<2;++j) {
			const int a=grid[i][j],b=grid[i+1][j],c=grid[i][j+1],d=grid[i+1][j+1];
			auto add_outward=[&](int q0,int q1,int q2) {
				const auto& p0=soup.vertices[q0]; const auto& p1=soup.vertices[q1]; const auto& p2=soup.vertices[q2];
				const std::array<double,3> normal{{(p1[1]-p0[1])*(p2[2]-p0[2])-(p1[2]-p0[2])*(p2[1]-p0[1]),(p1[2]-p0[2])*(p2[0]-p0[0])-(p1[0]-p0[0])*(p2[2]-p0[2]),(p1[0]-p0[0])*(p2[1]-p0[1])-(p1[1]-p0[1])*(p2[0]-p0[0])}};
				const std::array<double,3> radial{{(p0[0]+p1[0]+p2[0])/3.-.5,(p0[1]+p1[1]+p2[1])/3.-.5,(p0[2]+p1[2]+p2[2])/3.-.5}};
				if(normal[0]*radial[0]+normal[1]*radial[1]+normal[2]*radial[2]<0.) std::swap(q1,q2);
				soup.triangles.push_back(Face(q0,q1,q2,label));
			};
			add_outward(a,b,d); add_outward(a,d,c);
		}
	};
	add_face(7,[](int x,int y){return std::array<int,3>{{x,y,2}};});
	add_face(9,[](int x,int y){return std::array<int,3>{{x,y,0}};});
	add_face(1,[](int y,int z){return std::array<int,3>{{0,y,z}};});
	add_face(2,[](int y,int z){return std::array<int,3>{{2,y,z}};});
	add_face(9,[](int x,int z){return std::array<int,3>{{x,0,z}};});
	add_face(9,[](int x,int z){return std::array<int,3>{{x,2,z}};});
	return iga::PrescribedSurfaceMotion({{0.0,soup},{1.0,soup}}).Evaluate(1.0,0.0,1.0);
}

iga::DistributedSurfaceInterface Structure(const std::string& ref)
{
	iga::DistributedSurfaceInterface v; v.id={"structure","membrane","patch"}; v.subsystem_id="membrane";
	v.boundary_labels={7}; v.reference_mesh_identity_sha256=ref;
	v.provides={iga::SurfaceFieldQuantity::Displacement,iga::SurfaceFieldQuantity::Velocity};
	v.requires={iga::SurfaceFieldQuantity::TractionOnStructure}; return v;
}
iga::DistributedSurfaceInterface Fluid(const std::string& ref)
{
	iga::DistributedSurfaceInterface v; v.id={"fluid","immersed","patch"}; v.subsystem_id="immersed";
	v.boundary_labels={7}; v.reference_mesh_identity_sha256=ref;
	v.provides={iga::SurfaceFieldQuantity::TractionOnStructure};
	v.requires={iga::SurfaceFieldQuantity::Displacement,iga::SurfaceFieldQuantity::Velocity}; return v;
}
iga::MaterialSurfacePatchMap Map(const iga::MaterialSurfaceKinematics& initial)
{
	iga::DistributedSurfaceLayout layout; layout.global_node_count=9; layout.partition_count=1; layout.partition_rank=0;
	layout.owned_global_node_ids={11,12,13,14,15,16,17,18,19};
	for(std::size_t i=0;i<9;++i) layout.reference_positions.push_back({layout.owned_global_node_ids[i],initial.ReferenceMaterialVerticesM()[i]});
	for(std::uint64_t x=0;x<2;++x) for(std::uint64_t y=0;y<2;++y) { const auto a=11+3*x+y,b=a+3,c=a+1,d=b+1; layout.reference_triangles.push_back({{a,b,d}}); layout.reference_triangles.push_back({{a,d,c}}); }
	layout.owned_reference_lumped_areas_m2={.03,.045,.015,.045,.09,.045,.015,.045,.03};
	std::vector<iga::MaterialSurfacePatchMap::GlobalToSourceVertex> vertices; for(std::size_t i=0;i<9;++i) vertices.push_back({layout.owned_global_node_ids[i],static_cast<std::uint32_t>(i)});
	const std::vector<std::uint32_t> triangles{0,1,2,3,4,5,6,7}; const std::vector<std::uint64_t> clamps{11,12,13,14,16,17,18,19};
	layout.reference_mesh_identity_sha256=iga::MaterialSurfacePatchMap::BuildReferenceIdentitySha256(initial,7,vertices,layout.reference_triangles,triangles,clamps);
	layout.layout_identity_sha256=iga::BuildDistributedSurfaceLayoutIdentitySha256(layout);
	return iga::MaterialSurfacePatchMap::Create(Structure(layout.reference_mesh_identity_sha256),layout,initial,7,vertices,triangles,clamps);
}
iga::MovingImmersedTransientFlowOptions Options()
{
	iga::MovingImmersedTransientFlowOptions v; v.grid={{{0,0,0}},{{1,1,1}},{{2,2,2}}};
	v.grid={{{0,0,0}},{{1,1,1}},{{3,3,3}}}; v.geometry.volume.max_depth=3; v.geometry.volume.max_nodes=v.geometry.volume.max_leaves=v.geometry.volume.max_points=1000000;
	v.flow.parameters={1.,1.,1.}; v.flow.wall_labels={7,9}; v.flow.ports={{"inlet",1,iga::ImmersedFlowPortControlMode::Pressure,.10},{"outlet",2,iga::ImmersedFlowPortControlMode::Pressure,0.}}; v.flow.ksp_relative_tolerance=1e-12; v.flow.lu_pivot_shift=1e-20;
	return v;
}
iga::SurfaceFieldStamp Stamp(const iga::DistributedSurfaceLayout& l, std::uint64_t iteration)
{
	iga::SurfaceFieldStamp s; s.time_s=2.; s.step=1; s.coupling_iteration=iteration;
	s.reference_mesh_identity_sha256=l.reference_mesh_identity_sha256; s.layout_identity_sha256=l.layout_identity_sha256;
	s.partition_identity_sha256=iga::BuildDistributedSurfacePartitionIdentitySha256(l);
	iga::Sha256 h; h.Append("structure-state",15); s.producer_state_identity_sha256=h.Hex(); return s;
}
iga::SurfaceFieldStampEnvelope Envelope(const iga::DistributedSurfaceLayout& l, std::uint64_t iteration)
{ auto s=Stamp(l,iteration); return iga::MakeSurfaceFieldStampEnvelope(s); }
iga::SurfaceKinematics Kinematics(const iga::MaterialSurfacePatchMap& map, std::uint64_t iteration)
{
	iga::SurfaceKinematics k; k.interface=map.Interface().id; k.stamp=Stamp(map.Layout(),iteration);
	k.displacement_m.assign(9,{{0,0,0}}); k.velocity_m_per_s.assign(9,{{0,0,0}});
	k.displacement_m[4]={{0,0,.006}}; k.velocity_m_per_s[4]={{0,0,.006}}; return k;
}
}

static_assert(!std::is_copy_constructible<iga::MovingImmersedTransientFlowFsiRuntime>::value, "single owner");
static_assert(!std::is_move_constructible<iga::MovingImmersedTransientFlowFsiRuntime>::value, "issued capabilities must not move");

int main(int argc, char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr); int status=0;
	try {
		const auto initial=Initial(); const auto map=Map(initial); const auto& layout=map.Layout(); const auto fluid=Fluid(map.ReferenceIdentitySha256());
		const iga::FsiCouplingEdge edge("wall",fluid.id,map.Interface().id,iga::FsiCouplingLaw::FluidStructureTractionKinematics);
		auto wrong_fluid=fluid; wrong_fluid.id.interface_id="wrong";
		Reject([&]{ iga::MovingImmersedTransientFlowFsiRuntime("fluid","immersed",edge,wrong_fluid,map.Interface(),layout,layout,initial,map,Options()); });
		auto wrong_structure=map.Interface(); wrong_structure.id.subsystem_id="wrong";
		Reject([&]{ iga::MovingImmersedTransientFlowFsiRuntime("fluid","immersed",edge,fluid,wrong_structure,layout,layout,initial,map,Options()); });
		iga::MovingImmersedTransientFlowFsiRuntime runtime("fluid","immersed",edge,fluid,map.Interface(),layout,layout,initial,map,Options());
		Reject([&]{ runtime.GetSurfaceTraction("patch"); });
		Reject([&]{ runtime.GetCommittedSurfaceTraction("patch",Stamp(layout,0)); });
		Reject([&]{ (void)runtime.ConservationDiagnostics(); });
		runtime.BeginMacroStep({1,1.,1.}); const auto k=Kinematics(map,0); const auto e=Envelope(layout,0);
		Reject([&]{ runtime.SetSurfaceKinematics("patch",k); });
		runtime.BeginCouplingIteration(0,k.stamp,e);
		auto stale=k; stale.stamp.coupling_iteration=1;
		Reject([&]{ runtime.SetSurfaceKinematics("patch",stale); });
		auto foreign=k; foreign.interface=fluid.id;
		Reject([&]{ runtime.SetSurfaceKinematics("patch",foreign); });
		runtime.SetSurfaceKinematics("patch",k);
		Reject([&]{ runtime.SetSurfaceKinematics("patch",k); });
		Reject([&]{ runtime.SetSurfaceKinematics("other",k); });
		runtime.SolveFluidTrial();
		const auto flow=runtime.TrialFlowDiagnostics();
		Require(flow.converged,"moving channel trial did not converge");
		Require(flow.identity_history_nodes==flow.active_nodes && flow.missing_history_nodes==0
			&& flow.committed_history_nodes+flow.extended_history_nodes==flow.identity_history_nodes
			&& !flow.history_hash_sha256.empty() && !flow.moving_map_identity_sha256.empty(),"moving channel state transfer coverage is incomplete");
		const auto wall=flow.wall_penalty.by_boundary_id.find(7);
		Require(wall!=flow.wall_penalty.by_boundary_id.end() && wall->second.selected_points>0
			&& wall->second.selected_area_m2>0.0 && flow.wall_penalty.maximum_eta>0.0
			&& runtime.MaxTrialSurfaceRelativeVelocityNormForTesting()>0.0,"moving channel Nitsche wall contribution or speed is absent");
		const auto traction=runtime.GetSurfaceTraction("patch"); iga::ValidateSurfaceFieldStampMatchesEnvelope(traction.stamp,e,layout);
		Require(!traction.projection_identity_sha256.empty(),"moving channel did not publish traction identity");
		const auto target=runtime.TrialMaterialKinematicsSnapshot();
		Require(target.SourceVerticesM()[4][2]>initial.SourceVerticesM()[4][2],"free patch center did not move");
		const auto conservation=runtime.ConservationDiagnostics();
		Require(std::isfinite(conservation.total_fluid_surface_outward_flow_m3_s)
			&& std::isfinite(conservation.normalized_moving_mass_defect)
			&& conservation.fluid_surface_outward_flow_by_boundary_label_m3_s.count(1)==1
			&& conservation.fluid_surface_outward_flow_by_boundary_label_m3_s.count(2)==1,"moving channel conservation diagnostics are incomplete");
		Require(conservation.normalized_moving_mass_defect<0.03
			&& std::isfinite(conservation.normalized_wall_relative_leakage)
			&& std::isfinite(conservation.normalized_discrete_moving_wall_continuity_defect),"moving channel conservation bounds are not credible");
		const auto audit=runtime.TrialSurfaceTractionDiagnostics();
		Require(Norm(audit.nodal_resultant_n)>0. || Norm(audit.nodal_moment_n_m)>0.,"moving channel produced trivial patch traction");
		for(int q=0;q<3;++q) Require(std::isfinite(audit.quadrature_resultant_n[q])&&std::isfinite(audit.nodal_resultant_n[q])&&std::isfinite(audit.quadrature_moment_n_m[q])&&std::isfinite(audit.nodal_moment_n_m[q])&&std::abs(audit.quadrature_resultant_n[q]-audit.nodal_resultant_n[q])<=1e-11&&std::abs(audit.quadrature_moment_n_m[q]-audit.nodal_moment_n_m[q])<=1e-11,"patch traction resultant/moment audit failed");
		const auto original_trial=TakeTrial(runtime);
		const auto committed_before_faults=TakeCommitted(runtime);
		const auto geometry=runtime.CommittedGeometry().GeometryIdentitySha256();
		runtime.RejectCouplingIteration();
		Require(runtime.CommittedGeometry().GeometryIdentitySha256()==geometry&&runtime.CommittedGlobalState().Index()==0,"reject changed committed fluid state");
		CheckNoTrialLeakage(runtime); CheckCommitted(runtime,committed_before_faults);
		runtime.BeginCouplingIteration(0,k.stamp,e); runtime.SetSurfaceKinematics("patch",k); runtime.SolveFluidTrial();
		const auto ordinary_retry=TakeTrial(runtime);
		Require(SameTrial(ordinary_retry,original_trial),"ordinary reject/retry changed exact fluid FSI publication");
		runtime.AbortStep();
		Require(runtime.CommittedGeometry().GeometryIdentitySha256()==geometry&&runtime.CommittedGlobalState().Index()==0,"abort changed committed fluid state");
		CheckNoTrialLeakage(runtime); CheckCommitted(runtime,committed_before_faults);
		runtime.BeginMacroStep({1,1.,1.}); runtime.BeginCouplingIteration(0,k.stamp,e); runtime.SetSurfaceKinematics("patch",k);
		runtime.FailNextLatePublicationForTesting(); Reject([&]{ runtime.SolveFluidTrial(); });
		Require(runtime.Lifecycle().Phase()==iga::FsiTrialPhase::InputReady && runtime.CommittedGeometry().GeometryIdentitySha256()==geometry
			&& runtime.CommittedGlobalState().Index()==0,"late publication failure changed lifecycle or numerical commit");
		CheckNoTrialLeakage(runtime); CheckCommitted(runtime,committed_before_faults);
		runtime.SolveFluidTrial();
		const auto late_retry=TakeTrial(runtime);
		Require(SameTrial(late_retry,original_trial),"late-publication failure retry changed exact fluid FSI publication");
		runtime.RejectCouplingIteration();
		CheckNoTrialLeakage(runtime); CheckCommitted(runtime,committed_before_faults);
		runtime.BeginCouplingIteration(0,k.stamp,e); runtime.SetSurfaceKinematics("patch",k); runtime.SolveFluidTrial();
		const auto prepare_trial=TakeTrial(runtime);
		Require(SameTrial(prepare_trial,original_trial),"late-failure reject/retry changed exact fluid FSI publication");
		runtime.FailNextPrepareForTesting(); Reject([&]{ runtime.PrepareCommitStep(); });
		CheckNoTrialLeakage(runtime); CheckCommitted(runtime,committed_before_faults);
		runtime.BeginCouplingIteration(0,k.stamp,e); runtime.SetSurfaceKinematics("patch",k); runtime.SolveFluidTrial();
		const auto prepare_retry=TakeTrial(runtime);
		Require(SameTrial(prepare_retry,original_trial),"prepare failure retry changed exact fluid FSI publication");
		runtime.PrepareCommitStep(); Require(runtime.CommittedGeometry().GeometryIdentitySha256()==geometry,"prepare changed committed geometry"); runtime.FinalizeCommitStep();
		const auto committed=runtime.GetCommittedSurfaceTraction("patch",prepare_retry.traction.stamp);
		Require(SameTraction(committed,prepare_retry.traction)&&runtime.CommittedGlobalState().Index()==1,"finalize did not expose exact committed traction");
		Require(SameTractionDiagnostics(runtime.CommittedSurfaceTractionDiagnostics(),prepare_retry.traction_diagnostics),"finalize did not expose exact committed traction diagnostics");
		const auto committed_conservation=runtime.ConservationDiagnostics();
		Require(SameConservation(committed_conservation,original_trial.conservation)
			&&runtime.CommittedCompositionIdentitySha256()==original_trial.composition,"commit did not retain transition conservation/composition");
		std::cout<<"channel patch_nodes="<<layout.global_node_count<<" patch_triangles="<<layout.reference_triangles.size()<<" converged=true open_flow="<<committed_conservation.open_port_outward_flow_m3_s<<" moving_mass="<<committed_conservation.normalized_moving_mass_defect<<" resultant="<<Norm(audit.nodal_resultant_n)<<"\n";
	} catch (const std::exception& error) { std::cerr<<error.what()<<"\n"; status=1; }
	PetscFinalize(); return status;
}
