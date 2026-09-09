#include "CompliantChannelFsiFixture.hpp"
#include "MovingImmersedTransientFlowFsiRuntime.hpp"
#include "PretensionedMembraneFsiRuntime.hpp"
#include "StrongFluidStructureCoupling.hpp"
#include "ReferenceStateOutput.hpp"

#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>

namespace {

void Require(bool value, const char* message)
{ if (!value) throw std::runtime_error(message); }
double Norm(const std::array<double,3>& value)
{ return std::sqrt(value[0]*value[0]+value[1]*value[1]+value[2]*value[2]); }
bool Finite(double value) { return std::isfinite(value); }

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr); int status=0;
	int rank=0, ranks=1;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank); MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	iga::CurrentPhaseProfile().EnableFromEnvironment();
	try {
		auto reference = iga::test::ReferenceOutputFromArguments(argc, argv, "fsi", ranks);
		const auto initial=iga::compliant_channel_fixture::InitialMaterial();
		const auto map=iga::compliant_channel_fixture::PatchMap(initial);
		const auto& layout=map.Layout();
		const auto fluid_surface=iga::compliant_channel_fixture::Fluid(map.ReferenceIdentitySha256());
		const iga::FsiCouplingEdge edge("compliant-channel-wall",fluid_surface.id,map.Interface().id,
			iga::FsiCouplingLaw::FluidStructureTractionKinematics);
		// SI material: positive inertia/damping, 20 N/m pretension and a modest
		// 5 N/m^3 foundation give measurable (sub-mm) feedback in this .6 m patch.
		const auto membrane_material=iga::compliant_channel_fixture::MembraneMaterial();
		const auto clamped=iga::compliant_channel_fixture::ClampedGlobalNodeIds();
		iga::MovingImmersedTransientFlowFsiRuntime fluid("fluid","immersed",edge,fluid_surface,
			map.Interface(),layout,layout,initial,map,iga::compliant_channel_fixture::FlowOptions());
		iga::PretensionedMembraneFsiRuntime structure("structure","membrane",edge,map.Interface(),
			fluid_surface,layout,layout,membrane_material,clamped);
		const auto coupling_options=iga::compliant_channel_fixture::CouplingOptions();
		iga::StrongFluidStructureCoupling<iga::MovingImmersedTransientFlowFsiRuntime,
			iga::PretensionedMembraneFsiRuntime> coordinator(fluid,structure,edge,layout,coupling_options);
		const auto result=coordinator.Execute({1,1.0,.05});

		Require(result.converged && result.status=="converged" && result.iterations>=2,
			"real compliant-channel strong coordinator did not converge after two-way iterations");
		Require(result.history.size()==result.iterations && result.history.front().area_weighted_rms_residual_m>0.0,
			"predictor response is trivial or result history is incomplete");
		Require(result.history.front().area_weighted_rms_residual_m>result.history.front().convergence_threshold_m,
			"predictor passed a loose convergence tolerance");
		Require(result.history.back().area_weighted_rms_residual_m<=result.history.back().convergence_threshold_m,
			"accepted membrane residual does not meet the exact strong-coupling threshold");
		Require(result.history.back().area_weighted_rms_residual_m<result.history.front().area_weighted_rms_residual_m,
			"strong coupling did not reduce the membrane residual");
		bool aitken_applied=false;
		std::uint64_t total_newton=0,total_ksp=0;
		for(const auto& d:result.history) {
			Require(Finite(d.area_weighted_rms_residual_m)&&Finite(d.convergence_threshold_m),"nonfinite strong residual diagnostic");
			if(d.aitken_proposal_applied) {
				aitken_applied=true; Require(d.relaxation && *d.relaxation>=.1 && *d.relaxation<=.8,
					"Aitken factor escaped configured bounds");
			}
			Require(d.fluid_nonlinear_iterations && d.fluid_ksp_iterations && d.fluid_residual_norm && d.fluid_linear_relative_residual,
				"production coordinator omitted fluid iteration diagnostics");
			Require(Finite(*d.fluid_residual_norm)&&Finite(*d.fluid_linear_relative_residual),"nonfinite trial Newton/KSP diagnostics");
			total_newton+=*d.fluid_nonlinear_iterations; total_ksp+=*d.fluid_ksp_iterations;
			std::cout<<"fsi_iteration="<<d.iteration<<" residual_rms_m="<<d.area_weighted_rms_residual_m
				<<" threshold_m="<<d.convergence_threshold_m<<" aitken="
				<<(d.relaxation ? *d.relaxation : 0.0)<<" applied="<<d.aitken_proposal_applied
				<<" newton="<<*d.fluid_nonlinear_iterations<<" ksp="<<*d.fluid_ksp_iterations<<"\n";
		}
		Require(aitken_applied,"strong iteration never applied weighted Aitken relaxation");

		const auto membrane_state=structure.CommittedStateSnapshot();
		Require(membrane_state.displacement_m.size()==9 && membrane_state.velocity_m_per_s.size()==9,
			"committed membrane layout changed");
		Require(membrane_state.displacement_m[4]>0.0 && membrane_state.velocity_m_per_s[4]>0.0,
			"accepted free-center membrane response is not outward");
		for(std::size_t i=0;i<membrane_state.displacement_m.size();++i) if(i!=4)
			Require(membrane_state.displacement_m[i]==0.0 && membrane_state.velocity_m_per_s[i]==0.0,
				"clamped membrane perimeter moved");
		const auto& flow=fluid.CommittedFlowDiagnostics();
		Require(flow.converged && !flow.trial_active && flow.committed && flow.identity_history_nodes==flow.active_nodes
			&& flow.missing_history_nodes==0 && !flow.history_hash_sha256.empty() && !flow.layout_hash_sha256.empty(),
			"committed fluid solve/state history/layout identity is incomplete");
		Require(Finite(flow.residual_norm)&&Finite(flow.true_linear_relative_residual)&&flow.nonlinear_iterations>0,
			"committed fluid Newton/KSP diagnostics are invalid");
		const auto wall=flow.wall_penalty.by_boundary_id.find(7);
		Require(wall!=flow.wall_penalty.by_boundary_id.end() && wall->second.selected_points>0
			&& wall->second.selected_area_m2>0.0 && flow.wall_penalty.maximum_eta>0.0,
			"moving-wall Nitsche diagnostics are incomplete");
		// The trial is finalized, so obtain the wall speed from the committed
		// kinematics rather than retaining a trial-only diagnostic.
		const auto full=fluid.CommittedMaterialKinematicsSnapshot();
		Require(std::abs(full.SourceVerticesM()[4][2]-initial.SourceVerticesM()[4][2])>0.0,
			"committed fluid geometry did not receive membrane motion");
		bool moving_wall=false;
		for(std::uint32_t triangle=0;triangle<full.Surface().Triangles().size();++triangle)
			if(full.Surface().Triangles()[triangle].boundary_id==7 && Norm(full.WallVelocity(triangle,{{1./3.,1./3.,1./3.}}))>0.0) moving_wall=true;
		Require(moving_wall,"committed top patch has zero wall velocity");
		const auto traction=fluid.CommittedSurfaceTractionSnapshot();
		Require(traction && !traction->projection_identity_sha256.empty(),"committed fluid traction identity is absent");
		const auto audit=fluid.CommittedSurfaceTractionDiagnostics();
		// Label 7 is the cube's top patch, whose fixture orientation is outward +z.
		Require(audit.nodal_resultant_n[2]>0.0,"coupled top-patch normal load is not outward");
		for(int q=0;q<3;++q) Require(Finite(audit.quadrature_resultant_n[q])&&Finite(audit.nodal_resultant_n[q])
			&&Finite(audit.quadrature_moment_n_m[q])&&Finite(audit.nodal_moment_n_m[q])
			&&std::abs(audit.quadrature_resultant_n[q]-audit.nodal_resultant_n[q])<=1.e-11
			&&std::abs(audit.quadrature_moment_n_m[q]-audit.nodal_moment_n_m[q])<=1.e-11,
			"fluid traction resultant/moment projection audit failed");
		const auto conservation=fluid.ConservationDiagnostics();
		// The accepted two-way response retains the PR8.3b conservation gates.
		Require(Finite(conservation.normalized_moving_mass_defect)&&Finite(conservation.normalized_wall_relative_leakage)
			&&Finite(conservation.normalized_discrete_moving_wall_continuity_defect)
			&&conservation.normalized_moving_mass_defect<.03
			&&conservation.normalized_wall_relative_leakage<.03
			&&conservation.normalized_discrete_moving_wall_continuity_defect<1.e-8,
			"moving-domain conservation exceeds compliant-channel bounds");
		Require(fluid.CommittedGlobalState().Index()==1 && fluid.CommittedGlobalState().TimeS()==1.05
			&&fluid.MovingDiagnostics().idle && !fluid.MovingDiagnostics().trial_active
			&&fluid.CommittedCompositionIdentitySha256()==result.final_committed_fluid_composition_identity_sha256
			&&structure.CommittedStateIdentitySha256()==result.final_committed_structure_state_identity_sha256,
			"fluid/structure did not finalize the same committed FSI epoch");
		std::cout<<"compliant_channel_fsi converged=true iterations="<<result.iterations
			<<" center_displacement_m="<<membrane_state.displacement_m[4]
			<<" center_velocity_m_s="<<membrane_state.velocity_m_per_s[4]
			<<" total_newton="<<total_newton<<" total_ksp="<<total_ksp
			<<" open_flow_m3_s="<<conservation.open_port_outward_flow_m3_s
			<<" moving_mass="<<conservation.normalized_moving_mass_defect
			<<" wall_leakage="<<conservation.normalized_wall_relative_leakage
			<<" continuity="<<conservation.normalized_discrete_moving_wall_continuity_defect
			<<" traction_resultant_n="<<Norm(audit.nodal_resultant_n)<<"\n";
		if (reference) {
			iga::PhaseScope output_phase(iga::ProfilePhase::Output);
			const auto& state = fluid.CommittedGlobalState();
			reference->Add("fluid_velocity", "m/s", state.NodeIds(), 3,
				[&](std::size_t row, std::size_t column) { return state.Coefficients().at(row).at(column); });
			reference->Add("fluid_pressure", "pa", state.NodeIds(), 1,
				[&](std::size_t row, std::size_t) { return state.Coefficients().at(row).at(3); });
			if (!state.PortIds().empty()) reference->Add("controller_pressure", "pa", state.PortIds(), 1,
				[&](std::size_t row, std::size_t) { return state.PortMultipliers().at(row); });
			if (state.HasGaugeMultiplier()) reference->Add("gauge_multiplier", "1/s", std::vector<int>{0}, 1,
				[&](std::size_t, std::size_t) { return state.GaugeMultiplier(); });
			reference->Add("membrane_displacement", "m", layout.owned_global_node_ids, 1,
				[&](std::size_t row, std::size_t) { return membrane_state.displacement_m.at(row); });
			reference->Add("membrane_velocity", "m/s", layout.owned_global_node_ids, 1,
				[&](std::size_t row, std::size_t) { return membrane_state.velocity_m_per_s.at(row); });
			reference->Add("surface_traction", "pa", layout.owned_global_node_ids, 3,
				[&](std::size_t row, std::size_t column) { return traction->traction_on_structure_pa.at(row).at(column); });
			reference->Add("surface_force", "n", layout.owned_global_node_ids, 3,
				[&](std::size_t row, std::size_t column) { return traction->consistent_nodal_force_n.at(row).at(column); });
			std::vector<std::uint64_t> vertices(full.SourceVerticesM().size());
			for (std::size_t i = 0; i < vertices.size(); ++i) vertices[i] = i;
			reference->Add("material_position", "m", vertices, 3,
				[&](std::size_t row, std::size_t column) { return full.SourceVerticesM().at(row).at(column); });
			reference->Add("material_displacement", "m", vertices, 3,
				[&](std::size_t row, std::size_t column) {
					return full.SourceVerticesM().at(row).at(column)-full.ReferenceMaterialVerticesM().at(row).at(column);
				});
			reference->Add("material_velocity", "m/s", vertices, 3,
				[&](std::size_t row, std::size_t column) { return full.SourceVertexVelocitiesMPerS().at(row).at(column); });
			reference->Finish(state.TimeS(), state.Index());
		}
	} catch(const std::exception& error) { std::cerr<<error.what()<<"\n"; status=1; }
	iga::CurrentPhaseProfile().Write(std::cout,rank,ranks,status);
	PetscFinalize(); return status;
}
