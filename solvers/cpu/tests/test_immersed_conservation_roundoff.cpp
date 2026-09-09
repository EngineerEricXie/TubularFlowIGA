#include "MovingImmersedTransientFlowRuntime.hpp"
#include "PrescribedSurfaceMotion.hpp"
#include <iostream>

namespace {
void Require(bool condition, const char* message)
{
	if (!condition) throw std::runtime_error(message);
}
}

void MovingConservation()
{
	iga::RawSurfaceSoup soup;
	soup.vertices = {{{{0,0,0}},{{1,0,0}},{{1,1,0}},{{0,1,0}},{{0,0,1}},{{1,0,1}},{{1,1,1}},{{0,1,1}}}};
	const std::array<std::array<std::int64_t,3>,12> triangles{{{{0,2,1}},{{0,3,2}},{{4,5,6}},{{4,6,7}},
		{{0,1,5}},{{0,5,4}},{{1,2,6}},{{1,6,5}},{{2,3,7}},{{2,7,6}},{{3,0,4}},{{3,4,7}}}};
	for (std::size_t i = 0; i < triangles.size(); ++i) {
		iga::RawSurfaceTriangle triangle; triangle.indices = triangles[i]; triangle.boundary_id = i < 2 ? 8 : i < 4 ? 9 : 7;
		soup.triangles.push_back(triangle);
	}
	iga::PrescribedSurfaceMotion motion({{0,soup},{1,soup}});
	iga::MovingImmersedTransientFlowOptions options;
	options.grid = {{{0,0,0}},{{1,1,1}},{{4,1,1}}}; options.geometry.volume.max_depth = 2;
	options.geometry.volume_storage = iga::CutCellVolumeQuadratureStorageMode::Compact;
	options.flow.parameters = {1,1,0}; options.flow.wall_labels = {7};
	options.flow.ports = {{"outlet",9,iga::ImmersedFlowPortControlMode::FlowRate,.0001},
		{"inlet",8,iga::ImmersedFlowPortControlMode::FlowRate,-.0001}};
	options.flow.nonlinear_absolute_tolerance = 1e-14; options.flow.nonlinear_relative_tolerance = 1e-11;
	options.flow.ksp_relative_tolerance = 1e-11; options.flow.lu_pivot_shift = 1e-12;
	iga::MovingImmersedTransientFlowRuntime runtime(motion.Evaluate(0,0,1),options);
	runtime.BeginTrial(motion.Evaluate(.125,0,.125),1,.125);
	Require(runtime.SolveTrial(),"moving owner stationary solve failed");
	const auto conservation = runtime.ConservationDiagnostics();
	Require(std::abs(conservation.open_port_outward_flow_m3_s) < 1e-12,"physical port balance failed");
	Require(std::abs(conservation.divergence_theorem_defect_m3_s) < 1e-12,"physical divergence balance failed");
	const auto reject = [&] {
		bool rejected = false;
		try { runtime.ConservationDiagnostics(); } catch (const std::logic_error&) { rejected = true; }
		Require(rejected,"injected finite continuity defect accepted");
	};
	{ iga::ImmersedTransientFlowRuntime::FiniteDiscreteMovingWallContinuityDefectScopeForTesting injected; reject(); }
	{ iga::MovingImmersedTransientFlowRuntime::FiniteDiscreteMovingWallContinuityDefectScopeForTesting injected; reject(); }
	runtime.PrepareCommit(); runtime.FinalizeCommit();
	Require(runtime.ConservationDiagnostics().total_fluid_surface_outward_flow_m3_s == conservation.total_fluid_surface_outward_flow_m3_s,
		"committed conservation record changed");
}

int main(int argc,char** argv)
{
	if (PetscInitialize(&argc,&argv,nullptr,nullptr)) return 1;
	int status = 0;
	try {
		using iga::immersed_transient_detail::MovingWallContinuityIdentityReconciles;
		using iga::immersed_transient_detail::MovingWallContinuityIdentityRoundoffTolerance;
		// The total visits a wall contribution between equal/opposite port
		// contributions. The independently accumulated port sum is exactly zero.
		double total = 0, port = 0, wall = 0, absolute = 0;
		const std::array<double,3> contributions{{1e-4,1e-20,-1e-4}};
		for (std::size_t i = 0; i < contributions.size(); ++i) {
			total += contributions[i]; absolute += std::abs(contributions[i]);
			if (i == 1) wall += contributions[i]; else port += contributions[i];
		}
		Require(!MovingWallContinuityIdentityReconciles(port,wall,0,port,wall,total),"fixture does not expose lost cancellation scale");
		Require(MovingWallContinuityIdentityReconciles(port,wall,0,port,wall,total,absolute,3),"raw quadrature cancellation rejected");
		const double tolerance = MovingWallContinuityIdentityRoundoffTolerance(port,wall,0,port,wall,total,absolute,3);
		Require(tolerance > 0 && tolerance < 1e-16,"roundoff tolerance became a physical gate");
		Require(!MovingWallContinuityIdentityReconciles(port,wall,0,port,wall,total+2*tolerance,absolute,3),"inconsistent flux accepted");
		Require(MovingWallContinuityIdentityReconciles(1e12,-1e12,-1e12,0,0,1e-2),"material aggregate cancellation rejected");
		Require(!MovingWallContinuityIdentityReconciles(1e12,-1e12,-1e12,0,0,1e-1),"material aggregate inconsistency accepted");
		Require(!MovingWallContinuityIdentityReconciles(0,0,0,0,0,0,-1,3),"negative absolute sum accepted");
		Require(!MovingWallContinuityIdentityReconciles(0,0,0,0,0,0,std::numeric_limits<double>::infinity(),3),"infinite absolute sum accepted");
		Require(!MovingWallContinuityIdentityReconciles(0,0,0,0,0,0,absolute,std::numeric_limits<std::uint64_t>::max()),"unbounded accumulation accepted");
		MovingConservation();
		std::cout << "immersed_conservation_roundoff=passed defect=" << port+wall-total << " tolerance=" << tolerance << '\n';
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n'; status = 1;
	}
	PetscFinalize(); return status;
}
