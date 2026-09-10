#include "MovingImmersedTransientFlowRuntime.hpp"
#include "PrescribedSurfaceMotion.hpp"
#include <iostream>

namespace {
iga::RawSurfaceSoup Soup()
{
	iga::RawSurfaceSoup soup;
	soup.vertices = {{{{0,0,0}},{{1,0,0}},{{1,1,0}},{{0,1,0}},{{0,0,1}},{{1,0,1}},{{1,1,1}},{{0,1,1}}}};
	const std::array<std::array<std::int64_t,3>,12> triangles{{{{0,2,1}},{{0,3,2}},{{4,5,6}},{{4,6,7}},
		{{0,1,5}},{{0,5,4}},{{1,2,6}},{{1,6,5}},{{2,3,7}},{{2,7,6}},{{3,0,4}},{{3,4,7}}}};
	for (std::size_t i = 0; i < triangles.size(); ++i) {
		iga::RawSurfaceTriangle triangle; triangle.indices = triangles[i]; triangle.boundary_id = i < 2 ? 8 : i < 4 ? 9 : 7;
		soup.triangles.push_back(triangle);
	}
	return soup;
}
void Require(bool value,const char* message)
{
	if(!value)throw std::runtime_error(message);
}
void Run(iga::CutCellVolumeQuadratureStorageMode storage)
{
	auto source=Soup();for(auto& x:source.vertices)for(double& value:x)value+=.1;
	auto target=source;for(auto& x:target.vertices)x[0]+=.74;
	iga::PrescribedSurfaceMotion motion({{0,source},{.25,target}});
	iga::MovingImmersedTransientFlowOptions options;
	options.grid={{{0,0,0}},{{2.4,1.2,1.2}},{{8,3,3}}};
	options.geometry.volume.max_depth=2;options.geometry.volume_storage=storage;options.geometry.volume_fitting.emplace();
	options.extension_layers=2;options.flow.parameters={1,.1,0};options.flow.wall_labels={7};
	options.flow.ports={{"inlet",8,iga::ImmersedFlowPortControlMode::FlowRate,0},{"outlet",9,iga::ImmersedFlowPortControlMode::FlowRate,0}};
	iga::MovingImmersedTransientFlowRuntime runtime(motion.Evaluate(0,0,.125),options);
	const auto& layout=runtime.CommittedLayout();std::vector<std::array<double,4>> fields(layout.NodeIds().size(),{{2.96,0,0,0}});
	runtime.InitializeCommittedGlobalState(iga::ImmersedGlobalFlowState(0,0,layout,fields,{0,0},true,0));
	runtime.BeginTrial(motion.Evaluate(.125,0,.125),1,.125);
	std::size_t outside=0;double maximum=0;
	const auto& geometry=runtime.TrialGeometry();
	std::vector<std::array<double,3>> before;
	for(std::uint64_t cell=0;cell<geometry.Domain().Cells().size();++cell) {
		const auto& d=geometry.Volume().Cell(cell).diagnostics;
		before.push_back({{d.lower_reference_volume,d.estimated_reference_volume,d.upper_reference_volume}});
		const double violation=std::max(d.lower_reference_volume-d.estimated_reference_volume,d.estimated_reference_volume-d.upper_reference_volume);
		if(violation>0) { ++outside;maximum=std::max(maximum,violation); }
	}
	runtime.Assemble();
	const auto& d=runtime.TrialDiagnostics().wall_penalty;
	Require(d.fraction_ordering_tolerance>0&&std::isfinite(d.fraction_ordering_tolerance),"wall fraction ordering allowance missing");
	Require(d.fraction_lower<=d.fraction_estimate+d.fraction_ordering_tolerance
		&&d.fraction_estimate<=d.fraction_upper+d.fraction_ordering_tolerance,"aggregated fractions exceed ordering allowance");
	for(std::uint64_t cell=0;cell<geometry.Domain().Cells().size();++cell) {
		const auto& d=geometry.Volume().Cell(cell).diagnostics;
		Require(before[cell]==std::array<double,3>{{d.lower_reference_volume,d.estimated_reference_volume,d.upper_reference_volume}},"wall diagnostics changed raw quadrature bounds or estimate");
	}
	std::cout<<"fitted_wall_diagnostics storage="<<static_cast<int>(storage)<<" outside_cells="<<outside
		<<" maximum_interval_gap="<<maximum<<" aggregate_allowance="<<d.fraction_ordering_tolerance<<" passed\n";
}
}
int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);int status=0;
	try { Run(iga::CutCellVolumeQuadratureStorageMode::Expanded);Run(iga::CutCellVolumeQuadratureStorageMode::Compact); }
	catch(const std::exception& error) { std::cerr<<error.what()<<'\n';status=1; }
	PetscFinalize();return status;
}
