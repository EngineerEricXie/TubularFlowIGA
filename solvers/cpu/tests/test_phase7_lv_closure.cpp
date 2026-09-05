#include "IdealizedLeftVentricleFixture.hpp"
#include "MovingImmersedFlowSnapshotCapture.hpp"
#include "MovingImmersedFlowSnapshotPublisher.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs=std::filesystem;
namespace {
constexpr double kFlowM3S=1.0e-6; // low/moderate-Re numerical inlet target
constexpr double kControllerTolerance=1.0e-8;
constexpr double kControllerAcceptanceRelativeTolerance=1.0e-6;
constexpr double kLinearTolerance=1.0e-8;
constexpr double kWallRelativeTolerance=.10;
constexpr double kWallRelativeReleaseTolerance=.08;
constexpr double kConservationTolerance=.03;
constexpr std::uint64_t kCycleSteps=16;
void Require(bool value,const char* text){if(!value)throw std::runtime_error(text);}

double CanonicalTime(std::uint64_t index)
{
	Require(index<=kCycleSteps,"LV canonical time index is outside the cycle");
	return iga::IdealizedLeftVentricleFixture::PeriodS*static_cast<double>(index)/static_cast<double>(kCycleSteps);
}
struct TimeLatticeStep { double source,target,dt; };
TimeLatticeStep NextTimeLatticeStep(std::uint64_t index,double source)
{
	Require(index>0&&index<=kCycleSteps,"LV time-lattice transition index is invalid");
	const double expected_source=CanonicalTime(index-1),target=CanonicalTime(index),dt=target-source;
	Require(source==expected_source,"LV committed time is not the canonical predecessor");
	Require(std::isfinite(source)&&std::isfinite(target)&&std::isfinite(dt)&&dt>0.,"LV time-lattice step is not finite and positive");
	Require(iga::CheckedTransientTargetTime(source,dt)==target,"LV time-lattice target is not the exact transient successor");
	const double nominal_dt=iga::IdealizedLeftVentricleFixture::PeriodS/static_cast<double>(kCycleSteps);
	const double ulp_tolerance=8.*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(source),std::abs(target),std::abs(nominal_dt)});
	Require(std::abs(dt-nominal_dt)<=ulp_tolerance,"LV time-lattice step differs from nominal dt by more than ULP noise");
	return {source,target,dt};
}
iga::PrescribedSurfaceMotion::Evaluation LatticeEvaluation(const iga::PrescribedSurfaceMotion& motion,std::uint64_t index)
{
	Require(index<=kCycleSteps,"LV lattice evaluation index is outside the cycle");
	if(index==0)return motion.Evaluate(CanonicalTime(0),CanonicalTime(0),CanonicalTime(1));
	const auto step=NextTimeLatticeStep(index,CanonicalTime(index-1));
	return motion.Evaluate(step.target,step.source,step.target);
}

iga::MovingCutGeometryOptions GeometryOptions()
{
	iga::MovingCutGeometryOptions value;
	value.volume.max_depth=2;
	value.volume.empty_rule_rescue_max_depth=9;
	value.volume.max_nodes=value.volume.max_leaves=value.volume.max_points=value.volume.max_records=value.volume.max_logical_points=2000000;
	value.volume.max_retained_bytes=200000000;
	return value;
}
iga::MovingImmersedTransientFlowOptions Options()
{
	iga::MovingImmersedTransientFlowOptions value;
	value.grid={{{-.033,-.033,-.007}},{{.031,.031,.067}},{{6,6,7}}};
	value.geometry=GeometryOptions(); value.extension_layers=3;
	value.extension.max_unknowns=4096; value.extension.max_dense_bytes=128u*1024u*1024u;
	value.flow.parameters={1050.0,0.012,0.0}; // kg/m3, Pa s: numerical Re is intentionally modest
	value.flow.wall_labels={0};
	value.flow.wall_inertial_gamma0=12.0;
	value.flow.ports={{"inlet",1,iga::ImmersedFlowPortControlMode::FlowRate,-kFlowM3S},{"outlet",2,iga::ImmersedFlowPortControlMode::Pressure,0.0}};
	value.flow.ksp_relative_tolerance=1e-10; value.flow.nonlinear_relative_tolerance=1e-9; value.flow.nonlinear_absolute_tolerance=1e-11;
	value.flow.flow_controller_relative_tolerance=kControllerTolerance; value.flow.flow_controller_absolute_tolerance_m3_s=1e-13; value.flow.flow_controller_reference_flow_m3_s=kFlowM3S;
	return value;
}
std::array<unsigned,3> NodeIndex(std::int32_t id,const iga::CubicCartesianGridSpec& grid)
{const unsigned nx=grid.cells[0]+3,ny=grid.cells[1]+3;return {{unsigned(id)%nx,(unsigned(id)/nx)%ny,unsigned(id)/(nx*ny)}};}
double Knot(const iga::CubicCartesianGridSpec& grid,int d,unsigned i)
{if(i<=3)return grid.lower_m[d];if(i>=grid.cells[d]+3)return grid.upper_m[d];return grid.lower_m[d]+(grid.upper_m[d]-grid.lower_m[d])*(i-3)/grid.cells[d];}
std::vector<std::array<double,4>> RotationalSeed(const iga::ImmersedActiveLayout& layout,const iga::CubicCartesianGridSpec& grid)
{
	std::vector<std::array<double,4>> result;result.reserve(layout.NodeIds().size()); constexpr double omega=.35;
	for(const auto id:layout.NodeIds()){const auto p=NodeIndex(id,grid);const double x=(Knot(grid,0,p[0]+1)+Knot(grid,0,p[0]+2)+Knot(grid,0,p[0]+3))/3.;const double y=(Knot(grid,1,p[1]+1)+Knot(grid,1,p[1]+2)+Knot(grid,1,p[1]+3))/3.;result.push_back({{-omega*y,omega*x,0.,0.}});}return result;
}
iga::ImmersedGlobalFlowState InitialRotationalState(const iga::MovingImmersedTransientFlowRuntime& runtime,const iga::CubicCartesianGridSpec& grid)
{
	const auto& layout=runtime.CommittedLayout();const auto& baseline=runtime.CommittedGlobalState();
	Require(layout.Valid(),"LV initial layout is invalid");
	Require(layout.GeometryIdentity()==runtime.CommittedGeometry().GeometryIdentitySha256(),"LV initial layout geometry identity is stale");
	Require(baseline.GeometryIdentity()==layout.GeometryIdentity()&&baseline.TimeS()==runtime.CommittedGeometry().Evaluation().EvaluatedTimeS()&&baseline.Index()==0,"LV initial baseline geometry/time/index differs from committed epoch");
	Require(baseline.NodeIds()==layout.NodeIds()&&baseline.PortIds()==layout.PortIds()&&baseline.HasGaugeMultiplier()==layout.HasGaugeRow(),"LV initial baseline does not exactly cover committed layout");
	std::vector<std::uint64_t> expected_controller_ids;
	for(const auto& port:runtime.ConfiguredPorts())if(port.control_mode==iga::ImmersedFlowPortControlMode::FlowRate)expected_controller_ids.push_back(static_cast<std::uint64_t>(port.boundary_label));
	std::sort(expected_controller_ids.begin(),expected_controller_ids.end());
	Require(expected_controller_ids==layout.PortIds(),"LV initial flow-controller IDs/order differ from committed layout");
	auto coefficients=RotationalSeed(layout,grid);
	std::vector<double> controller_multipliers(layout.PortIds().size(),0.0);
	Require(coefficients.size()==layout.NodeIds().size()&&controller_multipliers.size()==layout.PortIds().size(),"LV initial field/controller sizes do not match committed layout");
	return iga::ImmersedGlobalFlowState(baseline.TimeS(),baseline.Index(),layout,std::move(coefficients),std::move(controller_multipliers),layout.HasGaugeRow(),0.0);
}
iga::MovingImmersedFlowSnapshotPublicationIdentity Identity(const iga::MovingImmersedTransientFlowRuntime& r)
{return {r.CommittedGeometry().GeometryIdentitySha256(),r.CommittedGeometry().PublicationIdentitySha256(),r.CommittedLayout().HashSha256(),r.CommittedGlobalState().HashSha256(),r.CommittedDiagnostics().committed_state_hash_sha256};}
struct WallQuadratureMetrics {
	double material_rms_m_per_s=0.,relative_rms_m_per_s=0.,normal_relative_rms_m_per_s=0.,tangential_relative_rms_m_per_s=0.,area_m2=0.;
};
WallQuadratureMetrics IntegrateWallQuadrature(const iga::MovingImmersedTransientFlowRuntime& runtime)
{
	const auto& geometry=runtime.CommittedGeometry();const auto& layout=runtime.CommittedLayout();const auto& state=runtime.CommittedGlobalState();long double material_l2=0.,relative_l2=0.,normal_l2=0.,tangential_l2=0.,area=0.;
	for(std::uint64_t id=0;id<geometry.Surface().Cells().size();++id){const auto element=geometry.Domain().Background().MaterializeElement(id);const auto& q=geometry.Surface().UsableRule(geometry.Domain(),id);const auto& p=geometry.Surface().UsableProvenance(geometry.Domain(),id);for(std::size_t i=0;i<q.Points().size();++i){const auto& point=q.Points()[i];if(point.boundary_id!=0)continue;const auto basis=iga::EvaluateBasis(element,point.parametric[0],point.parametric[1],point.parametric[2],false);std::array<double,3> u{};for(std::size_t a=0;a<element.connectivity.size();++a){const auto& coefficient=state.Coefficients()[layout.LocalNode(element.connectivity[a])];for(std::size_t d=0;d<3;++d)u[d]+=basis.value[a]*coefficient[d];}const auto w=geometry.Evaluation().WallVelocity(p[i].canonical_triangle,p[i].canonical_barycentric);const std::array<double,3> difference{{u[0]-w[0],u[1]-w[1],u[2]-w[2]}};const double relative_norm=std::sqrt(difference[0]*difference[0]+difference[1]*difference[1]+difference[2]*difference[2]),normal=difference[0]*point.normal[0]+difference[1]*point.normal[1]+difference[2]*point.normal[2],tangential_squared=std::max(0.,relative_norm*relative_norm-normal*normal),material_squared=w[0]*w[0]+w[1]*w[1]+w[2]*w[2];Require(std::isfinite(point.weight)&&point.weight>0.&&std::isfinite(relative_norm)&&std::isfinite(normal)&&std::isfinite(tangential_squared)&&std::isfinite(material_squared),"LV wall quadrature value is invalid");area+=point.weight;material_l2+=point.weight*material_squared;relative_l2+=point.weight*relative_norm*relative_norm;normal_l2+=point.weight*normal*normal;tangential_l2+=point.weight*tangential_squared;}}
	Require(std::isfinite(static_cast<double>(area))&&area>0.,"LV wall quadrature area is invalid");const double inverse_area=1./static_cast<double>(area);WallQuadratureMetrics result;result.area_m2=static_cast<double>(area);result.material_rms_m_per_s=std::sqrt(static_cast<double>(material_l2)*inverse_area);result.relative_rms_m_per_s=std::sqrt(static_cast<double>(relative_l2)*inverse_area);result.normal_relative_rms_m_per_s=std::sqrt(static_cast<double>(normal_l2)*inverse_area);result.tangential_relative_rms_m_per_s=std::sqrt(static_cast<double>(tangential_l2)*inverse_area);return result;
}
void GeometryOnly(const iga::PrescribedSurfaceMotion& motion)
{
	const auto& frames=motion.Frames(); Require(frames.size()==3,"LV fixture frame count changed");
	Require(iga::IdealizedLeftVentricleFixture::BitwiseEqualVertices(frames.front().surface,frames.back().surface),"LV source vertices are not exactly periodic");
	const auto ed=iga::ClosedTriangulatedSurface::Build(frames.front().surface),es=iga::ClosedTriangulatedSurface::Build(frames[1].surface),final=iga::ClosedTriangulatedSurface::Build(frames.back().surface);
	Require(ed.Diagnostics().component_count==1&&ed.Diagnostics().boundary_id_histogram.at(0)==84&&ed.Diagnostics().boundary_id_histogram.at(1)==6&&ed.Diagnostics().boundary_id_histogram.at(2)==6,"LV topology or labels changed");
	Require(ed.Diagnostics().canonical_sha256==final.Diagnostics().canonical_sha256,"LV source canonical hash is not periodic");
	const double ratio=es.Diagnostics().volume_m3/ed.Diagnostics().volume_m3;
	Require(std::abs(ratio-iga::IdealizedLeftVentricleFixture::EndSystolicVolumeRatio)<=5e-13*iga::IdealizedLeftVentricleFixture::EndSystolicVolumeRatio,"LV ES/ED triangulated-volume ratio changed");
}
enum class RunMode { FullCycle, GeometryOnly, GeometrySweep, OneStep, TimeLattice };
struct CommandLine { RunMode mode=RunMode::FullCycle; fs::path output_root; bool retain_output=false; };
CommandLine ParseCommandLine(int argc,char** argv)
{
	const auto usage=[] { throw std::invalid_argument("usage: phase7_lv_closure_test [output-directory] | --geometry-only | --geometry-sweep | --one-step [output-directory] | --time-lattice"); };
	if(argc==1)return {};
	const std::string first(argv[1]);
	if(first=="--geometry-only"){if(argc!=2)usage();return {RunMode::GeometryOnly,{ },false};}
	if(first=="--geometry-sweep"){if(argc!=2)usage();return {RunMode::GeometrySweep,{ },false};}
	if(first=="--time-lattice"){if(argc!=2)usage();return {RunMode::TimeLattice,{ },false};}
	if(first=="--one-step"){if(argc==2)return {RunMode::OneStep,{ },false};if(argc==3)return {RunMode::OneStep,fs::path(argv[2]),true};usage();}
	if(first.rfind("--",0)==0||argc!=2)usage();
	return {RunMode::FullCycle,fs::path(argv[1]),true};
}
bool Active(iga::CellClassification value)
{return value==iga::CellClassification::Inside||value==iga::CellClassification::Cut;}
std::pair<std::size_t,std::size_t> ActiveNodeChanges(const std::vector<std::int32_t>& previous,const std::vector<std::int32_t>& current)
{
	std::size_t removals=0,additions=0,old_index=0,new_index=0;
	while(old_index<previous.size()||new_index<current.size()){
		if(new_index==current.size()||(old_index<previous.size()&&previous[old_index]<current[new_index])){++removals;++old_index;}
		else if(old_index==previous.size()||current[new_index]<previous[old_index]){++additions;++new_index;}
		else {++old_index;++new_index;}
	}
	return {removals,additions};
}
void CheckGeometryProbeEpoch(std::uint64_t index,const iga::MovingCutGeometry& geometry,long long active_node_delta,std::size_t active_node_removals,std::size_t active_node_additions,double build_s)
{
	const auto& d=geometry.Diagnostics(); const auto& volume=d.volume;
	Require(geometry.Surface().Usable()&&geometry.Ghost().Usable(),"LV geometry probe catalog is unusable");
	Require(geometry.Domain().Diagnostics().ambiguous_count==0&&volume.ambiguous_samples==0&&volume.predicate_ambiguities==0,"LV geometry probe has an ambiguity");
	Require(d.catalog_lower_physical_volume_m3<=d.closed_surface_physical_volume_m3&&d.closed_surface_physical_volume_m3<=d.catalog_upper_physical_volume_m3,"LV geometry probe closed volume is outside its bracket");
	for(const auto& cell:geometry.Volume().Cells()){
		Require(cell.usable,"LV geometry probe volume cell is unusable");
		Require(std::isfinite(cell.diagnostics.estimated_reference_volume)&&cell.diagnostics.estimated_reference_volume>=0.&&cell.diagnostics.estimated_reference_volume<=1.,"LV geometry probe estimate is outside [0,1]");
		if(cell.classification==iga::CellClassification::Cut&&cell.diagnostics.logical_output_points==0)
			Require(cell.diagnostics.unresolved_reference_volume==0.,"LV geometry probe has an unresolved empty Cut cell");
	}
	if(index==2){
		for(const std::uint64_t id:{146,156})Require(geometry.Volume().Cell(id).diagnostics.rescue_effective_depth==9,
			"LV geometry probe rescue depth changed");
	}
	std::size_t removals=0,additions=0;
	for(const auto& transition:d.transition_counts){if(Active(transition.first.first)&&!Active(transition.first.second))removals+=transition.second;if(!Active(transition.first.first)&&Active(transition.first.second))additions+=transition.second;}
	std::cout.precision(17);
	std::cout<<"phase7_geometry_epoch index="<<index<<" time="<<d.time_s<<" build_s="<<build_s<<" logical_points="<<volume.logical_output_points<<" active_nodes="<<d.active_nodes<<" active_node_delta="<<active_node_delta<<" active_node_removals="<<active_node_removals<<" active_node_additions="<<active_node_additions<<" active_cells="<<d.active_cells<<" transitions="<<d.transitions.size()<<" removals="<<removals<<" additions="<<additions<<" classification_ambiguities="<<geometry.Domain().Diagnostics().ambiguous_count<<" volume_ambiguities="<<volume.ambiguous_samples<<" predicate_ambiguities="<<volume.predicate_ambiguities<<" certified_reference="<<volume.certified_reference_volume<<" unresolved_reference="<<volume.unresolved_reference_volume<<" estimated_reference="<<volume.estimated_reference_volume<<" usable=true bracket=["<<d.catalog_lower_physical_volume_m3<<','<<d.catalog_upper_physical_volume_m3<<"] closed_volume="<<d.closed_surface_physical_volume_m3<<std::endl;
}
void GeometrySweep(const iga::PrescribedSurfaceMotion& motion)
{
	constexpr std::uint64_t steps=kCycleSteps; auto options=Options();
	std::unique_ptr<iga::MovingCutGeometry> previous; std::vector<std::int32_t> previous_nodes; double ed_volume=0.,es_volume=0.; std::size_t previous_active_nodes=0,worst_points=0; bool contraction_removals=false,expansion_additions=false;
	const auto started=std::chrono::steady_clock::now();
	for(std::uint64_t index=0;index<=steps;++index){
		auto evaluation=LatticeEvaluation(motion,index);
		const auto build_started=std::chrono::steady_clock::now(); auto geometry=iga::MovingCutGeometry::Build(options.grid,std::move(evaluation),options.geometry,previous.get()); const double build_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-build_started).count();
		const auto& d=geometry->Diagnostics(); const long long active_node_delta=index ? static_cast<long long>(d.active_nodes)-static_cast<long long>(previous_active_nodes) : 0; const auto node_changes=ActiveNodeChanges(previous_nodes,geometry->Ghost().ActiveNodes()); CheckGeometryProbeEpoch(index,*geometry,active_node_delta,node_changes.first,node_changes.second,build_s); worst_points=std::max(worst_points,d.volume.logical_output_points);
		for(const auto& transition:d.transition_counts){if(Active(transition.first.first)&&!Active(transition.first.second)&&index<=8)contraction_removals=true;if(!Active(transition.first.first)&&Active(transition.first.second)&&index>8)expansion_additions=true;}
		if(index==0)ed_volume=d.closed_surface_physical_volume_m3;
		if(index==8)es_volume=d.closed_surface_physical_volume_m3;
		if(index==8){const auto es=LatticeEvaluation(motion,8);Require(geometry->Evaluation().EvaluatedTimeS()==CanonicalTime(8)&&geometry->Evaluation().SourceVerticesM()==motion.Frames()[1].surface.vertices&&geometry->Evaluation().Surface().CanonicalSha256()==es.Surface().CanonicalSha256(),"LV geometry probe ES landmark surface changed");}
		if(index==steps){const auto ed=LatticeEvaluation(motion,0),final=LatticeEvaluation(motion,steps);Require(geometry->Evaluation().EvaluatedTimeS()==iga::IdealizedLeftVentricleFixture::PeriodS&&geometry->Evaluation().SourceVerticesM()==motion.Frames().back().surface.vertices&&geometry->Evaluation().Surface().CanonicalSha256()==final.Surface().CanonicalSha256()&&final.SourceVerticesM()==ed.SourceVerticesM()&&final.Surface().CanonicalSha256()==ed.Surface().CanonicalSha256(),"LV geometry probe ED periodicity changed");Require(d.closed_surface_physical_volume_m3==ed_volume,"LV geometry probe ED closed volume is not periodic");}
		previous_active_nodes=d.active_nodes; previous_nodes=geometry->Ghost().ActiveNodes(); previous=std::move(geometry);
	}
	Require(contraction_removals&&expansion_additions,"LV geometry probe did not exercise contraction removals and expansion additions");
	Require(worst_points<=110000,"LV geometry probe logical-point reduction gate failed");
	{
		auto cap_options=options.geometry; cap_options.volume.empty_rule_rescue_max_depth=8;
		bool rejected=false;
		try { auto ignored=iga::MovingCutGeometry::Build(options.grid,LatticeEvaluation(motion,2),cap_options); (void)ignored; }
		catch(const std::exception& error) { rejected=std::string(error.what()).find("empty-rule rescue exhausted")!=std::string::npos; }
		Require(rejected,"LV geometry probe cap-eight rescue did not fail closed");
	}
	Require(std::abs(es_volume/ed_volume-iga::IdealizedLeftVentricleFixture::EndSystolicVolumeRatio)<=5e-13*iga::IdealizedLeftVentricleFixture::EndSystolicVolumeRatio,"LV geometry probe ES/ED closed-volume ratio changed");
	std::cout<<"phase7_geometry_summary epochs=17 depth=2 worst_logical_points="<<worst_points<<" ed_volume="<<ed_volume<<" es_volume="<<es_volume<<" elapsed_s="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<std::endl;
}
void TimeLatticeProbe(const iga::PrescribedSurfaceMotion& motion)
{
	const double nominal_dt=iga::IdealizedLeftVentricleFixture::PeriodS/static_cast<double>(kCycleSteps);
	double minimum_dt=std::numeric_limits<double>::infinity(),maximum_dt=0.;
	Require(CanonicalTime(0)==0.&&CanonicalTime(8)==iga::IdealizedLeftVentricleFixture::PeriodS/2.&&CanonicalTime(kCycleSteps)==iga::IdealizedLeftVentricleFixture::PeriodS,"LV canonical time lattice does not land exactly on ED/ES/ED");
	for(std::uint64_t index=1;index<=kCycleSteps;++index){const auto step=NextTimeLatticeStep(index,CanonicalTime(index-1));minimum_dt=std::min(minimum_dt,step.dt);maximum_dt=std::max(maximum_dt,step.dt);}
	const auto ed=LatticeEvaluation(motion,0),es=LatticeEvaluation(motion,8),final=LatticeEvaluation(motion,kCycleSteps);
	Require(es.EvaluatedTimeS()==iga::IdealizedLeftVentricleFixture::PeriodS/2.&&es.SourceVerticesM()==motion.Frames()[1].surface.vertices,"LV time-lattice ES frame is not exact");
	Require(final.EvaluatedTimeS()==iga::IdealizedLeftVentricleFixture::PeriodS&&final.SourceVerticesM()==ed.SourceVerticesM()&&final.Surface().CanonicalSha256()==ed.Surface().CanonicalSha256(),"LV time-lattice final ED frame is not exact");
	std::cout.precision(17);
	std::cout<<"phase7_time_lattice_summary steps="<<kCycleSteps<<" nominal_dt="<<nominal_dt<<" min_dt="<<minimum_dt<<" max_dt="<<maximum_dt<<" es_time="<<es.EvaluatedTimeS()<<" final_time="<<final.EvaluatedTimeS()<<std::endl;
}
void CheckSnapshot(const iga::MovingImmersedFlowSnapshot& s,const iga::MovingCutGeometry& g,bool require_port_flows_available)
{
	const auto& m=s.Metrics();const auto& d=g.Diagnostics();const auto& request=s.Request();
	Require(std::abs(m.quadrature_volume_m3-m.audited_volume_m3)<=1e-12*m.audited_volume_m3,"snapshot/catalog volume audit failed");
	const double tolerance=1e-12*d.closed_surface_physical_volume_m3;Require(d.catalog_lower_physical_volume_m3<=d.closed_surface_physical_volume_m3+tolerance&&d.catalog_upper_physical_volume_m3+tolerance>=d.closed_surface_physical_volume_m3,"catalog bracket excludes closed volume");
	for(const auto& p:s.Points())for(const double x:{p.pressure,p.speed_m_per_s,p.q_criterion_per_s2,p.enstrophy_density_per_s2})Require(std::isfinite(x),"snapshot field is nonfinite");
	for(const double x:{m.q_positive_volume_fraction,m.stagnant_volume_fraction,m.well_mixed_replacement_fraction_over_step})Require(std::isfinite(x)&&x>=0.&&x<=1.,"snapshot fraction is invalid");
	Require(request.port_flows_available==require_port_flows_available,"snapshot port-flow availability differs from committed-capture semantics");
	if(!request.port_flows_available){
		Require(!m.endpoint_turnover_available&&!m.endpoint_turnover_time_s&&m.endpoint_turnover_rate_per_s==0.&&m.well_mixed_replacement_fraction_over_step==0.,"unmeasured ports fabricated endpoint turnover metrics");
		return;
	}
	Require(m.endpoint_turnover_available,"available retained port flows did not retain turnover availability");
	double qin=0.;for(const auto& port:request.port_flows)qin+=std::max(0.,-port.outward_flow_m3_s);
	Require(std::abs(m.inlet_flow_m3_s-qin)<=1e-12*std::max(1.,qin),"snapshot inlet flow differs from retained port measurements");
	if(qin>0.){
		Require(m.endpoint_turnover_time_s.has_value(),"positive retained inlet flow omitted turnover time");
		Require(std::abs(m.endpoint_turnover_rate_per_s*(*m.endpoint_turnover_time_s)-1.)<=1e-12,"turnover rate/time identity failed");
		const double replacement=-std::expm1(-qin*request.dt_s/m.quadrature_volume_m3);Require(std::abs(m.well_mixed_replacement_fraction_over_step-replacement)<=1e-12,"well-mixed replacement formula changed");
	} else Require(!m.endpoint_turnover_time_s&&m.endpoint_turnover_rate_per_s==0.&&m.well_mixed_replacement_fraction_over_step==0.,"known-zero retained inlet flow has noncanonical turnover metrics");
}
void Print(std::uint64_t index,const iga::MovingImmersedTransientFlowRuntime& r,const iga::MovingImmersedFlowSnapshot& s,const iga::MovingImmersedFlowSnapshotPublication& p)
{
	const auto& d=r.CommittedDiagnostics();const auto& m=s.Metrics();const auto& request=s.Request();std::cout.precision(17);
	std::cout<<"phase7_epoch index="<<index<<" time="<<r.CommittedGlobalState().TimeS()<<" nodes="<<d.active_nodes<<" newton="<<d.nonlinear_iterations<<" ksp="<<d.ksp_iterations<<" residual="<<d.residual_norm<<" volume="<<m.quadrature_volume_m3<<" flow_in="<<m.inlet_flow_m3_s<<" eta_mu=["<<d.wall_penalty.minimum_eta_mu<<','<<d.wall_penalty.maximum_eta_mu<<"] eta_t=["<<d.wall_penalty.minimum_eta_t<<','<<d.wall_penalty.maximum_eta_t<<"] eta=["<<d.wall_penalty.minimum_eta<<','<<d.wall_penalty.maximum_eta<<"] eta_t_fraction=["<<d.wall_penalty.minimum_eta_t_fraction<<','<<d.wall_penalty.maximum_eta_t_fraction<<']';
	if(index==0){Require(!request.transition_available&&request.dt_s==0.,"LV epoch zero unexpectedly has a committed transition");std::cout<<" transition_available=false conservation=unavailable";}
	else {Require(request.transition_available&&request.dt_s>0.,"LV committed epoch is missing its transition");const auto c=r.ConservationDiagnostics();std::cout<<" transition_available=true discrete_continuity="<<c.normalized_discrete_moving_wall_continuity_defect<<" reynolds="<<c.normalized_reynolds_defect<<" moving_mass="<<c.normalized_moving_mass_defect<<" divergence_quadrature_defect="<<c.normalized_divergence_theorem_defect<<" wall_leak="<<c.normalized_wall_relative_leakage;}
	std::cout<<" enstrophy="<<m.enstrophy_integral_m3_per_s2<<" q_volume="<<m.q_positive_volume_m3<<" snapshot_hash="<<s.ContentHashSha256()<<" file_hash="<<p.metrics_sha256<<std::endl;
}
struct WallGateReport { WallQuadratureMetrics quadrature; double ratio=0.; bool public_pass=false,release_pass=false; };
WallGateReport PrintWallGate(std::uint64_t index,const iga::MovingImmersedTransientFlowRuntime& runtime,const iga::MovingImmersedFlowSnapshotOptions& snapshot_options)
{
	const auto quadrature=IntegrateWallQuadrature(runtime);const auto& d=runtime.CommittedDiagnostics();const double scale=std::max(snapshot_options.stagnant_speed_threshold_m_per_s,quadrature.material_rms_m_per_s),ratio=quadrature.relative_rms_m_per_s/scale;const bool public_pass=ratio<=kWallRelativeTolerance,release_pass=ratio<=kWallRelativeReleaseTolerance;std::cout.precision(17);
	std::cout<<"phase7_wall_gate index="<<index<<" time="<<runtime.CommittedGlobalState().TimeS()<<" material_wall_rms="<<quadrature.material_rms_m_per_s<<" relative_full_rms="<<quadrature.relative_rms_m_per_s<<" relative_normal_rms="<<quadrature.normal_relative_rms_m_per_s<<" relative_tangential_rms="<<quadrature.tangential_relative_rms_m_per_s<<" wall_area="<<quadrature.area_m2<<" scale="<<scale<<" ratio="<<ratio<<" public_pass="<<(public_pass?"true":"false")<<" release_pass="<<(release_pass?"true":"false")<<" eta_mu=["<<d.wall_penalty.minimum_eta_mu<<','<<d.wall_penalty.maximum_eta_mu<<"] eta_t=["<<d.wall_penalty.minimum_eta_t<<','<<d.wall_penalty.maximum_eta_t<<"] eta=["<<d.wall_penalty.minimum_eta<<','<<d.wall_penalty.maximum_eta<<"] eta_mu_h_n_over_mu=["<<d.wall_penalty.minimum_eta_mu_h_n_over_mu<<','<<d.wall_penalty.maximum_eta_mu_h_n_over_mu<<"] eta_t_dt_over_rho_h_n=["<<d.wall_penalty.minimum_eta_t_dt_over_rho_h_n<<','<<d.wall_penalty.maximum_eta_t_dt_over_rho_h_n<<"] eta_mu_fraction=["<<d.wall_penalty.minimum_eta_mu_fraction<<','<<d.wall_penalty.maximum_eta_mu_fraction<<"] eta_t_fraction=["<<d.wall_penalty.minimum_eta_t_fraction<<','<<d.wall_penalty.maximum_eta_t_fraction<<']'<<std::endl;
	return {quadrature,ratio,public_pass,release_pass};
}
bool SameBits(double left,double right)
{return std::memcmp(&left,&right,sizeof(double))==0;}
bool SamePorts(const std::vector<iga::ImmersedTransientFlowDiagnostics::Port>& left,const std::vector<iga::ImmersedTransientFlowDiagnostics::Port>& right)
{
	if(left.size()!=right.size())return false;
	for(std::size_t i=0;i<left.size();++i){const auto& a=left[i];const auto& b=right[i];if(a.id!=b.id||a.boundary_label!=b.boundary_label||a.control_mode!=b.control_mode||!SameBits(a.target,b.target)||!SameBits(a.multiplier,b.multiplier)||!SameBits(a.controller_error,b.controller_error)||a.multiplier_row!=b.multiplier_row||!SameBits(a.measurement.area_m2,b.measurement.area_m2)||!SameBits(a.measurement.outward_flow_m3_s,b.measurement.outward_flow_m3_s)||!SameBits(a.measurement.mean_pressure_pa,b.measurement.mean_pressure_pa)||!SameBits(a.measurement.mean_normal_traction_pa,b.measurement.mean_normal_traction_pa)||!SameBits(a.measurement.mean_velocity_squared_m2_s2,b.measurement.mean_velocity_squared_m2_s2)||a.measurement_valid!=b.measurement_valid)return false;}
	return true;
}
struct FlowRateObservation { std::string id; int boundary_label; double target=0.,outward_flow=0.,controller_error=0.; };
double ControllerGate(const iga::MovingImmersedTransientFlowOptions& options,double target)
{return std::max(options.flow.flow_controller_absolute_tolerance_m3_s,std::abs(target)*kControllerAcceptanceRelativeTolerance);}
bool SameFlowRateObservations(const std::vector<FlowRateObservation>& left,const std::vector<FlowRateObservation>& right)
{
	if(left.size()!=right.size())return false;
	for(std::size_t i=0;i<left.size();++i)if(left[i].id!=right[i].id||left[i].boundary_label!=right[i].boundary_label||!SameBits(left[i].target,right[i].target)||!SameBits(left[i].outward_flow,right[i].outward_flow)||!SameBits(left[i].controller_error,right[i].controller_error))return false;
	return true;
}
std::vector<FlowRateObservation> CheckTrialFlowRatePorts(const iga::MovingImmersedTransientFlowRuntime& runtime,const iga::MovingImmersedTransientFlowOptions& options,double& max_controller)
{
	const auto conservation=runtime.ConservationDiagnostics(); std::vector<FlowRateObservation> result;
	for(const auto& configured:runtime.ConfiguredPorts())if(configured.control_mode==iga::ImmersedFlowPortControlMode::FlowRate){
		const auto found=conservation.fluid_surface_outward_flow_by_boundary_label_m3_s.find(configured.boundary_label);
		Require(found!=conservation.fluid_surface_outward_flow_by_boundary_label_m3_s.end()&&std::isfinite(configured.value)&&std::isfinite(found->second),"LV trial flow-rate conservation measurement is unavailable");
		// Surface flux is outward by contract, so a negative inlet target is target-Q_outward.
		const double error=configured.value-found->second; Require(std::abs(error)<=ControllerGate(options,configured.value),"LV flow controller gate failed");
		result.push_back({configured.id,configured.boundary_label,configured.value,found->second,error}); max_controller=std::max(max_controller,std::abs(error));
	}
	return result;
}
void CheckCommittedFlowRatePorts(const iga::MovingImmersedTransientFlowRuntime& runtime,const iga::MovingImmersedTransientFlowOptions& options,const std::vector<FlowRateObservation>& trial_observations)
{
	constexpr double summation_order_factor=128.; const auto& ports=runtime.CommittedDiagnostics().ports;
	Require(trial_observations.size()==static_cast<std::size_t>(std::count_if(runtime.ConfiguredPorts().begin(),runtime.ConfiguredPorts().end(),[](const auto& port){return port.control_mode==iga::ImmersedFlowPortControlMode::FlowRate;})),"LV committed flow-rate observation count is incomplete");
	for(const auto& trial:trial_observations){
		const auto found=std::find_if(ports.begin(),ports.end(),[&](const auto& port){return port.id==trial.id;});
		Require(found!=ports.end()&&found->measurement_valid&&found->boundary_label==trial.boundary_label&&found->control_mode==iga::ImmersedFlowPortControlMode::FlowRate&&SameBits(found->target,trial.target)&&std::isfinite(found->measurement.outward_flow_m3_s)&&std::isfinite(found->controller_error),"LV committed flow-rate port record is invalid");
		// The port and conservation reductions may traverse surface quadrature in different orders.
		const double flux_tolerance=summation_order_factor*std::numeric_limits<double>::epsilon()*std::max({std::abs(trial.target),std::abs(trial.outward_flow),std::abs(found->measurement.outward_flow_m3_s)});
		Require(std::abs(found->measurement.outward_flow_m3_s-trial.outward_flow)<=flux_tolerance,"LV committed port flow differs from precommit conservation measurement");
		Require(std::abs(found->controller_error)<=ControllerGate(options,found->target),"LV committed flow controller gate failed");
	}
}
std::vector<FlowRateObservation> CheckTrialNumerics(const iga::MovingImmersedTransientFlowRuntime& runtime,const iga::MovingImmersedTransientFlowOptions& options,double& max_linear,double& max_controller)
{
	const auto& d=runtime.TrialDiagnostics(); Require(d.converged,"LV nonlinear convergence flag is false");
	for(const auto& n:d.newton_steps){if(n.ksp_iterations)Require(n.ksp_reason>0,"LV KSP did not converge");Require(n.linear_relative_residual<=kLinearTolerance,"LV true linear residual gate failed");max_linear=std::max(max_linear,n.linear_relative_residual);}
	const auto ports=CheckTrialFlowRatePorts(runtime,options,max_controller);
	const auto c=runtime.ConservationDiagnostics(); Require(std::isfinite(c.normalized_divergence_theorem_defect)&&c.normalized_discrete_moving_wall_continuity_defect<=1e-8&&c.normalized_reynolds_defect<=kConservationTolerance&&c.normalized_moving_mass_defect<=kConservationTolerance&&c.normalized_wall_relative_leakage<=1e-2,"LV moving conservation gate failed");
	return ports;
}
struct DivergenceQuadratureAudit { double uncertainty_m3=0.,integral_m3_s=0.,rdiv_m3_s=0.; };
DivergenceQuadratureAudit AuditDivergenceQuadrature(const iga::MovingImmersedTransientFlowRuntime& runtime,const iga::MovingImmersedTransientFlowOptions& options,std::uint32_t depth)
{
	auto geometry_options=options.geometry;geometry_options.volume.max_depth=depth;geometry_options.volume.empty_rule_rescue_max_depth=9;
	auto geometry=iga::MovingCutGeometry::Build(options.grid,runtime.TrialGeometry().Evaluation(),geometry_options);
	const auto& volume=geometry->Volume();const auto& domain=geometry->Domain();const auto& layout=runtime.TrialLayout();const auto state=runtime.TrialState();
	Require(geometry->Diagnostics().volume.ambiguous_samples==0&&geometry->Diagnostics().volume.predicate_ambiguities==0&&domain.Diagnostics().ambiguous_count==0,"LV divergence audit has an ambiguity");
	DivergenceQuadratureAudit result;result.uncertainty_m3=geometry->Diagnostics().catalog_upper_physical_volume_m3-geometry->Diagnostics().catalog_lower_physical_volume_m3;
	Require(std::isfinite(result.uncertainty_m3)&&result.uncertainty_m3>=0.,"LV divergence audit uncertainty is invalid");
	for(std::uint64_t cell=0;cell<domain.Cells().size();++cell){
		Require(volume.Cell(cell).usable,"LV divergence audit volume cell is unusable");
		if(domain.Cells()[cell].classification==iga::CellClassification::Outside)continue;
		const auto element=domain.Background().MaterializeElement(cell);
		auto add=[&](const iga::VolumeQuadraturePoint& point){const auto basis=iga::EvaluateBasis(element,point.parametric[0],point.parametric[1],point.parametric[2],false);double divergence=0.;for(std::size_t a=0;a<element.connectivity.size();++a){const auto local=layout.LocalNode(element.connectivity[a]);for(int q=0;q<3;++q)divergence+=PetscRealPart(state[4*local+q])*basis.gradient[a][q];}result.integral_m3_s+=point.weight*basis.raw_determinant*divergence;};
		if(volume.StorageMode()==iga::CutCellVolumeQuadratureStorageMode::Expanded)for(const auto& point:volume.UsableRule(domain,cell).Points())add(point);else iga::ForEachVolumePoint(volume.UsableCompactRule(domain,cell),add);
	}
	Require(std::isfinite(result.integral_m3_s),"LV divergence audit integral is nonfinite");
	result.rdiv_m3_s=result.integral_m3_s-runtime.ConservationDiagnostics().total_fluid_surface_outward_flow_m3_s;
	Require(std::isfinite(result.rdiv_m3_s),"LV divergence audit residual is nonfinite");return result;
}
void CheckDivergenceQuadratureRefinementAudit(const iga::MovingImmersedTransientFlowRuntime& runtime,const iga::MovingImmersedTransientFlowOptions& options)
{
	const std::array<DivergenceQuadratureAudit,3> audit{{AuditDivergenceQuadrature(runtime,options,1),AuditDivergenceQuadrature(runtime,options,2),AuditDivergenceQuadrature(runtime,options,3)}};const auto production=runtime.ConservationDiagnostics();
	Require(audit[2].uncertainty_m3<audit[1].uncertainty_m3&&audit[1].uncertainty_m3<audit[0].uncertainty_m3,"LV divergence audit uncertainty did not strictly refine");
	Require(std::abs(audit[2].rdiv_m3_s)<std::abs(audit[1].rdiv_m3_s)&&std::abs(audit[1].rdiv_m3_s)<std::abs(audit[0].rdiv_m3_s)&&std::abs(audit[2].rdiv_m3_s)<=.25*std::abs(audit[0].rdiv_m3_s),"LV divergence audit residual did not strictly refine");
	Require(std::abs(audit[2].integral_m3_s-audit[1].integral_m3_s)<std::abs(audit[1].integral_m3_s-audit[0].integral_m3_s),"LV divergence audit endpoint integral did not refine");
	const double summation_tolerance=512.*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(production.endpoint_volume_divergence_m3_s),std::abs(production.divergence_theorem_defect_m3_s)})*std::max<std::size_t>(1,runtime.TrialGeometry().Diagnostics().volume.logical_output_points);
	Require(std::abs(audit[1].integral_m3_s-production.endpoint_volume_divergence_m3_s)<=summation_tolerance&&std::abs(audit[1].rdiv_m3_s-production.divergence_theorem_defect_m3_s)<=summation_tolerance,"LV depth-two divergence audit does not reproduce production endpoint quadrature");
	std::cout<<"phase7_divergence_quadrature_audit U1="<<audit[0].uncertainty_m3<<" U2="<<audit[1].uncertainty_m3<<" U3="<<audit[2].uncertainty_m3<<" I1="<<audit[0].integral_m3_s<<" I2="<<audit[1].integral_m3_s<<" I3="<<audit[2].integral_m3_s<<" R1="<<audit[0].rdiv_m3_s<<" R2="<<audit[1].rdiv_m3_s<<" R3="<<audit[2].rdiv_m3_s<<std::endl;
}
void OneStepProbe(const iga::PrescribedSurfaceMotion& motion,const CommandLine& command)
{
	const auto started=std::chrono::steady_clock::now(); const bool temporary=!command.retain_output; const fs::path root=temporary ? fs::temp_directory_path()/"tubularflowiga-phase7-lv-one-step" : command.output_root;
	if(temporary){std::error_code ignored;fs::remove_all(root,ignored);} fs::create_directories(root); const auto stem=root/"lv_one_step";
	auto options=Options(); const auto initial_started=std::chrono::steady_clock::now(); iga::MovingImmersedTransientFlowRuntime runtime(LatticeEvaluation(motion,0),options); runtime.InitializeCommittedGlobalState(InitialRotationalState(runtime,options.grid)); const double initial_runtime_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-initial_started).count();
	iga::MovingImmersedFlowSnapshotOptions snapshot_options;snapshot_options.stagnant_speed_threshold_m_per_s=1e-4;snapshot_options.maximum_points=2000000;snapshot_options.maximum_output_bytes=500000000;
	bool q_positive=false; const auto initial_snapshot_started=std::chrono::steady_clock::now();
	{
		auto initial=iga::BuildCommittedMovingImmersedFlowSnapshot(runtime,snapshot_options);CheckSnapshot(initial,runtime.CommittedGeometry(),false);auto publication=iga::MovingImmersedFlowSnapshotPublisher::Publish(initial,stem,Identity(runtime),std::nullopt);Print(0,runtime,initial,publication);q_positive=initial.Metrics().enstrophy_integral_m3_per_s2>0.&&initial.Metrics().q_positive_volume_m3>1024.*std::numeric_limits<double>::epsilon()*initial.Metrics().quadrature_volume_m3;
	}
	const double initial_snapshot_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-initial_snapshot_started).count();
	const auto lattice=NextTimeLatticeStep(1,runtime.CommittedGlobalState().TimeS()); Require(lattice.dt==.05,"LV one-step canonical dt is no longer exactly .05"); const auto begin_started=std::chrono::steady_clock::now(); runtime.BeginTrial(motion.Evaluate(lattice.target,lattice.source,lattice.target),1,lattice.dt); const double begin_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin_started).count(); const auto committed_ports=runtime.CommittedDiagnostics().ports;
	double max_linear=0.,max_controller=0.; const auto solve_started=std::chrono::steady_clock::now(); Require(runtime.SolveTrial(),"LV one-step solve did not converge"); const double solve_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-solve_started).count(); Require(SamePorts(committed_ports,runtime.CommittedDiagnostics().ports),"LV solved trial published committed port diagnostics"); const auto solved=runtime.TrialDiagnostics(); const auto solved_observations=CheckTrialNumerics(runtime,options,max_linear,max_controller); CheckDivergenceQuadratureRefinementAudit(runtime,options); const auto before=runtime.TrialState(); const auto map=runtime.Diagnostics().map_identity_sha256,input=solved.input_hash_sha256,attempt=solved.attempt_hash_sha256;
	const auto replay_started=std::chrono::steady_clock::now(); runtime.Rollback(); Require(runtime.SolveTrial(),"LV one-step rollback retry did not converge"); const double replay_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-replay_started).count(); Require(SamePorts(committed_ports,runtime.CommittedDiagnostics().ports),"LV replayed trial published committed port diagnostics"); const auto replay=runtime.TrialDiagnostics(); const auto replay_observations=CheckTrialNumerics(runtime,options,max_linear,max_controller); Require(map==runtime.Diagnostics().map_identity_sha256&&input==replay.input_hash_sha256&&attempt==replay.attempt_hash_sha256&&solved.solved_state_hash_sha256==replay.solved_state_hash_sha256&&before==runtime.TrialState()&&SameFlowRateObservations(solved_observations,replay_observations),"LV one-step rollback retry changed deterministic identities, state, or independent flow observation");
	const auto conservation=runtime.ConservationDiagnostics(); const auto commit_started=std::chrono::steady_clock::now(); runtime.PrepareCommit();runtime.FinalizeCommit(); const double commit_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-commit_started).count(); Require(runtime.Diagnostics().rollback_count==1&&runtime.Diagnostics().prepare_count==1&&runtime.Diagnostics().finalize_count==1,"LV one-step transactional lifecycle evidence changed"); CheckCommittedFlowRatePorts(runtime,options,replay_observations);
	double wall_ratio=0.; const auto snapshot_started=std::chrono::steady_clock::now();
	{
		auto snapshot=iga::BuildCommittedMovingImmersedFlowSnapshot(runtime,snapshot_options);const auto wall_gate=PrintWallGate(1,runtime,snapshot_options);CheckSnapshot(snapshot,runtime.CommittedGeometry(),true);wall_ratio=wall_gate.ratio;const double wall_metric_tolerance=512.*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(snapshot.Metrics().wall_relative_velocity_rms_m_per_s),std::abs(wall_gate.quadrature.relative_rms_m_per_s)});Require(std::abs(snapshot.Metrics().wall_relative_velocity_rms_m_per_s-wall_gate.quadrature.relative_rms_m_per_s)<=wall_metric_tolerance,"LV wall quadrature full RMS differs from snapshot metric");Require(wall_gate.public_pass,"LV wall-relative velocity gate failed");Require(wall_gate.release_pass,"LV wall-relative velocity release margin failed");q_positive=q_positive||(snapshot.Metrics().enstrophy_integral_m3_per_s2>0.&&snapshot.Metrics().q_positive_volume_m3>1024.*std::numeric_limits<double>::epsilon()*snapshot.Metrics().quadrature_volume_m3);const auto state_hash=runtime.CommittedGlobalState().HashSha256();iga::MovingImmersedFlowSnapshotPublisher::SetTestFault(iga::MovingImmersedFlowSnapshotPublicationFault::AfterTemporaryVtuWrite);bool failed=false;try{(void)iga::MovingImmersedFlowSnapshotPublisher::Publish(snapshot,stem,Identity(runtime),std::nullopt);}catch(const std::exception&){failed=true;}iga::MovingImmersedFlowSnapshotPublisher::SetTestFault(iga::MovingImmersedFlowSnapshotPublicationFault::None);Require(failed&&!fs::exists(iga::MovingImmersedFlowSnapshotPublisher::EpochDirectory(stem,1))&&runtime.CommittedGlobalState().HashSha256()==state_hash,"LV one-step injected publisher failure left an epoch or changed solver state");const auto publication=iga::MovingImmersedFlowSnapshotPublisher::Publish(snapshot,stem,Identity(runtime),std::nullopt);const auto retry=iga::MovingImmersedFlowSnapshotPublisher::Publish(snapshot,stem,Identity(runtime),std::nullopt);Require(publication.metrics_sha256==retry.metrics_sha256&&publication.vtu_sha256==retry.vtu_sha256,"LV one-step publisher retry is not idempotent");Print(1,runtime,snapshot,publication);
	}
	const double snapshot_publish_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-snapshot_started).count(); Require(q_positive,"LV one-step did not retain a physical Q-positive rotating state"); iga::MovingImmersedFlowSnapshotPublisher::RebuildCollection(stem); std::ifstream pvd(iga::MovingImmersedFlowSnapshotPublisher::CollectionPath(stem)); const std::string pvd_text((std::istreambuf_iterator<char>(pvd)),{}); std::size_t epochs=0,pos=0; while((pos=pvd_text.find("<DataSet timestep=",pos))!=std::string::npos){++epochs;++pos;} Require(epochs==2&&fs::is_regular_file(iga::MovingImmersedFlowSnapshotPublisher::EpochDirectory(stem,0)/"fields.vtu")&&fs::is_regular_file(iga::MovingImmersedFlowSnapshotPublisher::EpochDirectory(stem,1)/"fields.vtu"),"LV one-step PVD evidence is incomplete");
	std::cout<<"phase7_one_step_summary grid=6x6x7 depth=2 logical_points="<<runtime.CommittedGeometry().Diagnostics().volume.logical_output_points<<" max_linear="<<max_linear<<" max_controller="<<max_controller<<" max_discrete_continuity="<<conservation.normalized_discrete_moving_wall_continuity_defect<<" max_divergence_quadrature_defect="<<conservation.normalized_divergence_theorem_defect<<" max_reynolds="<<conservation.normalized_reynolds_defect<<" max_moving_mass="<<conservation.normalized_moving_mass_defect<<" max_wall_leakage="<<conservation.normalized_wall_relative_leakage<<" wall_ratio="<<wall_ratio<<" initial_runtime_s="<<initial_runtime_s<<" initial_snapshot_publish_s="<<initial_snapshot_s<<" begin_s="<<begin_s<<" solve_s="<<solve_s<<" rollback_replay_s="<<replay_s<<" commit_s="<<commit_s<<" snapshot_publish_s="<<snapshot_publish_s<<" elapsed_s="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<std::endl;
	if(temporary){std::error_code ignored;fs::remove_all(root,ignored);}
}
}

int main(int argc,char** argv)
{
	int status=0;bool petsc_initialized=false;
	try {
		const auto command=ParseCommandLine(argc,argv); const auto started=std::chrono::steady_clock::now(); const auto motion=iga::IdealizedLeftVentricleFixture::Motion(); GeometryOnly(motion);
		if(command.mode==RunMode::GeometryOnly){std::cout<<"phase7_geometry ed_es_ratio="<<iga::IdealizedLeftVentricleFixture::EndSystolicVolumeRatio<<" sectors=12 rings=4 labels=84,6,6\n";return 0;}
		if(command.mode==RunMode::GeometrySweep){GeometrySweep(motion);return 0;}
		if(command.mode==RunMode::TimeLattice){TimeLatticeProbe(motion);return 0;}
		int petsc_argc=1; char** petsc_argv=argv;
		if(PetscInitialize(&petsc_argc,&petsc_argv,nullptr,nullptr)!=0){std::cerr<<"phase7_lv_closure: PetscInitialize failed\n";return 1;}petsc_initialized=true;
		if(command.mode==RunMode::OneStep){OneStepProbe(motion,command);if(PetscFinalize()!=0){std::cerr<<"phase7_lv_closure: PetscFinalize failed\n";return 1;}petsc_initialized=false;return 0;}
		{
		const bool temporary=!command.retain_output; const fs::path root=temporary ? fs::temp_directory_path()/"tubularflowiga-phase7-lv-closure" : command.output_root; if(temporary){std::error_code ignored;fs::remove_all(root,ignored);} fs::create_directories(root); const auto stem=root/"lv_cycle";
		auto options=Options(); iga::MovingImmersedTransientFlowRuntime runtime(LatticeEvaluation(motion,0),options);
		runtime.InitializeCommittedGlobalState(InitialRotationalState(runtime,options.grid));
		iga::MovingImmersedFlowSnapshotOptions snapshot_options;snapshot_options.stagnant_speed_threshold_m_per_s=1e-4;snapshot_options.maximum_points=2000000;snapshot_options.maximum_output_bytes=500000000;
		bool q_positive=false;
		{
			auto initial=iga::BuildCommittedMovingImmersedFlowSnapshot(runtime,snapshot_options);CheckSnapshot(initial,runtime.CommittedGeometry(),false);auto initial_pub=iga::MovingImmersedFlowSnapshotPublisher::Publish(initial,stem,Identity(runtime),std::nullopt);Print(0,runtime,initial,initial_pub);q_positive=initial.Metrics().enstrophy_integral_m3_per_s2>0.&&initial.Metrics().q_positive_volume_m3>1024.*std::numeric_limits<double>::epsilon()*initial.Metrics().quadrature_volume_m3;
		}
		double max_linear=0.,max_controller=0.,max_wall_ratio=0.,max_continuity=0.,max_divergence_quadrature=0.,max_reynolds=0.,max_mass=0.,max_leak=0.;bool contraction_transition=false,expansion_transition=false;
		constexpr std::uint64_t steps=kCycleSteps;const double nominal_dt=iga::IdealizedLeftVentricleFixture::PeriodS/steps;double minimum_dt=std::numeric_limits<double>::infinity(),maximum_dt=0.;
		for(std::uint64_t step=1;step<=steps;++step){const auto lattice=NextTimeLatticeStep(step,runtime.CommittedGlobalState().TimeS());minimum_dt=std::min(minimum_dt,lattice.dt);maximum_dt=std::max(maximum_dt,lattice.dt);runtime.BeginTrial(motion.Evaluate(lattice.target,lattice.source,lattice.target),step,lattice.dt);const auto committed_ports=runtime.CommittedDiagnostics().ports;Require(runtime.SolveTrial(),"LV step did not converge");Require(SamePorts(committed_ports,runtime.CommittedDiagnostics().ports),"LV solved trial published committed port diagnostics");auto solved=runtime.TrialDiagnostics();auto trial_observations=CheckTrialNumerics(runtime,options,max_linear,max_controller);
			if(step==1){const auto before=runtime.TrialState();const auto map=runtime.Diagnostics().map_identity_sha256,input=solved.input_hash_sha256,attempt=solved.attempt_hash_sha256;const auto solved_observations=trial_observations;runtime.Rollback();Require(runtime.SolveTrial(),"LV rollback retry did not converge");Require(SamePorts(committed_ports,runtime.CommittedDiagnostics().ports),"LV replayed trial published committed port diagnostics");const auto replay=runtime.TrialDiagnostics();trial_observations=CheckTrialNumerics(runtime,options,max_linear,max_controller);Require(map==runtime.Diagnostics().map_identity_sha256&&input==replay.input_hash_sha256&&attempt==replay.attempt_hash_sha256&&solved.solved_state_hash_sha256==replay.solved_state_hash_sha256&&before==runtime.TrialState()&&SameFlowRateObservations(solved_observations,trial_observations),"LV numerical rollback retry changed deterministic identities, state, or independent flow observation");}
			const auto c=runtime.ConservationDiagnostics();Require(c.source_time_s==lattice.source&&c.target_time_s==lattice.target&&c.dt_s==lattice.dt,"LV conservation record does not retain the canonical lattice transition");max_continuity=std::max(max_continuity,c.normalized_discrete_moving_wall_continuity_defect);max_divergence_quadrature=std::max(max_divergence_quadrature,c.normalized_divergence_theorem_defect);max_reynolds=std::max(max_reynolds,c.normalized_reynolds_defect);max_mass=std::max(max_mass,c.normalized_moving_mass_defect);max_leak=std::max(max_leak,c.normalized_wall_relative_leakage);Require(std::isfinite(c.normalized_divergence_theorem_defect)&&c.normalized_discrete_moving_wall_continuity_defect<=1e-8&&c.normalized_reynolds_defect<=kConservationTolerance&&c.normalized_moving_mass_defect<=kConservationTolerance&&c.normalized_wall_relative_leakage<=1e-2,"LV moving conservation gate failed");
			const auto& trans=runtime.TrialGeometry().Diagnostics().transition_counts;for(const auto& item:trans)if(item.first.first!=item.first.second){if(step<=8)contraction_transition=true;else expansion_transition=true;}
			runtime.PrepareCommit();runtime.FinalizeCommit();Require(runtime.CommittedGlobalState().TimeS()==lattice.target,"LV committed state did not land on the canonical target time");if(step==8){const auto es=LatticeEvaluation(motion,8);const auto& actual=runtime.CommittedGeometry().Evaluation();Require(actual.EvaluatedTimeS()==iga::IdealizedLeftVentricleFixture::PeriodS/2.&&actual.SourceVerticesM()==motion.Frames()[1].surface.vertices&&actual.Surface().CanonicalSha256()==es.Surface().CanonicalSha256(),"LV committed ES landmark surface/hash changed");}if(step==steps){const auto ed=LatticeEvaluation(motion,0),final=LatticeEvaluation(motion,steps);const auto& actual=runtime.CommittedGeometry().Evaluation();Require(runtime.CommittedGlobalState().TimeS()==iga::IdealizedLeftVentricleFixture::PeriodS&&actual.EvaluatedTimeS()==iga::IdealizedLeftVentricleFixture::PeriodS&&actual.SourceVerticesM()==motion.Frames().back().surface.vertices&&actual.Surface().CanonicalSha256()==final.Surface().CanonicalSha256()&&final.SourceVerticesM()==ed.SourceVerticesM()&&final.Surface().CanonicalSha256()==ed.Surface().CanonicalSha256(),"LV committed final ED landmark surface/hash changed");}CheckCommittedFlowRatePorts(runtime,options,trial_observations);auto snapshot=iga::BuildCommittedMovingImmersedFlowSnapshot(runtime,snapshot_options);const auto wall_gate=PrintWallGate(step,runtime,snapshot_options);CheckSnapshot(snapshot,runtime.CommittedGeometry(),true);Require(snapshot.Request().dt_s==lattice.dt,"LV snapshot did not retain the canonical lattice dt");const double wall_metric_tolerance=512.*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(snapshot.Metrics().wall_relative_velocity_rms_m_per_s),std::abs(wall_gate.quadrature.relative_rms_m_per_s)});Require(std::abs(snapshot.Metrics().wall_relative_velocity_rms_m_per_s-wall_gate.quadrature.relative_rms_m_per_s)<=wall_metric_tolerance,"LV wall quadrature full RMS differs from snapshot metric");max_wall_ratio=std::max(max_wall_ratio,wall_gate.ratio);Require(wall_gate.public_pass,"LV wall-relative velocity gate failed");Require(wall_gate.release_pass,"LV wall-relative velocity release margin failed");
			if(step==2){const auto state_hash=runtime.CommittedGlobalState().HashSha256();iga::MovingImmersedFlowSnapshotPublisher::SetTestFault(iga::MovingImmersedFlowSnapshotPublicationFault::AfterTemporaryVtuWrite);bool failed=false;try{(void)iga::MovingImmersedFlowSnapshotPublisher::Publish(snapshot,stem,Identity(runtime),std::nullopt);}catch(const std::exception&){failed=true;}iga::MovingImmersedFlowSnapshotPublisher::SetTestFault(iga::MovingImmersedFlowSnapshotPublicationFault::None);Require(failed&&!fs::exists(iga::MovingImmersedFlowSnapshotPublisher::EpochDirectory(stem,step))&&runtime.CommittedGlobalState().HashSha256()==state_hash,"LV injected publisher failure left an epoch or changed solver state");}
			auto publication=iga::MovingImmersedFlowSnapshotPublisher::Publish(snapshot,stem,Identity(runtime),std::nullopt);auto retry=iga::MovingImmersedFlowSnapshotPublisher::Publish(snapshot,stem,Identity(runtime),std::nullopt);Require(publication.metrics_sha256==retry.metrics_sha256&&publication.vtu_sha256==retry.vtu_sha256,"LV publisher retry is not idempotent");Print(step,runtime,snapshot,publication);q_positive=q_positive||(snapshot.Metrics().enstrophy_integral_m3_per_s2>0.&&snapshot.Metrics().q_positive_volume_m3>1024.*std::numeric_limits<double>::epsilon()*snapshot.Metrics().quadrature_volume_m3);
		}
		Require(contraction_transition&&expansion_transition,"LV cycle did not exercise both active-layout transition directions");Require(q_positive,"LV cycle did not retain a physical Q-positive rotating state");const auto ed=LatticeEvaluation(motion,0),final=LatticeEvaluation(motion,steps);Require(final.SourceVerticesM()==ed.SourceVerticesM()&&final.Surface().CanonicalSha256()==ed.Surface().CanonicalSha256(),"LV evaluated periodicity changed");
		iga::MovingImmersedFlowSnapshotPublisher::RebuildCollection(stem);std::ifstream pvd(iga::MovingImmersedFlowSnapshotPublisher::CollectionPath(stem));const std::string pvd_text((std::istreambuf_iterator<char>(pvd)),{});std::size_t epochs=0,pos=0;while((pos=pvd_text.find("<DataSet timestep=",pos))!=std::string::npos){++epochs;++pos;}Require(epochs==steps+1,"LV PVD epoch count is not strict");for(std::uint64_t i=0;i<=steps;++i)Require(fs::is_regular_file(iga::MovingImmersedFlowSnapshotPublisher::EpochDirectory(stem,i)/"fields.vtu")&&fs::is_regular_file(iga::MovingImmersedFlowSnapshotPublisher::EpochDirectory(stem,i)/"metrics.json"),"LV PVD epoch index is incomplete");std::cout<<"phase7_summary grid=6x6x7 depth=2 rho=1050 mu=0.012 period=0.8 flow="<<kFlowM3S<<" nominal_dt="<<nominal_dt<<" min_dt="<<minimum_dt<<" max_dt="<<maximum_dt<<" max_linear="<<max_linear<<" max_controller="<<max_controller<<" max_wall_ratio="<<max_wall_ratio<<" max_discrete_continuity="<<max_continuity<<" max_divergence_quadrature_defect="<<max_divergence_quadrature<<" max_reynolds="<<max_reynolds<<" max_moving_mass="<<max_mass<<" max_wall_leakage="<<max_leak<<" elapsed_s="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<'\n'; if(temporary){std::error_code ignored;fs::remove_all(root,ignored);}
		}
	} catch(const std::exception& error){std::cerr<<"phase7_lv_closure: "<<error.what()<<'\n';status=1;}
	catch(...){std::cerr<<"phase7_lv_closure: non-standard exception\n";status=1;}
	if(petsc_initialized&&PetscFinalize()!=0){std::cerr<<"phase7_lv_closure: PetscFinalize failed\n";status=1;}
	return status;
}
