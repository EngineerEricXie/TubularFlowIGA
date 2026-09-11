#include "CompliantChannelFsiFixture.hpp"
#include "CheckedText.hpp"
#include "ExecutionResources.hpp"
#include "MovingImmersedTransientFlowFsiRuntime.hpp"
#include "PretensionedMembraneFsiRuntime.hpp"
#include "StrongFluidStructureCoupling.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <utility>
#include <vector>

namespace {

namespace fs = std::filesystem;

void Require(bool value, const char* message)
{ if (!value) throw std::runtime_error(message); }

double Dot(const std::array<double,3>& left, const std::array<double,3>& right)
{ return left[0]*right[0]+left[1]*right[1]+left[2]*right[2]; }

std::array<double,3> Cross(const std::array<double,3>& left, const std::array<double,3>& right)
{ return {{left[1]*right[2]-left[2]*right[1],left[2]*right[0]-left[0]*right[2],left[0]*right[1]-left[1]*right[0]}}; }

double Norm(const std::array<double,3>& value)
{ return std::sqrt(Dot(value,value)); }

void RequireFinite(double value, const char* message)
{ Require(std::isfinite(value),message); }

void RequireFinite(const std::array<double,3>& value, const char* message)
{ for (const double component : value) RequireFinite(component,message); }

std::string ReadAll(const fs::path& path)
{
	std::ifstream input(path,std::ios::binary);
	if (!input) throw std::runtime_error("cannot read exported ParaView file");
	return iga::ReadCheckedText(input);
}

void RequireContains(const std::string& text, const std::string& token, const char* message)
{ Require(text.find(token)!=std::string::npos,message); }

void Data(std::ostream& output, const char* type, const char* name, unsigned components,
	const std::string& values)
{
	output<<"        <DataArray type=\""<<type<<"\" Name=\""<<name<<"\"";
	if (components!=1) output<<" NumberOfComponents=\""<<components<<"\"";
	output<<" format=\"ascii\">\n          "<<values<<"\n        </DataArray>\n";
}

template <class Writer> std::string Values(Writer&& writer)
{
	std::ostringstream values;
	values.exceptions(std::ios::badbit | std::ios::failbit);
	values<<std::setprecision(17); writer(values); return values.str();
}

void WriteFluidVtu(const fs::path& path, const iga::MovingImmersedFlowSnapshot& snapshot)
{
	std::ofstream output(path,std::ios::binary|std::ios::trunc);
	if (!output) throw std::runtime_error("cannot create fluid ParaView VTU");
	const auto& points=snapshot.Points();
	output<<std::setprecision(17)<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n  <UnstructuredGrid>\n    <FieldData>\n";
	Data(output,"Float64","time_s",1,Values([&](auto& v){v<<snapshot.Request().time_s;}));
	Data(output,"String","field_units",1,"velocity_m_per_s; pressure_pa; vorticity_per_s; integration_weight_m3");
	output<<"    </FieldData>\n    <Piece NumberOfPoints=\""<<points.size()<<"\" NumberOfCells=\""<<points.size()<<"\">\n      <PointData>\n";
	Data(output,"Float64","velocity_m_per_s",3,Values([&](auto& v){for(const auto& p:points)for(double x:p.velocity_m_per_s)v<<x<<' ';}));
	Data(output,"Float64","pressure_pa",1,Values([&](auto& v){for(const auto& p:points)v<<p.pressure<<' ';}));
	Data(output,"Float64","vorticity_per_s",3,Values([&](auto& v){for(const auto& p:points)for(double x:p.vorticity_per_s)v<<x<<' ';}));
	Data(output,"Float64","speed_m_per_s",1,Values([&](auto& v){for(const auto& p:points)v<<p.speed_m_per_s<<' ';}));
	Data(output,"Float64","q_criterion_per_s2",1,Values([&](auto& v){for(const auto& p:points)v<<p.q_criterion_per_s2<<' ';}));
	Data(output,"Float64","enstrophy_density_per_s2",1,Values([&](auto& v){for(const auto& p:points)v<<p.enstrophy_density_per_s2<<' ';}));
	Data(output,"Float64","integration_weight_m3",1,Values([&](auto& v){for(const auto& p:points)v<<p.physical_integration_weight_m3<<' ';}));
	Data(output,"Float64","reference_parametric",3,Values([&](auto& v){for(const auto& p:points)for(double x:p.parametric)v<<x<<' ';}));
	Data(output,"UInt64","background_cell_id",1,Values([&](auto& v){for(const auto& p:points)v<<p.background_cell_id<<' ';}));
	Data(output,"UInt64","quadrature_ordinal",1,Values([&](auto& v){for(const auto& p:points)v<<p.cell_quadrature_ordinal<<' ';}));
	output<<"      </PointData>\n      <Points>\n";
	Data(output,"Float64","points",3,Values([&](auto& v){for(const auto& p:points)for(double x:p.physical_m)v<<x<<' ';}));
	output<<"      </Points>\n      <Cells>\n";
	Data(output,"Int64","connectivity",1,Values([&](auto& v){for(std::size_t i=0;i<points.size();++i)v<<i<<' ';}));
	Data(output,"Int64","offsets",1,Values([&](auto& v){for(std::size_t i=0;i<points.size();++i)v<<(i+1)<<' ';}));
	Data(output,"UInt8","types",1,Values([&](auto& v){for(std::size_t i=0;i<points.size();++i)v<<1<<' ';}));
	output<<"      </Cells>\n    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n";
	output.close(); if (!output) throw std::runtime_error("cannot write fluid ParaView VTU");
}

struct MembraneExport {
	std::vector<std::uint64_t> global_ids;
	std::vector<std::array<double,3>> reference_positions,positions,displacement,velocity,normals,traction,force;
	std::vector<double> normal_displacement,normal_velocity;
	std::vector<std::array<std::uint64_t,3>> triangles;
};

MembraneExport BuildMembraneExport(const iga::MaterialSurfacePatchMap& map,
	const iga::MaterialSurfaceKinematics& material, const iga::PretensionedMembraneFsiRuntime& structure,
	const iga::SurfaceTraction& fluid_traction)
{
	const auto& layout=map.Layout(); const auto state=structure.CommittedStateSnapshot();
	const auto snapshot=structure.StrongCouplingCommittedSnapshot();
	Require(layout.owned_global_node_ids.size()==layout.global_node_count && layout.global_node_count==9,
		"ParaView export requires the exact globally owned compliant membrane layout");
	Require(state.displacement_m.size()==layout.global_node_count && state.velocity_m_per_s.size()==layout.global_node_count
		&& snapshot.immutable_reference_normals.size()==layout.global_node_count
		&& fluid_traction.traction_on_structure_pa.size()==layout.global_node_count
		&& fluid_traction.consistent_nodal_force_n.size()==layout.global_node_count,
		"ParaView membrane fields do not align with the exact global layout");
	MembraneExport result; result.global_ids=layout.owned_global_node_ids;
	for(std::size_t local=0;local<layout.global_node_count;++local) {
		const auto global=layout.owned_global_node_ids[local];
		const auto source=map.SourceVertexForGlobalNode(global);
		const auto& reference=material.ReferenceMaterialVerticesM().at(source);
		const auto& fluid_current=material.SourceVerticesM().at(source);
		const auto& normal=snapshot.immutable_reference_normals[local];
		std::array<double,3> displacement{},velocity{};
		for(int component=0;component<3;++component) {
			displacement[component]=state.displacement_m[local]*normal[component];
			velocity[component]=state.velocity_m_per_s[local]*normal[component];
		}
		std::array<double,3> current{}; for(int component=0;component<3;++component) current[component]=reference[component]+displacement[component];
		RequireFinite(reference,"nonfinite membrane reference coordinate"); RequireFinite(fluid_current,"nonfinite committed fluid-wall coordinate"); RequireFinite(current,"nonfinite membrane current coordinate");
		RequireFinite(displacement,"nonfinite membrane displacement"); RequireFinite(velocity,"nonfinite membrane velocity");
		RequireFinite(normal,"nonfinite membrane normal"); RequireFinite(fluid_traction.traction_on_structure_pa[local],"nonfinite fluid traction");
		RequireFinite(fluid_traction.consistent_nodal_force_n[local],"nonfinite fluid nodal force");
		result.reference_positions.push_back(reference); result.positions.push_back(current); result.displacement.push_back(displacement);
		result.velocity.push_back(velocity); result.normals.push_back(normal); result.traction.push_back(fluid_traction.traction_on_structure_pa[local]);
		result.force.push_back(fluid_traction.consistent_nodal_force_n[local]); result.normal_displacement.push_back(Dot(displacement,normal));
		result.normal_velocity.push_back(Dot(velocity,normal));
	}
	for(const auto& global_triangle:layout.reference_triangles) {
		std::array<std::uint64_t,3> local{};
		for(int corner=0;corner<3;++corner) {
			const auto found=std::lower_bound(layout.owned_global_node_ids.begin(),layout.owned_global_node_ids.end(),global_triangle[corner]);
			Require(found!=layout.owned_global_node_ids.end() && *found==global_triangle[corner],"membrane triangle has an unknown global node");
			local[corner]=static_cast<std::uint64_t>(found-layout.owned_global_node_ids.begin());
		}
		const std::array<double,3> left{{result.positions[local[1]][0]-result.positions[local[0]][0],result.positions[local[1]][1]-result.positions[local[0]][1],result.positions[local[1]][2]-result.positions[local[0]][2]}};
		const std::array<double,3> right{{result.positions[local[2]][0]-result.positions[local[0]][0],result.positions[local[2]][1]-result.positions[local[0]][1],result.positions[local[2]][2]-result.positions[local[0]][2]}};
		Require(Cross(left,right)[2]>0.0,"membrane triangle lost the required outward +z orientation"); result.triangles.push_back(local);
	}
	Require(result.triangles.size()==8,"ParaView membrane triangle count changed"); return result;
}

void WriteMembraneVtu(const fs::path& path, const MembraneExport& membrane, double time_s)
{
	std::ofstream output(path,std::ios::binary|std::ios::trunc); if(!output) throw std::runtime_error("cannot create membrane ParaView VTU");
	output<<std::setprecision(17)<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n  <UnstructuredGrid>\n    <FieldData>\n";
	Data(output,"Float64","time_s",1,Values([&](auto& v){v<<time_s;}));
	Data(output,"String","field_units",1,"coordinates_m; displacement_m; velocity_m_per_s; traction_pa; nodal_force_n");
	output<<"    </FieldData>\n    <Piece NumberOfPoints=\""<<membrane.positions.size()<<"\" NumberOfCells=\""<<membrane.triangles.size()<<"\">\n      <PointData>\n";
	Data(output,"UInt64","global_node_id",1,Values([&](auto& v){for(auto x:membrane.global_ids)v<<x<<' ';}));
	Data(output,"Float64","reference_coordinates_m",3,Values([&](auto& v){for(const auto& x:membrane.reference_positions)for(double y:x)v<<y<<' ';}));
	Data(output,"Float64","displacement_m",3,Values([&](auto& v){for(const auto& x:membrane.displacement)for(double y:x)v<<y<<' ';}));
	Data(output,"Float64","velocity_m_per_s",3,Values([&](auto& v){for(const auto& x:membrane.velocity)for(double y:x)v<<y<<' ';}));
	Data(output,"Float64","reference_normal",3,Values([&](auto& v){for(const auto& x:membrane.normals)for(double y:x)v<<y<<' ';}));
	Data(output,"Float64","normal_displacement_m",1,Values([&](auto& v){for(double x:membrane.normal_displacement)v<<x<<' ';}));
	Data(output,"Float64","normal_velocity_m_per_s",1,Values([&](auto& v){for(double x:membrane.normal_velocity)v<<x<<' ';}));
	Data(output,"Float64","fluid_traction_on_structure_pa",3,Values([&](auto& v){for(const auto& x:membrane.traction)for(double y:x)v<<y<<' ';}));
	Data(output,"Float64","fluid_nodal_force_on_structure_n",3,Values([&](auto& v){for(const auto& x:membrane.force)for(double y:x)v<<y<<' ';}));
	output<<"      </PointData>\n      <Points>\n";
	Data(output,"Float64","points",3,Values([&](auto& v){for(const auto& x:membrane.positions)for(double y:x)v<<y<<' ';}));
	output<<"      </Points>\n      <Cells>\n";
	Data(output,"Int64","connectivity",1,Values([&](auto& v){for(const auto& x:membrane.triangles)for(auto y:x)v<<y<<' ';}));
	Data(output,"Int64","offsets",1,Values([&](auto& v){for(std::size_t i=0;i<membrane.triangles.size();++i)v<<3*(i+1)<<' ';}));
	Data(output,"UInt8","types",1,Values([&](auto& v){for(std::size_t i=0;i<membrane.triangles.size();++i)v<<5<<' ';}));
	output<<"      </Cells>\n    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n";
	output.close(); if(!output) throw std::runtime_error("cannot write membrane ParaView VTU");
}

void WriteCollection(const fs::path& root, double time_s)
{
	std::ofstream vtm(root/"t1.05"/"fsi.vtm",std::ios::binary|std::ios::trunc); if(!vtm) throw std::runtime_error("cannot create ParaView VTM");
	vtm<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\">\n  <vtkMultiBlockDataSet>\n    <Block index=\"0\" name=\"fluid\"><DataSet index=\"0\" file=\"fluid.vtu\"/></Block>\n    <Block index=\"1\" name=\"membrane\"><DataSet index=\"0\" file=\"membrane.vtu\"/></Block>\n  </vtkMultiBlockDataSet>\n</VTKFile>\n";
	vtm.close(); if(!vtm) throw std::runtime_error("cannot write ParaView VTM");
	std::ofstream pvd(root/"phase8_compliant_channel_fsi.pvd",std::ios::binary|std::ios::trunc); if(!pvd) throw std::runtime_error("cannot create ParaView PVD");
	pvd<<std::setprecision(17)<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n  <Collection>\n    <DataSet timestep=\""<<time_s<<"\" group=\"\" part=\"0\" file=\"t1.05/fsi.vtm\"/>\n  </Collection>\n</VTKFile>\n";
	pvd.close(); if(!pvd) throw std::runtime_error("cannot write ParaView PVD");
}

void WriteReadme(const fs::path& root, double time_s, const iga::StrongFluidStructureCouplingResult& result,
	const iga::MovingImmersedTransientFlowConservationDiagnostics& conservation,
	const iga::FluidSurfaceTractionDiagnostics& traction)
{
	std::ofstream output(root/"README.md",std::ios::binary|std::ios::trunc); if(!output) throw std::runtime_error("cannot create ParaView README");
	output<<std::setprecision(17)<<"# Phase 8 compliant-channel FSI ParaView export\n\nOpen `phase8_compliant_channel_fsi.pvd` in ParaView. It contains one solved final frame at t = "<<time_s<<" s.\n\n"
		<<"`t1.05/fsi.vtm` has two blocks: `fluid` and `membrane`. Fluid is the exact committed moving immersed-flow snapshot represented as VTK vertex cells at its volume quadrature samples. Its point fields are `velocity_m_per_s` [m/s], `pressure_pa` [Pa], `vorticity_per_s` [1/s], `speed_m_per_s` [m/s], `q_criterion_per_s2` [1/s^2], `enstrophy_density_per_s2` [1/s^2], `integration_weight_m3` [m^3], `reference_parametric`, `background_cell_id`, and `quadrature_ordinal`.\n\n"
		<<"Membrane coordinates are current displaced coordinates [m]. Its point fields are `global_node_id`, `reference_coordinates_m` [m], `displacement_m` [m], `velocity_m_per_s` [m/s], `reference_normal`, `normal_displacement_m` [m], `normal_velocity_m_per_s` [m/s], `fluid_traction_on_structure_pa` [Pa], and `fluid_nodal_force_on_structure_n` [N]. Connectivity is the exact global-node/triangle map (global IDs 11--19), with every triangle outward +z.\n\n"
		<<"Accepted FSI: iterations="<<result.iterations<<", moving_mass="<<conservation.normalized_moving_mass_defect<<", wall_leakage="<<conservation.normalized_wall_relative_leakage<<", continuity="<<conservation.normalized_discrete_moving_wall_continuity_defect<<", traction_resultant_n="<<Norm(traction.nodal_resultant_n)<<".\n";
	output.close(); if(!output) throw std::runtime_error("cannot write ParaView README");
}

void ValidateExport(const fs::path& root, std::size_t fluid_points, const MembraneExport& membrane)
{
	const auto pvd=root/"phase8_compliant_channel_fsi.pvd",vtm=root/"t1.05"/"fsi.vtm",fluid=root/"t1.05"/"fluid.vtu",structure=root/"t1.05"/"membrane.vtu";
	for(const auto& path:{pvd,vtm,fluid,structure,root/"README.md"}) Require(fs::is_regular_file(path),"ParaView export file is missing");
	const auto pvd_text=ReadAll(pvd),vtm_text=ReadAll(vtm),fluid_text=ReadAll(fluid),membrane_text=ReadAll(structure);
	RequireContains(pvd_text,"<VTKFile type=\"Collection\"","PVD XML type is invalid"); RequireContains(pvd_text,"file=\"t1.05/fsi.vtm\"","PVD relative VTM reference is invalid");
	Require(fs::is_regular_file(pvd.parent_path()/"t1.05/fsi.vtm"),"PVD relative VTM reference does not resolve");
	RequireContains(vtm_text,"<VTKFile type=\"vtkMultiBlockDataSet\"","VTM XML type is invalid"); RequireContains(vtm_text,"file=\"fluid.vtu\"","VTM fluid reference is invalid"); RequireContains(vtm_text,"file=\"membrane.vtu\"","VTM membrane reference is invalid");
	Require(fs::is_regular_file(vtm.parent_path()/"fluid.vtu") && fs::is_regular_file(vtm.parent_path()/"membrane.vtu"),"VTM relative dataset reference does not resolve");
	RequireContains(fluid_text,"NumberOfPoints=\""+std::to_string(fluid_points)+"\" NumberOfCells=\""+std::to_string(fluid_points)+"\"","fluid XML point/cell counts are invalid");
	RequireContains(membrane_text,"NumberOfPoints=\"9\" NumberOfCells=\"8\"","membrane XML point/cell counts are invalid");
	for(const char* field:{"velocity_m_per_s","pressure_pa"}) RequireContains(fluid_text,"Name=\""+std::string(field)+"\"","fluid XML field is missing");
	for(const char* field:{"reference_coordinates_m","displacement_m","velocity_m_per_s","normal_displacement_m","normal_velocity_m_per_s","fluid_traction_on_structure_pa","fluid_nodal_force_on_structure_n"}) RequireContains(membrane_text,"Name=\""+std::string(field)+"\"","membrane XML field is missing");
	Require(membrane.global_ids.size()==9 && membrane.triangles.size()==8,"membrane export mapping count is invalid");
}

} // namespace

int main(int argc, char** argv)
{
	fs::path output_directory; PetscInitialize(&argc,&argv,nullptr,nullptr); int status=0;
	try {
		iga::RequireExecutionResources(PETSC_COMM_WORLD, &std::cout);
		int ranks = 1;
		MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
		Require(ranks == 1, "phase8_compliant_channel_fsi_paraview requires one MPI rank; distributed FSI is not supported by this exporter");
		Require(argc==2,"usage: phase8_compliant_channel_fsi_paraview OUTPUT_DIRECTORY");
		const fs::path root=argv[1]; output_directory=root; Require(!root.empty(),"ParaView output directory is empty"); fs::create_directories(root/"t1.05");
		if(!fs::is_directory(root/"t1.05")) throw std::runtime_error("ParaView output directory is unavailable");
		std::error_code ignored; fs::remove(root/"export_failure.txt",ignored);
		const auto started=std::chrono::steady_clock::now();
		const auto initial=iga::compliant_channel_fixture::InitialMaterial(); const auto map=iga::compliant_channel_fixture::PatchMap(initial); const auto& layout=map.Layout();
		const auto fluid_surface=iga::compliant_channel_fixture::Fluid(map.ReferenceIdentitySha256());
		const iga::FsiCouplingEdge edge("compliant-channel-wall",fluid_surface.id,map.Interface().id,iga::FsiCouplingLaw::FluidStructureTractionKinematics);
		iga::MovingImmersedTransientFlowFsiRuntime fluid("fluid","immersed",edge,fluid_surface,map.Interface(),layout,layout,initial,map,iga::compliant_channel_fixture::FlowOptions());
		iga::PretensionedMembraneFsiRuntime structure("structure","membrane",edge,map.Interface(),fluid_surface,layout,layout,iga::compliant_channel_fixture::MembraneMaterial(),iga::compliant_channel_fixture::ClampedGlobalNodeIds());
		iga::StrongFluidStructureCoupling<iga::MovingImmersedTransientFlowFsiRuntime,iga::PretensionedMembraneFsiRuntime> coordinator(fluid,structure,edge,layout,iga::compliant_channel_fixture::CouplingOptions());
		const auto result=coordinator.Execute({1,1.0,.05});
		Require(result.converged && result.status=="converged" && result.iterations>=2 && result.history.size()==result.iterations,"approved compliant-channel FSI case did not converge");
		Require(result.history.front().area_weighted_rms_residual_m>result.history.front().convergence_threshold_m && result.history.back().area_weighted_rms_residual_m<=result.history.back().convergence_threshold_m && result.history.back().area_weighted_rms_residual_m<result.history.front().area_weighted_rms_residual_m,"approved strong-coupling residual gates failed");
		const auto& flow=fluid.CommittedFlowDiagnostics(); const auto conservation=fluid.ConservationDiagnostics(); const auto traction=fluid.CommittedSurfaceTractionSnapshot(); const auto traction_audit=fluid.CommittedSurfaceTractionDiagnostics();
		Require(flow.converged && flow.committed && !flow.trial_active && flow.identity_history_nodes==flow.active_nodes && flow.missing_history_nodes==0 && flow.nonlinear_iterations>0 && std::isfinite(flow.residual_norm) && std::isfinite(flow.true_linear_relative_residual),"approved committed fluid gate failed");
		Require(traction && !traction->projection_identity_sha256.empty() && traction_audit.nodal_resultant_n[2]>0.0,"committed fluid-on-structure traction is unavailable or misoriented");
		for(int component=0;component<3;++component) RequireFinite(traction_audit.nodal_resultant_n[component],"nonfinite traction resultant");
		Require(std::isfinite(conservation.normalized_moving_mass_defect) && std::isfinite(conservation.normalized_wall_relative_leakage) && std::isfinite(conservation.normalized_discrete_moving_wall_continuity_defect) && conservation.normalized_moving_mass_defect<.03 && conservation.normalized_wall_relative_leakage<.03 && conservation.normalized_discrete_moving_wall_continuity_defect<1.e-8,"approved moving-domain conservation gates failed");
		Require(fluid.CommittedGlobalState().Index()==1 && fluid.CommittedGlobalState().TimeS()==1.05 && fluid.MovingDiagnostics().idle && !fluid.MovingDiagnostics().trial_active && fluid.CommittedCompositionIdentitySha256()==result.final_committed_fluid_composition_identity_sha256 && structure.CommittedStateIdentitySha256()==result.final_committed_structure_state_identity_sha256,"approved FSI epoch finalization gate failed");
		iga::MovingImmersedFlowSnapshotOptions snapshot_options; snapshot_options.maximum_points=2000000; snapshot_options.maximum_output_bytes=500000000;
		const auto snapshot=fluid.CommittedFlowSnapshot(snapshot_options); Require(snapshot.Request().time_s==1.05 && snapshot.Request().index==1 && !snapshot.Points().empty(),"committed fluid snapshot is unavailable");
		for(const auto& point:snapshot.Points()) { RequireFinite(point.physical_m,"nonfinite fluid point coordinate"); RequireFinite(point.velocity_m_per_s,"nonfinite fluid velocity"); RequireFinite(point.pressure,"nonfinite fluid pressure"); RequireFinite(point.vorticity_per_s,"nonfinite fluid vorticity"); RequireFinite(point.speed_m_per_s,"nonfinite fluid speed"); RequireFinite(point.q_criterion_per_s2,"nonfinite fluid Q criterion"); RequireFinite(point.enstrophy_density_per_s2,"nonfinite fluid enstrophy"); RequireFinite(point.physical_integration_weight_m3,"nonfinite fluid integration weight"); }
		const auto membrane=BuildMembraneExport(map,fluid.CommittedMaterialKinematicsSnapshot(),structure,*traction);
		WriteFluidVtu(root/"t1.05"/"fluid.vtu",snapshot); WriteMembraneVtu(root/"t1.05"/"membrane.vtu",membrane,snapshot.Request().time_s); WriteCollection(root,snapshot.Request().time_s); WriteReadme(root,snapshot.Request().time_s,result,conservation,traction_audit); ValidateExport(root,snapshot.Points().size(),membrane);
		struct rusage usage{}; getrusage(RUSAGE_SELF,&usage); const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
		std::cout<<std::setprecision(17)<<"phase8_paraview_export converged=true iterations="<<result.iterations<<" time_s="<<snapshot.Request().time_s<<" fluid_points="<<snapshot.Points().size()<<" membrane_points="<<membrane.positions.size()<<" membrane_triangles="<<membrane.triangles.size()<<" center_displacement_m="<<structure.CommittedStateSnapshot().displacement_m[4]<<" center_velocity_m_s="<<structure.CommittedStateSnapshot().velocity_m_per_s[4]<<" moving_mass="<<conservation.normalized_moving_mass_defect<<" wall_leakage="<<conservation.normalized_wall_relative_leakage<<" continuity="<<conservation.normalized_discrete_moving_wall_continuity_defect<<" traction_resultant_n="<<Norm(traction_audit.nodal_resultant_n)<<" wall_time_s="<<elapsed<<" peak_rss_kib="<<usage.ru_maxrss<<" output="<<root.string()<<"\n";
	} catch(const std::exception& error) {
		if(!output_directory.empty()) { std::ofstream diagnostic(output_directory/"export_failure.txt",std::ios::binary|std::ios::trunc); if(diagnostic) diagnostic<<error.what()<<"\n"; }
		std::cerr<<error.what()<<"\n"; status=1;
	}
	PetscFinalize(); return status;
}
