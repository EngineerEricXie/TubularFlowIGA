#include "MaterialSurfacePatchMap.hpp"
#include <set>
#include <map>
#include "ExecutionResources.hpp"
#include "MovingImmersedTransientFlowFsiRuntime.hpp"
#include "PretensionedMembraneFsiRuntime.hpp"
#include "StrongFluidStructureCoupling.hpp"
#include "NativeFsiReferenceOutput.hpp"

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

void Data(std::ostream& output, const char* type, const char* name, unsigned components,
    const std::string& values)
{
    if (std::string(type)=="String") {
        output<<"        <Array type=\"String\" Name=\""<<name<<"\" NumberOfTuples=\"1\" format=\"ascii\">\n          ";
        for (unsigned char byte:values) output<<static_cast<unsigned>(byte)<<' ';
        output<<"0\n        </Array>\n";
        return;
    }
    output<<"        <DataArray type=\""<<type<<"\" Name=\""<<name<<"\"";
    if (components!=1) output<<" NumberOfComponents=\""<<components<<"\"";
    if (std::string(name)=="time_s" || std::string(name)=="display_displacement_scale") output<<" NumberOfTuples=\"1\"";
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
	Require(layout.owned_global_node_ids.size()==layout.global_node_count,
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
		Require(Norm(Cross(left,right))>0.0,"degenerate deformed membrane triangle"); result.triangles.push_back(local);
	}
	return result;
}

void WriteMembraneVtu(const fs::path& path, const MembraneExport& membrane, double time_s, double display_scale=1.)
{
	std::ofstream output(path,std::ios::binary|std::ios::trunc); if(!output) throw std::runtime_error("cannot create membrane ParaView VTU");
	output<<std::setprecision(17)<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n  <UnstructuredGrid>\n    <FieldData>\n";
	Data(output,"Float64","time_s",1,Values([&](auto& v){v<<time_s;}));
	Data(output,"Float64","display_displacement_scale",1,Values([&](auto& v){v<<display_scale;}));
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
	Data(output,"Float64","points",3,Values([&](auto& v){for(std::size_t i=0;i<membrane.positions.size();++i)for(int c=0;c<3;++c)v<<membrane.reference_positions[i][c]+display_scale*membrane.displacement[i][c]<<' ';}));
	output<<"      </Points>\n      <Cells>\n";
	Data(output,"Int64","connectivity",1,Values([&](auto& v){for(const auto& x:membrane.triangles)for(auto y:x)v<<y<<' ';}));
	Data(output,"Int64","offsets",1,Values([&](auto& v){for(std::size_t i=0;i<membrane.triangles.size();++i)v<<3*(i+1)<<' ';}));
	Data(output,"UInt8","types",1,Values([&](auto& v){for(std::size_t i=0;i<membrane.triangles.size();++i)v<<5<<' ';}));
	output<<"      </Cells>\n    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n";
	output.close(); if(!output) throw std::runtime_error("cannot write membrane ParaView VTU");
}


void WriteNativeReference(const fs::path& directory,
    const iga::MovingImmersedTransientFlowFsiRuntime& fluid,
    const iga::PretensionedMembraneFsiRuntime& structure,
    const iga::MaterialSurfacePatchMap& map,const iga::SurfaceTraction& traction,
    bool native_gates_passed)
{
    iga::NativeFsiReferenceOutput output(directory);
    const auto& state=fluid.CommittedGlobalState();
    output.Add("fluid_velocity","m/s",state.NodeIds(),3,
        [&](std::size_t r,std::size_t c){return state.Coefficients().at(r).at(c);});
    output.Add("fluid_pressure","pa",state.NodeIds(),1,
        [&](std::size_t r,std::size_t){return state.Coefficients().at(r).at(3);});
    if(!state.PortIds().empty()) output.Add("controller_pressure","pa",state.PortIds(),1,
        [&](std::size_t r,std::size_t){return state.PortMultipliers().at(r);});
    if(state.HasGaugeMultiplier()) output.Add("gauge_multiplier","1/s",std::vector<int>{0},1,
        [&](std::size_t,std::size_t){return state.GaugeMultiplier();});
    const auto membrane=structure.CommittedStateSnapshot(); const auto& ids=map.Layout().owned_global_node_ids;
    output.Add("membrane_displacement","m",ids,1,
        [&](std::size_t r,std::size_t){return membrane.displacement_m.at(r);});
    output.Add("membrane_velocity","m/s",ids,1,
        [&](std::size_t r,std::size_t){return membrane.velocity_m_per_s.at(r);});
    output.Add("surface_traction","pa",ids,3,
        [&](std::size_t r,std::size_t c){return traction.traction_on_structure_pa.at(r).at(c);});
    output.Add("surface_force","n",ids,3,
        [&](std::size_t r,std::size_t c){return traction.consistent_nodal_force_n.at(r).at(c);});
    const auto full=fluid.CommittedMaterialKinematicsSnapshot();
    std::vector<std::uint64_t> vertices(full.SourceVerticesM().size());
    for(std::size_t r=0;r<vertices.size();++r) vertices[r]=r;
    output.Add("material_position","m",vertices,3,
        [&](std::size_t r,std::size_t c){return full.SourceVerticesM().at(r).at(c);});
    output.Add("material_displacement","m",vertices,3,
        [&](std::size_t r,std::size_t c){return full.SourceVerticesM().at(r).at(c)-full.ReferenceMaterialVerticesM().at(r).at(c);});
    output.Add("material_velocity","m/s",vertices,3,
        [&](std::size_t r,std::size_t c){return full.SourceVertexVelocitiesMPerS().at(r).at(c);});
    const auto& diagnostics=fluid.CommittedFlowDiagnostics(); std::vector<int> port_ids;
    for(const auto& p:diagnostics.ports) { Require(p.measurement_valid,"native reference has unmeasured committed port"); port_ids.push_back(p.boundary_label); }
    if(!port_ids.empty()) {
        output.Add("port_flow","m^3/s",port_ids,1,
            [&](std::size_t r,std::size_t){return diagnostics.ports.at(r).measurement.outward_flow_m3_s;});
        output.Add("port_mean_pressure","pa",port_ids,1,
            [&](std::size_t r,std::size_t){return diagnostics.ports.at(r).measurement.mean_pressure_pa;});
        output.Add("port_mean_normal_traction","pa",port_ids,1,
            [&](std::size_t r,std::size_t){return diagnostics.ports.at(r).measurement.mean_normal_traction_pa;});
        output.Add("port_area","m^2",port_ids,1,
            [&](std::size_t r,std::size_t){return diagnostics.ports.at(r).measurement.area_m2;});
    }
    std::vector<int> iterations; for(const auto& n:diagnostics.newton_steps) iterations.push_back(n.iteration);
    if(!iterations.empty()) {
        output.Add("newton_residual","mixed",iterations,1,
            [&](std::size_t r,std::size_t){return diagnostics.newton_steps.at(r).residual_norm;});
        output.Add("newton_update","mixed",iterations,1,
            [&](std::size_t r,std::size_t){return diagnostics.newton_steps.at(r).update_norm;});
        output.Add("newton_damping","1",iterations,1,
            [&](std::size_t r,std::size_t){return diagnostics.newton_steps.at(r).damping;});
        output.Add("newton_ksp_iterations","1",iterations,1,
            [&](std::size_t r,std::size_t){return static_cast<double>(diagnostics.newton_steps.at(r).ksp_iterations);});
        output.Add("newton_ksp_reason","1",iterations,1,
            [&](std::size_t r,std::size_t){return static_cast<double>(diagnostics.newton_steps.at(r).ksp_reason);});
        output.Add("newton_linear_relative_residual","1",iterations,1,
            [&](std::size_t r,std::size_t){return diagnostics.newton_steps.at(r).linear_relative_residual;});
    }
    output.Finish(state.TimeS(),state.Index(),native_gates_passed);
}

iga::MaterialSurfacePatchMap ReadVessel(const fs::path& path)
{
    std::ifstream in(path); std::size_t nv=0,nt=0; in>>nv>>nt;
    Require(in && nv>0 && nv<10000 && nt>0 && nt<30000,"invalid vessel mesh header");
    std::vector<std::array<double,3>> points(nv),velocity(nv,{{0,0,0}});
    for(auto& p:points)for(auto& x:p)in>>x;
    std::vector<iga::RawSurfaceTriangle> faces(nt);
    for(auto& f:faces)in>>f.indices[0]>>f.indices[1]>>f.indices[2]>>f.boundary_id;
    Require(bool(in),"invalid vessel mesh records");
    const auto initial=iga::MaterialSurfaceKinematics::CreateFromSourceTopology(points,points,velocity,faces,0.,-1.,0.);
    std::set<std::uint64_t> ids,clamps;
    std::map<std::uint64_t,double> areas;
    iga::DistributedSurfaceLayout layout;layout.partition_count=1;layout.partition_rank=0;
    std::vector<std::uint32_t> source_triangles;
    for(std::size_t i=0;i<nt;++i)if(faces[i].boundary_id==7) {
        const auto& f=faces[i];std::array<std::uint64_t,3> tri;
        for(int c=0;c<3;++c){tri[c]=f.indices[c];ids.insert(tri[c]);}
        std::array<double,3> a{},b{};
        for(int c=0;c<3;++c){a[c]=points[tri[1]][c]-points[tri[0]][c];b[c]=points[tri[2]][c]-points[tri[0]][c];}
        const double area=Norm(Cross(a,b))*.5;
        for(auto id:tri)areas[id]+=area/3.;
        layout.reference_triangles.push_back(tri);source_triangles.push_back(i);
    }
    for(const auto& f:faces)if(f.boundary_id!=7)for(auto id:f.indices)if(ids.count(id))clamps.insert(id);
    std::vector<iga::MaterialSurfacePatchMap::GlobalToSourceVertex> mapping;
    for(auto id:ids){layout.owned_global_node_ids.push_back(id);layout.reference_positions.push_back({id,points[id]});layout.owned_reference_lumped_areas_m2.push_back(areas[id]);mapping.push_back({id,static_cast<std::uint32_t>(id)});}
    layout.global_node_count=ids.size();const std::vector<std::uint64_t> fixed(clamps.begin(),clamps.end());
    layout.reference_mesh_identity_sha256=iga::MaterialSurfacePatchMap::BuildReferenceIdentitySha256(initial,7,mapping,layout.reference_triangles,source_triangles,fixed);
    layout.layout_identity_sha256=iga::BuildDistributedSurfaceLayoutIdentitySha256(layout);
    iga::DistributedSurfaceInterface surface;surface.id={"structure","membrane","wall"};surface.subsystem_id="membrane";surface.boundary_labels={7};surface.reference_mesh_identity_sha256=layout.reference_mesh_identity_sha256;
    surface.provides={iga::SurfaceFieldQuantity::Displacement,iga::SurfaceFieldQuantity::Velocity};surface.requires={iga::SurfaceFieldQuantity::TractionOnStructure};
    std::cout<<"vessel vertices="<<nv<<" triangles="<<nt<<" membrane_nodes="<<ids.size()<<" clamps="<<clamps.size()<<std::endl;
    return iga::MaterialSurfacePatchMap::Create(surface,layout,initial,7,mapping,source_triangles,fixed);
}

void Collection(const fs::path& root,const std::vector<double>& times)
{
    std::ofstream out(root/"fsi.pvd");
    out<<"<?xml version=\"1.0\"?><VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\"><Collection>\n";
    for(std::size_t i=1;i<=times.size();++i) {
        const std::string directory="step-"+std::to_string(i);
        std::ofstream frame(root/directory/"fsi.vtm");
        frame<<"<?xml version=\"1.0\"?><VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\"><vtkMultiBlockDataSet>\n"
             <<"<DataSet index=\"0\" name=\"fluid\" file=\"fluid.vtu\"/>\n"
             <<"<DataSet index=\"1\" name=\"membrane\" file=\"membrane.vtu\"/>\n"
             <<"</vtkMultiBlockDataSet></VTKFile>\n";
        frame.close();Require(bool(frame),"cannot write FSI multiblock frame");
        out<<std::setprecision(17)<<"<DataSet timestep=\""<<times[i-1]<<"\" group=\"\" part=\"0\" file=\""<<directory<<"/fsi.vtm\"/>\n";
    }
    out<<"</Collection></VTKFile>\n";
    out.close();Require(bool(out),"cannot write FSI time collection");
    std::ofstream wall(root/"wall-display.pvd");
    wall<<"<?xml version=\"1.0\"?><VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\"><Collection>\n";
    for(std::size_t i=0;i<times.size();++i)
        wall<<std::setprecision(17)<<"<DataSet timestep=\""<<times[i]<<"\" file=\"step-"<<i+1<<"/membrane-display.vtu\"/>\n";
    wall<<"</Collection></VTKFile>\n";
    wall.close();Require(bool(wall),"cannot write display collection");
}
} // namespace

int main(int argc,char** argv)
{
    fs::path output_directory;
    int provided=0;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);PetscInitialize(&argc,&argv,nullptr,nullptr);int status=0;
    try {
        iga::CurrentPhaseProfile().EnableFromEnvironment();
        std::cout<<std::unitbuf;
        Require(argc>=3,"usage: bifurcation_fsi SURFACE.txt OUTPUT [-demo_steps N] [-demo_grid N] [-demo_dt S] [-demo_pressure PA] [-demo_wall_inertial_gamma G]");
        int ranks=0;MPI_Comm_size(PETSC_COMM_WORLD,&ranks);Require(ranks==1,"this bounded FSI demonstration uses one MPI rank");
        iga::RequireExecutionResources(PETSC_COMM_WORLD,&std::cout);
        PetscInt steps=8,grid=16,depth=3;PetscReal dt=.01,pressure=20.,wall_inertial_gamma=1.,amplitude=0.,period=1.,display_scale=100.;
        PetscBool visualization_only=PETSC_FALSE,native_reference=PETSC_FALSE;
        PetscOptionsGetInt(nullptr,nullptr,"-demo_steps",&steps,nullptr);PetscOptionsGetInt(nullptr,nullptr,"-demo_grid",&grid,nullptr);
        PetscOptionsGetReal(nullptr,nullptr,"-demo_dt",&dt,nullptr);PetscOptionsGetReal(nullptr,nullptr,"-demo_pressure",&pressure,nullptr);
        PetscOptionsGetReal(nullptr,nullptr,"-demo_wall_inertial_gamma",&wall_inertial_gamma,nullptr);
        PetscOptionsGetReal(nullptr,nullptr,"-demo_pulse_amplitude",&amplitude,nullptr);
        PetscOptionsGetReal(nullptr,nullptr,"-demo_period",&period,nullptr);
        PetscOptionsGetReal(nullptr,nullptr,"-demo_display_scale",&display_scale,nullptr);
        PetscOptionsGetInt(nullptr,nullptr,"-demo_cut_depth",&depth,nullptr);
        PetscOptionsGetBool(nullptr,nullptr,"-demo_visualization_only",&visualization_only,nullptr);
        PetscOptionsGetBool(nullptr,nullptr,"-demo_native_reference",&native_reference,nullptr);
        Require(steps>0&&grid>=6&&std::isfinite(dt)&&dt>0&&std::isfinite(pressure)&&pressure>0,"invalid demo parameters");
        Require(std::isfinite(amplitude)&&amplitude>=0&&std::isfinite(period)&&period>0&&std::isfinite(display_scale)&&display_scale>0&&depth>=1&&depth<=3,"invalid pulse/display parameters");
        const auto inlet_pressure=[=](double t){return pressure+amplitude*.5*(1.-std::cos(2.*std::acos(-1.)*t/period));};
        Require(std::isfinite(wall_inertial_gamma)&&wall_inertial_gamma>=0,"invalid wall inertial gamma");
        const fs::path root=argv[2];Require(!fs::exists(root),"output directory already exists");fs::create_directories(root);output_directory=root;
        const auto map=ReadVessel(argv[1]);const auto& layout=map.Layout();const auto initial=map.FullReference();
        auto interface=map.Interface();interface.id={"fluid","immersed","wall"};interface.subsystem_id="immersed";
        interface.provides={iga::SurfaceFieldQuantity::TractionOnStructure};interface.requires={iga::SurfaceFieldQuantity::Displacement,iga::SurfaceFieldQuantity::Velocity};
        const iga::FsiCouplingEdge edge("y-vessel-wall",interface.id,map.Interface().id,iga::FsiCouplingLaw::FluidStructureTractionKinematics);
        iga::MovingImmersedTransientFlowOptions flow;
        flow.grid={{{-.011,-.010,-.0035}},{{.011,.010,.0035}},{{static_cast<std::uint32_t>(grid),static_cast<std::uint32_t>(grid),static_cast<std::uint32_t>(std::max(4,static_cast<int>(grid)/3))}}};
        flow.geometry.volume.max_depth=depth;
        if(depth<3) flow.geometry.volume.empty_rule_rescue_max_depth=3;
        flow.geometry.volume.max_nodes=flow.geometry.volume.max_leaves=flow.geometry.volume.max_points=4000000;
        flow.geometry.volume_storage=iga::CutCellVolumeQuadratureStorageMode::Compact;
        flow.flow.parameters={1060.,.0035,1.};flow.flow.wall_labels={7};
        flow.flow.wall_inertial_gamma0=wall_inertial_gamma;
        std::cout<<"wall stabilization viscous_gamma="<<flow.flow.wall_gamma0<<" inertial_gamma="<<wall_inertial_gamma<<std::endl;
        flow.flow.ports={{"inlet",1,iga::ImmersedFlowPortControlMode::Pressure,pressure},{"daughter_lower",2,iga::ImmersedFlowPortControlMode::Pressure,0.},{"daughter_upper",3,iga::ImmersedFlowPortControlMode::Pressure,0.}};
        flow.flow.flow_controller_reference_flow_m3_s=1e-6;
        flow.flow.ksp_relative_tolerance=1e-12;flow.flow.lu_pivot_shift=1e-20;
        flow.flow.nonlinear_maximum_iterations=24;
        iga::MovingImmersedTransientFlowFsiRuntime fluid("fluid","immersed",edge,interface,map.Interface(),layout,layout,initial,map,flow);
        iga::PretensionedMembraneMaterial material{.3,30.,10.,2e7};
        iga::PretensionedMembraneFsiRuntime structure("structure","membrane",edge,map.Interface(),interface,layout,layout,material,map.ConfiguredClampedGlobalNodeIds());
        iga::StrongFluidStructureCouplingOptions controls;controls.maximum_iterations=32;controls.absolute_displacement_tolerance_m=1e-10;controls.relative_displacement_tolerance=1e-3;
        controls.aitken_controls.initial_relaxation=.5;controls.aitken_controls.minimum_relaxation=.1;controls.aitken_controls.maximum_relaxation=.8;
        iga::StrongFluidStructureCoupling<iga::MovingImmersedTransientFlowFsiRuntime,iga::PretensionedMembraneFsiRuntime> coordinator(fluid,structure,edge,layout,controls);
        std::ofstream history(root/"history.csv");history<<"step,time_s,iterations,max_displacement_m,mass_defect,wall_leakage,continuity,inlet_m3_s,lower_m3_s,upper_m3_s,normalization_scale_m3_s,moving_mass_defect_m3_s,wall_relative_leakage_m3_s,backward_euler_volume_rate_m3_s,total_fluid_surface_outward_flow_m3_s,total_material_surface_outward_flow_m3_s,inlet_pressure_pa,conservation_passed\n";
        std::ofstream coupling(root/"coupling.csv");coupling<<"step,iteration,rms_m,threshold_m\n";
        std::vector<double> accepted_times;
        int conservation_failures=0;
        for(int i=1;i<=steps;++i) {
            const double start_time=fluid.CommittedGlobalState().TimeS();
            std::cout<<"FSI step="<<i<<" start time="<<start_time<<std::endl;
            fluid.SetPortControlValue("inlet",inlet_pressure(start_time+dt));
            const auto result=coordinator.Execute({i,start_time,dt});
            for(const auto& h:result.history) { coupling<<std::setprecision(17)<<i<<','<<h.iteration<<','<<h.area_weighted_rms_residual_m<<','<<h.convergence_threshold_m<<'\n'; }
            coupling.flush();
            Require(result.converged,"Y-vessel strong coupling failed to converge");
            const auto c=fluid.ConservationDiagnostics();const auto traction=fluid.CommittedSurfaceTractionSnapshot();Require(traction.has_value(),"missing traction");
            const auto membrane=BuildMembraneExport(map,fluid.CommittedMaterialKinematicsSnapshot(),structure,*traction);
            double displacement=0;for(auto d:membrane.normal_displacement)displacement=std::max(displacement,std::abs(d));
            const bool conservation_passed=c.normalized_moving_mass_defect<.03&&c.normalized_wall_relative_leakage<.03&&c.normalized_discrete_moving_wall_continuity_defect<1e-8;
            Require(std::isfinite(c.normalized_moving_mass_defect)&&std::isfinite(c.normalized_wall_relative_leakage)&&std::isfinite(c.normalized_discrete_moving_wall_continuity_defect),"nonfinite conservation diagnostic");
            Require(std::isfinite(displacement)&&displacement>0&&displacement<.0001,"nontrivial small-displacement gate failed");
            if(!conservation_passed) {
                ++conservation_failures;
                Require(visualization_only,"moving conservation gate failed");
                std::cout<<"VISUALIZATION ONLY: conservation gate exceeded at step "<<i<<std::endl;
            }
            const auto directory=root/("step-"+std::to_string(i));fs::create_directory(directory);
            iga::MovingImmersedFlowSnapshotOptions snapshot_options;snapshot_options.maximum_points=4000000;snapshot_options.maximum_output_bytes=1000000000;
            WriteFluidVtu(directory/"fluid.vtu",fluid.CommittedFlowSnapshot(snapshot_options));WriteMembraneVtu(directory/"membrane.vtu",membrane,fluid.CommittedGlobalState().TimeS());
            WriteMembraneVtu(directory/"membrane-display.vtu",membrane,fluid.CommittedGlobalState().TimeS(),display_scale);
            accepted_times.push_back(fluid.CommittedGlobalState().TimeS());Collection(root,accepted_times);
            history<<std::setprecision(17)<<i<<','<<accepted_times.back()<<','<<result.iterations<<','<<displacement<<','<<c.normalized_moving_mass_defect<<','<<c.normalized_wall_relative_leakage<<','<<c.normalized_discrete_moving_wall_continuity_defect;
            for(int label:{1,2,3}) { history<<','<<c.fluid_surface_outward_flow_by_boundary_label_m3_s.at(label); }
            history<<','<<c.normalization_scale_m3_s<<','<<c.moving_mass_defect_m3_s<<','<<c.wall_relative_leakage_m3_s
                   <<','<<c.backward_euler_volume_rate_m3_s<<','<<c.total_fluid_surface_outward_flow_m3_s<<','<<c.total_material_surface_outward_flow_m3_s;
            history<<','<<inlet_pressure(start_time+dt)<<','<<(conservation_passed?1:0)<<'\n';history.flush();
            if(!visualization_only) Require(c.fluid_surface_outward_flow_by_boundary_label_m3_s.at(1)<0. && c.fluid_surface_outward_flow_by_boundary_label_m3_s.at(2)>0. && c.fluid_surface_outward_flow_by_boundary_label_m3_s.at(3)>0.,"expected inlet and two daughter outflows are absent");
            Require(std::isfinite(displacement)&&displacement>0&&displacement<.0001,"nontrivial small-displacement gate failed");

            if(native_reference) {
                iga::PhaseScope output_phase(iga::ProfilePhase::Output);
                const bool direction_passed=c.fluid_surface_outward_flow_by_boundary_label_m3_s.at(1)<0.
                    && c.fluid_surface_outward_flow_by_boundary_label_m3_s.at(2)>0.
                    && c.fluid_surface_outward_flow_by_boundary_label_m3_s.at(3)>0.;
                WriteNativeReference(directory/"native",fluid,structure,map,*traction,conservation_passed&&direction_passed);
            }
            std::cout<<"FSI accepted step="<<i<<" iterations="<<result.iterations<<" displacement_m="<<displacement<<" mass="<<c.normalized_moving_mass_defect<<" leakage="<<c.normalized_wall_relative_leakage<<" continuity="<<c.normalized_discrete_moving_wall_continuity_defect<<std::endl;
            iga::CurrentPhaseProfile().Write(std::cout,0,1,0);
        }
        std::ofstream manifest(root/"run.json");
        manifest<<std::setprecision(17)<<"{\"status\":\""<<(visualization_only?"visualization_only":"passed")<<"\",\"conservation_failed_steps\":"<<conservation_failures<<",\"pulse_amplitude_pa\":"<<amplitude<<",\"period_s\":"<<period<<",\"display_scale\":"<<display_scale<<",\"cut_depth\":"<<depth<<",\"steps\":"<<steps<<",\"dt_s\":"<<dt<<",\"grid_xy\":"<<grid<<",\"inlet_pressure_pa\":"<<pressure<<",\"density_kg_m3\":1060,\"viscosity_pa_s\":0.0035,\"wall_gamma0\":2,\"wall_inertial_gamma0\":"<<wall_inertial_gamma<<",\"model\":\"small-displacement pretensioned membrane\"}\n";
        manifest.close();Require(bool(manifest),"cannot write passed run manifest");
        std::cout<<"bifurcation_fsi completed steps="<<steps<<std::endl;
    }catch(const std::exception& e){
        if(!output_directory.empty()){std::ofstream failure(output_directory/"failure.txt");failure<<e.what()<<'\n';}
        std::cerr<<e.what()<<std::endl;status=1;
    }
    PetscFinalize();MPI_Finalize();return status;
}
