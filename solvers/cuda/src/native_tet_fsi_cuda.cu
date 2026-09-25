#define main native_tet_flow_cuda_standalone_main
#include "native_tet_flow_cuda.cu"
#undef main

#include "NativeTetAleConservation.hpp"
#include "NativeTetAleMeshMotion.hpp"
#include "NativeTetLinearElasticVisualization.hpp"
#include "NativeTetMatchingFsiInterface.hpp"
#include "NativeTetSolidStaticSolver.hpp"

namespace {

struct FsiCase
{
	std::filesystem::path fluid_mesh,solid_mesh;
	int interface_label=0,steps=0;
	double dt=0.,density=0.,viscosity=0.,young=0.,poisson=0.,solid_density=0.;
	std::map<int,std::array<double,3>> velocity_by_label;
	std::set<int> outlet_labels;
	std::vector<std::size_t> fixed_nodes;
};

FsiCase ReadFsiCase(const std::filesystem::path& path)
{
	using namespace iga::config_detail;
	std::ifstream stream(path);
	if(!stream)throw std::runtime_error("cannot open FEM CUDA FSI case");
	const auto parsed=JsonParser(iga::ReadCheckedText(stream)).Parse();
	const auto& root=RequireObject(parsed,"FEM CUDA FSI case");
	auto field=[](const auto& object,const char* key)->const JsonValue& {
		const auto* value=Find(object,key);
		if(!value)throw std::invalid_argument(std::string("FSI case requires ")+key);
		return *value;
	};
	auto number=[&](const auto& object,const char* key){
		return RequireNumber(field(object,key),key);};
	if(RequireInteger(field(root,"schema_version"),"schema_version")!=1)
		throw std::invalid_argument("unsupported FEM CUDA FSI schema");
	FsiCase result;
	result.fluid_mesh=path.parent_path()/RequireString(field(root,"fluid_mesh_file"),
		"fluid_mesh_file");
	result.solid_mesh=path.parent_path()/RequireString(field(root,"solid_mesh_file"),
		"solid_mesh_file");
	result.interface_label=RequireInteger(field(root,"interface_label"),"interface_label");
	const auto& fluid=RequireObject(field(root,"fluid"),"fluid");
	result.density=number(fluid,"density_kg_m3");
	result.viscosity=number(fluid,"dynamic_viscosity_pa_s");
	if(const auto* boundaries=Find(fluid,"boundary_velocity_by_label_m_s"))
		for(const auto& item:RequireObject(*boundaries,"boundary velocity")){
			const auto& vector=RequireArray(item.second,"boundary velocity");
			std::array<double,3> value{};
			for(int axis=0;axis<3;++axis)
				value[axis]=RequireNumber(vector[axis],"boundary velocity component");
			result.velocity_by_label.emplace(std::stoi(item.first),value);
		}
	if(const auto* labels=Find(fluid,"natural_boundary_labels"))
		for(const auto& label:RequireArray(*labels,"natural boundary labels"))
			result.outlet_labels.insert(RequireInteger(label,"natural boundary label"));
	const auto& solid=RequireObject(field(root,"solid"),"solid");
	result.young=number(solid,"young_modulus_pa");
	result.poisson=number(solid,"poisson_ratio");
	result.solid_density=number(solid,"density_kg_m3");
	for(const auto& node:RequireArray(field(solid,"fixed_nodes"),"fixed_nodes"))
		result.fixed_nodes.push_back(RequireInteger(node,"fixed node"));
	const auto& time=RequireObject(field(root,"time"),"time");
	result.dt=number(time,"dt_s");
	result.steps=RequireInteger(field(time,"steps"),"steps");
	return result;
}

iga::NativeTetMesh ReadTet(const std::filesystem::path& path)
{
	std::ifstream stream(path);
	if(!stream)throw std::runtime_error("cannot open FEM CUDA FSI mesh");
	return iga::ReadNativeTetMeshGmsh41(stream);
}

std::map<std::array<double,3>,std::uint32_t> InterfaceCoordinates(
	const iga::NativeTetMesh& mesh,int label)
{
	std::map<std::array<double,3>,std::uint32_t> result;
	for(const auto& face:mesh.boundary_triangles)
		if(face.boundary_label==label)
			for(const auto node:face.nodes)result.emplace(mesh.points[node],node);
	return result;
}

std::vector<double> InterfaceForces(const iga::NativeTetMesh& current,
	const iga::NativeTetMesh& solid,
	const std::map<std::uint32_t,std::uint32_t>& fluid_to_solid,int label,
	const std::vector<iga::NativeTetInterfaceTriangleTraction>& tractions)
{
	std::map<std::uint64_t,std::array<double,3>> by_id;
	for(const auto& traction:tractions)by_id.emplace(traction.triangle_id,
		traction.traction_on_structure_pa);
	std::vector<double> force(3*solid.points.size(),0.);
	for(const auto& face:current.boundary_triangles){
		if(face.boundary_label!=label)continue;
		const auto& a=current.points[face.nodes[0]];
		const auto& b=current.points[face.nodes[1]];
		const auto& c=current.points[face.nodes[2]];
		std::array<double,3> ab{{b[0]-a[0],b[1]-a[1],b[2]-a[2]}};
		std::array<double,3> ac{{c[0]-a[0],c[1]-a[1],c[2]-a[2]}};
		const double area=.5*iga::NativeTetInterfaceNorm(
			iga::NativeTetInterfaceCross(ab,ac));
		for(const auto fluid_node:face.nodes){
			const auto solid_node=fluid_to_solid.at(fluid_node);
			for(int axis=0;axis<3;++axis)
				force[3*solid_node+axis]+=area*by_id.at(face.id)[axis]/3.;
		}
	}
	return force;
}

iga::NativeTetMesh Deformed(const iga::NativeTetMesh& reference,
	const std::vector<std::array<double,3>>& displacement)
{
	auto result=reference;
	for(std::size_t node=0;node<result.points.size();++node)
		for(int axis=0;axis<3;++axis)
			result.points[node][axis]+=displacement[node][axis];
	return result;
}

void RunFsi(const std::filesystem::path& case_path,
	const std::filesystem::path& output)
{
	const auto scenario=ReadFsiCase(case_path);
	const auto reference=ReadTet(scenario.fluid_mesh);
	const auto solid=ReadTet(scenario.solid_mesh);
	const auto topology=iga::BuildNativeTaylorHoodTopology(reference);
	const auto rows=ElementRows(reference,topology);
	const int cells=static_cast<int>(reference.cells.size());
	const int velocity_nodes=static_cast<int>(reference.points.size()+topology.edges.size());
	const int dofs=3*velocity_nodes+static_cast<int>(reference.points.size());
	const auto fluid_interface=InterfaceCoordinates(reference,scenario.interface_label);
	const auto solid_interface=InterfaceCoordinates(solid,scenario.interface_label);
	if(fluid_interface.size()!=solid_interface.size())
		throw std::invalid_argument("FEM CUDA FSI interface sizes differ");
	std::map<std::uint32_t,std::uint32_t> fluid_to_solid;
	for(const auto& point:fluid_interface)
		fluid_to_solid.emplace(point.second,solid_interface.at(point.first));
	std::vector<double> state(dofs,0.),solid_displacement(3*solid.points.size(),0.);
	std::vector<std::array<double,3>> ale_displacement(reference.points.size()),
		previous_ale=ale_displacement;
	std::vector<int> constrained(dofs,0);
	std::vector<double> prescribed(dofs,0.),mesh_velocity(cells*12,0.);
	std::array<double,12> acceleration{};
	iga::cuda::DeviceBuffer<iga::cuda::NativeTetFlowQuadrature> device_quadrature(
		static_cast<std::size_t>(cells)*iga::cuda::kTetFlowQuadrature);
	iga::cuda::DeviceBuffer<iga::cuda::NativeTetFlowField> fields(
		static_cast<std::size_t>(cells)*iga::cuda::kTetFlowQuadrature);
	iga::cuda::DeviceBuffer<int> device_rows(rows.size()),device_constrained(dofs);
	iga::cuda::DeviceBuffer<double> device_state(dofs),device_previous(dofs);
	iga::cuda::DeviceBuffer<double> device_mesh_velocity(mesh_velocity.size());
	iga::cuda::DeviceBuffer<double> device_acceleration(acceleration.size());
	iga::cuda::DeviceBuffer<double> device_element_residual(cells*34);
	iga::cuda::DeviceBuffer<double> device_element_jacobian(cells*34*34);
	iga::cuda::DeviceBuffer<double> device_residual(dofs),device_jacobian(
		static_cast<std::size_t>(dofs)*dofs),device_rhs(dofs);
	iga::cuda::DeviceBuffer<double> device_prescribed(dofs),device_boundary_load(dofs);
	device_rows.CopyFromHost(rows.data(),rows.size());
	device_state.CopyFromHost(state.data(),state.size());
	device_acceleration.CopyFromHost(acceleration.data(),acceleration.size());
	DenseSolver solver(dofs,device_jacobian.data());
	std::map<std::size_t,double> solid_constraints;
	for(const auto node:scenario.fixed_nodes)
		for(int axis=0;axis<3;++axis)solid_constraints.emplace(3*node+axis,0.);
	iga::NativeTetSolidStaticOptions solid_options;
	solid_options.maximum_iterations=40;
	std::vector<std::pair<double,std::filesystem::path>> flow_series,wall_series;
	std::filesystem::create_directories(output);
	double assembly_seconds=0.,solve_seconds=0.;
	for(int step=1;step<=scenario.steps;++step){
		const auto current=Deformed(reference,ale_displacement);
		const auto quadrature=Quadrature(current);
		std::fill(constrained.begin(),constrained.end(),0);
		std::fill(prescribed.begin(),prescribed.end(),0.);
		for(const auto& boundary:topology.boundary_velocity_nodes){
			if(boundary.first==scenario.interface_label
				||scenario.outlet_labels.count(boundary.first))continue;
			const auto value=scenario.velocity_by_label.count(boundary.first)
				?scenario.velocity_by_label.at(boundary.first)
				:std::array<double,3>{{0.,0.,0.}};
			for(const auto node:boundary.second)
				for(int axis=0;axis<3;++axis){
					constrained[3*node+axis]=1;
					prescribed[3*node+axis]=value[axis];
				}
		}
		for(const auto node:topology.boundary_velocity_nodes.at(
			scenario.interface_label)){
			std::array<double,3> velocity{};
			if(node<reference.points.size())
				for(int axis=0;axis<3;++axis)velocity[axis]=
					(ale_displacement[node][axis]-previous_ale[node][axis])/scenario.dt;
			else{
				const auto edge=topology.edges[node-reference.points.size()];
				for(int axis=0;axis<3;++axis)velocity[axis]=.5*(
					(ale_displacement[edge[0]][axis]-previous_ale[edge[0]][axis])
					+(ale_displacement[edge[1]][axis]-previous_ale[edge[1]][axis]))
					/scenario.dt;
			}
			for(int axis=0;axis<3;++axis){
				constrained[3*node+axis]=1;prescribed[3*node+axis]=velocity[axis];
			}
		}
		std::map<int,double> pressure;
		for(const auto label:scenario.outlet_labels)pressure.emplace(label,0.);
		const auto boundary_load=BoundaryLoad(current,topology,pressure,dofs);
		for(int cell=0;cell<cells;++cell)
			for(int local=0;local<4;++local)
				for(int axis=0;axis<3;++axis)
					mesh_velocity[12*cell+3*local+axis]=
						(ale_displacement[reference.cells[cell].nodes[local]][axis]
						-previous_ale[reference.cells[cell].nodes[local]][axis])/scenario.dt;
		device_quadrature.CopyFromHost(quadrature.data(),quadrature.size());
		device_constrained.CopyFromHost(constrained.data(),constrained.size());
		device_prescribed.CopyFromHost(prescribed.data(),prescribed.size());
		device_boundary_load.CopyFromHost(boundary_load.data(),boundary_load.size());
		device_mesh_velocity.CopyFromHost(mesh_velocity.data(),mesh_velocity.size());
		iga::cuda::Check(cudaMemcpy(device_previous.data(),device_state.data(),
			static_cast<std::size_t>(dofs)*sizeof(double),cudaMemcpyDeviceToDevice),
			"copy previous FEM CUDA FSI state");
		std::vector<double> increment(dofs);
		bool converged=false;
		for(int iteration=1;iteration<=30;++iteration){
			const auto assembly_start=std::chrono::steady_clock::now();
			device_residual.Clear();device_jacobian.Clear();
			const int field_count=static_cast<int>(quadrature.size());
			const int entry_count=cells*iga::cuda::kTetFlowEntries;
			iga::cuda::EvaluateNativeTetFlowFields<<<(field_count+127)/128,128>>>(
				device_quadrature.data(),device_rows.data(),device_state.data(),
				device_previous.data(),device_mesh_velocity.data(),cells,fields.data());
			iga::cuda::AssembleNativeTetFlowElements<<<(entry_count+127)/128,128>>>(
				device_quadrature.data(),fields.data(),device_acceleration.data(),
				scenario.density,scenario.viscosity,scenario.dt,cells,
				device_element_residual.data(),device_element_jacobian.data());
			iga::cuda::ScatterNativeTetFlowElements<<<(entry_count+127)/128,128>>>(
				device_rows.data(),device_element_residual.data(),
				device_element_jacobian.data(),cells,dofs,
				device_residual.data(),device_jacobian.data());
			iga::cuda::ApplyNativeTetFlowBoundaryRows<<<
				(static_cast<std::size_t>(dofs)*dofs+127)/128,128>>>(
				device_constrained.data(),dofs,device_jacobian.data());
			iga::cuda::NativeTetFlowNewtonRightHandSide<<<(dofs+127)/128,128>>>(
				device_residual.data(),device_boundary_load.data(),device_state.data(),
				device_constrained.data(),device_prescribed.data(),dofs,device_rhs.data());
			iga::cuda::CheckKernel("FEM CUDA FSI assembly");
			iga::cuda::Check(cudaDeviceSynchronize(),"FEM CUDA FSI assembly sync");
			assembly_seconds+=std::chrono::duration<double>(
				std::chrono::steady_clock::now()-assembly_start).count();
			const auto solve_start=std::chrono::steady_clock::now();
			solver.Solve(dofs,device_jacobian.data(),device_rhs.data());
			device_rhs.CopyToHost(increment.data(),increment.size());
			solve_seconds+=std::chrono::duration<double>(
				std::chrono::steady_clock::now()-solve_start).count();
			iga::cuda::UpdateNativeTetFlowState<<<(dofs+127)/128,128>>>(
				device_state.data(),device_rhs.data(),dofs);
			const double maximum=*std::max_element(increment.begin(),increment.end(),
				[](double a,double b){return std::abs(a)<std::abs(b);});
			if(std::abs(maximum)<1e-8){converged=true;break;}
		}
		if(!converged)throw std::runtime_error("FEM CUDA FSI flow did not converge");
		device_state.CopyToHost(state.data(),state.size());
		const auto tractions=iga::EvaluateNativeTetFluidTractionOnMatchingStructure(
			current,topology,state,scenario.viscosity,scenario.interface_label);
		const auto load=InterfaceForces(current,solid,fluid_to_solid,
			scenario.interface_label,tractions);
		const auto old_solid=solid_displacement;
		const auto solid_result=iga::SolveNativeTetSolidStatic(solid,
			{scenario.young,scenario.poisson,scenario.solid_density},load,
			solid_constraints,solid_displacement,solid_options);
		solid_displacement=solid_result.displacement_m;
		std::map<std::uint64_t,std::array<double,3>> traction_by_face;
		for(const auto& traction:tractions)
			traction_by_face.emplace(traction.triangle_id,
				traction.traction_on_structure_pa);
		std::array<double,3> interface_force{};
		double interface_power=0.;
		for(const auto& face:current.boundary_triangles){
			if(face.boundary_label!=scenario.interface_label)continue;
			const auto& a=current.points[face.nodes[0]];
			const auto& b=current.points[face.nodes[1]];
			const auto& c=current.points[face.nodes[2]];
			const std::array<double,3> ab{{b[0]-a[0],b[1]-a[1],b[2]-a[2]}};
			const std::array<double,3> ac{{c[0]-a[0],c[1]-a[1],c[2]-a[2]}};
			const double area=.5*iga::NativeTetInterfaceNorm(
				iga::NativeTetInterfaceCross(ab,ac));
			for(int axis=0;axis<3;++axis){
				const double force=area*traction_by_face.at(face.id)[axis];
				interface_force[axis]+=force;
				double velocity=0.;
				for(const auto node:face.nodes){
					const auto solid_node=fluid_to_solid.at(node);
					velocity+=(solid_displacement[3*solid_node+axis]
						-old_solid[3*solid_node+axis])/(3.*scenario.dt);
				}
				interface_power+=force*velocity;
			}
		}
		std::map<std::uint32_t,std::array<double,3>> ale_boundary;
		for(const auto& boundary:reference.boundary_triangles)
			for(const auto node:boundary.nodes)ale_boundary.try_emplace(node,
				std::array<double,3>{{0.,0.,0.}});
		for(const auto& item:fluid_to_solid)
			for(int axis=0;axis<3;++axis)
				ale_boundary[item.first][axis]=solid_displacement[3*item.second+axis];
		previous_ale=ale_displacement;
		ale_displacement=iga::SolveNativeTetAleHarmonicMotion(reference,
			ale_boundary).displacement_m;
		const auto next=Deformed(reference,ale_displacement);
		std::vector<std::array<double,3>> next_mesh_velocity(reference.points.size());
		for(std::size_t node=0;node<reference.points.size();++node)
			for(int axis=0;axis<3;++axis)
				next_mesh_velocity[node][axis]=(next.points[node][axis]
					-current.points[node][axis])/scenario.dt;
		std::vector<std::array<double,3>> vertex_velocity(reference.points.size());
		for(std::size_t node=0;node<reference.points.size();++node)
			for(int axis=0;axis<3;++axis)vertex_velocity[node][axis]=state[3*node+axis];
		const auto conservation=iga::EvaluateNativeTetAleConservation(current,next,
			next_mesh_velocity,vertex_velocity,vertex_velocity,scenario.dt);
		iga::NativeTetLinearElasticResult solid_field;
		solid_field.displacement_m.resize(solid.points.size());
		for(std::size_t node=0;node<solid.points.size();++node)
			for(int axis=0;axis<3;++axis)
				solid_field.displacement_m[node][axis]=solid_displacement[3*node+axis];
		const auto flow_path=output/("flow_step_"+std::to_string(step)+".vtu");
		const auto wall_path=output/("wall_step_"+std::to_string(step)+".vtu");
		iga::WriteVtuPartition(flow_path,iga::BuildNativeTetHydraulicVtkPartition(
			reference,current,state,scenario.viscosity,0,1),step*scenario.dt);
		iga::WriteVtuPartition(wall_path,iga::BuildNativeTetLinearElasticVtkPartition(
			solid,solid_field,0,1),step*scenario.dt);
		flow_series.push_back({step*scenario.dt,flow_path});
		wall_series.push_back({step*scenario.dt,wall_path});
		iga::WritePvd(output/"flow.pvd",flow_series);
		iga::WritePvd(output/"wall.pvd",wall_series);
		double maximum_displacement=0.;
		for(const auto value:solid_displacement)
			maximum_displacement=std::max(maximum_displacement,std::abs(value));
		std::cout<<std::setprecision(17)<<"native_tet_fsi_cuda_step step="<<step
			<<" gcl_residual_m3_s="<<conservation.gcl_residual_m3_s
			<<" moving_balance_residual_m3_s="
				<<conservation.moving_domain_balance_residual_m3_s
			<<" wall_max_displacement_m="<<maximum_displacement
			<<" interface_force_n="<<iga::NativeTetInterfaceNorm(interface_force)
			<<" interface_power_w="<<interface_power
			<<" solid_residual_n="<<solid_result.final_free_residual_n<<'\n';
	}
	const auto usage=[](){struct rusage value{};getrusage(RUSAGE_SELF,&value);
		return value.ru_maxrss;}();
	std::cout<<"native_tet_fsi_cuda_complete steps="<<scenario.steps
		<<" assembly_s="<<assembly_seconds<<" solve_s="<<solve_seconds
		<<" host_peak_rss_kib="<<usage
		<<" cuda_peak_bytes="<<iga::cuda::DeviceAllocationCounter::Peak()<<'\n';
}

} // namespace

int main(int argc,char** argv)
{
	try{
		if(argc!=3)throw std::invalid_argument(
			"usage: native_tet_fsi_cuda CASE.json OUTPUT_DIR");
		RunFsi(argv[1],argv[2]);
		return 0;
	}catch(const std::exception& error){
		std::cerr<<"native_tet_fsi_cuda: "<<error.what()<<'\n';
		return 1;
	}
}
