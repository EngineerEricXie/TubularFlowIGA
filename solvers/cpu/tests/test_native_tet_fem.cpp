#include "NativeTetFem.hpp"
#include "NativeTetBoundaryFlow.hpp"
#include "NativeTetHydraulicVisualization.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

bool Close(double first, double second)
{
	return std::abs(first-second) <= 1.0e-13*std::max({1.0, std::abs(first), std::abs(second)});
}

std::string Mesh(bool omit_last_boundary = false)
{
	return std::string(
		"$MeshFormat\n4.1 0 8\n$EndMeshFormat\n"
		"$PhysicalNames\n5\n"
		"2 1 \"boundary_label_0\"\n2 2 \"boundary_label_1\"\n"
		"2 3 \"boundary_label_2\"\n2 4 \"boundary_label_0\"\n3 1 \"fluid\"\n"
		"$EndPhysicalNames\n"
		"$Entities\n0 0 4 1\n"
		"1 0 0 0 1 1 1 1 1 0\n2 0 0 0 1 1 1 1 2 0\n"
		"3 0 0 0 1 1 1 1 3 0\n4 0 0 0 1 1 1 1 4 0\n"
		"1 0 0 0 1 1 1 1 1 0\n$EndEntities\n"
		"$Nodes\n1 4 1 4\n3 1 0 4\n1\n2\n3\n4\n"
		"0 0 0\n1 0 0\n0 1 0\n0 0 1\n$EndNodes\n"
		"$Elements\n")
		+(omit_last_boundary ? "4 4 1 4\n" : "5 5 1 5\n")
		+"2 1 2 1\n1 1 3 2\n"
		"2 2 2 1\n2 1 2 4\n"
		"2 3 2 1\n3 1 4 3\n"
		+(omit_last_boundary ? std::string() : "2 4 2 1\n4 2 3 4\n")
		+"3 1 4 1\n5 1 2 3 4\n$EndElements\n";
}

template <class Function>
void Reject(Function&& function)
{
	bool rejected = false;
	try { function(); }
	catch (const std::runtime_error&) { rejected = true; }
	assert(rejected);
}

} // namespace

int main(int argc, char** argv)
{
	if (argc == 2) {
		std::ifstream input(argv[1]);
		if (!input) throw std::runtime_error("cannot open requested native tetrahedral mesh");
		const auto mesh = iga::ReadNativeTetMeshGmsh41(input);
		const auto topology = iga::BuildNativeTaylorHoodTopology(mesh);
		std::cout << "native tetrahedral FEM external mesh passed nodes=" << mesh.points.size()
			<< " edges=" << topology.edges.size() << " cells=" << mesh.cells.size() << '\n';
		return 0;
	}
	if (argc != 1) throw std::runtime_error("usage: native_tet_fem_test [mesh.msh]");
	std::istringstream input(Mesh());
	const auto mesh = iga::ReadNativeTetMeshGmsh41(input);
	assert(mesh.points.size() == 4);
	assert(mesh.cells.size() == 1);
	assert(mesh.boundary_triangles.size() == 4);
	assert(mesh.boundary_triangles[0].boundary_label == 0);
	assert(mesh.boundary_triangles[1].boundary_label == 1);
	assert(mesh.boundary_triangles[2].boundary_label == 2);

	const auto topology = iga::BuildNativeTaylorHoodTopology(mesh);
	const auto& cell = mesh.cells.front();
	assert(topology.edges.size() == 6);
	assert(topology.cell_velocity_nodes.size() == 1);
	for (std::uint32_t i = 0; i < 10; ++i)
		assert(topology.cell_velocity_nodes[0][i] == i);
	assert(topology.boundary_pressure_nodes.at(1).size() == 3);
	assert(topology.boundary_velocity_nodes.at(1).size() == 6);
	assert(topology.boundary_velocity_nodes.at(0).size() == 9);
	{
		const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();
		std::vector<double> state(3*velocity_nodes+mesh.points.size(),0.);
		for(std::size_t node=0;node<velocity_nodes;++node){
			std::array<double,3> point{};
			if(node<mesh.points.size())point=mesh.points[node];
			else{
				const auto& edge=topology.edges[node-mesh.points.size()];
				for(int axis=0;axis<3;++axis)
					point[axis]=0.5*(mesh.points[edge[0]][axis]
						+mesh.points[edge[1]][axis]);
			}
			state[3*node]=2.*point[0];
			state[3*node+1]=3.*point[1];
			state[3*node+2]=-point[2];
		}
		for(std::size_t node=0;node<mesh.points.size();++node)
			state[3*velocity_nodes+node]=5.;
		const auto piece=iga::BuildNativeTetHydraulicVtkPartition(
			mesh,mesh,state,2.,0,1);
		assert(piece.cell_arrays.size()==1);
		assert(piece.cell_arrays[0].name=="cauchy_stress_pa");
		assert(piece.cell_arrays[0].components==9);
		const std::array<double,9> expected{{3.,0.,0.,0.,7.,0.,0.,0.,-9.}};
		for(std::size_t entry=0;entry<expected.size();++entry)
			assert(Close(piece.cell_arrays[0].values[entry],expected[entry]));
	}

	for (const auto coordinate : {std::array<double,3>{{0.1,0.2,0.3}},
		std::array<double,3>{{0.25,0.25,0.25}}, std::array<double,3>{{0,0,0}}}) {
		const auto basis = iga::EvaluateNativeTaylorHoodBasis(
			coordinate[0], coordinate[1], coordinate[2]);
		double velocity_sum = 0.0, pressure_sum = 0.0;
		std::array<double,3> gradient_sum{{0,0,0}};
		for (std::size_t i = 0; i < basis.velocity.size(); ++i) {
			velocity_sum += basis.velocity[i];
			for (int component = 0; component < 3; ++component)
				gradient_sum[component] += basis.velocity_gradient_reference[i][component];
		}
		for (const auto value : basis.pressure) pressure_sum += value;
		assert(Close(velocity_sum, 1.0));
		assert(Close(pressure_sum, 1.0));
		for (const auto value : gradient_sum) assert(Close(value, 0.0));
	}
	const auto vertex = iga::EvaluateNativeTaylorHoodBasis(0,0,0);
	assert(Close(vertex.velocity[0], 1.0));
	for (std::size_t i = 1; i < vertex.velocity.size(); ++i)
		assert(Close(vertex.velocity[i], 0.0));
	const auto geometry = iga::EvaluateNativeTetGeometry(mesh, mesh.cells.front());
	assert(Close(geometry.determinant, 1.0));
	const std::array<std::array<double,3>,4> expected_gradients{{{{-1,-1,-1}},
		{{1,0,0}},{{0,1,0}},{{0,0,1}}}};
	for (std::size_t i = 0; i < 4; ++i)
		for (int component = 0; component < 3; ++component)
			assert(Close(geometry.barycentric_gradients[i][component],
				expected_gradients[i][component]));
	double reference_volume = 0.0, linear_moment = 0.0;
	for (const auto& point : iga::NativeTetDegreeFiveQuadrature()) {
		reference_volume += point.weight;
		linear_moment += point.weight*point.reference[0];
	}
	assert(Close(reference_volume, 1.0/6.0));
	assert(Close(linear_moment, 1.0/24.0));

	using System = iga::NativeTaylorHoodElementSystem;
	std::array<double, System::dofs> zero{};
	const iga::NativeNavierStokesParameters parameters{1.7, 0.13};
	const auto zero_system = iga::BuildNativeTaylorHoodNavierStokesElement(
		mesh, mesh.cells.front(), zero, parameters);
	for (const auto value : zero_system.residual) assert(Close(value, 0.0));
	auto forced_parameters=parameters;
	forced_parameters.body_acceleration_m_s2[0]={{1.0,2.0,3.0,4.0}};
	const auto forced_system=iga::BuildNativeTaylorHoodNavierStokesElement(
		mesh,mesh.cells.front(),zero,forced_parameters);
	double integrated_x_force=0.0;
	for(std::size_t node=0;node<10;++node)integrated_x_force+=forced_system.residual[3*node];
	assert(Close(integrated_x_force,-parameters.density*3.25/6.0));
	for(std::size_t node=0;node<10;++node){
		assert(Close(forced_system.residual[3*node+1],0.0));
		assert(Close(forced_system.residual[3*node+2],0.0));
	}
	Reject([&] {
		auto invalid=parameters;
		invalid.body_acceleration_m_s2[2][3]=std::numeric_limits<double>::quiet_NaN();
		(void)iga::BuildNativeTaylorHoodNavierStokesElement(mesh,cell,zero,invalid);
	});
	std::array<double, System::dofs> state{};
	for (std::size_t i = 0; i < state.size(); ++i)
		state[i] = 0.01*std::sin(static_cast<double>(i+1));
	const auto system = iga::BuildNativeTaylorHoodNavierStokesElement(
		mesh, mesh.cells.front(), state, parameters);
	const double epsilon = 1.0e-7;
	double difference_squared = 0.0, reference_squared = 0.0;
	for (std::size_t column = 0; column < state.size(); ++column) {
		auto plus = state, minus = state;
		plus[column] += epsilon;
		minus[column] -= epsilon;
		const auto plus_system = iga::BuildNativeTaylorHoodNavierStokesElement(
			mesh, mesh.cells.front(), plus, parameters);
		const auto minus_system = iga::BuildNativeTaylorHoodNavierStokesElement(
			mesh, mesh.cells.front(), minus, parameters);
		for (std::size_t row = 0; row < state.size(); ++row) {
			const double finite_difference =
				(plus_system.residual[row]-minus_system.residual[row])/(2.0*epsilon);
			const double analytic = system.jacobian[row*System::dofs+column];
			difference_squared += (analytic-finite_difference)*(analytic-finite_difference);
			reference_squared += finite_difference*finite_difference;
		}
	}
	std::array<double, System::dofs> affine_state{};
	const std::array<std::array<int,2>,6> local_edges{{{{0,1}},{{0,2}},{{0,3}},
		{{1,2}},{{1,3}},{{2,3}}}};
	for (std::size_t node = 0; node < 4; ++node)
		affine_state[3*node] = mesh.points[cell.nodes[node]][0];
	for (std::size_t edge = 0; edge < local_edges.size(); ++edge)
		affine_state[3*(4+edge)] = 0.5*(mesh.points[cell.nodes[local_edges[edge][0]]][0]
			+mesh.points[cell.nodes[local_edges[edge][1]]][0]);
	std::array<std::array<double,3>,4> translating_grid{};
	for (auto& velocity : translating_grid) velocity = {{0.4,0,0}};
	auto shifted_state = affine_state;
	for (std::size_t node = 0; node < 10; ++node) shifted_state[3*node] -= 0.4;
	const auto ale_system = iga::BuildNativeTaylorHoodAleNavierStokesElement(mesh,cell,
		affine_state,translating_grid,{1000.0,0.004});
	const auto shifted_fixed_system = iga::BuildNativeTaylorHoodNavierStokesElement(mesh,cell,
		shifted_state,{1000.0,0.004});
	for (std::size_t row = 0; row < ale_system.residual.size(); ++row)
		assert(std::abs(ale_system.residual[row]-shifted_fixed_system.residual[row])<1e-11);
	for (std::size_t entry = 0; entry < ale_system.jacobian.size(); ++entry)
		assert(std::abs(ale_system.jacobian[entry]-shifted_fixed_system.jacobian[entry])<1e-11);
	std::array<double,System::dofs> uniform_state{};
	std::array<std::array<double,3>,4> uniform_grid{};
	for (std::size_t node = 0; node < 10; ++node) {
		uniform_state[3*node]=0.3;uniform_state[3*node+1]=-0.2;
		uniform_state[3*node+2]=0.1;
	}
	for (auto& velocity : uniform_grid) velocity={{0.3,-0.2,0.1}};
	const auto uniform_system=iga::BuildNativeTaylorHoodAleNavierStokesElement(
		mesh,cell,uniform_state,uniform_grid,{1000.0,0.004});
	for (const auto value : uniform_system.residual) assert(std::abs(value)<1e-12);
	std::vector<double> uniform_global_state(3*10+4,0.0);
	for(std::size_t node=0;node<10;++node) {
		uniform_global_state[3*node]=3.0;
		uniform_global_state[3*node+1]=-4.0;
	}
	const auto velocity_qoi=iga::EvaluateNativeTetVolumeVelocity(
		mesh,topology,uniform_global_state,2.0);
	assert(Close(velocity_qoi.volume_m3,1.0/6.0));
	assert(Close(velocity_qoi.mean_speed_m_s,5.0));
	assert(Close(velocity_qoi.rms_speed_m_s,5.0));
	assert(Close(velocity_qoi.kinetic_energy_j,25.0/6.0));
	std::vector<double> pressure_state(3*10+4,0.0);
	for(std::size_t node=0;node<4;++node) pressure_state[3*10+node]=2.0;
	const auto pressure_tractions=iga::EvaluateNativeTetBoundaryTractions(
		mesh,topology,pressure_state,0.5);
	for(const auto& labelled:pressure_tractions) {
		assert(Close(labelled.second.mean_normal_traction_pa,-2.0));
		assert(Close(labelled.second.rms_normal_traction_pa,2.0));
		assert(Close(labelled.second.mean_total_traction_pa,2.0));
		assert(Close(labelled.second.rms_total_traction_pa,2.0));
		assert(Close(labelled.second.mean_tangential_traction_pa,0.0));
		assert(Close(labelled.second.rms_tangential_traction_pa,0.0));
		for(const double value:labelled.second.integrated_tangential_traction_n)
			assert(Close(value,0.0));
	}
	const auto filtered_tractions=iga::EvaluateNativeTetBoundaryTractions(
		mesh,topology,pressure_state,0.5,[](int label,const std::array<double,3>&) {
			return label==1;
		});
	assert(filtered_tractions.size()==1&&filtered_tractions.count(1)==1);
	assert(Close(filtered_tractions.at(1).mean_total_traction_pa,2.0));
	std::vector<double> affine_global_state(3*10+4,0.0);
	for(std::size_t node=0;node<10;++node) {
		std::array<double,3> coordinate{};
		if(node<mesh.points.size()) coordinate=mesh.points[node];
		else {
			const auto edge=topology.edges[node-mesh.points.size()];
			for(int component=0;component<3;++component)
				coordinate[component]=0.5*(mesh.points[edge[0]][component]
					+mesh.points[edge[1]][component]);
		}
		affine_global_state[3*node]=coordinate[0]+2.0*coordinate[1]+3.0*coordinate[2];
		affine_global_state[3*node+1]=-coordinate[0]+coordinate[2];
		affine_global_state[3*node+2]=2.0*coordinate[1];
	}
	const auto recovered=iga::RecoverNativeTetVelocityGradients(
		mesh,topology,affine_global_state);
	const std::array<iga::NativeTetVelocityGradient,1> expected{{{{
		{{1,2,3}},{{-1,0,1}},{{0,2,0}}}}}};
	for(const auto& gradient:recovered)
		for(int component=0;component<3;++component)
			for(int derivative=0;derivative<3;++derivative)
				assert(Close(gradient[component][derivative],
					expected[0][component][derivative]));
	assert(std::sqrt(difference_squared/reference_squared) < 2.0e-9);

	Reject([] {
		std::istringstream incomplete(Mesh(true));
		(void)iga::ReadNativeTetMeshGmsh41(incomplete);
	});
	Reject([] {
		std::string binary = Mesh();
		binary.replace(binary.find("4.1 0 8"), 7, "4.1 1 8");
		std::istringstream invalid(binary);
		(void)iga::ReadNativeTetMeshGmsh41(invalid);
	});
	std::cout << "native tetrahedral FEM mesh/topology/basis tests passed\n";
	return 0;
}
