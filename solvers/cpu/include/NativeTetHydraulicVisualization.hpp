#ifndef IGA_NATIVE_TET_HYDRAULIC_VISUALIZATION_HPP
#define IGA_NATIVE_TET_HYDRAULIC_VISUALIZATION_HPP

#include "NativeTetFem.hpp"
#include "PartitionedVtkOutput.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace iga {

// Native local P2 edges: 01,02,03,12,13,23. VTK quadratic tetra (type 24)
// expects 01,12,20,03,13,23; see vtkQuadraticTetra class documentation.
// https://vtk.org/doc/nightly/release/8.1/html/classvtkQuadraticTetra.html
inline VtkPartition BuildNativeTetHydraulicVtkPartition(const NativeTetMesh& reference,
	const NativeTetMesh& mesh,const std::vector<double>& state,
	double dynamic_viscosity_pa_s,int rank,int ranks)
{
	if(rank<0||ranks<1||rank>=ranks)
		throw std::invalid_argument("native tetra visualization rank is invalid");
	if(!std::isfinite(dynamic_viscosity_pa_s)||dynamic_viscosity_pa_s<=0.)
		throw std::invalid_argument("native tetra visualization viscosity is invalid");
	if(reference.points.size()!=mesh.points.size()
		||reference.cells.size()!=mesh.cells.size())
		throw std::invalid_argument("native tetra visualization reference topology differs");
	for(std::size_t cell=0;cell<mesh.cells.size();++cell)
		if(reference.cells[cell].id!=mesh.cells[cell].id
			||reference.cells[cell].nodes!=mesh.cells[cell].nodes)
			throw std::invalid_argument("native tetra visualization reference cell differs");
	const auto topology=BuildNativeTaylorHoodTopology(mesh);
	const auto vertex_nodes=mesh.points.size();
	const auto velocity_nodes=vertex_nodes+topology.edges.size();
	if(state.size()!=3*velocity_nodes+vertex_nodes)
		throw std::invalid_argument("native tetra visualization field shape differs");
	for(const double value:state)
		if(!std::isfinite(value))
			throw std::invalid_argument("native tetra visualization field is nonfinite");
	const std::array<std::size_t,10> vtk_local{{0,1,2,3,4,7,5,6,8,9}};
	std::map<std::uint32_t,std::int64_t> local_point;
	VtkPartition piece;
	VtkPointArray velocity{"velocity_m_s",3,{}};
	VtkPointArray pressure{"pressure_pa",1,{}};
	VtkPointArray reference_position{"reference_position_m",3,{}};
	VtkPointArray displacement{"displacement_m",3,{}};
	VtkPointArray stress{"cauchy_stress_pa",9,{}};
	for(std::size_t cell=static_cast<std::size_t>(rank);cell<mesh.cells.size();
		cell+=static_cast<std::size_t>(ranks)){
		const auto geometry=EvaluateNativeTetGeometry(mesh,mesh.cells[cell]);
		const auto basis=EvaluateNativeTaylorHoodPhysicalBasis(geometry,
			{{0.25,0.25,0.25}});
		std::array<std::array<double,3>,3> gradient{};
		for(std::size_t local=0;local<10;++local){
			const auto global=topology.cell_velocity_nodes[cell][local];
			for(int component=0;component<3;++component)
				for(int axis=0;axis<3;++axis)
					gradient[component][axis]+=state[3*global+component]
						*basis.velocity_gradients[local][axis];
		}
		double cell_pressure=0.;
		for(std::size_t local=0;local<4;++local)
			cell_pressure+=basis.pressure[local]
				*state[3*velocity_nodes+mesh.cells[cell].nodes[local]];
		for(int row=0;row<3;++row)
			for(int column=0;column<3;++column)
				stress.values.push_back(dynamic_viscosity_pa_s
					*(gradient[row][column]+gradient[column][row])
					-(row==column?cell_pressure:0.));
		for(const auto local:vtk_local){
			const auto global=topology.cell_velocity_nodes[cell][local];
			auto inserted=local_point.emplace(global,
				static_cast<std::int64_t>(local_point.size()));
			if(inserted.second){
				piece.point_ids.push_back(global);
				std::array<double,3> coordinate{};
				std::array<double,3> reference_coordinate{};
				double point_pressure=0.;
				if(global<vertex_nodes){
					coordinate=mesh.points[global];
					reference_coordinate=reference.points[global];
					point_pressure=state[3*velocity_nodes+global];
				}else{
					const auto& edge=topology.edges[global-vertex_nodes];
					for(int axis=0;axis<3;++axis){
						coordinate[axis]=0.5*(mesh.points[edge[0]][axis]
							+mesh.points[edge[1]][axis]);
						reference_coordinate[axis]=0.5*(reference.points[edge[0]][axis]
							+reference.points[edge[1]][axis]);
					}
					point_pressure=0.5*(state[3*velocity_nodes+edge[0]]
						+state[3*velocity_nodes+edge[1]]);
				}
				for(int axis=0;axis<3;++axis){
					piece.grid.points.push_back(coordinate[axis]);
					velocity.values.push_back(state[3*global+axis]);
					reference_position.values.push_back(reference_coordinate[axis]);
					displacement.values.push_back(coordinate[axis]
						-reference_coordinate[axis]);
				}
				pressure.values.push_back(point_pressure);
			}
			piece.grid.connectivity.push_back(inserted.first->second);
		}
		piece.grid.offsets.push_back(
			static_cast<std::int64_t>(piece.grid.connectivity.size()));
		piece.grid.types.push_back(24);
		if(mesh.cells[cell].id>static_cast<std::uint64_t>(
			std::numeric_limits<std::int64_t>::max()))
			throw std::invalid_argument("native tetra visualization cell id exceeds range");
		piece.cell_ids.push_back(static_cast<std::int64_t>(mesh.cells[cell].id));
	}
	piece.point_arrays.push_back(std::move(velocity));
	piece.point_arrays.push_back(std::move(pressure));
	piece.point_arrays.push_back(std::move(reference_position));
	piece.point_arrays.push_back(std::move(displacement));
	piece.cell_arrays.push_back(std::move(stress));
	ValidateVtkPartition(piece,0.);
	return piece;
}

} // namespace iga

#endif
