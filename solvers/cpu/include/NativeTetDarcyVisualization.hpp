#ifndef IGA_NATIVE_TET_DARCY_VISUALIZATION_HPP
#define IGA_NATIVE_TET_DARCY_VISUALIZATION_HPP

#include "NativeTetDarcyPetsc.hpp"
#include "PartitionedVtkOutput.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace iga {

inline VtkPartition BuildNativeTetDarcyVtkPartition(const NativeTetMesh& mesh,
	const NativeTetDarcyResult& solution,const std::vector<double>& mobility,
	const std::vector<double>& source,int rank,int ranks)
{
	if(rank<0||ranks<1||rank>=ranks
		||solution.pressure_pa.size()!=mesh.points.size()
		||solution.cell_flux_m_s.size()!=mesh.cells.size()
		||solution.conservative_face_flow_m3_s.size()!=mesh.cells.size()
		||mobility.size()!=mesh.cells.size()||source.size()!=mesh.cells.size())
		throw std::invalid_argument("native Darcy visualization shape or rank is invalid");
	VtkPartition piece;
	VtkPointArray pressure{"pressure_pa",1,{}};
	VtkPointArray flux{"darcy_flux_m_s",3,{}};
	VtkPointArray recovered_flux{"darcy_rt0_centroid_flux_m_s",3,{}};
	VtkPointArray face_flow{"conservative_outward_face_flow_m3_s",4,{}};
	VtkPointArray cell_mobility{"mobility_m2_pa_s",1,{}};
	VtkPointArray cell_source{"source_s_inv",1,{}};
	std::map<std::uint32_t,std::int64_t> local_point;
	for(std::size_t index=static_cast<std::size_t>(rank);index<mesh.cells.size();
		index+=static_cast<std::size_t>(ranks)){
		const auto& cell=mesh.cells[index];
		EvaluateNativeTetGeometry(mesh,cell);
		for(const auto global:cell.nodes){
			const auto inserted=local_point.emplace(global,
				static_cast<std::int64_t>(local_point.size()));
			if(inserted.second){
				piece.point_ids.push_back(global);
				pressure.values.push_back(solution.pressure_pa.at(global));
				for(int axis=0;axis<3;++axis)
					piece.grid.points.push_back(mesh.points.at(global)[axis]);
			}
			piece.grid.connectivity.push_back(inserted.first->second);
		}
		piece.grid.offsets.push_back(
			static_cast<std::int64_t>(piece.grid.connectivity.size()));
		piece.grid.types.push_back(10);
		if(cell.id>static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()))
			throw std::invalid_argument("native Darcy visualization cell id exceeds range");
		piece.cell_ids.push_back(static_cast<std::int64_t>(cell.id));
		for(int axis=0;axis<3;++axis)
			flux.values.push_back(solution.cell_flux_m_s[index][axis]);
		std::array<double,3> centroid{};
		for(const auto node:cell.nodes)
			for(int axis=0;axis<3;++axis)
				centroid[axis]+=mesh.points.at(node)[axis]/4.;
		const auto recovered=EvaluateNativeTetDarcyRt0Flux(mesh,solution,index,centroid);
		for(int axis=0;axis<3;++axis)
			recovered_flux.values.push_back(recovered[axis]);
		for(int local=0;local<4;++local)
			face_flow.values.push_back(
				solution.conservative_face_flow_m3_s[index][local]);
		cell_mobility.values.push_back(mobility[index]);
		cell_source.values.push_back(source[index]);
	}
	piece.point_arrays.push_back(std::move(pressure));
	piece.cell_arrays.push_back(std::move(flux));
	piece.cell_arrays.push_back(std::move(recovered_flux));
	piece.cell_arrays.push_back(std::move(face_flow));
	piece.cell_arrays.push_back(std::move(cell_mobility));
	piece.cell_arrays.push_back(std::move(cell_source));
	ValidateVtkPartition(piece,0.);
	return piece;
}

} // namespace iga

#endif
