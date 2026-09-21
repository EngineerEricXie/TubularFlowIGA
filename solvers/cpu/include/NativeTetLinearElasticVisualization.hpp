#ifndef IGA_NATIVE_TET_LINEAR_ELASTIC_VISUALIZATION_HPP
#define IGA_NATIVE_TET_LINEAR_ELASTIC_VISUALIZATION_HPP

#include "NativeTetLinearElasticPetsc.hpp"
#include "PartitionedVtkOutput.hpp"

#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>

namespace iga {

inline VtkPartition BuildNativeTetLinearElasticVtkPartition(
	const NativeTetMesh& mesh,const NativeTetLinearElasticResult& result,
	int rank,int ranks)
{
	if(rank<0||ranks<1||rank>=ranks||
		result.displacement_m.size()!=mesh.points.size())
		throw std::invalid_argument("solid VTK mesh, field or rank invalid");
	VtkPartition piece;
	VtkPointArray displacement{"displacement_m",3,{}};
	VtkPointArray reference{"reference_position_m",3,{}};
	std::map<std::uint32_t,std::int64_t> local_point;
	for(std::size_t index=static_cast<std::size_t>(rank);index<mesh.cells.size();
		index+=static_cast<std::size_t>(ranks)){
		const auto& cell=mesh.cells[index];
		for(auto node:cell.nodes){
			const auto inserted=local_point.emplace(node,
				static_cast<std::int64_t>(local_point.size()));
			if(inserted.second){
				piece.point_ids.push_back(node);
				for(int axis=0;axis<3;++axis){
					const double x=mesh.points.at(node)[axis];
					const double u=result.displacement_m[node][axis];
					piece.grid.points.push_back(x+u);
					displacement.values.push_back(u);
					reference.values.push_back(x);
				}
			}
			piece.grid.connectivity.push_back(inserted.first->second);
		}
		piece.grid.offsets.push_back(
			static_cast<std::int64_t>(piece.grid.connectivity.size()));
		piece.grid.types.push_back(10);
		if(cell.id>static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()))
			throw std::invalid_argument("solid cell id exceeds VTK range");
		piece.cell_ids.push_back(static_cast<std::int64_t>(cell.id));
	}
	piece.point_arrays.push_back(std::move(displacement));
	piece.point_arrays.push_back(std::move(reference));
	ValidateVtkPartition(piece,0.);
	return piece;
}

} // namespace iga

#endif
