#ifndef IGA_NATIVE_TET_SPECIES_VISUALIZATION_HPP
#define IGA_NATIVE_TET_SPECIES_VISUALIZATION_HPP

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

// Rank-owned linear tetrahedra expose the native P1 concentration directly.
inline VtkPartition BuildNativeTetSpeciesVtkPartition(const NativeTetMesh& reference,
	const NativeTetMesh& current,const std::vector<double>& concentration_mol_m3,
	int rank,int ranks)
{
	if(rank<0||ranks<1||rank>=ranks)
		throw std::invalid_argument("native species visualization rank is invalid");
	if(reference.points.size()!=current.points.size()
		||reference.cells.size()!=current.cells.size()
		||concentration_mol_m3.size()!=current.points.size())
		throw std::invalid_argument("native species visualization shape differs");
	for(std::size_t cell=0;cell<current.cells.size();++cell)
		if(reference.cells[cell].id!=current.cells[cell].id
			||reference.cells[cell].nodes!=current.cells[cell].nodes)
			throw std::invalid_argument("native species visualization topology differs");
	VtkPartition piece;
	VtkPointArray concentration{"concentration_mol_m3",1,{}};
	VtkPointArray reference_position{"reference_position_m",3,{}};
	VtkPointArray displacement{"displacement_m",3,{}};
	std::map<std::uint32_t,std::int64_t> local_point;
	for(std::size_t cell=static_cast<std::size_t>(rank);cell<current.cells.size();
		cell+=static_cast<std::size_t>(ranks)){
		EvaluateNativeTetGeometry(current,current.cells[cell]);
		for(const auto global:current.cells[cell].nodes){
			const auto inserted=local_point.emplace(global,
				static_cast<std::int64_t>(local_point.size()));
			if(inserted.second){
				piece.point_ids.push_back(global);
				const double scalar=concentration_mol_m3[global];
				if(!std::isfinite(scalar))
					throw std::invalid_argument("native species visualization concentration is nonfinite");
				concentration.values.push_back(scalar);
				for(int axis=0;axis<3;++axis){
					piece.grid.points.push_back(current.points[global][axis]);
					reference_position.values.push_back(reference.points[global][axis]);
					displacement.values.push_back(current.points[global][axis]
						-reference.points[global][axis]);
				}
			}
			piece.grid.connectivity.push_back(inserted.first->second);
		}
		piece.grid.offsets.push_back(
			static_cast<std::int64_t>(piece.grid.connectivity.size()));
		piece.grid.types.push_back(10);
		if(current.cells[cell].id>static_cast<std::uint64_t>(
			std::numeric_limits<std::int64_t>::max()))
			throw std::invalid_argument("native species visualization cell id exceeds range");
		piece.cell_ids.push_back(static_cast<std::int64_t>(current.cells[cell].id));
	}
	piece.point_arrays.push_back(std::move(concentration));
	piece.point_arrays.push_back(std::move(reference_position));
	piece.point_arrays.push_back(std::move(displacement));
	ValidateVtkPartition(piece,0.);
	return piece;
}

} // namespace iga

#endif
