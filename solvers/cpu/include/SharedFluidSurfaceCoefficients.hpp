#ifndef IGA_SHARED_FLUID_SURFACE_COEFFICIENTS_HPP
#define IGA_SHARED_FLUID_SURFACE_COEFFICIENTS_HPP

#include "FluidSurfaceTraction.hpp"
#include "DistributedPointIdentity.hpp"

namespace iga {
// Validate exact shared coefficient copies without replicating the global
// field. Identical repeated nodes are valid; both local and remote conflicts
// fail collectively. This certifies agreement, not ownership or time freshness.
inline void ValidateSharedFluidSurfaceCoefficients(MPI_Comm comm,
	const CubicCartesianBackground& background,const std::vector<FluidSurfaceElementState>& state,
	PointIdentityLimits limits={})
{
	using namespace point_identity_detail;
	int ranks=0;MPI_Comm_size(comm,&ranks);
	std::vector<std::string> buckets;
	std::string geometry_identity;
	CollectiveLocalStage(comm,"shared fluid coefficient input",[&] {
		if(!limits.max_local_occurrences||state.size()>limits.max_local_occurrences)
			throw std::invalid_argument("fluid coefficient record limit");
		Sha256 hash;distributed_surface_detail::AppendString(hash,"SharedFluidSurfaceCoefficients/grid/v1");
		for(int axis=0;axis<3;++axis) {
			hash.AppendNormalizedDouble(background.Spec().lower_m[axis]);hash.AppendNormalizedDouble(background.Spec().upper_m[axis]);
			hash.AppendLittleEndian64(background.Spec().cells[axis]);
		}
		geometry_identity=hash.Hex();
		std::map<std::uint64_t,std::array<std::uint64_t,4>> coefficients;
		std::size_t records=0;
		for(const auto& item:state) {
			if(item.nodal_state.size()>limits.max_local_occurrences-records)throw std::invalid_argument("fluid coefficient occurrence limit");
			records+=item.nodal_state.size();
			const auto element=background.MaterializeElement(item.cell_id);
			if(item.nodal_state.size()!=element.connectivity.size())throw std::invalid_argument("fluid coefficient count mismatch");
			for(std::size_t node=0;node<element.connectivity.size();++node) {
				std::array<std::uint64_t,4> bits{};
				for(int component=0;component<4;++component) {
					const double value=item.nodal_state[node][component];
					if(!std::isfinite(value))throw std::invalid_argument("nonfinite shared fluid coefficient");
					std::memcpy(&bits[component],&value,sizeof(value));
				}
				const auto inserted=coefficients.emplace(element.connectivity[node],bits);
				if(!inserted.second&&inserted.first->second!=bits)throw std::invalid_argument("conflicting local fluid coefficient copies");
			}
		}
		buckets.resize(ranks);std::size_t used=0;
		for(const auto& item:coefficients) {
			Charge(used,40,limits.max_wire_bytes);auto& packet=buckets[item.first%ranks];Put(packet,item.first);
			for(auto value:item.second)Put(packet,value);
		}
	});
	RequireCollectiveSameText(comm,"shared fluid coefficient grid",geometry_identity);
	const auto incoming=Exchange(comm,buckets,limits.max_wire_bytes);
	CollectiveLocalStage(comm,"shared fluid coefficient agreement",[&] {
		std::map<std::uint64_t,std::array<std::uint64_t,4>> coefficients;
		for(int peer=0;peer<ranks;++peer) {
			auto packet=incoming.Peer(peer);
			while(!packet.empty()) {
				const auto id=Get(packet);std::array<std::uint64_t,4> bits{};
				for(auto& value:bits)value=Get(packet);
				const auto inserted=coefficients.emplace(id,bits);
				if(!inserted.second&&inserted.first->second!=bits)throw std::runtime_error("conflicting remote fluid coefficient copies");
			}
		}
	});
}
} // namespace iga
#endif
