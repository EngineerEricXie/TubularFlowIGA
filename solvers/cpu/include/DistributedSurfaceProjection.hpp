#ifndef IGA_DISTRIBUTED_SURFACE_PROJECTION_HPP
#define IGA_DISTRIBUTED_SURFACE_PROJECTION_HPP

#include "FluidSurfaceTraction.hpp"
#include "SurfaceOwnershipValidation.hpp"
#include "OwnedPointValues.hpp"

namespace iga {
struct SurfaceMassEntry {
	std::uint64_t row_node=0,column_node=0;
	double value_m2=0.;
};

// Bounded single-owner consistent projection. Only the chosen owner allocates
// the dense matrix and solves it. Callers certify quadrature authority and the
// matching force/mass epoch before invoking this algebraic operation.
inline std::vector<std::array<double,3>> ProjectDistributedSurfaceTraction(MPI_Comm comm,
	const SurfaceInterfaceRef& reference,const DistributedSurfaceLayout& layout,
	const std::vector<std::array<double,3>>& owned_force_n,
	const std::vector<SurfaceMassEntry>& local_mass,int projection_owner,
	std::size_t maximum_nodes=4096,PointIdentityLimits limits={})
{
	using namespace point_identity_detail;
	int rank=0,ranks=0;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&ranks);
	RequireCollectiveSameInt(comm,"surface projection owner",projection_owner);
	std::vector<std::string> buckets;
	std::vector<std::uint64_t> root_ids;
	std::vector<double> force_tuples;
	CollectiveLocalStage(comm,"surface projection input",[&] {
		ValidateDistributedSurfaceLayout(layout);
		if(projection_owner<0||projection_owner>=ranks||!maximum_nodes||layout.global_node_count>maximum_nodes
			||layout.global_node_count>std::numeric_limits<std::size_t>::max()/layout.global_node_count
			||owned_force_n.size()!=layout.owned_global_node_ids.size())
			throw std::invalid_argument("invalid bounded surface projection dimensions or owner");
		if(!limits.max_local_occurrences||local_mass.size()>limits.max_local_occurrences)
			throw std::invalid_argument("surface projection mass record limit");
		buckets.resize(ranks);std::size_t used=0;
		for(const auto& entry:local_mass) {
			if(!std::isfinite(entry.value_m2))throw std::invalid_argument("nonfinite surface mass entry");
			Charge(used,24,limits.max_wire_bytes);
			auto& packet=buckets[projection_owner];Put(packet,entry.row_node);Put(packet,entry.column_node);
			std::uint64_t bits=0;std::memcpy(&bits,&entry.value_m2,sizeof(bits));Put(packet,bits);
		}
		for(const auto& value:owned_force_n)for(double component:value)force_tuples.push_back(component);
		if(rank==projection_owner)for(const auto& node:layout.reference_positions)root_ids.push_back(node.global_node_id);
	});
	ValidateSurfacePublicationOwnership(reference,layout,comm);
	const auto incoming=Exchange(comm,buckets,limits.max_wire_bytes);
	const auto forces=FetchOwnedPointValues(comm,layout.owned_global_node_ids,force_tuples,root_ids,3,limits);
	std::vector<double> projected;
	CollectiveLocalStage(comm,"surface consistent projection solve",[&] {
		if(rank!=projection_owner)return;
		const auto nodes=root_ids.size();
		std::vector<double> mass(nodes*nodes,0.);
		const auto index=[&](std::uint64_t id) {
			const auto found=std::lower_bound(root_ids.begin(),root_ids.end(),id);
			if(found==root_ids.end()||*found!=id)throw std::invalid_argument("unknown surface mass node");
			return static_cast<std::size_t>(found-root_ids.begin());
		};
		for(int peer=0;peer<ranks;++peer) {
			auto packet=incoming.Peer(peer);
			while(!packet.empty()) {
				const auto row=Get(packet),column=Get(packet),bits=Get(packet);
				double value=0.;std::memcpy(&value,&bits,sizeof(bits));
				auto& target=mass[index(row)*nodes+index(column)];target+=value;
				if(!std::isfinite(target))throw std::overflow_error("surface mass accumulation overflow");
			}
		}
		for(std::size_t row=0;row<nodes;++row)for(std::size_t col=0;col<row;++col) {
			const auto a=mass[row*nodes+col],b=mass[col*nodes+row];
			if(std::abs(a-b)>128.*std::numeric_limits<double>::epsilon()*std::max(std::abs(a),std::abs(b)))
				throw std::runtime_error("asymmetric surface consistent mass");
		}
		projected=fluid_surface_traction_detail::SolveConsistentMass(std::move(mass),forces,nodes);
	});
	const auto tuples=FetchOwnedPointValues(comm,root_ids,projected,layout.owned_global_node_ids,3,limits);
	std::vector<std::array<double,3>> result;
	CollectiveLocalStage(comm,"surface projected owner values",[&] {
		result.resize(layout.owned_global_node_ids.size());
		for(std::size_t row=0;row<result.size();++row)for(int axis=0;axis<3;++axis)result[row][axis]=tuples[3*row+axis];
	});
	return result;
}
} // namespace iga
#endif
