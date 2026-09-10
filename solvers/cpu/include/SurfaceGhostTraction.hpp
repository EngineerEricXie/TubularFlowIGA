#ifndef IGA_SURFACE_GHOST_TRACTION_HPP
#define IGA_SURFACE_GHOST_TRACTION_HPP

#include "SurfaceOwnershipValidation.hpp"
#include "OwnedPointValues.hpp"

namespace iga {

// Requested-node data, deliberately distinct from an owned publication.
struct SurfaceGhostTraction {
	std::vector<std::uint64_t> requested_node_ids;
	std::vector<std::array<double,3>> traction_on_structure_pa,consistent_nodal_force_n;
};

// Every owner supplies its independently expected stamp. Producer-state IDs
// may be partition-local; interface, geometry, time, step and iteration must
// agree globally. The caller must obtain expected stamps from its transaction
// authority, not merely copy an untrusted incoming publication's stamp. The
// independently expected projection identity is also partition-local.
inline SurfaceGhostTraction FetchSurfaceGhostTraction(MPI_Comm comm,
	const DistributedSurfaceLayout& layout,const SurfaceTraction& publication,
	const SurfaceInterfaceRef& expected_interface,const SurfaceFieldStamp& expected_stamp,
	const std::string& expected_projection_identity,const std::vector<std::uint64_t>& requested,PointIdentityLimits limits={})
{
	std::string epoch;
	std::vector<double> tuples;
	CollectiveLocalStage(comm,"surface traction ghost publication binding",[&] {
		ValidateSurfaceInterfaceRef(expected_interface);
		ValidateSurfaceFieldStamp(expected_stamp,layout);
		ValidateSurfaceTraction(publication,layout);
		if(!IsLowercaseSha256(expected_projection_identity)||publication.projection_identity_sha256!=expected_projection_identity)
			throw std::runtime_error("surface traction ghost projection differs from expected identity");
		if(!(publication.interface==expected_interface)
			||BuildSurfaceFieldStampIdentitySha256(publication.stamp,layout)
				!=BuildSurfaceFieldStampIdentitySha256(expected_stamp,layout))
			throw std::runtime_error("surface traction ghost publication differs from expected owner stamp");
		if(!limits.max_local_occurrences||layout.owned_global_node_ids.size()>limits.max_local_occurrences
			||requested.size()>limits.max_local_occurrences-layout.owned_global_node_ids.size()
			||layout.owned_global_node_ids.size()>tuples.max_size()/6)
			throw std::invalid_argument("surface traction ghost tuples exceed local record limit");
		Sha256 hash;
		for(const auto* text:{&expected_interface.domain_id,&expected_interface.subsystem_id,&expected_interface.interface_id,
			&expected_stamp.reference_mesh_identity_sha256,&expected_stamp.layout_identity_sha256})
			distributed_surface_detail::AppendString(hash,*text);
		hash.AppendNormalizedDouble(expected_stamp.time_s);hash.AppendLittleEndian64(expected_stamp.step);
		hash.AppendLittleEndian64(expected_stamp.coupling_iteration);epoch=hash.Hex();
		tuples.reserve(6*layout.owned_global_node_ids.size());
		for(std::size_t row=0;row<layout.owned_global_node_ids.size();++row) {
			for(double value:publication.traction_on_structure_pa[row])tuples.push_back(value);
			for(double value:publication.consistent_nodal_force_n[row])tuples.push_back(value);
		}
	});
	RequireCollectiveSameText(comm,"surface traction ghost common epoch",epoch);
	ValidateSurfacePublicationOwnership(expected_interface,layout,comm);
	const auto fetched=FetchOwnedPointValues(comm,layout.owned_global_node_ids,tuples,requested,6,limits);
	SurfaceGhostTraction result;
	CollectiveLocalStage(comm,"surface traction ghost result",[&] {
		result.requested_node_ids=requested;
		result.traction_on_structure_pa.resize(requested.size());result.consistent_nodal_force_n.resize(requested.size());
		for(std::size_t row=0;row<requested.size();++row)for(int component=0;component<3;++component) {
			result.traction_on_structure_pa[row][component]=fetched[6*row+component];
			result.consistent_nodal_force_n[row][component]=fetched[6*row+3+component];
		}
	});
	return result;
}
} // namespace iga
#endif
