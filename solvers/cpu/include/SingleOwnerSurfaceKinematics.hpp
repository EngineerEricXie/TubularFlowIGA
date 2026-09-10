#ifndef IGA_SINGLE_OWNER_SURFACE_KINEMATICS_HPP
#define IGA_SINGLE_OWNER_SURFACE_KINEMATICS_HPP

#include "SurfaceOwnershipValidation.hpp"
#include "OwnedPointValues.hpp"
#include "FsiDomainRuntime.hpp"

namespace iga {
// The numerical owner supplies a complete single-partition publication and an
// independently captured identity from its solver/lifecycle. Other ranks pass
// null source pointers. Result stamps use the destination publication layout.
inline SurfaceKinematics DistributeSingleOwnerSurfaceKinematics(MPI_Comm comm,
	const DistributedSurfaceLayout& partition,int owner,const DistributedSurfaceLayout* source_layout,
	const SurfaceKinematics* source,const std::string& expected_source_identity,
	const SurfaceInterfaceRef& expected_interface,const FsiTrialContext& expected_context,
	PointIdentityLimits limits={})
{
	int rank=0,ranks=0;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&ranks);
	RequireCollectiveSameInt(comm,"surface kinematics numerical owner",owner);
	std::string common,producer(64,'0');
	std::vector<std::uint64_t> owned;
	std::vector<double> tuples;
	CollectiveLocalStage(comm,"single owner kinematics binding",[&] {
		ValidateDistributedSurfaceLayout(partition);ValidateSurfaceInterfaceRef(expected_interface);ValidateFsiTrialContext(expected_context);
		if(owner<0||owner>=ranks)throw std::invalid_argument("invalid kinematics numerical owner");
		Sha256 hash;
		for(const auto* value:{&expected_interface.domain_id,&expected_interface.subsystem_id,&expected_interface.interface_id})
			distributed_surface_detail::AppendString(hash,*value);
		hash.AppendLittleEndian64(expected_context.step);hash.AppendNormalizedDouble(expected_context.start_time_s);
		hash.AppendNormalizedDouble(expected_context.dt_s);hash.AppendLittleEndian64(expected_context.coupling_iteration);common=hash.Hex();
		if(rank!=owner) {
			if(source||source_layout||!expected_source_identity.empty())throw std::invalid_argument("kinematics source supplied by nonowner");
			return;
		}
		if(!source||!source_layout||!IsLowercaseSha256(expected_source_identity))throw std::invalid_argument("missing numerical kinematics authority");
		ValidateSurfaceKinematics(*source,*source_layout);
		if(source_layout->partition_count!=1||source_layout->partition_rank!=0
			||!(source->interface==expected_interface)||source->stamp.time_s!=expected_context.EndTime()
			||source->stamp.step!=expected_context.step||source->stamp.coupling_iteration!=expected_context.coupling_iteration
			||BuildSurfaceKinematicsIdentitySha256(*source,*source_layout)!=expected_source_identity)
			throw std::runtime_error("numerical kinematics differs from expected source or context");
		auto reference=*source_layout;reference.partition_count=partition.partition_count;
		if(BuildDistributedSurfaceLayoutIdentitySha256(reference)!=partition.layout_identity_sha256)
			throw std::runtime_error("kinematics source and destination reference geometry differ");
		owned=source_layout->owned_global_node_ids;producer=source->stamp.producer_state_identity_sha256;
		if(owned.size()>tuples.max_size()/6)throw std::overflow_error("numerical kinematics tuple overflow");
		tuples.reserve(6*owned.size());
		for(std::size_t row=0;row<owned.size();++row) {
			for(double value:source->displacement_m[row])tuples.push_back(value);
			for(double value:source->velocity_m_per_s[row])tuples.push_back(value);
		}
	});
	RequireCollectiveSameText(comm,"surface kinematics destination context",common);
	ValidateSurfacePublicationOwnership(expected_interface,partition,comm);
	MPI_Bcast(producer.data(),64,MPI_CHAR,owner,comm);
	const auto values=FetchOwnedPointValues(comm,owned,tuples,partition.owned_global_node_ids,6,limits);
	SurfaceKinematics result;
	CollectiveLocalStage(comm,"surface kinematics destination publication",[&] {
		result.interface=expected_interface;result.stamp.time_s=expected_context.EndTime();result.stamp.step=expected_context.step;
		result.stamp.coupling_iteration=expected_context.coupling_iteration;
		result.stamp.reference_mesh_identity_sha256=partition.reference_mesh_identity_sha256;
		result.stamp.layout_identity_sha256=partition.layout_identity_sha256;
		result.stamp.partition_identity_sha256=BuildDistributedSurfacePartitionIdentitySha256(partition);
		result.stamp.producer_state_identity_sha256=producer;
		result.displacement_m.resize(partition.owned_global_node_ids.size());result.velocity_m_per_s.resize(result.displacement_m.size());
		for(std::size_t row=0;row<result.displacement_m.size();++row)for(int axis=0;axis<3;++axis) {
			result.displacement_m[row][axis]=values[6*row+axis];result.velocity_m_per_s[row][axis]=values[6*row+3+axis];
		}
		ValidateSurfaceKinematics(result,partition);
	});
	return result;
}
} // namespace iga
#endif
