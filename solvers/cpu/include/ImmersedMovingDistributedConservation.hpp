#ifndef IGA_IMMERSED_MOVING_DISTRIBUTED_CONSERVATION_HPP
#define IGA_IMMERSED_MOVING_DISTRIBUTED_CONSERVATION_HPP

#include "ImmersedTransientFlowRuntime.hpp"

namespace iga {
struct ImmersedMovingDistributedConservation {
	ImmersedTransientFlowConservationDiagnostics endpoint;
	std::string source_geometry_identity_sha256,target_geometry_identity_sha256;
	std::string source_publication_identity_sha256,target_publication_identity_sha256;
	double source_time_s=0,target_time_s=0,dt_s=0;
	std::uint64_t source_index=0,target_index=0;
	double source_audited_volume_m3=0,target_audited_volume_m3=0;
	double backward_euler_volume_rate_m3_s=0,reynolds_defect_m3_s=0,moving_mass_defect_m3_s=0;
	double normalization_scale_m3_s=1;
	double normalized_divergence_theorem_defect=0,normalized_reynolds_defect=0;
	double normalized_moving_mass_defect=0,normalized_wall_relative_leakage=0;
};

} // namespace iga
#endif
