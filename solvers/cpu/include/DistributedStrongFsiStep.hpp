#ifndef IGA_DISTRIBUTED_STRONG_FSI_STEP_HPP
#define IGA_DISTRIBUTED_STRONG_FSI_STEP_HPP

#include "DistributedFsiCommitCoordinator.hpp"
#include "DistributedFsiConvergence.hpp"
#include "DistributedWeightedAitken.hpp"
#include "MaterialSurfacePatchMap.hpp"
#include "StrongFluidStructureCoupling.hpp"
#include <map>

namespace iga {

struct DistributedStrongFsiIteration {
	DistributedFsiConvergence displacement;
	double maximum_velocity_residual_times_dt_m=0.;
	double relaxation=0.;
};
struct DistributedStrongFsiResult {
	std::vector<DistributedStrongFsiIteration> history;
	bool converged=false;
};

// The supplied material is the independently expected accepted fluid epoch.
// Geometry remains replicated, but iterations, residuals and publications hold
// only owned surface values. Membrane normals use immutable reference facets.
template<class FluidRuntime,class StructureRuntime>
DistributedStrongFsiResult SolveDistributedStrongFsiStep(MPI_Comm comm,
	FluidRuntime& fluid,StructureRuntime& structure,const MaterialSurfacePatchMap& map,
	const MaterialSurfaceKinematics& committed,FsiTrialContext context,
	const StrongFluidStructureCouplingOptions& options={})
{
	const auto& layout=map.Layout();
	std::vector<std::array<double,3>> normals,accepted_displacement;
	std::vector<bool> clamped;
	std::vector<double> weights;
	SurfaceKinematics current;
	std::string identity,material_identity;
	DistributedStrongFsiResult result;
	CollectiveLocalStage(comm,"distributed strong FSI preflight",[&] {
		int rank=0,ranks=0;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&ranks);
		for(MPI_Comm runtime_comm:{fluid.Communicator(),structure.Communicator()}) {
			int comparison=MPI_UNEQUAL;MPI_Comm_compare(comm,runtime_comm,&comparison);
			if(comparison!=MPI_IDENT&&comparison!=MPI_CONGRUENT)
				throw std::invalid_argument("strong FSI runtime communicator mismatch");
		}
		if(layout.partition_rank!=static_cast<std::uint64_t>(rank)
			|| layout.partition_count!=static_cast<std::uint64_t>(ranks))
			throw std::invalid_argument("strong FSI surface partition communicator mismatch");
		ValidateFsiTrialContext(context);committed.Validate();
		if(context.coupling_iteration!=0 || !options.maximum_iterations
			|| committed.EvaluatedTimeS()!=context.start_time_s || committed.StepEndS()!=context.start_time_s)
			throw std::invalid_argument("strong FSI requires initial iteration and matching accepted clock");
		if(committed.MaterialIdentitySha256()!=map.FullReference().MaterialIdentitySha256()
			|| committed.TopologyIdentitySha256()!=map.FullReference().TopologyIdentitySha256())
			throw std::invalid_argument("strong FSI accepted material differs from patch reference");
		for(double value:{options.absolute_displacement_tolerance_m,options.relative_displacement_tolerance})
			if(!std::isfinite(value)||value<0)throw std::invalid_argument("invalid strong FSI tolerance");
		if(!std::isfinite(options.reference_displacement_scale_m)||options.reference_displacement_scale_m<=0)
			throw std::invalid_argument("invalid strong FSI displacement scale");
		std::map<std::uint64_t,std::array<double,3>> positions,sums;
		for(const auto& position:layout.reference_positions)positions.emplace(position.global_node_id,position.position_m);
		for(const auto& triangle:layout.reference_triangles) {
			const auto& a=positions.at(triangle[0]);const auto& b=positions.at(triangle[1]);const auto& c=positions.at(triangle[2]);
			std::array<double,3> x{},y{},cross{};
			for(int axis=0;axis<3;++axis) { x[axis]=b[axis]-a[axis];y[axis]=c[axis]-a[axis]; }
			for(int axis=0;axis<3;++axis)cross[axis]=x[(axis+1)%3]*y[(axis+2)%3]-x[(axis+2)%3]*y[(axis+1)%3];
			for(auto node:triangle)for(int axis=0;axis<3;++axis)sums[node][axis]+=cross[axis];
		}
		current.interface=map.Interface().id;
		current.stamp={context.EndTime(),context.step,0,layout.reference_mesh_identity_sha256,
			layout.layout_identity_sha256,BuildDistributedSurfacePartitionIdentitySha256(layout),{}};
		for(auto id:layout.owned_global_node_ids) {
			auto normal=sums.at(id);const double length=std::hypot(normal[0],std::hypot(normal[1],normal[2]));
			if(!std::isfinite(length)||length<=0)throw std::invalid_argument("invalid FSI reference normal");
			for(auto& value:normal)value/=length;
			normals.push_back(normal);
			const auto source=map.SourceVertexForGlobalNode(id);
			const bool clamp=std::binary_search(map.ConfiguredClampedGlobalNodeIds().begin(),map.ConfiguredClampedGlobalNodeIds().end(),id);
			clamped.push_back(clamp);
			std::array<double,3> displacement{},velocity{};
			for(int axis=0;axis<3;++axis) {
				displacement[axis]=committed.SourceVerticesM()[source][axis]-committed.ReferenceMaterialVerticesM()[source][axis];
				velocity[axis]=clamp?0.:committed.SourceVertexVelocitiesMPerS()[source][axis];
			}
			accepted_displacement.push_back(displacement);
			for(int axis=0;axis<3;++axis)displacement[axis]=clamp?0.:displacement[axis]+context.dt_s*velocity[axis];
			current.displacement_m.push_back(displacement);current.velocity_m_per_s.push_back(velocity);
		}
		material_identity=committed.ContentIdentitySha256();
		weights=layout.owned_reference_lumped_areas_m2;
		Sha256 hash;distributed_surface_detail::AppendString(hash,material_identity);
		hash.AppendLittleEndian64(context.step);hash.AppendNormalizedDouble(context.start_time_s);hash.AppendNormalizedDouble(context.dt_s);
		hash.AppendLittleEndian64(options.maximum_iterations);hash.AppendLittleEndian64(options.fail_before_finalize_for_testing);
		for(double value:{options.absolute_displacement_tolerance_m,options.relative_displacement_tolerance,options.reference_displacement_scale_m})hash.AppendNormalizedDouble(value);
		identity=hash.Hex();
	});
	RequireCollectiveSameText(comm,"distributed strong FSI step agreement",identity);
	DistributedWeightedAitken aitken(comm,current.stamp.partition_identity_sha256,
		std::move(weights),options.aitken_controls);
	const auto stamp=[&](SurfaceKinematics& publication,std::uint64_t iteration) {
		publication.stamp.coupling_iteration=iteration;
		Sha256 hash;distributed_surface_detail::AppendString(hash,"DistributedStrongFsiIterate/v1");
		distributed_surface_detail::AppendString(hash,identity);
		distributed_surface_detail::AppendString(hash,publication.stamp.partition_identity_sha256);
		hash.AppendLittleEndian64(iteration);
		for(const auto* values:{&publication.displacement_m,&publication.velocity_m_per_s})
			for(const auto& vector:*values)for(double value:vector)hash.AppendNormalizedDouble(value);
		publication.stamp.producer_state_identity_sha256=hash.Hex();ValidateSurfaceKinematics(publication,layout);
	};
	try {
		for(std::uint64_t iteration=0;iteration<options.maximum_iterations;++iteration) {
			context.coupling_iteration=iteration;
			CollectiveLocalStage(comm,"distributed strong FSI iterate stamp",[&] { stamp(current,iteration); });
			fluid.SolveTrial(context,current,current.stamp,material_identity);
			const auto& traction=fluid.TrialTraction();
			structure.SolveTrial(context,traction,traction.stamp,traction.projection_identity_sha256);
			std::vector<double> raw,current_scalar,residual;
			double local_velocity=0.,global_velocity=0.;
			CollectiveLocalStage(comm,"distributed strong FSI residual",[&] {
				const auto& publication=structure.TrialKinematics();ValidateSurfaceKinematics(publication,layout);
				if(!(publication.interface==map.Interface().id)||publication.stamp.step!=context.step
					|| publication.stamp.time_s!=context.EndTime()||publication.stamp.coupling_iteration!=iteration)
					throw std::invalid_argument("strong FSI structure publication epoch mismatch");
				for(std::size_t node=0;node<normals.size();++node) {
					double a=0,b=0;
					for(int axis=0;axis<3;++axis) {
						a+=publication.displacement_m[node][axis]*normals[node][axis];b+=current.displacement_m[node][axis]*normals[node][axis];
						local_velocity=std::max(local_velocity,std::abs(publication.velocity_m_per_s[node][axis]-current.velocity_m_per_s[node][axis])*context.dt_s);
					}
					raw.push_back(a);current_scalar.push_back(b);residual.push_back(a-b);
				}
				if(!std::isfinite(local_velocity))throw std::overflow_error("strong FSI velocity residual overflow");
			});
			MPI_Allreduce(&local_velocity,&global_velocity,1,MPI_DOUBLE,MPI_MAX,comm);
			DistributedStrongFsiIteration diagnostic;
			diagnostic.displacement=EvaluateDistributedFsiConvergence(comm,layout.owned_reference_lumped_areas_m2,
				residual,raw,current_scalar,options.reference_displacement_scale_m,
				options.absolute_displacement_tolerance_m,options.relative_displacement_tolerance);
			diagnostic.maximum_velocity_residual_times_dt_m=global_velocity;
			if(diagnostic.displacement.converged && global_velocity<=diagnostic.displacement.convergence_threshold_m) {
				CollectiveLocalStage(comm,"distributed strong FSI result preparation",[&] {
					result.history.push_back(diagnostic);result.converged=true;
					if(options.fail_before_finalize_for_testing)throw std::runtime_error("injected strong FSI precommit failure");
				});
				DistributedFsiCommitCoordinator::Commit(comm,context,fluid,structure);
				return result;
			}
			fluid.AbortTrial();structure.AbortTrial();
			const auto proposal=aitken.Propose(current_scalar,residual,options.reference_displacement_scale_m);
			SurfaceKinematics next;
			CollectiveLocalStage(comm,"distributed strong FSI relaxed publication",[&] {
				next=current;
				for(std::size_t node=0;node<normals.size();++node)for(int axis=0;axis<3;++axis) {
					const double value=clamped[node]?0.:proposal.next[node]*normals[node][axis];
					next.displacement_m[node][axis]=value;
					// Fluid wall velocity follows its accepted geometry, including
					// the small interface mismatch left by the previous FSI solve.
					next.velocity_m_per_s[node][axis]=clamped[node]?0.:(value-accepted_displacement[node][axis])/context.dt_s;
				}
				stamp(next,iteration+1);diagnostic.relaxation=proposal.relaxation_factor;
				result.history.push_back(diagnostic);
			});
			aitken.AcceptApplied(proposal,residual,proposal.relaxation_factor);
			using std::swap;swap(current,next);
		}
		throw std::runtime_error("distributed strong FSI iteration did not converge");
	} catch(...) { fluid.AbortTrial();structure.AbortTrial();throw; }
}

} // namespace iga
#endif
