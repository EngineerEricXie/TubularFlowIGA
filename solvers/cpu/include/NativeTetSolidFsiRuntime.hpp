#ifndef IGA_NATIVE_TET_SOLID_FSI_RUNTIME_HPP
#define IGA_NATIVE_TET_SOLID_FSI_RUNTIME_HPP

// Structure-side adapter joining the native tetrahedral hyperelastic solver to
// the existing field-valued strong-coupling transaction.
#include "FsiDomainRuntime.hpp"
#include "NativeTetSolidStaticSolver.hpp"
#include "StrongFluidStructureCoupling.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace iga {

class NativeTetSolidFsiRuntime final : public FsiStructureDomainRuntime {
public:
	NativeTetSolidFsiRuntime(const NativeTetSolidFsiRuntime&) = delete;
	NativeTetSolidFsiRuntime& operator=(const NativeTetSolidFsiRuntime&) = delete;
	NativeTetSolidFsiRuntime(NativeTetSolidFsiRuntime&&) = delete;
	NativeTetSolidFsiRuntime& operator=(NativeTetSolidFsiRuntime&&) = delete;

	NativeTetSolidFsiRuntime(std::string domain_id, std::string subsystem_id,
		FsiCouplingEdge edge, DistributedSurfaceInterface structure_surface,
		DistributedSurfaceInterface fluid_surface, DistributedSurfaceLayout structure_layout,
		DistributedSurfaceLayout fluid_layout, NativeTetMesh mesh, int interface_label,
		NativeTetSolidMaterial material, std::map<std::size_t, double> constraints,
		NativeTetSolidStaticOptions options = {}, std::vector<double> initial_displacement_m = {},
		std::vector<double> initial_velocity_m_per_s = {})
		: domain_id_(std::move(domain_id)), subsystem_id_(std::move(subsystem_id)),
		  edge_(std::move(edge)), catalog_{std::move(structure_surface)},
		  fluid_surface_(std::move(fluid_surface)), structure_layout_(std::move(structure_layout)),
		  fluid_layout_(std::move(fluid_layout)), mesh_(std::move(mesh)),
		  interface_label_(interface_label), material_(material), constraints_(std::move(constraints)),
		  options_(options), lifecycle_(domain_id_, subsystem_id_, edge_, catalog_.front(),
			fluid_surface_, structure_layout_, fluid_layout_)
	{
		if (!(catalog_.front().id == edge_.structure) || !(fluid_surface_.id == edge_.fluid))
			throw std::runtime_error("native tetrahedral solid FSI surfaces do not match the edge");
		if (catalog_.front().id.domain_id != domain_id_
			|| catalog_.front().id.subsystem_id != subsystem_id_)
			throw std::runtime_error("native tetrahedral solid FSI local endpoint is inconsistent");
		if (structure_layout_.partition_count != 1 || structure_layout_.partition_rank != 0)
			throw std::runtime_error("native tetrahedral solid FSI adapter is currently single-partition");
		if (catalog_.front().boundary_labels != std::vector<std::int64_t>{interface_label_})
			throw std::runtime_error("native tetrahedral solid FSI label differs from its surface contract");
		ValidateInterfaceTopology();
		const std::size_t dofs = 3*mesh_.points.size();
		if (initial_displacement_m.empty()) initial_displacement_m.assign(dofs, 0.0);
		if (initial_velocity_m_per_s.empty()) initial_velocity_m_per_s.assign(dofs, 0.0);
		if (initial_displacement_m.size() != dofs)
			throw std::runtime_error("native tetrahedral solid FSI initial state size mismatch");
		if (initial_velocity_m_per_s.size() != dofs)
			throw std::runtime_error("native tetrahedral solid FSI initial velocity size mismatch");
		for (double value : initial_displacement_m)
			if (!std::isfinite(value)) throw std::runtime_error("native tetrahedral solid FSI initial state is nonfinite");
		for (double value : initial_velocity_m_per_s)
			if (!std::isfinite(value)) throw std::runtime_error("native tetrahedral solid FSI initial velocity is nonfinite");
		for (const auto& constraint : constraints_) {
			if (constraint.first >= dofs || !std::isfinite(constraint.second))
				throw std::runtime_error("native tetrahedral solid FSI constraint is invalid");
			if (initial_displacement_m[constraint.first] != constraint.second)
				throw std::runtime_error("native tetrahedral solid FSI initial state violates a constraint");
			if (initial_velocity_m_per_s[constraint.first] != 0.0)
				throw std::runtime_error("native tetrahedral solid FSI initial velocity violates a constraint");
		}
		for (auto id : structure_layout_.owned_global_node_ids)
			for (int component = 0; component < 3; ++component) {
				const auto found = constraints_.find(3*static_cast<std::size_t>(id)+component);
				if (found != constraints_.end() && found->second != 0.0)
					throw std::runtime_error("native tetrahedral solid FSI interface supports only zero prescribed constraints");
			}
		committed_displacement_m_ = std::move(initial_displacement_m);
		committed_velocity_m_per_s_ = std::move(initial_velocity_m_per_s);
		reference_normals_ = BuildReferenceNormals();
		model_identity_sha256_ = BuildModelIdentity();
		committed_state_identity_sha256_ = BuildStateIdentity(committed_displacement_m_,
			committed_velocity_m_per_s_);
	}

	const std::vector<DistributedSurfaceInterface>& SurfaceInterfaces() const noexcept override
	{ return catalog_; }
	const FsiTrialLifecycle& Lifecycle() const noexcept { return lifecycle_; }
	const FsiCouplingEdge& StrongCouplingEdge() const noexcept { return edge_; }
	const DistributedSurfaceLayout& StrongCouplingLayout() const noexcept { return structure_layout_; }
	const std::vector<double>& CommittedDisplacementM() const noexcept { return committed_displacement_m_; }
	const std::vector<double>& CommittedVelocityMPerS() const noexcept { return committed_velocity_m_per_s_; }
	const std::string& ModelIdentitySha256() const noexcept { return model_identity_sha256_; }
	const std::string& CommittedStateIdentitySha256() const noexcept
	{ return committed_state_identity_sha256_; }

	StrongCouplingStructureSnapshot StrongCouplingCommittedSnapshot() const
	{
		StrongCouplingStructureSnapshot result;
		const auto& ids = structure_layout_.owned_global_node_ids;
		result.displacement_m.resize(ids.size()); result.velocity_m_per_s.resize(ids.size());
		result.immutable_reference_normals = reference_normals_; result.clamped.resize(ids.size());
		for (std::size_t row = 0; row < ids.size(); ++row) {
			const std::size_t node = static_cast<std::size_t>(ids[row]);
			for (int component = 0; component < 3; ++component) {
				result.displacement_m[row] += committed_displacement_m_[3*node+component]*reference_normals_[row][component];
				result.velocity_m_per_s[row] += committed_velocity_m_per_s_[3*node+component]*reference_normals_[row][component];
			}
			result.clamped[row] = IsFullyZeroConstrained(node);
		}
		result.committed_state_identity_sha256 = committed_state_identity_sha256_;
		result.model_identity_sha256 = model_identity_sha256_;
		return result;
	}

	void BeginMacroStep(const DomainStepContext& step)
	{
		if (trial_ || prepared_) throw std::runtime_error("native tetrahedral solid FSI retains stale trial state");
		lifecycle_.BeginStep(step);
	}
	void BeginStep(const DomainStepContext& step) { BeginMacroStep(step); }
	void BeginCouplingIteration(std::uint64_t iteration,
		const SurfaceFieldStamp& expected_traction, const SurfaceFieldStampEnvelope& output_envelope)
	{
		if (trial_ || prepared_) throw std::runtime_error("native tetrahedral solid FSI retains stale trial state");
		SurfaceFieldStamp input_copy = expected_traction;
		SurfaceFieldStampEnvelope output_copy = output_envelope;
		lifecycle_.BeginIteration(iteration, input_copy, output_copy);
		using std::swap;
		swap(expected_traction_, input_copy); swap(output_envelope_, output_copy);
	}
	void BeginIteration(std::uint64_t iteration, const SurfaceFieldStamp& expected_traction,
		const SurfaceFieldStampEnvelope& output_envelope)
	{ BeginCouplingIteration(iteration, expected_traction, output_envelope); }

	void SetSurfaceTraction(const std::string& interface_id, const SurfaceTraction& traction) override
	{
		if (interface_id != edge_.structure.interface_id)
			throw std::runtime_error("unknown native tetrahedral solid FSI interface");
		SurfaceTraction candidate = traction;
		ValidateFsiStructureTractionInput(edge_, candidate, structure_layout_, expected_traction_);
		lifecycle_.MarkInput(candidate.interface, candidate.stamp); traction_ = std::move(candidate);
	}

	void SolveSolidTrial()
	{
		lifecycle_.RequireSolveAllowed();
		if (!traction_) throw std::runtime_error("native tetrahedral solid FSI lacks an accepted traction publication");
		std::vector<double> external_force(3*mesh_.points.size(), 0.0);
		for (std::size_t row = 0; row < structure_layout_.owned_global_node_ids.size(); ++row) {
			const std::size_t node = static_cast<std::size_t>(structure_layout_.owned_global_node_ids[row]);
			for (int component = 0; component < 3; ++component)
				external_force[3*node+component] += traction_->consistent_nodal_force_n[row][component];
		}
		NativeTetSolidFsiTrial candidate;
		candidate.result = SolveNativeTetSolidStatic(mesh_, material_, external_force, constraints_,
			committed_displacement_m_, options_);
		candidate.velocity_m_per_s.resize(candidate.result.displacement_m.size());
		const double dt = lifecycle_.Context().dt_s;
		for (std::size_t dof = 0; dof < candidate.velocity_m_per_s.size(); ++dof)
			candidate.velocity_m_per_s[dof] =
				(candidate.result.displacement_m[dof]-committed_displacement_m_[dof])/dt;
		candidate.kinematics = BuildKinematics(candidate.result.displacement_m,
			candidate.velocity_m_per_s, *traction_);
		try { lifecycle_.MarkSolved(candidate.kinematics.interface, candidate.kinematics.stamp); }
		catch (...) { throw std::runtime_error("native tetrahedral solid FSI output stamp is inconsistent"); }
		trial_.emplace(std::move(candidate));
	}
	// Compatibility name used by the existing generic strong coordinator.
	void SolveMembraneTrial() { SolveSolidTrial(); }
	void SolveTrial() { SolveSolidTrial(); }

	SurfaceKinematics GetSurfaceKinematics(const std::string& interface_id) const override
	{
		if (interface_id != edge_.structure.interface_id)
			throw std::runtime_error("unknown native tetrahedral solid FSI interface");
		if (!trial_) throw std::runtime_error("native tetrahedral solid FSI has no solved trial");
		lifecycle_.RequireTrialOutput(trial_->kinematics.interface, trial_->kinematics.stamp);
		return trial_->kinematics;
	}
	SurfaceKinematics GetCommittedSurfaceKinematics(const std::string& interface_id,
		const SurfaceFieldStamp& stamp) const
	{
		if (interface_id != edge_.structure.interface_id)
			throw std::runtime_error("unknown native tetrahedral solid FSI interface");
		lifecycle_.RequireCommittedOutput(catalog_.front().id, stamp); return committed_kinematics_;
	}
	StrongCouplingAcceptanceDiagnostics TrialCouplingAcceptanceDiagnostics() const
	{
		if(!trial_)throw std::runtime_error("native tetrahedral solid FSI diagnostics require a trial");
		return {trial_->result.final_free_residual_n,0.0};
	}

	void RejectCouplingIteration()
	{
		if (!lifecycle_.HasActiveStep())
			throw std::runtime_error("native tetrahedral solid FSI reject requires an active step");
		DiscardTrialNoexcept(); lifecycle_.RejectIteration();
	}
	void RollbackTrial() { RejectCouplingIteration(); }
	void PrepareCommitStep()
	{
		if (!trial_) throw std::runtime_error("native tetrahedral solid FSI prepare requires a trial");
		lifecycle_.RequireTrialOutput(trial_->kinematics.interface, trial_->kinematics.stamp);
		NativeTetSolidFsiPrepared candidate;
		candidate.displacement_m = trial_->result.displacement_m;
		candidate.velocity_m_per_s = trial_->velocity_m_per_s;
		candidate.kinematics = trial_->kinematics;
		candidate.state_identity_sha256 = BuildStateIdentity(candidate.displacement_m,
			candidate.velocity_m_per_s);
		try { lifecycle_.PrepareCommit(); prepared_.emplace(std::move(candidate)); }
		catch (...) {
			trial_.reset(); traction_.reset(); prepared_.reset();
			if (lifecycle_.HasActiveStep()) lifecycle_.RejectIteration();
			throw;
		}
	}
	void PrepareCommit() { PrepareCommitStep(); }
	void FinalizeCommitStep()
	{ CoordinatorRequireFinalizeAllowed(); CoordinatorFinalizeCommitNoexcept(); }
	void FinalizeCommit() { FinalizeCommitStep(); }
	void AbortStep() noexcept { DiscardTrialNoexcept(); lifecycle_.AbortStep(); }

	std::string CoordinatorPreparedCommittedKinematicsIdentitySha256() const
	{
		if (!prepared_) throw std::runtime_error("native tetrahedral solid FSI prepared kinematics are unavailable");
		return BuildSurfaceKinematicsIdentitySha256(prepared_->kinematics, structure_layout_);
	}
	std::string CoordinatorPreparedCommittedStateIdentitySha256() const
	{
		if (!prepared_) throw std::runtime_error("native tetrahedral solid FSI prepared state is unavailable");
		return prepared_->state_identity_sha256;
	}

private:
	friend class StrongFluidStructureCouplingAccess;
	struct NativeTetSolidFsiTrial {
		NativeTetSolidStaticResult result; std::vector<double> velocity_m_per_s; SurfaceKinematics kinematics;
	};
	struct NativeTetSolidFsiPrepared {
		std::vector<double> displacement_m, velocity_m_per_s; SurfaceKinematics kinematics;
		std::string state_identity_sha256;
	};
	static_assert(std::is_nothrow_move_constructible<NativeTetSolidFsiPrepared>::value,
		"prepared native tetrahedral solid state must move without throwing");
	void CoordinatorRequireFinalizeAllowed() const
	{
		if (!trial_ || !prepared_)
			throw std::runtime_error("native tetrahedral solid FSI finalize requires prepared numerical state");
		lifecycle_.RequireFinalizeAllowed();
	}
	void CoordinatorFinalizeCommitNoexcept() noexcept
	{
		using std::swap;
		swap(committed_displacement_m_, prepared_->displacement_m);
		swap(committed_velocity_m_per_s_, prepared_->velocity_m_per_s);
		swap(committed_kinematics_, prepared_->kinematics);
		swap(committed_state_identity_sha256_, prepared_->state_identity_sha256);
		lifecycle_.FinalizePreparedCommitNoexcept();
		trial_.reset(); prepared_.reset(); traction_.reset();
	}
	void DiscardTrialNoexcept() noexcept { trial_.reset(); prepared_.reset(); traction_.reset(); }
	bool IsFullyZeroConstrained(std::size_t node) const
	{
		for (int component = 0; component < 3; ++component) {
			const auto found = constraints_.find(3*node+component);
			if (found == constraints_.end() || found->second != 0.0) return false;
		}
		return true;
	}
	void ValidateInterfaceTopology() const
	{
		std::set<std::array<std::uint64_t, 3>> mesh_faces, layout_faces;
		for (const auto& triangle : mesh_.boundary_triangles) if (triangle.boundary_label == interface_label_) {
			std::array<std::uint64_t, 3> face{{triangle.nodes[0], triangle.nodes[1], triangle.nodes[2]}};
			std::sort(face.begin(), face.end());
			if (!mesh_faces.insert(face).second)
				throw std::runtime_error("native tetrahedral solid FSI has duplicate labelled mesh faces");
		}
		for (auto face : structure_layout_.reference_triangles) {
			std::sort(face.begin(), face.end());
			if (!layout_faces.insert(face).second)
				throw std::runtime_error("native tetrahedral solid FSI layout has duplicate interface faces");
		}
		if (mesh_faces.empty() || mesh_faces != layout_faces)
			throw std::runtime_error("native tetrahedral solid FSI layout does not exactly match its labelled faces");
		for (const auto& position : structure_layout_.reference_positions) {
			if (position.global_node_id >= mesh_.points.size())
				throw std::runtime_error("native tetrahedral solid FSI surface node is absent from the volume mesh");
			if (mesh_.points[static_cast<std::size_t>(position.global_node_id)] != position.position_m)
				throw std::runtime_error("native tetrahedral solid FSI reference position differs from the volume mesh");
		}
	}
	std::vector<std::array<double, 3>> BuildReferenceNormals() const
	{
		std::map<std::uint64_t, std::array<double, 3>> accumulated;
		for (const auto& triangle : structure_layout_.reference_triangles) {
			const auto& a=mesh_.points[triangle[0]];const auto& b=mesh_.points[triangle[1]];
			const auto& c=mesh_.points[triangle[2]];
			const std::array<double,3> ab{{b[0]-a[0],b[1]-a[1],b[2]-a[2]}};
			const std::array<double,3> ac{{c[0]-a[0],c[1]-a[1],c[2]-a[2]}};
			const std::array<double,3> cross{{ab[1]*ac[2]-ab[2]*ac[1],ab[2]*ac[0]-ab[0]*ac[2],ab[0]*ac[1]-ab[1]*ac[0]}};
			for (auto node : triangle) for (int component=0;component<3;++component)
				accumulated[node][component]+=cross[component];
		}
		std::vector<std::array<double,3>> result;
		for (auto node : structure_layout_.owned_global_node_ids) {
			auto normal=accumulated.at(node);
			const double norm=std::sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
			if (!(norm>0.0)||!std::isfinite(norm)) throw std::runtime_error("native tetrahedral solid FSI reference normal is invalid");
			for(double& value:normal)value/=norm;
			result.push_back(normal);
		}
		return result;
	}
	std::string BuildModelIdentity() const
	{
		Sha256 hash;distributed_surface_detail::AppendString(hash,"NativeTetSolidFsiModel/v1");
		distributed_surface_detail::AppendString(hash,BuildFsiCouplingEdgeIdentitySha256(edge_));
		distributed_surface_detail::AppendString(hash,structure_layout_.layout_identity_sha256);
		hash.AppendLittleEndian32(static_cast<std::uint32_t>(interface_label_));
		for(const auto& point:mesh_.points)for(double value:point)hash.AppendNormalizedDouble(value);
		for(const auto& cell:mesh_.cells)for(auto node:cell.nodes)hash.AppendLittleEndian32(node);
		for(const auto& triangle:mesh_.boundary_triangles){
			hash.AppendLittleEndian64(triangle.id);
			for(auto node:triangle.nodes)hash.AppendLittleEndian32(node);
			hash.AppendLittleEndian32(static_cast<std::uint32_t>(triangle.boundary_label));
		}
		hash.AppendNormalizedDouble(material_.young_modulus_pa);hash.AppendNormalizedDouble(material_.poisson_ratio);
		hash.AppendNormalizedDouble(material_.density_kg_m3);
		hash.AppendLittleEndian32(static_cast<std::uint32_t>(options_.maximum_iterations));
		hash.AppendNormalizedDouble(options_.relative_tolerance);
		hash.AppendNormalizedDouble(options_.absolute_tolerance_n);
		hash.AppendNormalizedDouble(options_.minimum_line_search);
		hash.AppendLittleEndian64(static_cast<std::uint64_t>(options_.maximum_dofs));
		for(const auto& constraint:constraints_){hash.AppendLittleEndian64(constraint.first);hash.AppendNormalizedDouble(constraint.second);}
		return hash.Hex();
	}
	std::string BuildStateIdentity(const std::vector<double>& displacement,
		const std::vector<double>& velocity) const
	{
		Sha256 hash;distributed_surface_detail::AppendString(hash,"NativeTetSolidFsiState/v1");
		distributed_surface_detail::AppendString(hash,model_identity_sha256_.empty()?BuildModelIdentity():model_identity_sha256_);
		for(double value:displacement)hash.AppendNormalizedDouble(value);
		for(double value:velocity)hash.AppendNormalizedDouble(value);
		return hash.Hex();
	}
	SurfaceKinematics BuildKinematics(const std::vector<double>& displacement,
		const std::vector<double>& velocity,const SurfaceTraction& traction) const
	{
		SurfaceKinematics result;result.interface=catalog_.front().id;
		result.stamp={output_envelope_.time_s,output_envelope_.step,output_envelope_.coupling_iteration,
			output_envelope_.reference_mesh_identity_sha256,output_envelope_.layout_identity_sha256,
			output_envelope_.partition_identity_sha256,{}};
		Sha256 hash;distributed_surface_detail::AppendString(hash,"NativeTetSolidFsiTrial/v1");
		distributed_surface_detail::AppendString(hash,model_identity_sha256_);
		distributed_surface_detail::AppendString(hash,committed_state_identity_sha256_);
		distributed_surface_detail::AppendString(hash,BuildSurfaceTractionIdentitySha256(traction,structure_layout_));
		for(double value:displacement)hash.AppendNormalizedDouble(value);
		for(double value:velocity)hash.AppendNormalizedDouble(value);
		result.stamp.producer_state_identity_sha256=hash.Hex();
		for(auto id:structure_layout_.owned_global_node_ids){std::array<double,3> u{},v{};
			for(int component=0;component<3;++component){u[component]=displacement[3*id+component];v[component]=velocity[3*id+component];}
			result.displacement_m.push_back(u);result.velocity_m_per_s.push_back(v);}
		ValidateSurfaceKinematics(result,structure_layout_);return result;
	}

	std::string domain_id_,subsystem_id_;FsiCouplingEdge edge_;
	std::vector<DistributedSurfaceInterface> catalog_;DistributedSurfaceInterface fluid_surface_;
	DistributedSurfaceLayout structure_layout_,fluid_layout_;NativeTetMesh mesh_;int interface_label_=-1;
	NativeTetSolidMaterial material_;std::map<std::size_t,double> constraints_;NativeTetSolidStaticOptions options_;
	FsiTrialLifecycle lifecycle_;std::vector<std::array<double,3>> reference_normals_;
	std::vector<double> committed_displacement_m_,committed_velocity_m_per_s_;
	std::string model_identity_sha256_,committed_state_identity_sha256_;
	SurfaceFieldStamp expected_traction_;SurfaceFieldStampEnvelope output_envelope_;
	std::optional<SurfaceTraction> traction_;std::optional<NativeTetSolidFsiTrial> trial_;
	std::optional<NativeTetSolidFsiPrepared> prepared_;SurfaceKinematics committed_kinematics_;
};

static_assert(std::is_nothrow_swappable<std::vector<double>>::value,
	"native tetrahedral solid committed vectors must swap without throwing");
static_assert(std::is_nothrow_swappable<SurfaceKinematics>::value,
	"native tetrahedral solid publication must swap without throwing");

} // namespace iga
#endif
