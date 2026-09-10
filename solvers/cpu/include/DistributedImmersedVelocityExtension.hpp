#ifndef IGA_DISTRIBUTED_IMMERSED_VELOCITY_EXTENSION_HPP
#define IGA_DISTRIBUTED_IMMERSED_VELOCITY_EXTENSION_HPP

#include "ImmersedVelocityExtension.hpp"
#include "OwnedPointValues.hpp"
#include "PetscSolverOptions.hpp"
#include "RuntimeCleanup.hpp"
#include "RuntimeConstruction.hpp"
#include "ExecutionResources.hpp"
#include "PetscReadArray.hpp"
#include <petscksp.h>
#include <map>

namespace iga {
// Geometry/traces are currently replicated. Numerical matrix rows, anchor
// tuples and solved unknowns are distributed. All operations and Close are
// collective; source fields must be independently certified as committed.
class DistributedImmersedVelocityExtension
{
public:
	DistributedImmersedVelocityExtension(MPI_Comm comm,const MovingCutGeometry& old_geometry,
		const ImmersedActiveLayout& old_layout,const MovingCutGeometry& new_geometry,
		const ImmersedActiveLayout& new_layout,std::uint32_t layers,
		ImmersedVelocityExtensionOptions options={},PetscOptions source_options=nullptr)
		: comm_(comm),source_options_(source_options)
	{
		int rank=0,ranks=1;MPI_Comm_rank(comm_,&rank);MPI_Comm_size(comm_,&ranks);
		std::string identity;
		CollectiveLocalStage(comm_,"distributed extension topology",[&] {
			RequirePetscRealDouble();
			topology_.InitializeTopology(old_geometry,old_layout,new_geometry,new_layout,layers,options);
			target_ids_=new_layout.NodeIds();
			Sha256 hash;
			for(const auto& value:{topology_.HashOperator(),old_layout.HashSha256(),new_layout.HashSha256()})
				hash.Append(value.data(),value.size());
			hash.AppendNormalizedDouble(topology_.old_time_s_);hash.AppendNormalizedDouble(topology_.new_time_s_);
			hash.AppendNormalizedDouble(options.residual_tolerance);identity=hash.Hex();
			const auto count=topology_.unknown_nodes_.size();
			if(count>static_cast<std::size_t>(std::numeric_limits<PetscInt>::max()))throw std::length_error("extension row count overflow");
			begin_=static_cast<PetscInt>(count*static_cast<std::size_t>(rank)/ranks);
			end_=static_cast<PetscInt>(count*static_cast<std::size_t>(rank+1)/ranks);
			rows_.resize(end_-begin_);
			for(std::size_t i=0;i<topology_.anchor_nodes_.size();++i)
				if(i%ranks==static_cast<std::size_t>(rank))needed_anchors_.push_back(topology_.anchor_nodes_[i]);
			for(std::size_t f=0;f<topology_.faces_.size();++f) {
				const auto& face=topology_.faces_[f];bool needed=false;
				for(auto id:face.nodes){const auto row=UnknownRow(id);needed=needed||(row>=begin_&&row<end_);}
				if(!needed)continue;
				faces_.push_back(f);
				for(auto id:face.nodes)if(UnknownRow(id)<0)needed_anchors_.push_back(id);
				for(std::size_t q=0;q<face.alpha.size();++q)for(std::size_t a=0;a<face.nodes.size();++a) {
					const auto row=UnknownRow(face.nodes[a]);if(row<begin_||row>=end_)continue;
					for(std::size_t b=0;b<face.nodes.size();++b) {
						const auto col=UnknownRow(face.nodes[b]);if(col<0)continue;
						auto& entry=rows_[row-begin_][col];
						entry=ImmersedVelocityExtension::AddProduct(entry,face.alpha[q],face.jumps[q][std::min(a,b)],face.jumps[q][std::max(a,b)],"extension sparse entry overflow");
					}
				}
			}
			std::sort(needed_anchors_.begin(),needed_anchors_.end());
			needed_anchors_.erase(std::unique(needed_anchors_.begin(),needed_anchors_.end()),needed_anchors_.end());
		});
		RequireCollectiveSameText(comm_,"distributed extension topology agreement",identity);
		try { BuildMatrix(); } catch(...) { Release();throw; }
	}
	~DistributedImmersedVelocityExtension() { Release(); }
	DistributedImmersedVelocityExtension(const DistributedImmersedVelocityExtension&)=delete;
	DistributedImmersedVelocityExtension& operator=(const DistributedImmersedVelocityExtension&)=delete;
	void RequireTargetLayout(MPI_Comm communicator,const ImmersedActiveLayout& layout) const
	{
		int relation=MPI_UNEQUAL;
		if(closed_||layout.HashSha256()!=topology_.new_layout_identity_
			||MPI_Comm_compare(comm_,communicator,&relation)!=MPI_SUCCESS
			||(relation!=MPI_IDENT&&relation!=MPI_CONGRUENT))
			throw std::invalid_argument("extension target layout or communicator differs");
	}

	// Extract only locally owned node rows from the runtime's committed Vec.
	// Controller/gauge rows do not enter the extension. The runtime remains
	// responsible for certifying that this Vec is its committed state.
	std::vector<double> ExtendCommitted(Vec committed,const ImmersedActiveLayout& source_layout,
		double source_time,std::uint64_t source_index,double target_time,
		std::uint64_t target_index,double dt,const std::vector<std::uint64_t>& target_ids,
		PointIdentityLimits limits={})
	{
		std::vector<std::uint64_t> owned;
		std::vector<double> coefficients;
		std::string clock;
		CollectiveLocalStage(comm_,"extension committed state extraction",[&] {
			if(closed_||!committed)throw std::invalid_argument("extension committed source is unavailable");
			if(source_layout.HashSha256()!=topology_.old_layout_identity_)
				throw std::invalid_argument("extension committed source layout differs");
			if(source_time!=topology_.old_time_s_||target_time!=topology_.new_time_s_
				||!std::isfinite(dt)||!(dt>0.)||target_time!=CheckedTransientTargetTime(source_time,dt)
				||source_index==std::numeric_limits<std::uint64_t>::max()||target_index!=source_index+1)
				throw std::invalid_argument("extension committed clock differs");
			MPI_Comm source_comm=MPI_COMM_NULL;int relation=MPI_UNEQUAL;
			Local(PetscObjectGetComm(reinterpret_cast<PetscObject>(committed),&source_comm));
			if(MPI_Comm_compare(comm_,source_comm,&relation)!=MPI_SUCCESS
				||(relation!=MPI_IDENT&&relation!=MPI_CONGRUENT))
				throw std::invalid_argument("extension committed communicator differs");
			PetscInt total=0,first=0,last=0;
			Local(VecGetSize(committed,&total));Local(VecGetOwnershipRange(committed,&first,&last));
			if(total<0||static_cast<std::size_t>(total)!=source_layout.Rows())
				throw std::invalid_argument("extension committed row count differs");
			const auto physical=static_cast<PetscInt>(source_layout.NodeFieldRows());
			const auto node_first=std::min(first,physical),node_last=std::min(last,physical);
			if(node_first%4||node_last%4)throw std::invalid_argument("extension committed ownership splits a node");
			PetscReadArray view;view.Acquire(committed);
			for(PetscInt row=node_first;row<node_last;row+=4) {
				owned.push_back(source_layout.NodeIds()[row/4]);
				for(PetscInt component=0;component<4;++component)
					coefficients.push_back(PetscRealPart(view.Data()[row-first+component]));
			}
			view.Restore();
			Sha256 hash;hash.AppendNormalizedDouble(source_time);hash.AppendLittleEndian64(source_index);
			hash.AppendNormalizedDouble(target_time);hash.AppendLittleEndian64(target_index);hash.AppendNormalizedDouble(dt);
			clock=hash.Hex();
		});
		RequireCollectiveSameText(comm_,"extension committed clock agreement",clock);
		return Extend(owned,coefficients,target_ids,limits);
	}

	// Stage coefficients in the target Vec's local row order, without writing
	// either Vec. Controllers follow stable IDs; the target gauge is canonical
	// positive zero, matching the serial moving-runtime seed policy.
	std::vector<PetscScalar> BuildOwnedTargetSeed(Vec committed,const ImmersedActiveLayout& source_layout,
		Vec target_template,const ImmersedActiveLayout& target_layout,double source_time,
		std::uint64_t source_index,double target_time,std::uint64_t target_index,double dt,
		PointIdentityLimits limits={})
	{
		PetscInt first=0,last=0,node_first=0,node_last=0;
		std::vector<std::uint64_t> nodes,ports;
		CollectiveLocalStage(comm_,"extension target seed layout",[&] {
			RequireTargetLayout(comm_,target_layout);
			if(!target_template)throw std::invalid_argument("extension target Vec is null");
			MPI_Comm target_comm=MPI_COMM_NULL;int relation=MPI_UNEQUAL;
			Local(PetscObjectGetComm(reinterpret_cast<PetscObject>(target_template),&target_comm));
			if(MPI_Comm_compare(comm_,target_comm,&relation)!=MPI_SUCCESS
				||(relation!=MPI_IDENT&&relation!=MPI_CONGRUENT))throw std::invalid_argument("extension target Vec communicator differs");
			PetscInt rows=0;Local(VecGetSize(target_template,&rows));Local(VecGetOwnershipRange(target_template,&first,&last));
			if(rows<0||static_cast<std::size_t>(rows)!=target_layout.Rows())throw std::invalid_argument("extension target Vec rows differ");
			const auto physical=static_cast<PetscInt>(target_layout.NodeFieldRows());
			node_first=std::min(first,physical);node_last=std::min(last,physical);
			if(node_first%4||node_last%4)throw std::invalid_argument("extension target Vec splits a node");
			for(PetscInt row=node_first;row<node_last;row+=4)nodes.push_back(target_layout.NodeIds()[row/4]);
			for(auto id:target_layout.PortIds()) {
				const auto row=target_layout.ControllerRow(id);
				if(row>=static_cast<std::size_t>(first)&&row<static_cast<std::size_t>(last))ports.push_back(id);
			}
		});
		const auto coefficients=ExtendCommitted(committed,source_layout,source_time,source_index,target_time,target_index,dt,nodes,limits);
		std::vector<std::uint64_t> owned_ports;std::vector<double> port_values;
		CollectiveLocalStage(comm_,"extension source controllers",[&] {
			PetscInt source_first=0,source_last=0;Local(VecGetOwnershipRange(committed,&source_first,&source_last));
			PetscReadArray view;view.Acquire(committed);
			for(auto id:source_layout.PortIds()) {
				const auto row=source_layout.ControllerRow(id);
				if(row>=static_cast<std::size_t>(source_first)&&row<static_cast<std::size_t>(source_last)) {
					owned_ports.push_back(id);port_values.push_back(PetscRealPart(view.Data()[row-source_first]));
				}
			}
			view.Restore();
		});
		const auto controllers=FetchOwnedPointValues(comm_,owned_ports,port_values,ports,1,limits);
		std::vector<PetscScalar> seed;
		CollectiveLocalStage(comm_,"extension target seed candidate",[&] {
			seed.assign(last-first,0.);
			std::copy(coefficients.begin(),coefficients.end(),seed.begin());
			for(std::size_t i=0;i<ports.size();++i)seed[target_layout.ControllerRow(ports[i])-first]=controllers[i];
		});
		return seed;
	}

	// Four coefficients per node: velocity xyz and scalar pressure warm start.
	// Returns only requested target IDs, in request order, including duplicates.
	std::vector<double> Extend(const std::vector<std::uint64_t>& owned_ids,
		const std::vector<double>& owned_coefficients,const std::vector<std::uint64_t>& target_ids,
		PointIdentityLimits limits={})
	{
		CollectiveLocalStage(comm_,"distributed extension field preflight",[&] {
			if(closed_)throw std::logic_error("distributed extension is closed");
			for(auto id:owned_ids)if(id>static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max())
				||!std::binary_search(topology_.anchor_nodes_.begin(),topology_.anchor_nodes_.end(),static_cast<std::int32_t>(id)))
				throw std::invalid_argument("extension source is not an anchor");
			for(auto id:target_ids)if(id>static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max())
				||!std::binary_search(target_ids_.begin(),target_ids_.end(),static_cast<std::int32_t>(id)))
				throw std::invalid_argument("extension target is outside active layout");
		});
		const auto anchors=FetchOwnedPointValues(comm_,owned_ids,owned_coefficients,needed_anchors_,4,limits);
		std::vector<std::uint64_t> output_ids;std::vector<double> output;
		std::array<double,4> residuals{};
		CollectiveLocalStage(comm_,"distributed extension output staging",[&] {
			output_ids=owned_ids;output=owned_coefficients;
			for(PetscInt row=begin_;row<end_;++row)output_ids.push_back(topology_.unknown_nodes_[row]);
			output.resize(output_ids.size()*4);
		});
		if(matrix_)for(int component=0;component<4;++component) {
			std::vector<double> rhs;
			CollectiveLocalStage(comm_,"distributed extension RHS",[&] {
				rhs.assign(end_-begin_,0.);
				for(auto index:faces_) {
					const auto& face=topology_.faces_[index];
					for(std::size_t q=0;q<face.alpha.size();++q) {
						long double trace=0.;
						for(std::size_t a=0;a<face.nodes.size();++a)if(UnknownRow(face.nodes[a])<0) {
							const auto pos=std::lower_bound(needed_anchors_.begin(),needed_anchors_.end(),face.nodes[a])-needed_anchors_.begin();
							trace+=static_cast<long double>(face.jumps[q][a])*anchors[4*pos+component];
						}
						const double value=ImmersedVelocityExtension::FiniteDouble(trace,"extension anchor trace overflow");
						for(std::size_t a=0;a<face.nodes.size();++a) {
							const auto row=UnknownRow(face.nodes[a]);if(row<begin_||row>=end_)continue;
							rhs[row-begin_]=ImmersedVelocityExtension::AddProduct(rhs[row-begin_],-face.alpha[q],face.jumps[q][a],value,"extension RHS overflow");
						}
					}
				}
				PetscScalar* values=nullptr;Local(VecGetArray(rhs_,&values));
				std::copy(rhs.begin(),rhs.end(),values);Local(VecRestoreArray(rhs_,&values));
			});
			Check("extension scaled RHS",VecPointwiseMult(work_,rhs_,scaling_));
			solver_options_->Call("extension solve",[&] { return KSPSolve(solver_,work_,solution_); });
			CollectiveLocalStage(comm_,"extension convergence reason",[&] {
				KSPConvergedReason reason;Local(KSPGetConvergedReason(solver_,&reason));
				if(reason<=0)throw std::runtime_error("extension KSP did not converge");
			});
			Check("extension unscale solution",VecPointwiseMult(solution_,solution_,scaling_));
			Check("extension original matrix action",MatMult(matrix_,solution_,work_));
			Check("extension original residual",VecAXPY(work_,-1.,rhs_));
			PetscReal rn=0,xn=0,bn=0;
			Check("extension residual norm",VecNorm(work_,NORM_2,&rn));
			Check("extension solution norm",VecNorm(solution_,NORM_2,&xn));
			Check("extension RHS norm",VecNorm(rhs_,NORM_2,&bn));
			CollectiveLocalStage(comm_,"extension physical residual gate",[&] {
				const long double denominator=static_cast<long double>(matrix_norm_)*xn+bn;
				if(!std::isfinite(rn)||!std::isfinite(denominator))throw std::runtime_error("nonfinite extension residual");
				residuals[component]=denominator>0.?static_cast<double>(rn/denominator):(rn==0.?0.:std::numeric_limits<double>::infinity());
				if(residuals[component]>topology_.options_.residual_tolerance)throw std::runtime_error("extension physical residual gate failed");
				const PetscScalar* values=nullptr;Local(VecGetArrayRead(solution_,&values));
				for(PetscInt i=0;i<end_-begin_;++i)output[4*(owned_ids.size()+i)+component]=PetscRealPart(values[i]);
				Local(VecRestoreArrayRead(solution_,&values));
			});
		}
		auto result=FetchOwnedPointValues(comm_,output_ids,output,target_ids,4,limits);
		residuals_=residuals;return result;
	}
	const std::array<double,4>& Residuals() const noexcept { return residuals_; }
	PetscInt OwnedRows() const noexcept { return end_-begin_; }
	void Close() { Release();cleanup_.Check(comm_,"distributed extension close"); }
private:
	PetscInt UnknownRow(std::int32_t id) const
	{
		const auto& ids=topology_.unknown_nodes_;const auto it=std::lower_bound(ids.begin(),ids.end(),id);
		return it==ids.end()||*it!=id?-1:static_cast<PetscInt>(it-ids.begin());
	}
	static void Local(PetscErrorCode error) { if(error)throw std::runtime_error("extension local PETSc operation failed"); }
	void Check(const char* stage,PetscErrorCode code) const { RequireCollectivePetscSuccess(comm_,stage,code); }
	void BuildMatrix()
	{
		const auto count=static_cast<PetscInt>(topology_.unknown_nodes_.size());if(!count)return;
		std::vector<PetscInt> diagonal,off;
		CollectiveLocalStage(comm_,"extension sparse preallocation",[&] {
			diagonal.resize(rows_.size());off.resize(rows_.size());
			for(std::size_t i=0;i<rows_.size();++i)for(const auto& entry:rows_[i])
				++(entry.first>=begin_&&entry.first<end_?diagonal[i]:off[i]);
		});
		Check("extension matrix create",MatCreateAIJ(comm_,end_-begin_,end_-begin_,count,count,0,diagonal.data(),0,off.data(),&matrix_));
		Check("extension strict allocation",MatSetOption(matrix_,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE));
		CollectiveLocalStage(comm_,"extension owned row assembly",[&] {
			for(std::size_t i=0;i<rows_.size();++i)for(const auto& entry:rows_[i])Local(MatSetValue(matrix_,begin_+i,entry.first,entry.second,INSERT_VALUES));
		});
		Check("extension assembly begin",MatAssemblyBegin(matrix_,MAT_FINAL_ASSEMBLY));
		Check("extension assembly end",MatAssemblyEnd(matrix_,MAT_FINAL_ASSEMBLY));
		Check("extension matrix norm",MatNorm(matrix_,NORM_INFINITY,&matrix_norm_));
		Check("extension vector create",MatCreateVecs(matrix_,&solution_,&rhs_));
		Check("extension work create",VecDuplicate(rhs_,&work_));
		Check("extension scaling create",VecDuplicate(rhs_,&scaling_));
		Check("extension diagonal",MatGetDiagonal(matrix_,scaling_));
		CollectiveLocalStage(comm_,"extension positive diagonal",[&] {
			PetscScalar* values=nullptr;Local(VecGetArray(scaling_,&values));
			bool valid=true;
			for(PetscInt i=0;i<end_-begin_;++i) {
				const double d=PetscRealPart(values[i]);valid=valid&&std::isfinite(d)&&d>0.;
				values[i]=d>0.?1./std::sqrt(d):0.;
			}
			Local(VecRestoreArray(scaling_,&values));if(!valid)throw std::runtime_error("extension nonpositive diagonal");
		});
		Check("extension scaled matrix copy",MatDuplicate(matrix_,MAT_COPY_VALUES,&scaled_));
		Check("extension diagonal scaling",MatDiagonalScale(scaled_,scaling_,scaling_));
		Check("extension solver create",KSPCreate(comm_,&solver_));
		Check("extension solver type",KSPSetType(solver_,KSPPREONLY));
		PC pc=nullptr;Check("extension PC lookup",KSPGetPC(solver_,&pc));
		Check("extension PC type",PCSetType(pc,PCLU));
		Check("extension LU backend",PCFactorSetMatSolverType(pc,MATSOLVERMUMPS));
		Check("extension solver operator",KSPSetOperators(solver_,scaled_,scaled_));
		solver_options_=AllocateCollectiveRuntime<PetscSolverOptions>(comm_,comm_,"immersed_extension_",source_options_,PetscOptionEntries{},std::set<std::string>{},"immersed_extension_",false);
		solver_options_->Attach(solver_);
		solver_options_->Attach(scaled_);
		solver_options_->Call("extension solver options",[&] { return KSPSetFromOptions(solver_); });
		solver_options_->RecordUsed();
	}
	void Release() noexcept
	{
		closed_=true;
		cleanup_.Observe("extension KSPDestroy",KSPDestroy(&solver_));
		for(auto* value:{&solution_,&rhs_,&work_,&scaling_})cleanup_.Observe("extension VecDestroy",VecDestroy(value));
		cleanup_.Observe("extension MatDestroy scaled",MatDestroy(&scaled_));cleanup_.Observe("extension MatDestroy original",MatDestroy(&matrix_));
	}
	MPI_Comm comm_;ImmersedVelocityExtension topology_;std::vector<std::int32_t> target_ids_;
	PetscOptions source_options_=nullptr;
	PetscInt begin_=0,end_=0;std::vector<std::map<PetscInt,double>> rows_;std::vector<std::size_t> faces_;
	std::vector<std::uint64_t> needed_anchors_;std::array<double,4> residuals_{};
	Mat matrix_=nullptr,scaled_=nullptr;Vec solution_=nullptr,rhs_=nullptr,work_=nullptr,scaling_=nullptr;
	KSP solver_=nullptr;PetscReal matrix_norm_=0.;std::unique_ptr<PetscSolverOptions> solver_options_;
	RuntimeCleanupResult cleanup_;bool closed_=false;
};
} // namespace iga
#endif
