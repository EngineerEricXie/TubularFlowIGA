#ifndef IGA_NATIVE_TET_SOLID_STATIC_SOLVER_HPP
#define IGA_NATIVE_TET_SOLID_STATIC_SOLVER_HPP

#include "NativeTetHyperelasticSolid.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>

namespace iga {

struct NativeTetSolidStaticOptions
{
	int maximum_iterations=30;
	double relative_tolerance=1.0e-8;
	double absolute_tolerance_n=1.0e-10;
	double minimum_line_search=1.0/1024.0;
	std::size_t maximum_dofs=1200;
};

struct NativeTetSolidStaticResult
{
	std::vector<double> displacement_m;
	std::vector<double> reaction_n;
	int iterations=0;
	double initial_free_residual_n=0.0;
	double final_free_residual_n=0.0;
	double strain_energy_j=0.0;
	double external_load_potential_j=0.0;
	double minimum_deformation_jacobian=0.0;
	bool converged=false;
};

namespace native_tet_solid_static_detail {

struct Assembly
{
	std::vector<double> residual;
	std::vector<double> tangent;
	double energy=0.0;
	double minimum_jacobian=std::numeric_limits<double>::infinity();
};

inline void ValidateOptions(const NativeTetSolidStaticOptions& options)
{
	if(options.maximum_iterations<=0||!std::isfinite(options.relative_tolerance)
		||!(options.relative_tolerance>0.0)||!std::isfinite(options.absolute_tolerance_n)
		||!(options.absolute_tolerance_n>0.0)||!std::isfinite(options.minimum_line_search)
		||!(options.minimum_line_search>0.0)||options.minimum_line_search>1.0
		||options.maximum_dofs==0)
		throw std::invalid_argument("native tetrahedral solid static options are invalid");
}

inline Assembly Assemble(const NativeTetMesh& mesh,const std::vector<double>& displacement,
	const std::vector<double>& external_force,const NativeTetSolidMaterial& material)
{
	const std::size_t dofs=3*mesh.points.size();
	Assembly result;result.residual.assign(dofs,0.0);result.tangent.assign(dofs*dofs,0.0);
	for(const auto& cell:mesh.cells){
		std::array<double,12> local{};
		for(std::size_t a=0;a<4;++a)for(int component=0;component<3;++component)
			local[3*a+component]=displacement[3*cell.nodes[a]+component];
		const auto element=BuildNativeTetHyperelasticSolidElement(mesh,cell,local,material);
		result.energy+=element.strain_energy_j;
		result.minimum_jacobian=std::min(result.minimum_jacobian,element.deformation_jacobian);
		for(std::size_t a=0;a<4;++a)for(int i=0;i<3;++i){
			const std::size_t row=3*cell.nodes[a]+i,local_row=3*a+i;
			result.residual[row]+=element.residual_n[local_row];
			for(std::size_t b=0;b<4;++b)for(int k=0;k<3;++k){
				const std::size_t column=3*cell.nodes[b]+k,local_column=3*b+k;
				result.tangent[row*dofs+column]+=element.tangent_n_m[local_row*12+local_column];
			}
		}
	}
	for(std::size_t dof=0;dof<dofs;++dof)result.residual[dof]-=external_force[dof];
	return result;
}

inline double FreeNorm(const std::vector<double>& residual,const std::vector<std::size_t>& free)
{
	long double sum=0.0L;for(auto dof:free)sum+=static_cast<long double>(residual[dof])*residual[dof];
	return std::sqrt(static_cast<double>(sum));
}

inline std::vector<double> SolveDense(std::vector<double> matrix,std::vector<double> rhs)
{
	const std::size_t n=rhs.size();
	for(std::size_t column=0;column<n;++column){
		std::size_t pivot=column;double magnitude=std::abs(matrix[column*n+column]);
		for(std::size_t row=column+1;row<n;++row)if(std::abs(matrix[row*n+column])>magnitude){pivot=row;magnitude=std::abs(matrix[row*n+column]);}
		if(!std::isfinite(magnitude)||magnitude<=1e-14)
			throw std::runtime_error("native tetrahedral solid static tangent is singular");
		if(pivot!=column){for(std::size_t j=column;j<n;++j)std::swap(matrix[column*n+j],matrix[pivot*n+j]);std::swap(rhs[column],rhs[pivot]);}
		for(std::size_t row=column+1;row<n;++row){
			const double factor=matrix[row*n+column]/matrix[column*n+column];
			for(std::size_t j=column;j<n;++j)matrix[row*n+j]-=factor*matrix[column*n+j];
			rhs[row]-=factor*rhs[column];
		}
	}
	for(std::size_t i=n;i-->0;){
		for(std::size_t j=i+1;j<n;++j)rhs[i]-=matrix[i*n+j]*rhs[j];
		rhs[i]/=matrix[i*n+i];if(!std::isfinite(rhs[i]))throw std::runtime_error("native tetrahedral solid static update is nonfinite");
	}
	return rhs;
}

} // namespace native_tet_solid_static_detail

inline NativeTetSolidStaticResult SolveNativeTetSolidStatic(const NativeTetMesh& mesh,
	const NativeTetSolidMaterial& material,const std::vector<double>& external_force_n,
	const std::map<std::size_t,double>& prescribed_displacement_m,
	std::vector<double> initial_displacement_m={},
	const NativeTetSolidStaticOptions& options={})
{
	using namespace native_tet_solid_static_detail;
	ValidateNativeTetSolidMaterial(material);ValidateOptions(options);
	const std::size_t dofs=3*mesh.points.size();
	if(mesh.points.empty()||mesh.cells.empty()||dofs>options.maximum_dofs)
		throw std::invalid_argument("native tetrahedral solid static mesh is empty or exceeds the DOF cap");
	if(external_force_n.size()!=dofs)throw std::invalid_argument("native tetrahedral solid external-force size mismatch");
	for(double value:external_force_n)if(!std::isfinite(value))throw std::invalid_argument("native tetrahedral solid external force is nonfinite");
	if(initial_displacement_m.empty())initial_displacement_m.assign(dofs,0.0);
	if(initial_displacement_m.size()!=dofs)throw std::invalid_argument("native tetrahedral solid initial-state size mismatch");
	for(const auto& item:prescribed_displacement_m)if(item.first>=dofs||!std::isfinite(item.second))
		throw std::invalid_argument("native tetrahedral solid prescribed displacement is invalid");
	std::vector<std::size_t> free;free.reserve(dofs-prescribed_displacement_m.size());
	for(std::size_t dof=0;dof<dofs;++dof)if(!prescribed_displacement_m.count(dof))free.push_back(dof);
	if(free.empty())throw std::invalid_argument("native tetrahedral solid static solve has no free DOF");
	auto displacement=std::move(initial_displacement_m);
	for(const auto& item:prescribed_displacement_m)displacement[item.first]=item.second;
	auto assembly=Assemble(mesh,displacement,external_force_n,material);
	NativeTetSolidStaticResult result;result.initial_free_residual_n=FreeNorm(assembly.residual,free);
	const double threshold=options.absolute_tolerance_n+options.relative_tolerance*result.initial_free_residual_n;
	for(int iteration=0;iteration<=options.maximum_iterations;++iteration){
		const double norm=FreeNorm(assembly.residual,free);
		if(norm<=threshold){result.iterations=iteration;result.converged=true;break;}
		if(iteration==options.maximum_iterations)break;
		const std::size_t n=free.size();std::vector<double> reduced(n*n),rhs(n);
		for(std::size_t i=0;i<n;++i){rhs[i]=-assembly.residual[free[i]];for(std::size_t j=0;j<n;++j)reduced[i*n+j]=assembly.tangent[free[i]*dofs+free[j]];}
		const auto update=SolveDense(std::move(reduced),std::move(rhs));
		bool accepted=false;
		for(double scale=1.0;scale>=options.minimum_line_search;scale*=0.5){
			auto trial=displacement;for(std::size_t i=0;i<n;++i)trial[free[i]]+=scale*update[i];
			try{
				auto candidate=Assemble(mesh,trial,external_force_n,material);
				if(FreeNorm(candidate.residual,free)<norm){displacement=std::move(trial);assembly=std::move(candidate);accepted=true;break;}
			}catch(const std::runtime_error&){}
		}
		if(!accepted)throw std::runtime_error("native tetrahedral solid static line search failed");
	}
	if(!result.converged)throw std::runtime_error("native tetrahedral solid static Newton did not converge");
	result.displacement_m=displacement;result.reaction_n=assembly.residual;
	result.final_free_residual_n=FreeNorm(assembly.residual,free);result.strain_energy_j=assembly.energy;
	result.minimum_deformation_jacobian=assembly.minimum_jacobian;
	for(std::size_t dof=0;dof<dofs;++dof)result.external_load_potential_j-=external_force_n[dof]*displacement[dof];
	return result;
}

} // namespace iga

#endif
