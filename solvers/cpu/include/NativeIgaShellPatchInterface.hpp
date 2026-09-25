#ifndef IGA_NATIVE_IGA_SHELL_PATCH_INTERFACE_HPP
#define IGA_NATIVE_IGA_SHELL_PATCH_INTERFACE_HPP

#include "NativeIgaKirchhoffLove.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace iga {

struct NativeIgaShellInterfaceDiagnostics
{
	double maximum_reference_gap_m=0.0;
	double maximum_current_gap_m=0.0;
	double maximum_reference_tangent_jump=0.0;
	double maximum_current_tangent_jump=0.0;
	double minimum_reference_normal_dot=1.0;
	double minimum_current_normal_dot=1.0;
	double maximum_membrane_resultant_jump_n_m=0.0;
	double maximum_bending_moment_jump_n=0.0;
};

struct NativeIgaShellLinearConstraint
{
	std::vector<std::pair<std::size_t,double>> terms;
	double right_hand_side=0.0;
};

inline std::vector<NativeIgaShellLinearConstraint> BuildNativeIgaShellCompatibleUConstraints(
	const NativeNurbsShellPatch& left,const NativeNurbsShellPatch& right)
{
	ValidateNativeNurbsShellPatch(left);ValidateNativeNurbsShellPatch(right);
	if(left.controls_v!=right.controls_v||left.degree_v!=right.degree_v
		||left.knots_v!=right.knots_v)
		throw std::invalid_argument("native shell compatible interface requires matching seam bases");
	if(left.controls_u<2||right.controls_u<2)
		throw std::invalid_argument("native shell compatible interface lacks adjacent control rows");
	const double left_span=left.knots_u[left.controls_u]-left.knots_u[left.controls_u-1];
	const double right_span=right.knots_u[static_cast<std::size_t>(right.degree_u+1)]
		-right.knots_u[static_cast<std::size_t>(right.degree_u)];
	if(!(left_span>0.0)||!(right_span>0.0))
		throw std::invalid_argument("native shell compatible interface has a zero end span");
	const double left_scale=left.degree_u/left_span,right_scale=right.degree_u/right_span;
	const std::size_t right_offset=3*left.control_points.size();
	std::vector<NativeIgaShellLinearConstraint> result;
	result.reserve(6*left.controls_v);
	for(std::size_t j=0;j<left.controls_v;++j)
		for(std::size_t component=0;component<3;++component) {
			const std::size_t left_seam=3*(j*left.controls_u+left.controls_u-1)+component;
			const std::size_t left_inner=3*(j*left.controls_u+left.controls_u-2)+component;
			const std::size_t right_seam=right_offset+3*(j*right.controls_u)+component;
			const std::size_t right_inner=right_offset+3*(j*right.controls_u+1)+component;
			result.push_back({{{left_seam,-1.0},{right_seam,1.0}},0.0});
			result.push_back({{{left_seam,-left_scale},{left_inner,left_scale},
				{right_seam,-right_scale},{right_inner,right_scale}},0.0});
		}
	return result;
}

inline double EvaluateNativeIgaShellConstraint(const NativeIgaShellLinearConstraint& constraint,
	const std::vector<double>& values)
{
	double result=-constraint.right_hand_side;
	for(const auto& term:constraint.terms) {
		if(term.first>=values.size()) throw std::invalid_argument("native shell constraint index is out of range");
		result+=term.second*values[term.first];
	}
	return result;
}

struct NativeIgaShellInterfaceTransformation
{
	std::size_t local_dofs=0;
	std::size_t independent_dofs=0;
	// Row-major local = matrix * independent.
	std::vector<double> matrix;
};

inline NativeIgaShellInterfaceTransformation BuildNativeIgaShellCompatibleUTransformation(
	const NativeNurbsShellPatch& left,const NativeNurbsShellPatch& right)
{
	// Reuse all structural checks and obtain the end-span scaling.
	BuildNativeIgaShellCompatibleUConstraints(left,right);
	const std::size_t left_dofs=3*left.control_points.size();
	const std::size_t right_dofs=3*right.control_points.size();
	const std::size_t eliminated=6*right.controls_v;
	NativeIgaShellInterfaceTransformation result;
	result.local_dofs=left_dofs+right_dofs;
	result.independent_dofs=result.local_dofs-eliminated;
	result.matrix.assign(result.local_dofs*result.independent_dofs,0.0);
	for(std::size_t dof=0;dof<left_dofs;++dof)
		result.matrix[dof*result.independent_dofs+dof]=1.0;
	const double left_span=left.knots_u[left.controls_u]-left.knots_u[left.controls_u-1];
	const double right_span=right.knots_u[static_cast<std::size_t>(right.degree_u+1)]
		-right.knots_u[static_cast<std::size_t>(right.degree_u)];
	const double ratio=(left.degree_u/left_span)/(right.degree_u/right_span);
	std::size_t next=left_dofs;
	for(std::size_t j=0;j<right.controls_v;++j)
		for(std::size_t i=0;i<right.controls_u;++i)
			for(std::size_t component=0;component<3;++component) {
				const std::size_t row=left_dofs+3*(j*right.controls_u+i)+component;
				const std::size_t left_seam=3*(j*left.controls_u+left.controls_u-1)+component;
				const std::size_t left_inner=3*(j*left.controls_u+left.controls_u-2)+component;
				if(i==0) result.matrix[row*result.independent_dofs+left_seam]=1.0;
				else if(i==1) {
					result.matrix[row*result.independent_dofs+left_seam]=1.0+ratio;
					result.matrix[row*result.independent_dofs+left_inner]=-ratio;
				} else result.matrix[row*result.independent_dofs+next++]=1.0;
			}
	if(next!=result.independent_dofs)
		throw std::logic_error("native shell interface transformation size mismatch");
	return result;
}

inline std::vector<double> ApplyNativeIgaShellInterfaceTransformation(
	const NativeIgaShellInterfaceTransformation& transformation,
	const std::vector<double>& independent)
{
	if(independent.size()!=transformation.independent_dofs)
		throw std::invalid_argument("native shell independent interface vector size mismatch");
	std::vector<double> local(transformation.local_dofs,0.0);
	for(std::size_t row=0;row<transformation.local_dofs;++row)
		for(std::size_t column=0;column<transformation.independent_dofs;++column)
			local[row]+=transformation.matrix[row*transformation.independent_dofs+column]
				*independent[column];
	return local;
}

inline std::vector<double> ReduceNativeIgaShellInterfaceResidual(
	const NativeIgaShellInterfaceTransformation& transformation,
	const std::vector<double>& local_residual)
{
	if(local_residual.size()!=transformation.local_dofs)
		throw std::invalid_argument("native shell local residual size mismatch");
	std::vector<double> result(transformation.independent_dofs,0.0);
	for(std::size_t row=0;row<transformation.local_dofs;++row)
		for(std::size_t column=0;column<transformation.independent_dofs;++column)
			result[column]+=transformation.matrix[row*transformation.independent_dofs+column]
				*local_residual[row];
	return result;
}

inline std::vector<double> ReduceNativeIgaShellInterfaceTangent(
	const NativeIgaShellInterfaceTransformation& transformation,
	const std::vector<double>& local_tangent)
{
	if(local_tangent.size()!=transformation.local_dofs*transformation.local_dofs)
		throw std::invalid_argument("native shell local tangent size mismatch");
	const std::size_t local=transformation.local_dofs,global=transformation.independent_dofs;
	std::vector<double> temporary(local*global,0.0),result(global*global,0.0);
	for(std::size_t i=0;i<local;++i) for(std::size_t j=0;j<local;++j)
		for(std::size_t b=0;b<global;++b)
			temporary[i*global+b]+=local_tangent[i*local+j]
				*transformation.matrix[j*global+b];
	for(std::size_t a=0;a<global;++a) for(std::size_t i=0;i<local;++i)
		for(std::size_t b=0;b<global;++b)
			result[a*global+b]+=transformation.matrix[i*global+a]*temporary[i*global+b];
	return result;
}

inline double NativeShellVectorDistance(const std::array<double,3>& a,
	const std::array<double,3>& b)
{
	double sum=0.0;for(int i=0;i<3;++i) { const double d=a[i]-b[i];sum+=d*d; }
	return std::sqrt(sum);
}

inline double NativeShellRelativeVectorJump(const std::array<double,3>& a,
	const std::array<double,3>& b)
{
	const double difference=NativeShellVectorDistance(a,b);
	return difference/std::max({1.0e-30,std::sqrt(NativeShellDot(a,a)),std::sqrt(NativeShellDot(b,b))});
}

inline NativeIgaShellInterfaceDiagnostics CheckNativeIgaShellCompatibleUInterface(
	const NativeNurbsShellPatch& left_reference,const NativeNurbsShellPatch& right_reference,
	const NativeNurbsShellPatch& left_current,const NativeNurbsShellPatch& right_current,
	const NativeKirchhoffLoveMaterial& material,double geometry_tolerance_m=1e-10,
	double tangent_tolerance=1e-9,double resultant_tolerance=1e-7)
{
	for(const auto* patch:{&left_reference,&right_reference,&left_current,&right_current})
		ValidateNativeNurbsShellPatch(*patch);
	if(!(geometry_tolerance_m>0.0)||!(tangent_tolerance>0.0)||!(resultant_tolerance>0.0))
		throw std::invalid_argument("native shell interface tolerances must be positive");
	const double left_u=left_reference.knots_u[left_reference.controls_u];
	const double right_u=right_reference.knots_u[static_cast<std::size_t>(right_reference.degree_u)];
	const double left_current_u=left_current.knots_u[left_current.controls_u];
	const double right_current_u=right_current.knots_u[static_cast<std::size_t>(right_current.degree_u)];
	const double v0=left_reference.knots_v[static_cast<std::size_t>(left_reference.degree_v)];
	const double v1=left_reference.knots_v[left_reference.controls_v];
	if(right_reference.knots_v[static_cast<std::size_t>(right_reference.degree_v)]!=v0
		||right_reference.knots_v[right_reference.controls_v]!=v1
		||left_current.knots_v[static_cast<std::size_t>(left_current.degree_v)]!=v0
		||right_current.knots_v[static_cast<std::size_t>(right_current.degree_v)]!=v0
		||left_current.knots_v[left_current.controls_v]!=v1
		||right_current.knots_v[right_current.controls_v]!=v1)
		throw std::invalid_argument("native shell interface v parameter domains differ");
	NativeIgaShellInterfaceDiagnostics result;
	for(const double fraction:{0.0,0.1127016653792583,0.5,0.8872983346207417,1.0}) {
		const double v=v0+fraction*(v1-v0);
		const auto lr=EvaluateNativeNurbsShellSurface(left_reference,left_u,v);
		const auto rr=EvaluateNativeNurbsShellSurface(right_reference,right_u,v);
		const auto lc=EvaluateNativeNurbsShellSurface(left_current,left_current_u,v);
		const auto rc=EvaluateNativeNurbsShellSurface(right_current,right_current_u,v);
		result.maximum_reference_gap_m=std::max(result.maximum_reference_gap_m,
			NativeShellVectorDistance(lr.position,rr.position));
		result.maximum_current_gap_m=std::max(result.maximum_current_gap_m,
			NativeShellVectorDistance(lc.position,rc.position));
		result.maximum_reference_tangent_jump=std::max(result.maximum_reference_tangent_jump,
			std::max(NativeShellRelativeVectorJump(lr.tangent_u,rr.tangent_u),
				NativeShellRelativeVectorJump(lr.tangent_v,rr.tangent_v)));
		result.maximum_current_tangent_jump=std::max(result.maximum_current_tangent_jump,
			std::max(NativeShellRelativeVectorJump(lc.tangent_u,rc.tangent_u),
				NativeShellRelativeVectorJump(lc.tangent_v,rc.tangent_v)));
		result.minimum_reference_normal_dot=std::min(result.minimum_reference_normal_dot,
			NativeShellDot(lr.unit_normal,rr.unit_normal));
		result.minimum_current_normal_dot=std::min(result.minimum_current_normal_dot,
			NativeShellDot(lc.unit_normal,rc.unit_normal));
		const auto left_response=EvaluateNativeKirchhoffLovePoint(lr,lc,material);
		const auto right_response=EvaluateNativeKirchhoffLovePoint(rr,rc,material);
		for(int component=0;component<2;++component) {
			result.maximum_membrane_resultant_jump_n_m=std::max(
				result.maximum_membrane_resultant_jump_n_m,std::abs(
					left_response.membrane_resultant_contravariant[0][component]
					-right_response.membrane_resultant_contravariant[0][component]));
			result.maximum_bending_moment_jump_n=std::max(
				result.maximum_bending_moment_jump_n,std::abs(
					left_response.bending_moment_contravariant[0][component]
					-right_response.bending_moment_contravariant[0][component]));
		}
	}
	if(result.maximum_reference_gap_m>geometry_tolerance_m
		||result.maximum_current_gap_m>geometry_tolerance_m)
		throw std::runtime_error("native shell interface has a geometric gap");
	if(result.maximum_reference_tangent_jump>tangent_tolerance
		||result.maximum_current_tangent_jump>tangent_tolerance
		||result.minimum_reference_normal_dot<1.0-tangent_tolerance
		||result.minimum_current_normal_dot<1.0-tangent_tolerance)
		throw std::runtime_error("native shell interface is not C1-compatible");
	if(result.maximum_membrane_resultant_jump_n_m>resultant_tolerance
		||result.maximum_bending_moment_jump_n>resultant_tolerance)
		throw std::runtime_error("native shell interface force or moment is discontinuous");
	return result;
}

} // namespace iga

#endif
