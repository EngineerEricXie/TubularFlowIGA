#ifndef IGA_NATIVE_TET_ALE_KINEMATICS_HPP
#define IGA_NATIVE_TET_ALE_KINEMATICS_HPP

#include "NativeTetFem.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace iga {

class NativeTetAleFsiRuntime;

inline double NativeAleDeterminant(const std::array<double,3>& a,
	const std::array<double,3>& b,const std::array<double,3>& c,const std::array<double,3>& d)
{
	const std::array<double,3> u{{b[0]-a[0],b[1]-a[1],b[2]-a[2]}};
	const std::array<double,3> v{{c[0]-a[0],c[1]-a[1],c[2]-a[2]}};
	const std::array<double,3> w{{d[0]-a[0],d[1]-a[1],d[2]-a[2]}};
	return u[0]*(v[1]*w[2]-v[2]*w[1])-u[1]*(v[0]*w[2]-v[2]*w[0])
		+u[2]*(v[0]*w[1]-v[1]*w[0]);
}

inline double NativeAleCornerScaledJacobian(const std::array<double,3>& origin,
	const std::array<double,3>& a,const std::array<double,3>& b,const std::array<double,3>& c)
{
	std::array<double,3> u{},v{},w{};
	double nu=0.0,nv=0.0,nw=0.0;
	for(int i=0;i<3;++i) { u[i]=a[i]-origin[i];v[i]=b[i]-origin[i];w[i]=c[i]-origin[i];
		nu+=u[i]*u[i];nv+=v[i]*v[i];nw+=w[i]*w[i]; }
	const double denominator=std::sqrt(nu*nv*nw);
	return denominator>0.0?std::abs(NativeAleDeterminant({{0,0,0}},u,v,w))/denominator:0.0;
}

inline double NativeAleTetScaledJacobian(const std::array<double,3>& a,
	const std::array<double,3>& b,const std::array<double,3>& c,const std::array<double,3>& d)
{
	return std::min({NativeAleCornerScaledJacobian(a,b,c,d),
		NativeAleCornerScaledJacobian(b,a,d,c),NativeAleCornerScaledJacobian(c,d,a,b),
		NativeAleCornerScaledJacobian(d,c,b,a)});
}

struct NativeTetAleQuality
{
	double minimum_determinant_ratio=1.0;
	double minimum_scaled_jacobian=1.0;
};

struct NativeTetAleTrial
{
	double committed_time_s=0.0,next_time_s=0.0,dt_s=0.0;
	std::vector<std::array<double,3>> displacement_m,current_points_m,mesh_velocity_m_s;
	NativeTetAleQuality quality;
};

inline NativeTetMesh BuildNativeTetAleCurrentMesh(const NativeTetMesh& reference_mesh,
	const NativeTetAleTrial& trial)
{
	if (trial.current_points_m.size() != reference_mesh.points.size()
		|| trial.mesh_velocity_m_s.size() != reference_mesh.points.size()
		|| !(trial.dt_s > 0.0))
		throw std::invalid_argument("native ALE trial is incompatible with the reference mesh");
	NativeTetMesh current = reference_mesh;
	current.points = trial.current_points_m;
	return current;
}

class NativeTetAleKinematics
{
public:
	NativeTetAleKinematics(std::vector<std::array<double,3>> reference,
		std::vector<NativeTetCell> cells,double initial_time_s=0.0,
		std::vector<std::array<double,3>> initial_displacement_m={})
		:reference_(std::move(reference)),cells_(std::move(cells)),committed_time_s_(initial_time_s)
	{
		if(reference_.empty()||cells_.empty()||!std::isfinite(initial_time_s))
			throw std::invalid_argument("native ALE reference mesh or time is invalid");
		if(initial_displacement_m.empty())initial_displacement_m.assign(reference_.size(),{{0,0,0}});
		if(initial_displacement_m.size()!=reference_.size())
			throw std::invalid_argument("native ALE initial displacement size mismatch");
		for(const auto& displacement:initial_displacement_m)for(double value:displacement)
			if(!std::isfinite(value))throw std::invalid_argument("native ALE initial displacement is nonfinite");
		committed_displacement_=std::move(initial_displacement_m);
		for(const auto& cell:cells_) {
			for(const auto node:cell.nodes) if(node>=reference_.size())
				throw std::invalid_argument("native ALE cell index is out of range");
			if(!(NativeAleDeterminant(reference_[cell.nodes[0]],reference_[cell.nodes[1]],
				reference_[cell.nodes[2]],reference_[cell.nodes[3]])>0.0))
				throw std::invalid_argument("native ALE reference tetrahedron is inverted");
			std::array<std::array<double,3>,4> current;
			for(std::size_t local=0;local<4;++local)for(int component=0;component<3;++component)
				current[local][component]=reference_[cell.nodes[local]][component]
					+committed_displacement_[cell.nodes[local]][component];
			if(!(NativeAleDeterminant(current[0],current[1],current[2],current[3])>0.0))
				throw std::invalid_argument("native ALE initial displacement inverts a tetrahedron");
		}
	}

	const NativeTetAleTrial& BeginTrial(const std::vector<std::array<double,3>>& displacement,
		double next_time_s,double minimum_ratio=0.05,double minimum_scaled_jacobian=1e-3)
	{
		if(has_prepared_||displacement.size()!=reference_.size()||!std::isfinite(next_time_s)
			||!(next_time_s>committed_time_s_)||!(minimum_ratio>0.0)
			||!(minimum_scaled_jacobian>0.0))
			throw std::invalid_argument("native ALE trial input is invalid");
		has_trial_=false;trial_=NativeTetAleTrial{};
		NativeTetAleTrial candidate;candidate.committed_time_s=committed_time_s_;
		candidate.next_time_s=next_time_s;candidate.dt_s=next_time_s-committed_time_s_;
		candidate.displacement_m=displacement;candidate.current_points_m.resize(reference_.size());
		candidate.mesh_velocity_m_s.resize(reference_.size());
		for(std::size_t node=0;node<reference_.size();++node)
			for(int component=0;component<3;++component) {
				if(!std::isfinite(displacement[node][component]))
					throw std::invalid_argument("native ALE displacement is nonfinite");
				candidate.current_points_m[node][component]=reference_[node][component]
					+displacement[node][component];
				candidate.mesh_velocity_m_s[node][component]=(displacement[node][component]
					-committed_displacement_[node][component])/candidate.dt_s;
			}
		candidate.quality.minimum_determinant_ratio=std::numeric_limits<double>::infinity();
		candidate.quality.minimum_scaled_jacobian=std::numeric_limits<double>::infinity();
		for(const auto& cell:cells_) {
			const double reference_det=NativeAleDeterminant(reference_[cell.nodes[0]],
				reference_[cell.nodes[1]],reference_[cell.nodes[2]],reference_[cell.nodes[3]]);
			const auto& x0=candidate.current_points_m[cell.nodes[0]];
			const auto& x1=candidate.current_points_m[cell.nodes[1]];
			const auto& x2=candidate.current_points_m[cell.nodes[2]];
			const auto& x3=candidate.current_points_m[cell.nodes[3]];
			const double current_det=NativeAleDeterminant(x0,x1,x2,x3);
			if(!(current_det>0.0)||!std::isfinite(current_det))
				throw std::runtime_error("native ALE trial has a nonpositive current determinant");
			candidate.quality.minimum_determinant_ratio=std::min(
				candidate.quality.minimum_determinant_ratio,current_det/reference_det);
			candidate.quality.minimum_scaled_jacobian=std::min(
				candidate.quality.minimum_scaled_jacobian,NativeAleTetScaledJacobian(x0,x1,x2,x3));
		}
		if(candidate.quality.minimum_determinant_ratio<minimum_ratio
			||candidate.quality.minimum_scaled_jacobian<minimum_scaled_jacobian)
			throw std::runtime_error("native ALE trial failed a mesh-quality gate");
		trial_=std::move(candidate);has_trial_=true;return trial_;
	}

	void CommitTrial()
	{
		PrepareTrialCommit();
		FinalizePreparedTrialNoexcept();
	}
	void PrepareTrialCommit()
	{
		if(!has_trial_||has_prepared_)
			throw std::logic_error("native ALE prepare requires one unprepared trial");
		prepared_displacement_=trial_.displacement_m;
		prepared_time_s_=trial_.next_time_s;has_prepared_=true;
	}
	void RequireFinalizePreparedTrialAllowed() const
	{
		if(!has_trial_||!has_prepared_)
			throw std::logic_error("native ALE finalize requires a prepared trial");
	}
	void RejectTrial() noexcept { has_prepared_=false;has_trial_=false; }
	double CommittedTime() const { return committed_time_s_; }
	const std::vector<std::array<double,3>>& CommittedDisplacement() const
	{
		return committed_displacement_;
	}
	bool HasTrial() const { return has_trial_; }

private:
	friend class NativeTetAleFsiRuntime;
	void FinalizePreparedTrialNoexcept() noexcept
	{
		using std::swap;swap(committed_displacement_,prepared_displacement_);
		committed_time_s_=prepared_time_s_;has_prepared_=false;has_trial_=false;
	}
	std::vector<std::array<double,3>> reference_,committed_displacement_,prepared_displacement_;
	std::vector<NativeTetCell> cells_;
	double committed_time_s_=0.0;
	NativeTetAleTrial trial_;
	double prepared_time_s_=0.0;
	bool has_trial_=false,has_prepared_=false;
};

} // namespace iga

#endif
