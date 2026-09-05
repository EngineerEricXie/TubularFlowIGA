#ifndef IGA_IDEALIZED_LEFT_VENTRICLE_FIXTURE_HPP
#define IGA_IDEALIZED_LEFT_VENTRICLE_FIXTURE_HPP

// Deterministic, idealized (not patient-specific) closed chamber used by the
// prescribed-motion regression.  The basal cap is intentionally partitioned
// into two connected half-disks so it exercises mixed flow/pressure ports.
#include "PrescribedSurfaceMotion.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <vector>

namespace iga {

struct IdealizedLeftVentricleFixture {
	static constexpr std::uint32_t AzimuthalSectors = 12;
	static constexpr std::uint32_t AxialRings = 4;
	static constexpr double EndDiastolicLengthM = 0.060;
	static constexpr double EndDiastolicRadiusM = 0.025;
	static constexpr double PeriodS = 0.80;
	static constexpr double EndSystolicVolumeRatio = 0.7695; // .90 * .90 * .95

	static RawSurfaceSoup EndDiastolicSurface()
	{
		RawSurfaceSoup result;
		result.vertices.reserve(1+AxialRings*AzimuthalSectors+1);
		result.vertices.push_back({{0.0,0.0,0.0}}); // apex
		const double pi=std::acos(-1.0);
		for (std::uint32_t ring=1; ring<=AxialRings; ++ring) {
			const double z=EndDiastolicLengthM*static_cast<double>(ring)/AxialRings;
			const double r=EndDiastolicRadiusM*std::sin(0.5*pi*z/EndDiastolicLengthM);
			for (std::uint32_t sector=0; sector<AzimuthalSectors; ++sector) {
				const double theta=2.0*pi*static_cast<double>(sector)/AzimuthalSectors;
				result.vertices.push_back({{r*std::cos(theta),r*std::sin(theta),z}});
			}
		}
		const std::int64_t basal_center=static_cast<std::int64_t>(result.vertices.size());
		result.vertices.push_back({{0.0,0.0,EndDiastolicLengthM}});
		auto ring_id=[](std::uint32_t ring,std::uint32_t sector) {
			return static_cast<std::int64_t>(1+(ring-1)*AzimuthalSectors+sector%AzimuthalSectors);
		};
		for (std::uint32_t sector=0; sector<AzimuthalSectors; ++sector) {
			const auto next=(sector+1)%AzimuthalSectors;
			// Reverse of the natural parameter ordering: outward radial normal.
			result.triangles.push_back({{{0,ring_id(1,next),ring_id(1,sector)}},0});
		}
		for (std::uint32_t ring=1; ring<AxialRings; ++ring)
			for (std::uint32_t sector=0; sector<AzimuthalSectors; ++sector) {
				const auto next=(sector+1)%AzimuthalSectors;
				const auto a=ring_id(ring,sector), b=ring_id(ring,next);
				const auto c=ring_id(ring+1,sector), d=ring_id(ring+1,next);
				result.triangles.push_back({{{a,d,c}},0});
				result.triangles.push_back({{{a,b,d}},0});
			}
		for (std::uint32_t sector=0; sector<AzimuthalSectors; ++sector) {
			const auto next=(sector+1)%AzimuthalSectors;
			const int label=sector<AzimuthalSectors/2 ? 1 : 2;
			result.triangles.push_back({{{basal_center,ring_id(AxialRings,sector),ring_id(AxialRings,next)}},label});
		}
		return result;
	}

	static RawSurfaceSoup EndSystolicSurface()
	{
		auto result=EndDiastolicSurface();
		for (auto& vertex:result.vertices) { vertex[0]*=.90; vertex[1]*=.90; vertex[2]=EndDiastolicLengthM+.95*(vertex[2]-EndDiastolicLengthM); }
		return result;
	}

	static PrescribedSurfaceMotion Motion()
	{
		PrescribedSurfaceMotionOptions options;
		options.require_containment_in_fixed_bounds=true;
		options.fixed_bounds_m.minimum={{-0.032,-0.032,-0.006}};
		options.fixed_bounds_m.maximum={{0.032,0.032,0.066}};
		options.maximum_displacement_m=0.006;
		options.maximum_velocity_m_per_s=0.04;
		options.extension_band_m=0.012;
		return PrescribedSurfaceMotion({{0.0,EndDiastolicSurface()},{0.5*PeriodS,EndSystolicSurface()},{PeriodS,EndDiastolicSurface()}},options);
	}

	static bool BitwiseEqualVertices(const RawSurfaceSoup& a,const RawSurfaceSoup& b)
	{
		if(a.vertices.size()!=b.vertices.size())return false;
		for(std::size_t i=0;i<a.vertices.size();++i)
			for(std::size_t d=0;d<3;++d) if(std::memcmp(&a.vertices[i][d],&b.vertices[i][d],sizeof(double))!=0)return false;
		if(a.triangles.size()!=b.triangles.size())return false;
		for(std::size_t i=0;i<a.triangles.size();++i)
			if(a.triangles[i].indices!=b.triangles[i].indices||a.triangles[i].boundary_id!=b.triangles[i].boundary_id)return false;
		return true;
	}

	static double TriangulatedVolume(const RawSurfaceSoup& soup)
	{
		return ClosedTriangulatedSurface::Build(soup).Diagnostics().volume_m3;
	}
};

} // namespace iga

#endif
