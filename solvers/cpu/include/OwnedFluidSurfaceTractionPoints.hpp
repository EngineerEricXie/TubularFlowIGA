#ifndef IGA_OWNED_FLUID_SURFACE_TRACTION_POINTS_HPP
#define IGA_OWNED_FLUID_SURFACE_TRACTION_POINTS_HPP

#include "DistributedSurfaceTractionAssembly.hpp"
#include "SharedFluidSurfaceCoefficients.hpp"
#include "FsiDomainRuntime.hpp"

namespace iga {
// Local extraction only. The MPI caller coordinates exceptions, certifies cell
// coverage and state/geometry epoch, and validates shared IGA coefficients.
inline std::vector<SurfaceCellTractionPoints> BuildOwnedFluidSurfaceTractionPoints(
	const CartesianDomainClassification& domain,const ImmersedSurfaceQuadratureCatalog& catalog,
	const MaterialSurfaceKinematics& material,const MaterialSurfacePatchMap& map,
	double viscosity,const std::vector<std::uint64_t>& owned_cells,
	const std::vector<FluidSurfaceElementState>& state)
{
	material.Validate();
	if(!std::isfinite(viscosity)||!(viscosity>0.)
		||material.MaterialIdentitySha256()!=map.FullReference().MaterialIdentitySha256()
		||material.TopologyIdentitySha256()!=map.FullReference().TopologyIdentitySha256()
		||material.Surface().CanonicalSha256()!=domain.SurfaceCanonicalHash()
		||catalog.SurfaceCanonicalHash()!=domain.SurfaceCanonicalHash())
		throw std::invalid_argument("owned fluid traction geometry or viscosity mismatch");
	std::map<std::uint64_t,const FluidSurfaceElementState*> states;
	for(const auto& item:state) {
		if(!states.emplace(item.cell_id,&item).second)throw std::invalid_argument("duplicate fluid traction cell state");
		for(const auto& value:item.nodal_state)for(double component:value)
			if(!std::isfinite(component))throw std::invalid_argument("nonfinite fluid traction coefficient");
	}
	std::map<std::size_t,std::size_t> layout_for_canonical;
	const auto current_mapping = map.LayoutTrianglesByCanonical(material);
	for (std::size_t canonical = 0; canonical < current_mapping.size(); ++canonical)
		if (current_mapping[canonical] != std::numeric_limits<std::size_t>::max())
			layout_for_canonical.emplace(canonical, current_mapping[canonical]);
	std::set<std::uint64_t> seen;
	std::vector<SurfaceCellTractionPoints> result;result.reserve(owned_cells.size());
	std::size_t used_states=0;
	for(auto id:owned_cells) {
		if(id>=domain.Cells().size()||!seen.insert(id).second)throw std::invalid_argument("invalid owned fluid traction cell");
		SurfaceCellTractionPoints cell;cell.cell_id=id;
		if(domain.Cells()[id].classification==CellClassification::Cut) {
			const auto& provenance=catalog.UsableProvenance(domain,id);
			const auto& rule=catalog.UsableRule(domain,id);
			const bool retained=std::any_of(provenance.begin(),provenance.end(),[&](const auto& p){return layout_for_canonical.count(p.canonical_triangle)!=0;});
			if(retained) {
				const auto found=states.find(id);
				if(found==states.end())throw std::invalid_argument("missing retained fluid traction state");
				++used_states;const auto& coefficients=found->second->nodal_state;
				const auto element=domain.Background().MaterializeElement(id);
				if(coefficients.size()!=element.connectivity.size())throw std::invalid_argument("fluid traction coefficient count mismatch");
				for(std::size_t q=0;q<rule.Points().size();++q) {
					const auto& p=rule.Points()[q];const auto& origin=provenance.at(q);
					const auto selected=layout_for_canonical.find(origin.canonical_triangle);
					if(selected==layout_for_canonical.end())continue;
					const auto& triangle=material.CanonicalTriangleProvenance().at(origin.canonical_triangle);
					if(triangle.source_triangle!=map.SourceTriangleForLayoutTriangle(selected->second)
						||triangle.boundary_id!=map.PatchLabel()||p.boundary_id!=static_cast<std::int32_t>(map.PatchLabel()))
						throw std::invalid_argument("owned traction point provenance mismatch");
					SurfaceP1TractionPoint point;point.barycentric=origin.canonical_barycentric;point.weight_m2=p.weight;point.physical_m=p.physical;
					for(int corner=0;corner<3;++corner) {
						const auto source=triangle.source_vertex_indices[triangle.canonical_corner_to_source_corner[corner]];
						bool mapped=false;
						for(auto node:map.Layout().reference_triangles[selected->second])if(map.SourceVertexForGlobalNode(node)==source) {
							point.node_ids[corner]=node;mapped=true;break;
						}
						if(!mapped)throw std::invalid_argument("owned traction point has no material node mapping");
					}
					const auto basis=EvaluateBasis(element,p.parametric[0],p.parametric[1],p.parametric[2],false);
					double pressure=0.;std::array<std::array<double,3>,3> gradient{};
					for(std::size_t a=0;a<coefficients.size();++a) {
						pressure+=coefficients[a][3]*basis.value[a];
						for(int velocity=0;velocity<3;++velocity)for(int direction=0;direction<3;++direction)
							gradient[velocity][direction]+=coefficients[a][velocity]*basis.gradient[a][direction];
					}
					point.traction_pa=FluidOnStructureCauchyTraction(pressure,gradient,viscosity,p.normal);
					cell.points.push_back(point);
				}
			}
		}
		result.push_back(std::move(cell));
	}
	if(used_states!=states.size())throw std::invalid_argument("extra fluid traction state outside retained owned cells");
	return result;
}
// Collective entry point for distributed callers; local extraction remains
// available for callers that already hold an independently verified state.
inline std::vector<SurfaceCellTractionPoints> BuildDistributedFluidSurfaceTractionPoints(
	MPI_Comm comm,const CartesianDomainClassification& domain,const ImmersedSurfaceQuadratureCatalog& catalog,
	const MaterialSurfaceKinematics& material,const MaterialSurfacePatchMap& map,double viscosity,
	const std::vector<std::uint64_t>& owned_cells,const std::vector<FluidSurfaceElementState>& state,
	PointIdentityLimits limits={})
{
	std::string common_identity;
	CollectiveLocalStage(comm,"fluid traction common material binding",[&] {
		material.Validate();
		if(!std::isfinite(viscosity)||!(viscosity>0.))throw std::invalid_argument("invalid fluid traction viscosity");
		Sha256 hash;distributed_surface_detail::AppendString(hash,"DistributedFluidSurfaceTraction/common/v1");
		for(const auto* identity:{&material.ContentIdentitySha256(),&material.MaterialIdentitySha256(),
			&material.TopologyIdentitySha256(),&map.ReferenceIdentitySha256(),&domain.SurfaceCanonicalHash(),&catalog.SurfaceCanonicalHash()})
			distributed_surface_detail::AppendString(hash,*identity);
		hash.AppendNormalizedDouble(viscosity);hash.AppendLittleEndian32(map.PatchLabel());
		hash.AppendLittleEndian64(map.LayoutTriangleToSourceTriangles().size());
		for(auto triangle:map.LayoutTriangleToSourceTriangles())hash.AppendLittleEndian32(triangle);
		const auto& options=catalog.Options();
		hash.AppendLittleEndian64(options.max_candidates);hash.AppendLittleEndian64(options.max_fragments);
		hash.AppendLittleEndian64(options.max_points);hash.AppendLittleEndian64(options.max_exact_limbs);
		common_identity=hash.Hex();
	});
	RequireCollectiveSameText(comm,"fluid traction material and quadrature settings",common_identity);
	ValidateSharedFluidSurfaceCoefficients(comm,domain.Background(),state,limits);
	std::vector<SurfaceCellTractionPoints> result;
	CollectiveLocalStage(comm,"owned fluid surface traction extraction",[&] {
		result=BuildOwnedFluidSurfaceTractionPoints(domain,catalog,material,map,viscosity,owned_cells,state);
	});
	return result;
}
// Expected identities must be captured independently by the transaction
// authority before receiving the state under validation. The fluid identity
// uses the established producer-state/v3 hash with this rank's retained state.
inline std::vector<SurfaceCellTractionPoints> BuildTrialFluidSurfaceTractionPoints(
	MPI_Comm comm,const CartesianDomainClassification& domain,const ImmersedSurfaceQuadratureCatalog& catalog,
	const MaterialSurfaceKinematics& material,const MaterialSurfacePatchMap& map,double viscosity,
	const std::vector<std::uint64_t>& owned_cells,const std::vector<FluidSurfaceElementState>& state,
	const FsiTrialContext& expected_context,const std::string& expected_material_identity,
	const std::string& expected_local_fluid_identity,PointIdentityLimits limits={})
{
	std::string context_identity;
	CollectiveLocalStage(comm,"fluid traction expected trial binding",[&] {
		ValidateFsiTrialContext(expected_context);material.Validate();
		if(!IsLowercaseSha256(expected_material_identity)||!IsLowercaseSha256(expected_local_fluid_identity)
			||material.ContentIdentitySha256()!=expected_material_identity
			||material.EvaluatedTimeS()!=expected_context.EndTime()
			||material.StepStartS()!=expected_context.start_time_s||material.StepEndS()!=expected_context.EndTime())
			throw std::runtime_error("fluid traction material differs from expected trial");
		const auto actual=BuildFluidSurfaceTractionStateIdentitySha256(domain,catalog,material,map,viscosity,state);
		if(actual!=expected_local_fluid_identity)throw std::runtime_error("fluid traction state differs from expected producer identity");
		Sha256 hash;distributed_surface_detail::AppendString(hash,"DistributedFluidSurfaceTraction/trial/v1");
		hash.AppendLittleEndian64(expected_context.step);hash.AppendNormalizedDouble(expected_context.start_time_s);
		hash.AppendNormalizedDouble(expected_context.dt_s);hash.AppendLittleEndian64(expected_context.coupling_iteration);
		distributed_surface_detail::AppendString(hash,expected_material_identity);context_identity=hash.Hex();
	});
	RequireCollectiveSameText(comm,"fluid traction expected common trial",context_identity);
	return BuildDistributedFluidSurfaceTractionPoints(comm,domain,catalog,material,map,viscosity,owned_cells,state,limits);
}
} // namespace iga
#endif
