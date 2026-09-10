#ifndef IGA_PARALLEL_BEZIER_VISUALIZATION_HPP
#define IGA_PARALLEL_BEZIER_VISUALIZATION_HPP

#include "BezierVisualization.hpp"
#include "DistributedPointValues.hpp"
#include "ParallelVtkOutput.hpp"

namespace iga {

// Value(node, component) is LOCAL access to already available owned/ghost
// control values. It must not perform MPI. Only supplied owned elements are
// traversed. Cross-element overlap certification remains an upstream duty.
template<class Value>
VtkPartition BuildParallelBezierPartition(MPI_Comm comm,const std::vector<Element>& elements,
	std::uint64_t global_elements,const std::vector<VtkArraySchema>& fields,Value value,
	PointIdentityLimits limits={})
{
	int rank=0,ranks=0;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&ranks);
	std::vector<PointIdentityOccurrence> occurrences;std::vector<std::int64_t> ids;
	std::vector<double> tuples;std::vector<std::uint64_t> counts;
	std::string identity;int components=3;
	std::array<double,3> minimum{{INFINITY,INFINITY,INFINITY}},maximum{{-INFINITY,-INFINITY,-INFINITY}};
	CollectiveLocalStage(comm,"parallel Bezier local extraction",[&] {
		partitioned_vtk_detail::ValidateSchema(fields);
		Sha256 hash;hash.AppendLittleEndian64(global_elements);hash.AppendLittleEndian64(fields.size());
		for(const auto& field:fields) {
			if(field.components>std::numeric_limits<int>::max()-components)throw std::overflow_error("Bezier field component count overflows");
			components+=field.components;hash.AppendLittleEndian64(field.name.size());hash.Append(field.name.data(),field.name.size());hash.AppendLittleEndian32(field.components);
		}
		identity=hash.Hex();counts.resize(ranks);
		if(elements.size()>limits.max_local_occurrences/kBezierPointCount)throw std::runtime_error("Bezier occurrence cap reached");
		const auto points=elements.size()*kBezierPointCount;
		if(points>tuples.max_size()/static_cast<std::size_t>(components))throw std::overflow_error("Bezier tuple count overflows");
		occurrences.reserve(points);ids.reserve(points);tuples.reserve(points*static_cast<std::size_t>(components));
		for(const auto& element:elements) {
			if(element.id>=global_elements||element.owner!=rank)throw std::invalid_argument("invalid owned Bezier element");
			for(std::size_t point=0;point<kBezierPointCount;++point) {
				const auto signature=BuildBezierPointSignature(element,point);
				const auto id=BezierPointOccurrenceId(element.id,point);
				occurrences.push_back({id,EncodeBezierPointSignature(signature)});ids.push_back(id);
				for(unsigned q=0;q<3;++q) {
					const auto coordinate=element.bezier_points[point][q];
					if(!std::isfinite(coordinate))throw std::invalid_argument("nonfinite Bezier coordinate");
					minimum[q]=std::min(minimum[q],coordinate);maximum[q]=std::max(maximum[q],coordinate);tuples.push_back(coordinate);
				}
				for(int component=0;component<components-3;++component) {
					double extracted=0;
					for(const auto& entry:signature.coefficients)extracted+=entry.second*value(entry.first,component);
					if(!std::isfinite(extracted))throw std::invalid_argument("nonfinite Bezier field");
					tuples.push_back(extracted);
				}
			}
		}
	});
	RequireCollectiveSameText(comm,"parallel Bezier schema/domain agreement",identity);
	const auto local_count=static_cast<std::uint64_t>(elements.size());
	MPI_Allgather(&local_count,1,MPI_UINT64_T,counts.data(),1,MPI_UINT64_T,comm);
	CollectiveLocalStage(comm,"parallel Bezier element count",[&] {
		std::uint64_t sum=0;
		for(auto count:counts) {
			if(sum>global_elements||count>global_elements-sum)throw std::runtime_error("Bezier element coverage count exceeds domain");
			sum+=count;
		}
		if(sum!=global_elements)throw std::runtime_error("Bezier element coverage is incomplete");
	});
	// Occurrence uniqueness plus in-range cell ids and exact total cell count
	// proves each global element appears once, without gathering all cell ids.
	const auto representatives=ResolveDistributedPointIdentities(comm,occurrences,limits);
	const auto canonical=ExchangePointRepresentativeValues(comm,ids,representatives,tuples,components,limits);
	std::array<double,3> global_min{},global_max{};
	MPI_Allreduce(minimum.data(),global_min.data(),3,MPI_DOUBLE,MPI_MIN,comm);
	MPI_Allreduce(maximum.data(),global_max.data(),3,MPI_DOUBLE,MPI_MAX,comm);
	VtkPartition result;
	CollectiveLocalStage(comm,"parallel Bezier piece construction",[&] {
		double diagonal=0,max_coordinate=0;
		if(global_elements)for(unsigned q=0;q<3;++q) {
			diagonal=std::hypot(diagonal,global_max[q]-global_min[q]);
			max_coordinate=std::max({max_coordinate,std::abs(global_min[q]),std::abs(global_max[q])});
		}
		const auto tolerance=std::max(1e-11*std::max(diagonal,1e-12),2e-6*std::max(1.,max_coordinate));
		if(!std::isfinite(tolerance))throw std::runtime_error("Bezier coordinate tolerance overflows");
		for(const auto& field:fields)result.point_arrays.push_back({field.name,field.components,{}});
		result.cell_arrays={{"HigherOrderDegrees",3,{}},{"owner",1,{}}};
		std::map<std::int64_t,std::int64_t> local_points;
		std::vector<std::int64_t> connectivity(ids.size());
		for(std::size_t i=0;i<ids.size();++i) {
			const auto offset=i*static_cast<std::size_t>(components);
			double distance=0;
			for(unsigned q=0;q<3;++q)distance=std::hypot(distance,tuples[offset+q]-canonical[offset+q]);
			if(distance>tolerance)throw std::runtime_error("matching Bezier signatures have inconsistent shared coordinates");
			const auto id=representatives[i].occurrence;const auto found=local_points.find(id);
			if(found!=local_points.end()) { connectivity[i]=found->second;continue; }
			const auto local_id=static_cast<std::int64_t>(result.point_ids.size());local_points.emplace(id,local_id);connectivity[i]=local_id;
			result.point_ids.push_back(id);
			for(unsigned q=0;q<3;++q)result.grid.points.push_back(canonical[offset+q]);
			std::size_t field_offset=offset+3;
			for(auto& field:result.point_arrays)for(int q=0;q<field.components;++q)field.values.push_back(canonical[field_offset++]);
		}
		constexpr std::array<double,4> gauss{{.06943184420297371,.33000947820757187,.6699905217924281,.9305681557970262}};
		const auto order=detail::VtkCubicHexTensorIndices();
		for(std::size_t cell=0;cell<elements.size();++cell) {
			Element geometry;std::set<std::int64_t> unique;
			for(std::size_t point=0;point<64;++point) {
				const auto local=connectivity[cell*64+point];
				if(!unique.insert(local).second)throw std::runtime_error("collapsed Bezier element after signature matching");
				for(unsigned q=0;q<3;++q)geometry.bezier_points[point][q]=result.grid.points[3*static_cast<std::size_t>(local)+q];
			}
			for(double u:gauss)for(double v:gauss)for(double w:gauss) {
				std::array<std::array<double,3>,3> jacobian{};detail::MapBezier(geometry,{{u,v,w}},&jacobian);
				const auto determinant=detail::Determinant(jacobian);
				if(!std::isfinite(determinant)||!(determinant>0))throw std::runtime_error("nonpositive parallel Bezier Jacobian");
			}
			for(int point:order)result.grid.connectivity.push_back(connectivity[cell*64+point]);
			result.grid.offsets.push_back(static_cast<std::int64_t>(result.grid.connectivity.size()));result.grid.types.push_back(kVtkBezierHexahedron);
			result.cell_ids.push_back(static_cast<std::int64_t>(elements[cell].id));
			result.cell_arrays[0].values.insert(result.cell_arrays[0].values.end(),{3,3,3});result.cell_arrays[1].values.push_back(rank);
		}
		ValidateVtkPartition(result,0);
	});
	return result;
}
} // namespace iga
#endif
