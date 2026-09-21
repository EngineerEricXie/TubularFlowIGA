#ifndef IGA_NATIVE_TET_MOVING_SPECIES_TRANSPORT_HPP
#define IGA_NATIVE_TET_MOVING_SPECIES_TRANSPORT_HPP

#include "NativeTetFem.hpp"
#include "NativeTetRt0Flux.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

struct NativeTetWallExchange
{
	double transfer_coefficient_m_s=0.;
	double external_concentration_mol_m3=0.;
};

// Concentration is mol/m^3. Boundary flux and source budgets are mol/s.
// The conservative backward-Euler inventory uses each step's own geometry.
struct NativeTetMovingSpeciesStep
{
	std::vector<double> concentration_mol_m3;
	double previous_inventory_mol=0.;
	double current_inventory_mol=0.;
	double outward_advective_flux_mol_s=0.;
	std::map<int,double> outward_advective_flux_by_label_mol_s;
	std::map<int,double> outward_relative_flow_by_label_m3_s;
	std::map<int,double> positive_relative_flow_by_label_m3_s;
	std::map<int,double> negative_relative_flow_by_label_m3_s;
	double source_mol_s=0.;
	double reaction_sink_mol_s=0.;
	double outward_wall_exchange_mol_s=0.;
	std::map<int,double> outward_wall_exchange_by_label_mol_s;
	double balance_defect_mol_s=0.;
};

struct NativeTetMovingSpeciesAssembly
{
	std::size_t nodes=0;
	std::map<std::pair<std::uint32_t,std::uint32_t>,double> matrix;
	std::vector<double> rhs;
	std::vector<double> outward_coefficients;
	std::map<int,std::vector<double>> outward_coefficients_by_label;
	double incoming_flux_mol_s=0.;
	std::map<int,double> incoming_flux_by_label_mol_s;
	std::map<int,double> outward_relative_flow_by_label_m3_s;
	std::map<int,double> positive_relative_flow_by_label_m3_s;
	std::map<int,double> negative_relative_flow_by_label_m3_s;
	double source_mol_s=0.;
	double first_order_decay_rate_s_inv=0.;
	std::map<int,std::vector<double>> wall_exchange_coefficients_by_label;
	std::map<int,double> wall_exchange_constant_by_label_mol_s;
};

namespace native_tet_moving_species_detail {

using Face=std::array<std::uint32_t,3>;

using SparseMatrix=std::map<std::pair<std::uint32_t,std::uint32_t>,double>;

inline std::array<double,16> ElementBlock(const SparseMatrix& matrix,
	const std::array<std::uint32_t,4>& nodes)
{
	std::array<double,16> block{};
	for(std::size_t row=0;row<4;++row)
		for(std::size_t column=0;column<4;++column){
			const auto entry=matrix.find({nodes[row],nodes[column]});
			if(entry!=matrix.end())block[4*row+column]=entry->second;
		}
	return block;
}

// A symmetric graph Laplacian removes positive off-diagonal entries from each
// owned cell/face contribution. Its row and column sums vanish, preserving
// constant states and the assembled inventory balance across MPI partitions.
inline void AddMonotoneGraphDiffusion(SparseMatrix& matrix,
	const std::array<std::uint32_t,4>& nodes,
	const std::array<double,16>& before)
{
	const auto after=ElementBlock(matrix,nodes);
	for(std::size_t row=0;row<4;++row)
		for(std::size_t column=row+1;column<4;++column){
			const double forward=after[4*row+column]-before[4*row+column];
			const double reverse=after[4*column+row]-before[4*column+row];
			const double diffusion=std::max({0.,forward,reverse});
			matrix[{nodes[row],nodes[column]}]-=diffusion;
			matrix[{nodes[column],nodes[row]}]-=diffusion;
			matrix[{nodes[row],nodes[row]}]+=diffusion;
			matrix[{nodes[column],nodes[column]}]+=diffusion;
		}
}

struct FaceOwner
{
	std::size_t cell=0;
	std::size_t opposite=0;
	unsigned uses=0;
};

inline std::map<Face,FaceOwner> BoundaryOwners(const NativeTetMesh& mesh)
{
	std::map<Face,FaceOwner> owners;
	for(std::size_t cell=0;cell<mesh.cells.size();++cell)
		for(std::size_t opposite=0;opposite<4;++opposite){
			Face face{};std::size_t index=0;
			for(std::size_t local=0;local<4;++local)
				if(local!=opposite)face[index++]=mesh.cells[cell].nodes[local];
			std::sort(face.begin(),face.end());
			auto& owner=owners[face];
			if(owner.uses++==0){owner.cell=cell;owner.opposite=opposite;}
			if(owner.uses>2)
				throw std::invalid_argument("native moving species mesh is nonmanifold");
		}
	return owners;
}

inline std::vector<double> SolveDense(std::vector<double> matrix,
	std::vector<double> rhs)
{
	const std::size_t size=rhs.size();
	if(matrix.size()!=size*size||size==0)
		throw std::invalid_argument("native moving species dense system size is invalid");
	double largest=0.;
	for(const double value:matrix)largest=std::max(largest,std::abs(value));
	for(std::size_t column=0;column<size;++column){
		std::size_t pivot=column;
		for(std::size_t row=column+1;row<size;++row)
			if(std::abs(matrix[row*size+column])>
				std::abs(matrix[pivot*size+column]))pivot=row;
		if(!(std::abs(matrix[pivot*size+column])>1e-14*largest))
			throw std::runtime_error("native moving species dense system is singular");
		if(pivot!=column){
			for(std::size_t entry=column;entry<size;++entry)
				std::swap(matrix[column*size+entry],matrix[pivot*size+entry]);
			std::swap(rhs[column],rhs[pivot]);
		}
		for(std::size_t row=column+1;row<size;++row){
			const double factor=matrix[row*size+column]/matrix[column*size+column];
			matrix[row*size+column]=0.;
			for(std::size_t entry=column+1;entry<size;++entry)
				matrix[row*size+entry]-=factor*matrix[column*size+entry];
			rhs[row]-=factor*rhs[column];
		}
	}
	std::vector<double> result(size,0.);
	for(std::size_t reverse=size;reverse>0;--reverse){
		const std::size_t row=reverse-1;double value=rhs[row];
		for(std::size_t column=row+1;column<size;++column)
			value-=matrix[row*size+column]*result[column];
		result[row]=value/matrix[row*size+row];
		if(!std::isfinite(result[row]))
			throw std::runtime_error("native moving species solution is nonfinite");
	}
	return result;
}

inline double Inventory(const NativeTetMesh& mesh,const std::vector<double>& concentration)
{
	double inventory=0.;
	for(const auto& cell:mesh.cells){
		const double volume=EvaluateNativeTetGeometry(mesh,cell).determinant/6.;
		for(const auto node:cell.nodes)inventory+=volume*concentration[node]/4.;
	}
	return inventory;
}

} // namespace native_tet_moving_species_detail

// Native weak-form assembly. Cells and boundary faces are partitioned by
// owner rank; each caller contributes only its own sparse algebra entries.
inline NativeTetMovingSpeciesAssembly AssembleNativeTetMovingSpeciesStep(
	const NativeTetMesh& previous_mesh,const NativeTetMesh& current_mesh,
	const std::vector<std::array<double,3>>& fluid_velocity_nodes_m_s,
	const std::vector<std::array<double,3>>& mesh_velocity_nodes_m_s,
	const std::vector<double>& previous_concentration_mol_m3,
	const std::map<int,double>& inflow_concentration_mol_m3,
	double diffusivity_m2_s,double source_mol_m3_s,double dt_s,
	int rank=0,int ranks=1,bool monotone=false,
	double first_order_decay_rate_s_inv=0.,
	const std::map<int,NativeTetWallExchange>& wall_exchange={},
	const std::vector<std::array<double,4>>& rt0_face_flow_m3_s={},
	const std::vector<double>& cell_source_mol_m3_s={})
{
	const std::size_t nodes=current_mesh.points.size();
	if(nodes==0||previous_mesh.points.size()!=nodes
		||previous_mesh.cells.size()!=current_mesh.cells.size()
		||previous_mesh.boundary_triangles.size()!=current_mesh.boundary_triangles.size()
		||rank<0||ranks<=0||rank>=ranks
		||mesh_velocity_nodes_m_s.size()!=nodes
		||previous_concentration_mol_m3.size()!=nodes
		||!(dt_s>0.)||!std::isfinite(dt_s)
		||!(diffusivity_m2_s>=0.)||!std::isfinite(diffusivity_m2_s)
		||!std::isfinite(source_mol_m3_s)
		||(!cell_source_mol_m3_s.empty()
			&&cell_source_mol_m3_s.size()!=current_mesh.cells.size())
		||!(first_order_decay_rate_s_inv>=0.)
		||!std::isfinite(first_order_decay_rate_s_inv))
		throw std::invalid_argument("native moving species assembly input is invalid");
	for(const double source:cell_source_mol_m3_s)
		if(!std::isfinite(source))
			throw std::invalid_argument("native moving species cell source is nonfinite");
	for(std::size_t cell=0;cell<current_mesh.cells.size();++cell)
		if(previous_mesh.cells[cell].nodes!=current_mesh.cells[cell].nodes)
			throw std::invalid_argument("native moving species topology changed during the step");
	for(std::size_t cell=0;cell<current_mesh.cells.size();++cell){
		EvaluateNativeTetGeometry(previous_mesh,previous_mesh.cells[cell]);
		EvaluateNativeTetGeometry(current_mesh,current_mesh.cells[cell]);
	}
	for(std::size_t face=0;face<current_mesh.boundary_triangles.size();++face)
		if(previous_mesh.boundary_triangles[face].nodes
			!=current_mesh.boundary_triangles[face].nodes
			||previous_mesh.boundary_triangles[face].boundary_label
				!=current_mesh.boundary_triangles[face].boundary_label)
			throw std::invalid_argument("native moving species boundary changed during the step");
	for(const double concentration:previous_concentration_mol_m3)
		if(!std::isfinite(concentration)||concentration<0.)
			throw std::invalid_argument("native moving species previous concentration is invalid");
	for(std::size_t node=0;node<nodes;++node)
		for(int component=0;component<3;++component){
			const double velocity=mesh_velocity_nodes_m_s[node][component];
			const double expected=(current_mesh.points[node][component]
				-previous_mesh.points[node][component])/dt_s;
			const double tolerance=256.*std::numeric_limits<double>::epsilon()
				*std::max({1.,std::abs(velocity),std::abs(expected)});
			if(!std::isfinite(velocity)||!std::isfinite(expected)
				||std::abs(velocity-expected)>tolerance)
				throw std::invalid_argument("native moving species mesh velocity disagrees with geometry");
		}
	for(const auto& inlet:inflow_concentration_mol_m3)
		if(!std::isfinite(inlet.second)||inlet.second<0.)
			throw std::invalid_argument("native moving species inflow concentration is invalid");
	for(const auto& item:wall_exchange)
		if(item.first<0||!(item.second.transfer_coefficient_m_s>=0.)
			||!std::isfinite(item.second.transfer_coefficient_m_s)
			||!(item.second.external_concentration_mol_m3>=0.)
			||!std::isfinite(item.second.external_concentration_mol_m3))
			throw std::invalid_argument("native moving species wall exchange is invalid");
	const auto topology=BuildNativeTaylorHoodTopology(current_mesh);
	const bool use_rt0=!rt0_face_flow_m3_s.empty();
	if((use_rt0&&(!fluid_velocity_nodes_m_s.empty()
		||rt0_face_flow_m3_s.size()!=current_mesh.cells.size()))
		||(!use_rt0&&fluid_velocity_nodes_m_s.size()!=nodes+topology.edges.size()))
		throw std::invalid_argument("native moving species P2 velocity size is invalid");
	for(const auto& velocity:fluid_velocity_nodes_m_s)
		for(const double component:velocity)
			if(!std::isfinite(component))
				throw std::invalid_argument("native moving species fluid velocity is nonfinite");
	double global_speed_bound=0.;
	for(const auto& velocity:fluid_velocity_nodes_m_s)
		global_speed_bound=std::max(global_speed_bound,std::hypot(
			velocity[0],velocity[1],velocity[2]));
	if(use_rt0)
		for(std::size_t cell_index=0;cell_index<current_mesh.cells.size();++cell_index){
			for(std::size_t local=0;local<4;++local){
				const auto value=EvaluateNativeTetRt0Flux(current_mesh,
					rt0_face_flow_m3_s,cell_index,
					current_mesh.points.at(current_mesh.cells[cell_index].nodes[local]));
				global_speed_bound=std::max(global_speed_bound,std::hypot(
					value[0],value[1],value[2]));
			}
		}
	for(const auto& velocity:mesh_velocity_nodes_m_s)
		global_speed_bound=std::max(global_speed_bound,std::hypot(
			velocity[0],velocity[1],velocity[2]));
	NativeTetMovingSpeciesAssembly assembly;
	assembly.nodes=nodes;
	assembly.first_order_decay_rate_s_inv=first_order_decay_rate_s_inv;
	assembly.rhs.assign(nodes,0.);
	assembly.outward_coefficients.assign(nodes,0.);
	for(const auto& triangle:current_mesh.boundary_triangles){
		assembly.outward_coefficients_by_label.try_emplace(
			triangle.boundary_label,nodes,0.);
		assembly.incoming_flux_by_label_mol_s.try_emplace(
			triangle.boundary_label,0.);
		assembly.outward_relative_flow_by_label_m3_s.try_emplace(
			triangle.boundary_label,0.);
		assembly.positive_relative_flow_by_label_m3_s.try_emplace(
			triangle.boundary_label,0.);
		assembly.negative_relative_flow_by_label_m3_s.try_emplace(
			triangle.boundary_label,0.);
	}
	for(const auto& item:wall_exchange){
		if(!assembly.outward_coefficients_by_label.count(item.first))
			throw std::invalid_argument("native moving species wall exchange label is missing");
		assembly.wall_exchange_coefficients_by_label.emplace(item.first,
			std::vector<double>(nodes,0.));
		assembly.wall_exchange_constant_by_label_mol_s.emplace(item.first,0.);
	}
	auto& matrix=assembly.matrix;
	auto& rhs=assembly.rhs;
	auto& outward_coefficients=assembly.outward_coefficients;
	auto& incoming_flux_mol_s=assembly.incoming_flux_mol_s;
	auto& source_mol_s=assembly.source_mol_s;
	for(std::size_t cell_index=0;cell_index<current_mesh.cells.size();++cell_index){
		if(cell_index%static_cast<std::size_t>(ranks)
			!=static_cast<std::size_t>(rank))continue;
		const auto& cell=current_mesh.cells[cell_index];
		const auto before=monotone
			?native_tet_moving_species_detail::ElementBlock(matrix,cell.nodes)
			:std::array<double,16>{};
		const double previous_volume=EvaluateNativeTetGeometry(previous_mesh,
			previous_mesh.cells[cell_index]).determinant/6.;
		const auto geometry=EvaluateNativeTetGeometry(current_mesh,cell);
		const double current_volume=geometry.determinant/6.;
		const double cell_source=source_mol_m3_s+
			(cell_source_mol_m3_s.empty()?0.:cell_source_mol_m3_s[cell_index]);
		source_mol_s+=cell_source*current_volume;
		for(std::size_t row=0;row<4;++row){
			const auto global_row=cell.nodes[row];
			rhs[global_row]+=cell_source*current_volume/4.;
			for(std::size_t column=0;column<4;++column){
				const auto global_column=cell.nodes[column];
				const double factor=row==column?2.:1.;
				matrix[{global_row,global_column}]
					+=current_volume*factor/(20.*dt_s);
				matrix[{global_row,global_column}]
					+=first_order_decay_rate_s_inv*current_volume*factor/20.;
				rhs[global_row]+=previous_volume*factor/(20.*dt_s)
					*previous_concentration_mol_m3[global_column];
				const auto& a=geometry.barycentric_gradients[row];
				const auto& b=geometry.barycentric_gradients[column];
				matrix[{global_row,global_column}]+=diffusivity_m2_s
					*current_volume*(a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
			}
		}
		for(const auto& point:NativeTetDegreeFiveQuadrature()){
			const auto basis=EvaluateNativeTaylorHoodBasis(point.reference[0],
				point.reference[1],point.reference[2]);
			std::array<double,3> relative_velocity{};
			if(use_rt0){
				std::array<double,3> position{};
				for(std::size_t local=0;local<4;++local)
					for(int component=0;component<3;++component)
						position[component]+=basis.pressure[local]
							*current_mesh.points.at(cell.nodes[local])[component];
				relative_velocity=EvaluateNativeTetRt0Flux(current_mesh,
					rt0_face_flow_m3_s,cell_index,position);
			}else for(std::size_t local=0;local<10;++local){
				const auto& value=fluid_velocity_nodes_m_s[
					topology.cell_velocity_nodes[cell_index][local]];
				for(int component=0;component<3;++component)
					relative_velocity[component]+=basis.velocity[local]*value[component];
			}
			for(std::size_t local=0;local<4;++local){
				const auto& value=mesh_velocity_nodes_m_s[cell.nodes[local]];
				for(int component=0;component<3;++component)
					relative_velocity[component]-=basis.pressure[local]*value[component];
			}
			for(std::size_t row=0;row<4;++row){
				const auto& gradient=geometry.barycentric_gradients[row];
				const double direction=gradient[0]*relative_velocity[0]
					+gradient[1]*relative_velocity[1]
					+gradient[2]*relative_velocity[2];
				for(std::size_t column=0;column<4;++column)
					matrix[{cell.nodes[row],cell.nodes[column]}]
						-=geometry.determinant*point.weight*direction
							*basis.pressure[column];
			}
		}
		if(monotone)native_tet_moving_species_detail::AddMonotoneGraphDiffusion(
			matrix,cell.nodes,before);
	}
	const auto owners=native_tet_moving_species_detail::BoundaryOwners(current_mesh);
	std::map<native_tet_moving_species_detail::Face,unsigned> labelled;
	// Positive degree-five triangle quadrature, with weights summing to one.
	const std::array<std::array<double,4>,7> face_quadrature{{
		{{1./3.,1./3.,1./3.,0.225}},
		{{0.059715871789770,0.470142064105115,0.470142064105115,0.132394152788506}},
		{{0.470142064105115,0.059715871789770,0.470142064105115,0.132394152788506}},
		{{0.470142064105115,0.470142064105115,0.059715871789770,0.132394152788506}},
		{{0.797426985353087,0.101286507323456,0.101286507323456,0.125939180544827}},
		{{0.101286507323456,0.797426985353087,0.101286507323456,0.125939180544827}},
		{{0.101286507323456,0.101286507323456,0.797426985353087,0.125939180544827}}}};
	for(const auto& triangle:current_mesh.boundary_triangles){
		auto face=triangle.nodes;std::sort(face.begin(),face.end());
		const auto owner=owners.find(face);
		if(owner==owners.end()||owner->second.uses!=1||++labelled[face]!=1)
			throw std::invalid_argument("native moving species boundary face is invalid");
		const auto cell_index=owner->second.cell;
		const bool owned=cell_index%static_cast<std::size_t>(ranks)
			==static_cast<std::size_t>(rank);
		const auto& cell=current_mesh.cells[cell_index];
		const auto before=monotone&&owned
			?native_tet_moving_species_detail::ElementBlock(matrix,cell.nodes)
			:std::array<double,16>{};
		std::array<std::size_t,3> local_nodes{};
		for(std::size_t vertex=0;vertex<3;++vertex){
			const auto found=std::find(cell.nodes.begin(),cell.nodes.end(),
				triangle.nodes[vertex]);
			if(found==cell.nodes.end())
				throw std::invalid_argument("native moving species boundary owner is invalid");
			local_nodes[vertex]=static_cast<std::size_t>(found-cell.nodes.begin());
		}
		const auto& x0=current_mesh.points[triangle.nodes[0]];
		const auto& x1=current_mesh.points[triangle.nodes[1]];
		const auto& x2=current_mesh.points[triangle.nodes[2]];
		std::array<double,3> a{},b{},area_vector{};
		for(int component=0;component<3;++component){
			a[component]=x1[component]-x0[component];
			b[component]=x2[component]-x0[component];
		}
		area_vector={{0.5*(a[1]*b[2]-a[2]*b[1]),
			0.5*(a[2]*b[0]-a[0]*b[2]),0.5*(a[0]*b[1]-a[1]*b[0])}};
		const auto& interior=current_mesh.points[
			cell.nodes[owner->second.opposite]];
		double inward=0.,area_squared=0.;
		for(int component=0;component<3;++component){
			inward+=area_vector[component]*(interior[component]-x0[component]);
			area_squared+=area_vector[component]*area_vector[component];
		}
		if(!(area_squared>0.)||!std::isfinite(area_squared))
			throw std::invalid_argument("native moving species boundary area is invalid");
		if(inward>0.)for(auto& component:area_vector)component=-component;
		const auto exchange=wall_exchange.find(triangle.boundary_label);
		const double face_area=std::sqrt(area_squared);
		for(const auto& sample:face_quadrature){
			std::array<double,4> lambda{};
			for(std::size_t vertex=0;vertex<3;++vertex)
				lambda[local_nodes[vertex]]=sample[vertex];
			const auto basis=EvaluateNativeTaylorHoodBasis(lambda[1],lambda[2],lambda[3]);
			std::array<double,3> relative_velocity{};
			double speed_bound=0.;
			if(use_rt0){
				std::array<double,3> position{};
				for(std::size_t local=0;local<4;++local)
					for(int component=0;component<3;++component)
						position[component]+=lambda[local]
							*current_mesh.points.at(cell.nodes[local])[component];
				relative_velocity=EvaluateNativeTetRt0Flux(current_mesh,
					rt0_face_flow_m3_s,cell_index,position);
				speed_bound=std::hypot(relative_velocity[0],relative_velocity[1],
					relative_velocity[2]);
			}else for(std::size_t local=0;local<10;++local){
				const auto& value=fluid_velocity_nodes_m_s[
					topology.cell_velocity_nodes[cell_index][local]];
				speed_bound+=std::abs(basis.velocity[local])
					*std::sqrt(value[0]*value[0]+value[1]*value[1]+value[2]*value[2]);
				for(int component=0;component<3;++component)
					relative_velocity[component]+=basis.velocity[local]*value[component];
			}
			for(std::size_t local=0;local<4;++local){
				const auto& value=mesh_velocity_nodes_m_s[cell.nodes[local]];
				speed_bound+=lambda[local]
					*std::sqrt(value[0]*value[0]+value[1]*value[1]+value[2]*value[2]);
				for(int component=0;component<3;++component)
					relative_velocity[component]-=lambda[local]*value[component];
			}
			double normal_flux=0.;
			for(int component=0;component<3;++component)
				normal_flux+=relative_velocity[component]*area_vector[component];
			normal_flux*=sample[3];
			const double roundoff=64.*std::numeric_limits<double>::epsilon()
				*std::sqrt(area_squared)*sample[3]
				*std::max(speed_bound,global_speed_bound);
			if(std::abs(normal_flux)<=roundoff)normal_flux=0.;
			if(exchange!=wall_exchange.end()&&normal_flux!=0.)
				throw std::invalid_argument("native moving species wall exchange face has relative flow");
			const auto prescribed=inflow_concentration_mol_m3.find(
				triangle.boundary_label);
			if(normal_flux<0.&&prescribed==inflow_concentration_mol_m3.end()){
				std::ostringstream message;
				message<<"native moving species inflow on boundary label "
					<<triangle.boundary_label<<" has no concentration (quadrature flow "
					<<std::scientific<<normal_flux<<" m^3/s, roundoff "
					<<roundoff<<" m^3/s)";
				throw std::invalid_argument(message.str());
			}
			if(!owned)continue;
			if(exchange!=wall_exchange.end()){
				const double transfer=exchange->second.transfer_coefficient_m_s
					*face_area*sample[3];
				const double external=exchange->second.external_concentration_mol_m3;
				assembly.wall_exchange_constant_by_label_mol_s.at(
					triangle.boundary_label)-=transfer*external;
				for(std::size_t row=0;row<4;++row){
					rhs[cell.nodes[row]]+=transfer*external*lambda[row];
					for(std::size_t column=0;column<4;++column){
						matrix[{cell.nodes[row],cell.nodes[column]}]
							+=transfer*lambda[row]*lambda[column];
						assembly.wall_exchange_coefficients_by_label.at(
							triangle.boundary_label)[cell.nodes[column]]
							+=transfer*lambda[row]*lambda[column];
					}
				}
			}
			assembly.outward_relative_flow_by_label_m3_s.at(
				triangle.boundary_label)+=normal_flux;
			if(normal_flux>=0.)
				assembly.positive_relative_flow_by_label_m3_s.at(
					triangle.boundary_label)+=normal_flux;
			else assembly.negative_relative_flow_by_label_m3_s.at(
				triangle.boundary_label)+=normal_flux;
			if(normal_flux>=0.){
				for(std::size_t row=0;row<4;++row)
					for(std::size_t column=0;column<4;++column){
						const auto contribution=lambda[row]*lambda[column]*normal_flux;
						matrix[{cell.nodes[row],cell.nodes[column]}]+=contribution;
						outward_coefficients[cell.nodes[column]]+=contribution;
						assembly.outward_coefficients_by_label.at(
							triangle.boundary_label)[cell.nodes[column]]+=contribution;
					}
			}else{
				incoming_flux_mol_s+=normal_flux*prescribed->second;
				assembly.incoming_flux_by_label_mol_s.at(
					triangle.boundary_label)+=normal_flux*prescribed->second;
				for(std::size_t row=0;row<4;++row)
					rhs[cell.nodes[row]]-=lambda[row]*normal_flux*prescribed->second;
			}
		}
		if(monotone&&owned)
			native_tet_moving_species_detail::AddMonotoneGraphDiffusion(
				matrix,cell.nodes,before);
	}
	for(const auto& item:owners)
		if(item.second.uses==1&&labelled[item.first]!=1)
			throw std::invalid_argument("native moving species boundary is unlabelled");
	return assembly;
}

inline NativeTetMovingSpeciesStep FinishNativeTetMovingSpeciesStep(
	const NativeTetMesh& previous_mesh,const NativeTetMesh& current_mesh,
	const std::vector<double>& previous_concentration_mol_m3,
	std::vector<double> current_concentration_mol_m3,
	const NativeTetMovingSpeciesAssembly& assembly,double dt_s,
	bool require_nonnegative=false)
{
	if(current_concentration_mol_m3.size()!=current_mesh.points.size()
		||assembly.outward_coefficients.size()!=current_mesh.points.size()
		||!(dt_s>0.))
		throw std::invalid_argument("native moving species result size is invalid");
	for(const double value:current_concentration_mol_m3)
		if(!std::isfinite(value))
			throw std::runtime_error("native moving species result is nonfinite");
	if(require_nonnegative){
		double scale=1.;
		for(const double value:current_concentration_mol_m3)
			scale=std::max(scale,std::abs(value));
		for(const double value:current_concentration_mol_m3)
			if(value<-1e-12*scale)
				throw std::runtime_error("native moving species monotone result is negative");
	}
	NativeTetMovingSpeciesStep result;
	result.concentration_mol_m3=std::move(current_concentration_mol_m3);
	result.previous_inventory_mol=native_tet_moving_species_detail::Inventory(
		previous_mesh,previous_concentration_mol_m3);
	result.current_inventory_mol=native_tet_moving_species_detail::Inventory(
		current_mesh,result.concentration_mol_m3);
	result.outward_advective_flux_mol_s=assembly.incoming_flux_mol_s;
	result.outward_relative_flow_by_label_m3_s=
		assembly.outward_relative_flow_by_label_m3_s;
	result.positive_relative_flow_by_label_m3_s=
		assembly.positive_relative_flow_by_label_m3_s;
	result.negative_relative_flow_by_label_m3_s=
		assembly.negative_relative_flow_by_label_m3_s;
	for(const auto& item:assembly.outward_coefficients_by_label){
		if(item.second.size()!=assembly.nodes)
			throw std::invalid_argument("native moving species labelled flux size is invalid");
		double labelled_flux=assembly.incoming_flux_by_label_mol_s.at(item.first);
		for(std::size_t node=0;node<assembly.nodes;++node)
			labelled_flux+=item.second[node]*result.concentration_mol_m3[node];
		result.outward_advective_flux_by_label_mol_s.emplace(item.first,labelled_flux);
	}
	for(std::size_t node=0;node<assembly.nodes;++node)
		result.outward_advective_flux_mol_s+=assembly.outward_coefficients[node]
			*result.concentration_mol_m3[node];
	result.source_mol_s=assembly.source_mol_s;
	result.reaction_sink_mol_s=assembly.first_order_decay_rate_s_inv
		*result.current_inventory_mol;
	for(const auto& item:assembly.wall_exchange_coefficients_by_label){
		if(item.second.size()!=assembly.nodes)
			throw std::invalid_argument("native moving species wall exchange size is invalid");
		double flux=assembly.wall_exchange_constant_by_label_mol_s.at(item.first);
		for(std::size_t node=0;node<assembly.nodes;++node)
			flux+=item.second[node]*result.concentration_mol_m3[node];
		result.outward_wall_exchange_by_label_mol_s.emplace(item.first,flux);
		result.outward_wall_exchange_mol_s+=flux;
	}
	result.balance_defect_mol_s=(result.current_inventory_mol
		-result.previous_inventory_mol)/dt_s
		+result.outward_advective_flux_mol_s+result.reaction_sink_mol_s
		+result.outward_wall_exchange_mol_s
		-result.source_mol_s;
	return result;
}

// Bounded dense reference algebra; the PETSc path consumes the same sparse
// project-owned assembly without materializing this global dense matrix.
inline NativeTetMovingSpeciesStep SolveNativeTetMovingSpeciesDenseStep(
	const NativeTetMesh& previous_mesh,const NativeTetMesh& current_mesh,
	const std::vector<std::array<double,3>>& fluid_velocity_nodes_m_s,
	const std::vector<std::array<double,3>>& mesh_velocity_nodes_m_s,
	const std::vector<double>& previous_concentration_mol_m3,
	const std::map<int,double>& inflow_concentration_mol_m3,
	double diffusivity_m2_s,double source_mol_m3_s,double dt_s,
	bool monotone=false,double first_order_decay_rate_s_inv=0.,
	const std::map<int,NativeTetWallExchange>& wall_exchange={},
	const std::vector<std::array<double,4>>& rt0_face_flow_m3_s={},
	const std::vector<double>& cell_source_mol_m3_s={})
{
	if(current_mesh.points.size()>128)
		throw std::invalid_argument("native moving species dense limit is 128 vertices");
	auto assembly=AssembleNativeTetMovingSpeciesStep(previous_mesh,current_mesh,
		fluid_velocity_nodes_m_s,mesh_velocity_nodes_m_s,
		previous_concentration_mol_m3,inflow_concentration_mol_m3,
		diffusivity_m2_s,source_mol_m3_s,dt_s,0,1,monotone,
		first_order_decay_rate_s_inv,wall_exchange,rt0_face_flow_m3_s,
		cell_source_mol_m3_s);
	std::vector<double> matrix(assembly.nodes*assembly.nodes,0.);
	for(const auto& entry:assembly.matrix)
		matrix[static_cast<std::size_t>(entry.first.first)*assembly.nodes
			+entry.first.second]=entry.second;
	auto solution=native_tet_moving_species_detail::SolveDense(
		std::move(matrix),assembly.rhs);
	return FinishNativeTetMovingSpeciesStep(previous_mesh,current_mesh,
		previous_concentration_mol_m3,std::move(solution),assembly,dt_s,
		monotone);
}

// Adapter for the native Taylor-Hood flow state, whose pressure tail is not
// part of the scalar transport velocity. No external FEM field object enters.
inline NativeTetMovingSpeciesStep SolveNativeTetMovingSpeciesDenseStep(
	const NativeTetMesh& previous_mesh,const NativeTetMesh& current_mesh,
	const std::vector<double>& native_flow_state,
	const std::vector<std::array<double,3>>& mesh_velocity_nodes_m_s,
	const std::vector<double>& previous_concentration_mol_m3,
	const std::map<int,double>& inflow_concentration_mol_m3,
	double diffusivity_m2_s,double source_mol_m3_s,double dt_s,
	bool monotone=false,double first_order_decay_rate_s_inv=0.,
	const std::map<int,NativeTetWallExchange>& wall_exchange={})
{
	const auto topology=BuildNativeTaylorHoodTopology(current_mesh);
	const std::size_t velocity_nodes=current_mesh.points.size()+topology.edges.size();
	if(native_flow_state.size()!=3*velocity_nodes+current_mesh.points.size())
		throw std::invalid_argument("native moving species flow state size is invalid");
	std::vector<std::array<double,3>> velocity(velocity_nodes);
	for(std::size_t node=0;node<velocity_nodes;++node)
		for(int component=0;component<3;++component)
			velocity[node][component]=native_flow_state[3*node+component];
	return SolveNativeTetMovingSpeciesDenseStep(previous_mesh,current_mesh,
		velocity,mesh_velocity_nodes_m_s,previous_concentration_mol_m3,
		inflow_concentration_mol_m3,diffusivity_m2_s,source_mol_m3_s,dt_s,
		monotone,first_order_decay_rate_s_inv,wall_exchange);
}

} // namespace iga

#endif
