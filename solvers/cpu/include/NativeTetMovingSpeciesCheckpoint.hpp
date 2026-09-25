#ifndef IGA_NATIVE_TET_MOVING_SPECIES_CHECKPOINT_HPP
#define IGA_NATIVE_TET_MOVING_SPECIES_CHECKPOINT_HPP

#include "CheckpointMetadataCodec.hpp"
#include "NativeTetMovingSpeciesPetscRuntime.hpp"
#include "Sha256.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace iga {

struct NativeTetMovingSpeciesCommittedState
{
	NativeTetMesh reference_mesh,current_mesh;
	std::vector<double> concentration_mol_m3;
	double diffusivity_m2_s=0.;
	std::uint64_t accepted_steps=0;
	double time_s=0.;
	std::string model_identity_sha256;
};

struct NativeTetMovingSpeciesCheckpoint
{
	std::uint64_t accepted_steps=0;
	double time_s=0.;
	std::string model_identity_sha256,state_identity_sha256;
	std::vector<std::array<double,3>> current_points_m;
	std::vector<double> concentration_mol_m3;
};

namespace native_tet_moving_species_checkpoint_detail {

constexpr std::uint64_t maximum_nodes=250000;

inline bool IsSha256(const std::string& value)
{
	if(value.size()!=64)return false;
	for(const char c:value)
		if(!((c>='0'&&c<='9')||(c>='a'&&c<='f')))return false;
	return true;
}

inline void AppendTag(Sha256& hash,const char* value,std::size_t size)
{
	hash.AppendLittleEndian64(size);hash.Append(value,size);
}

inline std::string StateIdentity(const NativeTetMovingSpeciesCheckpoint& value)
{
	Sha256 hash;constexpr char tag[]="NativeTetMovingSpeciesState/v1";
	AppendTag(hash,tag,sizeof(tag)-1);
	AppendTag(hash,value.model_identity_sha256.data(),
		value.model_identity_sha256.size());
	hash.AppendLittleEndian64(value.accepted_steps);
	hash.AppendNormalizedDouble(value.time_s);
	hash.AppendLittleEndian64(value.current_points_m.size());
	for(const auto& point:value.current_points_m)
		for(const double component:point)hash.AppendNormalizedDouble(component);
	hash.AppendLittleEndian64(value.concentration_mol_m3.size());
	for(const double concentration:value.concentration_mol_m3)
		hash.AppendNormalizedDouble(concentration);
	return hash.Hex();
}

inline void Validate(const NativeTetMovingSpeciesCheckpoint& value)
{
	using checkpoint_metadata::Require;
	Require(value.current_points_m.size()>0
		&&value.current_points_m.size()<=maximum_nodes
		&&value.current_points_m.size()==value.concentration_mol_m3.size(),
		"native moving species checkpoint shape is invalid");
	Require(std::isfinite(value.time_s)&&value.time_s>=0.
		&&((value.accepted_steps==0)==(value.time_s==0.)),
		"native moving species checkpoint clock is invalid");
	Require(IsSha256(value.model_identity_sha256)
		&&IsSha256(value.state_identity_sha256),
		"native moving species checkpoint identity is invalid");
	for(const auto& point:value.current_points_m)
		for(const double component:point)
			Require(std::isfinite(component),"native moving species checkpoint point is nonfinite");
	for(const double concentration:value.concentration_mol_m3)
		Require(std::isfinite(concentration)&&concentration>=0.,
			"native moving species checkpoint concentration is invalid");
	Require(value.state_identity_sha256==StateIdentity(value),
		"native moving species checkpoint state identity mismatch");
}

inline std::string Digest(std::string_view bytes)
{
	Sha256 hash;hash.Append(bytes.data(),bytes.size());return hash.Hex();
}

inline bool SameTopology(const NativeTetMesh& current,const NativeTetMesh& reference)
{
	if(current.points.size()!=reference.points.size()
		||current.cells.size()!=reference.cells.size()
		||current.boundary_triangles.size()!=reference.boundary_triangles.size())
		return false;
	for(std::size_t index=0;index<reference.cells.size();++index)
		if(current.cells[index].id!=reference.cells[index].id
			||current.cells[index].nodes!=reference.cells[index].nodes)return false;
	for(std::size_t index=0;index<reference.boundary_triangles.size();++index){
		const auto& face=current.boundary_triangles[index];
		const auto& expected=reference.boundary_triangles[index];
		if(face.id!=expected.id||face.nodes!=expected.nodes
			||face.boundary_label!=expected.boundary_label)return false;
	}
	return true;
}

} // namespace native_tet_moving_species_checkpoint_detail

inline std::string NativeTetMovingSpeciesModelIdentitySha256(
	const NativeTetMesh& reference_mesh,double diffusivity_m2_s)
{
	if(reference_mesh.points.empty()||reference_mesh.cells.empty()
		||!(diffusivity_m2_s>=0.)||!std::isfinite(diffusivity_m2_s))
		throw std::invalid_argument("native moving species model is invalid");
	Sha256 hash;constexpr char tag[]="NativeTetMovingSpeciesModel/v1";
	native_tet_moving_species_checkpoint_detail::AppendTag(hash,tag,sizeof(tag)-1);
	hash.AppendNormalizedDouble(diffusivity_m2_s);
	hash.AppendLittleEndian64(reference_mesh.points.size());
	for(const auto& point:reference_mesh.points)
		for(const double component:point)hash.AppendNormalizedDouble(component);
	hash.AppendLittleEndian64(reference_mesh.cells.size());
	for(const auto& cell:reference_mesh.cells){
		hash.AppendLittleEndian64(cell.id);
		for(const auto node:cell.nodes)hash.AppendLittleEndian32(node);
		EvaluateNativeTetGeometry(reference_mesh,cell);
	}
	hash.AppendLittleEndian64(reference_mesh.boundary_triangles.size());
	for(const auto& triangle:reference_mesh.boundary_triangles){
		hash.AppendLittleEndian64(triangle.id);
		for(const auto node:triangle.nodes)hash.AppendLittleEndian32(node);
		hash.AppendLittleEndian32(static_cast<std::uint32_t>(triangle.boundary_label));
	}
	return hash.Hex();
}

inline NativeTetMovingSpeciesCommittedState InitializeNativeTetMovingSpeciesState(
	const NativeTetMesh& reference_mesh,std::vector<double> concentration_mol_m3,
	double diffusivity_m2_s)
{
	if(concentration_mol_m3.size()!=reference_mesh.points.size())
		throw std::invalid_argument("native moving species initial state size is invalid");
	for(const double value:concentration_mol_m3)
		if(!std::isfinite(value)||value<0.)
			throw std::invalid_argument("native moving species initial concentration is invalid");
	NativeTetMovingSpeciesCommittedState state;
	state.reference_mesh=reference_mesh;
	state.current_mesh=reference_mesh;
	state.concentration_mol_m3=std::move(concentration_mol_m3);
	state.diffusivity_m2_s=diffusivity_m2_s;
	state.model_identity_sha256=NativeTetMovingSpeciesModelIdentitySha256(
		reference_mesh,diffusivity_m2_s);
	return state;
}

// A failed solve or budget gate leaves the committed object unchanged.
inline NativeTetMovingSpeciesPetscResult AdvanceNativeTetMovingSpeciesState(
	NativeTetMovingSpeciesCommittedState& state,const NativeTetMesh& trial_mesh,
	const std::vector<std::array<double,3>>& fluid_velocity_nodes_m_s,
	const std::vector<std::array<double,3>>& mesh_velocity_nodes_m_s,
	const std::map<int,double>& inflow_concentration_mol_m3,
	double source_mol_m3_s,double dt_s)
{
	if(state.model_identity_sha256!=NativeTetMovingSpeciesModelIdentitySha256(
		state.reference_mesh,state.diffusivity_m2_s)
		||!std::isfinite(state.time_s)||state.time_s<0.
		||(state.accepted_steps==0)!=(state.time_s==0.)
		||!native_tet_moving_species_checkpoint_detail::SameTopology(
			state.current_mesh,state.reference_mesh)
		||!native_tet_moving_species_checkpoint_detail::SameTopology(
			trial_mesh,state.reference_mesh)
		||state.accepted_steps==std::numeric_limits<std::uint64_t>::max())
		throw std::invalid_argument("native moving species committed state is invalid");
	auto result=SolveNativeTetMovingSpeciesPetscStep(state.current_mesh,trial_mesh,
		fluid_velocity_nodes_m_s,mesh_velocity_nodes_m_s,
		state.concentration_mol_m3,inflow_concentration_mol_m3,
		state.diffusivity_m2_s,source_mol_m3_s,dt_s);
	const double scale=std::max({1e-30,
		std::abs(result.step.previous_inventory_mol/dt_s),
		std::abs(result.step.outward_advective_flux_mol_s),
		std::abs(result.step.source_mol_s)});
	if(!std::isfinite(result.step.balance_defect_mol_s)
		||std::abs(result.step.balance_defect_mol_s)>1e-12*scale+1e-14)
		throw std::runtime_error("native moving species accepted-step mass balance failed");
	for(const double value:result.step.concentration_mol_m3)
		if(!std::isfinite(value)||value<0.)
			throw std::runtime_error("native moving species accepted-step concentration is invalid");
	NativeTetMovingSpeciesCommittedState candidate=state;
	candidate.current_mesh=trial_mesh;
	candidate.concentration_mol_m3=result.step.concentration_mol_m3;
	candidate.time_s+=dt_s;
	if(!std::isfinite(candidate.time_s))
		throw std::runtime_error("native moving species accepted-step time is invalid");
	++candidate.accepted_steps;
	state=std::move(candidate);
	return result;
}

inline NativeTetMovingSpeciesCheckpoint CaptureNativeTetMovingSpeciesCheckpoint(
	const NativeTetMovingSpeciesCommittedState& state)
{
	if(state.model_identity_sha256!=NativeTetMovingSpeciesModelIdentitySha256(
		state.reference_mesh,state.diffusivity_m2_s))
		throw std::invalid_argument("native moving species checkpoint model identity is stale");
	if(!native_tet_moving_species_checkpoint_detail::SameTopology(
		state.current_mesh,state.reference_mesh))
		throw std::invalid_argument("native moving species checkpoint topology changed");
	for(const auto& cell:state.current_mesh.cells)
		EvaluateNativeTetGeometry(state.current_mesh,cell);
	NativeTetMovingSpeciesCheckpoint value;
	value.accepted_steps=state.accepted_steps;
	value.time_s=state.time_s;
	value.model_identity_sha256=state.model_identity_sha256;
	value.current_points_m=state.current_mesh.points;
	value.concentration_mol_m3=state.concentration_mol_m3;
	value.state_identity_sha256=
		native_tet_moving_species_checkpoint_detail::StateIdentity(value);
	native_tet_moving_species_checkpoint_detail::Validate(value);
	return value;
}

inline std::string SerializeNativeTetMovingSpeciesCheckpoint(
	const NativeTetMovingSpeciesCheckpoint& value)
{
	using namespace native_tet_moving_species_checkpoint_detail;
	Validate(value);
	checkpoint_metadata::Writer output;
	output.Text("IGA_NATIVE_MOVING_SPECIES");output.Unsigned(1);
	output.Unsigned(value.accepted_steps);output.Real(value.time_s);
	output.Text(value.model_identity_sha256);
	output.Text(value.state_identity_sha256);
	output.Unsigned(value.current_points_m.size());
	for(const auto& point:value.current_points_m)
		for(const double component:point)output.Real(component);
	output.Unsigned(value.concentration_mol_m3.size());
	for(const double concentration:value.concentration_mol_m3)
		output.Real(concentration);
	auto bytes=output.Bytes();
	checkpoint_metadata::Writer checksum;checksum.Text(Digest(bytes));
	checkpoint_metadata::Require(bytes.size()+checksum.Bytes().size()
		<=checkpoint_metadata::maximum_bytes,
		"native moving species checkpoint exceeds size limit");
	bytes+=checksum.Bytes();return bytes;
}

inline NativeTetMovingSpeciesCheckpoint ParseNativeTetMovingSpeciesCheckpoint(
	std::string_view bytes)
{
	using namespace native_tet_moving_species_checkpoint_detail;
	using checkpoint_metadata::Require;
	constexpr std::size_t checksum_bytes=72;
	Require(bytes.size()>=checksum_bytes&&bytes.size()<=checkpoint_metadata::maximum_bytes,
		"native moving species checkpoint size is invalid");
	const auto payload=bytes.substr(0,bytes.size()-checksum_bytes);
	checkpoint_metadata::Reader checksum(bytes.substr(bytes.size()-checksum_bytes));
	const auto expected=checksum.Text();checksum.Finish();
	Require(expected==Digest(payload),"native moving species checkpoint checksum mismatch");
	checkpoint_metadata::Reader input(payload);
	Require(input.Text()=="IGA_NATIVE_MOVING_SPECIES"&&input.Unsigned()==1,
		"unsupported native moving species checkpoint schema");
	NativeTetMovingSpeciesCheckpoint value;
	value.accepted_steps=input.Unsigned();value.time_s=input.Real();
	value.model_identity_sha256=input.Text();
	value.state_identity_sha256=input.Text();
	const auto points=input.Unsigned();
	Require(points<=maximum_nodes,"native moving species checkpoint points exceed limit");
	value.current_points_m.resize(static_cast<std::size_t>(points));
	for(auto& point:value.current_points_m)
		for(double& component:point)component=input.Real();
	const auto concentrations=input.Unsigned();
	Require(concentrations==points,
		"native moving species checkpoint concentration count is invalid");
	value.concentration_mol_m3.resize(static_cast<std::size_t>(concentrations));
	for(double& concentration:value.concentration_mol_m3)
		concentration=input.Real();
	input.Finish();Validate(value);return value;
}

inline NativeTetMovingSpeciesCommittedState RestoreNativeTetMovingSpeciesCheckpoint(
	const NativeTetMovingSpeciesCheckpoint& value,
	const NativeTetMesh& reference_mesh,double diffusivity_m2_s)
{
	native_tet_moving_species_checkpoint_detail::Validate(value);
	auto state=InitializeNativeTetMovingSpeciesState(reference_mesh,
		value.concentration_mol_m3,diffusivity_m2_s);
	if(value.model_identity_sha256!=state.model_identity_sha256
		||value.current_points_m.size()!=reference_mesh.points.size())
		throw std::invalid_argument("native moving species checkpoint model mismatch");
	state.current_mesh.points=value.current_points_m;
	for(const auto& cell:state.current_mesh.cells)
		EvaluateNativeTetGeometry(state.current_mesh,cell);
	state.accepted_steps=value.accepted_steps;
	state.time_s=value.time_s;
	return state;
}

} // namespace iga

#endif
