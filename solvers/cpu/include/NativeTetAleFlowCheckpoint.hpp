#ifndef IGA_NATIVE_TET_ALE_FLOW_CHECKPOINT_HPP
#define IGA_NATIVE_TET_ALE_FLOW_CHECKPOINT_HPP

#include "CheckpointMetadataCodec.hpp"
#include "NativeTetAlePetscRuntime.hpp"
#include "Sha256.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace iga {

struct NativeTetAleFlowCheckpoint
{
	std::string model_identity_sha256,state_identity_sha256;
	std::uint64_t accepted_steps=0;
	double time_s=0.;
	std::vector<std::array<double,3>> current_points_m;
	std::vector<double> flow_state;
};

namespace native_tet_ale_flow_checkpoint_detail {

constexpr std::uint64_t maximum_nodes=250000;
constexpr std::uint64_t maximum_flow_values=1000000;

inline void AppendText(Sha256& hash,const std::string& value)
{
	hash.AppendLittleEndian64(value.size());hash.Append(value.data(),value.size());
}

inline bool IsSha256(const std::string& value)
{
	if(value.size()!=64)return false;
	for(const char c:value)
		if(!((c>='0'&&c<='9')||(c>='a'&&c<='f')))return false;
	return true;
}

inline std::string StateIdentity(const NativeTetAleFlowCheckpoint& value)
{
	Sha256 hash;AppendText(hash,"NativeTetAleFlowState/v1");
	AppendText(hash,value.model_identity_sha256);
	hash.AppendLittleEndian64(value.accepted_steps);
	hash.AppendNormalizedDouble(value.time_s);
	hash.AppendLittleEndian64(value.current_points_m.size());
	for(const auto& point:value.current_points_m)
		for(const double component:point)hash.AppendNormalizedDouble(component);
	hash.AppendLittleEndian64(value.flow_state.size());
	for(const double component:value.flow_state)hash.AppendNormalizedDouble(component);
	return hash.Hex();
}

inline void Validate(const NativeTetAleFlowCheckpoint& value)
{
	using checkpoint_metadata::Require;
	Require(value.accepted_steps>0&&value.accepted_steps
		<=static_cast<std::uint64_t>(std::numeric_limits<int>::max())
		&&std::isfinite(value.time_s)&&value.time_s>0.,
		"native tetra flow checkpoint clock is invalid");
	Require(value.current_points_m.size()>0
		&&value.current_points_m.size()<=maximum_nodes
		&&value.flow_state.size()>0
		&&value.flow_state.size()<=maximum_flow_values,
		"native tetra flow checkpoint shape is invalid");
	Require(IsSha256(value.model_identity_sha256)
		&&IsSha256(value.state_identity_sha256),
		"native tetra flow checkpoint identity is invalid");
	for(const auto& point:value.current_points_m)
		for(const double component:point)
			Require(std::isfinite(component),
				"native tetra flow checkpoint point is nonfinite");
	for(const double component:value.flow_state)
		Require(std::isfinite(component),
			"native tetra flow checkpoint state is nonfinite");
	Require(value.state_identity_sha256==StateIdentity(value),
		"native tetra flow checkpoint state identity mismatch");
}

inline std::string Digest(std::string_view bytes)
{
	Sha256 hash;hash.Append(bytes.data(),bytes.size());return hash.Hex();
}

} // namespace native_tet_ale_flow_checkpoint_detail

inline std::string NativeTetAleFlowModelIdentitySha256(
	const std::string& domain_id,const NativeTetMesh& reference_mesh,
	const std::vector<CouplingPort>& ports,const NativeNavierStokesParameters& parameters,
	const std::set<int>& wall_labels,const std::string& motion_model_id,
	double nonlinear_tolerance,std::size_t maximum_iterations)
{
	if(domain_id.empty()||reference_mesh.points.empty()||reference_mesh.cells.empty()
		||!(parameters.density>0.)||!std::isfinite(parameters.density)
		||!(parameters.dynamic_viscosity>0.)
		||!std::isfinite(parameters.dynamic_viscosity)
		||!(nonlinear_tolerance>0.)||!std::isfinite(nonlinear_tolerance)
		||maximum_iterations==0)
		throw std::invalid_argument("native tetra flow checkpoint model is invalid");
	Sha256 hash;
	using native_tet_ale_flow_checkpoint_detail::AppendText;
	AppendText(hash,"NativeTetAleFlowModel/v1");
	AppendText(hash,domain_id);AppendText(hash,motion_model_id);
	hash.AppendNormalizedDouble(parameters.density);
	hash.AppendNormalizedDouble(parameters.dynamic_viscosity);
	for(const auto& acceleration:parameters.body_acceleration_m_s2)
		for(const double coefficient:acceleration){
			if(!std::isfinite(coefficient))
				throw std::invalid_argument("native tetra flow body acceleration is nonfinite");
			hash.AppendNormalizedDouble(coefficient);
		}
	hash.AppendNormalizedDouble(nonlinear_tolerance);
	hash.AppendLittleEndian64(maximum_iterations);
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
	for(const auto& face:reference_mesh.boundary_triangles){
		hash.AppendLittleEndian64(face.id);
		for(const auto node:face.nodes)hash.AppendLittleEndian32(node);
		hash.AppendLittleEndian32(static_cast<std::uint32_t>(face.boundary_label));
	}
	hash.AppendLittleEndian64(ports.size());
	for(const auto& port:ports){
		AppendText(hash,port.id);AppendText(hash,port.subsystem_id);
		AppendText(hash,port.locator_kind);
		AppendText(hash,port.locator);
		hash.AppendLittleEndian32(
			static_cast<std::uint32_t>(port.orientation.native_to_outward_sign));
		hash.AppendLittleEndian64(port.provides.size());
		for(const auto quantity:port.provides)
			hash.AppendLittleEndian32(static_cast<std::uint32_t>(quantity));
		hash.AppendLittleEndian64(port.requires.size());
		for(const auto quantity:port.requires)
			hash.AppendLittleEndian32(static_cast<std::uint32_t>(quantity));
		hash.AppendLittleEndian64(port.species.size());
		for(const auto& species:port.species)AppendText(hash,species);
	}
	hash.AppendLittleEndian64(wall_labels.size());
	for(const int label:wall_labels)
		hash.AppendLittleEndian32(static_cast<std::uint32_t>(label));
	return hash.Hex();
}

inline std::string SerializeNativeTetAleFlowCheckpoint(
	const NativeTetAleFlowCheckpoint& value)
{
	using namespace native_tet_ale_flow_checkpoint_detail;
	Validate(value);
	checkpoint_metadata::Writer output;
	output.Text("IGA_NATIVE_TET_ALE_FLOW");output.Unsigned(1);
	output.Text(value.model_identity_sha256);output.Text(value.state_identity_sha256);
	output.Unsigned(value.accepted_steps);output.Real(value.time_s);
	output.Unsigned(value.current_points_m.size());
	for(const auto& point:value.current_points_m)
		for(const double component:point)output.Real(component);
	output.Unsigned(value.flow_state.size());
	for(const double component:value.flow_state)output.Real(component);
	auto bytes=output.Bytes();
	checkpoint_metadata::Writer checksum;checksum.Text(Digest(bytes));
	checkpoint_metadata::Require(bytes.size()+checksum.Bytes().size()
		<=checkpoint_metadata::maximum_bytes,
		"native tetra flow checkpoint exceeds size limit");
	bytes+=checksum.Bytes();return bytes;
}

inline NativeTetAleFlowCheckpoint ParseNativeTetAleFlowCheckpoint(
	std::string_view bytes)
{
	using namespace native_tet_ale_flow_checkpoint_detail;
	using checkpoint_metadata::Require;
	constexpr std::size_t checksum_bytes=72;
	Require(bytes.size()>=checksum_bytes&&bytes.size()<=checkpoint_metadata::maximum_bytes,
		"native tetra flow checkpoint size is invalid");
	const auto payload=bytes.substr(0,bytes.size()-checksum_bytes);
	checkpoint_metadata::Reader checksum(bytes.substr(bytes.size()-checksum_bytes));
	const auto expected=checksum.Text();checksum.Finish();
	Require(expected==Digest(payload),"native tetra flow checkpoint checksum mismatch");
	checkpoint_metadata::Reader input(payload);
	Require(input.Text()=="IGA_NATIVE_TET_ALE_FLOW"&&input.Unsigned()==1,
		"unsupported native tetra flow checkpoint schema");
	NativeTetAleFlowCheckpoint value;
	value.model_identity_sha256=input.Text();value.state_identity_sha256=input.Text();
	value.accepted_steps=input.Unsigned();value.time_s=input.Real();
	const auto nodes=input.Unsigned();
	Require(nodes<=maximum_nodes,"native tetra flow checkpoint nodes exceed limit");
	value.current_points_m.resize(static_cast<std::size_t>(nodes));
	for(auto& point:value.current_points_m)
		for(double& component:point)component=input.Real();
	const auto scalars=input.Unsigned();
	Require(scalars<=maximum_flow_values,
		"native tetra flow checkpoint state exceeds limit");
	value.flow_state.resize(static_cast<std::size_t>(scalars));
	for(double& component:value.flow_state)component=input.Real();
	input.Finish();Validate(value);return value;
}

} // namespace iga

#endif
