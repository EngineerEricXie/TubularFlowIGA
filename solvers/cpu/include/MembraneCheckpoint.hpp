#ifndef IGA_MEMBRANE_CHECKPOINT_HPP
#define IGA_MEMBRANE_CHECKPOINT_HPP

#include "PretensionedMembrane.hpp"
#include "CheckpointMetadataCodec.hpp"

namespace iga {

// Numerical owner payload. The runtime/bundle authority supplies the accepted
// clock; publication stamps and atomic restoration belong to that outer layer.
inline std::string SerializeMembraneCheckpoint(const PretensionedMembrane& membrane,
	std::uint64_t accepted_step,double accepted_time,std::size_t maximum_nodes=4096)
{
	using checkpoint_metadata::Require;
	const auto& ids=membrane.Layout().owned_global_node_ids;
	const auto& state=membrane.CommittedState();
	Require(accepted_step>0&&std::isfinite(accepted_time)&&accepted_time>=0,"invalid membrane accepted clock");
	Require(!ids.empty()&&ids.size()<=maximum_nodes,"membrane checkpoint exceeds node limit");
	Require(membrane.CheckpointStateIdentitySha256(state)==membrane.CommittedStateIdentitySha256(),"membrane accepted identity mismatch");
	checkpoint_metadata::Writer output;output.Text("IGA_MEMBRANE_STATE/1");
	output.Text(membrane.ModelIdentitySha256());output.Text(membrane.CommittedStateIdentitySha256());
	output.Unsigned(accepted_step);output.Real(accepted_time);output.Unsigned(ids.size());
	for(std::size_t node=0;node<ids.size();++node) {
		output.Unsigned(ids[node]);output.Real(state.displacement_m[node]);output.Real(state.velocity_m_per_s[node]);
	}
	return output.Bytes();
}

inline PretensionedMembraneState ParseMembraneCheckpoint(std::string_view bytes,
	const PretensionedMembrane& model,const std::string& expected_state_identity,
	std::uint64_t expected_step,double expected_time,std::size_t maximum_nodes=4096)
{
	using checkpoint_metadata::Require;
	Require(expected_step>0&&std::isfinite(expected_time)&&expected_time>=0,"invalid expected membrane clock");
	checkpoint_metadata::Reader input(bytes);
	Require(input.Text()=="IGA_MEMBRANE_STATE/1","unsupported membrane checkpoint version");
	Require(input.Text()==model.ModelIdentitySha256(),"membrane checkpoint model differs");
	const auto identity=input.Text();Require(identity==expected_state_identity,"membrane checkpoint differs from expected state");
	Require(input.Unsigned()==expected_step&&input.Real()==expected_time,"membrane checkpoint epoch differs");
	const auto count=input.Unsigned();const auto& ids=model.Layout().owned_global_node_ids;
	Require(count>0&&count==ids.size()&&count<=maximum_nodes&&count<=bytes.size()/24,"invalid membrane checkpoint count");
	PretensionedMembraneState state;state.displacement_m.resize(count);state.velocity_m_per_s.resize(count);
	for(std::size_t node=0;node<count;++node) {
		Require(input.Unsigned()==ids[node],"membrane checkpoint stable node IDs differ");
		state.displacement_m[node]=input.Real();state.velocity_m_per_s[node]=input.Real();
	}
	input.Finish();
	Require(model.CheckpointStateIdentitySha256(state)==identity,"membrane checkpoint state content differs");
	return state;
}

} // namespace iga
#endif
