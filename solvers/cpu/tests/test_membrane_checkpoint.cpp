#include "MembraneCheckpoint.hpp"
#define main MembraneRegressionMain
#include "test_pretensioned_membrane.cpp"
#undef main

int main()
{
	try {
		const iga::PretensionedMembraneMaterial material{2.,.3,4.,1.};
		auto original=MakeSquare(2,material);const auto& layout=original.Layout();
		auto context=Context(layout,0.,.1,1,0,"checkpoint-first");
		std::vector<double> load(9,0.);load[4]=.2;
		auto first=original.SolveTrial(context,Traction(layout,context.expected_traction_stamp,load));
		original.FinalizeTrial(original.PrepareTrial(first));
		const auto bytes=iga::SerializeMembraneCheckpoint(original,1,.1);
		const auto identity=original.CommittedStateIdentitySha256();
		auto model=MakeSquare(2,material);
		const auto state=iga::ParseMembraneCheckpoint(bytes,model,identity,1,.1);
		auto restored=MakeSquare(2,material,state);
		assert(restored.CommittedStateIdentitySha256()==identity);
		assert(iga::SerializeMembraneCheckpoint(restored,1,.1)==bytes);
		const auto unchanged=model.CommittedStateIdentitySha256();
		Reject([&] { (void)iga::ParseMembraneCheckpoint(bytes,model,HashText("wrong-state"),1,.1); });
		Reject([&] { (void)iga::ParseMembraneCheckpoint(bytes,model,identity,2,.1); });
		Reject([&] { (void)iga::ParseMembraneCheckpoint(bytes,model,identity,1,.2); });
		Reject([&] { (void)iga::ParseMembraneCheckpoint(bytes,model,identity,1,.1,1); });
		Reject([&] { (void)iga::SerializeMembraneCheckpoint(original,1,.1,1); });
		for(std::size_t count:{std::size_t(0),bytes.size()/2,bytes.size()-1})
			Reject([&] { (void)iga::ParseMembraneCheckpoint(std::string_view(bytes).substr(0,count),model,identity,1,.1); });
		Reject([&] { (void)iga::ParseMembraneCheckpoint(bytes+"x",model,identity,1,.1); });
		auto wrong_model=MakeSquare(2,{2.,.3,5.,1.});
		Reject([&] { (void)iga::ParseMembraneCheckpoint(bytes,wrong_model,identity,1,.1); });
		auto corrupt=bytes;corrupt.back()^=1;
		Reject([&] { (void)iga::ParseMembraneCheckpoint(corrupt,model,identity,1,.1); });
		const auto body=bytes.size()-24*layout.owned_global_node_ids.size();
		auto wrong_node=bytes;wrong_node[body]^=1;
		Reject([&] { (void)iga::ParseMembraneCheckpoint(wrong_node,model,identity,1,.1); });
		for(bool infinite:{false,true}) {
			auto invalid=bytes;
			for(std::size_t byte=0;byte<8;++byte)invalid[body+8+byte]=0;
			// The first node is clamped: reject both finite displacement 1.0
			// and positive infinity before exposing any restore candidate.
			invalid[body+14]=static_cast<char>(0xf0);
			invalid[body+15]=static_cast<char>(infinite?0x7f:0x3f);
			Reject([&] { (void)iga::ParseMembraneCheckpoint(invalid,model,identity,1,.1); });
		}
		assert(model.CommittedStateIdentitySha256()==unchanged);
		context=Context(layout,.1,.15,2,1,"checkpoint-second");load[4]=-.1;
		const auto traction=Traction(layout,context.expected_traction_stamp,load);
		const auto continued=original.SolveTrial(context,traction),resumed=restored.SolveTrial(context,traction);
		assert(continued.state.displacement_m==resumed.state.displacement_m);
		assert(continued.state.velocity_m_per_s==resumed.state.velocity_m_per_s);
		original.FinalizeTrial(original.PrepareTrial(continued));restored.FinalizeTrial(restored.PrepareTrial(resumed));
		assert(original.CommittedStateIdentitySha256()==restored.CommittedStateIdentitySha256());
		std::cout<<"membrane checkpoint byte roundtrip, model/epoch/corruption rejection, unchanged candidate and exact next-step continuation passed\n";
	} catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
	return 0;
}
