#include "MaterialSurfaceCheckpoint.hpp"
#include "CompliantChannelFsiFixture.hpp"
#include <iostream>

namespace {
void Require(bool value) { if(!value)throw std::runtime_error("material checkpoint regression failed"); }
template<class Function> void Reject(Function&& function)
{
	bool rejected=false;try { function(); } catch(const std::exception&) { rejected=true; }
	Require(rejected);
}
}

int main()
{
	try {
		const auto initial=iga::compliant_channel_fixture::InitialMaterial();
		auto positions=initial.SourceVerticesM();auto velocities=initial.SourceVertexVelocitiesMPerS();
		positions[4][2]+=.006;velocities[4][2]=.006;
		std::vector<iga::RawSurfaceTriangle> topology;
		for(const auto& source:initial.SourceTriangles()) {
			iga::RawSurfaceTriangle triangle;triangle.boundary_id=source.boundary_id;
			for(int axis=0;axis<3;++axis)triangle.indices[axis]=source.source_vertex_indices[axis];
			topology.push_back(triangle);
		}
		const auto target=iga::MaterialSurfaceKinematics::CreateFromSourceTopology(initial.ReferenceMaterialVerticesM(),
			positions,velocities,topology,2.,1.,2.);
		for(const auto* material:{&initial,&target}) {
			const auto bytes=iga::SerializeMaterialSurfaceCheckpoint(*material);
			const auto restored=iga::ParseMaterialSurfaceCheckpoint(bytes,material->ContentIdentitySha256());
			Require(restored.ContentIdentitySha256()==material->ContentIdentitySha256());
			Require(restored.MaterialIdentitySha256()==material->MaterialIdentitySha256());
			Require(restored.TopologyIdentitySha256()==material->TopologyIdentitySha256());
			Require(restored.SourceTriangles()==material->SourceTriangles());
			Require(restored.CanonicalTriangleProvenance()==material->CanonicalTriangleProvenance());
			Require(iga::SerializeMaterialSurfaceCheckpoint(restored)==bytes);
			for(std::size_t count:{std::size_t(0),std::size_t(8),bytes.size()/2,bytes.size()-1})
				Reject([&] { (void)iga::ParseMaterialSurfaceCheckpoint(std::string_view(bytes).substr(0,count),material->ContentIdentitySha256()); });
			Reject([&] { (void)iga::ParseMaterialSurfaceCheckpoint(bytes+"x",material->ContentIdentitySha256()); });
			Reject([&] { (void)iga::ParseMaterialSurfaceCheckpoint(bytes,std::string(64,'0')); });
			Reject([&] { (void)iga::ParseMaterialSurfaceCheckpoint(bytes,material->ContentIdentitySha256(),{1,1}); });
			Reject([&] { (void)iga::SerializeMaterialSurfaceCheckpoint(*material,{1,1}); });
			auto corrupt=bytes;corrupt[corrupt.size()-8]^=1;
			Reject([&] { (void)iga::ParseMaterialSurfaceCheckpoint(corrupt,material->ContentIdentitySha256()); });
			iga::checkpoint_metadata::Writer forged;forged.Text("IGA_MATERIAL_SURFACE/1");forged.Text(material->ContentIdentitySha256());
			forged.Real(2.);forged.Real(1.);forged.Real(2.);forged.Unsigned(100000);forged.Unsigned(200000);
			Reject([&] { (void)iga::ParseMaterialSurfaceCheckpoint(forged.Bytes(),material->ContentIdentitySha256()); });
			auto nonfinite=bytes;const auto offset=forged.Bytes().size();
			for(std::size_t byte=0;byte<8;++byte)nonfinite[offset+byte]=0;
			nonfinite[offset+6]=static_cast<char>(0xf0);nonfinite[offset+7]=static_cast<char>(0x7f);
			Reject([&] { (void)iga::ParseMaterialSurfaceCheckpoint(nonfinite,material->ContentIdentitySha256()); });
		}
		std::cout<<"material checkpoint reference/deformed byte roundtrip, provenance, corruption and caps passed\n";
	} catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
	return 0;
}
