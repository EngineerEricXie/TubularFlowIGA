#define main BezierVisualizationRegressionMain
#include "test_bezier_visualization.cpp"
#undef main
#include <iostream>
#include <map>

int main()
{
	try {
		const auto first=MakeElement(0,0,0,true),second=MakeElement(1,1,0,true);
		std::map<std::string,std::int64_t> canonical;
		for(const auto& element:{first,second})for(std::size_t point=0;point<64;++point) {
			const auto signature=iga::BuildBezierPointSignature(element,point);
			const auto key=iga::EncodeBezierPointSignature(signature);
			const auto id=iga::BezierPointOccurrenceId(element.id,point);
			const auto found=canonical.find(key);
			if(found==canonical.end())canonical.emplace(key,id);else found->second=std::min(found->second,id);
		}
		assert(canonical.size()==112);
		// Reverse element order and extraction rows: the elected occurrence ids
		// and complete keys must not depend on input traversal or local row order.
		std::map<std::string,std::int64_t> reordered;
		for(auto element:{second,first}) {
			std::reverse(element.connectivity.begin(),element.connectivity.end());
			std::reverse(element.extraction.begin(),element.extraction.end());
			for(std::size_t point=0;point<64;++point) {
				const auto key=iga::EncodeBezierPointSignature(iga::BuildBezierPointSignature(element,point));
				const auto id=iga::BezierPointOccurrenceId(element.id,point);
				const auto found=reordered.find(key);
				if(found==reordered.end())reordered.emplace(key,id);else found->second=std::min(found->second,id);
			}
		}
		assert(reordered==canonical);
		const auto coincident=MakeElement(2,0,1000,false);
		assert(iga::EncodeBezierPointSignature(iga::BuildBezierPointSignature(first,0))
			!=iga::EncodeBezierPointSignature(iga::BuildBezierPointSignature(coincident,0)));
		const auto key=iga::EncodeBezierPointSignature(iga::BuildBezierPointSignature(first,0));
		assert(key.size()==20&&static_cast<unsigned char>(key[0])==1);
		const std::uint64_t coefficient=1000000000000ULL;
		for(unsigned byte=0;byte<8;++byte)assert(static_cast<unsigned char>(key[12+byte])==((coefficient>>(8*byte))&255u));
		int rejected=0;
		const auto reject=[&](auto function) { bool failed=false;try { function(); }catch(const std::exception&) { failed=true; }assert(failed);++rejected; };
		reject([&]{iga::BuildBezierPointSignature(first,64);});
		reject([&]{auto bad=first;bad.connectivity.pop_back();iga::BuildBezierPointSignature(bad,0);});
		reject([&]{auto bad=first;bad.extraction[0][0]=0;iga::BuildBezierPointSignature(bad,0);});
		reject([&]{auto bad=first;bad.extraction[0][0]=std::numeric_limits<double>::quiet_NaN();iga::BuildBezierPointSignature(bad,0);});
		reject([&]{iga::BezierPointOccurrenceId(std::numeric_limits<std::uint64_t>::max(),0);});
		reject([&]{iga::BezierPointOccurrenceId(0,64);});
		reject([&]{iga::EncodeBezierPointSignature({});});
		reject([&]{iga::BezierPointSignature bad;bad.key={{2,1},{1,1}};iga::EncodeBezierPointSignature(bad);});
		iga::BezierPointSignature negative;negative.key={{3,-2}};
		const auto signed_key=iga::EncodeBezierPointSignature(negative);
		assert(static_cast<unsigned char>(signed_key[12])==254);
		for(unsigned byte=13;byte<20;++byte)assert(static_cast<unsigned char>(signed_key[byte])==255);
		assert(iga::BezierPointOccurrenceId(static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max())/64,63)==std::numeric_limits<std::int64_t>::max());
		std::cout<<"bezier_point_signature_test: PASS 112 identities, 16 shared points, traversal parity, Int64 limit, "<<rejected<<" rejections\n";
	} catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
}
