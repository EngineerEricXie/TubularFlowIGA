#include "DistributedPointIdentity.hpp"
#include "DistributedPointValues.hpp"
#define main SerialBezierRegressionMain
#include "test_bezier_visualization.cpp"
#undef main
#include <iostream>
#include <map>

namespace {
void RequireIdentity(bool value,const char* message)
{
	if(!value)throw std::runtime_error(message);
}
void RunIdentities(MPI_Comm comm)
{
	int rank=0,ranks=0;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&ranks);
	std::vector<iga::PointIdentityOccurrence> global;
	for(const auto& element:{MakeElement(0,0,0,true),MakeElement(1,1,0,true)})
		for(std::size_t point=0;point<64;++point)global.push_back({iga::BezierPointOccurrenceId(element.id,point),
			iga::EncodeBezierPointSignature(iga::BuildBezierPointSignature(element,point))});
	// Binary keys, embedded NULs, and Int64 ids beyond exact double integers.
	global.push_back({9007199254740993LL,std::string("a\0b",3)});
	global.push_back({9007199254740994LL,std::string("a\0b",3)});
	global.push_back({9007199254740995LL,std::string("a\0c",3)});
	for(int mode=0;mode<3;++mode) {
		const auto owner=[&](std::size_t i) { return mode==0?0:mode==1?static_cast<int>(i%ranks):ranks-1-static_cast<int>(i%ranks); };
		std::map<std::string,std::pair<std::int64_t,int>> oracle;
		std::vector<iga::PointIdentityOccurrence> local;
		for(std::size_t i=0;i<global.size();++i) {
			const auto& occurrence=global[i];const auto found=oracle.find(occurrence.key);
			if(found==oracle.end()||occurrence.occurrence<found->second.first)oracle[occurrence.key]={occurrence.occurrence,owner(i)};
			if(owner(i)==rank)local.push_back(occurrence);
		}
		std::reverse(local.begin(),local.end());
		const auto result=iga::ResolveDistributedPointIdentities(comm,local);
		std::vector<std::int64_t> ids;std::vector<double> values;
		for(const auto& occurrence:local) {
			ids.push_back(occurrence.occurrence);
			values.push_back(static_cast<double>(occurrence.occurrence));
			values.push_back(-0.0);values.push_back(static_cast<double>(rank)+.25);
		}
		const auto transferred=iga::ExchangePointRepresentativeValues(comm,ids,result,values,3);
		iga::CollectiveLocalStage(comm,"identity oracle comparison",[&] {
			RequireIdentity(oracle.size()==114&&result.size()==local.size(),"incorrect unique identity count");
			for(std::size_t i=0;i<local.size();++i) {
				const auto expected=oracle.at(local[i].key);
				RequireIdentity(result[i].occurrence==expected.first&&result[i].rank==expected.second,"partition-dependent point representative");
				RequireIdentity(transferred[3*i]==static_cast<double>(expected.first)
					&&transferred[3*i+1]==0&&std::signbit(transferred[3*i+1])
					&&transferred[3*i+2]==expected.second+.25,"incorrect representative values or lost signed zero");
			}
		});
	}
	const auto empty=iga::ResolveDistributedPointIdentities(comm,{});
	iga::CollectiveLocalStage(comm,"empty identity test",[&] { RequireIdentity(empty.empty(),"all-empty identity result is nonempty"); });
	int failures=0;
	for(const auto& mode:std::vector<std::string>{"duplicate","empty-key","wire","local-cap","receive-cap"}) {
		if(mode=="receive-cap"&&ranks==1)continue;
		std::vector<iga::PointIdentityOccurrence> local{{rank,"shared"}};iga::PointIdentityLimits limits;
		if(mode=="duplicate"&&rank==ranks-1)local.push_back({0,"different-key"});
		if(mode=="empty-key"&&rank==ranks-1)local[0].key.clear();
		if(mode=="wire") { local[0].key=std::string(100,'x');limits.max_wire_bytes=64; }
		if(mode=="local-cap"&&rank==ranks-1)limits.max_local_occurrences=0;
		if(mode=="receive-cap")limits.max_wire_bytes=30;
		bool rejected=false;std::string message;
		try { iga::ResolveDistributedPointIdentities(comm,local,limits); }catch(const std::exception& error) { rejected=true;message=error.what(); }
		iga::CollectiveLocalStage(comm,"identity expected rejection",[&] {
			RequireIdentity(rejected,"invalid identity exchange accepted");
			const auto stage=mode=="duplicate"?"occurrence uniqueness":mode=="wire"?"key routing":mode=="receive-cap"?"receive allocation":"input validation";
			RequireIdentity(message.find(stage)!=std::string::npos,"identity failure occurred at an unexpected stage");
		});
		++failures;
		const auto retry=iga::ResolveDistributedPointIdentities(comm,{{rank,"shared"}});
		iga::CollectiveLocalStage(comm,"identity retry",[&] { RequireIdentity(retry.size()==1&&retry[0].occurrence==0&&retry[0].rank==0,"identity retry failed"); });
	}
	RequireIdentity(iga::ExchangePointRepresentativeValues(comm,{}, {}, {},3).empty(),"nonempty all-empty value exchange");
	int value_failures=0;
	for(const auto& mode:std::vector<std::string>{"owner","missing","nonfinite","components","receive-cap","response-cap"}) {
		if(ranks==1&&(mode=="components"||mode=="receive-cap"))continue;
		std::vector<iga::PointIdentityRepresentative> representative{{0,0}};
		std::vector<double> values{rank+.25,rank+.5};int components=2;iga::PointIdentityLimits limits;
		if(rank==ranks-1) {
			if(mode=="owner")representative[0].rank=ranks;
			if(mode=="missing")representative[0].occurrence=999;
			if(mode=="nonfinite")values[0]=std::numeric_limits<double>::quiet_NaN();
			if(mode=="components") { components=3;values.push_back(0); }
		}
		if(mode=="receive-cap")limits.max_wire_bytes=24;
		if(mode=="response-cap")limits.max_wire_bytes=20*static_cast<std::size_t>(ranks);
		const auto saved=values;bool rejected=false;std::string message;
		try { iga::ExchangePointRepresentativeValues(comm,{rank},representative,values,components,limits); }
		catch(const std::exception& error) { rejected=true;message=error.what(); }
		iga::CollectiveLocalStage(comm,"point value expected rejection",[&] {
			RequireIdentity(rejected,"invalid representative values accepted");
			const auto stage=(mode=="missing"||mode=="response-cap")?"response construction":mode=="components"?"component agreement":mode=="receive-cap"?"receive allocation":"input validation";
			RequireIdentity(message.find(stage)!=std::string::npos,"point values failed at unexpected stage");
			RequireIdentity(std::memcmp(saved.data(),values.data(),saved.size()*sizeof(double))==0,"failed exchange changed source values");
		});
		++value_failures;
		const auto retried=iga::ExchangePointRepresentativeValues(comm,{rank},{{0,0}},{rank+.25,rank+.5},2);
		iga::CollectiveLocalStage(comm,"point value retry",[&] {
			RequireIdentity(retried==std::vector<double>({.25,.5}),"representative values retry failed");
		});
	}
	if(rank==0)std::cout<<"distributed_point_values ranks="<<ranks<<" failures="<<value_failures<<" passed\n";
	if(rank==0)std::cout<<"distributed_point_identity ranks="<<ranks<<" layouts=3 keys=114 failures="<<failures<<" passed\n";
}
}
int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int rank=0,ranks=0;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	try {
		RequireIdentity(ranks==3,"require three ranks");RunIdentities(PETSC_COMM_WORLD);
		MPI_Comm group=MPI_COMM_NULL;MPI_Comm_split(PETSC_COMM_WORLD,rank==0?0:1,rank,&group);
		RunIdentities(group);MPI_Comm_free(&group);
	} catch(const std::exception& error) { std::cerr<<"rank "<<rank<<": "<<error.what()<<'\n';MPI_Abort(PETSC_COMM_WORLD,1);return 1; }
	PetscFinalize();return 0;
}
