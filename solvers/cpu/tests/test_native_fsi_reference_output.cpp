#include "NativeFsiReferenceOutput.hpp"
#include <cassert>
#include <limits>

template<class Work> void Reject(Work work) {
    bool rejected=false; try { work(); } catch(const std::exception&) { rejected=true; } assert(rejected);
}
int main(int argc,char** argv) {
    assert(argc==2); const std::filesystem::path root=argv[1];
    std::filesystem::create_directory(root);
    const std::vector<int> ids{2,7};
    iga::NativeFsiReferenceOutput good(root/"good");
    good.Add("fluid_velocity","m/s",ids,3,[](std::size_t r,std::size_t c){return .1*(r+1)*(c+1);});
    good.Finish(.05,1,false);
    Reject([&]{good.Finish(.05,1,true);});
    Reject([&]{good.Add("pressure","pa",ids,1,[](auto,auto){return 0.;});});
    Reject([&]{iga::NativeFsiReferenceOutput duplicate(root/"good");});
    iga::NativeFsiReferenceOutput positive(root/"positive");
    positive.Add("fluid_pressure","pa",ids,1,[](auto,auto){return 0.;}); positive.Finish(.05,1,true);
    iga::NativeFsiReferenceOutput unsorted(root/"unsorted");
    Reject([&]{unsorted.Add("pressure","pa",std::vector<int>{7,2},1,[](auto,auto){return 0.;});});
    Reject([&]{unsorted.Finish(.05,1,true);});
    assert(!std::filesystem::exists(root/"unsorted/manifest.json"));
    iga::NativeFsiReferenceOutput nonfinite(root/"nonfinite");
    Reject([&]{nonfinite.Add("pressure","pa",ids,1,[](auto,auto){return std::numeric_limits<double>::quiet_NaN();});});
    Reject([&]{nonfinite.Finish(.05,1,true);});
    assert(!std::filesystem::exists(root/"nonfinite/manifest.json"));
    iga::NativeFsiReferenceOutput negative(root/"negative");
    Reject([&]{negative.Add("pressure","pa",std::vector<int>{-1},1,[](auto,auto){return 0.;});});
    iga::NativeFsiReferenceOutput duplicate_field(root/"duplicate_field");
    duplicate_field.Add("pressure","pa",ids,1,[](auto,auto){return 0.;});
    Reject([&]{duplicate_field.Add("pressure","pa",ids,1,[](auto,auto){return 0.;});});
    Reject([&]{duplicate_field.Finish(.05,1,true);});
    iga::NativeFsiReferenceOutput incomplete(root/"incomplete");
    Reject([&]{incomplete.Finish(.05,1,true);});
    return 0;
}
