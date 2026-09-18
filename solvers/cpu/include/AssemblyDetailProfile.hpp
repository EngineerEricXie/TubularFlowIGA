#ifndef IGA_ASSEMBLY_DETAIL_PROFILE_HPP
#define IGA_ASSEMBLY_DETAIL_PROFILE_HPP
#include <chrono>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <map>
#include <string>
#include <utility>
namespace iga {
inline bool DetailProfileEnabled() noexcept {
    const char* v=std::getenv("IGA_PROFILE_DETAIL");
    return v && v[0]=='1' && v[1]=='\0';
}
inline double DetailNow() noexcept {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
struct AssemblyCallContext {
    const char* reason="diagnostic";
    long long coupling_iteration=-1,newton_iteration=-1;
    double damping=0.;
};
inline AssemblyCallContext& CurrentAssemblyCallContext() {
    static thread_local AssemblyCallContext context; return context;
}
class AssemblyCallScope {
public:
    explicit AssemblyCallScope(const char* reason,long long newton=-1,double damping=0.,long long coupling=-1)
        : saved_(CurrentAssemblyCallContext()) {
        auto& c=CurrentAssemblyCallContext(); c.reason=reason; c.newton_iteration=newton; c.damping=damping;
        if(coupling>=0) c.coupling_iteration=coupling;
    }
    ~AssemblyCallScope() { CurrentAssemblyCallContext()=saved_; }
private: AssemblyCallContext saved_;
};
class DetailTimer {
public:
    explicit DetailTimer(double* value) noexcept : value_(value),start_(value?DetailNow():0.) {}
    DetailTimer(double& value,bool enabled) noexcept : DetailTimer(enabled?&value:nullptr) {}
    DetailTimer(const DetailTimer&)=delete;
    ~DetailTimer() { if(value_) *value_+=DetailNow()-start_; }
private: double* value_; double start_;
};
// Caller-only publication; instrumentation I/O failures never change solve semantics.
class DetailRecord {
public:
    explicit DetailRecord(const char* kind) : enabled_(DetailProfileEnabled()),kind_(kind),
        start_(enabled_?DetailNow():0.),exceptions_(std::uncaught_exceptions()) {}
    ~DetailRecord() { Finish(); }
    bool Enabled() const noexcept { return enabled_; }
    DetailTimer Time(const char* key) { return DetailTimer(enabled_?&seconds_[key]:nullptr); }
    void Number(const char* key,double value) { if(enabled_) numbers_[key]=value; }
    void String(const char* key,const std::string& value) { if(enabled_) strings_[key]=value; }
    void Finish() noexcept {
        if(!enabled_ || finished_) return;
        finished_=true;
        try {
            const double elapsed=DetailNow()-start_;
            static thread_local std::ofstream out("assembly-detail.jsonl",std::ios::app);
            const auto& c=CurrentAssemblyCallContext();
            out<<std::setprecision(17)<<"{\"schema_version\":1,\"kind\":\""<<kind_
                <<"\",\"status\":\""<<(std::uncaught_exceptions()>exceptions_?"failed":"completed")
                <<"\",\"reason\":\""<<c.reason<<"\",\"coupling_iteration\":"<<c.coupling_iteration
                <<",\"newton_iteration\":"<<c.newton_iteration<<",\"damping\":"<<c.damping<<",\"wall_s\":"<<elapsed;
            for(const auto& x:numbers_) out<<",\""<<x.first<<"\":"<<x.second;
            for(const auto& x:strings_) out<<",\""<<x.first<<"\":\""<<x.second<<'"';
            out<<",\"timers_s\":{"; bool first=true;
            for(const auto& x:seconds_) { if(!first) out<<','; first=false; out<<'"'<<x.first<<"\":"<<x.second; }
            out<<"}}\n";
        } catch(...) {}
    }
private:
    bool enabled_,finished_=false; const char* kind_; double start_; int exceptions_;
    std::map<std::string,double> numbers_,seconds_;
    std::map<std::string,std::string> strings_;
};
template<class T,class... Args>
T DetailConstruct(DetailRecord& record,const char* stage,Args&&... args) {
    auto timer=record.Time(stage); return T(std::forward<Args>(args)...);
}
template<class Work>
auto DetailEvaluate(DetailRecord& record,const char* stage,Work work) {
    auto timer=record.Time(stage); return work();
}
} // namespace iga
#endif
