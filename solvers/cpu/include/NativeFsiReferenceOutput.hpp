#ifndef IGA_NATIVE_FSI_REFERENCE_OUTPUT_HPP
#define IGA_NATIVE_FSI_REFERENCE_OUTPUT_HPP

#include "Sha256.hpp"
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <locale>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace iga {

// Post-acceptance native coefficients, not a checkpoint. Uses the validated
// ReferenceStateOutput ID/precision/hash protocol, with CSV and truthful
// per-step gates (visualization-only steps must not claim a physics pass).
class NativeFsiReferenceOutput {
public:
    explicit NativeFsiReferenceOutput(std::filesystem::path directory)
        : directory_(std::move(directory)) {
        if(directory_.empty() || !std::filesystem::create_directory(directory_))
            throw std::runtime_error("native reference directory empty or already exists");
    }
    template<class Ids,class ValueAt>
    void Add(const std::string& name,const std::string& units,const Ids& ids,
             std::size_t columns,ValueAt value_at) {
        if(finished_ || failed_) throw std::logic_error("native reference closed or failed");
        failed_=true;
        if(name.empty() || name.find_first_not_of("abcdefghijklmnopqrstuvwxyz0123456789_")!=std::string::npos
           || units.empty() || units.find_first_not_of("abcdefghijklmnopqrstuvwxyz0123456789_/*^- ")!=std::string::npos
           || !names_.insert(name).second || ids.empty() || !columns)
            throw std::invalid_argument("invalid native reference metadata");
        std::ostringstream text; text.imbue(std::locale::classic());
        text<<std::scientific<<std::setprecision(17)<<"global_id";
        for(std::size_t c=0;c<columns;++c) text<<",c"<<c;
        text<<'\n';
        std::uint64_t previous=0;
        for(std::size_t row=0;row<ids.size();++row) {
            if constexpr(std::is_signed<typename Ids::value_type>::value)
                if(ids[row]<0) throw std::invalid_argument("negative native reference ID");
            const auto id=static_cast<std::uint64_t>(ids[row]);
            if(row && id<=previous) throw std::invalid_argument("native reference IDs not strictly increasing");
            previous=id; text<<id;
            for(std::size_t c=0;c<columns;++c) {
                const double value=value_at(row,c);
                if(!std::isfinite(value)) throw std::invalid_argument("nonfinite native reference value");
                text<<','<<value;
            }
            text<<'\n';
        }
        const auto contents=text.str(); Write(directory_/(name+".csv"),contents);
        Sha256 hash; hash.Append(contents.data(),contents.size());
        fields_.push_back({name,units,hash.Hex(),ids.size(),columns}); failed_=false;
    }
    void Finish(double time_s,std::uint64_t step,bool native_gates_passed) {
        if(finished_ || failed_ || fields_.empty() || !std::isfinite(time_s) || time_s<0.)
            throw std::logic_error("cannot publish incomplete native reference");
        failed_=true;
        std::ostringstream text; text.imbue(std::locale::classic());
        text<<std::setprecision(17)<<"{\n  \"schema_version\":1,\n  \"kind\":\"native_fsi_reference_state\",\n"
            <<"  \"case\":\"y-bifurcation\",\n  \"time_s\":"<<time_s<<",\n  \"step\":"<<step
            <<",\n  \"native_gates_passed\":"<<(native_gates_passed?"true":"false")<<",\n  \"checkpoint\":false,\n  \"fields\":[";
        for(std::size_t i=0;i<fields_.size();++i) {
            const auto& f=fields_[i];
            text<<(i?",":"")<<"\n    {\"name\":\""<<f.name<<"\",\"file\":\""<<f.name
                <<".csv\",\"units\":\""<<f.units<<"\",\"rows\":"<<f.rows<<",\"columns\":"<<f.columns
                <<",\"sha256\":\""<<f.sha256<<"\"}";
        }
        text<<"\n  ]\n}\n"; Write(directory_/"manifest.tmp",text.str());
        std::filesystem::rename(directory_/"manifest.tmp",directory_/"manifest.json");
        finished_=true; failed_=false;
    }
private:
    struct Field { std::string name,units,sha256; std::size_t rows,columns; };
    static void Write(const std::filesystem::path& path,const std::string& contents) {
        std::ofstream out(path,std::ios::binary);
        if(!out) throw std::runtime_error("cannot open native reference output");
        out.write(contents.data(),static_cast<std::streamsize>(contents.size())); out.close();
        if(!out) throw std::runtime_error("cannot complete native reference output");
    }
    std::filesystem::path directory_; std::set<std::string> names_; std::vector<Field> fields_;
    bool finished_=false,failed_=false;
};
} // namespace iga
#endif
