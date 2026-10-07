#ifndef IGA_NATIVE_TET_TEMPORAL_OUTPUT_HPP
#define IGA_NATIVE_TET_TEMPORAL_OUTPUT_HPP

#include "ParallelVtkOutput.hpp"
#include "TemporalUnstructuredVtkHdf.hpp"
#include <memory>

namespace iga {

inline bool NativeTetUseVtkHdf(const std::string& format)
{
	if(format=="auto"||format=="vtkhdf")return true;
	if(format=="pvtu")return false;
	throw std::invalid_argument("native visualization format must be auto, vtkhdf or pvtu");
}

// Native Tet mesh/state are already replicated. Build the complete grid only
// on root for HDF5; legacy PVTU builds owned cells on every rank. No field gather.
class NativeTetTemporalOutput
{
public:
	NativeTetTemporalOutput(MPI_Comm comm,std::filesystem::path directory,
		std::string stem,const std::string& format,bool moving=false,
		std::int64_t resume_steps=0,double resume_time=0.)
		: comm_(comm),directory_(std::move(directory)),stem_(std::move(stem)),
		  hdf_(NativeTetUseVtkHdf(format)),moving_(moving),resume_steps_(resume_steps),resume_time_(resume_time)
	{
		MPI_Comm_rank(comm_,&rank_);MPI_Comm_size(comm_,&ranks_);
		RequireCollectiveSameText(comm_,"native output configuration",
			directory_.string()+"\n"+stem_+"\n"+format+"\n"
			+std::to_string(moving_)+"\n"+std::to_string(resume_steps_));
		CollectiveLocalStage(comm_,"native PVTU series resume",[&]{
			if(hdf_||directory_.empty()||!resume_steps_)return;
			const auto path=directory_/(stem_+".pvd");
			if(!std::filesystem::exists(path))return;
			series_=ReadParallelVtkSeries(path);
			if(series_.size()>static_cast<std::uint64_t>(resume_steps_)
				||std::abs(series_.back().first-resume_time_)>1e-12)
				throw std::runtime_error("native PVTU checkpoint time or step count differs");
		});
	}
	template<class Builder>
	void Append(double time,const std::filesystem::path& legacy_directory,Builder&& build)
	{
		if(directory_.empty())return;
		std::ostringstream clock;clock<<std::setprecision(17)<<time;
		RequireCollectiveSameText(comm_,"native output time",clock.str());
		VtkPartition piece;
		CollectiveLocalStage(comm_,"native temporal field output",[&]{
			if(closed_)throw std::runtime_error("native output is closed");
			if(hdf_){
				if(rank_!=0)return;
				piece=build(0,1);
				if(!writer_){
					const auto path=directory_/(stem_+".vtkhdf");
					const auto count=std::filesystem::exists(path)?resume_steps_:0;
					writer_=std::make_unique<TemporalUnstructuredVtkHdfWriter>(
						path,piece,moving_,4,count,count?0:resume_steps_);
					if(count&&std::abs(writer_->LastTime()-resume_time_)>1e-12)
						throw std::runtime_error("native VTKHDF checkpoint time differs");
				}
				writer_->Append(time,piece);
			}else piece=build(rank_,ranks_);
		});
		if(!hdf_){
			const auto snapshot=directory_/legacy_directory;
			WriteParallelVtkSnapshot(comm_,snapshot,piece,time);
			series_.push_back({time,snapshot/"snapshot.pvtu"});
			WriteParallelVtkSeries(comm_,directory_/(stem_+".pvd"),series_);
		}
	}
	void Close()
	{
		CollectiveLocalStage(comm_,"native temporal output close",[&]{
			closed_=true;
			if(writer_)writer_->Close();
		});
	}
private:
	MPI_Comm comm_;
	std::filesystem::path directory_;
	std::string stem_;
	bool hdf_=true,moving_=false,closed_=false;
	std::int64_t resume_steps_=0;
	double resume_time_=0.;
	int rank_=0,ranks_=1;
	std::unique_ptr<TemporalUnstructuredVtkHdfWriter> writer_;
	std::vector<std::pair<double,std::filesystem::path>> series_;
};
} // namespace iga
#endif
