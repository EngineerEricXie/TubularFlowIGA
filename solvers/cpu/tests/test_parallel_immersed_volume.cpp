#include "CompliantChannelFsiFixture.hpp"
#include <atomic>
#include <iostream>
#include <thread>

namespace {
void Check(bool value, const char* message)
{
	if(!value) throw std::runtime_error(message);
}

template<class Work> void Reject(Work work, const char* message)
{
	bool rejected=false;
	try { work(); }
	catch(const std::exception& error) { rejected=std::string(error.what()).find(message)!=std::string::npos; }
	Check(rejected,"missing expected failure diagnostic");
}
}

int main(int argc, char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int status=0;
	try {
		const auto caller=std::this_thread::get_id();
		auto options=iga::compliant_channel_fixture::FlowOptions();
		options.geometry.volume.max_depth=1;
		const auto initial=iga::compliant_channel_fixture::InitialMaterial();
		std::vector<int> thread_counts{1};
#ifdef _OPENMP
		thread_counts.insert(thread_counts.end(),{2,4});
#endif
		for(const auto storage:{iga::CutCellVolumeQuadratureStorageMode::Expanded,iga::CutCellVolumeQuadratureStorageMode::Compact}) {
		options.geometry.volume_storage=storage;
		const auto geometry=iga::MovingCutGeometry::Build(options.grid,initial,options.geometry);
		options.flow.body_force=[caller](const std::array<double,3>& p) {
			Check(std::this_thread::get_id()==caller,"body force callback escaped prepare thread");
			return std::array<double,3>{{.01*p[0],-.02*p[1],.03*p[2]}};
		};
		std::vector<PetscScalar> reference_residual, reference_action;
		std::string reference_input;
		for(const int threads:thread_counts) {
			Check(::setenv("IGA_ASSEMBLY_THREADS",std::to_string(threads).c_str(),1)==0,"setenv failed");
			Check(::setenv("IGA_ASSEMBLY_BATCH_SIZE","8",1)==0,"setenv failed");
			iga::ImmersedTransientFlowRuntime runtime(*geometry,options.flow);
			runtime.SetCommittedGlobalState(iga::ImmersedGlobalFlowState(0.0,0,runtime.Layout(),
				std::vector<std::array<double,4>>(runtime.Layout().NodeIds().size()),{},false,0.0));
			const auto committed=runtime.CommittedState();
			const auto committed_hash=runtime.CommittedGlobalState().HashSha256();
			runtime.BeginTrial(1.0,1,1.0);
			const auto seed=runtime.TrialState();
			const auto input=runtime.Diagnostics().input_hash_sha256;
			std::vector<PetscScalar> trial(seed.size()), direction(seed.size());
			for(std::size_t i=0;i<trial.size();++i) {
				trial[i]=.003*(static_cast<int>((7*i+3)%19)-9);
				direction[i]=.02*(static_cast<int>((11*i+1)%23)-11);
			}
			std::atomic<int> worker_calls{0}, completed{0};
			runtime.SetVolumeProbeForTesting([&](std::uint64_t) {
				++completed;
				if(std::this_thread::get_id()!=caller) ++worker_calls;
			});
			runtime.SetTrialState(trial);
			runtime.Assemble();
			const auto residual=runtime.AssembledNegativeResidual();
			const auto action=runtime.AssembledJacobianAction(direction);
			Check(completed==27,"fixture did not compute all 27 volume cells");
			Check(threads==1 ? worker_calls==0 : worker_calls>0,"actual volume worker execution differs");
			if(threads==1) { reference_residual=residual; reference_action=action; reference_input=input; }
			else Check(residual==reference_residual && action==reference_action && input==reference_input,
				"parallel residual, Jacobian action or physical input identity changed");
			// Fail after real integration in the second batch: the first batch
			// has already been scattered, and later retry must clear its values.
			completed=0;
			runtime.SetVolumeProbeForTesting([&](std::uint64_t cell) {
				++completed;
				if(cell==10 || cell==12) throw std::runtime_error("volume cell "+std::to_string(cell));
			});
			Reject([&] { runtime.SolveTrial(); },"volume cell 10");
			Check(completed==16,"failed batch did not join all workers or stop at its boundary");
			Check(runtime.TrialState()==seed && runtime.CommittedState()==committed
				&& runtime.CommittedGlobalState().HashSha256()==committed_hash
				&& runtime.Diagnostics().trial_active && !runtime.Diagnostics().converged,
				"worker failure did not preserve committed state and roll back trial");
			runtime.SetVolumeProbeForTesting({});
			runtime.SetTrialState(trial); runtime.Assemble();
			Check(runtime.AssembledNegativeResidual()==residual && runtime.AssembledJacobianAction(direction)==action,
				"retry retained a failed assembly contribution");
			// Also exercise abort followed by a new trial after pending insertion.
			runtime.SetVolumeProbeForTesting([](std::uint64_t cell) {
				if(cell==10) throw std::runtime_error("abort pending volume");
			});
			Reject([&] { runtime.Assemble(); },"abort pending volume");
			runtime.AbortTrial();
			runtime.SetVolumeProbeForTesting({});
			runtime.BeginTrial(1.0,1,1.0); runtime.SetTrialState(trial); runtime.Assemble();
			Check(runtime.AssembledNegativeResidual()==residual && runtime.AssembledJacobianAction(direction)==action,
				"abort/new trial retained a failed assembly contribution");
			bool rejected_caller=false;
			std::thread foreign([&] {
				try { runtime.Assemble(); }
				catch(const std::logic_error&) { rejected_caller=true; }
			});
			foreign.join();
			Check(rejected_caller,"foreign caller reached PETSc assembly");
			runtime.AbortTrial();
			std::cout << "parallel_immersed_volume storage=" << (storage==iga::CutCellVolumeQuadratureStorageMode::Expanded ? "expanded" : "compact") << " threads=" << threads << " cells=27 exact_fields=1 rollback_retry=1 passed\n";
		}
		}
	} catch(const std::exception& error) { std::cerr << error.what() << '\n'; status=1; }
	PetscFinalize();
	return status;
}
