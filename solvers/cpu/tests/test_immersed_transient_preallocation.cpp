#ifdef IGA_PREALLOCATION_BASELINE_HEADER
#include IGA_PREALLOCATION_BASELINE_HEADER
#endif
#include "CompliantChannelFsiFixture.hpp"
#include <cassert>
#include <iostream>

int main(int argc,char** argv)
{
    PetscInitialize(&argc,&argv,nullptr,nullptr);
    int status=0;
    try {
        const auto initial=iga::compliant_channel_fixture::InitialMaterial();
        const auto setup=iga::compliant_channel_fixture::FlowOptions();
        const auto geometry=iga::MovingCutGeometry::Build(setup.grid,initial,setup.geometry);
        for(int controlled=0;controlled<2;++controlled) {
            auto options=setup.flow;
            if(controlled) {
                options.ports={{"inlet",1,iga::ImmersedFlowPortControlMode::FlowRate,-.01},
                    {"outlet",2,iga::ImmersedFlowPortControlMode::FlowRate,.01}};
            }
            iga::ImmersedTransientFlowRuntime runtime(*geometry,options);
            runtime.SetCommittedGlobalState(iga::ImmersedGlobalFlowState(0.,0,runtime.Layout(),
                std::vector<std::array<double,4>>(runtime.Layout().NodeIds().size()),
                controlled?std::vector<double>{0.,0.}:std::vector<double>{},controlled!=0,0.));
            runtime.BeginTrial(1.,1,1.);
            iga::Sha256 digest;
            for(int trial=0;trial<2;++trial) {
                // A zero initial velocity must not leave holes needed later
                // by convection, ghost terms, flow controllers or the gauge.
                std::vector<PetscScalar> state(runtime.Layout().Rows(),0.);
                if(trial)for(std::size_t i=0;i<state.size();++i)state[i]=.01*std::sin(.13*(i+1));
                runtime.SetTrialState(state);runtime.Assemble();
#ifndef IGA_PREALLOCATION_BASELINE_HEADER
                const auto storage=runtime.JacobianStorageInfo();
                assert(storage.mallocs==0. && storage.nz_used>0. && storage.nz_unneeded==0.);
#endif
                std::vector<PetscScalar> direction(state.size());
                for(std::size_t i=0;i<direction.size();++i)direction[i]=std::cos(.17*(i+1));
                for(auto value:runtime.AssembledNegativeResidual())digest.AppendNormalizedDouble(PetscRealPart(value));
                for(auto value:runtime.AssembledJacobianAction(direction))digest.AppendNormalizedDouble(PetscRealPart(value));
            }
            runtime.AbortTrial();
            std::cout<<"controlled="<<controlled<<" residual_and_action_sha256="<<digest.Hex()<<'\n';
        }
    } catch(const std::exception& error) {std::cerr<<error.what()<<'\n';status=1;}
    PetscFinalize();return status;
}
