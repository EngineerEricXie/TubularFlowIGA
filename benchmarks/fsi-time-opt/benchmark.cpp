// Same fixed Y-vessel operator workload as archived job 45981715.
#define main heartbeat_application_main
#include "driver.cpp"
#undef main
int main(int argc, char** argv)
{
    int provided = MPI_THREAD_SINGLE;
    if (MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided) != MPI_SUCCESS) return 2;
    if (provided < MPI_THREAD_FUNNELED) { MPI_Finalize(); return 2; }
    PetscInitialize(&argc,&argv,nullptr,nullptr);
    int status=0;
    try {
        std::cout << std::unitbuf;
        iga::CurrentPhaseProfile().EnableFromEnvironment();
        Require(argc>1,"surface input required");
        const auto map=ReadVessel(argv[1]);
        iga::MovingImmersedTransientFlowOptions options;
        options.grid={{{-.011,-.010,-.0035}},{{.011,.010,.0035}},{{16,16,5}}};
        options.geometry.volume.max_depth=1;
        options.geometry.volume.empty_rule_rescue_max_depth=3;
        options.geometry.volume.max_nodes=options.geometry.volume.max_leaves=options.geometry.volume.max_points=4000000;
        options.geometry.volume_storage=iga::CutCellVolumeQuadratureStorageMode::Compact;
        options.flow.parameters={1060.,.0035,1.}; options.flow.wall_labels={7};
        options.flow.wall_inertial_gamma0=1.;
        options.flow.ports={{"inlet",1,iga::ImmersedFlowPortControlMode::Pressure,20.},{"daughter_lower",2,iga::ImmersedFlowPortControlMode::Pressure,0.},{"daughter_upper",3,iga::ImmersedFlowPortControlMode::Pressure,0.}};
        options.flow.flow_controller_reference_flow_m3_s=1e-6;
        const auto geometry=iga::MovingCutGeometry::Build(options.grid,map.FullReference(),options.geometry);
        iga::ImmersedTransientFlowRuntime runtime(*geometry,options.flow);
        runtime.SetCommittedGlobalState(iga::ImmersedGlobalFlowState(0.,0,runtime.Layout(),std::vector<std::array<double,4>>(runtime.Layout().NodeIds().size()),{},false,0.));
        runtime.BeginTrial(.05,1,.05);
        // One warmup, then three samples per state, in one process (not six independent states).
        for(int sample=-1;sample<6;++sample) {
            const int state_id=sample<0?0:sample%2;
            std::vector<PetscScalar> state(runtime.Layout().Rows()),direction(state.size());
            for(std::size_t i=0;i<state.size();++i) {
                state[i]=state_id==0?0.:.01*std::sin(.13*(i+1));
                direction[i]=std::cos(.17*(i+1));
            }
            runtime.SetTrialState(state);
            const auto begin=std::chrono::steady_clock::now(); runtime.Assemble();
            const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
            iga::Sha256 hash;
            for(auto value:runtime.AssembledNegativeResidual()) hash.AppendNormalizedDouble(PetscRealPart(value));
            for(auto value:runtime.AssembledJacobianAction(direction)) hash.AppendNormalizedDouble(PetscRealPart(value));
            std::cout<<"benchmark_sample trial="<<sample<<" state="<<state_id<<" assembly_s="<<std::setprecision(17)<<elapsed<<" sha256="<<hash.Hex()<<" mallocs="<<runtime.JacobianStorageInfo().mallocs<<'\n';
        }
        runtime.AbortTrial(); iga::CurrentPhaseProfile().Write(std::cout,0,1,0);
    } catch(const std::exception& error) { std::cerr<<error.what()<<'\n'; status=1; }
    PetscFinalize(); MPI_Finalize(); return status;
}
