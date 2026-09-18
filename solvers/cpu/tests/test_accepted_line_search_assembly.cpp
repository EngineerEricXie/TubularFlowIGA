// Additional branch fixtures use the existing depth-3 branch geometry and
// original Options()/fields; they do not replace the original depth-4 suite.
#pragma push_macro("main")
#undef main
#define main p1b_original_immersed_suite_main
#include "test_immersed_transient_flow.cpp"
#undef main
#pragma pop_macro("main")

namespace {
struct ReuseResult {
    std::vector<PetscScalar> state, residual, action;
    iga::ImmersedTransientFlowDiagnostics diagnostics;
};
bool SamePhysical(const ReuseResult& a,const ReuseResult& b) {
    return a.state==b.state && a.residual==b.residual && a.action==b.action
        && SameNewtonRecords(a.diagnostics.newton_steps,b.diagnostics.newton_steps)
        && SamePortsBitwise(a.diagnostics.ports,b.diagnostics.ports)
        && SameBits(a.diagnostics.residual_norm,b.diagnostics.residual_norm)
        && SameBits(a.diagnostics.true_linear_relative_residual,b.diagnostics.true_linear_relative_residual)
        && SameBits(a.diagnostics.pressure_measure,b.diagnostics.pressure_measure)
        && SameBits(a.diagnostics.pressure_gauge_defect,b.diagnostics.pressure_gauge_defect)
        && a.diagnostics.nonlinear_iterations==b.diagnostics.nonlinear_iterations
        && a.diagnostics.ksp_iterations==b.diagnostics.ksp_iterations
        && a.diagnostics.ksp_reason==b.diagnostics.ksp_reason
        && a.diagnostics.input_hash_sha256==b.diagnostics.input_hash_sha256
        && a.diagnostics.solved_state_hash_sha256==b.diagnostics.solved_state_hash_sha256;
}
ReuseResult Capture(iga::ImmersedTransientFlowRuntime& runtime,bool capture_action=true) {
    ReuseResult result{runtime.TrialState(),runtime.AssembledNegativeResidual(),{},runtime.Diagnostics()};
    if(capture_action) result.action=runtime.AssembledJacobianAction(Direction(runtime));
    return result;
}
ReuseResult CheckCase(const iga::MovingCutGeometry& geometry,bool reuse,const std::string& name,bool residual_only=false) {
    auto options=Options(); options.reuse_accepted_line_search_assembly=reuse;
    options.line_search_residual_only=residual_only;
    if(name=="one-step") options.nonlinear_relative_tolerance=.99;
    iga::ImmersedTransientFlowRuntime runtime(geometry,options);
    const bool zero=name=="zero";
    runtime.SetCommittedGlobalState(iga::ImmersedGlobalFlowState(0.0,0,runtime.Layout(),
        zero?std::vector<std::array<double,4>>(runtime.Layout().NodeIds().size()):NonconstantFields(runtime.Layout()),
        zero?std::vector<double>{0.,0.}:std::vector<double>{.019,-.023},true,zero?0.:.031));
    runtime.BeginTrial(1.,1,1.); const auto seed=runtime.TrialState();
    const auto committed=runtime.CommittedState();
    bool hook_used=false;
    if(name=="half" || name=="nan-first" || name=="failure" || name=="nan-all") {
        runtime.SetLineSearchNormProbeForTesting([&](PetscInt it,double damping,double& candidate,double norm) {
            if(name=="failure") { candidate=norm; hook_used=true; }
            else if(name=="nan-all") { candidate=std::numeric_limits<double>::quiet_NaN(); hook_used=true; }
            else if(it==0 && damping==1.) { candidate=name=="half"?norm:std::numeric_limits<double>::quiet_NaN(); hook_used=true; }
        });
    }
    if(name=="throw-ready") runtime.SetNewtonBoundaryProbeForTesting([&](PetscInt) {
        hook_used=true; throw std::runtime_error("injected ready-boundary failure");
    });
    const bool failure=name=="failure" || name=="nan-all" || name=="throw-ready";
    if(failure) {
        RejectWithMessage([&] { runtime.SolveTrial(); },name=="throw-ready"?
            "injected ready-boundary failure":"immersed transient backtracking failed");
        assert(hook_used && runtime.TrialState()==seed && runtime.CommittedState()==committed);
        assert(runtime.Diagnostics().trial_active && !runtime.Diagnostics().converged);
        runtime.SetLineSearchNormProbeForTesting({}); runtime.SetNewtonBoundaryProbeForTesting({});
    }
    assert(runtime.SolveTrial()); auto result=Capture(runtime,!residual_only);
    if(zero) assert(result.diagnostics.newton_steps.empty() && result.diagnostics.ksp_iterations==0 && result.diagnostics.attempt_assembly_count==1);
    else if(name=="one-step") assert(result.diagnostics.newton_steps.size()==1);
    else assert(result.diagnostics.newton_steps.size()>1);
    if(name=="half" || name=="nan-first") assert(hook_used && result.diagnostics.newton_steps.front().damping==.5);
    // A residual-only accepted convergence intentionally leaves J stale. Build
    // a fresh full operator at the same state, prove R exact, and use that J for
    // the baseline/candidate comparison.
    runtime.Assemble(); assert(runtime.AssembledNegativeResidual()==result.residual);
    const auto fresh_action=runtime.AssembledJacobianAction(Direction(runtime));
    if(!residual_only) assert(fresh_action==result.action); result.action=fresh_action;
    runtime.Rollback(); assert(runtime.TrialState()==seed);
	const auto cache_hits_before_fault=runtime.Diagnostics().volume_basis_cache_hits;
	const auto cache_misses_before_fault=runtime.Diagnostics().volume_basis_cache_misses;
    { iga::ImmersedTransientFlowRuntime::FirstAssemblyFailureScopeForTesting fault;
      RejectWithMessage([&] { runtime.SolveTrial(); },"injected first immersed transient assembly failure");
      assert(fault.Consumed() && runtime.TrialState()==seed
		  && runtime.Diagnostics().volume_basis_cache_hits==cache_hits_before_fault
		  && runtime.Diagnostics().volume_basis_cache_misses==cache_misses_before_fault); }
    assert(runtime.SolveTrial()); auto retry=Capture(runtime,!residual_only); runtime.Assemble();
	assert(runtime.Diagnostics().volume_basis_cache_enabled
		&& runtime.Diagnostics().volume_basis_cache_hits>cache_hits_before_fault
		&& runtime.Diagnostics().volume_basis_cache_misses==cache_misses_before_fault);
    assert(runtime.AssembledNegativeResidual()==retry.residual); const auto retry_fresh_action=runtime.AssembledJacobianAction(Direction(runtime));
    if(!residual_only) assert(retry_fresh_action==retry.action); retry.action=retry_fresh_action;
    assert(SamePhysical(result,retry) && result.diagnostics.attempt_assembly_count==retry.diagnostics.attempt_assembly_count);
    // A state setter after convergence cannot retain readiness across calls.
    auto invalid=runtime.TrialState(); invalid.front()=std::numeric_limits<double>::quiet_NaN();
    Reject([&] { runtime.SetTrialState(invalid); }); assert(runtime.TrialState()==retry.state);
    runtime.AbortTrial(); assert(runtime.CommittedState()==committed && runtime.Diagnostics().idle);
    std::cout << "accepted_full_case " << name << " reuse=" << reuse
        << " steps=" << result.diagnostics.newton_steps.size()
        << " assemblies=" << result.diagnostics.attempt_assembly_count << " exact_retry_operator=1 passed\n";
    return result;
}
}
int main(int argc,char** argv) {
    std::string selected;
    bool p3=false;
    for(int i=1;i<argc;++i) if(std::string(argv[i])=="--reuse-case") {
        if(++i==argc) return 2;
        selected=argv[i];
    }
    for(int i=1;i<argc;++i) if(std::string(argv[i])=="--p3-residual-only") p3=true;
    if(!selected.empty() && selected!="zero" && selected!="one-step" && selected!="normal"
        && selected!="half" && selected!="nan-first" && selected!="failure"
        && selected!="nan-all" && selected!="throw-ready") return 2;
    PetscInitialize(&argc,&argv,nullptr,nullptr); int status=0;
    iga::CurrentPhaseProfile().EnableFromEnvironment(); std::cout << std::unitbuf;
    try {
        const auto soup=Cube(); iga::PrescribedSurfaceMotion motion({{0.,soup},{1.,soup}});
        iga::MovingCutGeometryOptions options; options.volume.max_depth=p3?2:3;
        options.volume.max_nodes=500000; options.volume.max_leaves=500000; options.volume.max_points=3000000;
        const iga::CubicCartesianGridSpec grid{{{0,0,0}},{{1,1,1}},{{3,3,3}}};
        const auto geometry=iga::MovingCutGeometry::Build(grid,motion.Evaluate(1.,0.,1.),options);
        for(const std::string name:{"zero","one-step","normal","half","nan-first","failure","nan-all","throw-ready"}) {
            if(!selected.empty() && selected!=name) continue;
            if(p3) {
                const auto full=CheckCase(*geometry,true,name,false),residual=CheckCase(*geometry,true,name,true);
                assert(SamePhysical(full,residual));
                const auto updates=residual.diagnostics.newton_steps.size();
                assert(residual.diagnostics.attempt_full_assembly_count==(updates?updates:1));
                assert(residual.diagnostics.attempt_residual_only_count>=updates);
                assert(full.diagnostics.attempt_residual_only_count==0 && (updates==0
                    ? full.diagnostics.attempt_full_assembly_count==residual.diagnostics.attempt_full_assembly_count
                    : full.diagnostics.attempt_full_assembly_count>residual.diagnostics.attempt_full_assembly_count));
                std::cout << "residual_only_pair " << name << " full_before=" << full.diagnostics.attempt_full_assembly_count
                    << " full_after=" << residual.diagnostics.attempt_full_assembly_count << " residual_calls=" << residual.diagnostics.attempt_residual_only_count << " bitwise=1 passed\n";
                continue;
            }
            const auto off=CheckCase(*geometry,false,name),on=CheckCase(*geometry,true,name);
            assert(SamePhysical(off,on));
            const auto saved=off.diagnostics.attempt_assembly_count-on.diagnostics.attempt_assembly_count;
            assert(saved==(off.diagnostics.newton_steps.empty()?0:off.diagnostics.newton_steps.size()-1));
            std::cout << "accepted_full_pair " << name << " removed=" << saved << " bitwise=1 passed\n";
        }
    } catch(const std::exception& error) { std::cerr << error.what() << '\n'; status=1; }
    int rank=0,ranks=1; MPI_Comm_rank(PETSC_COMM_WORLD,&rank); MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
    iga::CurrentPhaseProfile().Write(std::cout,rank,ranks,status); PetscFinalize(); return status;
}
