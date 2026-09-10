#define main MovingRuntimeRegressionMain
#include "test_moving_immersed_transient_flow.cpp"
#undef main

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int status = 0;
	try {
		PetscOptionsSetValue(nullptr, "-ksp_type", "unavailable_root_solver");
		PetscOptionsSetValue(nullptr, "-immersed_transient_ksp_type", "fgmres");
		PetscOptionsSetValue(nullptr, "-DOMAIN_RIGHT_FLOW_KSP_TYPE", "gmres");
		PetscBool aligned = PETSC_FALSE;
		Check(PetscOptionsHasName(nullptr, nullptr, "-test_solver_aligned", &aligned)==0, "aligned fixture option lookup failed");
		const auto surface = aligned ? Cube(0.,1.) : Cube(.15,.75);
		iga::PrescribedSurfaceMotion motion({{0., surface}, {2., surface}});
		auto left_options = Options(1e-5);
		// Small forced fields need a nonlinear residual tighter than the default
		// absolute stopping criterion to compare distinct Krylov strategies.
		left_options.flow.nonlinear_absolute_tolerance = 1e-14;
		left_options.flow.nonlinear_relative_tolerance = 1e-11;
		if (aligned) left_options.grid = {{{0,0,0}},{{1,1,1}},{{4,1,1}}};
		left_options.flow.solver_options_prefix = "domain_left_flow_";
		auto right_options = left_options; right_options.flow.solver_options_prefix = "domain_right_flow_";
		MovingImmersedTransientFlowRuntime left(motion.Evaluate(0.,0.,1.), left_options);
		MovingImmersedTransientFlowRuntime right(motion.Evaluate(0.,0.,1.), right_options);
		PetscOptionsSetValue(nullptr, "-immersed_transient_ksp_type", "unavailable_family_solver");
		PetscOptionsSetValue(nullptr, "-domain_right_flow_ksp_type", "unavailable_domain_solver");
		for (std::uint64_t step = 1; step <= 2; ++step) {
			left.BeginTrial(motion.Evaluate(step,step-1,step), step, 1.);
			right.BeginTrial(motion.Evaluate(step,step-1,step), step, 1.);
			Check(left.SolverConfiguration().ksp=="fgmres" && right.SolverConfiguration().ksp=="gmres", "epoch recaptured changed global options");
			Check(left.SolveTrial() && right.SolveTrial(), "configured moving solve failed");
			Check(left.TrialDiagnostics().input_hash_sha256!=right.TrialDiagnostics().input_hash_sha256, "input identity omitted solver options");
			const auto a = left.TrialState(), b = right.TrialState();
			Check(a.size()==b.size(), "configured moving layouts differ");
			double error = 0., norm = 0.;
			for (std::size_t i = 0; i < a.size(); ++i) { error += std::pow(PetscRealPart(a[i]-b[i]),2); norm += std::pow(PetscRealPart(a[i]),2); }
			const double relative = std::sqrt(error/std::max(norm,1e-300));
			std::cout << "moving_solver_comparison step=" << step << " relative_l2=" << relative << " reference_l2=" << std::sqrt(norm) << " difference_l2=" << std::sqrt(error)
				<< " left_residual=" << left.TrialDiagnostics().residual_norm << " right_residual=" << right.TrialDiagnostics().residual_norm << std::endl;
			Check(relative<1e-8, "moving independent solver fields differ");
			left.PrepareCommit(); right.PrepareCommit(); left.FinalizeCommit(); right.FinalizeCommit();
			std::cout << "moving_solver_options step=" << step << " relative_l2=" << relative << " passed\n";
		}
	} catch (const std::exception& error) { std::cerr << error.what() << '\n'; status = 1; }
	PetscFinalize(); return status;
}
