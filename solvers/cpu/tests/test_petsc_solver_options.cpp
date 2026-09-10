#include "PetscSolverOptions.hpp"
#include "ExecutionResources.hpp"
#include <iostream>

namespace {
void Check(PetscErrorCode code)
{
	if (code) throw std::runtime_error("test PETSc operation failed");
}
void Require(bool condition, const char* message)
{
	if (!condition) throw std::runtime_error(message);
}
PetscErrorCode CountErrors(MPI_Comm, int, const char*, const char*, PetscErrorCode code, PetscErrorType, const char*, void* context)
{
	++*static_cast<int*>(context);
	return code;
}
struct Options {
	Options() { Check(PetscOptionsCreate(&value)); }
	~Options() { PetscOptionsDestroy(&value); }
	void Set(const char* name, const char* value_text) { Check(PetscOptionsSetValue(value, name, value_text)); }
	PetscOptions value = nullptr;
};
void EnableView(Options& source, const std::string& family = {})
{
	PetscBool enabled = PETSC_FALSE;
	Check(PetscOptionsHasName(nullptr, nullptr, "-test_solver_view", &enabled));
	if (!enabled) return;
	for (const char* option : {"ksp_view", "ksp_converged_reason", "sub_ksp_converged_reason",
		"fieldsplit_first_ksp_converged_reason", "fieldsplit_second_ksp_converged_reason"})
		source.Set(("-"+family+option).c_str(), nullptr);
}
struct Linear {
	~Linear() { KSPDestroy(&solver); VecDestroy(&solution); VecDestroy(&rhs); MatDestroy(&matrix); }
	KSP solver = nullptr; Vec solution = nullptr, rhs = nullptr; Mat matrix = nullptr;
};
void Solve(MPI_Comm comm, iga::PetscSolverOptions& options, const char* type, const char* sub_pc, bool split = false, bool schur = false)
{
	int size = 0; MPI_Comm_size(comm, &size);
	Linear objects;
	Check(MatCreateAIJ(comm, 2, 2, 2*size, 2*size, 2, nullptr, 2, nullptr, &objects.matrix));
	PetscInt begin = 0, end = 0; Check(MatGetOwnershipRange(objects.matrix, &begin, &end));
	for (PetscInt row = begin; row < end; ++row) {
		Check(MatSetValue(objects.matrix, row, row, 4.0, INSERT_VALUES));
		if (row) Check(MatSetValue(objects.matrix, row, row-1, -1.0, INSERT_VALUES));
		if (row+1 < 2*size) Check(MatSetValue(objects.matrix, row, row+1, -1.0, INSERT_VALUES));
	}
	Check(MatAssemblyBegin(objects.matrix, MAT_FINAL_ASSEMBLY)); Check(MatAssemblyEnd(objects.matrix, MAT_FINAL_ASSEMBLY));
	Check(MatCreateVecs(objects.matrix, &objects.solution, &objects.rhs));
	Check(VecSet(objects.solution, 1.0)); Check(MatMult(objects.matrix, objects.solution, objects.rhs)); Check(VecSet(objects.solution, 0.0));
	Check(KSPCreate(comm, &objects.solver));
	// Deliberately create the PC before attachment, as production runtimes do.
	PC pc = nullptr; Check(KSPGetPC(objects.solver, &pc)); Check(PCSetType(pc, PCNONE));
	options.Attach(objects.solver);
	if (split) {
		Check(PCSetType(pc, PCFIELDSPLIT)); Check(PCFieldSplitSetBlockSize(pc, 2));
		const PetscInt first = 0, second = 1;
		Check(PCFieldSplitSetFields(pc, "first", 1, &first, &first));
		Check(PCFieldSplitSetFields(pc, "second", 1, &second, &second));
		Check(PCFieldSplitSetType(pc, schur ? PC_COMPOSITE_SCHUR : PC_COMPOSITE_ADDITIVE));
		if (schur) { Check(PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_FULL)); Check(PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_A11, nullptr)); }
	}
	options.Call("test solver options", [&] { return KSPSetFromOptions(objects.solver); });
	const char* actual = nullptr; Check(KSPGetType(objects.solver, &actual)); Require(std::string(actual)==type, "incorrect independent KSP override");
	Check(KSPSetOperators(objects.solver, objects.matrix, objects.matrix));
	iga::RequireKspFactorBackend(objects.solver, objects.matrix, comm);
	options.Call("test solve", [&] { return KSPSolve(objects.solver, objects.rhs, objects.solution); });
	KSPConvergedReason reason; PetscInt iterations = 0;
	Check(KSPGetConvergedReason(objects.solver, &reason)); Check(KSPGetIterationNumber(objects.solver, &iterations));
	Require(reason>0 && iterations>0, "solver did not converge");
	PetscInt local = 0; KSP* children = nullptr;
	if (split) Check(PCFieldSplitGetSubKSP(pc, &local, &children));
	else Check(PCBJacobiGetSubKSP(pc, &local, nullptr, &children));
	Require(local==(split ? 2 : 1), "unexpected child KSP layout");
	PC child = nullptr; Check(KSPGetPC(children[0], &child)); Check(PCGetType(child, &actual));
	if (std::string(actual)!=sub_pc) {
		const char* child_prefix = nullptr; PetscOptions child_options = nullptr;
		Check(KSPGetOptionsPrefix(children[0], &child_prefix)); Check(PetscObjectGetOptions(reinterpret_cast<PetscObject>(children[0]), &child_options));
		throw std::runtime_error(std::string("nested PC prefix inheritance failed: actual=")+actual+" expected="+sub_pc+" prefix="+(child_prefix ? child_prefix : "null")+" database="+(child_options==options.Database() ? "snapshot" : (child_options ? "other" : "global")));
	}
	if (split) {
		Check(KSPGetPC(children[1], &child)); Check(PCGetType(child, &actual));
		if (std::string(actual)!="lu") { const char* prefix = nullptr; Check(KSPGetOptionsPrefix(children[1], &prefix)); throw std::runtime_error(std::string("fieldsplit second child override failed: actual=")+actual+" prefix="+(prefix ? prefix : "null")); }
		MatSolverType backend = nullptr; Check(PCFactorGetMatSolverType(child, &backend));
		Require(backend && std::string(backend)=="mumps", "fieldsplit backend override failed");
		Check(PetscFree(children));
	}
	Check(VecShift(objects.solution, -1.0)); PetscReal norm = 0.0; Check(VecNorm(objects.solution, NORM_2, &norm));
	Require(norm<1e-10, "independently configured field differs from exact solution");
	options.RecordUsed();
	int rank = 0; MPI_Comm_rank(comm, &rank);
	std::cout << "solver_options rank=" << rank << " ranks=" << size << " prefix=" << options.Prefix()
		<< " ksp=" << type << " pc=" << (split ? "fieldsplit" : "bjacobi") << " sub_pc=" << sub_pc << " iterations=" << iterations << " error=" << norm << '\n';
}
void RunFamily(MPI_Comm comm)
{
	int rank = 0, size = 0; MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &size);
	Options source; EnableView(source, "immersed_static_");
	// Immersed families have never inherited the root KSP/PC defaults.
	source.Set("-ksp_type", "unavailable_solver"); source.Set("-pc_type", "none");
	source.Set("-IMMERSED_STATIC_ksp_type", "gmres"); source.Set("-immersed_static_ksp_rtol", "1e-12");
	source.Set("-immersed_static_pc_type", "bjacobi"); source.Set("-immersed_static_sub_ksp_type", "preonly");
	source.Set("-immersed_static_sub_pc_type", "lu");
	source.Set("-DOMAIN_RIGHT_FLOW_KSP_TYPE", "fgmres"); source.Set("-domain_right_flow_sub_pc_type", "jacobi");
	const auto before = iga::CapturePetscOptions(source.value);
	iga::PetscSolverOptions left(comm, "DOMAIN_LEFT_FLOW_", source.value, {}, {}, "immersed_static_", false);
	iga::PetscSolverOptions right(comm, "domain_right_flow_", source.value, {}, {}, "immersed_static_", false);
	source.Set("-immersed_static_ksp_type", "unavailable_solver");
	Solve(comm, left, "gmres", "lu"); Solve(comm, right, "fgmres", "jacobi");
	source.Set("-immersed_static_ksp_type", "gmres");
	Require(before==iga::CapturePetscOptions(source.value), "family snapshot changed source values");
	PetscBool used = PETSC_FALSE;
	Check(PetscOptionsUsed(source.value, "immersed_static_sub_pc_type", &used)); Require(used, "family child usage was not propagated");
	Check(PetscOptionsUsed(source.value, "domain_right_flow_ksp_type", &used)); Require(used, "domain usage was not propagated");
	Check(PetscOptionsUsed(source.value, "ksp_type", &used)); Require(!used, "family consumed ignored root KSP");
	int rejected = 0;
	try { iga::PetscSolverOptions bad(comm, "domain_left_flow_", source.value, {}, {}, rank==size-1 ? "bad" : "immersed_static_", false); }
	catch (const std::exception&) { ++rejected; }
	if (size>1) {
		try { iga::PetscSolverOptions bad(comm, "domain_left_flow_", source.value, {}, {}, "immersed_static_", rank==size-1); }
		catch (const std::exception&) { ++rejected; }
		try { iga::PetscSolverOptions bad(comm, "domain_left_flow_", source.value, {}, {}, rank==size-1 ? "immersed_transient_" : "immersed_static_", false); }
		catch (const std::exception&) { ++rejected; }
	}
	Require(rejected==(size>1 ? 3 : 1), "family policy disagreement was not rejected collectively");
}
void Run(MPI_Comm comm)
{
	int rank = 0, size = 0; MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &size);
	Options source; EnableView(source);
	source.Set("-ksp_type", "gmres"); source.Set("-ksp_rtol", "1e-12"); source.Set("-pc_type", "bjacobi");
	source.Set("-sub_ksp_type", "preonly"); source.Set("-sub_pc_type", "lu");
	source.Set("-left_ksp_type", "cg"); source.Set("-right_sub_pc_type", "jacobi");
	source.Set(("-"+std::string(506, 'z')).c_str(), "unrelated long option");
	source.Set("-unused_text", "spaces \"quotes\"\n-ksp_type cg "); source.Set("-flag", nullptr);
	const auto before = iga::CapturePetscOptions(source.value);
	PetscInt unused_before = 0, unused_after = 0; Check(PetscOptionsAllUsed(source.value, &unused_before));
	iga::PetscSolverOptions left(comm, "left_", source.value, {{"-ksp_type", std::string("richardson")}});
	iga::PetscSolverOptions right(comm, "right_", source.value);
	Check(PetscOptionsAllUsed(source.value, &unused_after)); Require(unused_before==unused_after, "snapshot consumed source options");
	// Values are copied exactly, not reparsed from an ambiguous display string.
	const char* text = nullptr; PetscBool present = PETSC_FALSE;
	Check(PetscOptionsFindPair(left.Database(), "left_", "-unused_text", &text, &present));
	Require(present && std::string(text)=="spaces \"quotes\"\n-ksp_type cg ", "snapshot lost exact option bytes");
	// Snapshot isolation: changes to the source do not retroactively reconfigure owners.
	source.Set("-ksp_type", "richardson");
	Solve(comm, right, "gmres", "jacobi"); Solve(comm, left, "cg", "lu");
	source.Set("-ksp_type", "gmres");
	Require(before==iga::CapturePetscOptions(source.value), "snapshot changed source values");
	PetscBool used = PETSC_FALSE; Check(PetscOptionsUsed(source.value, "left_ksp_type", &used)); Require(used, "prefixed option usage was not propagated");
	Check(PetscOptionsUsed(source.value, "sub_pc_type", &used)); Require(used, "nested inherited usage was not propagated");
	Check(PetscOptionsUsed(source.value, "unused_text", &used)); Require(used, "explicit snapshot lookup usage was not propagated");
	int rejections = 0;
	for (const auto& prefix : {std::string("-bad_"), std::string("no_tail"), std::string(129,'a')+"_"}) {
		try { iga::PetscSolverOptions bad(comm, rank==size-1 ? prefix : "valid_", source.value); }
		catch (const std::exception&) { ++rejections; }
	}
	if (size>1) {
		try { iga::PetscSolverOptions bad(comm, rank==size-1 ? "different_" : "same_", source.value); }
		catch (const std::exception&) { ++rejections; }
		if (rank==size-1) source.Set("-sub_pc_type", "jacobi");
		try { iga::PetscSolverOptions bad(comm, "different_input_", source.value); }
		catch (const std::exception&) { ++rejections; }
		source.Set("-sub_pc_type", "lu");
	}
	Require(rejections==(size>1 ? 5 : 3), "invalid options were not rejected collectively");
	iga::PetscSolverOptions retry(comm, "retry_", source.value); Solve(comm, retry, "gmres", "lu");
	const auto global_before = iga::CapturePetscOptions(nullptr);
	for (int mode = 0; mode < 3; ++mode) {
		int rejected = 0;
		try {
			retry.Call("injected options scope failure", [&]() -> PetscErrorCode {
				if (mode==0) return rank==size-1 ? PETSC_ERR_USER : 0;
				if (mode==1 && rank==size-1) throw std::runtime_error("local callback failure");
				if (mode==2) retry.Call("nested same-owner call", [] { return PetscErrorCode{0}; });
				return 0;
			});
		} catch (const std::exception&) { rejected = 1; }
		Require(rejected==1, "scope failure did not reach every rank");
		Require(iga::CapturePetscOptions(nullptr)==global_before, "failure left a private option database active");
	}
	// Exercise an actual PETSc error with the default caller handler installed.
	// PETSc's default traceback handler can abort non-root ranks here.
	{
		Linear invalid;
		Check(KSPCreate(comm, &invalid.solver)); retry.Attach(invalid.solver);
		int rejected = 0;
		try { retry.Call("invalid solver type", [&] { return KSPSetType(invalid.solver, "unavailable_solver"); }); }
		catch (const std::exception&) { rejected = 1; }
		Require(rejected==1, "actual PETSc error did not reach every rank");
		Require(iga::CapturePetscOptions(nullptr)==global_before, "PETSc error left private options active");
	}
	int caller_errors = 0;
	Check(PetscPushErrorHandler(CountErrors, &caller_errors));
	for (int mode = 0; mode < 3; ++mode) {
		try {
			retry.Call("caller handler isolation", [&]() -> PetscErrorCode {
				if (mode==1) return PetscError(comm, __LINE__, "Run", __FILE__, PETSC_ERR_USER, PETSC_ERROR_INITIAL, "scope error");
				if (mode==2) throw std::runtime_error("scope exception");
				return 0;
			});
		} catch (const std::exception&) { Require(mode!=0, "healthy scope failed"); }
		Require(caller_errors==mode, "private scope invoked the caller handler");
		const auto code = PetscError(comm, __LINE__, "Run", __FILE__, PETSC_ERR_USER, PETSC_ERROR_INITIAL, "caller handler probe");
		Require(code==PETSC_ERR_USER && caller_errors==mode+1, "caller error handler was not restored");
	}
	Check(PetscPopErrorHandler());
	retry.Call("healthy scope retry", [] { return PetscErrorCode{0}; });
	Require(iga::CapturePetscOptions(nullptr)==global_before, "retry changed the global option database");
	source.Set("-split_pc_type", "fieldsplit"); source.Set("-fieldsplit_first_ksp_type", "preonly");
	source.Set("-fieldsplit_first_pc_type", "jacobi"); source.Set("-fieldsplit_second_ksp_type", "preonly");
	source.Set("-split_fieldsplit_second_pc_type", "lu"); source.Set("-split_fieldsplit_second_pc_factor_mat_solver_type", "mumps");
	iga::PetscSolverOptions fields(comm, "split_", source.value); Solve(comm, fields, "gmres", "jacobi", true);
	Solve(comm, fields, "gmres", "jacobi", true, true);
}
}
int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int rank = 0; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	try {
		Run(PETSC_COMM_WORLD); RunFamily(PETSC_COMM_WORLD);
		MPI_Comm group = MPI_COMM_NULL; MPI_Comm_split(PETSC_COMM_WORLD, rank==0 ? 0 : 1, rank, &group);
		Run(group); RunFamily(group); MPI_Comm_free(&group);
	} catch (const std::exception& error) { std::cerr << "rank " << rank << ": " << error.what() << '\n'; MPI_Abort(PETSC_COMM_WORLD, 2); }
	PetscFinalize();
}
