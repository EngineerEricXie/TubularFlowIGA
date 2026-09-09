#define main implicit_failure_reference_main
#include "test_one_d_implicit_failure.cpp"
#undef main

namespace {
void RunOptions(MPI_Comm comm, const char* swc)
{
	int rank = 0, ranks = 0; MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &ranks);
	iga::OneDFlowSystemDefinition flow;
	flow.name = "flow"; flow.model = iga::OneDFlowModel::Compliant;
	flow.scheme = iga::OneDFlowScheme::ImplicitPetsc; flow.dynamic_viscosity = 0.004; flow.density = 1060;
	const auto network = iga::ReadOneDNetwork(swc, 1.0, 3, flow.dynamic_viscosity);
	iga::OneDFlowState initial;
	for (const int node : network.outlet_nodes) {
		iga::OneDOutletState outlet; outlet.node = node; outlet.kind = iga::OneDOutletKind::Pressure;
		initial.outlets.push_back(outlet);
	}
	iga::InitializeCompliantOneDFromRigid(network, flow, initial, 1e-9, 1e-3);
	Option("-ksp_type", "preonly"); Option("-pc_type", "lu"); Option("-pc_factor_mat_solver_type", "mumps");
	Option("-left_ksp_type", "fgmres"); Option("-right_ksp_type", "gmres");
	Option("-left_ksp_rtol", "1e-12"); Option("-right_ksp_rtol", "1e-12");
	Option("-left_snes_type", "newtonls"); Option("-right_snes_type", "newtontr");
	iga::OneDPetscSolverContext left(comm, "left_"), right(comm, "right_");
	const auto original = iga::CapturePetscOptions(nullptr);
	for (const auto formulation : {iga::OneDImplicitFormulation::PressureNetwork,
		iga::OneDImplicitFormulation::LinearizedAQ, iga::OneDImplicitFormulation::NonlinearAQ,
		iga::OneDImplicitFormulation::ImplicitPde}) {
		flow.formulation = formulation;
		auto reference = initial, first = initial, second = initial;
		iga::AdvanceImplicitOneD(network, flow, reference, 2e-9, 1e-3, comm);
		iga::AdvanceImplicitOneD(network, flow, first, 2e-9, 1e-3, comm, &left);
		iga::AdvanceImplicitOneD(network, flow, second, 2e-9, 1e-3, comm, &right);
		Compare(reference, first); Compare(reference, second);
		Require(left.linear.prefix=="left_" && right.linear.prefix=="right_", "wrong prefix");
		Require(left.linear.ksp=="fgmres" && right.linear.ksp=="gmres", "wrong independent KSP");
		Require(left.linear.factor_backend=="mumps" && right.linear.factor_backend=="mumps", "wrong backend");
		Require(left.linear.last_reason>0 && right.linear.last_reason>0, "linear convergence missing");
		if (formulation==iga::OneDImplicitFormulation::NonlinearAQ || formulation==iga::OneDImplicitFormulation::ImplicitPde) {
			Require(left.nonlinear_type=="newtonls" && right.nonlinear_type=="newtontr", "wrong independent SNES");
			Require(left.nonlinear_reason>0 && right.nonlinear_reason>0, "nonlinear convergence missing");
		}
		for (const auto& invalid : {std::string("ksp_type"), std::string("snes_type"), std::string("pc_factor_mat_solver_type")}) {
			if (invalid=="snes_type" && !left.nonlinear) continue;
			if (invalid=="pc_factor_mat_solver_type" && ranks==1) continue;
			const auto key = "-bad_"+invalid;
			Option(key.c_str(), invalid=="pc_factor_mat_solver_type" ? "petsc" : "not_an_iga_solver");
			iga::OneDPetscSolverContext bad(comm, "bad_");
			Option(key.c_str(), nullptr);
			auto current = initial;
			const auto before = iga::PackOneDCheckpointState(current, {}, network);
			bool rejected = false;
			try { iga::AdvanceImplicitOneD(network, flow, current, 2e-9, 1e-3, comm, &bad); }
			catch (const std::exception&) { rejected = true; }
			Require(rejected, "invalid scoped solver was accepted");
			Require(before==iga::PackOneDCheckpointState(current, {}, network), "failed solve published fields");
			iga::AdvanceImplicitOneD(network, flow, current, 2e-9, 1e-3, comm, &left); Compare(first, current);
		}
		Require(original==iga::CapturePetscOptions(nullptr), "source options changed");
		if (rank==0) { std::cout << "one_d_options ranks=" << ranks << " formulation=" << static_cast<int>(formulation) << " fields=passed retry=passed\n";
			iga::WriteOneDPetscSolverConfiguration(std::cout, left, 1); iga::WriteOneDPetscSolverConfiguration(std::cout, right, 1); }
	}
}
}
int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int rank = 0; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	try {
		Require(argc==2, "provide fixture SWC"); RunOptions(PETSC_COMM_WORLD, argv[1]);
		MPI_Comm group = MPI_COMM_NULL; MPI_Comm_split(PETSC_COMM_WORLD, rank==0 ? 0 : 1, rank, &group);
		RunOptions(group, argv[1]); MPI_Comm_free(&group);
	} catch (const std::exception& error) { std::cerr << "rank " << rank << ": " << error.what() << '\n'; MPI_Abort(PETSC_COMM_WORLD, 2); }
	PetscFinalize();
}
