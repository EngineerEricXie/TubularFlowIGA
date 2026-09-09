#define IGA_MULTIDOMAIN_NO_MAIN
#define IGA_COUPLED_CHECKPOINT_TESTING
#define IGA_ASSET_INPUT_TESTING
#include "../src/iga_1d_3d_bifurcation.cpp"

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int rank = 0; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	const char* requested = std::getenv("IGA_NATIVE_CHECKPOINT_FAULT");
	const std::string fault = requested ? requested : "";
	iga::CurrentPhaseProfile().EnableFromEnvironment();
	bool input_signal = false;
	iga::AssetReadProbeForTesting() = [&] {
		if (fault == "signal-input" && rank == 0 && !input_signal) { input_signal = true; ::raise(SIGUSR1); }
	};
	iga::coupled_checkpoint_detail::before_shard_payload = [&](const auto& epoch, const auto& spec, int fd) {
		if (epoch.accepted_steps != 4) return;
		if (fault == "stream" && rank == 0 && spec.id.find(".flow.rank-") != std::string::npos) {
			const char partial[5] = {1, 2, 3, 4, 5};
			if (::write(fd, partial, sizeof(partial)) != static_cast<ssize_t>(sizeof(partial))) ::_exit(88);
			std::cerr << "native_checkpoint_crash stage=stream bytes=5 step=4" << std::endl; ::_exit(86);
		}
		if (fault == "local-write" && rank == 1) throw std::runtime_error("injected rank-1 native checkpoint write failure");
	};
	iga::coupled_checkpoint_detail::before_manifest_publication = [&](const auto& epoch) {
		if (fault == "signal" && epoch.accepted_steps == 3 && rank == 0) ::raise(SIGUSR1);
		if (fault == "manifest" && epoch.accepted_steps == 4 && rank == 0) {
			std::cerr << "native_checkpoint_crash stage=before-manifest step=4" << std::endl; ::_exit(87);
		}
	};
	int status = iga::RunMultidomainFlow(argc, argv, PETSC_COMM_WORLD);
	int ranks = 1; MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
	try { iga::CollectiveLocalStage(PETSC_COMM_WORLD, "test graph profile", [&] { iga::CurrentPhaseProfile().Write(std::cout, rank, ranks, status); iga::FlushCheckedText(std::cout); }); }
	catch (const std::exception& error) { if (rank == 0) std::cerr << error.what() << '\n'; status = 1; }
	PetscFinalize(); return status;
}
