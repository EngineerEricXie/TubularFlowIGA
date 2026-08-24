template <int Fields>
int ConfiguredTransport(iga::Database& database, const fs::path& case_dir,
	const iga::SimulationConfiguration& configuration, const iga::CompiledLinearSystem& system,
	const fs::path& velocity_path, const fs::path* output_path)
{
	static_assert(Fields >= 1 && Fields <= 8, "supported CUDA field block width");
	if (system.fields.size() != Fields) throw std::runtime_error("internal CUDA field dispatch mismatch");
	const auto total_start = Clock::now();
	const auto labels = iga::ReadPointLabels((case_dir/"controlmesh.vtk").string(), database.header().nodes);
	const auto velocity_host = iga::ReadVelocity(velocity_path.string(), database.header().nodes);
	const auto boundaries = iga::ResolveScalarBoundaries(configuration, system, labels);
	const auto initial = iga::InitialScalarValues(configuration, system);

	std::vector<DeviceWeakTerm> host_terms;
	host_terms.reserve(system.terms.size());
	for (const auto& term : system.terms) {
		if (term.kind == iga::TermKind::Advection && term.velocity != "prescribed")
			throw std::runtime_error("CUDA advection currently requires velocity source 'prescribed'");
		host_terms.push_back({static_cast<int>(term.kind), static_cast<int>(term.equation),
			static_cast<int>(term.trial), term.coefficient});
	}
	std::vector<int> host_supg(Fields, 0);
	for (const auto& stabilization : system.stabilization) {
		if (stabilization.method != "supg" || stabilization.velocity != "prescribed")
			throw std::runtime_error("CUDA stabilization currently supports SUPG with prescribed velocity");
		host_supg[system.field_index.at(stabilization.equation)] = 1;
	}

	FlatMesh host(database);
	BlockPattern pattern_host(host);
	PrintMesh(host, &pattern_host);
	const auto preprocess_end = Clock::now();
	DeviceMesh mesh(host);
	DevicePattern pattern(pattern_host);
	ReferenceData reference;
	GeometryData geometry(mesh, reference, false);
	RequireGeometry(geometry);
	ElementTiles tiles(host);
	DeviceBuffer<DeviceWeakTerm> terms(host_terms.size());
	terms.CopyFromHost(host_terms.data(), host_terms.size());
	DeviceBuffer<int> supg(host_supg.size());
	supg.CopyFromHost(host_supg.data(), host_supg.size());
	DeviceBuffer<int> boundary_mask(boundaries.constrained.size());
	boundary_mask.CopyFromHost(boundaries.constrained.data(), boundaries.constrained.size());
	DeviceBuffer<double> boundary_value(boundaries.value.size());
	boundary_value.CopyFromHost(boundaries.value.data(), boundaries.value.size());
	auto velocity = Flatten(velocity_host);

	constexpr int linear_restart = 50;
	const std::size_t matrix_bytes = static_cast<std::size_t>(pattern.view().blocks)*Fields*Fields*sizeof(double);
	const std::size_t vector_bytes = static_cast<std::size_t>(host.nodes)*Fields*sizeof(double);
	CheckAvailableMemory(2*matrix_bytes+(12+linear_restart)*vector_bytes);
	BlockMatrix<Fields> left(pattern), previous(pattern);
	DeviceBuffer<double> current(host.nodes*Fields), next(host.nodes*Fields), rhs(host.nodes*Fields), forcing(host.nodes*Fields);
	left.Clear(); previous.Clear(); next.Clear(); rhs.Clear(); forcing.Clear();
	std::vector<double> host_current(host.nodes*Fields, 0.0);
	for (std::size_t node = 0; node < host.nodes; ++node)
		for (int field = 0; field < Fields; ++field) {
			const auto dof = node*Fields+field;
			host_current[dof] = boundaries.constrained[dof] ? boundaries.value[dof] : initial[field];
		}
	current.CopyFromHost(host_current.data(), host_current.size());

	const auto assembly_start = Clock::now();
	AssembleGenericTransportKernel<Fields><<<tiles.view().count,kPairTile>>>(
		mesh.view(), reference.view(), geometry.view(), tiles.view(), pattern.view(), velocity.data(),
		terms.data(), static_cast<int>(host_terms.size()), supg.data(), system.dt,
		left.values(), previous.values(), forcing.data());
	CheckKernel("AssembleGenericTransportKernel");
	const int node_blocks = (static_cast<int>(host.nodes)+255)/256;
	ApplyGenericBoundaryKernel<Fields><<<node_blocks,256>>>(
		pattern.view(), boundary_mask.data(), left.values(), previous.values());
	CheckKernel("ApplyGenericBoundaryKernel");
	Check(cudaDeviceSynchronize(), "configured transport assembly synchronize");
	const auto assembly_end = Clock::now();

	GmresWorkspace<Fields> workspace(pattern.view().nodes, linear_restart);
	long long total_iterations = 0;
	double peak_gpu_used = 0.0;
	unsigned int singular_blocks = 0;
	const auto solve_start = Clock::now();
	const int vector_size = static_cast<int>(host.nodes*Fields);
	const int vector_blocks = (vector_size+255)/256;
	for (int step = 0; step < system.steps; ++step) {
		BlockSpmvKernel<Fields><<<node_blocks,256>>>(pattern.view(), previous.values(), current.data(), rhs.data());
		CheckKernel("configured transport previous MatMult");
		AddVectorKernel<<<vector_blocks,256>>>(vector_size, forcing.data(), rhs.data());
		CheckKernel("configured transport volume source");
		SetGenericBoundaryVectorKernel<Fields><<<node_blocks,256>>>(
			static_cast<int>(host.nodes), boundary_mask.data(), boundary_value.data(), rhs.data());
		CheckKernel("configured transport boundary rhs");
		if (step == 0) next.Clear();
		else Check(cudaMemcpy(next.data(), current.data(), current.bytes(), cudaMemcpyDeviceToDevice),
			"configured transport warm start copy");
		const auto result = SolveGmres(left, workspace, rhs.data(), next.data(), 10000, 1e-8, false, 0, step == 0);
		if (!result.converged)
			throw std::runtime_error("configured transport GMRES failed at step "+std::to_string(step));
		total_iterations += result.iterations;
		peak_gpu_used = std::max(peak_gpu_used, result.device_used_gib);
		singular_blocks = std::max(singular_blocks, result.singular_diagonal_blocks);
		std::swap(current, next);
	}
	Check(cudaDeviceSynchronize(), "configured transport solve synchronize");
	const auto solve_end = Clock::now();
	BlasHandle blas;
	double norm = 0.0;
	Check(cublasDnrm2(blas, vector_size, current.data(), 1, &norm), "configured transport final norm");
	if (output_path) {
		std::vector<double> output(host.nodes*Fields);
		current.CopyToHost(output.data(), output.size());
		std::ofstream values(*output_path), metadata(output_path->string()+".fields");
		if (!values || !metadata) throw std::runtime_error("cannot create configured CUDA output");
		values << std::setprecision(17);
		for (const auto& name : system.fields) metadata << name << '\n';
		for (std::size_t node = 0; node < host.nodes; ++node) {
			values << node;
			for (int field = 0; field < Fields; ++field) values << ' ' << output[node*Fields+field];
			values << '\n';
		}
	}
	std::cout << "iga_cuda_solve system=" << system.name << " fields=" << Fields
		<< " steps=" << system.steps
		<< " preprocess_s=" << std::chrono::duration<double>(preprocess_end-total_start).count()
		<< " assembly_s=" << std::chrono::duration<double>(assembly_end-assembly_start).count()
		<< " solve_s=" << std::chrono::duration<double>(solve_end-solve_start).count()
		<< " total_iterations=" << total_iterations << " final_l2=" << norm
		<< " singular_diagonal_blocks=" << singular_blocks << " gpu_used_gib=" << peak_gpu_used << '\n';
	return 0;
}

int SolveConfigured(int argc, char** argv)
{
	if (argc < 4) throw std::runtime_error(
		"usage: iga_cuda solve DATABASE.ntiga CASE_DIR [SYSTEM] [OUTPUT] [VELOCITY]");
	iga::Database database(argv[2]);
	const fs::path case_dir(argv[3]);
	const auto configuration = iga::ReadSimulationConfiguration((case_dir/"simulation_config.json").string());
	const auto system_name = argc >= 5 ? std::string(argv[4]) : iga::FirstLinearTransportSystem(configuration);
	const auto system = iga::CompileLinearSystem(configuration, system_name);
	const fs::path velocity_path = argc >= 7 ? fs::path(argv[6]) : case_dir/"initial_velocityfield.txt";
	const fs::path output_storage = argc >= 6 ? fs::path(argv[5]) : fs::path();
	const fs::path* output = argc >= 6 ? &output_storage : nullptr;
	switch (system.fields.size()) {
	case 1: return ConfiguredTransport<1>(database, case_dir, configuration, system, velocity_path, output);
	case 2: return ConfiguredTransport<2>(database, case_dir, configuration, system, velocity_path, output);
	case 3: return ConfiguredTransport<3>(database, case_dir, configuration, system, velocity_path, output);
	case 4: return ConfiguredTransport<4>(database, case_dir, configuration, system, velocity_path, output);
	case 5: return ConfiguredTransport<5>(database, case_dir, configuration, system, velocity_path, output);
	case 6: return ConfiguredTransport<6>(database, case_dir, configuration, system, velocity_path, output);
	case 7: return ConfiguredTransport<7>(database, case_dir, configuration, system, velocity_path, output);
	case 8: return ConfiguredTransport<8>(database, case_dir, configuration, system, velocity_path, output);
	default: throw std::runtime_error("CUDA supports 1 through 8 coupled scalar fields per system");
	}
}
