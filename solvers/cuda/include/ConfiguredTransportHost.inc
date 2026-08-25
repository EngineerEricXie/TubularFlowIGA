struct CudaConfiguredTransportOptions {
	std::string system;
	fs::path output;
	fs::path velocity;
	fs::path checkpoint;
	fs::path restart;
	int output_every = 0;
	int checkpoint_every = 0;
	int stop_after_step = 0;
};

CudaConfiguredTransportOptions ParseCudaConfiguredTransportOptions(int argc, char** argv)
{
	CudaConfiguredTransportOptions options;
	int positional = 0;
	for (int i = 4; i < argc; ++i) {
		const std::string argument(argv[i]);
		if (argument.rfind("--", 0) != 0) {
			if (positional == 0) options.system = argument;
			else if (positional == 1) options.output = argument;
			else if (positional == 2) options.velocity = argument;
			else throw std::runtime_error("too many positional arguments");
			++positional;
			continue;
		}
		if (i+1 >= argc) throw std::runtime_error(argument+" requires a value");
		const std::string value(argv[++i]);
		if (argument == "--system") options.system = value;
		else if (argument == "--output") options.output = value;
		else if (argument == "--velocity") options.velocity = value;
		else if (argument == "--output-every")
			options.output_every = ParsePositiveInteger(value, argument);
		else if (argument == "--checkpoint") options.checkpoint = value;
		else if (argument == "--checkpoint-every")
			options.checkpoint_every = ParsePositiveInteger(value, argument);
		else if (argument == "--restart") options.restart = value;
		else if (argument == "--stop-after-step")
			options.stop_after_step = ParsePositiveInteger(value, argument);
		else throw std::runtime_error("unknown option: "+argument);
	}
	if (options.output_every > 0 && options.output.empty())
		throw std::runtime_error("--output-every requires --output");
	if (options.checkpoint_every > 0 && options.checkpoint.empty())
		throw std::runtime_error("--checkpoint-every requires --checkpoint");
	return options;
}

std::vector<double> ReadCudaTransportCheckpoint(const fs::path& prefix,
	const iga::TransportCheckpointMetadata& metadata)
{
	if (metadata.state_format != "raw_float64")
		throw std::runtime_error("CUDA transport restart requires raw_float64 state");
	fs::path path(metadata.state_file);
	if (path.is_relative())
		path = iga::TransportCheckpointMetadataPath(prefix).parent_path()/path;
	std::ifstream input(path, std::ios::binary);
	if (!input)
		throw std::runtime_error("cannot open CUDA transport checkpoint state: "+path.string());
	std::vector<double> values(
		static_cast<std::size_t>(metadata.nodes)*metadata.fields.size());
	input.read(reinterpret_cast<char*>(values.data()),
		static_cast<std::streamsize>(values.size()*sizeof(double)));
	char extra = 0;
	if (!input || input.read(&extra, 1))
		throw std::runtime_error(
			"CUDA transport checkpoint state is truncated or has trailing data");
	return values;
}

void WriteCudaTransportCheckpoint(const fs::path& prefix,
	iga::TransportCheckpointMetadata metadata, const DeviceBuffer<double>& state,
	std::size_t values)
{
	std::vector<double> host(values);
	state.CopyToHost(host.data(), host.size());
	const auto state_path = iga::TransportCheckpointStatePath(prefix);
	std::ofstream output(state_path, std::ios::binary);
	if (!output)
		throw std::runtime_error("cannot create CUDA transport checkpoint state: "+state_path.string());
	output.write(reinterpret_cast<const char*>(host.data()),
		static_cast<std::streamsize>(host.size()*sizeof(double)));
	if (!output)
		throw std::runtime_error("cannot write CUDA transport checkpoint state: "+state_path.string());
	output.close();
	if (!output)
		throw std::runtime_error("cannot close CUDA transport checkpoint state: "+state_path.string());
	metadata.state_file = state_path.filename().string();
	metadata.state_format = "raw_float64";
	iga::WriteTransportCheckpointMetadata(prefix, metadata);
}

template <int Fields>
void WriteCudaConfiguredOutput(const fs::path& path,
	const std::vector<std::string>& fields, std::size_t nodes,
	const DeviceBuffer<double>& state)
{
	std::vector<double> output(nodes*Fields);
	state.CopyToHost(output.data(), output.size());
	std::ofstream values(path), metadata(path.string()+".fields");
	if (!values || !metadata)
		throw std::runtime_error("cannot create configured CUDA output");
	values << std::setprecision(17);
	for (const auto& name : fields) metadata << name << '\n';
	for (std::size_t node = 0; node < nodes; ++node) {
		values << node;
		for (int field = 0; field < Fields; ++field)
			values << ' ' << output[node*Fields+field];
		values << '\n';
	}
	if (!values || !metadata)
		throw std::runtime_error("cannot write configured CUDA output");
}

template <int Fields>
int ConfiguredTransport(iga::Database& database, const fs::path& case_dir,
	const iga::SimulationConfiguration& configuration, const iga::CompiledLinearSystem& system,
	const CudaConfiguredTransportOptions& options)
{
	static_assert(Fields >= 1 && Fields <= 8, "supported CUDA field block width");
	if (system.fields.size() != Fields) throw std::runtime_error("internal CUDA field dispatch mismatch");
	const auto total_start = Clock::now();
	if (options.stop_after_step > system.steps)
		throw std::runtime_error("--stop-after-step exceeds configured transport steps");
	const auto run_end_step = options.stop_after_step > 0
		? options.stop_after_step : system.steps;
	const auto labels = iga::ReadPointLabels((case_dir/"controlmesh.vtk").string(), database.header().nodes);
	std::vector<std::array<double, 3>> prescribed_velocity;
	std::vector<iga::VelocitySnapshot> velocity_snapshots;
	const iga::VelocitySourceDefinition* velocity_source = nullptr;
	if (system.velocity_source == "prescribed")
		prescribed_velocity = iga::ReadVelocity((options.velocity.empty()
			? case_dir/"initial_velocityfield.txt" : options.velocity).string(),
			database.header().nodes);
	else {
		velocity_source = &iga::FindVelocitySource(configuration, system.velocity_source);
		velocity_snapshots = iga::ReadVelocityManifest(case_dir/velocity_source->manifest);
	}
	const auto initial_configuration = iga::MaterializeBoundaryWaveforms(configuration, case_dir, 0.0);
	const auto boundaries = iga::ResolveScalarBoundaries(initial_configuration, system, labels);
	const auto initial = iga::InitialScalarValues(configuration, system);

	std::vector<DeviceWeakTerm> host_terms;
	host_terms.reserve(system.terms.size());
	for (const auto& term : system.terms) {
		if (term.kind == iga::TermKind::Advection && term.velocity != system.velocity_source)
			throw std::runtime_error("CUDA transport received inconsistent velocity sources");
		host_terms.push_back({static_cast<int>(term.kind), static_cast<int>(term.equation),
			static_cast<int>(term.trial), term.coefficient});
	}
	std::vector<int> host_supg(Fields, 0);
	for (const auto& stabilization : system.stabilization) {
		if (stabilization.method != "supg" || stabilization.velocity != system.velocity_source)
			throw std::runtime_error("CUDA transport requires SUPG and advection to use the same velocity source");
		host_supg[system.field_index.at(stabilization.equation)] = 1;
	}

	FlatMesh host(database);
	std::vector<int> required_surface_labels;
	for (const auto& boundary : configuration.boundaries)
		for (const auto& condition : boundary.conditions)
			if (system.field_index.count(condition.field)
				&& (condition.kind == iga::FieldBoundaryKind::Flux
					|| condition.kind == iga::FieldBoundaryKind::Robin)
				&& std::find(required_surface_labels.begin(), required_surface_labels.end(),
					boundary.label) == required_surface_labels.end())
				required_surface_labels.push_back(boundary.label);
	std::vector<int> surface_element, surface_face, surface_first_pair;
	std::vector<double> surface_flux, surface_robin, surface_robin_source;
	std::vector<int> found_surface_labels;
	for (std::size_t element = 0; element < host.elements(); ++element)
		for (int face = 0; face < 6; ++face) {
			const int label = host.boundary_labels[element*6+static_cast<std::size_t>(face)];
			if (label < 0) continue;
			const auto boundary = std::find_if(configuration.boundaries.begin(), configuration.boundaries.end(),
				[label](const iga::NamedBoundaryDefinition& item) { return item.label == label; });
			if (boundary == configuration.boundaries.end())
				throw std::runtime_error("element boundary face has no simulation_config.json definition");
			std::array<double,Fields> flux{}, robin{}, robin_source{};
			bool active_surface = false;
			for (const auto& condition : boundary->conditions) {
				const auto field = system.field_index.find(condition.field);
				if (field == system.field_index.end()) continue;
				if (condition.kind == iga::FieldBoundaryKind::Flux) {
					flux[field->second] = condition.value[0];
					active_surface = true;
				} else if (condition.kind == iga::FieldBoundaryKind::Robin) {
					robin[field->second] = condition.coefficient;
					robin_source[field->second] = condition.coefficient*condition.exterior_value;
					active_surface = true;
				}
			}
			if (!active_surface) continue;
			if (std::find(found_surface_labels.begin(), found_surface_labels.end(), label)
				== found_surface_labels.end()) found_surface_labels.push_back(label);
			const int begin = host.element_offsets[element];
			const int nen = host.element_offsets[element+1]-begin;
			for (int pair = 0; pair < nen*nen; pair += kPairTile) {
				surface_element.push_back(static_cast<int>(element));
				surface_face.push_back(face);
				surface_first_pair.push_back(pair);
				surface_flux.insert(surface_flux.end(), flux.begin(), flux.end());
				surface_robin.insert(surface_robin.end(), robin.begin(), robin.end());
				surface_robin_source.insert(surface_robin_source.end(), robin_source.begin(), robin_source.end());
			}
		}
	for (const int label : required_surface_labels)
		if (std::find(found_surface_labels.begin(), found_surface_labels.end(), label)
			== found_surface_labels.end())
			throw std::runtime_error("configured scalar surface boundary label "
				+std::to_string(label)+" has no boundary faces in the .ntiga database; repack with iga_pack");
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
	DeviceBuffer<int> surface_elements(surface_element.size());
	DeviceBuffer<int> surface_faces(surface_face.size());
	DeviceBuffer<int> surface_pairs(surface_first_pair.size());
	DeviceBuffer<double> surface_flux_device(surface_flux.size());
	DeviceBuffer<double> surface_robin_device(surface_robin.size());
	DeviceBuffer<double> surface_robin_source_device(surface_robin_source.size());
	surface_elements.CopyFromHost(surface_element.data(), surface_element.size());
	surface_faces.CopyFromHost(surface_face.data(), surface_face.size());
	surface_pairs.CopyFromHost(surface_first_pair.data(), surface_first_pair.size());
	surface_flux_device.CopyFromHost(surface_flux.data(), surface_flux.size());
	surface_robin_device.CopyFromHost(surface_robin.data(), surface_robin.size());
	surface_robin_source_device.CopyFromHost(surface_robin_source.data(), surface_robin_source.size());
	const DeviceSurfaceView surfaces{static_cast<int>(surface_element.size()),
		surface_elements.data(), surface_faces.data(), surface_pairs.data(),
		surface_flux_device.data(), surface_robin_device.data(), surface_robin_source_device.data()};
	DeviceBuffer<double> velocity(host.nodes*3);

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
	int start_step = 0;
	if (!options.restart.empty()) {
		const auto metadata = iga::ReadTransportCheckpointMetadata(options.restart);
		iga::ValidateTransportCheckpoint(metadata, database.header().nodes,
			system.fields, system.name, system.velocity_source, system.steps, system.dt);
		const auto restart_state = ReadCudaTransportCheckpoint(options.restart, metadata);
		current.CopyFromHost(restart_state.data(), restart_state.size());
		start_step = metadata.completed_step;
		std::cout << "restart=" << options.restart.string()
			<< " completed_step=" << start_step
			<< " physical_time=" << metadata.physical_time << '\n';
	}

	const int node_blocks = (static_cast<int>(host.nodes)+255)/256;
	double assembly_seconds = 0.0;
	auto assemble_operators = [&](const std::vector<std::array<double, 3>>& host_velocity) {
		const auto start = Clock::now();
		velocity.CopyFromHost(reinterpret_cast<const double*>(host_velocity.data()), host.nodes*3);
		left.Clear(); previous.Clear(); forcing.Clear();
		AssembleGenericTransportKernel<Fields><<<tiles.view().count,kPairTile>>>(
			mesh.view(), reference.view(), geometry.view(), tiles.view(), pattern.view(), velocity.data(),
			terms.data(), static_cast<int>(host_terms.size()), supg.data(), system.dt,
			left.values(), previous.values(), forcing.data());
		CheckKernel("AssembleGenericTransportKernel");
		if (surfaces.count) {
			AssembleGenericSurfaceKernel<Fields><<<surfaces.count,kPairTile>>>(
				mesh.view(), reference.view(), surfaces, pattern.view(), system.dt,
				left.values(), forcing.data());
			CheckKernel("AssembleGenericSurfaceKernel");
		}
		ApplyGenericBoundaryKernel<Fields><<<node_blocks,256>>>(
			pattern.view(), boundary_mask.data(), left.values(), previous.values());
		CheckKernel("ApplyGenericBoundaryKernel");
		Check(cudaDeviceSynchronize(), "configured transport assembly synchronize");
		assembly_seconds += std::chrono::duration<double>(Clock::now()-start).count();
	};
	auto velocity_at = [&](double time) {
		const auto interpolation = iga::ResolveVelocityInterpolation(
			velocity_snapshots, time, velocity_source->out_of_range);
		auto lower = iga::ReadVelocity(
			(case_dir/velocity_snapshots[interpolation.lower].file).string(), database.header().nodes);
		if (interpolation.lower == interpolation.upper) return lower;
		const auto upper = iga::ReadVelocity(
			(case_dir/velocity_snapshots[interpolation.upper].file).string(), database.header().nodes);
		for (std::size_t node = 0; node < lower.size(); ++node)
			for (int component = 0; component < 3; ++component)
				lower[node][component] = (1.0-interpolation.upper_weight)*lower[node][component]
					+interpolation.upper_weight*upper[node][component];
		return lower;
	};
	if (!velocity_source) assemble_operators(prescribed_velocity);

	GmresWorkspace<Fields> workspace(pattern.view().nodes, linear_restart);
	long long total_iterations = 0;
	double peak_gpu_used = 0.0;
	unsigned int singular_blocks = 0;
	const auto solve_start = Clock::now();
	double linear_seconds = 0.0;
	const int vector_size = static_cast<int>(host.nodes*Fields);
	const int vector_blocks = (vector_size+255)/256;
	for (int step = start_step; step < run_end_step; ++step) {
		if (velocity_source) assemble_operators(velocity_at((step+1)*system.dt));
		const auto step_configuration = iga::MaterializeBoundaryWaveforms(
			configuration, case_dir, (step+1)*system.dt);
		const auto step_boundaries = iga::ResolveScalarBoundaries(step_configuration, system, labels);
		boundary_value.CopyFromHost(step_boundaries.value.data(), step_boundaries.value.size());
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
		const auto linear_start = Clock::now();
		const auto result = SolveGmres(left, workspace, rhs.data(), next.data(),
			10000, 1e-8, false, 0, velocity_source || step == start_step);
		linear_seconds += std::chrono::duration<double>(Clock::now()-linear_start).count();
		if (!result.converged)
			throw std::runtime_error("configured transport GMRES failed at step "+std::to_string(step));
		total_iterations += result.iterations;
		peak_gpu_used = std::max(peak_gpu_used, result.device_used_gib);
		singular_blocks = std::max(singular_blocks, result.singular_diagonal_blocks);
		std::swap(current, next);
		const auto completed_step = step+1;
		if (options.output_every > 0
			&& completed_step%options.output_every == 0)
			WriteCudaConfiguredOutput<Fields>(
				iga::TimeIndexedPath(options.output, completed_step),
				system.fields, host.nodes, current);
		if (!options.checkpoint.empty()
			&& (completed_step == run_end_step
				|| (options.checkpoint_every > 0
					&& completed_step%options.checkpoint_every == 0))) {
			iga::TransportCheckpointMetadata metadata;
			metadata.nodes = database.header().nodes;
			metadata.fields = system.fields;
			metadata.system = system.name;
			metadata.velocity_source = system.velocity_source;
			metadata.completed_step = completed_step;
			metadata.physical_time = completed_step*system.dt;
			metadata.dt = system.dt;
			WriteCudaTransportCheckpoint(
				options.checkpoint, metadata, current, host.nodes*Fields);
			std::cout << "checkpoint=" << options.checkpoint.string()
				<< " completed_step=" << completed_step << '\n';
		}
	}
	Check(cudaDeviceSynchronize(), "configured transport solve synchronize");
	const auto solve_end = Clock::now();
	BlasHandle blas;
	double norm = 0.0;
	Check(cublasDnrm2(blas, vector_size, current.data(), 1, &norm), "configured transport final norm");
	if (!options.output.empty())
		WriteCudaConfiguredOutput<Fields>(
			options.output, system.fields, host.nodes, current);
	std::cout << "iga_cuda_solve system=" << system.name << " fields=" << Fields
		<< " steps=" << system.steps << " run_end_step=" << run_end_step
		<< " preprocess_s=" << std::chrono::duration<double>(preprocess_end-total_start).count()
		<< " assembly_s=" << assembly_seconds
		<< " solve_s=" << linear_seconds
		<< " time_loop_s=" << std::chrono::duration<double>(solve_end-solve_start).count()
		<< " total_iterations=" << total_iterations << " final_l2=" << norm
		<< " singular_diagonal_blocks=" << singular_blocks << " gpu_used_gib=" << peak_gpu_used << '\n';
	return 0;
}

int SolveConfigured(int argc, char** argv)
{
	if (argc < 4) throw std::runtime_error(
		"usage: iga_cuda solve DATABASE.ntiga CASE_DIR [SYSTEM] [OUTPUT] [VELOCITY] "
		"[--system NAME] [--output PATH] [--velocity PATH] [--output-every N] "
		"[--checkpoint PREFIX] [--checkpoint-every N] [--restart PREFIX] "
		"[--stop-after-step N]");
	const auto options = ParseCudaConfiguredTransportOptions(argc, argv);
	iga::Database database(argv[2]);
	const fs::path case_dir(argv[3]);
	const auto configuration = iga::ReadSimulationConfiguration((case_dir/"simulation_config.json").string());
	const auto system_name = !options.system.empty() ? options.system
		: iga::FirstLinearTransportSystem(configuration);
	const auto system = iga::CompileLinearSystem(configuration, system_name);
	if (system.velocity_source != "prescribed" && !options.velocity.empty())
		throw std::runtime_error("a positional VELOCITY file cannot override a configured snapshot_series");
	switch (system.fields.size()) {
	case 1: return ConfiguredTransport<1>(database, case_dir, configuration, system, options);
	case 2: return ConfiguredTransport<2>(database, case_dir, configuration, system, options);
	case 3: return ConfiguredTransport<3>(database, case_dir, configuration, system, options);
	case 4: return ConfiguredTransport<4>(database, case_dir, configuration, system, options);
	case 5: return ConfiguredTransport<5>(database, case_dir, configuration, system, options);
	case 6: return ConfiguredTransport<6>(database, case_dir, configuration, system, options);
	case 7: return ConfiguredTransport<7>(database, case_dir, configuration, system, options);
	case 8: return ConfiguredTransport<8>(database, case_dir, configuration, system, options);
	default: throw std::runtime_error("CUDA supports 1 through 8 coupled scalar fields per system");
	}
}
