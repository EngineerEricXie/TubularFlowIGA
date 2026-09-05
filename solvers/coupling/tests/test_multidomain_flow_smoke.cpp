#define main BifurcationSmokeFixtureMain
#include "test_bifurcation_coupling_smoke.cpp"
#undef main

namespace {

std::string ReadFile(const fs::path& path)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot read multi-island configuration");
	return std::string((std::istreambuf_iterator<char>(input)),
		std::istreambuf_iterator<char>());
}

void RequireLogContains(const fs::path& path, const std::string& text)
{
	if (ReadFile(path).find(text) == std::string::npos)
		throw std::runtime_error("expected diagnostic is missing from "+path.string());
}

double DeclaredSpeciesAmountBound(double scale)
{
	return 1.0e-10+1.0e-6*std::max(1.0e-5, scale);
}

struct SpeciesResidualMaxima {
	double edge = 0.0;
	double domain = 0.0;
	double global = 0.0;
};

SpeciesResidualMaxima species_residual_maxima;
fs::path executable_directory;

void WriteOneDDomain(std::ostream& output, const std::string& id,
	const std::string& case_directory, bool coupled_root, bool coupled_terminal,
	bool trailing_comma)
{
	output << "    {\"id\":\"" << id << "\",\"dimension\":\"1d\","
		"\"kind\":\"network_flow\",\"case\":\"" << case_directory << "\","
		"\"inlet_policy\":\"" << (coupled_root ? "coupled_root" : "configured_open_loop")
		<< "\",\"ports\":[{\"id\":\"root\",\"locator_kind\":\"runtime_port\","
		"\"locator\":\"root\",\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],"
		"\"requires\":[" << (coupled_root ? "\"flow_rate\"" : "") << "]},"
		"{\"id\":\"terminal\",\"locator_kind\":\"runtime_port\","
		"\"locator\":\"outlet:2\",\"provides\":[\"area\",\"flow_rate\","
		"\"mean_pressure\"],\"requires\":["
		<< (coupled_terminal ? "\"mean_pressure\"" : "") << "]}]}"
		<< (trailing_comma ? ",\n" : "\n");
}

void WriteThreeDDomain(std::ostream& output, const std::string& id,
	const std::string& case_directory, const std::string& database,
	bool trailing_comma)
{
	output << "    {\"id\":\"" << id << "\",\"dimension\":\"3d\","
		"\"kind\":\"body_fitted_iga_flow\",\"case\":\"" << case_directory
		<< "\",\"database\":\"" << database << "\",\"ports\":["
		"{\"id\":\"inlet\",\"locator_kind\":\"boundary_label\",\"locator\":\"1\","
		"\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],"
		"\"requires\":[\"flow_rate\"]},"
		"{\"id\":\"outlet_a\",\"locator_kind\":\"boundary_label\",\"locator\":\"2\","
		"\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],"
		"\"requires\":[\"mean_pressure\"]},"
		"{\"id\":\"outlet_b\",\"locator_kind\":\"boundary_label\",\"locator\":\"3\","
		"\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],"
		"\"requires\":[\"mean_pressure\"]},"
		"{\"id\":\"wall\",\"locator_kind\":\"boundary_label\",\"locator\":\"0\","
		"\"provides\":[\"flow_rate\"],\"requires\":[]}]}"
		<< (trailing_comma ? ",\n" : "\n");
}

void WriteEdge(std::ostream& output, const std::string& id,
	const std::string& first_domain, const std::string& first_port,
	const std::string& second_domain, const std::string& second_port,
	double initial_pressure, bool trailing_comma)
{
	output << "    {\"id\":\"" << id << "\",\"a\":{\"domain\":\""
		<< first_domain << "\",\"port\":\"" << first_port << "\"},"
		"\"b\":{\"domain\":\"" << second_domain << "\",\"port\":\""
		<< second_port << "\"},\"mode\":\"pressure_flow\","
		"\"initial_pressure_pa\":" << Number(initial_pressure) << "}"
		<< (trailing_comma ? ",\n" : "\n");
}

void WriteMultiIslandGraph(const fs::path& root, const std::string& execution,
	bool permuted = false)
{
	const bool explicit_mode = execution == "explicit";
	std::ofstream output(root/"simulation_config.json");
	output << "{\n  \"schema_version\":5,\"time\":{\"dt\":" << Number(kDt)
		<< ",\"steps\":" << kSteps << "},\"start_domain\":\"source\",\n"
		<< "  \"execution\":{\"kind\":\"" << execution << "\","
		"\"maximum_iterations\":" << (explicit_mode ? 1 : 100)
		<< ",\"pressure_relative_tolerance\":1e-6,\"pressure_reference_pa\":1,"
		"\"flow_relative_tolerance\":1e-10,\"relaxation_factor\":0.5,"
		"\"minimum_relaxation\":0.05,\"maximum_relaxation\":0.999},\n"
		<< "  \"domains\":[\n";
	if (!permuted) {
		WriteOneDDomain(output, "source", "source", false, true, true);
		WriteThreeDDomain(output, "island_a", "three_a", "a.ntiga", true);
		WriteOneDDomain(output, "bridge", "bridge", true, true, true);
		WriteOneDDomain(output, "leaf_a", "leaf_a", true, false, true);
		WriteThreeDDomain(output, "island_b", "three_b", "b.ntiga", true);
		WriteOneDDomain(output, "leaf_b", "leaf_b", true, false, true);
		WriteOneDDomain(output, "leaf_c", "leaf_c", true, false, false);
	} else {
		WriteOneDDomain(output, "leaf_c", "leaf_c", true, false, true);
		WriteOneDDomain(output, "leaf_b", "leaf_b", true, false, true);
		WriteThreeDDomain(output, "island_b", "three_b", "b.ntiga", true);
		WriteOneDDomain(output, "leaf_a", "leaf_a", true, false, true);
		WriteOneDDomain(output, "bridge", "bridge", true, true, true);
		WriteThreeDDomain(output, "island_a", "three_a", "a.ntiga", true);
		WriteOneDDomain(output, "source", "source", false, true, false);
	}
	output << "  ],\n  \"couplings\":[\n";
	if (!permuted) {
		WriteEdge(output, "e0_source_a", "source", "terminal", "island_a", "inlet", 0.0, true);
		WriteEdge(output, "e1_a_bridge", "island_a", "outlet_a", "bridge", "root", 0.02, true);
		WriteEdge(output, "e2_a_leaf", "island_a", "outlet_b", "leaf_a", "root", 0.01, true);
		WriteEdge(output, "e3_bridge_b", "bridge", "terminal", "island_b", "inlet", 0.0, true);
		WriteEdge(output, "e4_b_leaf", "island_b", "outlet_a", "leaf_b", "root", 0.015, true);
		WriteEdge(output, "e5_b_leaf", "island_b", "outlet_b", "leaf_c", "root", 0.005, false);
	} else {
		WriteEdge(output, "e5_b_leaf", "leaf_c", "root", "island_b", "outlet_b", 0.005, true);
		WriteEdge(output, "e4_b_leaf", "leaf_b", "root", "island_b", "outlet_a", 0.015, true);
		WriteEdge(output, "e3_bridge_b", "island_b", "inlet", "bridge", "terminal", 0.0, true);
		WriteEdge(output, "e2_a_leaf", "leaf_a", "root", "island_a", "outlet_b", 0.01, true);
		WriteEdge(output, "e1_a_bridge", "bridge", "root", "island_a", "outlet_a", 0.02, true);
		WriteEdge(output, "e0_source_a", "island_a", "inlet", "source", "terminal", 0.0, false);
	}
	output << "  ]\n}\n";
	if (!output) throw std::runtime_error("cannot write multi-island graph");
}

void WriteImmersedCase(const fs::path& directory)
{
	std::ofstream surface(directory/"surface.vtp");
	surface << R"xml(<?xml version="1.0"?><VTKFile type="PolyData" version="1.0" byte_order="LittleEndian"><PolyData><Piece NumberOfPoints="8" NumberOfPolys="12"><Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">.1 .1 .1 .9 .1 .1 .9 .9 .1 .1 .9 .1 .1 .1 .9 .9 .1 .9 .9 .9 .9 .1 .9 .9</DataArray></Points><Polys><DataArray type="Int32" Name="connectivity" format="ascii">0 2 1 0 3 2 4 5 6 4 6 7 0 1 5 0 5 4 1 2 6 1 6 5 2 3 7 2 7 6 3 0 4 3 4 7</DataArray><DataArray type="Int32" Name="offsets" format="ascii">3 6 9 12 15 18 21 24 27 30 33 36</DataArray></Polys><CellData Scalars="boundary_id"><DataArray type="UInt32" Name="boundary_id" format="ascii">1 1 2 2 0 0 0 0 0 0 0 0</DataArray></CellData></Piece></PolyData></VTKFile>)xml";
	std::ofstream config(directory/"simulation_config.json");
	config << "{\"schema_version\":3,\"dimension\":\"3d\",\"fields\":[{\"name\":\"velocity\",\"kind\":\"vector3\"},{\"name\":\"pressure\",\"kind\":\"pressure\"}],\"time\":{\"dt\":" << Number(kDt) << ",\"steps\":" << kSteps << "},\"equation_systems\":[{\"name\":\"flow\",\"kind\":\"navier_stokes\",\"unknowns\":[\"velocity\",\"pressure\"],\"viscosity\":" << Number(kViscosity) << ",\"density\":" << Number(kDensity) << ",\"time_integration\":\"steady\"}],\"boundaries\":[]}";
	std::ofstream geometry(directory/"immersed_geometry.json");
	geometry << R"json({"surface":"surface.vtp","grid":{"lower_m":[0.012345678901234567,0.012345678901234567,0.012345678901234567],"upper_m":[0.98765432109876543,0.98765432109876543,0.98765432109876543],"cells":[3,3,3]},"volume_quadrature":{"storage":"compact","max_depth":4,"max_nodes":500000,"max_leaves":500000,"max_points":3000000,"max_records":3000000,"max_retained_bytes":30000000,"max_logical_points":3000000},"surface_quadrature":{"max_candidates":500000,"max_fragments":500000,"max_points":3000000,"max_exact_limbs":512},"ghost_penalty":{"gamma_u":0.01,"gamma_p":0.01,"max_faces":500000,"max_quadrature_points":3000000,"max_trace_entries":128},"wall_labels":[0],"runtime":{"lu_pivot_shift":1e-12,"ports":[{"id":"inlet","boundary_label":1,"control_mode":"flow_rate","value":-0.001},{"id":"outlet","boundary_label":2,"control_mode":"pressure","value":0}]}})json";
	if (!surface || !config || !geometry) throw std::runtime_error("cannot write immersed fixture");
}

void WriteImmersedGraph(const fs::path& root)
{
	std::ofstream output(root/"simulation_config.json", std::ios::trunc);
	output << "{\"schema_version\":5,\"time\":{\"dt\":" << Number(kDt)
		<< ",\"steps\":" << kSteps << "},\"start_domain\":\"source\",\"execution\":{\"kind\":\"explicit\",\"maximum_iterations\":1,\"pressure_relative_tolerance\":1e-6,\"pressure_reference_pa\":1,\"flow_relative_tolerance\":1e-10,\"relaxation_factor\":0.5,\"minimum_relaxation\":0.05,\"maximum_relaxation\":0.999},\"domains\":[";
	WriteOneDDomain(output, "source", "immersed_source", false, true, true);
	output << "{\"id\":\"immersed\",\"dimension\":\"3d\",\"kind\":\"three_d_immersed_flow\",\"case\":\"immersed\",\"ports\":[{\"id\":\"inlet\",\"locator_kind\":\"boundary_label\",\"locator\":\"1\",\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],\"requires\":[\"flow_rate\"]},{\"id\":\"outlet\",\"locator_kind\":\"boundary_label\",\"locator\":\"2\",\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],\"requires\":[\"mean_pressure\"]}]},";
	WriteOneDDomain(output, "sink", "immersed_sink", true, false, false);
	output << "],\"couplings\":[{\"id\":\"in\",\"a\":{\"domain\":\"source\",\"port\":\"terminal\"},\"b\":{\"domain\":\"immersed\",\"port\":\"inlet\"},\"mode\":\"pressure_flow\",\"initial_pressure_pa\":0},{\"id\":\"out\",\"a\":{\"domain\":\"immersed\",\"port\":\"outlet\"},\"b\":{\"domain\":\"sink\",\"port\":\"root\"},\"mode\":\"pressure_flow\",\"initial_pressure_pa\":0}]}";
	if (!output) throw std::runtime_error("cannot write immersed graph");
}

void ValidateImmersedManifest(const fs::path& output)
{
	const std::string manifest = ReadFile(output/"graph_binding_manifest.json");
	const std::string hash_key = "\"surface_hash\":\"";
	const auto hash_begin = manifest.find(hash_key);
	if (hash_begin == std::string::npos)
		throw std::runtime_error("immersed manifest is missing the surface hash");
	const auto hash_end = manifest.find('"', hash_begin+hash_key.size());
	if (hash_end == std::string::npos || hash_end-(hash_begin+hash_key.size()) != 64)
		throw std::runtime_error("immersed manifest surface hash is invalid");
	for (const auto& text : {"\"kind\":\"three_d_immersed_flow\"",
		"\"grid\":{\"lower_m\":[0.012345678901234567,0.012345678901234567,0.012345678901234567],\"upper_m\":[0.98765432109876539,0.98765432109876539,0.98765432109876539],\"cells\":[3,3,3]}",
		"\"catalog_audit\":{\"active_cells\":", "\"volume_points\":",
		"\"surface_points\":", "\"ghost_faces\":"})
		if (manifest.find(text) == std::string::npos)
			throw std::runtime_error("immersed manifest audit is incomplete");
}

void WriteOneDTransportCase(const fs::path& directory, double inlet_flow,
	const std::string& first_field, const std::string& second_field)
{
	std::ofstream network(directory/"tree.swc");
	network << std::setprecision(17) << "1 2 0 0 0 " << kRadius << " -1\n"
		<< "2 2 1 0 0 " << kRadius << " 1\n";
	std::ofstream config(directory/"simulation_config.json");
	config << "{\n  \"schema_version\":3,\"dimension\":\"1d\","
		"\"simulation_scope\":{\"mode\":\"flow_only\"},\n"
		<< "  \"coupling\":{\"scheme\":\"explicit_staggered\","
			"\"flow_epsilon_m3_s\":0.01},\n"
		<< "  \"geometry\":{\"kind\":\"swc_network\",\"file\":\"tree.swc\","
			"\"length_scale_to_m\":1},\n"
		<< "  \"fields\":[{\"name\":\"area\",\"kind\":\"scalar\"},"
			"{\"name\":\"flow_rate\",\"kind\":\"scalar\"},"
			"{\"name\":\"pressure\",\"kind\":\"pressure\"},"
			"{\"name\":\"" << first_field << "\",\"kind\":\"scalar\",\"initial_value\":1.25},"
			"{\"name\":\"" << second_field << "\",\"kind\":\"scalar\",\"initial_value\":3.5}],\n"
		<< "  \"time\":{\"dt\":" << Number(kDt) << ",\"steps\":" << kSteps
		<< ",\"output_every\":1},\n"
		<< "  \"temporal_functions\":[{\"name\":\"inlet_flow\",\"kind\":\"constant\","
			"\"units\":\"m3/s\",\"value\":" << Number(inlet_flow) << "}],\n"
		<< "  \"equation_systems\":["
			"{\"name\":\"flow\",\"kind\":\"network_flow_1d\","
			"\"unknowns\":[\"area\",\"flow_rate\",\"pressure\"],\"model\":\"rigid\","
			"\"scheme\":\"steady_poiseuille\",\"dynamic_viscosity\":" << Number(kViscosity)
		<< ",\"density\":" << Number(kDensity)
		<< ",\"discretization\":{\"cells_per_segment\":1}},"
			"{\"name\":\"transport\",\"kind\":\"network_transport_1d\","
			"\"unknowns\":[\"" << first_field << "\",\"" << second_field
		<< "\"],\"flow_system\":\"flow\",\"species\":["
			"{\"field\":\"" << first_field << "\",\"diffusivity\":0},"
			"{\"field\":\"" << second_field << "\",\"diffusivity\":0}]}],\n"
		<< "  \"boundaries\":["
			"{\"name\":\"inlet\",\"role\":\"inlet\",\"node_ids\":[1],\"conditions\":["
			"{\"field\":\"flow_rate\",\"type\":\"dirichlet\",\"quantity\":\"flow_rate\","
			"\"waveform\":\"inlet_flow\"},"
			"{\"field\":\"" << first_field << "\",\"type\":\"dirichlet\",\"value\":1.25},"
			"{\"field\":\"" << second_field << "\",\"type\":\"dirichlet\",\"value\":3.5}]},"
			"{\"name\":\"outlet\",\"role\":\"outlet\",\"node_ids\":[2],\"conditions\":["
			"{\"field\":\"pressure\",\"type\":\"pressure\",\"value\":0}]}]\n}\n";
	if (!network || !config) throw std::runtime_error("cannot write schema-v6 1D case");
}

void WriteThreeDTransportCase(const fs::path& directory, bool prescribed = true)
{
	WriteThreeDCase(directory);
	std::ofstream config(directory/"simulation_config.json", std::ios::trunc);
	config << "{\n  \"schema_version\":3,\"dimension\":\"3d\",\n"
		<< "  \"simulation_scope\":{\"mode\":\"flow_only\"},\n"
		<< "  \"coupling\":{\"scheme\":\"explicit_staggered\",\"flow_epsilon_m3_s\":1e-14,"
			"\"three_d_ports\":{\"inlet_label\":1,\"outlet_labels\":[2,3]}},\n"
		<< (prescribed ? "" : "  \"velocity_sources\":[{\"name\":\"velocity\","
			"\"kind\":\"snapshot_series\",\"manifest\":\"velocity_series.json\","
			"\"interpolation\":\"linear\",\"out_of_range\":\"hold\"}],\n")
		<< "  \"fields\":[{\"name\":\"velocity\",\"kind\":\"vector3\"},"
			"{\"name\":\"pressure\",\"kind\":\"pressure\"},"
			"{\"name\":\"three_red\",\"kind\":\"scalar\",\"initial_value\":1.25},"
			"{\"name\":\"three_blue\",\"kind\":\"scalar\",\"initial_value\":3.5}],\n"
		<< "  \"time\":{\"dt\":" << Number(kDt) << ",\"steps\":" << kSteps << "},\n"
		<< "  \"equation_systems\":["
			"{\"name\":\"flow\",\"kind\":\"navier_stokes\","
			"\"unknowns\":[\"velocity\",\"pressure\"],\"viscosity\":" << Number(kViscosity)
		<< ",\"density\":" << Number(kDensity)
		<< ",\"time_integration\":\"backward_euler\"},"
			"{\"name\":\"transport\",\"kind\":\"linear_transport\","
			"\"unknowns\":[\"three_red\",\"three_blue\"],\"terms\":["
			"{\"operator\":\"time_derivative\",\"equation\":\"three_red\"},"
			"{\"operator\":\"advection\",\"equation\":\"three_red\",\"velocity\":\""
		<< (prescribed ? "prescribed" : "velocity") << "\"},"
			"{\"operator\":\"time_derivative\",\"equation\":\"three_blue\"},"
			"{\"operator\":\"advection\",\"equation\":\"three_blue\",\"velocity\":\""
		<< (prescribed ? "prescribed" : "velocity") << "\"}]}],\n"
		<< "  \"boundaries\":[\n"
			"    {\"label\":0,\"name\":\"wall\",\"conditions\":["
			"{\"field\":\"velocity\",\"type\":\"dirichlet\",\"value\":[0,0,0]},"
			"{\"field\":\"three_red\",\"type\":\"no_flux\"},"
			"{\"field\":\"three_blue\",\"type\":\"no_flux\"}]},\n"
			"    {\"label\":1,\"name\":\"inlet\",\"conditions\":["
			"{\"field\":\"velocity\",\"type\":\"dirichlet\","
			"\"profile\":\"initial_velocityfield.txt\",\"scale\":1},"
			"{\"field\":\"three_red\",\"type\":\"dirichlet\",\"value\":1.25},"
			"{\"field\":\"three_blue\",\"type\":\"dirichlet\",\"value\":3.5}]},\n"
			"    {\"label\":2,\"name\":\"outlet_a\",\"conditions\":["
			"{\"field\":\"pressure\",\"type\":\"pressure_traction\",\"value\":0},"
			"{\"field\":\"three_red\",\"type\":\"advective_outflow\"},"
			"{\"field\":\"three_blue\",\"type\":\"advective_outflow\"}]},\n"
			"    {\"label\":3,\"name\":\"outlet_b\",\"conditions\":["
			"{\"field\":\"pressure\",\"type\":\"pressure_traction\",\"value\":0},"
			"{\"field\":\"three_red\",\"type\":\"advective_outflow\"},"
			"{\"field\":\"three_blue\",\"type\":\"advective_outflow\"}]}\n"
			"  ]\n}\n";
	if (!config) throw std::runtime_error("cannot write schema-v6 3D case");
}

void WriteSpeciesPort(std::ostream& output, const std::string& id,
	const std::string& locator_kind, const std::string& locator,
	const std::string& hydraulic_requirement)
{
	output << "{\"id\":\"" << id << "\",\"locator_kind\":\"" << locator_kind
		<< "\",\"locator\":\"" << locator << "\","
			"\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\","
			"\"species_concentration\",\"species_flux\"],\"requires\":["
		<< hydraulic_requirement
		<< (hydraulic_requirement.empty() ? "" : ",")
		<< "\"species_concentration\",\"species_flux\"],"
			"\"species\":[\"red_logical\",\"blue_logical\"]}";
}

void WriteSpeciesObservationPort(std::ostream& output, const std::string& id,
	const std::string& locator_kind, const std::string& locator,
	bool hydraulic_state = true)
{
	output << "{\"id\":\"" << id << "\",\"locator_kind\":\"" << locator_kind
		<< "\",\"locator\":\"" << locator << "\",\"provides\":[";
	if (hydraulic_state)
		output << "\"area\",\"flow_rate\",\"mean_pressure\",";
	else output << "\"flow_rate\",";
	output << "\"species_concentration\",\"species_flux\"],\"requires\":[],"
		"\"species\":[\"red_logical\",\"blue_logical\"]}";
}

void WriteSpeciesGraph(const fs::path& root, const std::string& database,
	bool permuted = false)
{
	std::ofstream output(root/"simulation_config.json", std::ios::trunc);
	output << "{\n  \"schema_version\":6,"
		"\"species\":[{\"id\":\"red_logical\",\"concentration_unit\":\"kg/m^3\"},"
		"{\"id\":\"blue_logical\",\"concentration_unit\":\"mol/m^3\"}],"
		"\"time\":{\"dt\":" << Number(kDt) << ",\"steps\":" << kSteps
		<< "},\"start_domain\":\"source\",\n"
			"  \"execution\":{\"kind\":\"explicit\","
			"\"species_routing\":{\"flow_switch_m3_s\":1e-6,"
			"\"flow_absolute_tolerance_m3_s\":1e-12,\"flow_relative_tolerance\":1e-8},"
			"\"species_amount_tolerances\":{"
			"\"red_logical\":{\"absolute_tolerance\":1e-10,\"reference_amount\":1e-5,"
			"\"relative_tolerance\":1e-6},"
			"\"blue_logical\":{\"absolute_tolerance\":1e-10,\"reference_amount\":1e-5,"
			"\"relative_tolerance\":1e-6}}},\n  \"domains\":[\n";
	const auto one_d = [&](const std::string& id, const std::string& case_directory,
		bool coupled, bool comma) {
		output << "    {\"id\":\"" << id << "\",\"dimension\":\"1d\","
			"\"kind\":\"network_flow\",\"case\":\"" << case_directory << "\","
			"\"inlet_policy\":\"" << (coupled ? "coupled_root" : "configured_open_loop")
			<< "\",\"species_bindings\":{\"red_logical\":\"" << id << "_red\","
			"\"blue_logical\":\"" << id << "_blue\"},\"ports\":[";
		if (coupled)
			WriteSpeciesPort(output, "root", "runtime_port", "root", "\"flow_rate\"");
		else WriteSpeciesObservationPort(output, "root", "runtime_port", "root");
		output << ',';
		if (!coupled)
			WriteSpeciesPort(output, "terminal", "runtime_port", "outlet:2",
				"\"mean_pressure\"");
		else WriteSpeciesObservationPort(output, "terminal", "runtime_port", "outlet:2");
		output << "]}" << (comma ? ",\n" : "\n");
	};
	const auto three_d = [&](bool comma) {
		output << "    {\"id\":\"junction\",\"dimension\":\"3d\","
			"\"kind\":\"body_fitted_iga_flow\",\"case\":\"junction\",\"database\":\""
			<< database << "\",\"species_bindings\":{\"red_logical\":\"three_red\","
			"\"blue_logical\":\"three_blue\"},\"ports\":[";
		WriteSpeciesPort(output, "inlet", "boundary_label", "1", "\"flow_rate\"");
		output << ',';
		WriteSpeciesPort(output, "outlet_a", "boundary_label", "2", "\"mean_pressure\"");
		output << ',';
		WriteSpeciesPort(output, "outlet_b", "boundary_label", "3", "\"mean_pressure\"");
		output << ',';
		WriteSpeciesObservationPort(output, "wall", "boundary_label", "0", false);
		output << "]}";
		output << (comma ? ",\n" : "\n");
	};
	if (!permuted) {
		one_d("source", "species_source", false, true);
		three_d(true);
		one_d("leaf_a", "species_leaf_a", true, true);
		one_d("leaf_b", "species_leaf_b", true, false);
	} else {
		one_d("leaf_b", "species_leaf_b", true, true);
		one_d("leaf_a", "species_leaf_a", true, true);
		three_d(true);
		one_d("source", "species_source", false, false);
	}
	output << "  ],\n  \"couplings\":[\n";
	const auto edge = [&](const std::string& id, const std::string& first_domain,
		const std::string& first_port, const std::string& second_domain,
		const std::string& second_port, double pressure, bool comma) {
		output << "    {\"id\":\"" << id << "\",\"a\":{\"domain\":\""
			<< first_domain << "\",\"port\":\"" << first_port << "\"},"
			"\"b\":{\"domain\":\"" << second_domain << "\",\"port\":\""
			<< second_port << "\"},\"mode\":\"pressure_flow\","
			"\"initial_pressure_pa\":" << Number(pressure)
			<< ",\"species\":[\"blue_logical\",\"red_logical\"]}"
			<< (comma ? ",\n" : "\n");
	};
	if (!permuted) {
		edge("edge_in", "source", "terminal", "junction", "inlet", 0.0, true);
		edge("edge_a", "junction", "outlet_a", "leaf_a", "root", 0.02, true);
		edge("edge_b", "junction", "outlet_b", "leaf_b", "root", 0.01, false);
	} else {
		edge("edge_b", "leaf_b", "root", "junction", "outlet_b", 0.01, true);
		edge("edge_a", "leaf_a", "root", "junction", "outlet_a", 0.02, true);
		edge("edge_in", "junction", "inlet", "source", "terminal", 0.0, false);
	}
	output << "  ]\n}\n";
	if (!output) throw std::runtime_error("cannot write schema-v6 graph");
}

int RunMultidomain(const fs::path& root, const fs::path& output,
	const std::string& launcher = {}, const std::string& prefix = {},
	bool transport = false)
{
	const auto log = root/(output.filename().string()+".log");
	const std::string linear_solver = transport
		? (launcher.empty() ? " -ksp_type gmres -pc_type lu -ksp_rtol 1e-12"
			: " -ksp_type gmres -pc_type bjacobi -sub_pc_type lu -ksp_rtol 1e-12")
		: (launcher.empty()
		? " -ksp_type preonly -pc_type lu"
		: " -ksp_type gmres -pc_type bjacobi -sub_pc_type lu -ksp_rtol 1e-12");
	return std::system((prefix+launcher+Quote(executable_directory/"iga_multidomain_flow")+" --graph-case "+Quote(root)
		+" --output-dir "+Quote(output)+linear_solver+">"
		+Quote(log)+" 2>&1").c_str());
}

int RunBifurcation(const fs::path& root, const fs::path& output)
{
	const auto log = root/(output.filename().string()+".log");
	return std::system((Quote(executable_directory/"iga_1d_3d_bifurcation")+" --graph-case "+Quote(root)
		+" --output-dir "+Quote(output)+" -ksp_type preonly -pc_type lu >"
		+Quote(log)+" 2>&1").c_str());
}

void ValidateMultidomain(const fs::path& output, const std::string& execution)
{
	if (!fs::is_regular_file(output/"graph_binding_manifest.json"))
		throw std::runtime_error("multi-island completion marker is missing");
	const auto steps = ReadCsv(output/"pressure_flow_steps.csv");
	const auto edges = ReadCsv(output/"pressure_flow_edges.csv");
	const auto ports = ReadCsv(output/"pressure_flow_ports.csv");
	const auto balances = ReadCsv(output/"pressure_flow_domain_balances.csv");
	const auto initialization = ReadCsv(output/"one_d_initialization.csv");
	if (steps.size() != kSteps || edges.size() != 6*kSteps
		|| ports.size() != 18*kSteps || balances.size() != 2*kSteps
		|| initialization.size() != 5)
		throw std::runtime_error("multi-island long-form row count is invalid");
	std::map<std::string, int> edge_counts;
	for (const auto& edge : edges) {
		++edge_counts[edge.at("edge_id")];
		if (std::abs(Value(edge, "flow_residual_m3_s")) > 1.0e-13
			|| Value(edge, "normalized_flow_residual") > 1.0e-10)
			throw std::runtime_error("multi-island edge conservation failed");
	}
	if (edge_counts.size() != 6)
		throw std::runtime_error("multi-island edge coverage is incomplete");
	std::map<std::string, int> balance_counts;
	for (const auto& row : balances) {
		++balance_counts[row.at("domain_id")];
		if (Value(row, "normalized_mass_imbalance") > 1.0e-3)
			throw std::runtime_error("multi-island 3D mass balance failed");
	}
	if (balance_counts != std::map<std::string, int>{{"island_a", kSteps},
		{"island_b", kSteps}})
		throw std::runtime_error("multi-island domain-balance coverage is invalid");
	for (const auto& row : steps) {
		const double external = std::abs(Value(row, "external_outward_flow_m3_s"));
		double scale = 0.0;
		const int step = static_cast<int>(Value(row, "step"));
		for (const auto& port : ports)
			if (static_cast<int>(Value(port, "step")) == step
				&& ((port.at("domain_id") == "source" && port.at("port_id") == "root")
					|| (port.at("domain_id").find("leaf_") == 0
						&& port.at("port_id") == "terminal")))
				scale += std::abs(Value(port, "outward_flow_m3_s"));
		if (external/std::max(1.0e-14, scale) > 1.0e-10)
			throw std::runtime_error("multi-island global external conservation failed");
	}
	if (execution == "aitken") {
		const auto iterations = ReadCsv(output/"pressure_flow_iterations.csv");
		if (iterations.size() <= edges.size())
			throw std::runtime_error("multi-island Aitken run did not exercise a retry");
	}
}

void ValidateSpeciesMultidomain(const fs::path& output)
{
	if (!fs::is_regular_file(output/"graph_binding_manifest.json"))
		throw std::runtime_error("schema-v6 completion marker is missing");
	const std::string marker = ReadFile(output/"graph_binding_manifest.json");
	for (const auto& text : {"\"schema_version\": 6", "\"species_mode\": true",
		"\"red_logical\"", "\"blue_logical\"", "\"concentration\":\"kg/m^3\"",
		"\"integrated_amount\":\"(kg/m^3)*m^3\"",
		"\"outward_rate\":\"(mol/m^3)*m^3/s\""})
		if (marker.find(text) == std::string::npos)
			throw std::runtime_error("schema-v6 species units/marker are incomplete");
	const auto hydraulic = ReadCsv(output/"species_hydraulic_steps.csv");
	const auto edges = ReadCsv(output/"species_edge_amounts.csv");
	const auto domains = ReadCsv(output/"species_domain_accounting.csv");
	const auto global = ReadCsv(output/"species_global_balance.csv");
	const auto ports = ReadCsv(output/"species_logical_ports.csv");
	if (hydraulic.size() != kSteps || edges.size() != 6*kSteps
		|| domains.size() != 8*kSteps || global.size() != 2*kSteps
		|| ports.size() != 20*kSteps)
		throw std::runtime_error("schema-v6 long-form row coverage is invalid");
	struct ExpectedEdge {
		const char* edge_id;
		const char* species_id;
		const char* first_domain;
		const char* first_port;
		const char* second_domain;
		const char* second_port;
	};
	const std::vector<ExpectedEdge> expected_edges{
		{"edge_a", "blue_logical", "junction", "outlet_a", "leaf_a", "root"},
		{"edge_a", "red_logical", "junction", "outlet_a", "leaf_a", "root"},
		{"edge_b", "blue_logical", "junction", "outlet_b", "leaf_b", "root"},
		{"edge_b", "red_logical", "junction", "outlet_b", "leaf_b", "root"},
		{"edge_in", "blue_logical", "junction", "inlet", "source", "terminal"},
		{"edge_in", "red_logical", "junction", "inlet", "source", "terminal"}};
	for (int step = 0; step < kSteps; ++step)
		for (std::size_t index = 0; index < expected_edges.size(); ++index) {
			const auto& row = edges.at(static_cast<std::size_t>(step)*expected_edges.size()+index);
			const auto& expected = expected_edges[index];
			if (row.at("edge_id") != expected.edge_id
				|| row.at("species_id") != expected.species_id
				|| row.at("first_domain_id") != expected.first_domain
				|| row.at("first_port_id") != expected.first_port
				|| row.at("second_domain_id") != expected.second_domain
				|| row.at("second_port_id") != expected.second_port)
				throw std::runtime_error("schema-v6 edge/species diagnostics are not canonical");
			const double amount_scale = std::max(std::abs(Value(row, "first_outward_amount")),
				std::abs(Value(row, "second_outward_amount")));
			species_residual_maxima.edge = std::max(species_residual_maxima.edge,
				std::abs(Value(row, "amount_residual")));
			if (std::abs(Value(row, "amount_residual")) > DeclaredSpeciesAmountBound(amount_scale))
				throw std::runtime_error("schema-v6 edge integrated amount gate failed");
		}
	for (const auto& row : domains) {
		const double m0 = Value(row, "M0_amount");
		const double m1 = Value(row, "M1_amount");
		const double source = Value(row, "source_amount");
		const double outward = Value(row, "total_outward_amount");
		const double recomputed = m1-m0+outward-source;
		species_residual_maxima.domain = std::max(species_residual_maxima.domain,
			std::abs(recomputed));
		if (std::abs(recomputed-Value(row, "recomputed_residual")) > 1.0e-18
			|| std::abs(recomputed) > DeclaredSpeciesAmountBound(
				std::max({std::abs(m1-m0), std::abs(source), std::abs(outward)})))
			throw std::runtime_error("schema-v6 domain integrated amount gate failed");
	}
	for (const auto& row : global) {
		const double residual = Value(row, "M1_amount")-Value(row, "M0_amount")
			+Value(row, "outward_amount")-Value(row, "source_amount");
		species_residual_maxima.global = std::max(species_residual_maxima.global,
			std::abs(residual));
		if (std::abs(residual-Value(row, "global_residual")) > 1.0e-18
			|| std::abs(residual) > DeclaredSpeciesAmountBound(
				Value(row, "gross_activity")))
			throw std::runtime_error("schema-v6 global integrated amount gate failed");
		if (!(Value(row, "gross_activity") > 0.0))
			throw std::runtime_error("schema-v6 global gross activity is not positive");
	}
	const auto concentration = [&](const std::string& domain, const std::string& port,
		const std::string& species, int step) {
		for (const auto& row : ports)
			if (row.at("domain_id") == domain && row.at("port_id") == port
				&& row.at("species_id") == species
				&& static_cast<int>(Value(row, "step")) == step)
				return Value(row, "concentration");
		throw std::runtime_error("schema-v6 logical transport port is missing");
	};
	for (int step = 1; step <= kSteps; ++step)
		for (const auto& species : {std::make_pair("red_logical", 1.25),
			std::make_pair("blue_logical", 3.5)}) {
			const double source = concentration("source", "terminal", species.first, step);
			for (const auto& receiver : {std::make_pair("junction", "inlet"),
				std::make_pair("leaf_a", "root"), std::make_pair("leaf_b", "root")})
				if (std::abs(concentration(receiver.first, receiver.second, species.first, step)
					-source) > 1.0e-10)
					throw std::runtime_error("schema-v6 receiver concentration differs from donor");
		}
	for (const auto& domain : domains) {
		double end_step_rate = 0.0;
		for (const auto& port : ports)
			if (port.at("domain_id") == domain.at("domain_id")
				&& port.at("species_id") == domain.at("species_id")
				&& Value(port, "step") == Value(domain, "step"))
				end_step_rate += Value(port, "end_step_outward_species_flux");
		if (std::abs(kDt*end_step_rate-Value(domain, "total_outward_amount")) > 1.0e-10)
			throw std::runtime_error(
				"schema-v6 constant-state rate/amount diagnostic is inconsistent");
	}
	for (const auto& species : {"red_logical", "blue_logical"})
		for (int step = 1; step <= kSteps; ++step) {
			double source_root_rate = 0.0;
			double leaf_terminal_rate = 0.0;
			double wall_rate = 0.0;
			bool source_root_seen = false;
			bool wall_seen = false;
			int leaf_terminal_count = 0;
			for (const auto& row : ports) {
				if (row.at("species_id") != species
					|| static_cast<int>(Value(row, "step")) != step) continue;
				if (row.at("domain_id") == "source" && row.at("port_id") == "root") {
					source_root_seen = true;
					source_root_rate = Value(row, "end_step_outward_species_flux");
				}
				if ((row.at("domain_id") == "leaf_a" || row.at("domain_id") == "leaf_b")
					&& row.at("port_id") == "terminal") {
					++leaf_terminal_count;
					leaf_terminal_rate += Value(row, "end_step_outward_species_flux");
				}
				if (row.at("domain_id") == "junction" && row.at("port_id") == "wall") {
					wall_seen = true;
					wall_rate = Value(row, "end_step_outward_species_flux");
				}
			}
			if (!source_root_seen || leaf_terminal_count != 2
				|| !wall_seen || !(source_root_rate < 0.0) || !(leaf_terminal_rate > 0.0)
				|| std::abs(wall_rate) > 1.0e-10)
				throw std::runtime_error(
					"schema-v6 external species observation ports are incomplete or unsigned");
		}
}

constexpr int kZeroDClockSteps = 8;
constexpr double kZeroDClockDt = 0.1;

void WriteZeroDModel(const fs::path& directory, bool source)
{
	std::ofstream output(directory/"zero_d_model.json");
	if (source) {
		output << R"json({"role":"source_reservoir","capacitance_m3_pa":0.1,
"resistance_pa_s_m3":100,"initial_pressure_pa":0.1,"prescribed_flow_m3_s":0.001})json";
	} else {
		output << R"json({"role":"terminal_rcr","proximal_resistance_pa_s_m3":0,
"distal_resistance_pa_s_m3":100,"capacitance_m3_pa":0.1,"distal_pressure_pa":0,
"initial_pressure_pa":0})json";
	}
	if (!output) throw std::runtime_error("cannot write 0D clock fixture model");
}

void WriteZeroDThreeDCase(const fs::path& directory)
{
	WriteThreeDCase(directory);
	std::string configuration = ReadFile(directory/"simulation_config.json");
	const std::string previous = "\"time\":{\"dt\":"+Number(kDt)
		+",\"steps\":"+std::to_string(kSteps)+"}";
	const std::string replacement = "\"time\":{\"dt\":"+Number(kZeroDClockDt)
		+",\"steps\":"+std::to_string(kZeroDClockSteps)+"}";
	const auto position = configuration.find(previous);
	if (position == std::string::npos)
		throw std::runtime_error("0D clock fixture cannot update 3D time grid");
	configuration.replace(position, previous.size(), replacement);
	std::ofstream output(directory/"simulation_config.json", std::ios::trunc);
	output << configuration;
	if (!output) throw std::runtime_error("cannot write 0D clock fixture 3D time grid");
}

void WriteZeroDClockGraph(const fs::path& root)
{
	std::ofstream output(root/"simulation_config.json", std::ios::trunc);
	output << "{\n  \"schema_version\":5,\"time\":{\"dt\":" << Number(kZeroDClockDt)
		<< ",\"steps\":" << kZeroDClockSteps << "},\"start_domain\":\"source_0d\",\n"
		<< "  \"execution\":{\"kind\":\"explicit\"},\n  \"domains\":[\n"
		<< "    {\"id\":\"source_0d\",\"dimension\":\"0d\",\"kind\":\"zero_d_flow\","
		<< "\"case\":\"zero_source\",\"zero_d_model\":\"zero_d_model.json\",\"ports\":[{"
		<< "\"id\":\"port\",\"locator_kind\":\"zero_d_port\",\"locator\":\"port\","
		<< "\"provides\":[\"mean_pressure\",\"flow_rate\"],\"requires\":[\"mean_pressure\"]}]},\n"
		<< "    {\"id\":\"junction\",\"dimension\":\"3d\",\"kind\":\"body_fitted_iga_flow\","
		<< "\"case\":\"zero_junction\",\"database\":\"zero_junction.ntiga\",\"ports\":["
		<< "{\"id\":\"inlet\",\"locator_kind\":\"boundary_label\",\"locator\":\"1\","
		<< "\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],\"requires\":[\"flow_rate\"]},"
		<< "{\"id\":\"outlet_a\",\"locator_kind\":\"boundary_label\",\"locator\":\"2\","
		<< "\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],\"requires\":[\"mean_pressure\"]},"
		<< "{\"id\":\"outlet_b\",\"locator_kind\":\"boundary_label\",\"locator\":\"3\","
		<< "\"provides\":[\"area\",\"flow_rate\",\"mean_pressure\"],\"requires\":[]},"
		<< "{\"id\":\"wall\",\"locator_kind\":\"boundary_label\",\"locator\":\"0\","
		<< "\"provides\":[\"flow_rate\"],\"requires\":[]}]},\n"
		<< "    {\"id\":\"terminal_0d\",\"dimension\":\"0d\",\"kind\":\"zero_d_flow\","
		<< "\"case\":\"zero_terminal\",\"zero_d_model\":\"zero_d_model.json\",\"ports\":[{"
		<< "\"id\":\"port\",\"locator_kind\":\"zero_d_port\",\"locator\":\"port\","
		<< "\"provides\":[\"mean_pressure\",\"flow_rate\"],\"requires\":[\"flow_rate\"]}]}\n"
		<< "  ],\n  \"couplings\":["
		<< "{\"id\":\"source_to_junction\",\"a\":{\"domain\":\"source_0d\",\"port\":\"port\"},"
		<< "\"b\":{\"domain\":\"junction\",\"port\":\"inlet\"},\"mode\":\"pressure_flow\",\"initial_pressure_pa\":0},"
		<< "{\"id\":\"junction_to_terminal\",\"a\":{\"domain\":\"junction\",\"port\":\"outlet_a\"},"
		<< "\"b\":{\"domain\":\"terminal_0d\",\"port\":\"port\"},\"mode\":\"pressure_flow\",\"initial_pressure_pa\":0}]\n}\n";
	if (!output) throw std::runtime_error("cannot write 0D clock graph");
}

void RequireManifestDigest(const std::string& manifest, const std::string& key)
{
	const auto start = manifest.find(key);
	if (start == std::string::npos)
		throw std::runtime_error("0D manifest is missing "+key);
	const auto value = start+key.size();
	const auto end = manifest.find('"', value);
	if (end == std::string::npos || end-value != 64)
		throw std::runtime_error("0D manifest digest is malformed");
}

std::size_t CountOccurrences(const std::string& text, const std::string& needle)
{
	std::size_t count = 0;
	for (std::size_t position = text.find(needle); position != std::string::npos;
		position = text.find(needle, position+needle.size()))
		++count;
	return count;
}

void ValidateZeroDClockRun(const fs::path& output)
{
	if (!fs::is_regular_file(output/"graph_binding_manifest.json"))
		throw std::runtime_error("0D clock fixture did not complete cleanly");
	const std::string manifest = ReadFile(output/"graph_binding_manifest.json");
	for (const auto& text : {"\"kind\":\"zero_d_flow\"", "\"model_role\":\"source_reservoir\"",
		"\"model_role\":\"terminal_rcr\""})
		if (manifest.find(text) == std::string::npos)
			throw std::runtime_error("0D manifest model role is incomplete");
	for (const auto& key : {"\"model_identity_sha256\":\"",
		"\"final_state_identity_sha256\":\"",
		"\"final_step_accounting_identity_sha256\":\""}) {
		if (CountOccurrences(manifest, key) != 2)
			throw std::runtime_error("0D manifest is missing a domain identity");
		RequireManifestDigest(manifest, key);
	}
	const auto steps = ReadCsv(output/"pressure_flow_steps.csv");
	const auto ports = ReadCsv(output/"pressure_flow_ports.csv");
	const auto history = ReadCsv(output/"zero_d_flow_history.csv");
	if (steps.size() != kZeroDClockSteps || history.size() != 2*kZeroDClockSteps)
		throw std::runtime_error("0D clock fixture long-form output is incomplete");
	int nullable_area_rows = 0;
	for (const auto& row : ports)
		if (row.at("domain_id") == "source_0d" || row.at("domain_id") == "terminal_0d") {
			if (!row.at("area_m2").empty())
				throw std::runtime_error("0D port output fabricated an area");
			++nullable_area_rows;
		}
	if (nullable_area_rows != 2*kZeroDClockSteps)
		throw std::runtime_error("0D port output omitted nullable-area rows");
	std::map<std::string, int> history_domains;
	for (const auto& row : history) {
		++history_domains[row.at("domain_id")];
		for (const auto& field : {"stored_pressure_pa", "initial_stored_volume_m3",
			"final_stored_volume_m3", "prescribed_source_amount_m3",
			"distal_sink_amount_m3", "outward_graph_port_amount_m3", "residual_m3"})
			if (!std::isfinite(Value(row, field)))
				throw std::runtime_error("0D history contains a non-finite value");
		const double scale = std::max({std::abs(Value(row, "final_stored_volume_m3")
			-Value(row, "initial_stored_volume_m3")),
			std::abs(Value(row, "prescribed_source_amount_m3")),
			std::abs(Value(row, "distal_sink_amount_m3")),
			std::abs(Value(row, "outward_graph_port_amount_m3"))});
		if (std::abs(Value(row, "residual_m3")) > 1.0e-12*std::max(1.0, scale))
			throw std::runtime_error("0D step balance is not conservative");
	}
	if (history_domains != std::map<std::string, int>{{"source_0d", kZeroDClockSteps},
		{"terminal_0d", kZeroDClockSteps}})
		throw std::runtime_error("0D history does not cover both model roles");
}

} // namespace

int main(int argc, char** argv)
{
	try {
		if (argc < 1) throw std::runtime_error("missing test executable path");
		executable_directory = fs::absolute(argv[0]).parent_path();
		const auto root = fs::temp_directory_path()
			/("tubularflowiga_multidomain_"+std::to_string(static_cast<long long>(getpid())));
		for (const auto& directory : {"three_a", "three_b", "source", "bridge",
			"leaf_a", "leaf_b", "leaf_c", "junction", "species_source",
			"species_leaf_a", "species_leaf_b", "immersed", "immersed_source",
			"immersed_sink", "zero_source", "zero_terminal", "zero_junction"})
			fs::create_directories(root/directory);
		// The eight 0.1 s steps cross the n*dt rounding divergence point.  The
		// body-fitted single-element case also exercises source-0D initialization
		// without turning this into a benchmark solve.
		WriteZeroDModel(root/"zero_source", true);
		WriteZeroDModel(root/"zero_terminal", false);
		WriteZeroDThreeDCase(root/"zero_junction");
		WriteDatabase(root/"zero_junction.ntiga", 1);
		WriteZeroDClockGraph(root);
		// The runner must reject a loaded model-role/port mismatch before it can
		// create a native 3D runtime or output directory.
		WriteZeroDModel(root/"zero_source", false);
		if (RunMultidomain(root, root/"zero_d_invalid_role") == 0
			|| fs::exists(root/"zero_d_invalid_role/graph_binding_manifest.json"))
			throw std::runtime_error("invalid schema-v5 0D model-role topology was accepted");
		RequireLogContains(root/"zero_d_invalid_role.log",
			"0D flow domain port does not match its model role");
		WriteZeroDModel(root/"zero_source", true);
		if (RunMultidomain(root, root/"zero_d_clock_one") != 0)
			throw std::runtime_error("schema-v5 0D exact-clock production run failed");
		ValidateZeroDClockRun(root/"zero_d_clock_one");
		if (RunMultidomain(root, root/"zero_d_clock_two") != 0)
			throw std::runtime_error("repeated schema-v5 0D production run failed");
		ValidateZeroDClockRun(root/"zero_d_clock_two");
		if (ReadFile(root/"zero_d_clock_one/zero_d_flow_history.csv")
			!= ReadFile(root/"zero_d_clock_two/zero_d_flow_history.csv"))
			throw std::runtime_error("schema-v5 0D history is not deterministic");
		if (std::getenv("TUBULARFLOWIGA_ZERO_D_CLOCK_FIXTURE_ONLY")) {
			fs::remove_all(root);
			std::cout << "schema-v5 0D exact-clock production fixture passed\n";
			return 0;
		}
		WriteOneDCase(root/"immersed_source", 1.0e-3);
		WriteOneDCase(root/"immersed_sink", 1.0e-3);
		WriteImmersedCase(root/"immersed");
		WriteImmersedGraph(root);
		if (RunMultidomain(root, root/"immersed_one") != 0)
			throw std::runtime_error("one-rank generic immersed multidomain run failed");
		ValidateImmersedManifest(root/"immersed_one");
		WriteThreeDCase(root/"three_a");
		WriteThreeDCase(root/"three_b");
		WriteOneDCase(root/"source", 1.0e-3);
		WriteOneDCase(root/"bridge", 5.5e-4);
		WriteOneDCase(root/"leaf_a", 4.5e-4);
		WriteOneDCase(root/"leaf_b", 3.0e-4);
		WriteOneDCase(root/"leaf_c", 2.5e-4);
		WriteDatabase(root/"a.ntiga", 1);
		WriteDatabase(root/"b.ntiga", 1);
		WriteMultiIslandGraph(root, "explicit");
		if (RunBifurcation(root, root/"bifurcation_rejected") == 0
			|| fs::exists(root/"bifurcation_rejected/graph_binding_manifest.json"))
			throw std::runtime_error("legacy bifurcation entry accepted a multi-island graph");
		if (RunMultidomain(root, root/"explicit_one") != 0)
			throw std::runtime_error("one-rank explicit multi-island run failed");
		ValidateMultidomain(root/"explicit_one", "explicit");
		WriteMultiIslandGraph(root, "explicit", true);
		if (RunMultidomain(root, root/"explicit_permuted") != 0)
			throw std::runtime_error("permuted explicit multi-island run failed");
		ValidateMultidomain(root/"explicit_permuted", "explicit");
		const auto canonical_edges = ReadCsv(root/"explicit_one/pressure_flow_edges.csv");
		const auto permuted_edges = ReadCsv(root/"explicit_permuted/pressure_flow_edges.csv");
		if (canonical_edges.size() != permuted_edges.size())
			throw std::runtime_error("permuted multi-island row count changed");
		for (std::size_t row = 0; row < canonical_edges.size(); ++row)
			if (canonical_edges[row].at("edge_id") != permuted_edges[row].at("edge_id")
				|| std::abs(Value(canonical_edges[row], "measured_pressure_pa")
					-Value(permuted_edges[row], "measured_pressure_pa")) > 1.0e-12)
				throw std::runtime_error("permuted multi-island execution changed");
		WriteMultiIslandGraph(root, "explicit");
		if (RunMultidomain(root, root/"failed", {},
			"TUBULARFLOWIGA_INJECT_BIFURCATION_FAILURE_STEP=1 ") == 0
			|| fs::exists(root/"failed/graph_binding_manifest.json"))
			throw std::runtime_error("failed multi-island run published a completion marker");
		WriteMultiIslandGraph(root, "aitken");
		if (RunMultidomain(root, root/"aitken_one") != 0)
			throw std::runtime_error("one-rank Aitken multi-island run failed");
		ValidateMultidomain(root/"aitken_one", "aitken");
		WriteDatabase(root/"a_two.ntiga", 2);
		WriteDatabase(root/"b_two.ntiga", 2);
		WriteMultiIslandGraph(root, "explicit");
		{
			std::string config = ReadFile(root/"simulation_config.json");
			const auto first = config.find("a.ntiga");
			const auto second = config.find("b.ntiga");
			if (first == std::string::npos || second == std::string::npos)
				throw std::runtime_error("multi-island database paths are missing");
			config.replace(first, 7, "a_two.ntiga");
			config.replace(config.find("b.ntiga"), 7, "b_two.ntiga");
			std::ofstream rewritten(root/"simulation_config.json", std::ios::trunc);
			rewritten << config;
		}
		if (RunMultidomain(root, root/"explicit_two", "mpiexec -np 2 ") != 0)
			throw std::runtime_error("two-rank explicit multi-island run failed");
		ValidateMultidomain(root/"explicit_two", "explicit");
		const auto one = ReadCsv(root/"explicit_one/pressure_flow_edges.csv");
		const auto two = ReadCsv(root/"explicit_two/pressure_flow_edges.csv");
		if (one.size() != two.size()) throw std::runtime_error("multi-island MPI row mismatch");
		for (std::size_t row = 0; row < one.size(); ++row) {
			if (one[row].at("edge_id") != two[row].at("edge_id"))
				throw std::runtime_error("multi-island MPI ordering mismatch");
			for (const auto& name : {"measured_pressure_pa", "first_outward_flow_m3_s",
				"second_outward_flow_m3_s"})
				if (std::abs(Value(one[row], name)-Value(two[row], name))
					> 1.0e-11*std::max({1.0, std::abs(Value(one[row], name)),
						std::abs(Value(two[row], name))}))
					throw std::runtime_error("multi-island MPI numerical mismatch");
		}
		WriteOneDTransportCase(root/"species_source", 1.0e-3,
			"source_red", "source_blue");
		WriteOneDTransportCase(root/"species_leaf_a", 6.0e-4,
			"leaf_a_red", "leaf_a_blue");
		WriteOneDTransportCase(root/"species_leaf_b", 4.0e-4,
			"leaf_b_red", "leaf_b_blue");
		WriteThreeDTransportCase(root/"junction");
		WriteDatabase(root/"junction.ntiga", 1);
		WriteSpeciesGraph(root, "junction.ntiga");
		if (RunBifurcation(root, root/"species_bifurcation_rejected") == 0
			|| fs::exists(root/"species_bifurcation_rejected/graph_binding_manifest.json"))
			throw std::runtime_error("legacy bifurcation entry accepted schema-v6 graph");
		if (RunMultidomain(root, root/"species_one", {}, {}, true) != 0)
			throw std::runtime_error("one-rank schema-v6 multidomain run failed");
		ValidateSpeciesMultidomain(root/"species_one");
		WriteSpeciesGraph(root, "junction.ntiga", true);
		if (RunMultidomain(root, root/"species_permuted", {}, {}, true) != 0)
			throw std::runtime_error("permuted schema-v6 multidomain run failed");
		ValidateSpeciesMultidomain(root/"species_permuted");
		const auto species_one = ReadCsv(root/"species_one/species_edge_amounts.csv");
		const auto species_permuted = ReadCsv(root/"species_permuted/species_edge_amounts.csv");
		if (species_one.size() != species_permuted.size())
			throw std::runtime_error("permuted schema-v6 edge/species row count changed");
		for (std::size_t row = 0; row < species_one.size(); ++row) {
			if (species_one[row].at("edge_id") != species_permuted[row].at("edge_id")
				|| species_one[row].at("species_id") != species_permuted[row].at("species_id")
				|| species_one[row].at("first_domain_id")
					!= species_permuted[row].at("first_domain_id")
				|| species_one[row].at("first_port_id")
					!= species_permuted[row].at("first_port_id")
				|| species_one[row].at("second_domain_id")
					!= species_permuted[row].at("second_domain_id")
				|| species_one[row].at("second_port_id")
					!= species_permuted[row].at("second_port_id")
				|| species_one[row].at("donor") != species_permuted[row].at("donor"))
				throw std::runtime_error("permuted schema-v6 edge/species ordering changed");
			for (const auto& name : {"first_outward_amount", "second_outward_amount",
				"amount_residual"})
				if (std::abs(Value(species_one[row], name)-Value(species_permuted[row], name))
					> 1.0e-10)
					throw std::runtime_error("permuted schema-v6 amounts changed");
		}
		WriteDatabase(root/"junction_two.ntiga", 2);
		WriteSpeciesGraph(root, "junction_two.ntiga");
		if (RunMultidomain(root, root/"species_failed", "mpiexec -np 2 ",
			"TUBULARFLOWIGA_INJECT_BIFURCATION_FAILURE_STEP=1 ", true) == 0
			|| fs::exists(root/"species_failed")
			|| fs::exists(root/"species_failed/graph_binding_manifest.json"))
			throw std::runtime_error("schema-v6 failure published accepted output");
		RequireLogContains(root/"species_failed.log",
			"injected bifurcation failure before commit");
		{
			std::string invalid = ReadFile(root/"simulation_config.json");
			const auto binding = invalid.find("\"source_red\"");
			if (binding == std::string::npos)
				throw std::runtime_error("schema-v6 source binding is missing from fixture");
			invalid.replace(binding, std::string("\"source_red\"").size(),
				"\"missing_native\"");
			std::ofstream rewritten(root/"simulation_config.json", std::ios::trunc);
			rewritten << invalid;
		}
		if (RunMultidomain(root, root/"species_wrong_binding", "mpiexec -np 2 ", {}, true) == 0
			|| fs::exists(root/"species_wrong_binding/graph_binding_manifest.json"))
			throw std::runtime_error("schema-v6 wrong species binding was accepted");
		RequireLogContains(root/"species_wrong_binding.log",
			"schema-v6 1D species binding names no transported native field");
		WriteOneDTransportCase(root/"species_source", 1.0e-3,
			"source_red", "source_blue");
		{
			std::string invalid = ReadFile(root/"species_source/simulation_config.json");
			const auto mode = invalid.find("\"mode\":\"flow_only\"");
			if (mode == std::string::npos)
				throw std::runtime_error("schema-v6 1D coupling fixture cannot be rewritten");
			invalid.replace(mode, std::string("\"mode\":\"flow_only\"").size(),
				"\"mode\":\"vca_replay\"");
			const auto scheme = invalid.find("\"scheme\":\"explicit_staggered\",");
			if (scheme == std::string::npos)
				throw std::runtime_error("schema-v6 1D coupling fixture cannot be rewritten");
			invalid.replace(scheme, std::string("\"scheme\":\"explicit_staggered\",").size(),
				"\"scheme\":\"explicit_staggered\",\"replay_file\":\"replay.json\",");
			std::ofstream rewritten(root/"species_source/simulation_config.json",
				std::ios::trunc);
			rewritten << invalid;
		}
		WriteSpeciesGraph(root, "junction_two.ntiga");
		if (RunMultidomain(root, root/"species_nonflow_one_d", "mpiexec -np 2 ", {}, true)
			== 0 || fs::exists(root/"species_nonflow_one_d/graph_binding_manifest.json"))
			throw std::runtime_error("schema-v6 non-flow-only 1D coupling was accepted");
		RequireLogContains(root/"species_nonflow_one_d.log",
			"schema-v6 multidomain transport requires one 1D flow system, one transport system, flow_only coupling");
		WriteOneDTransportCase(root/"species_source", 1.0e-3,
			"source_red", "source_blue");
		WriteSpeciesGraph(root, "junction_two.ntiga");
		WriteThreeDTransportCase(root/"junction", false);
		if (RunMultidomain(root, root/"species_bad_velocity", "mpiexec -np 2 ", {}, true) == 0
			|| fs::exists(root/"species_bad_velocity/graph_binding_manifest.json"))
			throw std::runtime_error("schema-v6 non-prescribed transport velocity was accepted");
		RequireLogContains(root/"species_bad_velocity.log",
			"schema-v6 multidomain transport requires prescribed 3D transport velocity");
		WriteThreeDTransportCase(root/"junction");
		WriteSpeciesGraph(root, "junction_two.ntiga");
		if (RunMultidomain(root, root/"species_two", "mpiexec -np 2 ", {}, true) != 0)
			throw std::runtime_error("two-rank schema-v6 multidomain run failed");
		ValidateSpeciesMultidomain(root/"species_two");
		const auto species_two = ReadCsv(root/"species_two/species_edge_amounts.csv");
		if (species_one.size() != species_two.size())
			throw std::runtime_error("schema-v6 MPI edge/species row count changed");
		for (std::size_t row = 0; row < species_one.size(); ++row) {
			if (species_one[row].at("edge_id") != species_two[row].at("edge_id")
				|| species_one[row].at("species_id") != species_two[row].at("species_id"))
				throw std::runtime_error("schema-v6 MPI edge/species ordering changed");
			for (const auto& name : {"first_outward_amount", "second_outward_amount",
				"amount_residual"})
				if (std::abs(Value(species_one[row], name)-Value(species_two[row], name))
					> 1.0e-9)
					throw std::runtime_error("schema-v6 MPI amounts changed");
		}
		fs::remove_all(root);
		std::cout << std::setprecision(17)
			<< "multidomain flow smoke test passed species_max_edge_residual="
			<< species_residual_maxima.edge << " species_max_domain_residual="
			<< species_residual_maxima.domain << " species_max_global_residual="
			<< species_residual_maxima.global << '\n';
		return 0;
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n';
		return 1;
	}
}
