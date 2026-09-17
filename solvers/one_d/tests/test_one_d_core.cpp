#include "OneDTransport.hpp"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace fs = std::filesystem;

namespace {

std::string Configuration(const std::string& model = "rigid",
	const std::string& scheme = "steady_poiseuille")
{
	return std::string(R"json({
  "schema_version": 3,
  "dimension": "1d",
  "geometry": {"kind": "swc_network", "file": "tree.swc", "length_scale_to_m": 1.0},
  "fields": [
    {"name": "area", "kind": "scalar"},
    {"name": "flow_rate", "kind": "scalar"},
    {"name": "pressure", "kind": "pressure"},
    {"name": "oxygen", "kind": "scalar", "initial_value": 0.2}
  ],
  "time": {"dt": 0.00001, "steps": 2, "output_every": 1},
  "temporal_functions": [
    {"name": "inlet_q", "kind": "constant", "units": "m3/s", "value": 1.0e-9},
    {"name": "pulse", "kind": "sinusoid", "units": "m3/s",
     "mean": 2.0, "amplitude": 3.0, "period": 1.0, "phase": 0.0},
    {"name": "measured", "kind": "periodic_table", "units": "m3/s",
     "file": "inlet_flow.csv", "period": 1.0, "interpolation": "linear"},
    {"name": "harmonics", "kind": "fourier", "units": "m3/s",
     "mean": 1.0, "period": 1.0, "cosine": [2.0], "sine": [3.0]}
  ],
  "equation_systems": [
    {"name": "flow", "kind": "network_flow_1d", "unknowns": ["area", "flow_rate", "pressure"],
     "model": ")json") + model + R"json(", "scheme": ")json" + scheme + R"json(",
     "dynamic_viscosity": 0.004, "density": 1060.0,
     "wall": {"model": "linear", "young_modulus": 1000000.0, "thickness_ratio": 0.1},
     "discretization": {"cells_per_segment": 2, "cfl": 0.75, "min_area_fraction": 0.2}},
    {"name": "transport", "kind": "network_transport_1d", "unknowns": ["oxygen"],
     "flow_system": "flow", "species": [
       {"field": "oxygen", "diffusivity": 1.0e-9, "reaction_rate": 0.0, "volume_source": 0.0}
     ]}
  ],
  "boundaries": [
    {"name": "inlet", "role": "inlet", "node_ids": [1], "conditions": [
      {"field": "flow_rate", "type": "dirichlet", "quantity": "flow_rate", "waveform": "inlet_q"},
      {"field": "oxygen", "type": "dirichlet", "value": 1.0}
    ]},
    {"name": "outlets", "role": "outlet", "conditions": [
      {"field": "pressure", "type": "pressure", "value": 0.0}
    ]},
    {"name": "wall", "role": "wall", "conditions": [
      {"field": "oxygen", "type": "no_flux"}
    ]}
  ]
})json";
}

std::string ZeroDConfiguration()
{
	return R"json({
  "schema_version": 3,
  "dimension": "0d",
  "geometry": {"kind": "swc_network", "file": "tree.swc", "length_scale_to_m": 1.0},
  "fields": [
    {"name": "area", "kind": "scalar"},
    {"name": "flow_rate", "kind": "scalar"},
    {"name": "pressure", "kind": "pressure"}
  ],
  "time": {"dt": 0.001, "steps": 10, "output_every": 1},
  "temporal_functions": [
    {"name": "inlet_q", "kind": "constant", "units": "m3/s", "value": 1.0e-9}
  ],
  "equation_systems": [{
    "name": "circuit", "kind": "network_flow_0d",
    "unknowns": ["area", "flow_rate", "pressure"],
    "model": "lumped", "scheme": "petsc", "formulation": "transient_rc",
    "dynamic_viscosity": 0.004, "density": 1060.0,
    "lumped_parameters": {
      "resistance_scale": 1.0, "compliance_scale": 0.0,
      "segment_resistance": {"2": 100000000.0},
      "segment_compliance": {"2": 1.0e-10}
    }
  }],
  "boundaries": [
    {"name": "inlet", "role": "inlet", "node_ids": [1], "conditions": [
      {"field": "flow_rate", "type": "dirichlet", "quantity": "flow_rate", "waveform": "inlet_q"}
    ]},
    {"name": "outlet", "role": "outlet", "conditions": [
      {"field": "pressure", "type": "windkessel_rc", "resistance": 1.0e9,
       "capacitance": 1.0e-10, "reference_pressure": 0.0, "initial_pressure": 0.0}
    ]}
  ]
})json";
}

void WriteTree(const fs::path& path)
{
	std::ofstream output(path);
	output << "1 2 0 0 0 0.001 -1\n"
		<< "2 2 0.01 0 0 0.001 1\n"
		<< "3 2 0.02 0.01 0 0.001 2\n"
		<< "4 2 0.02 -0.01 0 0.001 2\n";
}

} // namespace

int main()
{
	const auto configuration = iga::ParseOneDConfiguration(Configuration());
	assert(configuration.schema_version == 3);
	assert(configuration.flow_systems.size() == 1);
	assert(configuration.transport_systems.size() == 1);
	assert(configuration.temporal_functions.size() == 4);
	assert(!configuration.warnings.empty());
	assert(std::string(iga::NetworkFlowPhysicalDimension(
		configuration.flow_systems.front())) == "0d");
	const auto zero_d_configuration = iga::ParseOneDConfiguration(ZeroDConfiguration());
	assert(zero_d_configuration.dimension == "0d");
	assert(zero_d_configuration.flow_systems.front().model == iga::OneDFlowModel::Lumped);
	assert(zero_d_configuration.flow_systems.front().scheme == iga::OneDFlowScheme::Petsc);
	assert(zero_d_configuration.flow_systems.front().formulation
		== iga::OneDImplicitFormulation::TransientRc);
	assert(zero_d_configuration.flow_systems.front().lumped.compliance_scale == 0.0);
	assert(zero_d_configuration.boundaries.back().conditions.front().type == "windkessel_rc");
	assert(zero_d_configuration.warnings.empty());
	auto rlc_text = ZeroDConfiguration();
	const auto rc_name = rlc_text.find("\"transient_rc\"");
	assert(rc_name != std::string::npos);
	rlc_text.replace(rc_name, std::string("\"transient_rc\"").size(), "\"transient_rlc\"");
	const auto rlc_configuration = iga::ParseOneDConfiguration(rlc_text);
	assert(rlc_configuration.flow_systems.front().formulation
		== iga::OneDImplicitFormulation::TransientRlc);
	auto legacy_zero_d = ZeroDConfiguration();
	const auto petsc_name = legacy_zero_d.find("\"petsc\"");
	legacy_zero_d.replace(petsc_name, std::string("\"petsc\"").size(), "\"implicit_petsc\"");
	const auto transient_name = legacy_zero_d.find("\"transient_rc\"");
	legacy_zero_d.replace(transient_name, std::string("\"transient_rc\"").size(), "\"lumped_rc\"");
	const auto legacy_configuration = iga::ParseOneDConfiguration(legacy_zero_d);
	assert(legacy_configuration.flow_systems.front().scheme == iga::OneDFlowScheme::Petsc);
	assert(legacy_configuration.flow_systems.front().formulation
		== iga::OneDImplicitFormulation::TransientRc);
	assert(legacy_configuration.warnings.size() == 2);
	const auto waveform_directory = fs::temp_directory_path()/"tubularflowiga-one-d-waveforms";
	fs::create_directories(waveform_directory);
	{
		std::ofstream output(waveform_directory/"inlet_flow.csv");
		output << "time,value\n0,1\n0.5,3\n";
	}
	const iga::OneDWaveformEvaluator waveforms(configuration, waveform_directory);
	assert(std::abs(waveforms.Evaluate("inlet_q", 4.0)-1.0e-9) < 1.0e-20);
	assert(std::abs(waveforms.Evaluate("pulse", 0.25)-5.0) < 1.0e-14);
	assert(std::abs(waveforms.Evaluate("harmonics", 0.25)-4.0) < 1.0e-14);
	assert(std::abs(waveforms.Evaluate("measured", 0.25)-2.0) < 1.0e-14);
	fs::remove(waveform_directory/"inlet_flow.csv");
	assert(std::abs(waveforms.Evaluate("measured", 1.25)-2.0) < 1.0e-14);
	iga::OneDInletState velocity_inlet;
	velocity_inlet.quantity = "centerline_velocity";
	velocity_inlet.waveform = "pulse";
	assert(std::abs(iga::EvaluateOneDInlet(velocity_inlet, waveforms, 0.25, 0.2)-0.5) < 1.0e-14);
	auto invalid_waveform_configuration = configuration;
	for (auto& function : invalid_waveform_configuration.temporal_functions)
		if (function.name == "measured") function.file = "/tmp/outside.csv";
	bool waveform_path_rejected = false;
	try {
		const iga::OneDWaveformEvaluator invalid_waveforms(
			invalid_waveform_configuration, waveform_directory);
	}
	catch (const std::runtime_error&) { waveform_path_rejected = true; }
	assert(waveform_path_rejected);
	fs::remove(waveform_directory);
	auto obj_configuration_text = Configuration();
	const auto kind_position = obj_configuration_text.find("\"kind\": \"swc_network\"");
	const auto file_position = obj_configuration_text.find("\"file\": \"tree.swc\"");
	assert(kind_position != std::string::npos && file_position != std::string::npos);
	obj_configuration_text.replace(kind_position, std::string("\"kind\": \"swc_network\"").size(),
		"\"kind\": \"obj_network\"");
	const auto shifted_file_position = obj_configuration_text.find("\"file\": \"tree.swc\"");
	obj_configuration_text.replace(shifted_file_position, std::string("\"file\": \"tree.swc\"").size(),
		"\"file\": \"tree.obj\", \"root_node_id\": 1");
	const auto obj_configuration = iga::ParseOneDConfiguration(obj_configuration_text);
	assert(obj_configuration.geometry.kind == "obj_network");
	assert(obj_configuration.geometry.root_node_id == 1);

	const auto temporary = fs::temp_directory_path()/"tubularflowiga-one-d-core-test.swc";
	WriteTree(temporary);
	auto network = iga::ReadOneDNetwork(temporary, 1.0, 2, 0.004);
	fs::remove(temporary);
	assert(network.nodes.size() == 4);
	assert(network.segments.size() == 3);
	assert(network.outlet_nodes.size() == 2);
	assert(network.cells == 6);
	iga::ValidateOneDTopologyReferences(configuration, network);

	const auto temporary_obj = fs::temp_directory_path()/"tubularflowiga-one-d-core-test.obj";
	{
		std::ofstream output(temporary_obj);
		output << "v 0 0 0 0.001 0 0\n"
			<< "v 0.01 0 0 0.001 0 0\n"
			<< "v 0.02 0.01 0 0.001 0 0\n"
			<< "v 0.02 -0.01 0 0.001 0 0\n"
			<< "l 1 2\n"
			<< "l 2 3\n"
			<< "l 2 4\n";
	}
	const auto obj_network = iga::ReadOneDNetwork(temporary_obj, 1.0, 2, 0.004, 1);
	fs::remove(temporary_obj);
	assert(obj_network.root == 0);
	assert(obj_network.nodes.size() == 4);
	assert(obj_network.segments.size() == 3);
	assert(obj_network.outlet_nodes.size() == 2);
	const auto skeleton_output = fs::temp_directory_path()/"tubularflowiga-one-d-skeleton-output";
	fs::create_directories(skeleton_output);
	iga::WriteOneDSkeletonFiles(skeleton_output, obj_network, 1.0);
	assert(fs::is_regular_file(skeleton_output/"skeleton_normalized.swc"));
	assert(fs::is_regular_file(skeleton_output/"skeleton.vtp"));
	const auto normalized_network = iga::ReadOneDNetwork(
		skeleton_output/"skeleton_normalized.swc", 1.0, 2, 0.004);
	assert(normalized_network.nodes.size() == obj_network.nodes.size());
	{
		std::ifstream input(skeleton_output/"skeleton.vtp");
		std::ostringstream contents;
		contents << input.rdbuf();
		assert(contents.str().find("Name=\"role\"") != std::string::npos);
		assert(contents.str().find("Name=\"segment_id\"") != std::string::npos);
	}
	fs::remove(skeleton_output/"skeleton_normalized.swc");
	fs::remove(skeleton_output/"skeleton.vtp");
	fs::remove(skeleton_output);

	const auto& segment = network.segments.front();
	const double expected_resistance = 8.0*0.004*0.01/(iga::OneDPi*std::pow(0.001, 4.0));
	assert(std::abs(segment.resistance/expected_resistance-1.0) < 1.0e-14);

	iga::OneDFlowState flow_state;
	flow_state.outlets = iga::ResolveOneDOutlets(configuration, network);
	iga::SolveRigidOneD(network, configuration.flow_systems.front(), flow_state, 1.0e-9,
		configuration.time.dt);
	assert(flow_state.segment_flow.size() == 3);
	assert(std::abs(flow_state.segment_flow[1]-0.5e-9) < 1.0e-20);
	assert(std::abs(flow_state.segment_flow[2]-0.5e-9) < 1.0e-20);
	assert(flow_state.node_pressure[0] > flow_state.node_pressure[1]);
	iga::OneDFlowState compliant_initial_state;
	compliant_initial_state.outlets = iga::ResolveOneDOutlets(configuration, network);
	iga::InitializeCompliantOneDFromRigid(network, configuration.flow_systems.front(),
		compliant_initial_state, 1.0e-9, 0.0);
	assert(std::abs(compliant_initial_state.inlet_flow-1.0e-9) < 1.0e-20);
	iga::OneDFlowState unequal_pressure_state;
	unequal_pressure_state.outlets = iga::ResolveOneDOutlets(configuration, network);
	unequal_pressure_state.outlets[0].pressure = 1.0;
	const double branch_resistance = network.segments[1].resistance;
	iga::SolveRigidOneD(network, configuration.flow_systems.front(),
		unequal_pressure_state, 1.0e-9, configuration.time.dt);
	const double expected_first_branch_flow = 0.5e-9-0.5/branch_resistance;
	assert(std::abs(unequal_pressure_state.segment_flow[1]
		-expected_first_branch_flow) < 1.0e-20);
	assert(std::abs(unequal_pressure_state.segment_flow[2]
		-(1.0e-9-expected_first_branch_flow)) < 1.0e-20);
	assert(std::abs(unequal_pressure_state.node_pressure[2]-1.0) < 1.0e-12);
	assert(std::abs(unequal_pressure_state.node_pressure[3]) < 1.0e-12);

	const auto& wall = configuration.flow_systems.front().wall;
	const double expanded = segment.area0*1.1;
	const double pressure = iga::OneDPressureFromArea(expanded, segment.area0, segment.radius0, wall);
	assert(std::abs(iga::OneDAreaFromPressure(pressure, segment.area0, segment.radius0, wall)/expanded-1.0) < 1.0e-13);
	assert(iga::OneDWaveSpeed(segment.area0, segment.area0, segment.radius0, wall, 1060.0) > 0.0);
	iga::OneDWallDefinition olufsen = wall;
	olufsen.model = iga::OneDWallModel::Olufsen;
	const double olufsen_pressure = iga::OneDPressureFromArea(expanded,
		segment.area0, segment.radius0, olufsen);
	assert(std::abs(iga::OneDAreaFromPressure(olufsen_pressure, segment.area0,
		segment.radius0, olufsen)/expanded-1.0) < 1.0e-13);

	iga::OneDOutletState rcr;
	rcr.kind = iga::OneDOutletKind::WindkesselRcr;
	rcr.proximal_resistance = 1.0e8;
	rcr.distal_resistance = 1.0e9;
	rcr.capacitance = 1.0e-10;
	rcr.reference_pressure = 2.0;
	rcr.capacitor_pressure = 3.0;
	const double rcr_dt = 1.0e-3;
	const double rcr_flow = 1.0e-9;
	const double expected_capacitor = (3.0+rcr_dt*(rcr_flow
		+rcr.reference_pressure/rcr.distal_resistance)/rcr.capacitance)
		/(1.0+rcr_dt/(rcr.distal_resistance*rcr.capacitance));
	iga::AdvanceOutletState(rcr, rcr_flow, rcr_dt);
	assert(std::abs(rcr.capacitor_pressure-expected_capacitor) < 1.0e-14);
	assert(std::abs(rcr.pressure-(rcr.proximal_resistance*rcr_flow
		+expected_capacitor)) < 1.0e-14);

	iga::OneDFlowSystemDefinition table_flow = configuration.flow_systems.front();
	table_flow.junctions.loss_model = "table";
	table_flow.junctions.angle_table = {{0.0, 0.0}, {90.0, 0.9}};
	const double table_loss = iga::OneDJunctionLossCoefficient(table_flow, network, 1,
		network.segments[0], network.segments[1], 0.5);
	assert(table_loss > 0.0 && table_loss < 0.9);
	const iga::OneDConservativeState uniform{segment.area0, 1.0e-9};
	const auto physical_flux = iga::OneDFlux(uniform, segment,
		configuration.flow_systems.front());
	const auto numerical_flux = iga::OneDRusanovFlux(uniform, uniform, segment,
		configuration.flow_systems.front());
	assert(std::abs(physical_flux.area-numerical_flux.area) < 1.0e-30);
	assert(std::abs(physical_flux.flow-numerical_flux.flow) < 1.0e-30);

	auto transport = iga::InitializeOneDTransport(configuration,
		configuration.transport_systems.front(), network);
	iga::AdvanceOneDTransport(configuration, network, flow_state, transport, waveforms, 0.0,
		configuration.time.dt);
	assert(transport.species.size() == 1);
	for (const double value : transport.species.front().concentration)
		assert(std::isfinite(value));

	const auto explicit_configuration = iga::ParseOneDConfiguration(
		Configuration("compliant", "explicit_rusanov"));
	iga::OneDFlowState explicit_state;
	explicit_state.outlets = iga::ResolveOneDOutlets(explicit_configuration, network);
	iga::InitializeCompliantOneD(network, explicit_configuration.flow_systems.front(), explicit_state);
	iga::AdvanceExplicitOneD(network, explicit_configuration.flow_systems.front(),
		explicit_state, 1.0e-9, explicit_configuration.time.dt);
	for (std::size_t i = 0; i < explicit_state.area.size(); ++i) {
		assert(explicit_state.area[i] > 0.0);
		assert(std::isfinite(explicit_state.flow[i]));
	}

	bool rejected = false;
	try {
		iga::ParseOneDConfiguration(Configuration("rigid", "explicit_rusanov"));
	} catch (const std::runtime_error&) { rejected = true; }
	assert(rejected);

	const auto invalid = fs::temp_directory_path()/"tubularflowiga-one-d-invalid.swc";
	{
		std::ofstream output(invalid);
		output << "1 2 0 0 0 0.001 -1\n2 2 0 0 0 0.001 1\n";
	}
	rejected = false;
	try { iga::ReadOneDNetwork(invalid, 1.0, 1, 0.004); }
	catch (const std::runtime_error&) { rejected = true; }
	fs::remove(invalid);
	assert(rejected);

	const auto cyclic_obj = fs::temp_directory_path()/"tubularflowiga-one-d-cycle.obj";
	{
		std::ofstream output(cyclic_obj);
		output << "v 0 0 0 0.001 0 0\n"
			<< "v 1 0 0 0.0009 0 0\n"
			<< "v 2 1 0 0.0008 0 0\n"
			<< "v 2 -1 0 0.0008 0 0\n"
			<< "v 3 0 0 0.0007 0 0\n"
			<< "v 4 0 0 0.0006 0 0\n"
			<< "l 1 2\n"
			<< "l 2 3\n"
			<< "l 2 4\n"
			<< "l 3 5\n"
			<< "l 4 5\n"
			<< "l 5 6\n";
	}
	rejected = false;
	try { iga::ReadOneDNetwork(cyclic_obj, 1.0, 1, 0.004, 1); }
	catch (const std::runtime_error&) { rejected = true; }
	assert(rejected);
	const auto cyclic_network = iga::ReadOneDNetwork(
		cyclic_obj, 1.0, 1, 0.004, 1, true);
	fs::remove(cyclic_obj);
	assert(cyclic_network.has_cycles);
	assert(cyclic_network.segments.size() == 6);
	assert(cyclic_network.outlet_nodes.size() == 1);
	assert(cyclic_network.nodes[static_cast<std::size_t>(
		cyclic_network.outlet_nodes.front())].id == 6);

	std::cout << "one-dimensional core tests passed\n";
}
