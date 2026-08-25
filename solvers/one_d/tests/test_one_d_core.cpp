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
    {"name": "inlet_q", "kind": "constant", "units": "m3/s", "value": 1.0e-9}
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
	iga::AdvanceOneDTransport(configuration, network, flow_state, transport, ".", 0.0,
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

	std::cout << "one-dimensional core tests passed\n";
}
