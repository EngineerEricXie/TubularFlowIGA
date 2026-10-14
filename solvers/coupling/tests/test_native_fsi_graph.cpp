// Graph and configuration contract of the native tetrahedral FSI pair.
#include "MultidomainConfig.hpp"
#include "SimulationGraph.hpp"

#include <cassert>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

const std::string kIdentity(64, 'a');

iga::DistributedSurfaceLayout Layout()
{
	iga::DistributedSurfaceLayout layout;
	layout.reference_mesh_identity_sha256 = kIdentity;
	layout.global_node_count = 3;
	layout.partition_count = 1;
	layout.partition_rank = 0;
	layout.owned_global_node_ids = {0, 1, 2};
	layout.reference_positions = {{0, {{0.0, 0.0, 0.0}}}, {1, {{1.0, 0.0, 0.0}}}, {2, {{0.0, 1.0, 0.0}}}};
	layout.reference_triangles = {{{0, 1, 2}}};
	layout.owned_reference_lumped_areas_m2 = {1.0, 1.0, 1.0};
	layout.layout_identity_sha256 = iga::BuildDistributedSurfaceLayoutIdentitySha256(layout);
	return layout;
}

iga::DistributedSurfaceInterface Surface(const std::string& domain, bool fluid)
{
	iga::DistributedSurfaceInterface surface;
	surface.id = {domain, domain+"_runtime", "wall"};
	surface.subsystem_id = domain+"_runtime";
	surface.boundary_labels = {7};
	surface.reference_mesh_identity_sha256 = kIdentity;
	if (fluid) {
		surface.provides = {iga::SurfaceFieldQuantity::TractionOnStructure};
		surface.requires = {iga::SurfaceFieldQuantity::Displacement, iga::SurfaceFieldQuantity::Velocity};
	} else {
		surface.provides = {iga::SurfaceFieldQuantity::Displacement, iga::SurfaceFieldQuantity::Velocity};
		surface.requires = {iga::SurfaceFieldQuantity::TractionOnStructure};
	}
	return surface;
}

iga::CouplingPort LogicalPort(const std::string& domain_id)
{
	iga::CouplingPort port;
	port.id = "flow_port";
	port.subsystem_id = domain_id;
	port.locator_kind = "runtime_port";
	port.locator = "flow_port";
	port.provides = {iga::PortQuantity::Area, iga::PortQuantity::FlowRate, iga::PortQuantity::MeanPressure};
	port.requires = {iga::PortQuantity::MeanPressure};
	iga::ValidateCouplingPort(port);
	return port;
}

iga::DomainNode Domain(const std::string& id, iga::DomainKind kind, bool fluid)
{
	const auto surface = Surface(id, fluid);
	// Kinds other than the native FSI pair must still carry a scalar port.
	std::vector<iga::CouplingPort> ports;
	if (kind != iga::DomainKind::NativeTetAleFlow && kind != iga::DomainKind::NativeTetSolidStructure)
		ports.push_back(LogicalPort(id));
	return {id, kind, ports, {}, {surface}, {{surface.id, Layout()}}};
}

iga::SimulationGraph Pair(iga::DomainKind fluid_kind, iga::DomainKind structure_kind)
{
	auto fluid = Domain("fluid", fluid_kind, true);
	auto structure = Domain("structure", structure_kind, false);
	const iga::FsiCouplingEdge edge{"wall", fluid.surface_interfaces.front().id,
		structure.surface_interfaces.front().id, iga::FsiCouplingLaw::FluidStructureTractionKinematics};
	return iga::SimulationGraph({std::move(fluid), std::move(structure)}, {}, {}, {edge});
}

void RequireRejectedWith(const std::function<void()>& operation, const std::string& diagnostic)
{
	bool rejected = false;
	try { operation(); }
	catch (const std::runtime_error& error) {
		rejected = std::string(error.what()).find(diagnostic) != std::string::npos;
		if (!rejected) std::cerr << "unexpected diagnostic: " << error.what() << '\n';
	}
	assert(rejected);
}

const char* kConfiguration = R"({
  "schema_version": 5,
  "time": {"dt": 0.01, "steps": 2},
  "start_domain": "fluid",
  "execution": {"kind": "explicit"},
  "domains": [
    {"id": "fluid", "dimension": "3d", "kind": "native_tet_ale_flow", "case": "fluid", "ports": [],
     "surface_interfaces": [{"id": "wall", "boundary_labels": [7]}]},
    {"id": "solid", "dimension": "3d", "kind": "native_tet_solid_structure", "case": "solid", "ports": [],
     "surface_interfaces": [{"id": "wall", "boundary_labels": [7]}]}
  ],
  "couplings": [],
  "fsi_couplings": [{"id": "wall", "fluid": {"domain": "fluid", "interface": "wall"},
                     "structure": {"domain": "solid", "interface": "wall"}}]
})";

std::string Replace(std::string text, const std::string& from, const std::string& to)
{
	const auto position = text.find(from);
	assert(position != std::string::npos);
	return text.replace(position, from.size(), to);
}

} // namespace

int main()
{
	// Graph: only the native ALE fluid / tetrahedral solid pair is an FSI pair.
	{
		const auto graph = Pair(iga::DomainKind::NativeTetAleFlow, iga::DomainKind::NativeTetSolidStructure);
		assert(graph.FsiEdges().size() == 1);
		assert(graph.Domain("fluid").kind == iga::DomainKind::NativeTetAleFlow);
		assert(iga::IsFlowDomainKind(iga::DomainKind::NativeTetAleFlow));
		assert(iga::IsStructureDomainKind(iga::DomainKind::NativeTetSolidStructure));
		assert(iga::DomainTopologyDimensionOf(iga::DomainKind::NativeTetSolidStructure) == 3);
		assert(std::string(iga::DomainKindName(iga::DomainKind::NativeTetAleFlow)) == "native_tet_ale_flow");
	}
	RequireRejectedWith([] { Pair(iga::DomainKind::ThreeDBodyFittedFlow, iga::DomainKind::NativeTetSolidStructure); },
		"unsupported fluid/structure domain pair");
	RequireRejectedWith([] { Pair(iga::DomainKind::NativeTetAleFlow, iga::DomainKind::NativeTetAleFlow); },
		"unsupported fluid/structure domain pair");
	RequireRejectedWith([] {
		auto structure = Domain("structure", iga::DomainKind::NativeTetSolidStructure, false);
		structure.ports = {iga::CouplingPort{}};
		iga::SimulationGraph({std::move(structure)}, {}, {}, {});
	}, "structure domain cannot declare scalar ports");

	// Configuration: the JSON form of the same pair.
	{
		const auto configuration = iga::ParseMultidomainConfiguration(kConfiguration);
		assert(configuration.fsi_couplings.size() == 1);
		assert(configuration.fsi_couplings.front().fluid_domain_id == "fluid");
		assert(configuration.fsi_couplings.front().structure_interface_id == "wall");
		assert(configuration.domains.size() == 2);
		assert(configuration.domains[1].kind == iga::DomainKind::NativeTetSolidStructure);
		assert(configuration.domains[0].surface_interfaces.size() == 1);
		assert(configuration.domains[0].surface_interfaces.front().second.count(7));
		assert(configuration.graph.Edges().empty());
	}
	RequireRejectedWith([] { iga::ParseMultidomainConfiguration(Replace(kConfiguration,
		"\"native_tet_solid_structure\"", "\"native_tet_ale_flow\"")); }, "native_tet_ale_flow fluid to a");
	RequireRejectedWith([] { iga::ParseMultidomainConfiguration(Replace(kConfiguration,
		"\"fsi_couplings\"", "\"ignored\"")); }, "unknown");
	RequireRejectedWith([] { iga::ParseMultidomainConfiguration(Replace(kConfiguration,
		"\"interface\": \"wall\"},", "\"interface\": \"lid\"},")); },
		"undeclared surface interface");
	RequireRejectedWith([] { iga::ParseMultidomainConfiguration(Replace(kConfiguration,
		"\"case\": \"solid\", \"ports\": []", "\"case\": \"solid\", \"ports\": [{\"id\": \"p\", \"locator_kind\": \"runtime_port\","
		" \"locator\": \"p\", \"provides\": [\"area\", \"flow_rate\", \"mean_pressure\"], \"requires\": [\"mean_pressure\"]}]")); },
		"do not accept scalar ports");
	RequireRejectedWith([] { iga::ParseMultidomainConfiguration(Replace(kConfiguration,
		"\"kind\": \"native_tet_ale_flow\"", "\"kind\": \"body_fitted_iga_flow\", \"database\": \"x.ntiga\"")); },
		"surface_interfaces are only valid for native FSI domains");
	std::cout << "native_fsi_graph_test: PASS\n";
	return 0;
}
