#include "SimulationGraph.hpp"
#include "ZeroDFlowDomain.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <functional>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

template <class Function>
void RequireRejected(Function&& function)
{
	bool rejected = false;
	try { function(); }
	catch (const std::runtime_error&) { rejected = true; }
	assert(rejected);
}

iga::CouplingPort LogicalPort(const std::string& domain_id, const std::string& id,
	const std::string& locator, iga::PortQuantity accepted)
{
	iga::CouplingPort port;
	port.id = id;
	port.subsystem_id = domain_id;
	port.locator_kind = "runtime_port";
	port.locator = locator;
	port.provides = {iga::PortQuantity::Area, iga::PortQuantity::FlowRate,
		iga::PortQuantity::MeanPressure};
	port.requires = {accepted};
	iga::ValidateCouplingPort(port);
	return port;
}

iga::DomainNode Domain(const std::string& id, iga::DomainKind kind,
	std::vector<iga::CouplingPort> ports)
{
	return {id, kind, std::move(ports)};
}

std::string SurfaceHash(char character)
{
	return std::string(64, character);
}

iga::DistributedSurfaceLayout SurfaceLayout(const std::string& reference_mesh_identity,
	double area = 1.0)
{
	iga::DistributedSurfaceLayout layout;
	layout.reference_mesh_identity_sha256 = reference_mesh_identity;
	layout.global_node_count = 3;
	layout.partition_count = 1;
	layout.partition_rank = 0;
	layout.owned_global_node_ids = {0, 1, 2};
	layout.reference_positions = {{0, {{0.0, 0.0, 0.0}}},
		{1, {{1.0, 0.0, 0.0}}}, {2, {{0.0, 1.0, 0.0}}}};
	layout.reference_triangles = {{{0, 1, 2}}};
	layout.owned_reference_lumped_areas_m2 = {area, area, area};
	layout.layout_identity_sha256 = iga::BuildDistributedSurfaceLayoutIdentitySha256(layout);
	return layout;
}

iga::DistributedSurfaceInterface FsiSurface(const iga::SurfaceInterfaceRef& id,
	const std::string& reference_mesh_identity, bool fluid)
{
	iga::DistributedSurfaceInterface surface;
	surface.id = id;
	surface.subsystem_id = id.subsystem_id;
	surface.boundary_labels = {0, 7};
	surface.reference_mesh_identity_sha256 = reference_mesh_identity;
	if (fluid) {
		surface.provides = {iga::SurfaceFieldQuantity::TractionOnStructure};
		surface.requires = {iga::SurfaceFieldQuantity::Displacement,
			iga::SurfaceFieldQuantity::Velocity};
	} else {
		surface.provides = {iga::SurfaceFieldQuantity::Displacement,
			iga::SurfaceFieldQuantity::Velocity};
		surface.requires = {iga::SurfaceFieldQuantity::TractionOnStructure};
	}
	iga::ValidateDistributedSurfaceInterface(surface);
	return surface;
}

iga::DomainNode FsiDomain(const std::string& id, iga::DomainKind kind,
	const iga::SurfaceInterfaceRef& surface_id, const std::string& reference_mesh_identity,
	bool fluid, double area = 1.0)
{
	auto domain = Domain(id, kind, fluid
		? std::vector<iga::CouplingPort>{LogicalPort(id, "flow_port", "flow_port",
			iga::PortQuantity::MeanPressure)}
		: std::vector<iga::CouplingPort>{});
	auto surface = FsiSurface(surface_id, reference_mesh_identity, fluid);
	domain.surface_interfaces = {surface};
	domain.surface_layouts.emplace(surface.id, SurfaceLayout(reference_mesh_identity, area));
	return domain;
}

iga::CouplingPort SpeciesPort(const std::string& domain_id, const std::string& id,
	iga::PortQuantity hydraulic_input)
{
	auto port = LogicalPort(domain_id, id, id, hydraulic_input);
	port.provides.insert(iga::PortQuantity::SpeciesConcentration);
	port.provides.insert(iga::PortQuantity::SpeciesFlux);
	port.requires.insert(iga::PortQuantity::SpeciesConcentration);
	port.requires.insert(iga::PortQuantity::SpeciesFlux);
	port.species = {"drug_parent", "tracer_alpha"};
	iga::ValidateCouplingPort(port);
	return port;
}

iga::SimulationGraph StraightChain()
{
	auto upstream = Domain("upstream", iga::DomainKind::OneDFlow,
		{LogicalPort("upstream", "terminal", "outlet:2", iga::PortQuantity::MeanPressure)});
	auto three_d = Domain("roi", iga::DomainKind::ThreeDBodyFittedFlow,
		{LogicalPort("roi", "inlet", "boundary:1", iga::PortQuantity::FlowRate),
		 LogicalPort("roi", "outlet", "boundary:2", iga::PortQuantity::MeanPressure)});
	auto downstream = Domain("downstream", iga::DomainKind::OneDFlow,
		{LogicalPort("downstream", "root", "root", iga::PortQuantity::FlowRate)});
	iga::CouplingEdge left{"upstream_to_roi", {"upstream", "terminal"}, {"roi", "inlet"},
		iga::CouplingLaw::PressureFlow};
	iga::CouplingEdge right{"roi_to_downstream", {"roi", "outlet"}, {"downstream", "root"},
		iga::CouplingLaw::PressureFlow};
	return iga::SimulationGraph({std::move(downstream), std::move(three_d), std::move(upstream)},
		{std::move(right), std::move(left)});
}

iga::SimulationGraph Bifurcation(bool permuted = false)
{
	auto source = Domain("source", iga::DomainKind::OneDFlow,
		{LogicalPort("source", "terminal", "outlet:2", iga::PortQuantity::MeanPressure)});
	auto junction = Domain("junction", iga::DomainKind::ThreeDBodyFittedFlow,
		{LogicalPort("junction", "inlet", "boundary:1", iga::PortQuantity::FlowRate),
		 LogicalPort("junction", "right", "boundary:3", iga::PortQuantity::MeanPressure),
		 LogicalPort("junction", "left", "boundary:2", iga::PortQuantity::MeanPressure)});
	auto left = Domain("branch_a", iga::DomainKind::OneDFlow,
		{LogicalPort("branch_a", "root", "root", iga::PortQuantity::FlowRate)});
	auto right = Domain("branch_b", iga::DomainKind::OneDFlow,
		{LogicalPort("branch_b", "root", "root", iga::PortQuantity::FlowRate)});
	std::vector<iga::CouplingEdge> edges{
		{"c_right", {"junction", "right"}, {"branch_b", "root"},
			iga::CouplingLaw::PressureFlow},
		{"a_source", {"source", "terminal"}, {"junction", "inlet"},
			iga::CouplingLaw::PressureFlow},
		{"b_left", {"junction", "left"}, {"branch_a", "root"},
			iga::CouplingLaw::PressureFlow}};
	if (permuted) {
		std::reverse(edges.begin(), edges.end());
		for (auto& edge : edges) std::swap(edge.first, edge.second);
	}
	return iga::SimulationGraph({std::move(right), std::move(source), std::move(left),
		std::move(junction)}, std::move(edges));
}

} // namespace

int main()
{
	{
		auto first = Domain("one", iga::DomainKind::OneDFlow,
			{SpeciesPort("one", "outlet:2", iga::PortQuantity::MeanPressure)});
		first.species_bindings = {{"drug_parent", "native_drug"},
			{"tracer_alpha", "native_tracer"}};
		auto second = Domain("three", iga::DomainKind::ThreeDBodyFittedFlow,
			{SpeciesPort("three", "boundary:1", iga::PortQuantity::FlowRate)});
		second.species_bindings = {{"tracer_alpha", "scalar_a"},
			{"drug_parent", "scalar_b"}};
		const std::set<std::string> species{"tracer_alpha", "drug_parent"};
		const iga::SimulationGraph species_graph({first, second},
			{{"edge", {"one", "outlet:2"}, {"three", "boundary:1"},
				iga::CouplingLaw::PressureFlow, species}},
			{{"tracer_alpha", "kg/m^3"}, {"drug_parent", "mol/m^3"}});
		assert(species_graph.Species().begin()->first == "drug_parent");
		assert(species_graph.Edge("edge").species
			== std::set<std::string>({"drug_parent", "tracer_alpha"}));
		const iga::SimulationGraph permuted({second, first},
			{{"edge", {"three", "boundary:1"}, {"one", "outlet:2"},
				iga::CouplingLaw::PressureFlow, species}},
			{{"drug_parent", "mol/m^3"}, {"tracer_alpha", "kg/m^3"}});
		assert(permuted.Species() == species_graph.Species());
		RequireRejected([&first, &second, &species] {
			auto incomplete = second;
			incomplete.species_bindings.erase("drug_parent");
			iga::SimulationGraph({first, incomplete},
				{{"edge", {"one", "outlet:2"}, {"three", "boundary:1"},
					iga::CouplingLaw::PressureFlow, species}},
				{{"tracer_alpha", "kg/m^3"}, {"drug_parent", "mol/m^3"}});
		});
		RequireRejected([] {
			auto unattached = SpeciesPort("observation", "unused",
				iga::PortQuantity::MeanPressure);
			unattached.species = {"unbound_species"};
			iga::SimulationGraph({Domain("observation", iga::DomainKind::OneDFlow,
				{unattached})}, {}, {{"unbound_species", "mol/m^3"}});
		});
	}
	const auto graph = StraightChain();
	assert(graph.Domains().size() == 3);
	assert(graph.Edges().size() == 2);
	assert(graph.Domain("roi").kind == iga::DomainKind::ThreeDBodyFittedFlow);
	assert(graph.Port({"roi", "inlet"}).requires.count(iga::PortQuantity::FlowRate));
	assert(graph.Port({"roi", "inlet"}).provides.count(iga::PortQuantity::FlowRate));
	assert(graph.Edge("upstream_to_roi").law == iga::CouplingLaw::PressureFlow);
	assert(iga::DomainDimensionOf(iga::DomainKind::OneDFlow) == 1);
	assert(iga::DomainDimensionOf(iga::DomainKind::ThreeDBodyFittedFlow) == 3);
	assert(iga::DomainTopologyDimensionOf(iga::DomainKind::ZeroDFlow) == 0);
	assert(iga::DomainEmbeddingDimensionOf(iga::DomainKind::ZeroDFlow) == 0);
	assert(iga::IsFlowDomainKind(iga::DomainKind::ZeroDFlow));
	const auto plan = iga::MakeSequentialPlan(graph, "upstream");
	assert((plan.domain_ids == std::vector<std::string>{"upstream", "roi", "downstream"}));
	assert((plan.edge_ids == std::vector<std::string>{"upstream_to_roi", "roi_to_downstream"}));
	const auto reverse = iga::MakeSequentialPlan(graph, "downstream");
	assert((reverse.domain_ids == std::vector<std::string>{"downstream", "roi", "upstream"}));
	const auto pressure_flow_plan = iga::MakeSequentialPressureFlowPlan(graph, "upstream");
	assert(pressure_flow_plan.domain_ids == plan.domain_ids);
	{
		auto source_port = iga::MakeZeroDFlowPort("source", iga::ZeroDFlowRole::SourceReservoir);
		auto terminal_port = iga::MakeZeroDFlowPort("terminal", iga::ZeroDFlowRole::TerminalRcr);
		auto root = LogicalPort("network", "root", "root", iga::PortQuantity::FlowRate);
		auto outlet = LogicalPort("network", "outlet", "outlet:2", iga::PortQuantity::MeanPressure);
		const iga::SimulationGraph zero_d_chain({Domain("source", iga::DomainKind::ZeroDFlow,
			{source_port}), Domain("network", iga::DomainKind::OneDFlow, {root, outlet}),
			Domain("terminal", iga::DomainKind::ZeroDFlow, {terminal_port})},
			{{"source_to_network", {"source", "port"}, {"network", "root"},
				iga::CouplingLaw::PressureFlow},
			 {"network_to_terminal", {"network", "outlet"}, {"terminal", "port"},
				iga::CouplingLaw::PressureFlow}});
		assert((iga::MakeSequentialPlan(zero_d_chain, "source").domain_ids
			== std::vector<std::string>{"source", "network", "terminal"}));
		const iga::SimulationGraph zero_d_to_three_d({Domain("source", iga::DomainKind::ZeroDFlow,
			{source_port}), Domain("roi", iga::DomainKind::ThreeDBodyFittedFlow,
			{LogicalPort("roi", "inlet", "boundary:1", iga::PortQuantity::FlowRate)})},
			{{"source_to_roi", {"source", "port"}, {"roi", "inlet"},
				iga::CouplingLaw::PressureFlow}});
		assert((iga::MakeSequentialPlan(zero_d_to_three_d, "source").domain_ids
			== std::vector<std::string>{"source", "roi"}));
		RequireRejected([&] {
			iga::SimulationGraph({Domain("source", iga::DomainKind::ZeroDFlow, {source_port}),
				Domain("terminal", iga::DomainKind::ZeroDFlow, {terminal_port})},
				{{"invalid", {"source", "port"}, {"terminal", "port"},
					iga::CouplingLaw::PressureFlow}});
		});
	}
	const auto component_plan = iga::MakeAcyclicPressureFlowPlan(graph, "upstream");
	assert(component_plan.domain_order == plan.domain_ids);
	iga::ValidateConnectedGraph(graph, "upstream");
	RequireRejected([&graph] {
		(void)iga::MakeSequentialPressureFlowPlan(graph, "downstream");
	});
	{
		auto upstream = Domain("upstream", iga::DomainKind::OneDFlow,
			{LogicalPort("upstream", "terminal", "outlet:2", iga::PortQuantity::MeanPressure)});
		auto volume = Domain("volume", iga::DomainKind::ThreeDBodyFittedFlow,
			{LogicalPort("volume", "inlet", "1", iga::PortQuantity::FlowRate),
			 LogicalPort("volume", "outlet", "2", iga::PortQuantity::MeanPressure)});
		auto downstream = Domain("downstream", iga::DomainKind::OneDFlow,
			{LogicalPort("downstream", "root", "root", iga::PortQuantity::FlowRate)});
		const iga::SimulationGraph volume_chain({std::move(upstream), std::move(volume),
			std::move(downstream)}, {{"left", {"upstream", "terminal"}, {"volume", "inlet"},
				iga::CouplingLaw::PressureFlow}, {"right", {"volume", "outlet"},
				{"downstream", "root"}, iga::CouplingLaw::PressureFlow}});
		assert((iga::MakeSequentialPlan(volume_chain, "upstream").domain_ids
			== std::vector<std::string>{"upstream", "volume", "downstream"}));
		auto body = Domain("body", iga::DomainKind::ThreeDBodyFittedFlow,
			{LogicalPort("body", "outlet", "1", iga::PortQuantity::MeanPressure)});
		auto volume_neighbor = Domain("volume", iga::DomainKind::ThreeDBodyFittedFlow,
			{LogicalPort("volume", "inlet", "1", iga::PortQuantity::FlowRate)});
		const iga::SimulationGraph adjacent_three_d({std::move(body),
			std::move(volume_neighbor)}, {{"bad", {"body", "outlet"},
				{"volume", "inlet"}, iga::CouplingLaw::PressureFlow}});
		RequireRejected([&adjacent_three_d] {
			(void)iga::MakeSequentialPlan(adjacent_three_d, "body");
		});
		RequireRejected([&adjacent_three_d] {
			(void)iga::MakeAcyclicPressureFlowPlan(adjacent_three_d, "body");
		});
	}
	{
		std::vector<iga::DomainNode> copied_domains;
		for (const auto& domain : graph.Domains()) copied_domains.push_back(domain.second);
		auto copied_edges = graph.Edges();
		const iga::SimulationGraph immutable(copied_domains, copied_edges);
		copied_domains.front().id = "mutated";
		copied_edges.front().id = "mutated";
		assert(immutable.Domains().count("mutated") == 0);
		assert(immutable.Edge("roi_to_downstream").id == "roi_to_downstream");
	}
	{
		// No graph domain kind implements an FSI endpoint yet, so a well-formed
		// FSI edge between existing kinds is rejected by the pair check itself.
		const auto mesh_identity = SurfaceHash('a');
		const iga::SurfaceInterfaceRef fluid_surface{"fluid", "fluid_runtime", "wall"};
		const iga::SurfaceInterfaceRef structure_surface{"structure", "solid_runtime", "wall"};
		auto structure = FsiDomain("structure", iga::DomainKind::OneDFlow,
			structure_surface, mesh_identity, false);
		structure.ports = {LogicalPort("structure", "flow_port", "flow_port",
			iga::PortQuantity::MeanPressure)};
		bool rejected_pair = false;
		try {
			iga::SimulationGraph({FsiDomain("fluid", iga::DomainKind::ThreeDBodyFittedFlow,
				fluid_surface, mesh_identity, true), std::move(structure)}, {}, {},
				{{"fluid_structure", fluid_surface, structure_surface,
					iga::FsiCouplingLaw::FluidStructureTractionKinematics}});
		} catch (const std::runtime_error& error) {
			rejected_pair = std::string(error.what()).find(
				"unsupported fluid/structure domain pair") != std::string::npos;
		}
		assert(rejected_pair);
	}

	RequireRejected([&graph] { (void)graph.Domain("missing"); });
	RequireRejected([&graph] { (void)graph.Port({"roi", "missing"}); });
	RequireRejected([&graph] { (void)graph.Edge("missing"); });
	RequireRejected([&graph] { (void)iga::MakeSequentialPlan(graph, "missing"); });
	RequireRejected([&graph] { (void)iga::MakeSequentialPlan(graph, "roi"); });
	RequireRejected([] { iga::SimulationGraph({}, {}); });
	RequireRejected([] {
		iga::SimulationGraph({{"", iga::DomainKind::OneDFlow, {}}}, {});
	});
	RequireRejected([] {
		iga::SimulationGraph({{"empty", iga::DomainKind::OneDFlow, {}}}, {});
	});
	RequireRejected([] {
		auto port = LogicalPort("unknown", "port", "port", iga::PortQuantity::FlowRate);
		iga::SimulationGraph({Domain("unknown", static_cast<iga::DomainKind>(99), {port})}, {});
	});
	RequireRejected([] {
		auto port = LogicalPort("one", "terminal", "outlet:2", iga::PortQuantity::MeanPressure);
		iga::SimulationGraph({Domain("one", iga::DomainKind::OneDFlow, {port}),
			Domain("one", iga::DomainKind::OneDFlow, {port})}, {});
	});
	RequireRejected([] {
		auto port = LogicalPort("wrong", "terminal", "outlet:2", iga::PortQuantity::MeanPressure);
		iga::SimulationGraph({Domain("one", iga::DomainKind::OneDFlow, {port})}, {});
	});
	RequireRejected([] {
		auto first = LogicalPort("one", "a", "same", iga::PortQuantity::MeanPressure);
		auto second = LogicalPort("one", "b", "same", iga::PortQuantity::FlowRate);
		iga::SimulationGraph({Domain("one", iga::DomainKind::OneDFlow, {first, second})}, {});
	});
	RequireRejected([] {
		auto first = LogicalPort("one", "a", "a", iga::PortQuantity::MeanPressure);
		auto second = LogicalPort("two", "b", "b", iga::PortQuantity::FlowRate);
		iga::SimulationGraph({Domain("one", iga::DomainKind::OneDFlow, {first}),
			Domain("two", iga::DomainKind::ThreeDBodyFittedFlow, {second})},
			{{"bad", {"one", "missing"}, {"two", "b"}, iga::CouplingLaw::PressureFlow}});
	});
	RequireRejected([] {
		auto first = LogicalPort("one", "a", "a", iga::PortQuantity::MeanPressure);
		auto second = LogicalPort("two", "b", "b", iga::PortQuantity::FlowRate);
		iga::SimulationGraph({Domain("one", iga::DomainKind::OneDFlow, {first}),
			Domain("two", iga::DomainKind::ThreeDBodyFittedFlow, {second})},
			{{"", {"one", "a"}, {"two", "b"}, iga::CouplingLaw::PressureFlow}});
	});
	RequireRejected([] {
		auto first = LogicalPort("one", "a", "a", iga::PortQuantity::MeanPressure);
		auto second = LogicalPort("two", "b", "b", iga::PortQuantity::FlowRate);
		iga::SimulationGraph({Domain("one", iga::DomainKind::OneDFlow, {first}),
			Domain("two", iga::DomainKind::ThreeDBodyFittedFlow, {second})},
			{{"empty_endpoint", {"", "a"}, {"two", "b"}, iga::CouplingLaw::PressureFlow}});
	});
	RequireRejected([] {
		auto first = LogicalPort("one", "a", "a", iga::PortQuantity::MeanPressure);
		auto second = LogicalPort("two", "b", "b", iga::PortQuantity::MeanPressure);
		iga::SimulationGraph({Domain("one", iga::DomainKind::OneDFlow, {first}),
			Domain("two", iga::DomainKind::ThreeDBodyFittedFlow, {second})},
			{{"ambiguous", {"one", "a"}, {"two", "b"}, iga::CouplingLaw::PressureFlow}});
	});
	RequireRejected([] {
		auto first = LogicalPort("one", "a", "a", iga::PortQuantity::MeanPressure);
		auto second = LogicalPort("two", "b", "b", iga::PortQuantity::FlowRate);
		second.provides.erase(iga::PortQuantity::MeanPressure);
		iga::SimulationGraph({Domain("one", iga::DomainKind::OneDFlow, {first}),
			Domain("two", iga::DomainKind::ThreeDBodyFittedFlow, {second})},
			{{"missing_diagnostic", {"one", "a"}, {"two", "b"}, iga::CouplingLaw::PressureFlow}});
	});
	RequireRejected([] {
		auto first = LogicalPort("one", "a", "a", iga::PortQuantity::MeanPressure);
		auto second = LogicalPort("two", "b", "b", iga::PortQuantity::FlowRate);
		iga::SimulationGraph({Domain("one", iga::DomainKind::OneDFlow, {first}),
			Domain("two", iga::DomainKind::ThreeDBodyFittedFlow, {second})},
			{{"same", {"one", "a"}, {"two", "b"}, iga::CouplingLaw::PressureFlow},
			 {"same", {"one", "a"}, {"two", "b"}, iga::CouplingLaw::PressureFlow}});
	});
	RequireRejected([] {
		auto first = LogicalPort("one", "a", "a", iga::PortQuantity::MeanPressure);
		iga::SimulationGraph({Domain("one", iga::DomainKind::OneDFlow, {first})},
			{{"self", {"one", "a"}, {"one", "a"}, iga::CouplingLaw::PressureFlow}});
	});
	RequireRejected([] {
		auto first = LogicalPort("one", "a", "a", iga::PortQuantity::MeanPressure);
		auto second = LogicalPort("two", "b", "b", iga::PortQuantity::FlowRate);
		iga::SimulationGraph({Domain("one", iga::DomainKind::OneDFlow, {first}),
			Domain("two", iga::DomainKind::ThreeDBodyFittedFlow, {second})},
			{{"unknown_law", {"one", "a"}, {"two", "b"},
				static_cast<iga::CouplingLaw>(99)}});
	});
	RequireRejected([] {
		auto source = Domain("source", iga::DomainKind::OneDFlow,
			{LogicalPort("source", "out", "out", iga::PortQuantity::MeanPressure)});
		auto first = Domain("first", iga::DomainKind::ThreeDBodyFittedFlow,
			{LogicalPort("first", "in", "in", iga::PortQuantity::FlowRate)});
		auto second = Domain("second", iga::DomainKind::ThreeDBodyFittedFlow,
			{LogicalPort("second", "in", "in", iga::PortQuantity::FlowRate)});
		iga::SimulationGraph({std::move(source), std::move(first), std::move(second)},
			{{"first_edge", {"source", "out"}, {"first", "in"}, iga::CouplingLaw::PressureFlow},
			 {"second_edge", {"source", "out"}, {"second", "in"}, iga::CouplingLaw::PressureFlow}});
	});

	{
		const auto branch = Bifurcation();
		const auto branch_plan = iga::MakeAcyclicPressureFlowPlan(branch, "source");
		assert((branch_plan.domain_order == std::vector<std::string>{"source", "junction",
			"branch_a", "branch_b"}));
		assert(branch_plan.interfaces.size() == 3);
		assert(branch_plan.interfaces[0].edge_id == "a_source");
		assert(branch_plan.interfaces[1].edge_id == "b_left");
		assert(branch_plan.interfaces[2].edge_id == "c_right");
		assert(branch_plan.interfaces[1].flow_provider
			== iga::PortRef({"junction", "left"}));
		assert(branch_plan.interfaces[1].flow_receiver
			== iga::PortRef({"branch_a", "root"}));
		const auto permuted_plan = iga::MakeAcyclicPressureFlowPlan(Bifurcation(true), "source");
		assert(permuted_plan.domain_order == branch_plan.domain_order);
		for (std::size_t i = 0; i < branch_plan.interfaces.size(); ++i) {
			assert(permuted_plan.interfaces[i].edge_id == branch_plan.interfaces[i].edge_id);
			assert(permuted_plan.interfaces[i].flow_provider
				== branch_plan.interfaces[i].flow_provider);
		}
		RequireRejected([&branch] {
			(void)iga::MakeAcyclicPressureFlowPlan(branch, "branch_a");
		});
	}
	{
		auto center = Domain("center", iga::DomainKind::ThreeDBodyFittedFlow,
			{LogicalPort("center", "a", "a", iga::PortQuantity::FlowRate),
			 LogicalPort("center", "b", "b", iga::PortQuantity::FlowRate),
			 LogicalPort("center", "c", "c", iga::PortQuantity::FlowRate)});
		std::vector<iga::DomainNode> domains{std::move(center)};
		std::vector<iga::CouplingEdge> edges;
		for (const auto& id : {"a", "b", "c"}) {
			domains.push_back(Domain(id, iga::DomainKind::OneDFlow,
				{LogicalPort(id, "terminal", "terminal", iga::PortQuantity::MeanPressure)}));
			edges.push_back({std::string("edge_")+id, {"center", id}, {id, "terminal"},
				iga::CouplingLaw::PressureFlow});
		}
		const iga::SimulationGraph branch(std::move(domains), std::move(edges));
		assert(branch.Domains().size() == 4);
		iga::ValidateConnectedGraph(branch, "a");
		RequireRejected([&branch] { (void)iga::MakeSequentialPlan(branch, "a"); });
		RequireRejected([&branch] {
			(void)iga::MakeSequentialPressureFlowPlan(branch, "a");
		});
		RequireRejected([&branch] {
			(void)iga::MakeAcyclicPressureFlowPlan(branch, "a");
		});
	}
	{
		auto graph_with_island = StraightChain();
		std::vector<iga::DomainNode> domains;
		for (const auto& domain : graph_with_island.Domains()) domains.push_back(domain.second);
		domains.push_back(Domain("island", iga::DomainKind::ThreeDBodyFittedFlow,
			{LogicalPort("island", "unused", "unused", iga::PortQuantity::FlowRate)}));
		iga::SimulationGraph disconnected(std::move(domains), graph_with_island.Edges());
		RequireRejected([&disconnected] {
			iga::ValidateConnectedGraph(disconnected, "upstream");
		});
		RequireRejected([&disconnected] {
			(void)iga::MakeSequentialPlan(disconnected, "upstream");
		});
	}
	{
		auto left = LogicalPort("left", "a", "a", iga::PortQuantity::MeanPressure);
		auto right = LogicalPort("right", "b", "b", iga::PortQuantity::FlowRate);
		iga::SimulationGraph same_kind({Domain("left", iga::DomainKind::OneDFlow, {left}),
			Domain("right", iga::DomainKind::OneDFlow, {right})},
			{{"same_kind", {"left", "a"}, {"right", "b"}, iga::CouplingLaw::PressureFlow}});
		RequireRejected([&same_kind] { (void)iga::MakeSequentialPlan(same_kind, "left"); });
		RequireRejected([&same_kind] {
			(void)iga::MakeAcyclicPressureFlowPlan(same_kind, "left");
		});
	}
	{
		std::vector<iga::DomainNode> domains;
		for (int index = 0; index < 4; ++index) {
			const std::string id = "cycle_"+std::to_string(index);
			const bool one_d = index%2 == 0;
			domains.push_back(Domain(id, one_d ? iga::DomainKind::OneDFlow
				: iga::DomainKind::ThreeDBodyFittedFlow,
				{LogicalPort(id, "left", "left", one_d ? iga::PortQuantity::MeanPressure
					: iga::PortQuantity::FlowRate),
				 LogicalPort(id, "right", "right", one_d ? iga::PortQuantity::MeanPressure
					: iga::PortQuantity::FlowRate)}));
		}
		const iga::SimulationGraph cycle(std::move(domains), {
			{"cycle_edge_0", {"cycle_0", "right"}, {"cycle_1", "left"}, iga::CouplingLaw::PressureFlow},
			{"cycle_edge_1", {"cycle_1", "right"}, {"cycle_2", "left"}, iga::CouplingLaw::PressureFlow},
			{"cycle_edge_2", {"cycle_2", "right"}, {"cycle_3", "left"}, iga::CouplingLaw::PressureFlow},
			{"cycle_edge_3", {"cycle_3", "right"}, {"cycle_0", "left"}, iga::CouplingLaw::PressureFlow}});
		RequireRejected([&cycle] { (void)iga::MakeSequentialPlan(cycle, "cycle_0"); });
		RequireRejected([&cycle] {
			(void)iga::MakeAcyclicPressureFlowPlan(cycle, "cycle_0");
		});
	}
	return 0;
}
