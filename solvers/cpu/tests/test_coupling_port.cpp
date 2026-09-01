#include "CouplingPort.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

template <class Function>
void RequireRejected(Function&& function, const std::string& description)
{
	bool rejected = false;
	try {
		function();
	} catch (const std::runtime_error&) {
		rejected = true;
	}
	if (!rejected) throw std::runtime_error("expected rejection: "+description);
}

bool Close(double first, double second)
{
	return std::abs(first-second) <= 1.0e-12*std::max({1.0, std::abs(first), std::abs(second)});
}

} // namespace

int main()
{
	iga::PortBoundaryData zero_input;
	zero_input.time_s = 0.0;
	zero_input.outward_flow_m3_s = 0.0;
	zero_input.mean_pressure_pa = 0.0;
	zero_input.mean_normal_traction_pa = 0.0;
	zero_input.total_pressure_pa = 0.0;
	zero_input.concentration["tracer"] = 0.0;
	zero_input.outward_species_flux["tracer"] = 0.0;
	iga::ValidatePortBoundaryData(zero_input);
	assert(zero_input.outward_flow_m3_s.has_value());
	assert(*zero_input.outward_flow_m3_s == 0.0);

	iga::PortBoundaryData missing_input;
	iga::ValidatePortBoundaryData(missing_input);
	assert(!missing_input.outward_flow_m3_s.has_value());
	assert(!missing_input.mean_pressure_pa.has_value());

	iga::PortState state;
	state.time_s = 0.1;
	state.area_m2 = 2.0;
	state.outward_flow_m3_s = -3.0;
	state.mean_pressure_pa = 0.0;
	state.concentration["tracer"] = 0.0;
	state.outward_species_flux["tracer"] = -4.0;
	iga::ValidatePortState(state);

	RequireRejected([] {
		iga::PortState invalid;
		invalid.area_m2 = 0.0;
		iga::ValidatePortState(invalid);
	}, "zero area");
	RequireRejected([] {
		iga::PortState invalid;
		invalid.area_m2 = -1.0;
		iga::ValidatePortState(invalid);
	}, "negative area");
	RequireRejected([] {
		iga::PortState invalid;
		invalid.mean_pressure_pa = std::numeric_limits<double>::infinity();
		iga::ValidatePortState(invalid);
	}, "non-finite state value");
	RequireRejected([] {
		iga::PortBoundaryData invalid;
		invalid.outward_flow_m3_s = std::numeric_limits<double>::quiet_NaN();
		iga::ValidatePortBoundaryData(invalid);
	}, "non-finite boundary value");

	iga::PortOrientation root_orientation{-1};
	assert(Close(root_orientation.ToOutward(2.5), -2.5));
	assert(Close(root_orientation.ToNative(-2.5), 2.5));
	iga::PortOrientation three_d_orientation{1};
	assert(Close(three_d_orientation.ToOutward(-2.5), -2.5));
	RequireRejected([] { iga::PortOrientation{0}.Validate(); }, "zero orientation sign");
	RequireRejected([] { iga::PortOrientation{2}.Validate(); }, "non-unit orientation sign");

	iga::CouplingPort port;
	port.id = "distal";
	port.subsystem_id = "one_d";
	port.locator_kind = "network_node";
	port.locator = "42";
	port.orientation.native_to_outward_sign = 1;
	port.provides = {iga::PortQuantity::FlowRate, iga::PortQuantity::MeanPressure};
	port.requires = {iga::PortQuantity::SpeciesConcentration};
	port.species = {"tracer"};
	iga::ValidateCouplingPort(port);
	RequireRejected([&port] {
		auto invalid = port;
		invalid.id.clear();
		iga::ValidateCouplingPort(invalid);
	}, "empty port id");
	RequireRejected([&port] {
		auto invalid = port;
		invalid.subsystem_id.clear();
		iga::ValidateCouplingPort(invalid);
	}, "empty subsystem id");
	RequireRejected([&port] {
		auto invalid = port;
		invalid.locator_kind.clear();
		iga::ValidateCouplingPort(invalid);
	}, "empty locator kind");
	RequireRejected([&port] {
		auto invalid = port;
		invalid.locator.clear();
		iga::ValidateCouplingPort(invalid);
	}, "empty locator");
	auto bidirectional = port;
	bidirectional.provides.insert(iga::PortQuantity::SpeciesConcentration);
	iga::ValidateCouplingPort(bidirectional);
	assert(bidirectional.provides.count(iga::PortQuantity::SpeciesConcentration));
	assert(bidirectional.requires.count(iga::PortQuantity::SpeciesConcentration));
	RequireRejected([] {
		iga::ValidatePortQuantityDeclarations(
			{iga::PortQuantity::FlowRate, iga::PortQuantity::FlowRate}, {});
	}, "duplicate provided quantity declaration");
	iga::ValidatePortQuantityDeclarations({iga::PortQuantity::FlowRate},
		{iga::PortQuantity::FlowRate});
	RequireRejected([] {
		iga::ValidatePortQuantityDeclarations(
			{static_cast<iga::PortQuantity>(99)}, {});
	}, "unknown quantity declaration");
	RequireRejected([&port] {
		auto duplicate = port;
		iga::ValidateCouplingPorts({port, duplicate});
	}, "duplicate port id within one subsystem");
	auto same_id_other_subsystem = port;
	same_id_other_subsystem.subsystem_id = "three_d";
	iga::ValidateCouplingPorts({port, same_id_other_subsystem});
	RequireRejected([] {
		iga::PortBoundaryData invalid;
		invalid.total_pressure_pa = std::numeric_limits<double>::infinity();
		iga::ValidatePortBoundaryData(invalid);
	}, "non-finite total-pressure boundary value");
	RequireRejected([] {
		iga::PortState invalid;
		invalid.concentration[""] = 1.0;
		iga::ValidatePortState(invalid);
	}, "empty species name");

	iga::PortState first;
	first.outward_flow_m3_s = 1.0e-6;
	first.outward_species_flux["tracer"] = 2.0e-6;
	first.outward_species_flux["drug_parent"] = 3.0e-9;
	iga::PortState second;
	second.outward_flow_m3_s = -1.0e-6;
	second.outward_species_flux["tracer"] = -2.0e-6;
	second.outward_species_flux["drug_parent"] = -3.0e-9;
	const auto residual = iga::ConservativeEdgeResidual(first, second,
		{"drug_parent", "tracer"});
	assert(Close(residual.outward_flow_m3_s, 0.0));
	assert(Close(residual.outward_species_flux.at("tracer"), 0.0));
	assert(Close(residual.outward_species_flux.at("drug_parent"), 0.0));
	second.outward_flow_m3_s = -0.75e-6;
	const auto nonconservative = iga::ConservativeEdgeResidual(first, second);
	assert(Close(nonconservative.outward_flow_m3_s, 0.25e-6));
	RequireRejected([] {
		iga::PortState first;
		iga::PortState second;
		second.outward_flow_m3_s = 0.0;
		iga::ConservativeEdgeResidual(first, second);
	}, "missing edge flow");
	RequireRejected([] {
		iga::PortState first;
		first.outward_flow_m3_s = 0.0;
		first.outward_species_flux["tracer"] = 0.0;
		iga::PortState second;
		second.outward_flow_m3_s = 0.0;
		iga::ConservativeEdgeResidual(first, second);
	}, "mismatched edge species");
	const iga::SpeciesRoutingControls routing{1.0e-12, 1.0e-14, 1.0e-10};
	assert(iga::ResolveSpeciesDonor(2.0e-6, -2.0e-6, routing)
		== iga::SpeciesDonor::First);
	assert(iga::ResolveSpeciesDonor(-3.0e-6, 3.0e-6, routing)
		== iga::SpeciesDonor::Second);
	assert(iga::ResolveSpeciesDonor(0.5e-12, -0.5e-12, routing,
		iga::SpeciesDonor::Second) == iga::SpeciesDonor::Second);
	RequireRejected([&routing] {
		(void)iga::ResolveSpeciesDonor(0.0, 0.0, routing);
	}, "near-zero routing without committed ownership");
	RequireRejected([&routing] {
		(void)iga::ResolveSpeciesDonor(2.0e-6, 1.0e-6, routing);
	}, "same-sign species flows");
	RequireRejected([&routing] {
		(void)iga::ResolveSpeciesDonor(2.0e-6, -1.0e-6, routing);
	}, "nonconservative species flows");
	const auto registry = iga::MakeSpeciesRegistry({{"drug_parent", "mol/m^3"},
		{"tracer_alpha", "kg/m^3"}});
	assert(registry.begin()->first == "drug_parent");
	assert(iga::SpeciesFluxUnit(registry.at("tracer_alpha")) == "(kg/m^3)*m^3/s");
	RequireRejected([] {
		(void)iga::MakeSpeciesRegistry({{"oxygen", "mol/m^3"},
			{"oxygen", "kg/m^3"}});
	}, "duplicate logical species");

	std::cout << "coupling port tests passed\n";
}
