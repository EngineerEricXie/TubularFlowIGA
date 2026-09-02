#include "DomainRuntimeRegistry.hpp"
#include "ImmersedFlowCase.hpp"
#include "SimulationGraph.hpp"
#include "ThreeDImmersedFlowDomain.hpp"

#include <cassert>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

template <class Function> void Reject(Function&& function)
{
	bool rejected = false;
	try { function(); } catch (const std::exception&) { rejected = true; }
	assert(rejected);
}

iga::CouplingPort Port(const char* id, int label, iga::PortQuantity input)
{
	iga::CouplingPort port;
	port.id = id;
	port.subsystem_id = "immersed";
	port.locator_kind = "boundary_label";
	port.locator = std::to_string(label);
	port.provides = {iga::PortQuantity::Area, iga::PortQuantity::FlowRate,
		iga::PortQuantity::MeanPressure, iga::PortQuantity::MeanNormalTraction};
	port.requires = {input};
	return port;
}

void Write(const fs::path& path, const std::string& text)
{
	std::ofstream output(path);
	if (!output) throw std::runtime_error("cannot create test fixture");
	output << text;
}

void WriteCase(const fs::path& directory, const std::string& geometry)
{
	Write(directory/"surface.vtp", R"xml(<?xml version="1.0"?>
<VTKFile type="PolyData" version="1.0" byte_order="LittleEndian" header_type="UInt64"><PolyData><Piece NumberOfPoints="8" NumberOfPolys="12"><Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">.1 .1 .1 .9 .1 .1 .9 .9 .1 .1 .9 .1 .1 .1 .9 .9 .1 .9 .9 .9 .9 .1 .9 .9</DataArray></Points><Polys><DataArray type="Int32" Name="connectivity" format="ascii">0 2 1 0 3 2 4 5 6 4 6 7 0 1 5 0 5 4 1 2 6 1 6 5 2 3 7 2 7 6 3 0 4 3 4 7</DataArray><DataArray type="Int32" Name="offsets" format="ascii">3 6 9 12 15 18 21 24 27 30 33 36</DataArray></Polys><CellData Scalars="boundary_id"><DataArray type="UInt32" Name="boundary_id" format="ascii">1 1 2 2 0 0 0 0 0 0 0 0</DataArray></CellData></Piece></PolyData></VTKFile>)xml");
	Write(directory/"simulation_config.json", R"json({
"schema_version":3,"dimension":"3d","fields":[{"name":"velocity","kind":"vector3"},{"name":"pressure","kind":"pressure"}],"time":{"dt":1,"steps":1},"equation_systems":[{"name":"flow","kind":"navier_stokes","unknowns":["velocity","pressure"],"viscosity":1,"density":1,"time_integration":"steady"}],"boundaries":[]})json");
	Write(directory/"immersed_geometry.json", geometry);
}

void WriteDisconnectedClosedSurface(const fs::path& directory)
{
	Write(directory/"surface.vtp", R"xml(<?xml version="1.0"?>
<VTKFile type="PolyData" version="1.0" byte_order="LittleEndian"><PolyData><Piece NumberOfPoints="16" NumberOfPolys="24"><Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">.1 .1 .1 .9 .1 .1 .9 .9 .1 .1 .9 .1 .1 .1 .9 .9 .1 .9 .9 .9 .9 .1 .9 .9 .15 .15 .15 .35 .15 .15 .35 .35 .15 .15 .35 .15 .15 .15 .35 .35 .15 .35 .35 .35 .35 .15 .35 .35</DataArray></Points><Polys><DataArray type="Int32" Name="connectivity" format="ascii">0 2 1 0 3 2 4 5 6 4 6 7 0 1 5 0 5 4 1 2 6 1 6 5 2 3 7 2 7 6 3 0 4 3 4 7 8 10 9 8 11 10 12 13 14 12 14 15 8 9 13 8 13 12 9 10 14 9 14 13 10 11 15 10 15 14 11 8 12 11 12 15</DataArray><DataArray type="Int32" Name="offsets" format="ascii">3 6 9 12 15 18 21 24 27 30 33 36 39 42 45 48 51 54 57 60 63 66 69 72</DataArray></Polys><CellData Scalars="boundary_id"><DataArray type="UInt32" Name="boundary_id" format="ascii">1 1 2 2 0 0 0 0 0 0 0 0 1 1 2 2 0 0 0 0 0 0 0 0</DataArray></CellData></Piece></PolyData></VTKFile>)xml");
}

std::string ReplaceOnce(std::string value, const std::string& from, const std::string& to)
{
	const auto position = value.find(from);
	if (position == std::string::npos) throw std::runtime_error("fixture replacement is missing");
	value.replace(position, from.size(), to);
	return value;
}

void RewriteSimulation(const fs::path& directory, const std::string& from, const std::string& to)
{
	std::ifstream input(directory/"simulation_config.json");
	std::string text((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
	Write(directory/"simulation_config.json", ReplaceOnce(text, from, to));
}

std::string Geometry(const std::string& walls = "[0]")
{
	return "{\"surface\":\"surface.vtp\",\"grid\":{\"lower_m\":[0,0,0],\"upper_m\":[1,1,1],\"cells\":[3,3,3]},\"volume_quadrature\":{\"storage\":\"compact\",\"max_depth\":4,\"max_nodes\":500000,\"max_leaves\":500000,\"max_points\":3000000,\"max_records\":3000000,\"max_retained_bytes\":30000000,\"max_logical_points\":3000000},\"surface_quadrature\":{\"max_candidates\":500000,\"max_fragments\":500000,\"max_points\":3000000,\"max_exact_limbs\":512},\"ghost_penalty\":{\"gamma_u\":0.01,\"gamma_p\":0.01,\"max_faces\":500000,\"max_quadrature_points\":3000000,\"max_trace_entries\":128},\"wall_labels\":"+walls+",\"runtime\":{\"lu_pivot_shift\":1e-12,\"ports\":[{\"id\":\"inlet\",\"boundary_label\":1,\"control_mode\":\"flow_rate\",\"value\":-0.001},{\"id\":\"outlet\",\"boundary_label\":2,\"control_mode\":\"pressure\",\"value\":0}]}}";
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int status = 0;
	try {
		const auto nonce = std::chrono::steady_clock::now().time_since_epoch().count();
		const fs::path root = fs::temp_directory_path()/
			("tubularflowiga_immersed_case_"+std::to_string(nonce));
		fs::create_directories(root);
		const std::vector<iga::CouplingPort> ports = {
			Port("inlet", 1, iga::PortQuantity::FlowRate),
			Port("outlet", 2, iga::PortQuantity::MeanPressure)};
		WriteCase(root, Geometry());
		auto owner = iga::ImmersedFlowCase::Load(root, "immersed", ports, 1);
		assert(!owner->SurfaceHash().empty());
		assert(owner->Grid().cells[0] == 3);
		assert(owner->ClassificationDiagnostics().cut_count > 0);
		assert(owner->SurfaceDiagnostics().catalog_usable);
		assert(owner->RuntimeParameters().density == owner->Configuration().equation_systems[0].density);
		assert(owner->RuntimeParameters().dynamic_viscosity
			== owner->Configuration().equation_systems[0].viscosity);
		std::vector<std::unique_ptr<iga::CoupledDomainRuntime>> runtimes;
		runtimes.push_back(std::move(owner));
		iga::SimulationGraph graph({iga::DomainNode("immersed",
			iga::DomainKind::ThreeDImmersedFlow, ports)}, {});
		iga::DomainRuntimeRegistry registry(graph, std::move(runtimes));
		assert(registry.Runtime("immersed").Kind() == iga::DomainKind::ThreeDImmersedFlow);
		Reject([&] { (void)iga::ImmersedFlowCase::Load(root, "immersed", ports, 2); });
		WriteCase(root, Geometry("[]"));
		Reject([&] { (void)iga::ImmersedFlowCase::Load(root, "immersed", ports, 1); });
		WriteCase(root, Geometry("[0,1]"));
		Reject([&] { (void)iga::ImmersedFlowCase::Load(root, "immersed", ports, 1); });
		WriteCase(root, ReplaceOnce(Geometry(), "\"surface\":\"surface.vtp\"",
			"\"surface\":\"../surface.vtp\""));
		Reject([&] { (void)iga::ImmersedFlowCase::Load(root, "immersed", ports, 1); });
		WriteCase(root, Geometry());
		Write(root/"surface.vtp", R"xml(<?xml version="1.0"?><VTKFile type="PolyData" version="1.0" byte_order="LittleEndian"><PolyData><Piece NumberOfPoints="3" NumberOfPolys="1"><Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">0 0 0 1 0 0 0 1 0</DataArray></Points><Polys><DataArray type="Int32" Name="connectivity" format="ascii">0 1 2</DataArray><DataArray type="Int32" Name="offsets" format="ascii">3</DataArray></Polys><CellData Scalars="boundary_id"><DataArray type="UInt32" Name="boundary_id" format="ascii">0</DataArray></CellData></Piece></PolyData></VTKFile>)xml");
		Reject([&] { (void)iga::ImmersedFlowCase::Load(root, "immersed", ports, 1); });
		WriteCase(root, Geometry());
		WriteDisconnectedClosedSurface(root);
		Reject([&] { (void)iga::ImmersedFlowCase::Load(root, "immersed", ports, 1); });
		WriteCase(root, Geometry());
		RewriteSimulation(root, "\"time_integration\":\"steady\"",
			"\"time_integration\":\"backward_euler\"");
		Reject([&] { (void)iga::ImmersedFlowCase::Load(root, "immersed", ports, 1); });
		WriteCase(root, Geometry());
		RewriteSimulation(root, "{\"name\":\"pressure\",\"kind\":\"pressure\"}",
			"{\"name\":\"pressure\",\"kind\":\"pressure\"},{\"name\":\"species\",\"kind\":\"scalar\"}");
		Reject([&] { (void)iga::ImmersedFlowCase::Load(root, "immersed", ports, 1); });
		WriteCase(root, Geometry());
		RewriteSimulation(root, "\"boundaries\":[]", "\"boundaries\":[{\"label\":2,\"conditions\":[{\"field\":\"pressure\",\"type\":\"resistance\",\"resistance\":1}]}]");
		Reject([&] { (void)iga::ImmersedFlowCase::Load(root, "immersed", ports, 1); });
		WriteCase(root, Geometry());
		RewriteSimulation(root, "\"boundaries\":[]", "\"boundaries\":[{\"label\":1,\"conditions\":[{\"field\":\"velocity\",\"type\":\"dirichlet\",\"profile\":\"profile.txt\",\"value\":[0,0,0]}]}]");
		Reject([&] { (void)iga::ImmersedFlowCase::Load(root, "immersed", ports, 1); });
		WriteCase(root, Geometry());
		RewriteSimulation(root, "\"boundaries\":[]", "\"boundaries\":[{\"label\":2,\"conditions\":[{\"field\":\"pressure\",\"type\":\"pressure_traction\",\"value\":0,\"pressure_gauge\":true}]}]");
		Reject([&] { (void)iga::ImmersedFlowCase::Load(root, "immersed", ports, 1); });
		WriteCase(root, ReplaceOnce(Geometry(), "\"control_mode\":\"pressure\"",
			"\"control_mode\":\"total_pressure\""));
		Reject([&] { (void)iga::ImmersedFlowCase::Load(root, "immersed", ports, 1); });
		WriteCase(root, ReplaceOnce(Geometry(), "\"cells\":[3,3,3]", "\"cells\":[3.5,3,3]"));
		Reject([&] { (void)iga::ImmersedFlowCase::Load(root, "immersed", ports, 1); });
		WriteCase(root, ReplaceOnce(Geometry(), "\"max_depth\":4", "\"max_depth\":4294967296"));
		Reject([&] { (void)iga::ImmersedFlowCase::Load(root, "immersed", ports, 1); });
		WriteCase(root, ReplaceOnce(Geometry(), "\"gamma_u\":0.01", "\"gamma_u\":-0.01"));
		Reject([&] { (void)iga::ImmersedFlowCase::Load(root, "immersed", ports, 1); });
		WriteCase(root, ReplaceOnce(Geometry(), "\"lu_pivot_shift\":1e-12",
			"\"nonlinear_maximum_iterations\":1.5,\"ksp_maximum_iterations\":2147483648,\"lu_pivot_shift\":1e-12"));
		Reject([&] { (void)iga::ImmersedFlowCase::Load(root, "immersed", ports, 1); });
		WriteCase(root, ReplaceOnce(Geometry(), "\"id\":\"inlet\"", "\"id\":\"wrong\""));
		Reject([&] { (void)iga::ImmersedFlowCase::Load(root, "immersed", ports, 1); });
		WriteCase(root, ReplaceOnce(Geometry(), "\"boundary_label\":1", "\"boundary_label\":3"));
		Reject([&] { (void)iga::ImmersedFlowCase::Load(root, "immersed", ports, 1); });
		WriteCase(root, ReplaceOnce(Geometry(), "\"control_mode\":\"flow_rate\"",
			"\"control_mode\":\"pressure\""));
		Reject([&] { (void)iga::ImmersedFlowCase::Load(root, "immersed", ports, 1); });
		WriteCase(root, ReplaceOnce(Geometry(), "\"lu_pivot_shift\":1e-12",
			"\"density\":2,\"lu_pivot_shift\":1e-12"));
		Reject([&] { (void)iga::ImmersedFlowCase::Load(root, "immersed", ports, 1); });
		fs::remove_all(root);
	} catch (...) { status = 1; }
	PetscFinalize();
	return status;
}
