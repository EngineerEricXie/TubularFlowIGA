#include "OneDOutput.hpp"

#include <hdf5.h>

#include <cassert>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <vector>

namespace fs = std::filesystem;

void CheckTopology(const fs::path& directory, int subdivisions, bool cycle)
{
	const auto swc = directory/"bifurcation.swc";
	{
		std::ofstream output(swc);
		output << "1 2 0 0 0 0.001 -1\n"
			<< "2 2 0.01 0.01 0 0.001 1\n"
			<< "3 2 0.02 0.02 0 0.001 2\n"
			<< "4 2 0.02 0 0 0.001 2\n";
	}
	auto network = iga::ReadOneDNetwork(swc, 1.0, subdivisions, 0.004);
	if (cycle) {
		auto edge = network.segments.back();
		edge.index = static_cast<int>(network.segments.size());
		edge.parent = 2;
		edge.child = 3;
		edge.cell_offset = network.cells;
		network.cells += edge.cells;
		network.segments.push_back(edge);
	}
	const auto grid = iga::OneDVtkHdfGrid(network);
	iga::hdf_detail::ValidateUnstructuredGrid(grid);
	assert(grid.points.size() == network.nodes.size()+network.segments.size()*(subdivisions-1));
	assert(grid.types.size() == static_cast<std::size_t>(network.cells));
	for (std::size_t node = 0; node < network.nodes.size(); ++node)
		assert(grid.points[node] == network.nodes[node].position);
	std::vector<int> degree(grid.points.size(), 0);
	const auto centers = iga::OneDCellCenters(network);
	for (const auto& segment : network.segments) {
		const auto first = static_cast<std::size_t>(segment.cell_offset);
		assert(grid.connectivity[2*first] == segment.parent);
		assert(grid.connectivity[2*(first+segment.cells)-1] == segment.child);
		for (int local = 0; local < segment.cells; ++local) {
			const auto cell = first+local;
			const auto left = grid.connectivity[2*cell];
			const auto right = grid.connectivity[2*cell+1];
			++degree[left];
			++degree[right];
			assert(grid.types[cell] == 3);
			assert(grid.offsets[cell] == static_cast<std::int64_t>(2*cell));
			for (std::size_t axis = 0; axis < 3; ++axis)
				assert(std::abs((grid.points[left][axis]+grid.points[right][axis])/2
					-centers[cell][axis]) < 1e-14);
			if (local > 0) assert(grid.connectivity[2*cell-1] == left);
		}
	}
	assert(degree[1] == 3); // The real junction, not a cell midpoint.
	assert(degree[0] == 1);
	assert(degree[2] == (cycle ? 2 : 1));
	assert(degree[3] == (cycle ? 2 : 1));
	for (std::size_t point = network.nodes.size(); point < degree.size(); ++point)
		assert(degree[point] == 2);
	iga::OneDFlowState state;
	state.area.assign(network.cells, 1.0);
	state.flow.assign(network.cells, 2.0);
	state.pressure.resize(network.cells);
	std::map<std::string, std::vector<double>> derived;
	derived["oxygen"].resize(network.cells);
	for (int cell = 0; cell < network.cells; ++cell) {
		state.pressure[cell] = 100.0+cell;
		derived["oxygen"][cell] = 0.01*cell;
	}
	const auto stem = std::string(cycle ? "cycle_" : "branch_")+std::to_string(subdivisions);
	const auto path = directory/(stem+".vtkhdf");
	auto arrays = iga::OneDVisualizationArrays(state, {}, derived);
	{
		iga::TemporalVtkHdfPointWriter writer(path, grid, false, 4, iga::VtkHdfAssociation::Cells);
		writer.Append(0.0, arrays);
		writer.Append(0.1, arrays);
		writer.Close();
	}
	{
		iga::TemporalVtkHdfPointWriter writer(path, grid, true, 4, iga::VtkHdfAssociation::Cells);
		for (auto& array : arrays)
			if (array.name == "pressure")
				for (auto& value : array.values) value += 10.0;
		writer.Append(0.1, arrays); // Restart replaces the final stored step.
		writer.Append(0.2, arrays);
		writer.Close();
	}
	bool rejected = false;
	try { iga::TemporalVtkHdfPointWriter wrong_association(path, grid, true); }
	catch (const std::runtime_error&) { rejected = true; }
	assert(rejected);
	iga::WriteOneDVtp(directory/(stem+".vtp"), network, state, {}, derived);
}

int main()
{
	const auto directory = fs::temp_directory_path()/"tubularflowiga-one-d-output-test";
	fs::remove_all(directory);
	fs::create_directories(directory);
	for (const int subdivisions : {1, 3})
		for (const bool cycle : {false, true}) CheckTopology(directory, subdivisions, cycle);
	const auto swc = directory/"network.swc";
	{
		std::ofstream output(swc);
		output << "1 2 0 0 0 0.001 -1\n"
			<< "2 2 0.01 0 0 0.001 1\n"
			<< "3 2 0.02 0 0 0.001 2\n";
	}
	iga::OneDFlowSystemDefinition flow;
	flow.name = "flow";
	flow.model = iga::OneDFlowModel::Compliant;
	flow.scheme = iga::OneDFlowScheme::ImplicitPetsc;
	flow.formulation = iga::OneDImplicitFormulation::ImplicitPde;
	const auto network = iga::ReadOneDNetwork(swc, 1.0, 2, 0.004);
	iga::OneDFlowState state;
	state.area.assign(static_cast<std::size_t>(network.cells), 3.0e-6);
	state.flow.assign(static_cast<std::size_t>(network.cells), 1.0e-9);
	state.pressure.assign(static_cast<std::size_t>(network.cells), 100.0);
	state.node_pressure.assign(network.nodes.size(), 100.0);
	state.segment_flow.assign(network.segments.size(), 1.0e-9);
	state.inlet_flow = 1.0e-9;
	const std::vector<iga::OneDTransportState> transports;
	const std::map<std::string, std::vector<double>> derived;
	const auto vtkhdf_directory = directory/"vtkhdf";
	{
		iga::OneDOutputWriter writer(vtkhdf_directory, network, flow);
		writer.Write(0, 0.0, state, transports, derived);
		state.completed_step = 1;
		state.physical_time = 0.1;
		state.pressure.front() = 110.0;
		writer.Write(1, 0.1, state, transports, derived);
		writer.Finish(state, 0.0, 0.0, 0.0);
	}
	const auto vtkhdf = vtkhdf_directory/"profile_1d.vtkhdf";
	assert(fs::is_regular_file(vtkhdf));
	assert(!fs::exists(vtkhdf_directory/"profile_1d.pvd"));
	const auto file = H5Fopen(vtkhdf.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	assert(file >= 0);
	const auto steps = H5Gopen2(file, "/VTKHDF/Steps", H5P_DEFAULT);
	assert(steps >= 0);
	const auto count_attribute = H5Aopen(steps, "NSteps", H5P_DEFAULT);
	assert(count_attribute >= 0);
	std::int64_t count = 0;
	assert(H5Aread(count_attribute, H5T_NATIVE_INT64, &count) >= 0);
	assert(count == 2);
	H5Aclose(count_attribute);
	H5Gclose(steps);
	const auto pressure = H5Dopen2(file, "/VTKHDF/CellData/pressure", H5P_DEFAULT);
	assert(pressure >= 0);
	const auto pressure_space = H5Dget_space(pressure);
	assert(pressure_space >= 0);
	std::array<hsize_t, 1> dimensions{};
	assert(H5Sget_simple_extent_dims(pressure_space, dimensions.data(), nullptr) >= 0);
	assert(dimensions[0] == 2*static_cast<hsize_t>(network.cells));
	H5Sclose(pressure_space);
	H5Dclose(pressure);
	H5Fclose(file);

	// A single 0D edge must remain a full line, not an isolated midpoint.
	auto single = network;
	single.nodes.resize(2);
	single.segments.resize(1);
	single.segments[0].cells = 1;
	single.cells = 1;
	const auto single_grid = iga::OneDVtkHdfGrid(single);
	assert(single_grid.points.size() == 2 && single_grid.types.size() == 1);
	assert((single_grid.connectivity == std::vector<std::int64_t>{0, 1}));
	flow.model = iga::OneDFlowModel::Lumped;
	flow.formulation = iga::OneDImplicitFormulation::TransientRc;
	{
		iga::OneDOutputWriter writer(directory/"zero_d", network, flow);
		writer.Write(1, 0.1, state, transports, derived);
		writer.Finish(state, 0.0, 0.0, 0.0);
	}
	assert(fs::is_regular_file(directory/"zero_d/profile_0d.vtkhdf"));
	flow.model = iga::OneDFlowModel::Compliant;
	flow.formulation = iga::OneDImplicitFormulation::ImplicitPde;

	// Default point-associated callers remain compatible, including old files
	// written before the explicit DataAssociation metadata was introduced.
	const auto point_path = directory/"point_compatibility.vtkhdf";
	{
		iga::TemporalVtkHdfPointWriter writer(point_path, single_grid);
		writer.Append(0.0, {{"pressure", 1, {10.0, 20.0}}});
		writer.Close();
	}
	{
		const auto old = H5Fopen(point_path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
		const auto metadata = H5Gopen2(old, "/TubularFlowIGA", H5P_DEFAULT);
		assert(H5Adelete(metadata, "DataAssociation") >= 0);
		H5Gclose(metadata);
		H5Fclose(old);
		iga::TemporalVtkHdfPointWriter writer(point_path, single_grid, true);
		writer.Append(0.1, {{"pressure", 1, {30.0, 40.0}}});
		writer.Close();
	}

	const auto legacy_directory = directory/"vtp";
	{
		iga::OneDOutputWriter writer(legacy_directory, network, flow,
			iga::OneDVisualizationFormat::Vtp);
		writer.Write(1, 0.1, state, transports, derived);
		writer.Finish(state, 0.0, 0.0, 0.0);
	}
	assert(fs::is_regular_file(legacy_directory/"profile_1d.pvd"));
	assert(fs::is_regular_file(legacy_directory/"profile_1d_000001.vtp"));
	assert(!fs::exists(legacy_directory/"profile_1d.vtkhdf"));
	if (std::getenv("TUBULARFLOWIGA_KEEP_TEST_OUTPUT") == nullptr)
		fs::remove_all(directory);
}
