#include "OneDOutput.hpp"

#include <hdf5.h>

#include <cassert>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <vector>

namespace fs = std::filesystem;

int main()
{
	const auto directory = fs::temp_directory_path()/"tubularflowiga-one-d-output-test";
	fs::remove_all(directory);
	fs::create_directories(directory);
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
	const auto pressure = H5Dopen2(file, "/VTKHDF/PointData/pressure", H5P_DEFAULT);
	assert(pressure >= 0);
	const auto pressure_space = H5Dget_space(pressure);
	assert(pressure_space >= 0);
	std::array<hsize_t, 1> dimensions{};
	assert(H5Sget_simple_extent_dims(pressure_space, dimensions.data(), nullptr) >= 0);
	assert(dimensions[0] == 2*static_cast<hsize_t>(network.cells));
	H5Sclose(pressure_space);
	H5Dclose(pressure);
	H5Fclose(file);

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
