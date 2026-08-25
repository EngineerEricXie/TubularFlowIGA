#include "FlowCheckpoint.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

int main()
{
	iga::FlowCheckpointMetadata source;
	source.nodes = 123;
	source.completed_step = 7;
	source.physical_time = 0.7;
	source.dt = 0.1;
	source.density = 1060.0;
	source.viscosity = 0.0035;
	source.state_file = "flow.state";
	source.outlets.push_back({2, "windkessel_rcr", 3.0, 7.0, 5.0});
	const auto parsed = iga::ParseFlowCheckpointMetadata(
		iga::SerializeFlowCheckpointMetadata(source));
	assert(parsed.nodes == source.nodes);
	assert(parsed.completed_step == source.completed_step);
	assert(std::abs(parsed.physical_time-source.physical_time) < 1e-14);
	assert(parsed.state_file == source.state_file);
	assert(parsed.outlets.size() == 1);
	assert(parsed.outlets[0].label == 2);
	assert(parsed.outlets[0].capacitor_pressure == 5.0);
	iga::ValidateFlowCheckpoint(parsed, 123, 10, 0.1, 1060.0, 0.0035);
	assert(iga::TimeIndexedPath("velocity.txt", 12) == "velocity.step000012.txt");
	assert(iga::TimeIndexedPath("velocity", 3) == "velocity.step000003");
	bool rejected = false;
	try {
		iga::ValidateFlowCheckpoint(parsed, 124, 10, 0.1, 1060.0, 0.0035);
	} catch (const std::runtime_error&) {
		rejected = true;
	}
	assert(rejected);
	std::cout << "flow checkpoint metadata tests passed\n";
}
