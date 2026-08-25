#include "TransportCheckpoint.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

int main()
{
	iga::TransportCheckpointMetadata source;
	source.nodes = 321;
	source.fields = {"tracer", "oxygen"};
	source.system = "coupled transport";
	source.velocity_source = "flow_snapshots";
	source.completed_step = 4;
	source.physical_time = 0.2;
	source.dt = 0.05;
	source.state_file = "transport.state";
	const auto parsed = iga::ParseTransportCheckpointMetadata(
		iga::SerializeTransportCheckpointMetadata(source));
	assert(parsed.nodes == source.nodes);
	assert(parsed.fields == source.fields);
	assert(parsed.system == source.system);
	assert(parsed.velocity_source == source.velocity_source);
	assert(parsed.completed_step == source.completed_step);
	assert(std::abs(parsed.physical_time-source.physical_time) < 1e-14);
	iga::ValidateTransportCheckpoint(parsed, 321, source.fields,
		source.system, source.velocity_source, 10, 0.05);
	bool rejected = false;
	try {
		iga::ValidateTransportCheckpoint(parsed, 321, {"oxygen", "tracer"},
			source.system, source.velocity_source, 10, 0.05);
	} catch (const std::runtime_error&) {
		rejected = true;
	}
	assert(rejected);
	rejected = false;
	try {
		iga::ValidateTransportCheckpoint(parsed, 321, source.fields,
			source.system, "prescribed", 10, 0.05);
	} catch (const std::runtime_error&) {
		rejected = true;
	}
	assert(rejected);
	std::cout << "transport checkpoint metadata tests passed\n";
}
