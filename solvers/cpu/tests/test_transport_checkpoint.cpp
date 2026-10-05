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
	const auto schema_one = iga::SerializeTransportCheckpointMetadata(source);
	assert(schema_one.find("time_integration") == std::string::npos);
	assert(parsed.time_integration == "backward_euler" && parsed.history_file.empty());
	auto bdf2 = source;
	bdf2.schema_version = 2;
	bdf2.time_integration = "bdf2";
	bdf2.history_file = "transport.history";
	const auto parsed_bdf2 = iga::ParseTransportCheckpointMetadata(
		iga::SerializeTransportCheckpointMetadata(bdf2));
	assert(parsed_bdf2.schema_version == 2 && parsed_bdf2.time_integration == "bdf2");
	assert(parsed_bdf2.history_file == bdf2.history_file);
	auto rejects = [](iga::TransportCheckpointMetadata value) {
		try {
			(void)iga::ParseTransportCheckpointMetadata(
				iga::SerializeTransportCheckpointMetadata(value));
		} catch (const std::runtime_error&) {
			return true;
		}
		return false;
	};
	auto missing_history = bdf2;
	missing_history.history_file.clear();
	assert(rejects(missing_history));
	auto first_step = bdf2;
	first_step.completed_step = 0;
	first_step.physical_time = 0.0;
	assert(rejects(first_step));
	first_step.history_file.clear();
	assert(!rejects(first_step));
	auto backward_schema_two = bdf2;
	backward_schema_two.time_integration = "backward_euler";
	assert(rejects(backward_schema_two));
	std::cout << "transport checkpoint metadata tests passed\n";
}
