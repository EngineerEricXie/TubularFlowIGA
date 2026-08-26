#include "VcaCheckpoint.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

int main()
{
	iga::VcaCheckpointMetadata metadata;
	metadata.completed_step = 3;
	metadata.physical_time = 0.3;
	metadata.dt = 0.1;
	metadata.fields = {"oxygen", "glucose"};
	metadata.transport_state_file = "checkpoint.vca_transport.state";
	metadata.reservoir.volume_m3 = 1.0e-6;
	metadata.reservoir.temperature_c = 37.0;
	metadata.reservoir.hematocrit_percent = 35.0;
	metadata.reservoir.species = {{"oxygen", 0.2}, {"glucose", 5.0}};
	const auto parsed = iga::ParseVcaCheckpointMetadata(
		iga::SerializeVcaCheckpointMetadata(metadata));
	assert(parsed.completed_step == 3);
	assert(std::abs(parsed.physical_time-0.3) < 1.0e-14);
	assert(parsed.fields == metadata.fields);
	assert(std::abs(parsed.reservoir.species.at("glucose")-5.0) < 1.0e-14);
	iga::ValidateVcaCheckpoint(parsed, 3, 0.1, metadata.fields);
	bool rejected = false;
	try {
		iga::ValidateVcaCheckpoint(parsed, 2, 0.1, metadata.fields);
	} catch (const std::runtime_error&) {
		rejected = true;
	}
	assert(rejected);
	std::cout << "VCA checkpoint metadata tests passed\n";
}
