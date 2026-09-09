#include "GraphAcceptedHistory.hpp"
#include "CheckpointRecordStream.hpp"
#include <iostream>

namespace {
int rejected = 0;
void Require(bool value)
{
	if (!value) throw std::runtime_error("history checkpoint test failed");
}
template<class Work> void Reject(Work work)
{
	try { work(); } catch (const std::exception&) { ++rejected; return; }
	throw std::runtime_error("expected history rejection");
}
template<class Record> struct RepeatedRecords {
	const Record& record;
	std::size_t count;
	struct Iterator {
		const Record& record; std::size_t index;
		const Record& operator*() const
		{
			return record;
		}
		Iterator& operator++()
		{
			++index; return *this;
		}
		bool operator!=(const Iterator& other) const
		{
			return index != other.index;
		}
	};
	Iterator begin() const
	{
		return {record, 0};
	}
	Iterator end() const
	{
		return {record, count};
	}
};
}
int main()
{
	try {
		iga::GraphAcceptedFlowStep flow; flow.step = 3; flow.time_s = 0.3; flow.iterations = 2;
		flow.three_d_mass_m3_s = -0.0; flow.external_outward_flow_m3_s = -2;
		flow.three_d_balance = {{"body", {1e-9, 2}}}; flow.zero_d_states = {{"source", {2.5}}};
		flow.zero_d_accounting = {{"source", {1, 2, 3, 4, 5, -6}}};
		iga::PortState port; port.time_s = 0.3; port.area_m2 = 2; port.outward_flow_m3_s = -1;
		port.mean_pressure_pa = 4; port.mean_normal_traction_pa = -4; port.total_pressure_pa = 5;
		port.concentration = {{"tracer", 2.5}}; port.outward_species_flux = {{"tracer", -2.5}};
		flow.result.accepted_ports = {{{"source", "outlet"}, port}};
		for (int i = 0; i < 2; ++i) {
			iga::PressureFlowIterationState iteration; iteration.iteration = i; iteration.applied_relaxation = 0.5; iteration.converged = i == 1;
			iteration.edges.push_back({"edge", 1, 2, 3, 4, 5, 6, 7, 8}); flow.result.iterations.push_back(iteration);
		}
		iga::GraphAcceptedSpeciesStep species; species.step = 3; species.time_s = 0.3; species.hydraulic_iterations = 2;
		species.three_d_balance = flow.three_d_balance; species.transport_ports = flow.result.accepted_ports;
		species.result.hydraulic_iterations = flow.result.iterations; species.result.accepted_ports = flow.result.accepted_ports;
		species.result.transport_ports = flow.result.accepted_ports; species.result.transport_domain_order = {"source", "body"};
		species.result.donor_ownership = {{{"edge", "tracer"}, iga::SpeciesDonor::Second}};
		species.result.edge_amounts = {{"edge", "tracer", {"source", "outlet"}, {"body", "inlet"}, iga::SpeciesDonor::Second, -1, 1, 0, 0}};
		species.result.domain_balances = {{"body", "tracer", {1, 2, 3, {{"inlet", -1}, {"outlet", 3}}, 0}, 0, 0}};
		species.result.global_balances = {{"tracer", {"tracer", 1, 2, 3, 2, 6, 0, 0}}};
		const auto f = iga::SerializeGraphHistoryRecord(flow), s = iga::SerializeGraphHistoryRecord(species);
		Require(iga::SerializeGraphHistoryRecord(iga::ParseGraphFlowHistoryRecord(f)) == f);
		Require(iga::SerializeGraphHistoryRecord(iga::ParseGraphSpeciesHistoryRecord(s)) == s);
		for (std::size_t i = 0; i < f.size(); ++i) Reject([&] { iga::ParseGraphFlowHistoryRecord(std::string_view(f).substr(0, i)); });
		for (std::size_t i = 0; i < s.size(); ++i) Reject([&] { iga::ParseGraphSpeciesHistoryRecord(std::string_view(s).substr(0, i)); });
		Reject([&] { iga::ParseGraphFlowHistoryRecord(f+"x"); }); Reject([&] { iga::ParseGraphSpeciesHistoryRecord(s+"x"); });
		Reject([&] { iga::ParseGraphSpeciesHistoryRecord(f); }); Reject([&] { iga::ParseGraphFlowHistoryRecord(s); });
		Reject([&] { auto value = flow; value.iterations = 3; iga::ParseGraphFlowHistoryRecord(iga::SerializeGraphHistoryRecord(value)); });
		Reject([&] { auto value = species; value.hydraulic_iterations = 0; iga::SerializeGraphHistoryRecord(value); });
		Reject([&] { auto value = species; value.result.donor_ownership.begin()->second = static_cast<iga::SpeciesDonor>(123); iga::SerializeGraphHistoryRecord(value); });
		Reject([&] { auto value = flow; value.zero_d_states.clear(); iga::SerializeGraphHistoryRecord(value); });
		std::vector<iga::GraphAcceptedFlowStep> records{flow, flow}; std::string image;
		iga::WriteCheckpointRecords(records, [&](const void* p, std::size_t n) { image.append(static_cast<const char*>(p), n); });
		Require(image.size() == iga::CheckpointRecordBytes(records));
		for (std::size_t i = 0; i < image.size(); ++i) Reject([&] {
			iga::CheckpointRecordReader reader(2, [&](std::uint64_t, std::string_view bytes) { iga::ParseGraphFlowHistoryRecord(bytes); });
			reader.Consume(image.data(), i); reader.Finish();
		});
		std::uint64_t decoded = 0;
		iga::CheckpointRecordReader fragmented(2, [&](std::uint64_t index, std::string_view bytes) { Require(index == decoded++ && bytes == f); });
		for (const char& byte : image) fragmented.Consume(&byte, 1);
		fragmented.Finish(); Require(decoded == 2); Reject([&] { fragmented.Consume("x", 1); });
		for (std::uint64_t length : {std::uint64_t{0}, std::uint64_t{iga::checkpoint_metadata::maximum_bytes+1}, UINT64_MAX}) {
			iga::CheckpointRecordReader reader(1, [](std::uint64_t, std::string_view) { throw std::runtime_error("unexpected consumer"); });
			char bytes[8]; for (unsigned i = 0; i < 8; ++i) bytes[i] = static_cast<char>(length >> (8*i));
			Reject([&] { reader.Consume(bytes, 8); }); Reject([&] { reader.Consume(image.data(), image.size()); });
		}
		const RepeatedRecords<iga::GraphAcceptedSpeciesStep> large{species, 16000}; std::size_t maximum = 0; std::uint64_t consumed = 0;
		iga::CheckpointRecordReader reader(large.count, [&](std::uint64_t, std::string_view bytes) {
			Require(iga::SerializeGraphHistoryRecord(iga::ParseGraphSpeciesHistoryRecord(bytes)) == s); ++consumed;
		});
		iga::WriteCheckpointRecords(large, [&](const void* p, std::size_t n) { maximum = std::max(maximum, n); reader.Consume(p, n); }); reader.Finish();
		const auto large_bytes = iga::CheckpointRecordBytes(large); Require(large_bytes > iga::checkpoint_metadata::maximum_bytes && maximum <= 65536 && consumed == large.count);
		std::cout << "graph_checkpoint_history status=passed rejections=" << rejected << " large_bytes=" << large_bytes << " maximum_chunk=" << maximum << '\n';
	} catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
