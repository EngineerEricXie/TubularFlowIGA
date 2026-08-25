#ifndef IGA_ONE_D_OUTPUT_HPP
#define IGA_ONE_D_OUTPUT_HPP

#include "OneDTransport.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <vector>

namespace iga {

inline const char* OneDModelName(OneDFlowModel value)
{
	return value == OneDFlowModel::Rigid ? "rigid" : "compliant";
}

inline const char* OneDSchemeName(OneDFlowScheme value)
{
	if (value == OneDFlowScheme::SteadyPoiseuille) return "steady_poiseuille";
	if (value == OneDFlowScheme::ExplicitRusanov) return "explicit_rusanov";
	return "implicit_petsc";
}

inline const char* OneDFormulationName(OneDImplicitFormulation value)
{
	if (value == OneDImplicitFormulation::PressureNetwork) return "pressure_network";
	if (value == OneDImplicitFormulation::LinearizedAQ) return "linearized_aq";
	if (value == OneDImplicitFormulation::NonlinearAQ) return "nonlinear_aq";
	return "implicit_1d_pde";
}

inline std::string OneDEscapeJson(const std::string& value)
{
	std::string result;
	for (const char character : value) {
		if (character == '\\' || character == '"') result.push_back('\\');
		result.push_back(character);
	}
	return result;
}

inline std::string OneDEscapeXml(const std::string& value)
{
	std::string result;
	for (const char character : value) {
		if (character == '&') result += "&amp;";
		else if (character == '<') result += "&lt;";
		else if (character == '>') result += "&gt;";
		else if (character == '"') result += "&quot;";
		else result.push_back(character);
	}
	return result;
}

inline std::vector<std::array<double, 3>> OneDCellCenters(const OneDNetwork& network)
{
	std::vector<std::array<double, 3>> points(static_cast<std::size_t>(network.cells));
	for (const auto& segment : network.segments) {
		const auto& first = network.nodes[static_cast<std::size_t>(segment.parent)].position;
		const auto& last = network.nodes[static_cast<std::size_t>(segment.child)].position;
		for (int cell = 0; cell < segment.cells; ++cell) {
			const double fraction = (static_cast<double>(cell)+0.5)/segment.cells;
			auto& point = points[static_cast<std::size_t>(segment.cell_offset+cell)];
			for (int axis = 0; axis < 3; ++axis)
				point[static_cast<std::size_t>(axis)] = (1.0-fraction)*first[static_cast<std::size_t>(axis)]
					+fraction*last[static_cast<std::size_t>(axis)];
		}
	}
	return points;
}

inline void WriteOneDVtp(const std::filesystem::path& path,
	const OneDNetwork& network, const OneDFlowState& flow,
	const std::vector<OneDTransportState>& transports,
	const std::map<std::string, std::vector<double>>& derived)
{
	std::ofstream output(path);
	if (!output) throw std::runtime_error("cannot create 1d VTP output: " + path.string());
	output << std::setprecision(17);
	const auto points = OneDCellCenters(network);
	int lines = 0;
	for (const auto& segment : network.segments) if (segment.cells > 1) ++lines;
	output << "<?xml version=\"1.0\"?>\n"
		<< "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
		<< "  <PolyData>\n"
		<< "    <Piece NumberOfPoints=\"" << points.size() << "\" NumberOfVerts=\"0\" NumberOfLines=\""
		<< lines << "\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n"
		<< "      <PointData>\n";
	auto array = [&](const std::string& name, const std::vector<double>& values) {
		if (values.size() != points.size()) throw std::runtime_error("1d VTP field size mismatch");
		output << "        <DataArray type=\"Float64\" Name=\"" << OneDEscapeXml(name)
			<< "\" format=\"ascii\">\n          ";
		for (const double value : values) output << value << ' ';
		output << "\n        </DataArray>\n";
	};
	array("area", flow.area);
	array("flow_rate", flow.flow);
	array("pressure", flow.pressure);
	std::vector<double> velocity(flow.flow.size());
	for (std::size_t i = 0; i < velocity.size(); ++i) velocity[i] = flow.flow[i]/flow.area[i];
	array("velocity", velocity);
	for (const auto& transport : transports)
		for (const auto& species : transport.species) array(species.definition.field, species.concentration);
	for (const auto& field : derived) array(field.first, field.second);
	output << "      </PointData>\n"
		<< "      <Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n        ";
	for (const auto& point : points) output << point[0] << ' ' << point[1] << ' ' << point[2] << ' ';
	output << "\n      </DataArray></Points>\n"
		<< "      <Lines>\n"
		<< "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n          ";
	for (const auto& segment : network.segments) if (segment.cells > 1)
		for (int cell = 0; cell < segment.cells; ++cell) output << segment.cell_offset+cell << ' ';
	output << "\n        </DataArray>\n"
		<< "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n          ";
	int offset = 0;
	for (const auto& segment : network.segments) if (segment.cells > 1) {
		offset += segment.cells;
		output << offset << ' ';
	}
	output << "\n        </DataArray>\n"
		<< "      </Lines>\n"
		<< "    </Piece>\n  </PolyData>\n</VTKFile>\n";
	if (!output) throw std::runtime_error("cannot write 1d VTP output: " + path.string());
}

class OneDOutputWriter {
public:
	OneDOutputWriter(std::filesystem::path directory, const OneDNetwork& network,
		const OneDFlowSystemDefinition& system)
		: directory_(std::move(directory)), network_(network), system_(system)
	{
		std::filesystem::create_directories(directory_);
		flow_.open(directory_/"flow_timeseries.csv");
		branch_.open(directory_/"branch_timeseries.csv");
		profile_.open(directory_/"profile_1d.csv");
		species_.open(directory_/"species_profile_1d.csv");
		derived_.open(directory_/"derived_profile_1d.csv");
		if (!flow_ || !branch_ || !profile_ || !species_ || !derived_)
			throw std::runtime_error("cannot create 1d output files in " + directory_.string());
		flow_ << "time,inlet_flow_rate,outlet_flow_rate_sum,storage_rate,relative_continuity_residual,pressure_drop,min_area,max_area\n";
		branch_ << "time,parent_id,child_id,flow_rate,pressure_parent,pressure_child,segment_resistance\n";
		profile_ << "time,parent_id,child_id,cell,x,area,flow_rate,pressure,velocity\n";
		species_ << "time,parent_id,child_id,cell,x,species,concentration,species_flux\n";
		derived_ << "time,parent_id,child_id,cell,x,field,value\n";
		flow_ << std::setprecision(17);
		branch_ << std::setprecision(17);
		profile_ << std::setprecision(17);
		species_ << std::setprecision(17);
		derived_ << std::setprecision(17);
	}

	void Write(int step, double time, const OneDFlowState& state,
		const std::vector<OneDTransportState>& transports,
		const std::map<std::string, std::vector<double>>& derived)
	{
		const double inlet = state.inlet_flow;
		double outlet = 0.0;
		for (const auto& item : state.outlets) outlet += item.flow;
		double volume = 0.0;
		for (const auto& segment : network_.segments) {
			const double dx = segment.length/segment.cells;
			for (int cell = 0; cell < segment.cells; ++cell)
				volume += state.area[static_cast<std::size_t>(segment.cell_offset+cell)]*dx;
		}
		const double storage_rate = has_previous_volume_ && time > previous_time_
			? (volume-previous_volume_)/(time-previous_time_) : inlet-outlet;
		const double continuity = inlet-outlet-storage_rate;
		const double relative_continuity = std::abs(continuity)
			/std::max({std::abs(inlet), std::abs(outlet), std::abs(storage_rate), 1.0e-30});
		maximum_relative_continuity_residual_ = std::max(
			maximum_relative_continuity_residual_, relative_continuity);
		last_outlet_flow_ = outlet;
		previous_volume_ = volume;
		previous_time_ = time;
		has_previous_volume_ = true;
		double outlet_pressure = 0.0;
		if (!state.outlets.empty()) {
			for (const auto& item : state.outlets) outlet_pressure += item.pressure;
			outlet_pressure /= state.outlets.size();
		}
		const double root_pressure = state.node_pressure.empty() ? 0.0 : state.node_pressure[static_cast<std::size_t>(network_.root)];
		const auto area_bounds = std::minmax_element(state.area.begin(), state.area.end());
		flow_ << time << ',' << inlet << ',' << outlet << ',' << storage_rate << ',' << relative_continuity << ','
			<< root_pressure-outlet_pressure << ',' << *area_bounds.first << ',' << *area_bounds.second << '\n';
		for (const auto& segment : network_.segments) {
			const double parent_pressure = state.node_pressure.size() == network_.nodes.size()
				? state.node_pressure[static_cast<std::size_t>(segment.parent)]
				: state.pressure[static_cast<std::size_t>(segment.cell_offset)];
			const double child_pressure = state.node_pressure.size() == network_.nodes.size()
				? state.node_pressure[static_cast<std::size_t>(segment.child)]
				: state.pressure[static_cast<std::size_t>(segment.cell_offset+segment.cells-1)];
			branch_ << time << ',' << network_.nodes[static_cast<std::size_t>(segment.parent)].id << ','
				<< network_.nodes[static_cast<std::size_t>(segment.child)].id << ','
				<< state.segment_flow[static_cast<std::size_t>(segment.index)] << ','
				<< parent_pressure << ',' << child_pressure << ',' << segment.resistance << '\n';
			for (int cell = 0; cell < segment.cells; ++cell) {
				const auto index = static_cast<std::size_t>(segment.cell_offset+cell);
				const double x = (static_cast<double>(cell)+0.5)*segment.length/segment.cells;
				profile_ << time << ',' << network_.nodes[static_cast<std::size_t>(segment.parent)].id << ','
					<< network_.nodes[static_cast<std::size_t>(segment.child)].id << ',' << cell << ',' << x << ','
					<< state.area[index] << ',' << state.flow[index] << ',' << state.pressure[index] << ','
					<< state.flow[index]/state.area[index] << '\n';
				for (const auto& transport : transports)
					for (const auto& field : transport.species)
						species_ << time << ',' << network_.nodes[static_cast<std::size_t>(segment.parent)].id << ','
							<< network_.nodes[static_cast<std::size_t>(segment.child)].id << ',' << cell << ',' << x << ','
							<< field.definition.field << ',' << field.concentration[index] << ','
							<< state.flow[index]*field.concentration[index] << '\n';
				for (const auto& field : derived)
					derived_ << time << ',' << network_.nodes[static_cast<std::size_t>(segment.parent)].id << ','
						<< network_.nodes[static_cast<std::size_t>(segment.child)].id << ',' << cell << ',' << x << ','
						<< field.first << ',' << field.second[index] << '\n';
			}
		}
		std::ostringstream name;
		name << "profile_1d_" << std::setw(6) << std::setfill('0') << step << ".vtp";
		WriteOneDVtp(directory_/name.str(), network_, state, transports, derived);
		vtp_.push_back({time, name.str()});
	}

	void Finish(const OneDFlowState& state, double setup_seconds,
		double solve_seconds, double output_seconds)
	{
		std::ofstream pvd(directory_/"profile_1d.pvd");
		pvd << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n  <Collection>\n";
		for (const auto& item : vtp_)
			pvd << "    <DataSet timestep=\"" << std::setprecision(17) << item.first
				<< "\" group=\"\" part=\"0\" file=\"" << OneDEscapeXml(item.second) << "\"/>\n";
		pvd << "  </Collection>\n</VTKFile>\n";
		std::ofstream summary(directory_/"summary.json");
		struct rusage usage {};
		getrusage(RUSAGE_SELF, &usage);
		summary << std::setprecision(17)
			<< "{\n  \"schema_version\": 1,\n"
			<< "  \"system\": \"" << OneDEscapeJson(system_.name) << "\",\n"
			<< "  \"model\": \"" << OneDModelName(system_.model) << "\",\n"
			<< "  \"scheme\": \"" << OneDSchemeName(system_.scheme) << "\",\n"
			<< "  \"formulation\": \"" << OneDFormulationName(system_.formulation) << "\",\n"
			<< "  \"dynamic_viscosity_pa_s\": " << system_.dynamic_viscosity << ",\n"
			<< "  \"density_kg_m3\": " << system_.density << ",\n"
			<< "  \"nodes\": " << network_.nodes.size() << ",\n"
			<< "  \"segments\": " << network_.segments.size() << ",\n"
			<< "  \"cells\": " << network_.cells << ",\n"
			<< "  \"completed_step\": " << state.completed_step << ",\n"
			<< "  \"physical_time\": " << state.physical_time << ",\n"
			<< "  \"internal_substeps\": " << state.internal_substeps << ",\n"
			<< "  \"converged\": true,\n"
			<< "  \"final_inlet_flow_rate\": " << state.inlet_flow << ",\n"
			<< "  \"final_outlet_flow_rate_sum\": " << last_outlet_flow_ << ",\n"
			<< "  \"maximum_sampled_relative_continuity_residual\": "
			<< maximum_relative_continuity_residual_ << ",\n"
			<< "  \"setup_seconds\": " << setup_seconds << ",\n"
			<< "  \"solve_seconds\": " << solve_seconds << ",\n"
			<< "  \"output_seconds\": " << output_seconds << ",\n"
			<< "  \"peak_rss_kib\": " << usage.ru_maxrss << "\n}\n";
	}

private:
	std::filesystem::path directory_;
	const OneDNetwork& network_;
	const OneDFlowSystemDefinition& system_;
	std::ofstream flow_, branch_, profile_, species_, derived_;
	std::vector<std::pair<double, std::string>> vtp_;
	bool has_previous_volume_ = false;
	double previous_volume_ = 0.0;
	double previous_time_ = 0.0;
	double maximum_relative_continuity_residual_ = 0.0;
	double last_outlet_flow_ = 0.0;
};

inline void WriteOneDPhysiologyManifest(const std::filesystem::path& directory,
	const OneDConfiguration& configuration,
	const std::vector<OneDTransportState>& transports,
	const std::map<std::string, std::vector<double>>& derived)
{
	std::ofstream output(directory/"physiology_fields.json");
	if (!output) throw std::runtime_error("cannot create physiology_fields.json");
	std::set<std::string> solved;
	for (const auto& transport : transports)
		for (const auto& species : transport.species) solved.insert(species.definition.field);
	output << "{\n  \"enabled\": " << (configuration.physiology.enabled ? "true" : "false") << ",\n"
		<< "  \"fields\": {";
	bool first = true;
	for (const auto& name : solved) {
		output << (first ? "\n" : ",\n") << "    \"" << OneDEscapeJson(name) << "\": {\"status\": \"solved\"}";
		first = false;
	}
	for (const auto& name : configuration.physiology.derived_fields) {
		const bool available = derived.count(name) != 0;
		output << (first ? "\n" : ",\n") << "    \"" << OneDEscapeJson(name)
			<< "\": {\"status\": \"" << (available ? "derived" : "skipped") << "\"}";
		first = false;
	}
	if (!first) output << '\n';
	output << "  }\n}\n";
}

} // namespace iga

#endif
