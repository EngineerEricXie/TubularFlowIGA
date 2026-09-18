#ifndef IGA_ONE_D_OUTPUT_HPP
#define IGA_ONE_D_OUTPUT_HPP

#include "OneDTransport.hpp"
#include "TemporalVtkHdf.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <vector>

namespace iga {

enum class OneDVisualizationFormat { VtkHdf, Vtp };

inline OneDVisualizationFormat ParseOneDVisualizationFormat(const std::string& value)
{
	if (value == "auto" || value == "vtkhdf") return OneDVisualizationFormat::VtkHdf;
	if (value == "vtp") return OneDVisualizationFormat::Vtp;
	throw std::runtime_error("--visualization-format must be auto, vtkhdf, or vtp");
}

inline const char* OneDVisualizationFormatName(OneDVisualizationFormat value)
{
	return value == OneDVisualizationFormat::VtkHdf ? "vtkhdf" : "vtp";
}

inline const char* OneDModelName(OneDFlowModel value)
{
	if (value == OneDFlowModel::Rigid) return "rigid";
	if (value == OneDFlowModel::Compliant) return "compliant";
	return "lumped";
}

inline const char* OneDSchemeName(OneDFlowScheme value)
{
	if (value == OneDFlowScheme::SteadyPoiseuille) return "steady_poiseuille";
	if (value == OneDFlowScheme::RigidInertance) return "rigid_inertance";
	if (value == OneDFlowScheme::Petsc) return "petsc";
	if (value == OneDFlowScheme::ExplicitRusanov) return "explicit_rusanov";
	return "implicit_petsc";
}

inline const char* OneDFormulationName(OneDImplicitFormulation value)
{
	if (value == OneDImplicitFormulation::SteadyR) return "steady_r";
	if (value == OneDImplicitFormulation::TransientRc) return "transient_rc";
	if (value == OneDImplicitFormulation::TransientRlc) return "transient_rlc";
	if (value == OneDImplicitFormulation::NonlinearRlc) return "nonlinear_rlc";
	return "implicit_1d_pde";
}

inline const char* OneDOutletKindName(OneDOutletKind value)
{
	if (value == OneDOutletKind::Pressure) return "pressure";
	if (value == OneDOutletKind::Resistance) return "resistance";
	if (value == OneDOutletKind::WindkesselRc) return "windkessel_rc";
	return "windkessel_rcr";
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

inline std::vector<VtkPointArray> OneDVisualizationArrays(
	const OneDFlowState& flow, const std::vector<OneDTransportState>& transports,
	const std::map<std::string, std::vector<double>>& derived)
{
	std::vector<VtkPointArray> arrays;
	arrays.push_back({"area", 1, flow.area});
	arrays.push_back({"flow_rate", 1, flow.flow});
	arrays.push_back({"pressure", 1, flow.pressure});
	std::vector<double> velocity(flow.flow.size());
	for (std::size_t i = 0; i < velocity.size(); ++i)
		velocity[i] = flow.flow[i]/flow.area[i];
	arrays.push_back({"velocity", 1, std::move(velocity)});
	for (const auto& transport : transports)
		for (const auto& species : transport.species)
			arrays.push_back({species.definition.field, 1, species.concentration});
	for (const auto& field : derived)
		arrays.push_back({field.first, 1, field.second});
	return arrays;
}

inline VtkHdfUnstructuredGrid OneDVtkHdfGrid(const OneDNetwork& network)
{
	constexpr std::uint8_t vtk_vertex = 1;
	constexpr std::uint8_t vtk_line = 3;
	constexpr std::uint8_t vtk_poly_line = 4;
	VtkHdfUnstructuredGrid grid;
	grid.points = OneDCellCenters(network);
	grid.offsets.push_back(0);
	std::vector<std::vector<int>> node_cells(network.nodes.size());
	std::vector<bool> used(grid.points.size(), false);
	for (const auto& segment : network.segments) {
		node_cells[static_cast<std::size_t>(segment.parent)].push_back(segment.cell_offset);
		node_cells[static_cast<std::size_t>(segment.child)].push_back(
			segment.cell_offset+segment.cells-1);
		if (segment.cells <= 1) continue;
		for (int cell = 0; cell < segment.cells; ++cell) {
			const auto point = segment.cell_offset+cell;
			grid.connectivity.push_back(point);
			used[static_cast<std::size_t>(point)] = true;
		}
		grid.offsets.push_back(static_cast<std::int64_t>(grid.connectivity.size()));
		grid.types.push_back(segment.cells == 2 ? vtk_line : vtk_poly_line);
	}
	for (const auto& cells : node_cells)
		for (std::size_t i = 1; i < cells.size(); ++i) {
			grid.connectivity.push_back(cells.front());
			grid.connectivity.push_back(cells[i]);
			used[static_cast<std::size_t>(cells.front())] = true;
			used[static_cast<std::size_t>(cells[i])] = true;
			grid.offsets.push_back(static_cast<std::int64_t>(grid.connectivity.size()));
			grid.types.push_back(vtk_line);
		}
	for (std::size_t point = 0; point < used.size(); ++point)
		if (!used[point]) {
			grid.connectivity.push_back(static_cast<std::int64_t>(point));
			grid.offsets.push_back(static_cast<std::int64_t>(grid.connectivity.size()));
			grid.types.push_back(vtk_vertex);
		}
	std::vector<std::int32_t> segment_id(grid.points.size());
	std::vector<std::int32_t> parent_node_id(grid.points.size());
	std::vector<std::int32_t> child_node_id(grid.points.size());
	std::vector<std::int32_t> segment_cell(grid.points.size());
	std::vector<double> axial_position(grid.points.size());
	std::vector<double> reference_radius(grid.points.size());
	for (const auto& segment : network.segments)
		for (int cell = 0; cell < segment.cells; ++cell) {
			const auto point = static_cast<std::size_t>(segment.cell_offset+cell);
			segment_id[point] = segment.index;
			parent_node_id[point] = network.nodes[static_cast<std::size_t>(segment.parent)].id;
			child_node_id[point] = network.nodes[static_cast<std::size_t>(segment.child)].id;
			segment_cell[point] = cell;
			axial_position[point] = (static_cast<double>(cell)+0.5)
				*segment.length/segment.cells;
			reference_radius[point] = segment.radius0;
		}
	grid.point_int32.push_back({"segment_id", std::move(segment_id)});
	grid.point_int32.push_back({"parent_node_id", std::move(parent_node_id)});
	grid.point_int32.push_back({"child_node_id", std::move(child_node_id)});
	grid.point_int32.push_back({"segment_cell", std::move(segment_cell)});
	grid.point_double.push_back({"axial_position_m", std::move(axial_position)});
	grid.point_double.push_back({"reference_radius_m", std::move(reference_radius)});
	return grid;
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
	std::vector<std::vector<int>> node_cells(network.nodes.size());
	for (const auto& segment : network.segments) {
		node_cells[static_cast<std::size_t>(segment.parent)].push_back(segment.cell_offset);
		node_cells[static_cast<std::size_t>(segment.child)].push_back(
			segment.cell_offset+segment.cells-1);
	}
	int lines = 0;
	for (const auto& segment : network.segments) if (segment.cells > 1) ++lines;
	for (const auto& cells : node_cells)
		if (cells.size() > 1) lines += static_cast<int>(cells.size())-1;
	output << "<?xml version=\"1.0\"?>\n"
		<< "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
		<< "  <PolyData>\n"
		<< "    <Piece NumberOfPoints=\"" << points.size() << "\" NumberOfVerts=\"0\" NumberOfLines=\""
		<< lines << "\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n"
		<< "      <PointData>\n";
	for (const auto& array : OneDVisualizationArrays(flow, transports, derived)) {
		if (array.components != 1 || array.values.size() != points.size())
			throw std::runtime_error("1d VTP field size mismatch");
		output << "        <DataArray type=\"Float64\" Name=\"" << OneDEscapeXml(array.name)
			<< "\" format=\"ascii\">\n          ";
		for (const double value : array.values) output << value << ' ';
		output << "\n        </DataArray>\n";
	}
	output << "      </PointData>\n"
		<< "      <Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n        ";
	for (const auto& point : points) output << point[0] << ' ' << point[1] << ' ' << point[2] << ' ';
	output << "\n      </DataArray></Points>\n"
		<< "      <Lines>\n"
		<< "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n          ";
	for (const auto& segment : network.segments) if (segment.cells > 1)
		for (int cell = 0; cell < segment.cells; ++cell) output << segment.cell_offset+cell << ' ';
	for (const auto& cells : node_cells)
		for (std::size_t i = 1; i < cells.size(); ++i)
			output << cells.front() << ' ' << cells[i] << ' ';
	output << "\n        </DataArray>\n"
		<< "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n          ";
	int offset = 0;
	for (const auto& segment : network.segments) if (segment.cells > 1) {
		offset += segment.cells;
		output << offset << ' ';
	}
	for (const auto& cells : node_cells)
		for (std::size_t i = 1; i < cells.size(); ++i) {
			offset += 2;
			output << offset << ' ';
		}
	output << "\n        </DataArray>\n"
		<< "      </Lines>\n"
		<< "    </Piece>\n  </PolyData>\n</VTKFile>\n";
	output.close();
	if (!output) throw std::runtime_error("cannot write 1d VTP output: " + path.string());
}

inline void RemoveSupersededOneDVisualization(const std::filesystem::path& directory,
	const std::string& dimension, OneDVisualizationFormat format)
{
	if (!std::filesystem::exists(directory)) return;
	const std::string stem = "profile_"+dimension;
	if (format == OneDVisualizationFormat::Vtp) {
		std::filesystem::remove(directory/(stem+".vtkhdf"));
		return;
	}
	std::filesystem::remove(directory/(stem+".pvd"));
	const std::string prefix = stem+"_";
	for (const auto& entry : std::filesystem::directory_iterator(directory)) {
		if (!entry.is_regular_file()) continue;
		const auto name = entry.path().filename().string();
		if (entry.path().extension() == ".vtp" && name.rfind(prefix, 0) == 0)
			std::filesystem::remove(entry.path());
	}
}

class OneDOutputWriter {
public:
	OneDOutputWriter(std::filesystem::path directory, const OneDNetwork& network,
		const OneDFlowSystemDefinition& system,
		OneDVisualizationFormat visualization_format = OneDVisualizationFormat::VtkHdf,
		bool resume_visualization = false)
		: directory_(std::move(directory)), network_(network), system_(system),
		  visualization_format_(visualization_format)
	{
		std::filesystem::create_directories(directory_);
		const std::string dimension = NetworkFlowPhysicalDimension(system_);
		if (!resume_visualization)
			RemoveSupersededOneDVisualization(directory_, dimension, visualization_format_);
		if (visualization_format_ == OneDVisualizationFormat::VtkHdf)
			vtkhdf_ = std::make_unique<TemporalVtkHdfPointWriter>(
				directory_/("profile_"+dimension+".vtkhdf"),
				OneDVtkHdfGrid(network_), resume_visualization);
		flow_.open(directory_/"flow_timeseries.csv");
		outlet_.open(directory_/"outlet_timeseries.csv");
		node_.open(directory_/"node_timeseries.csv");
		branch_.open(directory_/"branch_timeseries.csv");
		profile_.open(directory_/("profile_"+dimension+".csv"));
		species_.open(directory_/("species_profile_"+dimension+".csv"));
		derived_.open(directory_/("derived_profile_"+dimension+".csv"));
		if (!flow_ || !outlet_ || !node_ || !branch_ || !profile_ || !species_ || !derived_)
			throw std::runtime_error("cannot create network output files in " + directory_.string());
		flow_ << "time,inlet_flow_rate,outlet_flow_rate_sum,storage_rate,relative_continuity_residual,distal_outlet_flow_rate_sum,outlet_capacitor_storage_rate,total_circuit_relative_residual,pressure_drop,min_area,max_area\n";
		outlet_ << "time,node_id,kind,inlet_flow,terminal_pressure,capacitor_pressure,distal_flow,capacitor_storage_rate,proximal_resistance,distal_resistance,capacitance,reference_pressure\n";
		node_ << "time,node_id,pressure,is_root,is_outlet\n";
		branch_ << "time,parent_id,child_id,flow_rate,pressure_parent,pressure_child,segment_resistance\n";
		profile_ << "time,parent_id,child_id,cell,x,area,flow_rate,pressure,velocity\n";
		species_ << "time,parent_id,child_id,cell,x,species,concentration,species_flux\n";
		derived_ << "time,parent_id,child_id,cell,x,field,value\n";
		flow_ << std::setprecision(17);
		outlet_ << std::setprecision(17);
		node_ << std::setprecision(17);
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
		const double selected_storage_rate = state.has_conservation_diagnostic
			? state.storage_rate : storage_rate;
		const double continuity = inlet-outlet-selected_storage_rate;
		const double relative_continuity = state.has_conservation_diagnostic
			? state.relative_continuity_residual : std::abs(continuity)
				/std::max({std::abs(inlet), std::abs(outlet),
					std::abs(selected_storage_rate), 1.0e-30});
		double distal_outlet = 0.0;
		double outlet_capacitor_storage = 0.0;
		for (const auto& item : state.outlets) {
			double distal_flow = item.flow;
			if (item.kind == OneDOutletKind::WindkesselRc
				|| item.kind == OneDOutletKind::WindkesselRcr)
				distal_flow = (item.capacitor_pressure-item.reference_pressure)
					/item.distal_resistance;
			distal_outlet += distal_flow;
			outlet_capacitor_storage += item.flow-distal_flow;
		}
		if (system_.model == OneDFlowModel::Lumped)
			for (std::size_t node = 0; node < network_.nodes.size(); ++node)
				node_ << time << ',' << network_.nodes[node].id << ','
					<< state.node_pressure[node] << ','
					<< (static_cast<int>(node) == network_.root ? 1 : 0) << ','
					<< (network_.nodes[node].children.empty() ? 1 : 0) << '\n';
		const double circuit_residual = inlet-distal_outlet-selected_storage_rate
			-outlet_capacitor_storage;
		const double relative_circuit_residual = std::abs(circuit_residual)
			/std::max({std::abs(inlet), std::abs(distal_outlet),
				std::abs(selected_storage_rate), std::abs(outlet_capacitor_storage), 1.0e-30});
		maximum_relative_continuity_residual_ = std::max(
			maximum_relative_continuity_residual_, relative_continuity);
		maximum_relative_circuit_residual_ = std::max(
			maximum_relative_circuit_residual_, relative_circuit_residual);
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
		flow_ << time << ',' << inlet << ',' << outlet << ',' << selected_storage_rate << ','
			<< relative_continuity << ',' << distal_outlet << ',' << outlet_capacitor_storage
			<< ',' << relative_circuit_residual << ',' << root_pressure-outlet_pressure << ','
			<< *area_bounds.first << ',' << *area_bounds.second << '\n';
		for (const auto& item : state.outlets) {
			double distal_flow = item.flow;
			double capacitor_storage_rate = 0.0;
			if (item.kind == OneDOutletKind::WindkesselRc
				|| item.kind == OneDOutletKind::WindkesselRcr) {
				distal_flow = (item.capacitor_pressure-item.reference_pressure)
					/item.distal_resistance;
				capacitor_storage_rate = item.flow-distal_flow;
			}
			outlet_ << time << ','
				<< network_.nodes[static_cast<std::size_t>(item.node)].id << ','
				<< OneDOutletKindName(item.kind) << ',' << item.flow << ',' << item.pressure
				<< ',' << item.capacitor_pressure << ',' << distal_flow << ','
				<< capacitor_storage_rate << ',' << item.proximal_resistance << ','
				<< item.distal_resistance << ',' << item.capacitance << ','
				<< item.reference_pressure << '\n';
		}
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
				<< parent_pressure << ',' << child_pressure << ','
				<< ((system_.formulation == OneDImplicitFormulation::SteadyR
					|| system_.formulation == OneDImplicitFormulation::TransientRc)
					? OneDLumpedSegmentResistance(network_, system_, segment)
					: segment.resistance) << '\n';
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
		if (vtkhdf_)
			vtkhdf_->Append(time, OneDVisualizationArrays(state, transports, derived));
		else {
			std::ostringstream name;
			name.exceptions(std::ios::badbit | std::ios::failbit);
			const std::string prefix = std::string("profile_")
				+NetworkFlowPhysicalDimension(system_)+"_";
			name << prefix << std::setw(6) << std::setfill('0') << step << ".vtp";
			WriteOneDVtp(directory_/name.str(), network_, state, transports, derived);
			vtp_.push_back({time, name.str()});
		}
		for (auto* stream : {&flow_, &outlet_, &node_, &branch_, &profile_, &species_, &derived_}) {
			stream->flush();
			if (!*stream)
				throw std::runtime_error("cannot write network CSV output in " + directory_.string());
		}
	}

	void Finish(const OneDFlowState& state, double setup_seconds,
		double solve_seconds, double output_seconds)
	{
		for (auto* stream : {&flow_, &outlet_, &node_, &branch_, &profile_, &species_, &derived_}) {
			stream->close();
			if (!*stream)
				throw std::runtime_error("cannot close network CSV output in " + directory_.string());
		}
		if (vtkhdf_) vtkhdf_->Close();
		if (visualization_format_ == OneDVisualizationFormat::Vtp) {
			const std::string collection = std::string("profile_")
				+NetworkFlowPhysicalDimension(system_)+".pvd";
			std::ofstream pvd(directory_/collection);
			pvd << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n  <Collection>\n";
			for (const auto& item : vtp_)
				pvd << "    <DataSet timestep=\"" << std::setprecision(17) << item.first
					<< "\" group=\"\" part=\"0\" file=\"" << OneDEscapeXml(item.second) << "\"/>\n";
			pvd << "  </Collection>\n</VTKFile>\n";
			pvd.close();
			if (!pvd)
				throw std::runtime_error("cannot write network PVD output in " + directory_.string());
		}
		std::ofstream summary(directory_/"summary.json");
		struct rusage usage {};
		getrusage(RUSAGE_SELF, &usage);
		summary << std::setprecision(17)
			<< "{\n  \"schema_version\": 1,\n"
			<< "  \"dimension\": \"" << NetworkFlowPhysicalDimension(system_) << "\",\n"
			<< "  \"system\": \"" << OneDEscapeJson(system_.name) << "\",\n"
			<< "  \"model\": \"" << OneDModelName(system_.model) << "\",\n"
			<< "  \"scheme\": \"" << OneDSchemeName(system_.scheme) << "\",\n"
			<< "  \"formulation\": \"" << OneDFormulationName(system_.formulation) << "\",\n"
			<< "  \"visualization_format\": \""
			<< OneDVisualizationFormatName(visualization_format_) << "\",\n"
			<< "  \"visualization_file\": \"profile_"
			<< NetworkFlowPhysicalDimension(system_)
			<< (visualization_format_ == OneDVisualizationFormat::VtkHdf
				? ".vtkhdf\",\n" : ".pvd\",\n")
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
			<< "  \"maximum_sampled_total_circuit_relative_residual\": "
			<< maximum_relative_circuit_residual_ << ",\n"
			<< "  \"setup_seconds\": " << setup_seconds << ",\n"
			<< "  \"solve_seconds\": " << solve_seconds << ",\n"
			<< "  \"output_seconds\": " << output_seconds << ",\n"
			<< "  \"peak_rss_kib\": " << usage.ru_maxrss << "\n}\n";
		summary.close();
		if (!summary)
			throw std::runtime_error("cannot write network summary in " + directory_.string());
	}

private:
	std::filesystem::path directory_;
	const OneDNetwork& network_;
	const OneDFlowSystemDefinition& system_;
	OneDVisualizationFormat visualization_format_ = OneDVisualizationFormat::VtkHdf;
	std::unique_ptr<TemporalVtkHdfPointWriter> vtkhdf_;
	std::ofstream flow_, outlet_, node_, branch_, profile_, species_, derived_;
	std::vector<std::pair<double, std::string>> vtp_;
	bool has_previous_volume_ = false;
	double previous_volume_ = 0.0;
	double previous_time_ = 0.0;
	double maximum_relative_continuity_residual_ = 0.0;
	double maximum_relative_circuit_residual_ = 0.0;
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
	output.close();
	if (!output) throw std::runtime_error("cannot write physiology_fields.json");
}

} // namespace iga

#endif
