#include "BoundaryFlow.hpp"
#include "CaseInput.hpp"
#include "IgaDatabase.hpp"
#include "VelocitySeries.hpp"
#include "WomersleyReference.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct FlowMetrics {
	std::map<int, double> boundary_flow;
	double net_flow = 0.0;
	double absolute_flow = 0.0;
	double relative_imbalance = 0.0;
	double volume_divergence = 0.0;
	double divergence_theorem_error = 0.0;
	double relative_divergence_theorem_error = 0.0;
};

FlowMetrics MeasureFlow(const std::vector<iga::Element>& elements,
	const std::vector<std::array<double, 3>>& velocity)
{
	FlowMetrics result;
	for (const auto& element : elements) {
		std::vector<std::array<double, 4>> nodal(element.connectivity.size());
		for (std::size_t a = 0; a < element.connectivity.size(); ++a) {
			if (element.connectivity[a] < 0)
				throw std::runtime_error("element connectivity is negative");
			const auto node = static_cast<std::size_t>(element.connectivity[a]);
			if (node >= velocity.size())
				throw std::runtime_error("element connectivity exceeds velocity field");
			nodal[a] = {velocity[node][0], velocity[node][1], velocity[node][2], 0.0};
		}
		for (std::size_t face = 0; face < element.boundary_labels.size(); ++face) {
			const auto label = element.boundary_labels[face];
			if (label < 0) continue;
			result.boundary_flow[label] +=
				iga::IntegrateBoundaryFlow(element, face, nodal);
		}
		result.volume_divergence += iga::IntegrateVolumeDivergence(element, nodal);
	}
	if (result.boundary_flow.empty())
		throw std::runtime_error(
			"database contains no packed boundary faces; repack with iga_pack");
	for (const auto& item : result.boundary_flow) {
		result.net_flow += item.second;
		result.absolute_flow += std::abs(item.second);
	}
	result.relative_imbalance = result.absolute_flow > 0.0
		? 2.0*std::abs(result.net_flow)/result.absolute_flow : 0.0;
	result.divergence_theorem_error = std::abs(
		result.net_flow-result.volume_divergence);
	const auto divergence_scale = std::max(result.absolute_flow,
		std::abs(result.volume_divergence));
	result.relative_divergence_theorem_error = divergence_scale > 0.0
		? result.divergence_theorem_error/divergence_scale : 0.0;
	return result;
}

void PrintFlowMetrics(const FlowMetrics& metrics, const std::string& prefix)
{
	for (const auto& item : metrics.boundary_flow)
		std::cout << prefix << "boundary_flow[" << item.first << "]="
			<< item.second << '\n';
	std::cout << prefix << "net_outward_flow=" << metrics.net_flow
		<< ' ' << prefix << "absolute_boundary_flow=" << metrics.absolute_flow
		<< ' ' << prefix << "relative_mass_imbalance="
		<< metrics.relative_imbalance
		<< ' ' << prefix << "volume_divergence=" << metrics.volume_divergence
		<< ' ' << prefix << "divergence_theorem_error="
		<< metrics.divergence_theorem_error
		<< ' ' << prefix << "relative_divergence_theorem_error="
		<< metrics.relative_divergence_theorem_error << '\n';
}

double RelativeVelocityL2(const std::vector<std::array<double, 3>>& reference,
	const std::vector<std::array<double, 3>>& current)
{
	if (reference.size() != current.size())
		throw std::runtime_error("cycle comparison velocity fields have different sizes");
	double difference_squared = 0.0;
	double reference_squared = 0.0;
	for (std::size_t node = 0; node < reference.size(); ++node)
		for (int component = 0; component < 3; ++component) {
			const auto difference = current[node][component]-reference[node][component];
			difference_squared += difference*difference;
			reference_squared += reference[node][component]*reference[node][component];
		}
	if (reference_squared == 0.0)
		return difference_squared == 0.0 ? 0.0
			: std::numeric_limits<double>::infinity();
	return std::sqrt(difference_squared/reference_squared);
}

double WomersleyVolumeRelativeL2(const std::vector<iga::Element>& elements,
	const std::vector<std::array<double, 3>>& velocity,
	const iga::WomersleyReferenceConfiguration& configuration, double physical_time)
{
	constexpr std::array<double, 4> points{{0.06943184420297371, 0.33000947820757187,
		0.6699905217924281, 0.9305681557970262}};
	constexpr std::array<double, 4> weights{{0.3478548451374539, 0.6521451548625461,
		0.6521451548625461, 0.3478548451374539}};
	double difference_squared = 0.0;
	double reference_squared = 0.0;
	for (const auto& element : elements) {
		std::vector<std::array<double, 3>> nodal(element.connectivity.size());
		for (std::size_t a = 0; a < element.connectivity.size(); ++a) {
			if (element.connectivity[a] < 0
				|| static_cast<std::size_t>(element.connectivity[a]) >= velocity.size())
				throw std::runtime_error(
					"element connectivity exceeds Womersley velocity field");
			nodal[a] = velocity[static_cast<std::size_t>(element.connectivity[a])];
		}
		for (std::size_t qz = 0; qz < 4; ++qz)
			for (std::size_t qy = 0; qy < 4; ++qy)
				for (std::size_t qx = 0; qx < 4; ++qx) {
					const auto basis = iga::EvaluateBoundaryBasis(
						element, points[qx], points[qy], points[qz]);
					std::array<double, 3> numerical{};
					for (std::size_t a = 0; a < nodal.size(); ++a)
						for (int component = 0; component < 3; ++component)
							numerical[component] += basis.value[a]*nodal[a][component];
					const auto reference = iga::WomersleyVelocity(
						configuration, basis.physical_coordinate, physical_time);
					const auto measure = weights[qx]*weights[qy]*weights[qz]
						*basis.determinant;
					for (int component = 0; component < 3; ++component) {
						const auto difference = numerical[component]-reference[component];
						difference_squared += measure*difference*difference;
						reference_squared += measure*reference[component]*reference[component];
					}
				}
	}
	if (reference_squared == 0.0)
		return difference_squared == 0.0 ? 0.0
			: std::numeric_limits<double>::infinity();
	return std::sqrt(difference_squared/reference_squared);
}

double ParsePositiveDouble(const std::string& text, const std::string& name)
{
	std::size_t used = 0;
	double value = 0.0;
	try {
		value = std::stod(text, &used);
	} catch (const std::exception&) {
		throw std::runtime_error(name+" requires a positive finite number");
	}
	if (used != text.size() || !std::isfinite(value) || !(value > 0.0))
		throw std::runtime_error(name+" requires a positive finite number");
	return value;
}

} // namespace

int main(int argc, char** argv)
{
	try {
		if (argc != 3 && argc != 5 && argc != 6 && argc != 7)
			throw std::runtime_error(
				"usage: iga_flow_validate DATABASE.ntiga VELOCITY.txt\n"
				"   or: iga_flow_validate DATABASE.ntiga --compare REFERENCE.txt CURRENT.txt\n"
				"   or: iga_flow_validate DATABASE.ntiga --manifest CASE_DIR MANIFEST.csv [PERIOD]\n"
				"   or: iga_flow_validate DATABASE.ntiga --compare-manifests "
				"REFERENCE_CASE_DIR REFERENCE.csv CURRENT_CASE_DIR CURRENT.csv\n"
				"   or: iga_flow_validate DATABASE.ntiga --womersley "
				"WOMERSLEY.json CURRENT_CASE_DIR CURRENT.csv");
		iga::Database database(argv[1]);
		std::vector<iga::Element> elements;
		elements.reserve(static_cast<std::size_t>(database.header().elements));
		std::uint64_t boundary_faces = 0;
		for (std::uint64_t index = 0; index < database.header().elements; ++index) {
			auto element = database.Load(index);
			for (const auto label : element.boundary_labels)
				if (label >= 0) ++boundary_faces;
			elements.push_back(std::move(element));
		}
		std::cout << std::setprecision(17);
		if (argc == 3) {
			const auto velocity = iga::ReadVelocity(argv[2], database.header().nodes);
			const auto metrics = MeasureFlow(elements, velocity);
			PrintFlowMetrics(metrics, "");
			std::cout << "boundary_faces=" << boundary_faces << '\n';
			return 0;
		}
		if (std::string(argv[2]) == "--compare") {
			if (argc != 5)
				throw std::runtime_error("--compare requires REFERENCE.txt CURRENT.txt");
			const auto reference = iga::ReadVelocity(argv[3], database.header().nodes);
			const auto current = iga::ReadVelocity(argv[4], database.header().nodes);
			PrintFlowMetrics(MeasureFlow(elements, reference), "reference.");
			PrintFlowMetrics(MeasureFlow(elements, current), "current.");
			std::cout << "relative_velocity_l2="
				<< RelativeVelocityL2(reference, current)
				<< " boundary_faces=" << boundary_faces << '\n';
			return 0;
		}
		if (std::string(argv[2]) == "--compare-manifests") {
			if (argc != 7)
				throw std::runtime_error(
					"--compare-manifests requires two case-directory/manifest pairs");
			const fs::path reference_dir(argv[3]);
			const fs::path current_dir(argv[5]);
			const auto reference_snapshots = iga::ReadVelocityManifest(
				reference_dir/fs::path(argv[4]));
			const auto current_snapshots = iga::ReadVelocityManifest(
				current_dir/fs::path(argv[6]));
			double maximum_error = 0.0;
			for (std::size_t index = 0; index < reference_snapshots.size(); ++index) {
				const auto reference_time = reference_snapshots[index].time;
				const auto current_snapshot = std::find_if(
					current_snapshots.begin(), current_snapshots.end(),
					[&](const iga::VelocitySnapshot& snapshot) {
						return std::abs(reference_time-snapshot.time)
							<= 1e-12*std::max(
								{1.0, std::abs(reference_time), std::abs(snapshot.time)});
					});
				if (current_snapshot == current_snapshots.end())
					throw std::runtime_error(
						"current velocity manifest has no snapshot at reference time "
						+std::to_string(reference_time));
				const auto reference = iga::ReadVelocity(
					(reference_dir/reference_snapshots[index].file).string(),
					database.header().nodes);
				const auto current = iga::ReadVelocity(
					(current_dir/current_snapshot->file).string(),
					database.header().nodes);
				const auto error = RelativeVelocityL2(reference, current);
				maximum_error = std::max(maximum_error, error);
				std::cout << "snapshot[" << index << "].time=" << reference_time
					<< " snapshot[" << index << "].relative_velocity_l2="
					<< error << '\n';
			}
			std::cout << "reference_snapshots=" << reference_snapshots.size()
				<< " current_snapshots=" << current_snapshots.size()
				<< " maximum_relative_velocity_l2=" << maximum_error << '\n';
			return 0;
		}
		if (std::string(argv[2]) == "--womersley") {
			if (argc != 6)
				throw std::runtime_error(
					"--womersley requires configuration, case directory, and manifest");
			const auto configuration = iga::ReadWomersleyReferenceConfiguration(argv[3]);
			const fs::path current_dir(argv[4]);
			const auto current_snapshots = iga::ReadVelocityManifest(
				current_dir/fs::path(argv[5]));
			double maximum_error = 0.0;
			for (std::size_t index = 0; index < configuration.sample_times.size(); ++index) {
				const auto reference_time = configuration.sample_times[index];
				const auto current_snapshot = std::find_if(
					current_snapshots.begin(), current_snapshots.end(),
					[&](const iga::VelocitySnapshot& snapshot) {
						return std::abs(reference_time-snapshot.time)
							<= 1e-12*std::max(
								{1.0, std::abs(reference_time), std::abs(snapshot.time)});
					});
				if (current_snapshot == current_snapshots.end())
					throw std::runtime_error(
						"current velocity manifest has no Womersley sample time "
						+std::to_string(reference_time));
				const auto current = iga::ReadVelocity(
					(current_dir/current_snapshot->file).string(), database.header().nodes);
				const auto error = WomersleyVolumeRelativeL2(
					elements, current, configuration, reference_time);
				maximum_error = std::max(maximum_error, error);
				std::cout << "snapshot[" << index << "].time=" << reference_time
					<< " snapshot[" << index << "].womersley_volume_relative_l2="
					<< error << '\n';
			}
			std::cout << "womersley_snapshots=" << configuration.sample_times.size()
				<< " current_snapshots=" << current_snapshots.size()
				<< " maximum_womersley_volume_relative_l2=" << maximum_error << '\n';
			return 0;
		}
		if (std::string(argv[2]) != "--manifest")
			throw std::runtime_error("time-series mode requires --manifest");
		const fs::path case_dir(argv[3]);
		const auto snapshots = iga::ReadVelocityManifest(case_dir/fs::path(argv[4]));
		const double period = argc == 6 ? ParsePositiveDouble(argv[5], "PERIOD") : 0.0;
		double maximum_imbalance = 0.0;
		double maximum_cycle_error = 0.0;
		std::size_t cycle_pairs = 0;
		for (std::size_t index = 0; index < snapshots.size(); ++index) {
			const auto velocity = iga::ReadVelocity(
				(case_dir/snapshots[index].file).string(), database.header().nodes);
			const auto metrics = MeasureFlow(elements, velocity);
			const auto prefix = "snapshot["+std::to_string(index)+"].";
			std::cout << prefix << "time=" << snapshots[index].time << '\n';
			PrintFlowMetrics(metrics, prefix);
			maximum_imbalance = std::max(maximum_imbalance, metrics.relative_imbalance);
			if (!(period > 0.0)) continue;
			const auto target = snapshots[index].time-period;
			const auto previous = std::find_if(snapshots.begin(), snapshots.begin()+index,
				[&](const iga::VelocitySnapshot& item) {
					return std::abs(item.time-target)
						<= 1e-12*std::max({1.0, std::abs(item.time), std::abs(target)});
				});
			if (previous == snapshots.begin()+index) continue;
			const auto reference = iga::ReadVelocity(
				(case_dir/previous->file).string(), database.header().nodes);
			const auto error = RelativeVelocityL2(reference, velocity);
			std::cout << prefix << "cycle_relative_velocity_l2=" << error << '\n';
			maximum_cycle_error = std::max(maximum_cycle_error, error);
			++cycle_pairs;
		}
		std::cout << "snapshots=" << snapshots.size()
			<< " boundary_faces=" << boundary_faces
			<< " maximum_relative_mass_imbalance=" << maximum_imbalance;
		if (period > 0.0)
			std::cout << " cycle_period=" << period
				<< " cycle_pairs=" << cycle_pairs
				<< " maximum_cycle_relative_velocity_l2=" << maximum_cycle_error;
		std::cout << '\n';
		if (period > 0.0 && cycle_pairs == 0)
			throw std::runtime_error(
				"manifest contains no snapshot pairs separated by PERIOD");
		return 0;
	} catch (const std::exception& error) {
		std::cerr << "iga_flow_validate: " << error.what() << '\n';
		return 1;
	}
}
