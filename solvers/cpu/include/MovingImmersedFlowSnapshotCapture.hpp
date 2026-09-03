#ifndef IGA_MOVING_IMMERSED_FLOW_SNAPSHOT_CAPTURE_HPP
#define IGA_MOVING_IMMERSED_FLOW_SNAPSHOT_CAPTURE_HPP

// Read-only conversion of the outer committed epoch into the pure snapshot
// request.  This deliberately accepts no trial or prepared publication.
#include "MovingImmersedFlowSnapshot.hpp"
#include "MovingImmersedTransientFlowRuntime.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace iga {

inline MovingImmersedFlowSnapshot BuildCommittedMovingImmersedFlowSnapshot(
	const MovingImmersedTransientFlowRuntime& runtime,
	MovingImmersedFlowSnapshotOptions options = {})
{
	if (!runtime.Idle())
		throw std::logic_error("moving immersed snapshot capture requires an idle committed runtime");
	const auto& geometry = runtime.CommittedGeometry();
	const auto& state = runtime.CommittedGlobalState();
	MovingImmersedFlowSnapshotRequest request;
	request.time_s = state.TimeS(); request.index = state.Index();
	for (const int label : runtime.ConfiguredWallLabels()) {
		if (label <= 0) throw std::logic_error("moving immersed snapshot capture wall label is invalid");
		request.wall_labels.push_back(static_cast<std::uint32_t>(label));
	}
	for (const auto& port : runtime.ConfiguredPorts()) {
		if (port.boundary_label <= 0) throw std::logic_error("moving immersed snapshot capture port label is invalid");
		request.port_labels.push_back(static_cast<std::uint32_t>(port.boundary_label));
	}
	std::sort(request.wall_labels.begin(), request.wall_labels.end());
	std::sort(request.port_labels.begin(), request.port_labels.end());
	if (state.Index() == 0) {
		const auto& ports = runtime.CommittedDiagnostics().ports;
		if (ports.size() != request.port_labels.size())
			throw std::logic_error("moving immersed snapshot capture configured port identity is inconsistent");
		bool available = true;
		for (const auto& port : ports)
			available = available && port.measurement_valid && std::isfinite(port.measurement.outward_flow_m3_s);
		if (available) {
			request.port_flows_available = true;
			for (const auto& port : ports) request.port_flows.push_back({static_cast<std::uint64_t>(port.boundary_label), port.measurement.outward_flow_m3_s});
			std::sort(request.port_flows.begin(), request.port_flows.end(), [](const auto& a, const auto& b) { return a.label < b.label; });
		}
	} else {
		const auto& outer = runtime.Diagnostics();
		const auto conservation = runtime.ConservationDiagnostics();
		const double volume = geometry.Diagnostics().closed_surface_physical_volume_m3;
		if (outer.transition_identity_sha256.empty()
			|| conservation.transition_identity_sha256 != outer.transition_identity_sha256
			|| outer.target_geometry_identity_sha256 != geometry.GeometryIdentitySha256()
			|| outer.target_publication_identity_sha256 != geometry.PublicationIdentitySha256()
			|| outer.target_time_s != state.TimeS() || outer.target_index != state.Index()
			|| outer.target_audited_volume_m3 != volume
			|| conservation.target_geometry_identity_sha256 != geometry.GeometryIdentitySha256()
			|| conservation.target_publication_identity_sha256 != geometry.PublicationIdentitySha256()
			|| conservation.target_time_s != state.TimeS() || conservation.target_index != state.Index()
			|| !std::isfinite(volume) || conservation.target_audited_volume_m3 != volume)
			throw std::logic_error("moving immersed snapshot capture retained conservation record is inconsistent");
		request.transition_available = true; request.dt_s = conservation.dt_s;
		request.port_flows_available = true;
		for (const auto label : request.port_labels) {
			const auto found = conservation.fluid_surface_outward_flow_by_boundary_label_m3_s.find(static_cast<int>(label));
			if (found == conservation.fluid_surface_outward_flow_by_boundary_label_m3_s.end() || !std::isfinite(found->second))
				throw std::logic_error("moving immersed snapshot capture retained port flow is unavailable");
			request.port_flows.push_back({label, found->second});
		}
	}
	return MovingImmersedFlowSnapshot::Build(geometry, runtime.CommittedLayout(), state, std::move(request), options);
}

} // namespace iga

#endif
