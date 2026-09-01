#!/bin/sh
set -eu

checker=${1:-./iga_config_check}

v3_one_d_output=$("$checker" ../../examples/one_d/rigid_straight/simulation_config.json)
printf '%s\n' "$v3_one_d_output" | grep -q 'schema_version=3 dimension=1d'

v3_three_d='{"schema_version":3,"dimension":"3d","fields":[{"name":"velocity","kind":"vector3"},{"name":"pressure","kind":"pressure"}],"time":{"dt":0.1,"steps":1},"equation_systems":[{"name":"flow","kind":"navier_stokes","unknowns":["velocity","pressure"],"viscosity":0.004,"density":1060,"time_integration":"steady"}],"boundaries":[]}'
v3_three_d_output=$(printf '%s\n' "$v3_three_d" | "$checker" -)
printf '%s\n' "$v3_three_d_output" | grep -q 'schema_version=3 dimension=3d'

v4_output=$("$checker" ../../examples/vascular_flow/straight_tube/simulation_config.json)
printf '%s\n' "$v4_output" | grep -q 'schema_version=4 dimension=3d'

v5='{"schema_version":5,"time":{"dt":0.01,"steps":2},"start_domain":"upstream","execution":{"kind":"explicit"},"domains":[{"id":"upstream","dimension":"1d","kind":"network_flow","case":"domains/upstream","inlet_policy":"configured_open_loop","ports":[{"id":"terminal","locator_kind":"runtime_port","locator":"outlet:2","provides":["area","flow_rate","mean_pressure"],"requires":["mean_pressure"]}]},{"id":"roi","dimension":"3d","kind":"body_fitted_iga_flow","case":"domains/roi","database":"domains/roi/roi.ntiga","ports":[{"id":"inlet","locator_kind":"boundary_label","locator":"1","provides":["area","flow_rate","mean_pressure"],"requires":["flow_rate"]}]}],"couplings":[{"id":"upstream_to_roi","a":{"domain":"upstream","port":"terminal"},"b":{"domain":"roi","port":"inlet"},"mode":"pressure_flow","initial_pressure_pa":0}]}'
v5_output=$(printf '%s\n' "$v5" | "$checker" -)
printf '%s\n' "$v5_output" | grep -q 'schema_version=5 mode=multidomain'
printf '%s\n' "$v5_output" | grep -q 'sequential_pressure_flow_runner_compatible=yes'

printf '%s\n' 'iga_config_check dispatch tests passed'
