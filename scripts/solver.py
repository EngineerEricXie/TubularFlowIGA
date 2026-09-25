#!/usr/bin/env python3
"""Run existing 0D, 1D, IGA and native tetra FEM cases from one command."""

import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]


def path_from(base, value):
	return (base / value).resolve()


def read_json(path):
	return json.loads(path.read_text(encoding="utf-8"))


def run(command, dry_run):
	print(" ".join(str(part) for part in command), flush=True)
	if not dry_run:
		subprocess.run([str(part) for part in command], check=True)


def legacy_domain(config_path, data, args):
	dimension = data.get("dimension")
	if dimension in ("0d", "1d"):
		return {"id": "main", "dimension": dimension,
			"solver": "network_" + dimension, "case_dir": str(config_path.parent)}
	if dimension == "3d":
		return {"id": "main", "dimension": "3d", "method": "IGA",
			"device": args.device, "case_dir": str(config_path.parent),
			"database": args.database, "physics": args.physics}
	if data.get("backend", "").startswith("native_tet_"):
		return {"id": "main", "dimension": "3d", "method": "FEM",
			"device": "CPU", "workflow_file": str(config_path)}
	if data.get("kind") == "idealized_cube_dual_tree_functional_only":
		return {"id": "main", "dimension": "3d", "method": "FEM",
			"device": "CPU", "backend": "idealized_cube_dual_tree",
			"case_file": str(config_path), "physics": ["flow", "darcy", "solid"]}
	if all(key in data.get("meshes", {}) for key in ("artery", "tissue", "vein")):
		return {"id": "main", "dimension": "3d", "method": "FEM",
			"device": "CPU", "backend": "native_tet_flow_darcy_chain",
			"case_file": str(config_path), "physics": ["flow", "darcy"]}
	if "fluid_mesh_file" in data and "solid_mesh_file" in data:
		return {"id": "main", "dimension": "3d", "method": "FEM",
			"device": "GPU" if args.backend == "native_tet_fsi_cuda" else "CPU",
			"backend": "native_tet_fsi_cuda" if args.backend == "native_tet_fsi_cuda"
				else "native_tet_fsi_steps",
			"case_file": str(config_path), "physics": ["flow", "solid"]}
	if "mesh_file" in data:
		return {"id": "main", "dimension": "3d", "method": "FEM",
			"device": "GPU" if args.backend in (
				"native_tet_flow_cuda", "native_tet_darcy_cuda",
				"native_tet_species_cuda") else "CPU",
			"case_file": str(config_path),
			"backend": args.backend}
	raise ValueError("unknown legacy config; specify a solver-v1 config")


def validate_domains(domains, connections):
	if not domains:
		raise ValueError("domains must contain at least one domain")
	ids = [domain["id"] for domain in domains]
	if len(set(ids)) != len(ids):
		raise ValueError("domain ids must be unique")
	by_id = {domain["id"]: domain for domain in domains}
	for domain in domains:
		if domain["dimension"] not in ("0d", "1d", "3d"):
			raise ValueError("dimension must be 0d, 1d or 3d")
		if domain["dimension"] == "3d":
			if domain.get("method") not in ("IGA", "FEM"):
				raise ValueError("3d method must be IGA or FEM")
			if domain.get("device", "CPU") not in ("CPU", "GPU"):
				raise ValueError("3d device must be CPU or GPU")
		else:
			if domain.get("device", "CPU") != "CPU" or "method" in domain:
				raise ValueError("0d/1d domains use the CPU network solver")
	for edge in connections:
		for end in ("from", "to"):
			domain_id, separator, port = edge[end].partition(":")
			if not separator or domain_id not in by_id or port not in by_id[domain_id].get("ports", {}):
				raise ValueError(f"unknown connection endpoint {edge[end]}")


def domain_case(base, domain):
	value = domain.get("case_file")
	if value is None:
		raise ValueError(f"{domain['id']} requires case_file")
	return path_from(base, value)


def runtime_case(case_dir, output, time, waveforms, every_steps, dry_run):
	if not time and not waveforms and not every_steps:
		return case_dir
	effective = output / "runtime_case"
	if dry_run:
		return effective
	effective.mkdir(parents=True, exist_ok=True)
	for asset in case_dir.iterdir():
		if asset.name != "simulation_config.json":
			link = effective / asset.name
			if not link.exists():
				link.symlink_to(asset)
	case = read_json(case_dir / "simulation_config.json")
	if time:
		if "dt_s" in time:
			case["time"]["dt"] = time["dt_s"]
		if "steps" in time:
			case["time"]["steps"] = time["steps"]
	if waveforms:
		case["temporal_functions"] = waveforms
	if every_steps:
		case["time"]["output_every"] = every_steps
	(effective / "simulation_config.json").write_text(
		json.dumps(case, indent=2)+"\n", encoding="utf-8")
	return effective


def native_hydraulic_graph(domains, connections, time):
	by_id = {domain["id"]: domain for domain in domains}
	vessels = [domain for domain in domains
		if domain["dimension"] == "3d" and domain["method"] == "FEM"]
	if len(vessels) != 1:
		raise ValueError("native hydraulic graph needs one FEM vessel")
	vessel = vessels[0]
	if vessel.get("device", "CPU") != "CPU":
		raise ValueError("this hydraulic graph route requires a CPU FEM vessel")
	if vessel.get("physics") not in (["flow"], ["flow", "species"]):
		raise ValueError("native hydraulic graph supports FEM flow with optional species")
	sources = [domain for domain in domains
		if domain.get("solver") in ("source", "network_1d")]
	terminals = [domain for domain in domains if domain.get("solver") == "rcr"]
	if len(sources) != 1 or len(terminals) != len(domains)-2:
		raise ValueError("native hydraulic graph needs one source and RCR terminals")
	source = sources[0]
	inlet = next((port for port, info in vessel["ports"].items()
		if info.get("role") == "inlet"), None)
	if inlet is None or len([info for info in vessel["ports"].values()
			if info.get("role") == "inlet"]) != 1:
		raise ValueError("native vessel needs one labelled inlet port")
	if len(connections) != len(domains)-1:
		raise ValueError("every source and terminal needs one connection")
	for edge in connections:
		allowed = [["pressure", "flow_rate"]]
		if source["dimension"] == "1d" and "species" in vessel.get("physics", []) \
				and edge["from"].startswith(source["id"]+":"):
			allowed.append(["pressure", "flow_rate", "species_concentration",
				"species_flux"])
		if edge.get("exchange") not in allowed \
				or edge.get("scheme", "sequential") != "sequential":
			raise ValueError("native hydraulic graph connection exchange differs")
	source_edges = [edge for edge in connections
		if edge["from"].startswith(source["id"]+":")
			and edge["to"] == vessel["id"]+":"+inlet]
	if len(source_edges) != 1:
		raise ValueError("source must connect to the FEM inlet")
	rcrs = []
	for terminal in terminals:
		matching = [edge for edge in connections
			if edge["to"].startswith(terminal["id"]+":")]
		if len(matching) != 1:
			raise ValueError("each RCR needs one vessel outlet")
		vessel_id, _, port_name = matching[0]["from"].partition(":")
		if vessel_id != vessel["id"] or vessel["ports"][port_name].get("role") != "outlet":
			raise ValueError("RCR must connect to a FEM outlet")
		rcrs.append({"boundary_label": vessel["ports"][port_name]["label"],
			**terminal["materials"]})
	if len({item["boundary_label"] for item in rcrs}) != len(rcrs):
		raise ValueError("FEM outlets must have distinct labels")
	case = {"schema_version": 2,
		"mesh_file": vessel["geometry"],
		"boundaries": {"inlet_label": vessel["ports"][inlet]["label"],
			"wall_labels": vessel["boundaries"]["wall_labels"]},
		"fluid": vessel["materials"],
		"terminal_rcrs": rcrs,
		"motion": vessel.get("motion", {"kind": "fixed"}),
		"coupling": vessel["coupling"],
		"time": {"dt_s": time["dt_s"], "steps": time["steps"]}}
	if source["dimension"] == "1d":
		port_name = source_edges[0]["from"].split(":", 1)[1]
		case["source_1d"] = {"case_dir": source["case_dir"],
			"outlet_node_id": source["ports"][port_name]["node_id"]}
		if source.get("species_bindings"):
			case["source_1d"]["species_bindings"] = source["species_bindings"]
		if source.get("system"):
			case["source_1d"]["system"] = source["system"]
	else:
		case["source"] = source["materials"]
	if "species" in vessel.get("physics", []):
		case["species"] = vessel["species"]
	return case


def native_flow_darcy_graph(base, domains, connections, time):
	by_id = {domain["id"]: domain for domain in domains}
	if len(domains) != 3 or any(domain["dimension"] != "3d"
			or domain["method"] != "FEM" or domain.get("device", "CPU") != "CPU"
			for domain in domains):
		raise ValueError("flow-Darcy graph needs three CPU FEM domains")
	artery, tissue, vein = (next((domain for domain in domains
		if domain.get("role") == role), None)
		for role in ("artery", "tissue", "vein"))
	if not all((artery, tissue, vein)):
		raise ValueError("flow-Darcy graph needs artery, tissue and vein roles")
	with_species = "species" in artery["physics"]
	if artery["physics"] != (["flow", "species"] if with_species else ["flow"]) \
			or tissue["physics"] != (["darcy", "species"] if with_species else ["darcy"]) \
			or vein["physics"] != (["flow", "species"] if with_species else ["flow"]):
		raise ValueError("flow-Darcy graph physics must match across its three domains")
	arterial, venous = [], []
	for edge in connections:
		if edge.get("exchange") != ["flow_rate"] or edge.get("scheme") != "sequential":
			raise ValueError("flow-Darcy graph uses sequential flow-rate exchange")
		first_domain, first_port = edge["from"].split(":", 1)
		second_domain, second_port = edge["to"].split(":", 1)
		first = by_id[first_domain]["ports"][first_port]
		second = by_id[second_domain]["ports"][second_port]
		if first_domain == artery["id"] and second_domain == tissue["id"]:
			arterial.append({"name": first_port, "artery_label": first["label"],
				"tissue_label": second["label"]})
		elif first_domain == tissue["id"] and second_domain == vein["id"]:
			venous.append({"tissue_label": first["label"],
				"vein_label": second["label"]})
		else:
			raise ValueError("flow-Darcy edge must connect artery to tissue or tissue to vein")
	if not arterial or not venous:
		raise ValueError("flow-Darcy graph needs both exchange directions")
	artery_inlet = next((info for info in artery["ports"].values()
		if info.get("role") == "inlet"), artery["ports"].get("inlet"))
	vein_outlet = next((info for info in vein["ports"].values()
		if info.get("role") == "outlet"), vein["ports"].get("outlet"))
	if artery_inlet is None or vein_outlet is None:
		raise ValueError("flow-Darcy graph needs artery inlet and vein outlet ports")
	case = {"schema_version": 1,
		"meshes": {role: str(path_from(base, domain["geometry"]))
			for role, domain in (("artery", artery), ("tissue", tissue),
				("vein", vein))},
		"artery": {"inlet_label": artery_inlet["label"],
			"wall_label": artery["boundaries"]["wall_label"],
			"inlet_velocity_m_s": artery["boundaries"]["inlet_velocity_m_s"]},
		"vein": {"outlet_label": vein_outlet["label"],
			"wall_label": vein["boundaries"]["wall_label"],
			"outlet_pressure_pa": vein["boundaries"]["outlet_pressure_pa"]},
		"materials": {**artery["materials"], **tissue["materials"]},
		"artery_to_tissue": arterial, "tissue_to_vein": venous}
	if with_species:
		case["species"] = artery["species"]
		case["time"] = {"dt_s": time["dt_s"], "steps": time["steps"]}
	return case


def native_gpu_hydraulic_graph(base, domains, connections, time):
	vessels = [domain for domain in domains if domain["dimension"] == "3d"]
	sources = [domain for domain in domains if domain.get("solver") == "source"]
	terminals = [domain for domain in domains if domain.get("solver") == "rcr"]
	if len(vessels) != 1 or len(sources) != 1 \
			or len(terminals) != len(domains)-2:
		raise ValueError("GPU hydraulic graph needs one source, one vessel, and RCR terminals")
	vessel, source = vessels[0], sources[0]
	if vessel.get("method") != "FEM" or vessel.get("device") != "GPU" \
			or vessel.get("backend") != "native_tet_flow_cuda" \
			or vessel.get("physics") != ["flow"]:
		raise ValueError("GPU hydraulic graph needs one native FEM CUDA flow vessel")
	inlet = next((name for name, port in vessel["ports"].items()
		if port.get("role") == "inlet"), None)
	if inlet is None:
		raise ValueError("GPU hydraulic vessel needs a labelled inlet")
	source_edges = [edge for edge in connections
		if edge["from"].startswith(source["id"]+":")
			and edge["to"] == vessel["id"]+":"+inlet]
	if len(source_edges) != 1 or len(connections) != len(domains)-1:
		raise ValueError("GPU source must connect once to the vessel inlet")
	outlets = []
	for terminal in terminals:
		matching = [edge for edge in connections
			if edge["to"].startswith(terminal["id"]+":")]
		if len(matching) != 1:
			raise ValueError("each GPU RCR needs one vessel outlet")
		domain_id, port_name = matching[0]["from"].split(":", 1)
		if domain_id != vessel["id"] or vessel["ports"][port_name].get("role") != "outlet":
			raise ValueError("GPU RCR must connect from a vessel outlet")
		outlets.append({"boundary_label": vessel["ports"][port_name]["label"],
			**terminal["materials"]})
	for edge in connections:
		if edge.get("exchange") != ["pressure", "flow_rate"] \
				or edge.get("scheme", "sequential") != "sequential":
			raise ValueError("GPU T7 ports use sequential pressure-flow exchange")
	case_path = domain_case(base, vessel)
	case = read_json(case_path)
	if any(case["time"].get(key) != time.get(key) for key in ("dt_s", "steps")):
		raise ValueError("GPU hydraulic graph outer time differs from its flow case")
	case["mesh_file"] = str(path_from(case_path.parent, case["mesh_file"]))
	case["t7_ports"] = {
		"inlet_boundary_label": vessel["ports"][inlet]["label"],
		"source": source["materials"],
		"outlets": outlets}
	return case


def run_1d_terminal_graph(base, domains, connections, time, waveforms,
		output, every_steps, dry_run, ranks, stop_after_step):
	network = next(domain for domain in domains if domain["dimension"] == "1d")
	terminals = {domain["id"]: domain for domain in domains
		if domain["dimension"] == "0d"}
	case_dir = path_from(base, network["case_dir"])
	case = read_json(case_dir / "simulation_config.json")
	boundaries = case["boundaries"]
	connected_terminals = set()
	connected_nodes = set()
	for edge in connections:
		first_id, first_port = edge["from"].split(":", 1)
		last_id, last_port = edge["to"].split(":", 1)
		if first_id != network["id"] or last_id not in terminals \
				or edge.get("exchange") != ["pressure", "flow_rate"] \
				or edge.get("scheme", "sequential") != "sequential":
			raise ValueError("1D terminal edges need network outlet to R/RC/RCR inlet")
		port = network["ports"][first_port]
		terminal = terminals[last_id]
		if port["role"] != "outlet" or terminal["ports"][last_port]["role"] != "inlet":
			raise ValueError("1D terminal port roles must be outlet to inlet")
		node_id = port["node_id"]
		if last_id in connected_terminals or node_id in connected_nodes:
			raise ValueError("each 1D terminal needs a distinct outlet node")
		connected_terminals.add(last_id)
		connected_nodes.add(node_id)
		for boundary in boundaries:
			if node_id in boundary.get("node_ids", []):
				boundary["node_ids"].remove(node_id)
		kind = terminal["solver"]
		material = terminal["materials"]
		condition = {"field": "pressure", "type": {
			"r": "resistance", "rc": "windkessel_rc",
			"rcr": "windkessel_rcr"}[kind]}
		if kind in ("r", "rc"):
			condition["resistance"] = material["resistance_pa_s_m3"]
		else:
			condition["proximal_resistance"] = material[
				"proximal_resistance_pa_s_m3"]
			condition["distal_resistance"] = material[
				"distal_resistance_pa_s_m3"]
		if kind in ("rc", "rcr"):
			condition["capacitance"] = material["capacitance_m3_pa"]
		condition["reference_pressure"] = material.get("reference_pressure_pa", 0.0)
		condition["initial_pressure"] = material.get("initial_pressure_pa", 0.0)
		boundaries.append({"name": last_id, "role": "outlet",
			"node_ids": [node_id], "conditions": [condition]})
	case["boundaries"] = [boundary for boundary in boundaries
		if boundary.get("node_ids", [1])]
	if connected_terminals != set(terminals):
		raise ValueError("every 0D terminal needs a 1D outlet connection")
	if time:
		if "dt_s" in time:
			case["time"]["dt"] = time["dt_s"]
		if "steps" in time:
			case["time"]["steps"] = time["steps"]
	if waveforms:
		case["temporal_functions"] = waveforms
	if every_steps:
		case["time"]["output_every"] = every_steps
	generated = output / "runtime_case"
	if not dry_run:
		generated.mkdir(parents=True, exist_ok=True)
		for asset in case_dir.iterdir():
			if asset.name != "simulation_config.json":
				link = generated / asset.name
				if not link.exists():
					link.symlink_to(asset)
		(generated / "simulation_config.json").write_text(
			json.dumps(case, indent=2)+"\n", encoding="utf-8")
	command = ["mpiexec", "-np", str(ranks), ROOT / "solvers/one_d/iga_1d",
		generated, "--output-dir", output / "fields"]
	if network.get("system"):
		command += ["--system", network["system"]]
	if stop_after_step:
		command += ["--stop-after-step", str(stop_after_step)]
	run(command + network.get("options", []), dry_run)


def run_graph(base, domains, connections, time, waveforms, every_steps,
		output, dry_run, ranks, args):
	if any(domain["dimension"] == "1d" for domain in domains) \
			and not any(domain["dimension"] == "3d" for domain in domains):
		if args.resume:
			raise ValueError("1D terminal graph restart uses the 1D backend options")
		if sum(domain["dimension"] == "1d" for domain in domains) != 1 \
				or any(domain["dimension"] not in ("0d", "1d")
					or (domain["dimension"] == "0d" and domain.get("solver")
						not in ("r", "rc", "rcr")) for domain in domains) \
				or len(connections) != len(domains)-1:
			raise ValueError("1D terminal graph needs one network and R/RC/RCR terminals")
		run_1d_terminal_graph(base, domains, connections, time, waveforms,
			output, every_steps, dry_run, ranks, args.stop_after_step)
		return
	gpu_vessels = [domain for domain in domains
		if domain["dimension"] == "3d" and domain.get("device") == "GPU"]
	if gpu_vessels and any(domain["dimension"] == "0d" for domain in domains):
		if ranks != 1 or waveforms or every_steps not in (None, 1) \
				or args.resume or args.stop_after_step:
			raise ValueError("GPU T7 hydraulic graph needs one rank and every-step output")
		case = native_gpu_hydraulic_graph(base, domains, connections, time)
		case_path = output / "case.json"
		if not dry_run:
			output.mkdir(parents=True, exist_ok=True)
			case_path.write_text(json.dumps(case, indent=2)+"\n", encoding="utf-8")
		run([ROOT / "solvers/cuda/native_tet_flow_cuda", case_path,
			output / "fields"], dry_run)
		return
	if all(domain["dimension"] == "3d" and domain.get("method") == "FEM"
			for domain in domains):
		gpu_drivers = [domain for domain in domains if domain.get("backend") in (
			"native_tet_flow_cuda", "native_tet_darcy_cuda")]
		gpu_flows = [domain for domain in domains
			if domain.get("backend") == "native_tet_flow_cuda"]
		gpu_darcys = [domain for domain in domains
			if domain.get("backend") == "native_tet_darcy_cuda"]
		gpu_species = [domain for domain in domains
			if domain.get("backend") == "native_tet_species_cuda"]
		if len(gpu_flows) == 1 and len(gpu_darcys) == 1 \
				and len(gpu_species) == len(domains)-2:
			if ranks != 1 or waveforms or every_steps not in (None, 1) \
					or args.resume or args.stop_after_step \
					or any(domain.get("device") != "GPU" for domain in domains):
				raise ValueError("GPU flow-Darcy graph needs one GPU and every-step output")
			flow, darcy = gpu_flows[0], gpu_darcys[0]
			if flow.get("physics") != ["flow"] or darcy.get("physics") != ["darcy"]:
				raise ValueError("GPU flow-Darcy graph physics differs")
			flow_edges = [edge for edge in connections
				if edge["from"].split(":", 1)[0] == flow["id"]
					and edge["to"].split(":", 1)[0] == darcy["id"]]
			species_edges = [edge for edge in connections
				if edge["from"].split(":", 1)[0] == darcy["id"]
					and any(edge["to"].split(":", 1)[0] == species["id"]
						for species in gpu_species)]
			if not flow_edges or len(connections) != len(flow_edges)+len(species_edges) \
					or len(species_edges) != len(gpu_species):
				raise ValueError("GPU flow-Darcy graph connections differ")
			ports = []
			for edge in flow_edges:
				flow_port = edge["from"].split(":", 1)[1]
				darcy_port = edge["to"].split(":", 1)[1]
				if edge.get("exchange") != ["flow_rate", "source"] \
						or edge.get("scheme", "sequential") != "sequential" \
						or "label" not in flow["ports"][flow_port] \
						or "cell_id_weights" not in darcy["ports"][darcy_port]:
					raise ValueError("GPU flow-Darcy source port differs")
				ports.append({"name": flow_port,
					"vessel_boundary_label": flow["ports"][flow_port]["label"],
					"tissue_cell_id_weights": darcy["ports"][darcy_port][
						"cell_id_weights"]})
			flow_case, darcy_case = domain_case(base, flow), domain_case(base, darcy)
			flow_data, darcy_data = read_json(flow_case), read_json(darcy_case)
			if flow_data["time"] != {key: time[key] for key in ("dt_s", "steps")}:
				raise ValueError("GPU flow-Darcy outer time differs from flow case")
			for species in gpu_species:
				matching = [edge for edge in species_edges
					if edge["to"].split(":", 1)[0] == species["id"]]
				species_case = domain_case(base, species)
				species_data = read_json(species_case)
				if species.get("physics") != ["species"] or len(matching) != 1 \
						or matching[0].get("exchange") != ["flow_rate"] \
						or matching[0].get("scheme", "sequential") != "sequential" \
						or path_from(species_case.parent, species_data["mesh_file"]) \
							!= path_from(darcy_case.parent, darcy_data["mesh_file"]) \
						or any(species_data["time"].get(key) != time.get(key)
							for key in ("dt_s", "steps")):
					raise ValueError("GPU Darcy-species connection, mesh, or time differs")
			flow_output, darcy_output = output / "flow", output / "darcy"
			source_map_path = output / "flow_darcy_source_map.json"
			source_map = {"schema_version": 1,
				"vessel_mesh_file": str(path_from(flow_case.parent, flow_data["mesh_file"])),
				"flow_velocity_file": str(flow_output / (
					"flow_velocity_step_"+str(time["steps"])+".bin")),
				"ports": ports}
			if not dry_run:
				output.mkdir(parents=True, exist_ok=True)
				source_map_path.write_text(json.dumps(source_map, indent=2)+"\n",
					encoding="utf-8")
			run([ROOT / "solvers/cuda/native_tet_flow_cuda", flow_case, flow_output],
				dry_run)
			run([ROOT / "solvers/cuda/native_tet_darcy_cuda", darcy_case,
				darcy_output, "--source-map", source_map_path], dry_run)
			for species in gpu_species:
				run([ROOT / "solvers/cuda/native_tet_species_cuda",
					domain_case(base, species), output / ("species_"+species["id"]),
					"--darcy-flux", darcy_output], dry_run)
			return
		if len(domains) > 2 and len(gpu_drivers) == 1 \
				and len(gpu_species) == len(domains)-1:
			if ranks != 1 or waveforms or every_steps not in (None, 1) \
					or args.resume or args.stop_after_step \
					or any(domain.get("device") != "GPU" for domain in domains):
				raise ValueError("GPU multi-species graph needs one GPU and every-step output")
			driver = gpu_drivers[0]
			driver_case = domain_case(base, driver)
			driver_data = read_json(driver_case)
			flow_driver = driver["backend"] == "native_tet_flow_cuda"
			exchange = ["velocity"] if flow_driver else ["flow_rate"]
			if driver.get("physics") != (["flow"] if flow_driver else ["darcy"]) \
					or len(connections) != len(gpu_species):
				raise ValueError("GPU multi-species graph needs one driver and one edge per species")
			for species in gpu_species:
				matching = [edge for edge in connections
					if edge["to"].startswith(species["id"]+":")]
				case = domain_case(base, species)
				case_data = read_json(case)
				if species.get("physics") != ["species"] or len(matching) != 1 \
						or not matching[0]["from"].startswith(driver["id"]+":") \
						or matching[0].get("exchange") != exchange \
						or matching[0].get("scheme", "sequential") != "sequential" \
						or path_from(case.parent, case_data["mesh_file"]) != path_from(
							driver_case.parent, driver_data["mesh_file"]) \
						or (flow_driver and case_data["time"] != driver_data["time"]) \
						or any(case_data["time"].get(key) != time.get(key)
							for key in ("dt_s", "steps")):
					raise ValueError("GPU species connection, mesh, or time differs")
			driver_output = output / ("flow" if flow_driver else "darcy")
			run([ROOT / "solvers/cuda" / driver["backend"],
				driver_case, driver_output], dry_run)
			for species in gpu_species:
				run([ROOT / "solvers/cuda/native_tet_species_cuda",
					domain_case(base, species), output / ("species_"+species["id"]),
					"--velocity-series" if flow_driver else "--darcy-flux",
					driver_output], dry_run)
			return
		if len(domains) == 2 and {domain.get("backend") for domain in domains} \
				== {"native_tet_darcy_cuda", "native_tet_species_cuda"}:
			if ranks != 1 or waveforms or every_steps not in (None, 1) \
					or args.resume or args.stop_after_step:
				raise ValueError("GPU Darcy-species graph needs one rank and every-step output")
			darcy = next(domain for domain in domains
				if domain["backend"] == "native_tet_darcy_cuda")
			species = next(domain for domain in domains
				if domain["backend"] == "native_tet_species_cuda")
			if any(domain.get("device") != "GPU" for domain in domains) \
					or darcy.get("physics") != ["darcy"] \
					or species.get("physics") != ["species"] \
					or len(connections) != 1 \
					or connections[0].get("exchange") != ["flow_rate"] \
					or connections[0].get("scheme", "sequential") != "sequential" \
					or not connections[0]["from"].startswith(darcy["id"]+":") \
					or not connections[0]["to"].startswith(species["id"]+":"):
				raise ValueError("GPU Darcy-species graph needs a sequential flow-rate connection")
			darcy_case = domain_case(base, darcy)
			species_case = domain_case(base, species)
			darcy_data, species_data = read_json(darcy_case), read_json(species_case)
			if path_from(darcy_case.parent, darcy_data["mesh_file"]) != path_from(
					species_case.parent, species_data["mesh_file"]) \
					or any(species_data["time"].get(key) != time.get(key)
						for key in ("dt_s", "steps")):
				raise ValueError("GPU Darcy and species must use the same mesh and species time")
			run([ROOT / "solvers/cuda/native_tet_darcy_cuda",
				darcy_case, output / "darcy"], dry_run)
			run([ROOT / "solvers/cuda/native_tet_species_cuda",
				species_case, output / "species", "--darcy-flux",
				output / "darcy"], dry_run)
			return
		if len(domains) == 2 and {domain.get("backend") for domain in domains} \
				== {"native_tet_flow_cuda", "native_tet_species_cuda"}:
			if ranks != 1 or waveforms or every_steps not in (None, 1) \
					or args.resume or args.stop_after_step:
				raise ValueError("GPU flow-species graph needs one rank and every-step output")
			flow = next(domain for domain in domains
				if domain["backend"] == "native_tet_flow_cuda")
			species = next(domain for domain in domains
				if domain["backend"] == "native_tet_species_cuda")
			if any(domain.get("device") != "GPU" for domain in domains) \
					or flow.get("physics") != ["flow"] \
					or species.get("physics") != ["species"] \
					or len(connections) != 1 \
					or connections[0].get("exchange") != ["velocity"] \
					or connections[0].get("scheme", "sequential") != "sequential" \
					or not connections[0]["from"].startswith(flow["id"]+":") \
					or not connections[0]["to"].startswith(species["id"]+":"):
				raise ValueError("GPU flow-species graph needs a sequential velocity connection")
			flow_case = domain_case(base, flow)
			species_case = domain_case(base, species)
			flow_data, species_data = read_json(flow_case), read_json(species_case)
			if path_from(flow_case.parent, flow_data["mesh_file"]) != path_from(
					species_case.parent, species_data["mesh_file"]) \
					or flow_data["time"] != species_data["time"] \
					or any(flow_data["time"].get(key) != time.get(key)
						for key in ("dt_s", "steps")):
				raise ValueError("GPU flow and species must use the same mesh and time")
			run([ROOT / "solvers/cuda/native_tet_flow_cuda",
				flow_case, output / "flow"], dry_run)
			run([ROOT / "solvers/cuda/native_tet_species_cuda",
				species_case, output / "species", "--velocity-series",
				output / "flow"], dry_run)
			return
		if args.resume or args.stop_after_step:
			raise ValueError("flow-Darcy graph does not have a checkpoint")
		case = native_flow_darcy_graph(base, domains, connections, time)
		case_path = output / "case.json"
		if not dry_run:
			output.mkdir(parents=True, exist_ok=True)
			case_path.write_text(json.dumps(case, indent=2)+"\n", encoding="utf-8")
		run(["mpiexec", "-np", str(ranks),
			ROOT / "solvers/cpu/native_tet_flow_darcy_chain",
			case_path, output / "fields"], dry_run)
		return
	case = native_hydraulic_graph(domains, connections, time)
	if case.get("source_1d"):
		if args.resume or args.stop_after_step:
			raise ValueError("1D/FEM graph checkpoint and partial run are unavailable")
		case["source_1d"]["case_dir"] = str(path_from(base,
			case["source_1d"]["case_dir"]))
	if case.get("species") and args.resume:
		raise ValueError("coupled flow/species restart is not available")
	case["mesh_file"] = str(path_from(base, case["mesh_file"]))
	case_path = output / "case.json"
	fields = output / "fields"
	checkpoints = output / "checkpoints"
	if not dry_run:
		output.mkdir(parents=True, exist_ok=True)
		case_path.write_text(json.dumps(case, indent=2)+"\n", encoding="utf-8")
	command = ["mpiexec", "-np", str(ranks),
		ROOT / "solvers/cpu/native_tet_hydraulic_graph", case_path]
	if not args.resume:
		run(command+["--check-input"], dry_run)
	if not case.get("species") and not case.get("source_1d"):
		command += ["--checkpoint-dir", checkpoints]
	command += ["--output-dir", fields]
	if args.resume:
		command += ["--restart-dir", checkpoints]
	if args.stop_after_step:
		command += ["--stop-after-step", str(args.stop_after_step)]
	run(command, dry_run)
	if not dry_run:
		write_pvd(fields, output / "series.pvd", time["dt_s"])


def write_pvd(fields, destination, dt):
	steps = sorted((int(item.name[5:]), item / "snapshot.pvtu")
		for item in fields.glob("step_*") if item.name[5:].isdigit())
	lines = ['<?xml version="1.0"?>',
		'<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">',
		'  <Collection>']
	for step, snapshot in steps:
		lines.append(f'    <DataSet timestep="{step*dt:.17g}" group="" part="0" file="{snapshot.relative_to(destination.parent)}"/>')
	lines += ['  </Collection>', '</VTKFile>']
	destination.write_text("\n".join(lines)+"\n", encoding="utf-8")


def run_one(base, domain, output, dry_run, ranks, args, time, waveforms, every_steps):
	dimension = domain["dimension"]
	device = domain.get("device", "CPU")
	if dimension in ("0d", "1d"):
		case_dir = path_from(base, domain["case_dir"])
		case_dir = runtime_case(case_dir, output, time, waveforms, every_steps, dry_run)
		binary = ROOT / "solvers/one_d" / ("iga_0d" if dimension == "0d" else "iga_1d")
		command = ["mpiexec", "-np", str(ranks), binary, case_dir,
			"--output-dir", output]
		if domain.get("system"):
			command += ["--system", domain["system"]]
		command += domain.get("options", [])
		run(command, dry_run)
		return
	if domain["method"] == "FEM":
		if domain.get("backend") == "native_tet_fsi_cuda":
			if device != "GPU" or ranks != 1 \
					or domain.get("physics", ["flow", "solid"]) != ["flow", "solid"]:
				raise ValueError("native FEM CUDA FSI needs one GPU flow-solid domain")
			if time.get("mode") not in (None, "transient") or waveforms \
					or every_steps not in (None, 1) or args.resume or args.stop_after_step:
				raise ValueError("native FEM CUDA FSI advances every transient step")
			case = domain_case(base, domain)
			case_time = read_json(case)["time"]
			if any(case_time.get(key) != time[key]
					for key in ("dt_s", "steps") if key in time):
				raise ValueError("FEM CUDA FSI outer time must match the referenced case")
			run([ROOT / "solvers/cuda/native_tet_fsi_cuda", case, output], dry_run)
			return
		if domain.get("backend") == "native_tet_species_cuda":
			if device != "GPU" or ranks != 1 or domain.get("physics", ["species"]) != ["species"]:
				raise ValueError("native FEM CUDA species needs one GPU species domain")
			if time.get("mode") not in (None, "transient") or waveforms \
					or every_steps not in (None, 1) or args.resume or args.stop_after_step:
				raise ValueError("native FEM CUDA species advances transient steps")
			case = domain_case(base, domain)
			case_time = read_json(case)["time"]
			if any(case_time.get(key) != time[key]
					for key in ("dt_s", "steps") if key in time):
				raise ValueError("FEM CUDA outer time must match the referenced species case")
			run([ROOT / "solvers/cuda/native_tet_species_cuda", case, output], dry_run)
			return
		if domain.get("backend") == "native_tet_darcy_cuda":
			if device != "GPU" or ranks != 1 or domain.get("physics", ["darcy"]) != ["darcy"]:
				raise ValueError("native FEM CUDA Darcy needs one GPU Darcy domain")
			if time.get("mode") not in (None, "steady") or waveforms or every_steps \
					or args.resume or args.stop_after_step:
				raise ValueError("native FEM CUDA Darcy is a steady solve")
			run([ROOT / "solvers/cuda/native_tet_darcy_cuda",
				domain_case(base, domain), output], dry_run)
			return
		if domain.get("backend") == "native_tet_flow_cuda":
			if device != "GPU" or ranks != 1 or domain.get("physics", ["flow"]) != ["flow"]:
				raise ValueError("native FEM CUDA flow needs one GPU flow domain")
			if time.get("mode") not in (None, "transient"):
				raise ValueError("native FEM CUDA flow advances transient steps")
			if waveforms or every_steps not in (None, 1) or args.resume \
					or args.stop_after_step:
				raise ValueError("native FEM CUDA flow writes every step without restart")
			case = domain_case(base, domain)
			case_time = read_json(case)["time"]
			if any(case_time.get(key) != time[key]
					for key in ("dt_s", "steps") if key in time):
				raise ValueError("FEM CUDA outer time must match the referenced case")
			run([ROOT / "solvers/cuda/native_tet_flow_cuda", case, output], dry_run)
			return
		if device != "CPU":
			raise ValueError("the selected FEM backend has no GPU route")
		if domain.get("backend") == "native_tet_fsi_steps":
			if ranks != 1 or waveforms or every_steps not in (None, 1):
				raise ValueError("native FEM FSI uses one CPU rank and every-step output")
			case = domain_case(base, domain)
			case_time = read_json(case)["time"]
			if any(case_time.get(key) != time[key]
					for key in ("dt_s", "steps") if key in time):
				raise ValueError("FEM FSI outer time must match the referenced case")
			run(["mpiexec", "-np", "1", ROOT / "solvers/cpu/native_tet_fsi_steps",
				case, output], dry_run)
			return
		if domain.get("backend") == "native_tet_flow_darcy_chain":
			case = domain_case(base, domain)
			case_data = read_json(case)
			with_species = bool(case_data.get("species"))
			if waveforms or every_steps or time.get("mode") not in (
				None, "transient" if with_species else "steady"):
				raise ValueError("flow-Darcy case time or waveform differs from its backend")
			if any(case_data.get("time", {}).get(key) != time[key]
					for key in ("dt_s", "steps") if key in time):
				raise ValueError("flow-Darcy outer time must match the referenced case")
			run(["mpiexec", "-np", str(ranks),
				ROOT / "solvers/cpu/native_tet_flow_darcy_chain",
				case, output], dry_run)
			return
		if domain.get("backend") == "idealized_cube_dual_tree":
			if waveforms or every_steps or time.get("mode") not in (None, "steady"):
				raise ValueError("idealized dual-tree flow uses its referenced steady case")
			case = domain_case(base, domain)
			species = None
			if domain.get("species_case_file"):
				species = path_from(base, domain["species_case_file"])
				flow_reference = path_from(species.parent,
					read_json(species)["flow_case_file"])
				if flow_reference != case:
					raise ValueError("dual-tree species case references a different flow case")
			if not dry_run:
				output.mkdir(parents=True, exist_ok=True)
			flow = output / "flow"
			command = [sys.executable,
				ROOT / "scripts/run_idealized_cube_dual_tree_fixed_flow.py",
				flow, "--case", case, "--ranks", str(ranks)]
			run(command, dry_run)
			if species:
				run([sys.executable,
					ROOT / "scripts/run_idealized_cube_dual_tree_oxygen.py",
					flow, output / "species", "--case", species,
					"--split-directory", flow / "submeshes",
					"--ranks", str(ranks)], dry_run)
			return
		if waveforms:
			raise ValueError("native FEM cases do not accept outer waveforms")
		if every_steps not in (None, 1):
			raise ValueError("native FEM writes every accepted step")
		workflow_file = domain.get("workflow_file")
		if workflow_file:
			workflow_path = path_from(base, workflow_file)
			workflow = read_json(workflow_path)
			if time and ("dt_s" in time or "steps" in time):
				case_time = read_json(path_from(workflow_path.parent,
					workflow["case_file"])).get("time", {})
				if any(case_time.get(key) != time[key]
						for key in ("dt_s", "steps") if key in time):
					raise ValueError("FEM workflow time must match its referenced case")
			for key in ("input_file", "case_file"):
				workflow[key] = str(path_from(workflow_path.parent, workflow[key]))
			workflow["output_directory"] = str(output)
			workflow["mpi_ranks"] = ranks
			if workflow.get("mesher", {}).get("executable"):
				workflow["mesher"]["executable"] = str(path_from(workflow_path.parent,
					workflow["mesher"]["executable"]))
			if dry_run:
				print(json.dumps(workflow, indent=2))
				run([sys.executable, ROOT / "scripts/run_native_tet_workflow.py",
					"<generated-workflow.json>"] + args.workflow_options, True)
				return
			with tempfile.TemporaryDirectory(prefix="tubular-solver-") as temp:
				effective = Path(temp) / "workflow.json"
				effective.write_text(json.dumps(workflow, sort_keys=True)+"\n", encoding="utf-8")
				run([sys.executable, ROOT / "scripts/run_native_tet_workflow.py", effective]
					+ args.workflow_options, False)
			if workflow["backend"] != "native_tet_p1_darcy_steady":
				case_time = read_json(Path(workflow["case_file"]))["time"]
				write_pvd(output / "fields", output / "series.pvd", case_time["dt_s"])
			return
		backend = domain.get("backend")
		if backend not in ("native_tet_hydraulic_graph", "native_tet_species_transport",
				"native_tet_darcy"):
			raise ValueError("FEM case_file requires a native_tet_* backend")
		case = domain_case(base, domain)
		case_data = read_json(case)
		binary = ROOT / "solvers/cpu" / backend
		if not dry_run:
			output.mkdir(parents=True, exist_ok=True)
		if time and "dt_s" in time:
			effective_data = case_data.copy()
			effective_data["mesh_file"] = str(path_from(case.parent,
				effective_data["mesh_file"]))
			effective_data["time"] = {"dt_s": time["dt_s"], "steps": time["steps"]}
			effective_case = output / "case.json"
			if not dry_run:
				effective_case.write_text(json.dumps(effective_data, indent=2)+"\n",
					encoding="utf-8")
			case = effective_case
		command = ["mpiexec", "-np", str(ranks), binary, case]
		if backend != "native_tet_darcy":
			run(command+["--check-input"], dry_run)
			if backend == "native_tet_hydraulic_graph" and not case_data.get("species"):
				command += ["--checkpoint-dir", output / "checkpoints"]
		else:
			run(command+["--check-input"], dry_run)
		run(command+["--output-dir", output / "fields"] + domain.get("options", []), dry_run)
		if not dry_run and backend != "native_tet_darcy":
			write_pvd(output / "fields", output / "series.pvd", read_json(case)["time"]["dt_s"])
		return
	if domain["method"] == "IGA":
		if device == "GPU" and ranks != 1:
			raise ValueError("IGA CUDA currently uses one GPU process")
		case_dir = path_from(base, domain["case_dir"])
		database = domain.get("database")
		if domain.get("prepare"):
			generated = output / "generated"
			prepare = [ROOT / "scripts/generate_case.sh", case_dir,
				"--output", generated, "--ranks", str(ranks)]
			if domain.get("template_directory"):
				prepare += ["--template-dir", path_from(base, domain["template_directory"])]
			run(prepare, dry_run)
			database = generated / "database" / f"{case_dir.name}-{ranks}.ntiga"
			case_dir = generated / "preprocessing"
		elif not database:
			raise ValueError("IGA requires a packed database path")
		database = path_from(base, database)
		case_dir = runtime_case(case_dir, output, time, waveforms, every_steps, dry_run)
		physics = domain.get("physics", ["flow"])
		if physics == ["flow"]:
			parallel_vtk = device == "CPU" and ranks > 1
			binary = ROOT / "solvers" / ("cuda" if device == "GPU" else "cpu") / (
				"iga_cuda" if device == "GPU" else "iga_navier_stokes")
			command = [binary]
			if device == "GPU":
				command += ["navier-stokes"]
			command += [database, case_dir, "--output",
				output / ("flow.pvtu" if parallel_vtk else "flow.vtu")]
			if parallel_vtk and "--visualization-format" not in domain.get("options", []):
				command += ["--visualization-format", "pvtu"]
			if every_steps:
				command += ["--output-every", str(every_steps)]
		elif physics == ["species"]:
			binary = ROOT / "solvers" / ("cuda" if device == "GPU" else "cpu") / (
				"iga_cuda" if device == "GPU" else "iga_solve")
			command = [binary]
			if device == "GPU":
				command += ["transport"]
			command += [database, case_dir]
		else:
			raise ValueError("IGA domain physics must select one existing solver")
		command += domain.get("options", [])
		if not dry_run:
			output.mkdir(parents=True, exist_ok=True)
		run((["mpiexec", "-np", str(ranks)] if device == "CPU" else [])+command,
			dry_run)
		return
	raise ValueError("unsupported domain")


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("configuration", type=Path)
	parser.add_argument("--dry-run", action="store_true")
	parser.add_argument("--database", help="database for direct legacy IGA config")
	parser.add_argument("--device", choices=("CPU", "GPU"), default="CPU")
	parser.add_argument("--physics", nargs="+", default=["flow"])
	parser.add_argument("--backend", help="backend for direct native FEM case")
	parser.add_argument("--resume", action="store_true", help="resume a FEM workflow")
	parser.add_argument("--stop-after-step", type=int)
	args = parser.parse_args()
	args.workflow_options = []
	if args.resume:
		args.workflow_options.append("--resume")
	if args.stop_after_step:
		args.workflow_options += ["--stop-after-step", str(args.stop_after_step)]
	configuration = args.configuration.resolve()
	data = read_json(configuration)
	if data.get("schema_version") == "solver-v1":
		domains = data["domains"]
		connections = data.get("connections", [])
		output = path_from(configuration.parent,
			data.get("output", {}).get("directory", "results"))
		ranks = data.get("resources", {}).get("mpi_ranks", 1)
		time = data.get("time", {})
		waveforms = data.get("waveforms", [])
		every_steps = data.get("output", {}).get("every_steps")
		if connections and not any(domain["dimension"] == "1d" for domain in domains) \
				and (waveforms or every_steps not in (None, 1)):
			raise ValueError("connected FEM cases do not accept outer waveforms or output frequency")
		validate_domains(domains, connections)
	else:
		domains = [legacy_domain(configuration, data, args)]
		connections = []
		output = configuration.parent / "results"
		ranks = 1
		time, waveforms, every_steps = {}, [], None
	if type(ranks) is not int or not 1 <= ranks <= 12:
		raise ValueError("resources.mpi_ranks must be between 1 and 12")
	if len(domains) == 1 and not connections:
		run_one(configuration.parent, domains[0], output, args.dry_run, ranks, args,
			time, waveforms, every_steps)
	elif len(domains) > 1 and connections:
		run_graph(configuration.parent, domains, connections,
			data["time"], waveforms, every_steps, output, args.dry_run, ranks, args)
	else:
		raise ValueError("multiple independent domains require separate configs")


if __name__ == "__main__":
	main()
