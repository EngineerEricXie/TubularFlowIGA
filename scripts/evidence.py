#!/usr/bin/env python3
"""Check evidence cards and record the test runs that support them.

An evidence card (benchmarks/*.json) is a contract: what was tested, the
expected values and tolerances, and what may not be claimed. Its validator
checks the card itself. The tests that back the card are rerun by its test
command, and each run writes a record with the commit, tree state, host, and
toolchain, so a result is tied to exact sources without listing file hashes
in the card.

  evidence.py list                       cards, validators, and test commands
  evidence.py check [CARD...]            run the validators (no solver runs)
  evidence.py check --require-run ...    also require a passing run of HEAD
  evidence.py run [CARD...]              run the test commands and record them

Runs are written to artifacts/evidence/<card>.json (log beside it). CARD is a
card name such as t7_native_hydraulic_cli_evidence; no CARD means all cards.
"""
import argparse
import datetime
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT/"artifacts/evidence"

# card -> (validator, test command used when the card has no test_command)
CARDS = {
	"supported_input_contract": ("validate_supported_input_contract.py", "make contract-audit"),
	"t1_geometry_inventory": ("validate_t1_geometry_inventory.py", "make t1-audit"),
	"t1_mesher_comparison": ("validate_t1_mesher_comparison.py", "make t1-audit"),
	"t2_fixed_flow_contract": ("validate_t2_fixed_flow_contract.py", "make t2-contract-audit"),
	"t2_native_functional_evidence": ("validate_t2_native_functional_evidence.py", "make t2-contract-audit"),
	"t2_native_bifurcation_evidence": ("validate_t2_native_bifurcation_evidence.py", "make t2-contract-audit"),
	"t2_native_spatial_evidence": ("validate_t2_native_spatial_evidence.py", "make t2-contract-audit"),
	"t2_native_temporal_evidence": ("validate_t2_native_temporal_evidence.py", "make t2-contract-audit"),
	"t3_native_tet_solid_contract": ("validate_t3_native_tet_solid_contract.py", "make t3-solid-contract-audit"),
	"t4_native_ale_contract": ("validate_t4_native_ale_contract.py", "make t4-contract-audit"),
	"t4_native_ale_reference_evidence": ("validate_t4_native_ale_evidence.py", "make t4-contract-audit"),
	"t5_native_matching_interface_contract": ("validate_t5_native_matching_interface_contract.py",
		"make t5-matching-interface-audit"),
	"t6_native_lv_qoi_evidence": ("validate_t6_native_lv_qoi_evidence.py", "make -C solvers/cpu lv-backflow-test"),
	"t6_native_tet_darcy_evidence": ("validate_t6_native_tet_darcy_evidence.py", None),
	"t6_native_wall_reservoir_evidence": ("validate_t6_native_wall_reservoir_evidence.py", None),
	"t7_native_0d_graph_evidence": ("validate_t7_native_0d_graph_evidence.py", None),
	"t7_native_hydraulic_cli_evidence": ("validate_t7_native_hydraulic_cli_evidence.py", None),
	"t7_native_moving_species_evidence": ("validate_t7_native_moving_species_evidence.py", None),
	"t9_native_species_cli_evidence": ("validate_t9_native_species_cli_evidence.py", None),
	"liver_roi_functional_evidence": ("validate_liver_roi_functional_evidence.py", "make t6-native-tet-darcy-test"),
}


def card_path(name):
	return ROOT/"benchmarks"/(name+".json")


def test_command(name):
	command = json.loads(card_path(name).read_text()).get("test_command") or CARDS[name][1]
	if not command:
		raise SystemExit(f"{name}: card has no test_command")
	return command


def sha256(path):
	return hashlib.sha256(path.read_bytes()).hexdigest()


def output(command):
	try:
		result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, timeout=30)
	except (OSError, subprocess.TimeoutExpired):
		return None
	text = (result.stdout or result.stderr).strip()
	return text.splitlines()[0] if text else None


def source_state():
	commit = output(["git", "rev-parse", "HEAD"])
	changes = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"],
		cwd=ROOT, capture_output=True, text=True).stdout.strip()
	return commit, bool(changes)


def petsc_version():
	header = Path(os.environ.get("PETSC_DIR", "/nonexistent"))/"include/petscversion.h"
	if not header.exists():
		return None
	text = header.read_text()
	parts = [re.search(rf"#define PETSC_VERSION_{key}\s+(\d+)", text) for key in ("MAJOR", "MINOR", "SUBMINOR")]
	if not all(parts):
		return None
	return ".".join(part.group(1) for part in parts)+" ("+os.environ.get("PETSC_ARCH", "")+")"


def toolchain():
	return {
		"cxx": output([os.environ.get("CXX", "c++"), "--version"]),
		"mpiexec": output(["mpiexec", "--version"]),
		"petsc": petsc_version(),
		"python": platform.python_version(),
	}


def run_cards(names):
	RUNS.mkdir(parents=True, exist_ok=True)
	groups = {}
	for name in names:
		groups.setdefault(test_command(name), []).append(name)
	failed = 0
	for command, members in groups.items():
		commit, dirty = source_state()
		log = RUNS/(members[0]+".log")
		print(f"+ {command}  ({', '.join(members)})", flush=True)
		started = time.time()
		with log.open("w") as stream:
			returncode = subprocess.run(command, shell=True, cwd=ROOT, stdout=stream,
				stderr=subprocess.STDOUT).returncode
		record = {
			"schema_version": 1,
			"test_command": command,
			"returncode": returncode,
			"passed": returncode == 0,
			"started_utc": datetime.datetime.fromtimestamp(started, datetime.timezone.utc).isoformat(timespec="seconds"),
			"elapsed_s": round(time.time()-started, 1),
			"commit": commit,
			"uncommitted_changes": dirty,
			"host": platform.node(),
			"platform": platform.platform(),
			"toolchain": toolchain(),
			"log": log.relative_to(ROOT).as_posix(),
		}
		for name in members:
			card = dict(record, card=card_path(name).relative_to(ROOT).as_posix(), card_sha256=sha256(card_path(name)))
			(RUNS/(name+".json")).write_text(json.dumps(card, indent=2)+"\n")
		print(f"  {'passed' if returncode == 0 else f'FAILED (exit {returncode})'} in {record['elapsed_s']} s; log {record['log']}")
		failed += returncode != 0
	return 1 if failed else 0


def check_cards(names, require_run):
	commit, dirty = source_state()
	failures = []
	for name in names:
		validator = ROOT/"scripts"/CARDS[name][0]
		result = subprocess.run([sys.executable, str(validator)], cwd=ROOT, capture_output=True, text=True)
		problem = None
		if result.returncode:
			lines = (result.stderr or result.stdout).strip().splitlines()
			problem = lines[-1] if lines else f"validator exited {result.returncode}"
		elif require_run:
			record_path = RUNS/(name+".json")
			record = json.loads(record_path.read_text()) if record_path.exists() else None
			if record is None:
				problem = "no recorded run; use evidence.py run "+name
			elif not record.get("passed"):
				problem = "recorded run failed; see "+record.get("log", "its log")
			elif record.get("commit") != commit or record.get("uncommitted_changes") or dirty:
				problem = "recorded run is not of the clean HEAD commit"
			elif record.get("card_sha256") != sha256(card_path(name)):
				problem = "card changed after the recorded run"
		print(f"{name}: {'PASS' if problem is None else 'FAIL: '+problem}")
		if problem:
			failures.append(name)
	print(f"{len(names)-len(failures)}/{len(names)} cards pass")
	return 1 if failures else 0


def main():
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("action", choices=("list", "check", "run"))
	parser.add_argument("cards", nargs="*", help="card names (default: all)")
	parser.add_argument("--require-run", action="store_true", help="check: require a passing run of HEAD")
	args = parser.parse_args()
	unknown = [name for name in args.cards if name not in CARDS]
	if unknown:
		parser.error("unknown card(s): "+", ".join(unknown)+"; see evidence.py list")
	names = args.cards or sorted(CARDS)
	if args.action == "list":
		for name in names:
			print(f"{name}\n  validator: scripts/{CARDS[name][0]}\n  tests:     {test_command(name)}")
		return 0
	if args.action == "run":
		return run_cards(names)
	return check_cards(names, args.require_run)


if __name__ == "__main__":
	sys.exit(main())
