#!/usr/bin/env bash

set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
ranks=${RANKS:-2}

if ! [[ $ranks =~ ^[0-9]+$ ]] || (( ranks < 2 )); then
	printf 'RANKS must be an integer of at least 2 because mpmetis rejects one partition\n' >&2
	exit 2
fi

if [[ $# -gt 1 ]]; then
	printf 'usage: %s [EMPTY_WORK_DIRECTORY]\n' "$0" >&2
	exit 2
fi
if [[ $# -eq 1 ]]; then
	work_dir=$1
	mkdir -p "$work_dir"
	if [[ -n $(find "$work_dir" -mindepth 1 -maxdepth 1 -print -quit) ]]; then
		printf 'work directory must be empty: %s\n' "$work_dir" >&2
		exit 2
	fi
else
	work_dir=$(mktemp -d "${TMPDIR:-/tmp}/tubularflowiga-smoke.XXXXXX")
fi

"$repo_dir/scripts/check_dependencies.sh" preprocessing
cp "$repo_dir"/examples/smoke/input/* "$work_dir"/

make -C "$repo_dir" mesh spline cpu
"$repo_dir/preprocessing/mesh/tubular_mesh" pipeline \
	"$work_dir" "$repo_dir/meshgeneration/template"
OMP_NUM_THREADS=${OMP_NUM_THREADS:-2} \
	"$repo_dir/preprocessing/spline/spline" "$work_dir/" --no-legacy-text
mpmetis "$work_dir/bzmeshinfo.txt" "$ranks"
database=$work_dir/smoke-$ranks.ntiga
"$repo_dir/solvers/cpu/iga_pack" "$work_dir" "$ranks" "$database"
"$repo_dir/solvers/cpu/iga_inspect" "$database"
"$repo_dir/solvers/cpu/iga_case_check" "$database" "$work_dir"

printf 'smoke preparation passed\ncase directory: %s\ndatabase: %s\nranks: %s\n' \
	"$work_dir" "$database" "$ranks"
"$repo_dir/solvers/cpu/iga_config_check" "$work_dir/simulation_config.json"
