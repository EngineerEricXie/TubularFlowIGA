#!/usr/bin/env bash

set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
ranks=${RANKS:-2}

if [[ $# -lt 1 || $# -gt 2 ]]; then
	printf 'usage: %s EXAMPLE [EMPTY_WORK_DIRECTORY]\n' "$0" >&2
	exit 2
fi
if ! [[ $ranks =~ ^[0-9]+$ ]] || (( ranks < 2 )); then
	printf 'RANKS must be an integer of at least 2 because mpmetis rejects one partition\n' >&2
	exit 2
fi

example=$1
if [[ ! $example =~ ^(neuron_transport|vascular_flow)/[A-Za-z0-9][A-Za-z0-9_-]*$ ]]; then
	printf 'example must be APPLICATION/CASE under examples/neuron_transport or examples/vascular_flow: %s\n' "$example" >&2
	exit 2
fi
case_name=${example##*/}
source_dir=$repo_dir/examples/$example
if [[ ! -d $source_dir ]]; then
	printf 'unknown example: %s\navailable examples:\n' "$example" >&2
	find "$repo_dir/examples/neuron_transport" "$repo_dir/examples/vascular_flow" \
		-mindepth 1 -maxdepth 1 -type d -printf '  %h/%f\n' \
		| sed "s|$repo_dir/examples/|  |" | sort >&2
	exit 2
fi

if [[ $# -eq 2 ]]; then
	work_dir=$2
else
	work_dir=$(mktemp -d "${TMPDIR:-/tmp}/tubularflowiga-${case_name}.XXXXXX")
fi

exec "$repo_dir/scripts/generate_case.sh" "$source_dir" \
	--output "$work_dir" --ranks "$ranks"
