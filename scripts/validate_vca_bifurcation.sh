#!/usr/bin/env bash

set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
relative_tolerance=${VCA_RELATIVE_TOLERANCE:-1e-6}
pressure_relative_tolerance=${VCA_PRESSURE_RELATIVE_TOLERANCE:-1e-4}
volume_balance_tolerance=${VCA_VOLUME_BALANCE_TOLERANCE:-1e-12}

if [[ $# -ne 1 ]]; then
	printf 'usage: %s EMPTY_WORK_DIRECTORY\n' "$0" >&2
	exit 2
fi

work_dir=$1
case_dir=$work_dir/preprocessing
database_dir=$work_dir/database
results_dir=$work_dir/results
solver=$repo_dir/solvers/cpu/iga_navier_stokes
packer=$repo_dir/solvers/cpu/iga_pack

if [[ ! -x $solver || ! -x $packer ]]; then
	printf 'build the PETSc CPU solver first: make cpu-petsc PETSC_DIR=... PETSC_ARCH=...\n' >&2
	exit 2
fi
for command in mpiexec awk cmp; do
	if ! command -v "$command" >/dev/null 2>&1; then
		printf 'required command is unavailable: %s\n' "$command" >&2
		exit 2
	fi
done

RunCase()
{
	local ranks=$1
	local database=$2
	local prefix=$3
	shift 3
	mkdir -p "$(dirname "$prefix")"
	mpiexec -np "$ranks" "$solver" "$database" "$case_dir" \
		--checkpoint "$prefix" --output "$(dirname "$prefix")/flow.txt" "$@"
}

RelativeL2()
{
	local reference=$1
	local candidate=$2
	local label=$3
	local tolerance=$4
	awk -v tolerance="$tolerance" -v label="$label" '
		NR == FNR {
			for (column = 1; column <= NF; ++column) reference_value[FNR, column] = $column
			reference_columns[FNR] = NF
			reference_rows = FNR
			next
		}
		NF != reference_columns[FNR] {
			printf "%s column count differs at row %d\\n", label, FNR > "/dev/stderr"
			exit 2
		}
		{
			for (column = 1; column <= NF; ++column) {
				delta = $column-reference_value[FNR, column]
				error += delta*delta
				norm += reference_value[FNR, column]*reference_value[FNR, column]
			}
		}
		END {
			if (FNR != reference_rows) {
				printf "%s row count differs\\n", label > "/dev/stderr"
				exit 2
			}
			relative = sqrt(error)/(sqrt(norm) > 0 ? sqrt(norm) : 1)
			printf "%s relative_l2=%.17g (tolerance %.17g)\\n", label, relative, tolerance
			if (relative > tolerance) exit 1
		}' "$reference" "$candidate"
}

ExtractNumber()
{
	local file=$1
	local key=$2
	sed -nE 's/.*"'"$key"'": *([-+0-9.eE]+).*/\1/p' "$file" | tail -n 1
}

CompareReservoirNumber()
{
	local key=$1
	local first=$2
	local second=$3
	local left
	local right
	left=$(ExtractNumber "$first" "$key")
	right=$(ExtractNumber "$second" "$key")
	if [[ -z $left || -z $right ]]; then
		printf 'cannot read reservoir %s from VCA checkpoint metadata\n' "$key" >&2
		exit 1
	fi
	awk -v key="$key" -v left="$left" -v right="$right" -v tolerance="$relative_tolerance" '
		BEGIN {
			left_magnitude = left < 0 ? -left : left
			right_magnitude = right < 0 ? -right : right
			if (left_magnitude == 0 && right_magnitude == 0) relative = 0
			else {
				difference = left-right
				if (difference < 0) difference = -difference
				relative = difference/(left_magnitude > right_magnitude ? left_magnitude : right_magnitude)
			}
			printf "reservoir %s relative_difference=%.17g (tolerance %.17g)\\n", key, relative, tolerance
			exit relative > tolerance
		}'
}

ValidateManifest()
{
	local manifest=$1
	local records=$2
	if [[ ! -s $manifest ]]; then
		printf 'missing VCA manifest: %s\n' "$manifest" >&2
		exit 1
	fi
	if [[ $(grep -c '"device_source_mol_s"' "$manifest") -ne $records ]]; then
		printf 'unexpected VCA balance record count in %s\n' "$manifest" >&2
		exit 1
	fi
	grep -Fq '"outlet_ids": [2, 3]' "$manifest"
	grep -Fq '"oxygen"' "$manifest"
	awk -v tolerance="$volume_balance_tolerance" -v manifest="$manifest" '
		/"volume_change_m3":/ {
			value = $0
			sub(/^.*"volume_change_m3": */, "", value)
			sub(/[,}].*$/, "", value)
			magnitude = value+0
			if (magnitude < 0) magnitude = -magnitude
			if (magnitude > tolerance) {
				printf "reservoir volume change %.17g exceeds tolerance %.17g in %s\\n", magnitude, tolerance, manifest > "/dev/stderr"
				exit 1
			}
		}' "$manifest"
	if grep -Eqi '(^|[^[:alpha:]])(nan|inf)([^[:alpha:]]|$)' "$manifest"; then
		printf 'non-finite VCA manifest value in %s\n' "$manifest" >&2
		exit 1
	fi
}

"$repo_dir/scripts/generate_case.sh" "$repo_dir/examples/vascular_flow/vca_bifurcation" \
	--output "$work_dir" --ranks 2
awk '{print 0}' "$case_dir/bzmeshinfo.txt.epart.2" > "$case_dir/bzmeshinfo.txt.epart.1"
"$packer" "$case_dir" 1 "$database_dir/vca_bifurcation-1.ntiga"

one_database=$database_dir/vca_bifurcation-1.ntiga
two_database=$database_dir/vca_bifurcation-2.ntiga

RunCase 1 "$one_database" "$results_dir/one/full/checkpoint"
RunCase 1 "$one_database" "$results_dir/one/split/checkpoint" --stop-after-step 2
RunCase 1 "$one_database" "$results_dir/one/resumed/checkpoint" \
	--restart "$results_dir/one/split/checkpoint"
RunCase 2 "$two_database" "$results_dir/two/full/checkpoint"
RunCase 2 "$two_database" "$results_dir/two/split/checkpoint" --stop-after-step 2
RunCase 2 "$two_database" "$results_dir/two/resumed/checkpoint" \
	--restart "$results_dir/two/split/checkpoint"

for ranks in one two; do
	cmp "$results_dir/$ranks/full/checkpoint.state" "$results_dir/$ranks/resumed/checkpoint.state"
	cmp "$results_dir/$ranks/full/checkpoint.vca_transport.state" \
		"$results_dir/$ranks/resumed/checkpoint.vca_transport.state"
done

RelativeL2 "$results_dir/one/full/flow.txt" "$results_dir/two/full/flow.txt" velocity "$relative_tolerance"
RelativeL2 "$results_dir/one/full/flow.txt.pressure" "$results_dir/two/full/flow.txt.pressure" pressure "$pressure_relative_tolerance"
for key in volume_m3 temperature_c hematocrit_percent oxygen; do
	CompareReservoirNumber "$key" "$results_dir/one/full/checkpoint.vca.json" \
		"$results_dir/two/full/checkpoint.vca.json"
done
for ranks in one two; do
	ValidateManifest "$results_dir/$ranks/full/coupling_manifest.json" 4
	ValidateManifest "$results_dir/$ranks/split/coupling_manifest.json" 2
	ValidateManifest "$results_dir/$ranks/resumed/coupling_manifest.json" 2
done

printf 'VCA bifurcation one-/two-rank validation passed: %s\n' "$work_dir"
