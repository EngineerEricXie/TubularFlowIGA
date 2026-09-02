#!/usr/bin/env bash

set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
ranks=2
output_argument=
clean=false
legacy_vtk=false
allow_preflight_failure=false

usage()
{
	printf 'usage: %s CASE_DIR [--output DIR] [--ranks N] [--clean] [--legacy-vtk] [--allow-preflight-failure]\n' "$0" >&2
}

if [[ $# -lt 1 ]]; then
	usage
	exit 2
fi
case_argument=$1
shift
while (( $# > 0 )); do
	case $1 in
		--output)
			if (( $# < 2 )); then usage; exit 2; fi
			output_argument=$2
			shift 2
			;;
		--ranks)
			if (( $# < 2 )); then usage; exit 2; fi
			ranks=$2
			shift 2
			;;
		--clean)
			clean=true
			shift
			;;
		--legacy-vtk)
			legacy_vtk=true
			shift
			;;
		--allow-preflight-failure)
			allow_preflight_failure=true
			shift
			;;
		*)
			printf 'unknown option: %s\n' "$1" >&2
			usage
			exit 2
			;;
	esac
done

if [[ ! -d $case_argument ]]; then
	printf 'case directory does not exist: %s\n' "$case_argument" >&2
	exit 2
fi
if ! [[ $ranks =~ ^[0-9]+$ ]] || (( ranks < 2 )); then
	printf 'ranks must be an integer of at least 2 because mpmetis rejects one partition\n' >&2
	exit 2
fi

case_dir=$(realpath "$case_argument")
case_name=${case_dir##*/}
if [[ ! $case_name =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ ]]; then
	printf 'case directory name must contain only letters, digits, underscores, or hyphens: %s\n' "$case_name" >&2
	exit 2
fi
if [[ ! -f $case_dir/simulation_config.json ]]; then
	printf 'case is missing simulation_config.json: %s\n' "$case_dir" >&2
	exit 2
fi

shopt -s nullglob
geometry_sources=("$case_dir"/*.swc "$case_dir"/*.obj)
shopt -u nullglob
if (( ${#geometry_sources[@]} == 0 )); then
	printf 'case has no SWC or OBJ geometry file: %s\n' "$case_dir" >&2
	exit 2
fi

if [[ -z $output_argument ]]; then
	output_argument=$case_dir/generated
fi
if [[ -L $output_argument ]]; then
	printf 'refusing to use a symbolic link as the generated-output directory: %s\n' "$output_argument" >&2
	exit 2
fi
output_dir=$(realpath -m "$output_argument")
if [[ $output_dir == / || $output_dir == "$repo_dir" || $output_dir == "$case_dir" ]]; then
	printf 'unsafe generated-output directory: %s\n' "$output_dir" >&2
	exit 2
fi

if [[ -e $output_dir && ! -d $output_dir ]]; then
	printf 'generated-output path is not a directory: %s\n' "$output_dir" >&2
	exit 2
fi
if [[ -d $output_dir ]]; then
	if [[ $clean != true && -n $(find "$output_dir" -mindepth 1 -maxdepth 1 -print -quit) ]]; then
		printf 'generated-output directory is not empty; rerun with --clean to replace it: %s\n' "$output_dir" >&2
		exit 2
	fi
fi

output_parent=${output_dir%/*}
output_base=${output_dir##*/}
mkdir -p "$output_parent"
staging_dir=$output_parent/.${output_base}.staging.$BASHPID
backup_dir=
if [[ -e $staging_dir ]]; then
	printf 'staging path already exists: %s\n' "$staging_dir" >&2
	exit 1
fi
cleanup()
{
	if [[ -n ${staging_dir:-} && -d $staging_dir ]]; then
		rm -rf -- "$staging_dir"
	fi
	if [[ -n ${backup_dir:-} && -e $backup_dir && ! -e $output_dir ]]; then
		mv "$backup_dir" "$output_dir"
	fi
}
trap cleanup EXIT

preprocessing_dir=$staging_dir/preprocessing
database_dir=$staging_dir/database
visualization_dir=$staging_dir/visualization
results_dir=$staging_dir/results
mkdir -p "$preprocessing_dir" "$database_dir" "$visualization_dir" "$results_dir"

"$repo_dir/scripts/check_dependencies.sh" preprocessing
"$repo_dir/scripts/check_dependencies.sh" vtkhdf
cp "$case_dir/simulation_config.json" "$preprocessing_dir/"
if [[ -f $case_dir/mesh_parameter.txt ]]; then
	cp "$case_dir/mesh_parameter.txt" "$preprocessing_dir/"
fi
cp "${geometry_sources[@]}" "$preprocessing_dir/"

make -C "$repo_dir" mesh spline cpu
pipeline_command=(
	"$repo_dir/preprocessing/mesh/tubular_mesh" pipeline
	"$preprocessing_dir" "$repo_dir/meshgeneration/template"
)
if [[ $allow_preflight_failure == true ]]; then
	pipeline_command+=(--allow-preflight-failure)
fi
"${pipeline_command[@]}"
for generated in skeleton_normalized.swc skeleton.vtp skeleton_smooth.swc \
	mesh_diagnostics.json skeleton_diagnostics.vtp controlmesh.vtk mesh_quality.json \
	initial_velocityfield.txt; do
	if [[ ! -s $preprocessing_dir/$generated ]]; then
		printf 'mesh preprocessing did not create required output: %s/%s\n' "$preprocessing_dir" "$generated" >&2
		exit 1
	fi
done

spline_command=(
	"$repo_dir/preprocessing/spline/spline" "$preprocessing_dir/" --no-legacy-text
)
if [[ $legacy_vtk == true ]]; then
	spline_command+=(--legacy-vtk)
fi
OMP_NUM_THREADS=${OMP_NUM_THREADS:-2} "${spline_command[@]}"
for generated in bzmeshinfo.txt spline_cache.igacache geometry_transform.json; do
	if [[ ! -s $preprocessing_dir/$generated ]]; then
		printf 'spline preprocessing did not create required output: %s/%s\n' "$preprocessing_dir" "$generated" >&2
		exit 1
	fi
done
if [[ $legacy_vtk == true ]]; then
	mv "$preprocessing_dir/bzmesh.vtk" "$visualization_dir/bzmesh.vtk"
fi

mpmetis "$preprocessing_dir/bzmeshinfo.txt" "$ranks"
database=$database_dir/$case_name-$ranks.ntiga
"$repo_dir/solvers/cpu/iga_pack" "$preprocessing_dir" "$ranks" "$database"
bezier_visualization=$visualization_dir/bzmesh.vtkhdf
"$repo_dir/solvers/cpu/iga_bezier_export" "$database" "$bezier_visualization"
for generated in bzmesh.vtkhdf bzmesh.bezier_geometry.json; do
	if [[ ! -s $visualization_dir/$generated ]]; then
		printf 'Bezier visualization export did not create required output: %s/%s\n' \
			"$visualization_dir" "$generated" >&2
		exit 1
	fi
done

"$repo_dir/solvers/cpu/iga_inspect" "$database"
"$repo_dir/solvers/cpu/iga_config_check" "$preprocessing_dir/simulation_config.json"
"$repo_dir/solvers/cpu/iga_case_check" "$database" "$preprocessing_dir"

mesh_points=$(sed -n 's/^[[:space:]]*"points":[[:space:]]*\([0-9][0-9]*\),*$/\1/p' "$preprocessing_dir/mesh_quality.json")
mesh_elements=$(sed -n 's/^[[:space:]]*"elements":[[:space:]]*\([0-9][0-9]*\),*$/\1/p' "$preprocessing_dir/mesh_quality.json")
minimum_scaled_jacobian=$(sed -n 's/^[[:space:]]*"minimum_scaled_jacobian":[[:space:]]*\([^,]*\),*$/\1/p' "$preprocessing_dir/mesh_quality.json")
geometry_preflight_valid=$(sed -n 's/^[[:space:]]*"valid":[[:space:]]*\(true\|false\),*$/\1/p' "$preprocessing_dir/mesh_diagnostics.json")
bezier_geometry_valid=$(sed -n 's/^[[:space:]]*"valid":[[:space:]]*\(true\|false\),*$/\1/p' "$visualization_dir/bzmesh.bezier_geometry.json")
bezier_unique_points=$(sed -n 's/^[[:space:]]*"unique_points":[[:space:]]*\([0-9][0-9]*\),*$/\1/p' "$visualization_dir/bzmesh.bezier_geometry.json")
for value in mesh_points mesh_elements minimum_scaled_jacobian geometry_preflight_valid \
	bezier_geometry_valid bezier_unique_points; do
	if [[ -z ${!value} ]]; then
		printf 'could not extract %s while writing manifest\n' "$value" >&2
		exit 1
	fi
done

{
	printf '{\n'
	printf '  "schema_version": 1,\n'
	printf '  "case_name": "%s",\n' "$case_name"
	printf '  "ranks": %s,\n' "$ranks"
	printf '  "geometry_preflight_valid": %s,\n' "$geometry_preflight_valid"
	printf '  "allow_preflight_failure": %s,\n' "$allow_preflight_failure"
	printf '  "mesh_points": %s,\n' "$mesh_points"
	printf '  "mesh_elements": %s,\n' "$mesh_elements"
	printf '  "minimum_scaled_jacobian": %s,\n' "$minimum_scaled_jacobian"
	printf '  "bezier_geometry_valid": %s,\n' "$bezier_geometry_valid"
	printf '  "bezier_unique_points": %s,\n' "$bezier_unique_points"
	printf '  "coordinate_spaces": {\n'
	printf '    "control_mesh": "source",\n'
	printf '    "database": "normalized",\n'
	printf '    "visualization": "source",\n'
	printf '    "transform": "preprocessing/geometry_transform.json"\n'
	printf '  },\n'
	printf '  "legacy_vtk": %s,\n' "$legacy_vtk"
	printf '  "paths": {\n'
	printf '    "runtime_case": "preprocessing",\n'
	printf '    "database": "database/%s-%s.ntiga",\n' "$case_name" "$ranks"
	printf '    "bezier_visualization": "visualization/bzmesh.vtkhdf",\n'
	printf '    "results": "results"\n'
	printf '  }\n'
	printf '}\n'
} > "$staging_dir/manifest.json"

if [[ -e $output_dir ]]; then
	backup_dir=$output_parent/.${output_base}.previous.$BASHPID
	mv "$output_dir" "$backup_dir"
fi
if mv "$staging_dir" "$output_dir"; then
	staging_dir=
else
	if [[ -n $backup_dir && -e $backup_dir ]]; then
		mv "$backup_dir" "$output_dir"
		backup_dir=
	fi
	exit 1
fi
if [[ -n $backup_dir ]]; then
	rm -rf -- "$backup_dir"
	backup_dir=
fi
trap - EXIT

final_preprocessing=$output_dir/preprocessing
final_database=$output_dir/database/$case_name-$ranks.ntiga
final_visualization=$output_dir/visualization/bzmesh.vtkhdf
final_results=$output_dir/results
printf 'case generation passed\ncase: %s\ngenerated root: %s\nruntime case: %s\ndatabase: %s\nBezier visualization: %s\nmanifest: %s\nresults directory: %s\nranks: %s\n' \
	"$case_dir" "$output_dir" "$final_preprocessing" "$final_database" \
	"$final_visualization" "$output_dir/manifest.json" "$final_results" "$ranks"
if grep -q '"neuron_transport"' "$case_dir/simulation_config.json"; then
	printf 'CPU neuron transport:\n  mpiexec -np %s %s/solvers/cpu/iga_solve %s %s --system neuron_transport --output %s/neuron-cpu.txt\n' \
		"$ranks" "$repo_dir" "$final_database" "$final_preprocessing" "$final_results"
	printf 'CUDA neuron transport:\n  %s/solvers/cuda/iga_cuda solve %s %s --system neuron_transport --output %s/neuron-cuda.txt\n' \
		"$repo_dir" "$final_database" "$final_preprocessing" "$final_results"
else
	printf 'CPU vascular flow:\n  mpiexec -np %s %s/solvers/cpu/iga_navier_stokes %s %s --output %s/velocity-cpu.txt\n' \
		"$ranks" "$repo_dir" "$final_database" "$final_preprocessing" "$final_results"
	printf 'CUDA vascular flow:\n  %s/solvers/cuda/iga_cuda navier-stokes %s %s --output %s/velocity-cuda.txt\n' \
		"$repo_dir" "$final_database" "$final_preprocessing" "$final_results"
fi
