#!/bin/bash
set -euo pipefail
while test $# -gt 0; do
    case "$1" in
        --root) root=$2; shift 2;;
        --variants) variants=$2; shift 2;;
        --case-list) cases=$2; shift 2;;
        *) exit 2;;
    esac
done
: "${LOCAL:?Slurm node-local LOCAL required}" "${SLURM_JOB_ID:?Slurm required}"
repo=/ocean/projects/mch260002p/thsieh1/TubularFlowIGA
work="$LOCAL/fsi-time-opt-$SLURM_JOB_ID"
out="$root/job-$SLURM_JOB_ID"
mkdir -p "$work" "$out"
started=$(date -u +%FT%TZ)
stage=staging
failure_status=failed_environment
finish() {
    result=$?
    trap - EXIT
    printf '{"schema_version":1,"job_id":"%s","status":"%s","exit_code":%s,"stage":"%s","started":"%s","ended":"%s"}\n' "$SLURM_JOB_ID" "$([ "$result" -eq 0 ] && echo passed || echo "$failure_status")" "$result" "$stage" "$started" "$(date -u +%FT%TZ)" > "$work/status.json"
    copy_status=0
    cp -a "$work/evidence/." "$out/" || copy_status=$?
    cp "$work/status.json" "$out/" || copy_status=$?
    if test "$copy_status" -ne 0; then echo "copy-back failed ($copy_status); local=$work" >&2; test "$result" -ne 0 || result=$copy_status; fi
    exit "$result"
}
mkdir -p "$work/evidence"
trap finish EXIT
cp "$cases" "$work/evidence/case.json"
sha256sum "$cases" "$0" > "$work/evidence/runner-input.sha256"
sha256sum -c "$root/surface.sha256"
cp "$root/surface.txt" "$work/surface.txt"
sha256sum "$work/surface.txt" > "$work/evidence/input.sha256"
{
    date -u +%FT%TZ; module list; mpicxx --version
    lscpu; command -v numactl >/dev/null && numactl --hardware || true
    scontrol show job "$SLURM_JOB_ID"
    grep Cpus_allowed_list /proc/self/status
    for name in OMP_NUM_THREADS OMP_THREAD_LIMIT OMP_DYNAMIC OMP_PROC_BIND OMP_PLACES IGA_ASSEMBLY_THREADS IGA_ASSEMBLY_BATCH_SIZE OPENBLAS_NUM_THREADS MKL_NUM_THREADS BLIS_NUM_THREADS PETSC_DIR PETSC_ARCH; do printf '%s=%s\n' "$name" "${!name:-}"; done
} > "$work/evidence/environment.txt" 2>&1
cp "$PETSC_DIR/$PETSC_ARCH/lib/petsc/conf/petscvariables" "$work/evidence/petscvariables.txt"
IFS=, read -r -a variant_list <<< "$variants"
for variant in "${variant_list[@]}"; do
    sha256sum -c "$root/$variant-source.sha256"
    mkdir -p "$work/$variant"
    tar -xf "$root/$variant-source.tar" -C "$work/$variant"
    cp "$root/$variant-source.sha256" "$work/evidence/$variant-source.sha256"
done
workload_variant=${variant_list[${#variant_list[@]}-1]}
workload_src="$work/$workload_variant"
sha256sum -c "$root/archived-workload.sha256"
tar -xOf "$root/archived-workload.tar.gz" ./driver.cpp > "$workload_src/benchmarks/fsi-time-opt/driver.cpp"
sha256sum "$workload_src/benchmarks/fsi-time-opt/driver.cpp" "$workload_src/benchmarks/fsi-time-opt/benchmark.cpp" > "$work/evidence/workload-source.sha256"
for variant in "${variant_list[@]}"; do
    stage="build-$variant"
    src="$work/$variant"
    include_flags=(-I"$src/solvers/cpu/include" -I"$src/include" -I"$src/solvers/cpu/tests")
    # PETSc make variables are expanded by its own make configuration, not guessed.
    make -C "$src/solvers/cpu" -j2 parallel_element_batch_test parallel_element_batch_openmp_test element_assembly_execution_openmp_test > "$work/evidence/$variant-build.stdout" 2> "$work/evidence/$variant-build.stderr"
    /usr/bin/time -v -o "$work/evidence/$variant-build.time" make -f "$workload_src/benchmarks/fsi-time-opt/runner.mk" SOURCE="$src" WORKLOAD_SOURCE="$workload_src" BUILD="$work/$variant" benchmark volume >> "$work/evidence/$variant-build.stdout" 2>> "$work/evidence/$variant-build.stderr"
    sha256sum "$work/$variant/benchmark" "$work/$variant/volume" > "$work/evidence/$variant-binary.sha256"
    ldd "$work/$variant/benchmark" > "$work/evidence/$variant-ldd.txt"
    cp "$work/$variant/benchmark" "$work/evidence/$variant-benchmark"
    cp "$work/$variant/volume" "$work/evidence/$variant-volume"
done
variant=${variant_list[${#variant_list[@]}-1]}
src="$work/$variant"
stage=W0
failure_status=failed_numerical
for target in parallel_element_batch_test parallel_element_batch_openmp_test element_assembly_execution_openmp_test; do
    /usr/bin/time -v -o "$work/evidence/W0-$target.time" "$src/solvers/cpu/$target" > "$work/evidence/W0-$target.stdout" 2> "$work/evidence/W0-$target.stderr"
done
"$src/solvers/cpu/element_assembly_execution_openmp_test" --mpi-single > "$work/evidence/W0-mpi-single.stdout" 2> "$work/evidence/W0-mpi-single.stderr"
mkdir -p "$work/evidence/W0-volume"
cd "$work/evidence/W0-volume"
IGA_PROFILE_DETAIL=1 /usr/bin/time -v -o timing.txt mpiexec -np 1 --map-by slot:PE="$OMP_NUM_THREADS" --bind-to core "$src/volume" > stdout.txt 2> stderr.txt
stage=W1
for pair in 0 1 2; do
    order=0,1,2; test "$pair" -ne 1 || order=2,1,0
    IFS=, read -r -a modes <<< "$order"
    for mode in "${modes[@]}"; do
        variant=${variant_list[0]}; detail=0; label=$variant
        if test "$mode" -ne 0; then variant=$workload_variant; label=$variant-off; fi
        if test "$mode" -eq 2; then detail=1; label=$variant-on; fi
        sample="$work/evidence/W1-pair-$pair-$label"
        mkdir -p "$sample"; cd "$sample"
        IGA_PROFILE_DETAIL=$detail /usr/bin/time -v -o timing.txt mpiexec -np 1 --map-by slot:PE="$OMP_NUM_THREADS" --bind-to core "$work/$variant/benchmark" "$work/surface.txt" > stdout.txt 2> stderr.txt
    done
done
stage=verification
failure_status=inconclusive
python3 "$workload_src/benchmarks/fsi-time-opt/verify_p0.py" "$work/evidence" "$variants"
stage=complete
