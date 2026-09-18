#!/bin/bash
# Batch-only screen: one immutable binary, fixed 8 threads, exact ordered outputs.
set -euo pipefail
while test $# -gt 0; do
    case "$1" in --root) root=$2; shift 2;; --variants) variant=$2; shift 2;; --case-list) cases=$2; shift 2;; *) exit 2;; esac
done
: "${root:?}" "${variant:?}" "${cases:?}" "${LOCAL:?}" "${SLURM_JOB_ID:?}"
work="$LOCAL/fsi-p2-$SLURM_JOB_ID"; out="$root/job-$SLURM_JOB_ID"; mkdir -p "$work/evidence" "$out" "$work/$variant"
stage=staging; failure=failed_environment
save() { printf '{"job_id":"%s","stage":"%s"}\n' "$SLURM_JOB_ID" "$stage" > "$work/evidence/stage.json"; cp -au "$work/evidence/." "$out/"; }
finish() { result=$?; trap - EXIT TERM USR1; printf '{"job_id":"%s","stage":"%s","status":"%s","exit_code":%s}\n' "$SLURM_JOB_ID" "$stage" "$([ "$result" -eq 0 ] && echo passed || echo "$failure")" "$result" > "$work/evidence/status.json"; cp -au "$work/evidence/." "$out/" || { test "$result" -ne 0 || result=1; }; exit "$result"; }
trap finish EXIT; trap 'failure=interrupted; exit 124' TERM USR1
sha256sum -c "$root/$variant-source.sha256"; tar -xf "$root/$variant-source.tar" -C "$work/$variant"; cp "$root/$variant-source.sha256" "$work/evidence/"
src="$work/$variant"; tools="$src/benchmarks/fsi-time-opt"
sha256sum -c "$root/surface.sha256"; cp "$root/surface.txt" "$work/surface.txt"
sha256sum -c "$root/archived-workload.sha256"; tar -xOf "$root/archived-workload.tar.gz" ./driver.cpp > "$tools/driver.cpp"
cp "$0" "$work/evidence/runner.sh"; cp "$cases" "$work/evidence/case.json"; sha256sum "$tools/driver.cpp" "$tools/benchmark.cpp" "$root/surface.txt" "$cases" > "$work/evidence/input.sha256"
{ date -u +%FT%TZ; module list; mpicxx --version; lscpu; scontrol show job "$SLURM_JOB_ID"; grep Cpus_allowed_list /proc/self/status; } > "$work/evidence/environment.txt" 2>&1
stage=build; save
make -C "$src/solvers/cpu" -j2 parallel_element_batch_test parallel_element_batch_openmp_test element_assembly_execution_openmp_test > "$work/evidence/build.stdout" 2> "$work/evidence/build.stderr"
/usr/bin/time -v -o "$work/evidence/build.time" make -f "$tools/runner.mk" SOURCE="$src" WORKLOAD_SOURCE="$src" BUILD="$src" benchmark volume >> "$work/evidence/build.stdout" 2>> "$work/evidence/build.stderr"
sha256sum "$src/benchmark" "$src/volume" > "$work/evidence/binary.sha256"; cp "$src/benchmark" "$src/volume" "$work/evidence/"
stage=W0; failure=failed_numerical; save
for target in parallel_element_batch_test parallel_element_batch_openmp_test element_assembly_execution_openmp_test; do "$src/solvers/cpu/$target" > "$work/evidence/$target.stdout" 2> "$work/evidence/$target.stderr"; done
"$src/solvers/cpu/element_assembly_execution_openmp_test" --mpi-single > "$work/evidence/mpi-single.stdout" 2> "$work/evidence/mpi-single.stderr"
stage=W1; save
for order in 8,16,32 32,16,8; do
    IFS=, read -r -a batches <<< "$order"
    for batch in "${batches[@]}"; do
        process="$work/evidence/W1-batch-$batch-$(test "$order" = 8,16,32 && echo A || echo B)"; mkdir -p "$process"; cd "$process"
        IGA_ASSEMBLY_BATCH_SIZE=$batch IGA_PROFILE_DETAIL=1 /usr/bin/time -v -o timing.txt mpiexec -np 1 --map-by slot:PE="$OMP_NUM_THREADS" --bind-to core "$src/benchmark" "$work/surface.txt" > stdout.txt 2> stderr.txt
    done
done
stage=selection; save
python3 "$tools/verify_p2.py" --select "$work/evidence"
selected=$(cat "$work/evidence/selected-batch.txt")
stage="W0-volume-batch-$selected"; save
sample="$work/evidence/W0-volume"; mkdir -p "$sample"; cd "$sample"
IGA_ASSEMBLY_BATCH_SIZE=$selected IGA_PROFILE_DETAIL=1 /usr/bin/time -v -o timing.txt mpiexec -np 1 --map-by slot:PE="$OMP_NUM_THREADS" --bind-to core "$src/volume" > stdout.txt 2> stderr.txt
stage=verification; save
python3 "$tools/verify_p2.py" "$work/evidence"
stage=complete
