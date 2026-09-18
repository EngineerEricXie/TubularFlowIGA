#!/bin/bash
# Small independent compile/export check only; no FSI simulation.
set -euo pipefail
module load anaconda3 openmpi/4.0.5-gcc10.2.0
export PETSC_DIR=/ocean/projects/mch260002p/thsieh1/TubularFlowIGA-hpc-acceptance/petsc
export PETSC_ARCH=arch-bridges2-acceptance
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OMP_DYNAMIC=FALSE
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
root=${1:?evidence root}; tools_variant=${2:?frozen native exporter workload}
resume=${3:-}
baseline=${4:-P1ABefore}; candidate=${5:-P1A}
: "${LOCAL:?Slurm LOCAL required}" "${SLURM_JOB_ID:?Slurm required}"
work="$LOCAL/w2-export-smoke-$SLURM_JOB_ID"; out="$root/job-$SLURM_JOB_ID"
mkdir -p "$work/evidence/binaries" "$out"
stage=staging
finish() {
    result=$?; trap - EXIT
    printf '{"job_id":"%s","stage":"%s","status":"%s","exit_code":%s,"simulation_executed":false}\n' "$SLURM_JOB_ID" "$stage" "$([ "$result" -eq 0 ] && echo passed || echo failed)" "$result" > "$work/evidence/status.json"
    copied=0; cp -a "$work/evidence/." "$out/" || copied=$?
    test "$copied" -eq 0 || { test "$result" -ne 0 || result=$copied; }
    exit "$result"
}
trap finish EXIT
for variant in "$baseline" "$candidate" "$tools_variant"; do
    sha256sum -c "$root/$variant-source.sha256"; mkdir -p "$work/$variant"
    tar -xf "$root/$variant-source.tar" -C "$work/$variant"
    cp "$root/$variant-source.sha256" "$work/evidence/$variant-source.sha256"
done
common="$work/$tools_variant"
cp "$0" "$work/evidence/runner.sh"
previous=
if test -n "$resume"; then
    test "$baseline,$candidate" = P1ABefore,P1A
    previous="$root/job-$resume"
    python3 - "$previous" "$work/evidence" "$tools_variant" <<'PY'
import hashlib,json,pathlib,sys
old,new=map(pathlib.Path,sys.argv[1:3]); tools=sys.argv[3]
status=json.loads((old/'status.json').read_text())
assert status['stage']=='heartbeat-usage-P1ABefore' and status['exit_code']==127
for variant in ('P1ABefore','P1A',tools):
    assert (old/f'{variant}-source.sha256').read_text().split()[0]==(new/f'{variant}-source.sha256').read_text().split()[0]
for line in (old/'binary.sha256').read_text().splitlines():
    checksum,path=line.split(); name='P1ABefore-heartbeat' if pathlib.Path(path).name=='heartbeat' else 'native_reference_test'
    assert hashlib.sha256((old/'binaries'/name).read_bytes()).hexdigest()==checksum
assert 'usage: bifurcation_fsi' in (old/'P1ABefore-usage.stderr').read_text()
print('native_export_and_baseline_compile_prefix_sources_binaries_usage verified')
PY
    printf '%s\n' "$resume" > "$work/evidence/reused-prefix-job.txt"
fi
{ date -u +%FT%TZ; module list; mpicxx --version; scontrol show job "$SLURM_JOB_ID"; } > "$work/evidence/environment.txt" 2>&1
stage=native-export-test-build
if test -n "$previous"; then
    cp "$previous/binaries/native_reference_test" "$work/native_reference_test"
    cp -a "$previous/native-tests" "$work/evidence/"
else
    make -f "$common/benchmarks/fsi-time-opt/heartbeat_case.mk" SOURCE="$work/$candidate" WORKLOAD_SOURCE="$common" BUILD="$work" native_reference_test > "$work/evidence/native-build.stdout" 2> "$work/evidence/native-build.stderr"
fi
stage=native-export-test
test -n "$previous" || "$work/native_reference_test" "$work/evidence/native-tests"
cp "$work/native_reference_test" "$work/evidence/binaries/native_reference_test"
sha256sum "$work/native_reference_test" > "$work/evidence/binary.sha256"
python3 - "$work/evidence/native-tests" <<'PY'
import csv,hashlib,json,math,pathlib,sys
root=pathlib.Path(sys.argv[1])
for directory,gate in [('good',False),('positive',True)]:
    manifest=json.loads((root/directory/'manifest.json').read_text())
    assert manifest['native_gates_passed'] is gate and manifest['checkpoint'] is False
    for field in manifest['fields']:
        path=root/directory/field['file']
        assert hashlib.sha256(path.read_bytes()).hexdigest()==field['sha256']
        rows=list(csv.reader(path.open()))
        assert len(rows)==field['rows']+1
        assert [int(row[0]) for row in rows[1:]]==[2,7]
        assert all(len(row)==field['columns']+1 for row in rows)
        assert all(math.isfinite(float(v)) for row in rows[1:] for v in row[1:])
for directory in ('unsorted','nonfinite','negative','duplicate_field','incomplete'):
    assert not (root/directory/'manifest.json').exists()
print('native_export_sorted_ids_finite_sha_truthful_gates_partial_failure passed')
PY
for variant in "$baseline" "$candidate"; do
    stage="heartbeat-compile-$variant"
    if test -n "$previous" && test "$variant" = P1ABefore; then
        cp "$previous/binaries/P1ABefore-heartbeat" "$work/$variant/heartbeat"
    else
        make -f "$common/benchmarks/fsi-time-opt/heartbeat_case.mk" SOURCE="$work/$variant" WORKLOAD_SOURCE="$common" BUILD="$work/$variant" heartbeat > "$work/evidence/$variant-build.stdout" 2> "$work/evidence/$variant-build.stderr"
    fi
    cp "$work/$variant/heartbeat" "$work/evidence/binaries/$variant-heartbeat"
    sha256sum "$work/$variant/heartbeat" >> "$work/evidence/binary.sha256"
    stage="heartbeat-usage-$variant"
    if mpiexec -np 1 --map-by slot:PE=1 --bind-to core "$work/$variant/heartbeat" > "$work/evidence/$variant-usage.stdout" 2> "$work/evidence/$variant-usage.stderr"; then
        echo "missing required driver arguments unexpectedly succeeded" >&2; exit 1
    else
        result=$?; test "$result" -eq 1
        grep -q 'usage: bifurcation_fsi' "$work/evidence/$variant-usage.stderr"
    fi
done
stage=complete
