#!/bin/bash
# Supply the allocation at submission: sbatch -A YOUR_PSC_PROJECT this-file.sbatch
#SBATCH -J iga-graph
#SBATCH -p RM
#SBATCH -N 2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH -t 02:00:00
#SBATCH --signal=B:USR1@300
#SBATCH --requeue
#SBATCH -o slurm/graph-%j.out

set -euo pipefail
# sbatch executes a spool copy; resolve the checkout from the submission directory.
repo_root=${IGA_REPO_ROOT:-${SLURM_SUBMIT_DIR:-}}
if [[ -z ${repo_root} ]]; then
	repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd)
fi
repo_root=$(cd -- "${repo_root}" && pwd)
if [[ ! -f ${repo_root}/scripts/hpc_scheduler_record.py ]]; then
	echo "submit from the repository root or set IGA_REPO_ROOT" >&2
	exit 2
fi
cd "${repo_root}"
module load anaconda3
source_case=${IGA_GRAPH_CASE:?Set IGA_GRAPH_CASE to a graph case on the shared filesystem}
shared_output=${IGA_OUTPUT_ROOT:?Set IGA_OUTPUT_ROOT to a new shared output root}
checkpoint_root=${IGA_CHECKPOINT_ROOT:?Set IGA_CHECKPOINT_ROOT on the shared filesystem}
: "${PETSC_DIR:?Set PETSC_DIR to the PETSc installation used for this build}"
PETSC_ARCH=${PETSC_ARCH:-}
module load openmpi/4.0.5-gcc10.2.0

attempt=${SLURM_RESTART_COUNT:-0}

attempt_root=${shared_output}/attempt-${attempt}
result=${attempt_root}/scheduler.json
command=("$0")
mkdir -p "${attempt_root}" "${checkpoint_root}"
record_failure()
{
	local status=$?
	trap - EXIT
	if (( status != 0 )); then
		python3 "${repo_root}/scripts/hpc_scheduler_record.py" \
			--output "${result}" --status failed --returncode "${status}" -- "${command[@]}" || true
	fi
	exit "${status}"
}
trap record_failure EXIT

if [[ ${SLURM_JOB_NUM_NODES:-0} -lt 2 || ${SLURM_NTASKS:-0} -lt 2 ]]; then
	echo "multinode_graph.sbatch requires at least two nodes and two tasks" >&2
	exit 2
fi
if [[ ! -x ${repo_root}/solvers/coupling/iga_multidomain_flow ]]; then
	echo "build solvers/coupling/iga_multidomain_flow before submission" >&2
	exit 2
fi
if [[ ${IGA_GRAPH_PREPARE:-0} == 1 && ${attempt} == 0 ]]; then
	fixture=${repo_root}/solvers/coupling/hpc_duct_solver_fixture
	if [[ ! -x ${fixture} ]]; then
		echo "build solvers/coupling/hpc_duct_solver_fixture before submission" >&2
		exit 2
	fi
	"${fixture}" "${source_case}" "${IGA_GRAPH_TRANSVERSE:-16}" \
		"${IGA_GRAPH_AXIAL:-64}" "${SLURM_NTASKS}" "${IGA_GRAPH_STEPS:-8}"
fi
if [[ ${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK}} != "${SLURM_CPUS_PER_TASK}" ]]; then
	echo "OMP_NUM_THREADS must equal SLURM_CPUS_PER_TASK" >&2
	exit 2
fi

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_THREAD_LIMIT=${SLURM_CPUS_PER_TASK}
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export IGA_PROFILE=1

test -d "${source_case}"
srun --nodes="${SLURM_JOB_NUM_NODES}" --ntasks="${SLURM_JOB_NUM_NODES}" \
	--ntasks-per-node=1 --cpus-per-task=1 hostname | sort -u >"${attempt_root}/nodes.txt"
if [[ $(wc -l <"${attempt_root}/nodes.txt") -lt 2 ]]; then
	echo "graph allocation did not expose two distinct hosts" >&2
	exit 2
fi

# Each node receives an identical private input tree. Outputs and checkpoints
# stay on the shared filesystem so a later allocation can resume them.
# PSC exports LOCAL in job steps; SLURM_TMPDIR is supported at other sites.
export SLURM_TMPDIR=${SLURM_TMPDIR:-${LOCAL:?node-local storage is unavailable}}
export IGA_STAGE_SOURCE=$(realpath "${source_case}")
srun --nodes="${SLURM_JOB_NUM_NODES}" --ntasks="${SLURM_JOB_NUM_NODES}" \
	--ntasks-per-node=1 --cpus-per-task=1 \
	bash -c 'mkdir -p "$SLURM_TMPDIR/tubularflowiga-graph"; rsync -a --delete "$IGA_STAGE_SOURCE/" "$SLURM_TMPDIR/tubularflowiga-graph/"'
staged_case=${SLURM_TMPDIR}/tubularflowiga-graph

python3 "${repo_root}/scripts/hpc_build_manifest.py" --component cpu \
	--dependency-status passed --output "${attempt_root}/build.json" \
	--binary "${repo_root}/solvers/coupling/iga_multidomain_flow"

restart=()
if find "${checkpoint_root}" -name manifest -type f -print -quit | grep -q .; then
	restart=(--restart-dir "${checkpoint_root}")
fi
command=(mpiexec --mca ess_base_forward_signals SIGUSR1 --map-by \
	"ppr:${SLURM_NTASKS_PER_NODE}:node:PE=${SLURM_CPUS_PER_TASK}" --bind-to core \
	--report-bindings -np "${SLURM_NTASKS}" \
	"${repo_root}/solvers/coupling/iga_multidomain_flow" \
	--graph-case "${staged_case}" --output-dir "${attempt_root}/results" \
	--checkpoint-dir "${checkpoint_root}" --checkpoint-every "${IGA_CHECKPOINT_EVERY:-10}" \
	"${restart[@]}")
python3 "${repo_root}/scripts/hpc_scheduler_record.py" --output "${result}" \
	--status starting -- "${command[@]}"

checkpoint_requested=0
launcher_pid=
timer_pid=
forward_checkpoint()
{
	checkpoint_requested=1
	if [[ -n ${launcher_pid} ]]; then kill -USR1 "${launcher_pid}"; fi
}
trap forward_checkpoint USR1
if [[ ${IGA_CHECKPOINT_TEST_DELAY:-0} != 0 && ${attempt} == 0 ]]; then
	batch_pid=$$
	(
		sleep "${IGA_CHECKPOINT_TEST_DELAY}"
		kill -USR1 "${batch_pid}"
	) &
	timer_pid=$!
fi
set +e
"${command[@]}" >"${attempt_root}/solver.stdout" 2>"${attempt_root}/solver.stderr" &
launcher_pid=$!
wait "${launcher_pid}"
status=$?
while (( checkpoint_requested && status >= 128 )) && kill -0 "${launcher_pid}" 2>/dev/null; do
	wait "${launcher_pid}"
	status=$?
done
set -e
if [[ -n ${timer_pid} ]]; then
	kill "${timer_pid}" 2>/dev/null || true
	wait "${timer_pid}" 2>/dev/null || true
fi
trap - USR1

record_args=(--output "${result}" --returncode "${status}")
if (( checkpoint_requested )); then record_args+=(--checkpoint-requested); fi
if (( status != 0 )); then
	python3 "${repo_root}/scripts/hpc_scheduler_record.py" "${record_args[@]}" \
		--status failed -- "${command[@]}"
	exit "${status}"
fi
if (( checkpoint_requested )); then
	python3 "${repo_root}/scripts/hpc_scheduler_record.py" "${record_args[@]}" \
		--status checkpointed -- "${command[@]}"
	if [[ ${IGA_REQUEUE_ON_SIGNAL:-1} == 1 ]]; then
		# Launcher startup may reset the flag; enable it after the accepted checkpoint.
		scontrol update "JobId=${SLURM_JOB_ID}" Requeue=1
		scontrol requeue "${SLURM_JOB_ID}"
		exit 0
	fi
else
	python3 "${repo_root}/scripts/hpc_scheduler_record.py" "${record_args[@]}" \
		--status passed -- "${command[@]}"
fi
