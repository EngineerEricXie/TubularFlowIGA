#!/bin/bash
#SBATCH -J iga-heartbeat
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=15200M
#SBATCH -t 04:00:00
set -euo pipefail
module load anaconda3 openmpi/4.0.5-gcc10.2.0
repo=${IGA_REPO_ROOT:-${SLURM_SUBMIT_DIR:?}}
: "${PETSC_DIR:?}" "${PETSC_ARCH:?}"
export CXX=mpicxx OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_THREAD_LIMIT=$OMP_NUM_THREADS OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores
export IGA_ASSEMBLY_THREADS=$OMP_NUM_THREADS IGA_PROFILE=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
root="$repo/.nsvms_diagnostics/heartbeat/job-$SLURM_JOB_ID"
mkdir -p "$root/source/solvers/cpu/src"
cd "$repo"
git diff > "$root/source.patch"
git rev-parse HEAD > "$root/source-commit.txt"
cp -a include "$root/source/"
cp -a solvers/cpu/include solvers/cpu/Makefile "$root/source/solvers/cpu/"
cp solvers/cpu/src/bifurcation_fsi.cpp "$root/source/solvers/cpu/src/"
cp -a examples/vascular_flow/bifurcation_fsi "$root/scripts"
cp solvers/cpu/slurm/bifurcation_heartbeat.sbatch "$root/submission.sbatch"
module list 2> "$root/modules.txt"
hostname > "$root/hostname.txt"
make -C "$root/source/solvers/cpu" bifurcation_fsi > "$root/build.stdout" 2> "$root/build.stderr"
python3 "$root/scripts/generate_surface.py" "$root/surface.txt" > "$root/geometry.stdout"
sha256sum "$root/source/solvers/cpu/bifurcation_fsi" "$root/surface.txt" > "$root/hashes.txt"
mpiexec -np 1 --map-by "slot:PE=$OMP_NUM_THREADS" --bind-to core \
 "$root/source/solvers/cpu/bifurcation_fsi" "$root/surface.txt" "$root/results" \
 -demo_visualization_only true -demo_steps 20 -demo_dt .05 -demo_period 1 \
 -demo_pressure 5 -demo_pulse_amplitude 30 -demo_grid 16 -demo_cut_depth 1 \
 -demo_display_scale 100 -demo_wall_inertial_gamma 1 \
 > "$root/solver.stdout" 2> "$root/solver.stderr"
python3 "$root/scripts/verify_heartbeat.py" "$root/results" > "$root/verification.stdout" 2> "$root/verification.stderr"
python3 "$root/scripts/render_heartbeat.py" "$root/results" > "$root/render.stdout" 2> "$root/render.stderr"
