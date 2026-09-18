#!/bin/bash
#SBATCH -J heartbeat-opt
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=15200M
#SBATCH -t 12:00:00
set -euo pipefail
module load anaconda3 openmpi/4.0.5-gcc10.2.0
export REPO=/ocean/projects/mch260002p/thsieh1/TubularFlowIGA
export PETSC_DIR=/ocean/projects/mch260002p/thsieh1/TubularFlowIGA-hpc-acceptance/petsc PETSC_ARCH=arch-bridges2-acceptance CXX=mpicxx
export OMP_NUM_THREADS=8 OMP_THREAD_LIMIT=8 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores
export IGA_ASSEMBLY_THREADS=8 IGA_PROFILE=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
previous="$REPO/.nsvms_diagnostics/heartbeat/job-45971412"
root="$REPO/benchmarks/heartbeat-optimized-$SLURM_JOB_ID"
work="${LOCAL:-/tmp}/heartbeat-optimized-$SLURM_JOB_ID"
mkdir -p "$root" "$work"
cp -a "$previous/source" "$work/"
cp -a "$previous/scripts" "$root/"
tar -xzf "$REPO/benchmarks/heartbeat-performance-45981715.tar.gz" -C "$work" ./candidate/include/ImmersedTransientFlowRuntime.hpp ./comparison.json
python3 - "$work/comparison.json" <<'PY'
import json,sys
assert json.load(open(sys.argv[1]))['status']=='passed'
PY
cp "$work/candidate/include/ImmersedTransientFlowRuntime.hpp" "$work/source/solvers/cpu/include/"
cp "$previous/surface.txt" "$root/surface.txt"
cp "$work/candidate/include/ImmersedTransientFlowRuntime.hpp" "$root/ImmersedTransientFlowRuntime.hpp"
printf '%s\n' "$previous/source" > "$root/base-source-snapshot.txt"
make -C "$work/source/solvers/cpu" bifurcation_fsi > "$root/build.stdout" 2> "$root/build.stderr"
cp "$work/source/solvers/cpu/bifurcation_fsi" "$root/bifurcation_fsi"
sha256sum "$root/bifurcation_fsi" "$root/ImmersedTransientFlowRuntime.hpp" "$root/surface.txt" > "$root/hashes.txt"
mpiexec -np 1 --map-by slot:PE=8 --bind-to core "$root/bifurcation_fsi" "$root/surface.txt" "$root/results" \
 -demo_visualization_only true -demo_steps 20 -demo_dt .05 -demo_period 1 \
 -demo_pressure 5 -demo_pulse_amplitude 30 -demo_grid 16 -demo_cut_depth 1 \
 -demo_display_scale 100 -demo_wall_inertial_gamma 1 \
 > "$root/solver.stdout" 2> "$root/solver.stderr"
python3 "$root/scripts/verify_heartbeat.py" "$root/results" > "$root/verification.stdout" 2> "$root/verification.stderr"
python3 "$root/scripts/render_heartbeat.py" "$root/results" > "$root/render.stdout" 2> "$root/render.stderr"
/jet/home/thsieh1/.conda/envs/research/bin/python "$root/scripts/verify_heartbeat_vtk.py" "$root/results" > "$root/vtk.stdout" 2> "$root/vtk.stderr"
tar -czf "$root/heartbeat-paraview.tar.gz" -C "$root" results
sha256sum "$root/heartbeat-paraview.tar.gz" > "$root/archive.sha256"
