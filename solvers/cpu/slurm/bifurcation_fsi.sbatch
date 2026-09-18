#!/bin/bash
#SBATCH -J iga-bifurcation-fsi
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=15200M
#SBATCH -t 24:00:00
set -euo pipefail
module load anaconda3 openmpi/4.0.5-gcc10.2.0
repo=${IGA_REPO_ROOT:-${SLURM_SUBMIT_DIR:?}}
root=${IGA_FSI_SHOWCASE_ROOT:?Set a shared showcase workspace}
: "${PETSC_DIR:?Set the compatible PETSc installation}"
export CXX=mpicxx
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_THREAD_LIMIT=$OMP_NUM_THREADS OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores
export IGA_ASSEMBLY_THREADS=$OMP_NUM_THREADS IGA_PROFILE=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
attempt="$root/job-$SLURM_JOB_ID"
mkdir -p "$root"
mkdir "$attempt"
record_exit() {
    code=$?
    trap - EXIT
    sacct --json -j "$SLURM_JOB_ID" > "$attempt/scheduler.json" 2> "$attempt/scheduler.stderr" || true
    exit "$code"
}
trap record_exit EXIT
cd "$repo"
git rev-parse HEAD > "$attempt/source-commit.txt"
git status --porcelain=v1 > "$attempt/source-status.txt"
git diff > "$attempt/source.patch"
cp solvers/cpu/src/bifurcation_fsi.cpp "$attempt/source.cpp"
for header in ImmersedTransientFlowRuntime MovingImmersedTransientFlowRuntime MovingImmersedTransientFlowFsiRuntime StrongFluidStructureCoupling ImmersedNitscheWall; do
    cp "solvers/cpu/include/$header.hpp" "$attempt/$header.hpp"
done
cp examples/vascular_flow/bifurcation_fsi/generate_surface.py "$attempt/generate_surface.py"
cp examples/vascular_flow/bifurcation_fsi/render.py "$attempt/render.py"
cp examples/vascular_flow/bifurcation_fsi/verify.py "$attempt/verify.py"
cp examples/vascular_flow/bifurcation_fsi/verify_vtk.py "$attempt/verify_vtk.py"
hostname > "$attempt/hostname.txt"
lscpu > "$attempt/hardware.txt"
module list 2> "$attempt/modules.txt"
make -C solvers/cpu bifurcation_fsi > "$attempt/build.stdout" 2> "$attempt/build.stderr"
cp solvers/cpu/bifurcation_fsi "$attempt/bifurcation_fsi"
ldd "$attempt/bifurcation_fsi" > "$attempt/linkage.txt"
python3 "$attempt/generate_surface.py" "$attempt/surface.txt" > "$attempt/geometry.stdout" 2> "$attempt/geometry.stderr"
sha256sum "$attempt/bifurcation_fsi" "$attempt/source.cpp" "$attempt/surface.txt" > "$attempt/hashes.txt"
/usr/bin/time -v mpiexec -np 1 --map-by "slot:PE=$OMP_NUM_THREADS" --bind-to core \
    "$attempt/bifurcation_fsi" "$attempt/surface.txt" "$attempt/results" \
    -demo_steps "${IGA_FSI_STEPS:-8}" -demo_grid "${IGA_FSI_GRID:-16}" \
    -demo_dt "${IGA_FSI_DT:-.01}" -demo_pressure "${IGA_FSI_PRESSURE:-20}" \
    -demo_wall_inertial_gamma "${IGA_FSI_WALL_INERTIAL_GAMMA:-1}" \
    -domain_fluid_flow_ksp_converged_reason \
    > "$attempt/solver.stdout" 2> "$attempt/solver.stderr"
python3 "$attempt/render.py" "$attempt/results" "$attempt/visuals" \
    > "$attempt/render.stdout" 2> "$attempt/render.stderr"
