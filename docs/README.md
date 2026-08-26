# Documentation index

Start with the root [README](../README.md) for capabilities and the shortest
vascular and neuron workflows. Then choose the document matching the task:

| Task | Document |
|---|---|
| Install Linux, WSL, PETSc, MPI, METIS, or CUDA dependencies | [DEPENDENCIES.md](DEPENDENCIES.md) |
| Follow a fresh-clone vascular and neuron walkthrough | [QUICKSTART.md](QUICKSTART.md) |
| Understand generated files and stage interfaces | [PIPELINE.md](PIPELINE.md) |
| Configure fields, operators, time stepping, checkpointing, and solver CLI | [PDE_CONFIGURATION.md](PDE_CONFIGURATION.md) |
| Configure native CPU 3D VCA coupling and review its validation scope | [THREE_D_VCA_IMPLEMENTATION_PLAN.md](THREE_D_VCA_IMPLEMENTATION_PLAN.md) |
| Configure wall, inlet, outlet, waveform, and Windkessel conditions | [BOUNDARY_CONDITIONS.md](BOUNDARY_CONDITIONS.md) |
| Run on PSC Bridges-2 | [BRIDGES2.md](BRIDGES2.md) |
| Review mesh correctness evidence | [MESH_CPP_VALIDATION.md](MESH_CPP_VALIDATION.md) |
| Review performance measurements | [BENCHMARKS.md](BENCHMARKS.md) |

Implementation-specific guides live next to their code:

- [Control-mesh generator](../preprocessing/mesh/README.md)
- [Spline and Bezier extraction](../preprocessing/spline/README.md)
- [MPI/PETSc CPU backend](../solvers/cpu/README.md)
- [CUDA backend](../solvers/cuda/README.md)

Runnable cases and their three-file input contract are indexed in
[examples/README.md](../examples/README.md).
