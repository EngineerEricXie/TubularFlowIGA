# User case inputs

Place each private or local 1D or 3D case in its own immediate subdirectory:

```text
Input/
└── MyCase/
    ├── simulation_config.json
    └── skeleton_initial.swc
```

The geometry may instead be a radius-annotated line OBJ named by
`simulation_config.json`. Do not place generated mesh, database, or simulation
files in the source directory manually. Three-dimensional cases use the full
mesh/database layout; native 1D cases write their solver products below
`MyCase/generated/results/`.

Configure machine and launch settings in the repository-level
[`execution.conf`](../execution.conf), then run one or more case names:

```bash
./scripts/run_cases.sh MyCase
```

With no case names, the runner processes every immediate subdirectory here.
Committed examples remain under [`examples/`](../examples/README.md); this
directory is for user-provided source cases.
