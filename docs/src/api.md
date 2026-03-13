**API and Types**

## Fixed-geometry models

- `DarcyModelMono`
- `DarcyModelDiph`
- `DarcyContinuity`, `DarcyMembrane`

Assembly/solve:

- `assemble_steady_mono!`, `assemble_unsteady_mono!`
- `assemble_steady_diph!`, `assemble_unsteady_diph!`
- `solve_steady!`, `solve_unsteady!`

Postprocessing/diagnostics:

- `recover_flux`, `recover_velocity`
- `mass_balance`, `compute_mass_balance`
- `boundary_discharge`, `interface_discharge`
- `integrated_source`, `integrated_well_rate`
- `face_mobility_values`

## Moving/free-boundary models

- `AbstractDarcyTracker`
- `HeightFunctionTracker`
- `MovingDarcyModelMono`
- `MovingDarcyModelDiph`
- `assemble_unsteady_moving!`
- `solve_unsteady_moving!`
- `interface_normal_velocity`
- `update_interface!`
- `interface_mass_residual`
- `interface_iteration_history`

Notes:

- Moving models are quasi-steady in this release.
- `storage != 0` is rejected for moving geometry.
- Moving two-phase models require `flux_jump = 0`.

Internal API

```@docs
PenguinDarcy._apply_box_bc_darcy!
```

Moving docstrings

```@docs
PenguinDarcy.InterfaceIterationHistory
PenguinDarcy.MovingStepDiagnostics
PenguinDarcy.MovingDarcySolution
PenguinDarcy.tracker_dimension_type
PenguinDarcy.tracker_state
PenguinDarcy.interface_positions
PenguinDarcy.rebuild_signed_distance_or_geometry
PenguinDarcy.rebuild_capacities
PenguinDarcy.interface_normals
PenguinDarcy.interface_curvature
PenguinDarcy.update_interface!
PenguinDarcy.predictor_interface_state
PenguinDarcy.check_graph_validity
PenguinDarcy.interface_iteration_history
```
