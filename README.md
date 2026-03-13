# PenguinDarcy.jl

[![In development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://PenguinxCutCell.github.io/PenguinDarcy.jl/dev)
![CI](https://github.com/PenguinxCutCell/PenguinDarcy.jl/actions/workflows/ci.yml/badge.svg)
![Coverage](https://codecov.io/gh/PenguinxCutCell/PenguinDarcy.jl/branch/main/graph/badge.svg)

`PenguinDarcy.jl` solves Darcy pressure equations on Cartesian cut-cell grids.

The package stays pressure-only (no mixed `(u,p)` saddle-point solve), and now includes both fixed-geometry and graph-based moving/free-boundary Darcy models.

## Feature matrix

| Capability | Status |
|---|---|
| Mono fixed geometry, steady | implemented |
| Mono fixed geometry, unsteady (`storage != 0`) | implemented |
| Two-phase fixed interface, steady | implemented |
| Two-phase fixed interface, unsteady (`storage != 0`) | implemented |
| Mono free-boundary (graph interface), quasi-steady | implemented |
| Two-phase free-boundary (graph interface), quasi-steady | implemented |
| Free-boundary backend: `HeightFunctionTracker` (`GlobalHeightFunctions.jl`) | implemented |
| Free-boundary backend: level-set | not implemented |
| Free-boundary backend: gVOF transport | not implemented |
| Free-boundary with topology changes / bubbles / droplets | not implemented |
| Moving two-phase with nonzero interfacial `flux_jump` | not implemented |
| Moving-geometry transient storage (`storage != 0`) | not implemented |

Important limitation of this first moving release: **graph interfaces only**.

## Governing models

Darcy velocity:

- `u = -Λ(∇p - ρg)`

Mono bulk equation:

- `∇·u = s + q_w` (quasi-steady moving mode)
- `S ∂t p + ∇·u = s + q_w` (fixed-geometry unsteady mode)

Two-phase interface conditions (moving mode):

- `[u·n] = 0`
- `V_n = (u₁·n)/ϕ_Γ = (u₂·n)/ϕ_Γ`
- `[p] = J_p` with optional capillary contribution `σ κ`

Graph kinematics uses the correct axis-aware form. For `x = h(y,t)`:

- `h_t = u_x - u_y h_y = V_n / (n·e_x)`

## Scope notes for moving models

- `HeightFunctionTracker` is the current tracker backend.
- `interface_porosity` in kinematics defaults to `1`, but is configurable.
- Surface tension is optional (`surface_tension = 0` by default).
- If moving models are built with nonzero `storage`, constructors throw a clear error.
- If moving two-phase models use nonzero `flux_jump`, constructors throw a clear error.

## Quick moving example (mono)

```julia
using PenguinDarcy
using CartesianGrids
using PenguinBCs

grid = CartesianGrid((0.0, 0.0), (1.0, 1.0), (81, 65))
y = collect(grid1d(grid, 2))
xf0 = fill(0.43, length(y)); xf0[end] = xf0[1]

tracker = HeightFunctionTracker(grid, xf0; axis=:x, periodic_transverse=true)
bc = BorderConditions(; left=Neumann(-0.2), right=Neumann(0.0), bottom=Periodic(), top=Periodic())

model = MovingDarcyModelMono(tracker, 1.0;
    source=(x,y,t)->0.0,
    bc_border=bc,
    p_ext=0.0,
)

sol = solve_unsteady_moving!(model, (0.0, 0.05); dt=0.01)
```

## Main API

Fixed geometry:

- `DarcyModelMono`, `DarcyModelDiph`
- `DarcyContinuity`, `DarcyMembrane`
- `assemble_steady_mono!`, `assemble_unsteady_mono!`
- `assemble_steady_diph!`, `assemble_unsteady_diph!`
- `solve_steady!`, `solve_unsteady!`

Moving/free-boundary:

- `AbstractDarcyTracker`, `HeightFunctionTracker`
- `MovingDarcyModelMono`, `MovingDarcyModelDiph`
- `assemble_unsteady_moving!`, `solve_unsteady_moving!`
- `interface_normal_velocity`, `update_interface!`
- `interface_mass_residual`, `interface_iteration_history`

Postprocessing/diagnostics:

- `recover_velocity`, `recover_flux`
- `mass_balance`, `compute_mass_balance`
- `boundary_discharge`, `interface_discharge`
- `integrated_source`, `integrated_well_rate`
- `face_mobility_values`

## Examples

Fixed-geometry examples:

- `examples/steady_1d_linear_pressure.jl`
- `examples/steady_2d_linear_gradient.jl`
- `examples/unsteady_1d_mms.jl`
- `examples/embedded_impermeable_circle.jl`
- `examples/gravity_hydrostatic_1d.jl`
- `examples/tensor_rotated_2d.jl`
- `examples/well_pair_2d.jl`
- `examples/two_domain_1d_piecewise.jl`
- `examples/two_domain_circle_2d.jl`

Moving/free-boundary examples:

- `examples/free_boundary_mono_planar_translation_2d.jl`
- `examples/free_boundary_mono_hydrostatic_equilibrium_2d.jl`
- `examples/free_boundary_diph_planar_translation_2d.jl`
- `examples/free_boundary_diph_muskat_perturbation_2d.jl`

## Roadmap

- Add level-set and gVOF tracker backends through the same `AbstractDarcyTracker` path.
- Add moving-geometry support for nonzero interfacial mass transfer and richer topology handling.
