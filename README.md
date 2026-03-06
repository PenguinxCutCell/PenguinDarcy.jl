# PenguinDarcy.jl

`PenguinDarcy.jl` solves steady and transient Darcy pressure equations on Cartesian cut-cell grids and reconstructs Darcy velocity from the discrete gradient operator.

It is intentionally a **scalar pressure** package (no saddle-point `(u,p)` backend), consistent with the Penguin ecosystem split.

## Scope

Current implemented scope:

- Monophasic Darcy (steady + unsteady)
- Explicit two-domain fixed-interface Darcy (steady + unsteady)
- Fixed geometry only
- Outer box BCs (`Dirichlet`, `Neumann`, `Robin`, `Periodic`)
- Gravity/body-force forcing
- Distributed sources + conservative wells
- Scalar/diagonal/full-tensor mobility
- Velocity/flux reconstruction as postprocessing

Not implemented:

- moving geometry
- mixed `(u,p)` saddle-point formulations

## Governing equations

Pressure form:

- `u = -Λ(∇p - ρg)`
- `S ∂t p + ∇·u = s + q_w`

Steady is recovered with `S = 0`.

Two-domain fixed-interface form uses one pressure equation per domain plus interface constraints.

## Variable mode: pressure vs head

`variable = :pressure` (default):
- solve in pressure form with explicit `ρg` forcing.

`variable = :head`:
- use the same internal scalar path with head-like form (`ρg` forcing omitted in flux map).

## Tensor mobility

Supported mobility inputs:

- scalar (`λ::Number` or callable)
- diagonal (`NTuple{N}` / vector of component mobilities)
- full tensor (`N×N` constant matrix or callable returning `N×N`)

Full tensors use cross-face interpolation (`cross_face_interp`) from `CartesianOperators.jl` so off-diagonal coupling is active.

## Wells and sources

- `source` contributes volumetrically via cell wet volume.
- Wells are conservative by construction.

Available well types:

- `PointWell(x0, rate; radius=..., phase=...)`
- `CellWell(mask_or_indices, rate; phase=...)`

## Interface models (two-domain)

- `DarcyContinuity(pressure_jump, flux_jump)` (default `DarcyContinuity(0,0)`)
- `DarcyMembrane(σ, flux_jump)`

Unknown layout follows Penguin conventions:

- mono: `[pω, pγ]`
- two-domain: `[p1ω, p1γ, p2ω, p2γ]`

## Quick mono steady example

```julia
using CartesianGeometry: geometric_moments, nan
using CartesianOperators
using PenguinBCs
using PenguinDarcy

grid = (range(0.0, 1.0; length=65),)
moms = geometric_moments((x, t=0.0) -> -1.0, grid, Float64, nan; method=:vofijul)
cap = assembled_capacity(moms; bc=0.0)
ops = DiffusionOps(cap; periodic=periodic_flags(BorderConditions(), 1))

bc = BorderConditions(; left=Dirichlet(1.0), right=Dirichlet(0.0))
model = DarcyModelMono(cap, ops, 2.0;
    source=(x,t)->0.0,
    storage=1.0,
    ρ=1.0,
    gravity=(0.0,),
    bc_border=bc,
)

sys = solve_steady!(model)
vel = recover_velocity(model, sys.x)
mb = mass_balance(model, sys.x)
```

## Quick two-domain example

```julia
model = DarcyModelDiph(
    cap1, ops1, λ1,
    cap2, ops2, λ2;
    source=((x...)->0.0, (x...)->0.0),
    storage=(S1, S2),
    ρ=(ρ1, ρ2),
    gravity=(0.0, -9.81),
    bc_border=bc,
    bc_interface=DarcyContinuity(0.0, 0.0),
)

sys = solve_steady!(model)
flux = recover_flux(model, sys.x)
mb = mass_balance(model, sys.x)
```

## Main API

- `DarcyModelMono`, `DarcyModelDiph`
- `DarcyContinuity`, `DarcyMembrane`
- `PointWell`, `CellWell`
- `assemble_steady_mono!`, `assemble_unsteady_mono!`
- `assemble_steady_diph!`, `assemble_unsteady_diph!`
- `solve_steady!`, `solve_unsteady!`
- `recover_velocity`, `recover_flux`
- `mass_balance`, `compute_mass_balance`
- `boundary_discharge`, `interface_discharge`
- `integrated_source`, `integrated_well_rate`
- `face_mobility_values`

## Examples

- `examples/steady_1d_linear_pressure.jl`
- `examples/steady_2d_linear_gradient.jl`
- `examples/unsteady_1d_mms.jl`
- `examples/embedded_impermeable_circle.jl`
- `examples/gravity_hydrostatic_1d.jl`
- `examples/tensor_rotated_2d.jl`
- `examples/well_pair_2d.jl`
- `examples/two_domain_1d_piecewise.jl`
- `examples/two_domain_circle_2d.jl`

Each example prints pressure extrema, integrated source/well rates, boundary discharges, interface discharge (when applicable), and global mass-balance residual.

## Roadmap

- `v0.1`: monophasic steady/unsteady fixed geometry
- `v0.2`: gravity + wells + tensor mobility
- `v0.3`: explicit two-domain fixed-interface Darcy
- later: moving geometry
- later: mixed `(u,p)` backend for Brinkman/Stokes-Darcy coupling
