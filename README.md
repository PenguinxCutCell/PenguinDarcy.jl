# PenguinDarcy.jl

`PenguinDarcy.jl` solves steady and transient Darcy pressure equations on Cartesian cut-cell grids and reconstructs Darcy velocity from the discrete gradient operator.

## Scope (v0.1)

- Monophasic Darcy only
- Steady and unsteady pressure solves
- Fixed geometry only (no moving body support)
- Outer box BCs (`Dirichlet`, `Neumann`, `Robin`, `Periodic`)
- Optional fixed embedded mono interface BC through `Robin`
- Postprocessed Darcy velocity/flux recovery

Not included in `v0.1`:

- diphasic/two-domain Darcy
- mixed `(u, p)` saddle-point solves
- moving geometry

## Governing equations

Steady:

`u = -λ∇p`, `∇·u = s`  ->  `-∇·(λ∇p) = s`

Unsteady:

`S ∂t p - ∇·(λ∇p) = s`

where:

- `p`: pressure
- `u`: Darcy velocity
- `λ`: mobility (`κ/μ`)
- `S`: storage coefficient
- `s`: source term

## BC conventions

- `Dirichlet`: imposed pressure `p`
- `Neumann`: imposed Darcy normal flux `u·n`
- `Robin`: `α p + β (u·n) = g`

For embedded impermeable boundaries, use `Robin(0, 1, 0)`.

## Quick steady example

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
model = DarcyModelMono(cap, ops, 2.0; source=(x, t) -> 0.0, bc_border=bc)

sys = solve_steady!(model)
vel = recover_velocity(model, sys.x)
mb = compute_mass_balance(model, sys.x)
```

## Quick unsteady example

```julia
model = DarcyModelMono(cap, ops, 1.0;
    source=(x, t) -> 0.0,
    storage=1.0,
    bc_border=bc,
)

u0 = zeros(cap.ntotal)
sol = solve_unsteady!(model, u0, (0.0, 0.5); dt=0.01, scheme=:CN)
```

## Public API

- `DarcyModelMono`
- `assemble_steady_mono!`
- `assemble_unsteady_mono!`
- `solve_steady!`
- `solve_unsteady!`
- `recover_velocity`
- `recover_flux`
- `compute_mass_balance`
- `boundary_discharge`
- `face_mobility_values`

## Roadmap

- `v0.1`: monophasic steady/unsteady fixed geometry
- later: gravity/head wrapper
- later: moving geometry
- later: explicit two-domain Darcy interface solver
- later: mixed `(u,p)` backend if needed for Brinkman/Stokes-Darcy coupling
