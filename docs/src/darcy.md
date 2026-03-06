**Darcy Model and PDE**

`PenguinDarcy.jl` solves steady and unsteady single-phase Darcy pressure equations on Cartesian cut-cell grids.

Steady:

`-∇·(λ∇p) = s`

Unsteady:

`S ∂t p - ∇·(λ∇p) = s`

where:

- `p`: pressure
- `λ`: mobility (`κ/μ`)
- `S`: storage
- `s`: source term

Boundary conventions

- `Dirichlet`: imposed pressure `p`
- `Neumann`: imposed normal Darcy flux `u·n`
- `Robin`: `α p + β (u·n) = g`

Embedded interface condition (optional)

- `bc_interface` supports `PenguinBCs.Robin`
- If `bc_interface = nothing`, the interface block is constrained with identity rows

Coefficient sampling

- `coeff_mode = :harmonic` (default), `:arithmetic`, `:face`, or `:cell`
- Face mobilities are produced by `face_mobility_values`

Recover Darcy velocity : `u = -λ∇p` with `recover_flux` or `recover_velocity` (alias)