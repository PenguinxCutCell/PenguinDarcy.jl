**API and Types**

Primary type

- `DarcyModelMono{N,T,LT,ST,ΣT,IT}`
  - Fields: `ops`, `cap`, `λ`, `source`, `storage`, `bc_border`, `bc_interface`, `layout`, `coeff_mode`
  - Constructor:
    - `DarcyModelMono(cap, ops, λ; source=..., storage=..., bc_border=..., bc_interface=..., layout=..., coeff_mode=...)`

Assembly

- `assemble_steady_mono!(sys, model, t)`
- `assemble_unsteady_mono!(sys, model, uⁿ, t, dt, scheme)`

Solvers

- `solve_steady!(model; t=0, method=:direct, kwargs...)`
- `solve_unsteady!(model, u0, tspan; dt, scheme=:BE|:CN|θ, method=:direct, save_history=true, kwargs...)`

Post-processing

- `recover_flux(model, state; t=0)`
- `recover_velocity(model, state; t=0)` (alias of `recover_flux`)
- `face_mobility_values(model; t=0)`
- `face_mobility_values(cap, λ, t, mode)`

Diagnostics

- `boundary_discharge(model, state, side; t=0)`
- `compute_mass_balance(model, state; t=0)`

Notes

- `u0` may be either the `ω` block (`length = ntotal`) or the full `ω+γ` state.
- Supported unsteady schemes are `:BE`, `:CN`, or numeric `θ`.
