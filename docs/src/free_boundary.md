**Free-Boundary Darcy / Hele-Shaw Models**

This page documents the first moving-interface release in `PenguinDarcy.jl`.

Implemented scope:

- one-phase free-boundary Darcy (graph interface)
- two-phase free-boundary Darcy / Muskat style (graph interface)
- quasi-steady pressure solve + implicit/damped interface update
- tracker backend: `HeightFunctionTracker` (powered by `GlobalHeightFunctions.jl`)

Not implemented in this release:

- topology changes (droplets/bubbles/multi-valued interface)
- level-set or gVOF tracker backends
- moving models with nonzero interfacial mass transfer (`flux_jump != 0`)
- moving models with nonzero pressure storage (`storage != 0`)

## Equations

Bulk Darcy law:

- `u = -Λ(∇p - ρg)`

One-phase moving mode (quasi-steady):

- `∇·u = s + q_w`
- free boundary: `p = p_ext + σ κ` (with `σ=0` by default)
- kinematics: `V_n = (u·n)/ϕ_Γ`

Two-phase moving mode (quasi-steady):

- `∇·u_i = s_i + q_{w,i}` in each phase
- interface conditions:
  - `[u·n] = 0`
  - `V_n = (u₁·n)/ϕ_Γ = (u₂·n)/ϕ_Γ`
  - `[p] = J_p + σ κ`

`ϕ_Γ` is controlled by `interface_porosity` (default `1`).

## Sign conventions and normals

The moving two-phase path reuses the same interface sign conventions as fixed-interface `DarcyContinuity`:

- pressure jump is `p1 - p2`
- moving mode currently requires `flux_jump = 0`

For graph interfaces, phase 1 is the `phi < 0` side of

- `phi(x) = x_axis - h(transverse)`

and the interface normal is `n = ∇phi / |∇phi|`.

## Graph kinematics

For `x = h(y,t)`:

- `h_t = u_x - u_y h_y`
- equivalent form used in code: `h_t = V_n / (n·e_x)`

The same axis-aware formula is used for `axis=:y` and `axis=:z`.

## Capillary jump

Curvature is computed directly from the graph representation.

- in 2D: graph curvature from first/second derivatives of `h`
- in 3D: graph mean-curvature expression from `h_y, h_z, h_yy, h_zz, h_yz`

A lightweight optional smoothing pass is available in `HeightFunctionTracker(curvature_smoothing=...)`.

## Nonlinear step algorithm

At each time step:

1. Predict interface state from tracker history.
2. Rebuild geometry/capacities from the predicted graph.
3. Assemble/solve Darcy pressure in the rebuilt geometry.
4. Reconstruct interface normal velocity.
5. Solve implicit graph update residual
   `R(ξ) = ξ - ξ^n - Δt G(ξ)`
   with damped Picard iterations.
6. Accept converged `ξ^{n+1}`, store iteration/mass diagnostics.

Convergence checks:

- `||R||∞ < tol_interface`
- `||δξ||∞ < tol_update`

## Public moving API

- `AbstractDarcyTracker`, `HeightFunctionTracker`
- `MovingDarcyModelMono`, `MovingDarcyModelDiph`
- `assemble_unsteady_moving!`, `solve_unsteady_moving!`
- `interface_normal_velocity`, `update_interface!`
- `interface_mass_residual`, `interface_iteration_history`

## Examples

- `examples/free_boundary_mono_planar_translation_2d.jl`
- `examples/free_boundary_mono_hydrostatic_equilibrium_2d.jl`
- `examples/free_boundary_diph_planar_translation_2d.jl`
- `examples/free_boundary_diph_muskat_perturbation_2d.jl`
