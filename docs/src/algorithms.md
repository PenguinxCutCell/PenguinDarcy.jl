**Numerical Algorithms**

1) Steady assembly (`assemble_steady_mono!`)

- Build weighted core operators (`K`, `C`, `J`, `L`) from `DiffusionOps` and face mobility values.
- Assemble mono block system in `(ω, γ)` unknown layout.
- Apply source vector and interface Robin terms.
- Apply outer box boundary conditions with Darcy flux semantics.
- Constrain inactive/halo rows with identity rows.

2) Unsteady assembly (`assemble_unsteady_mono!`)

- Start from steady assembly at `t + θΔt`.
- For `θ != 1`, add explicit correction term for the `ω` block.
- Add storage mass diagonal `M = (S .* V) / Δt` to `ω` rows.
- Add previous-state contribution `M * pⁿ` to RHS.

3) Time integration (`solve_unsteady!`)

- Advance from `t0` to `tend` with fixed nominal `dt` and a clipped final step.
- Solve each step through `PenguinSolverCore.solve!`.
- Return `(times, states, system, reused_constant_operator=false)`.

4) Moving free-boundary integration (`solve_unsteady_moving!`)

- Predict a graph-interface state from tracker history.
- Rebuild cut-cell geometry from the graph (`HeightFunctionTracker`).
- Reuse Darcy steady assembly/solve on rebuilt geometry.
- Recover interface normal velocity and solve implicit graph residual with damped Picard iterations.
- Store per-step diagnostics: interface residual history, mass residual, CFL-like indicator, and velocity-jump metric (two-phase).

Performance notes

- Canonical `(ω, γ)` layouts are assembled as direct block matrices.
- Non-canonical layouts use sparse block insertion helpers.
- Face coefficients support harmonic/arithmetic/cell/face modes.
