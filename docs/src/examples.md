**Examples and Tests**

Run examples from repository root:

```bash
julia --project=. examples/steady_1d_linear_pressure.jl
```

Available examples

- `examples/steady_1d_linear_pressure.jl`
  - 1D steady Darcy pressure benchmark.

- `examples/unsteady_1d_mms.jl`
  - 1D unsteady manufactured-solution style validation.

- `examples/steady_2d_linear_gradient.jl`
  - 2D steady linear-gradient pressure case.

- `examples/embedded_impermeable_circle.jl`
  - Embedded impermeable boundary using interface Robin condition.

Run tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Test files

- `test/test_steady_1d.jl`
- `test/test_unsteady_1d.jl`
- `test/test_steady_2d.jl`
- `test/test_embedded_boundary.jl`
