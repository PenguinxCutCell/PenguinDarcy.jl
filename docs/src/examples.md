**Examples and Tests**

Run examples from repository root:

```bash
julia --project=. examples/steady_1d_linear_pressure.jl
```

Fixed-geometry examples

- `examples/steady_1d_linear_pressure.jl`
- `examples/steady_2d_linear_gradient.jl`
- `examples/unsteady_1d_mms.jl`
- `examples/embedded_impermeable_circle.jl`
- `examples/gravity_hydrostatic_1d.jl`
- `examples/tensor_rotated_2d.jl`
- `examples/well_pair_2d.jl`
- `examples/two_domain_1d_piecewise.jl`
- `examples/two_domain_circle_2d.jl`

Moving/free-boundary examples

- `examples/free_boundary_mono_planar_translation_2d.jl`
- `examples/free_boundary_mono_hydrostatic_equilibrium_2d.jl`
- `examples/free_boundary_diph_planar_translation_2d.jl`
- `examples/free_boundary_diph_muskat_perturbation_2d.jl`

Run tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Moving/free-boundary tests

- `test/test_free_mono_planar_translation.jl`
- `test/test_free_mono_hydrostatic_equilibrium.jl`
- `test/test_free_diph_planar_translation.jl`
- `test/test_free_diph_hydrostatic_equilibrium.jl`
- `test/test_height_tracker_roundtrip.jl`
- `test/test_free_boundary_residual_history.jl`
- `test/test_free_boundary_api_errors.jl`
