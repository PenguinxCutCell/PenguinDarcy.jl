@testset "embedded impermeable boundary + no-interface reduction" begin
    nx, ny = 41, 41
    grid = (range(0.0, 1.0; length=nx), range(0.0, 1.0; length=ny))

    radius = 0.18
    center = (0.5, 0.5)
    # Negative outside disk -> active fluid is outside obstacle.
    body(x, y) = -(hypot(x - center[1], y - center[2]) - radius)

    moms = geometric_moments(body, grid, Float64, nan; method=:vofijul)
    cap = assembled_capacity(moms; bc=0.0)

    p_far(x, y, t=0.0) = 1.0 - x
    bc = BorderConditions(
        ; left=Dirichlet(p_far), right=Dirichlet(p_far),
        bottom=Dirichlet(p_far), top=Dirichlet(p_far),
    )
    ops = DiffusionOps(cap; periodic=periodic_flags(bc, 2))

    model = DarcyModelMono(
        cap,
        ops,
        1.0;
        source=(x, y, t) -> 0.0,
        bc_border=bc,
        bc_interface=Robin(0.0, 1.0, 0.0),
    )

    sys = solve_steady!(model)
    lay = model.layout.offsets
    pω = sys.x[lay.ω]
    pγ = sys.x[lay.γ]

    λf = face_mobility_values(model; t=0.0)
    Wλ = spdiagm(0 => model.ops.Winv.nzval .* λf)
    qγ = model.ops.H' * (Wλ * (model.ops.G * pω + model.ops.H * pγ))

    idxγ = interface_indices(cap)
    @test !isempty(idxγ)
    @test maximum(abs.(qγ[idxγ])) < 5e-6

    mb = compute_mass_balance(model, sys.x)
    @test abs(mb.imbalance) < 5e-6

    # No-interface reduction: Robin interface option must be inert when Γ == 0.
    capf = assembled_capacity(full_moments(grid); bc=0.0)
    opsf = DiffusionOps(capf; periodic=periodic_flags(bc, 2))

    model_none = DarcyModelMono(capf, opsf, 1.0; source=(x, y, t) -> 0.0, bc_border=bc)
    model_rob = DarcyModelMono(
        capf,
        opsf,
        1.0;
        source=(x, y, t) -> 0.0,
        bc_border=bc,
        bc_interface=Robin(0.0, 1.0, 0.0),
    )

    sys_none = solve_steady!(model_none)
    sys_rob = solve_steady!(model_rob)

    p_none = sys_none.x[model_none.layout.offsets.ω]
    p_rob = sys_rob.x[model_rob.layout.offsets.ω]
    @test maximum(abs.(p_none .- p_rob)) < 1e-12
end
