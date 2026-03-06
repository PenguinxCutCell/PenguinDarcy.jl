@testset "mono Darcy hydrostatic gravity equilibrium" begin
    ρ = 2.3
    g = 1.7
    λ = 1.9

    grid = (range(0.0, 1.0; length=81),)
    cap = assembled_capacity(full_moments(grid); bc=0.0)

    p_eq(x, t=0.0) = ρ * g * x
    bc = BorderConditions(; left=Dirichlet(p_eq), right=Dirichlet(p_eq))
    ops = DiffusionOps(cap; periodic=periodic_flags(bc, 1))

    model = DarcyModelMono(
        cap,
        ops,
        λ;
        source=(x, t) -> 0.0,
        ρ=ρ,
        gravity=(g,),
        bc_border=bc,
    )

    sys = solve_steady!(model)
    lay = model.layout.offsets
    pω = sys.x[lay.ω]

    idx = active_indices(cap)
    errp = sqrt(sum((pω[i] - p_eq(cap.C_ω[i][1]))^2 for i in idx) / length(idx))
    @test errp < 1e-10

    vel = recover_velocity(model, sys.x)
    LI = LinearIndices(cap.nnodes)
    idx_v = Int[]
    for I in CartesianIndices(cap.nnodes)
        i = LI[I]
        halo = any(d -> I[d] == cap.nnodes[d], 1:1)
        halo && continue
        if I[1] > 1 && cap.buf.V[i] > 0.0
            push!(idx_v, i)
        end
    end
    @test maximum(abs.(vel.x[idx_v])) < 2e-8

    mb = mass_balance(model, sys.x)
    @test abs(mb.source_integral) < 1e-12
    @test abs(mb.well_rate) < 1e-12
    @test abs(mb.imbalance) < 1e-8
end
