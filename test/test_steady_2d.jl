@testset "steady Darcy 2D linear gradient" begin
    grid = (range(0.0, 1.0; length=33), range(0.0, 1.0; length=33))
    cap = assembled_capacity(full_moments(grid); bc=0.0)

    a = 1.5
    b = -0.7
    λ = 2.2
    p_exact(x, y, t=0.0) = a * x + b * y

    bc = BorderConditions(
        ; left=Dirichlet(p_exact), right=Dirichlet(p_exact),
        bottom=Dirichlet(p_exact), top=Dirichlet(p_exact),
    )
    ops = DiffusionOps(cap; periodic=periodic_flags(bc, 2))
    model = DarcyModelMono(cap, ops, λ; source=(x, y, t) -> 0.0, bc_border=bc)
    sys = solve_steady!(model)

    lay = model.layout.offsets
    pω = sys.x[lay.ω]
    idx = active_indices(cap)
    err = sqrt(sum((pω[i] - p_exact(cap.C_ω[i]...))^2 for i in idx) / length(idx))
    @test err < 3e-10

    vel = recover_velocity(model, sys.x)
    ux_exact = -λ * a
    uy_exact = -λ * b
    LI = LinearIndices(cap.nnodes)
    idx_x = Int[]
    idx_y = Int[]
    for I in CartesianIndices(cap.nnodes)
        i = LI[I]
        halo = any(d -> I[d] == cap.nnodes[d], 1:2)
        halo && continue
        if cap.buf.V[i] > 0.0
            I[1] > 1 && push!(idx_x, i)
            I[2] > 1 && push!(idx_y, i)
        end
    end
    @test maximum(abs.(vel.x[idx_x] .- ux_exact)) < 2e-8
    @test maximum(abs.(vel.y[idx_y] .- uy_exact)) < 2e-8

    mb = compute_mass_balance(model, sys.x)
    @test abs(mb.source_integral) < 1e-12
    @test abs(mb.imbalance) < 1e-7
end
