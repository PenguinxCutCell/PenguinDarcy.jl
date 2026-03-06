@testset "steady Darcy 1D linear pressure" begin
    grid = (range(0.0, 1.0; length=65),)
    cap = assembled_capacity(full_moments(grid); bc=0.0)
    bc = BorderConditions(; left=Dirichlet(1.0), right=Dirichlet(0.0))
    ops = DiffusionOps(cap; periodic=periodic_flags(bc, 1))

    λ = 2.0
    model = DarcyModelMono(cap, ops, λ; source=(x, t) -> 0.0, bc_border=bc)
    sys = solve_steady!(model)

    lay = model.layout.offsets
    pω = sys.x[lay.ω]
    idx = active_indices(cap)
    err = sqrt(sum((pω[i] - (1.0 - cap.C_ω[i][1]))^2 for i in idx) / length(idx))
    @test err < 1e-10

    vel = recover_velocity(model, sys.x)
    vx = vel.x
    LI = LinearIndices(cap.nnodes)
    idx_v = Int[]
    for I in CartesianIndices(cap.nnodes)
        if I[1] > 1 && I[1] < cap.nnodes[1]
            i = LI[I]
            cap.buf.V[i] > 0.0 && push!(idx_v, i)
        end
    end
    vmaxerr = maximum(abs.(vx[idx_v] .- λ))
    @test vmaxerr < 2e-8

    qL = boundary_discharge(model, sys.x, :left)
    qR = boundary_discharge(model, sys.x, :right)
    @test qL ≈ -qR atol=1e-8

    mb = compute_mass_balance(model, sys.x)
    @test abs(mb.source_integral) < 1e-12
    @test abs(mb.imbalance) < 1e-8
end
