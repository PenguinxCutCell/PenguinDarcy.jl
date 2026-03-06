@testset "two-domain Darcy hydrostatic gravity equilibrium" begin
    L = 1.0
    ξ = 0.437
    λ1 = 1.3
    λ2 = 0.6
    ρ = 2.0
    g = 2.4

    p_eq(x, t=0.0) = ρ * g * x

    grid = (range(0.0, L; length=101),)
    moms1 = geometric_moments(x -> x - ξ, grid, Float64, nan; method=:vofijul)
    moms2 = geometric_moments(x -> -(x - ξ), grid, Float64, nan; method=:vofijul)

    cap1 = assembled_capacity(moms1; bc=0.0)
    cap2 = assembled_capacity(moms2; bc=0.0)

    bc = BorderConditions(; left=Dirichlet(p_eq), right=Dirichlet(p_eq))
    ops1 = DiffusionOps(cap1; periodic=periodic_flags(bc, 1))
    ops2 = DiffusionOps(cap2; periodic=periodic_flags(bc, 1))

    model = DarcyModelDiph(
        cap1,
        ops1,
        λ1,
        cap2,
        ops2,
        λ2;
        source=((x, t) -> (0.0, 0.0)),
        ρ=(ρ, ρ),
        gravity=(g,),
        bc_border=bc,
        bc_interface=DarcyContinuity(),
    )

    sys = solve_steady!(model)
    lay = model.layout.offsets
    p1ω = sys.x[lay.ω1]
    p2ω = sys.x[lay.ω2]

    idx1 = active_indices(cap1)
    idx2 = active_indices(cap2)
    e1 = sqrt(sum((p1ω[i] - p_eq(cap1.C_ω[i][1]))^2 for i in idx1) / length(idx1))
    e2 = sqrt(sum((p2ω[i] - p_eq(cap2.C_ω[i][1]))^2 for i in idx2) / length(idx2))
    @test e1 < 1e-9
    @test e2 < 1e-9

    flux = recover_flux(model, sys.x)
    LI = LinearIndices(cap1.nnodes)
    idxf1 = Int[]
    idxf2 = Int[]
    for I in CartesianIndices(cap1.nnodes)
        i = LI[I]
        halo = any(d -> I[d] == cap1.nnodes[d], 1:1)
        halo && continue
        if I[1] > 1 && cap1.buf.V[i] > 0.0
            push!(idxf1, i)
        end
        if I[1] > 1 && cap2.buf.V[i] > 0.0
            push!(idxf2, i)
        end
    end

    @test maximum(abs.(flux.phase1.x[idxf1])) < 2e-8
    @test maximum(abs.(flux.phase2.x[idxf2])) < 2e-8

    iq = interface_discharge(model, sys.x)
    @test abs(iq.balance) < 1e-8

    mb = mass_balance(model, sys.x)
    @test abs(mb.total.imbalance) < 1e-8
end
