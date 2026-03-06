@testset "two-domain Darcy steady 1D transmission" begin
    L = 1.0
    ξ = 0.37
    U0 = 1.0
    UL = 0.0
    λ1 = 1.5
    λ2 = 0.7

    q = (U0 - UL) / (ξ / λ1 + (L - ξ) / λ2)
    p1_exact(x) = U0 - (q / λ1) * x
    p2_exact(x) = UL + (q / λ2) * (L - x)

    function run_case(n)
        grid = (range(0.0, L; length=n),)
        moms1 = geometric_moments(x -> x - ξ, grid, Float64, nan; method=:vofijul)
        moms2 = geometric_moments(x -> -(x - ξ), grid, Float64, nan; method=:vofijul)

        cap1 = assembled_capacity(moms1; bc=0.0)
        cap2 = assembled_capacity(moms2; bc=0.0)

        bc = BorderConditions(; left=Dirichlet(U0), right=Dirichlet(UL))
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
            bc_border=bc,
            bc_interface=DarcyContinuity(),
        )

        sys = solve_steady!(model)
        lay = model.layout.offsets

        p1ω = sys.x[lay.ω1]
        p1γ = sys.x[lay.γ1]
        p2ω = sys.x[lay.ω2]
        p2γ = sys.x[lay.γ2]

        idx1 = active_indices(cap1)
        idx2 = active_indices(cap2)
        err1 = sqrt(sum((p1ω[i] - p1_exact(cap1.C_ω[i][1]))^2 for i in idx1) / length(idx1))
        err2 = sqrt(sum((p2ω[i] - p2_exact(cap2.C_ω[i][1]))^2 for i in idx2) / length(idx2))

        idxγ = interface_indices(cap1)
        @test !isempty(idxγ)

        vel = recover_flux(model, sys.x)
        q1 = model.ops1.H' * vel.phase1.faces
        q2 = model.ops2.H' * vel.phase2.faces

        jump_p = maximum(abs.(p1γ[idxγ] .- p2γ[idxγ]))
        cont_q = maximum(abs.(q1[idxγ] .+ q2[idxγ]))

        return err1, err2, step(grid[1]), jump_p, cont_q
    end

    errs1 = Float64[]
    errs2 = Float64[]
    pjump = Float64[]
    qjump = Float64[]

    for n in (41, 81, 161)
        e1, e2, _, jp, jq = run_case(n)
        push!(errs1, e1)
        push!(errs2, e2)
        push!(pjump, jp)
        push!(qjump, jq)
    end

    @test maximum(errs1) < 1e-10
    @test maximum(errs2) < 1e-10
    @test maximum(pjump) < 3e-2
    @test maximum(qjump) < 3e-2
end
