@testset "two-domain Darcy unsteady 1D BE/CN" begin
    L = 1.0
    ξ = 0.37
    λ1 = 1.2
    λ2 = 1.2
    S1 = λ1 * pi^2
    S2 = λ2 * pi^2

    p_exact(x, t) = exp(-t) * sin(pi * x)

    function run_error(n, dt, scheme, tend)
        grid = (range(0.0, L; length=n),)
        moms1 = geometric_moments(x -> x - ξ, grid, Float64, nan; method=:vofijul)
        moms2 = geometric_moments(x -> -(x - ξ), grid, Float64, nan; method=:vofijul)

        cap1 = assembled_capacity(moms1; bc=0.0)
        cap2 = assembled_capacity(moms2; bc=0.0)

        bc = BorderConditions(; left=Dirichlet(0.0), right=Dirichlet(0.0))
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
            storage=(S1, S2),
            bc_border=bc,
            bc_interface=DarcyContinuity(),
        )

        u01 = [p_exact(cap1.C_ω[i][1], 0.0) for i in 1:cap1.ntotal]
        u02 = [p_exact(cap2.C_ω[i][1], 0.0) for i in 1:cap2.ntotal]
        u0 = vcat(u01, u02)

        sol = solve_unsteady!(model, u0, (0.0, tend); dt=dt, scheme=scheme, save_history=false)

        lay = model.layout.offsets
        p1ω = sol.states[end][lay.ω1]
        p2ω = sol.states[end][lay.ω2]

        idx1 = active_indices(cap1)
        idx2 = active_indices(cap2)
        e1 = sqrt(sum((p1ω[i] - p_exact(cap1.C_ω[i][1], tend))^2 for i in idx1) / length(idx1))
        e2 = sqrt(sum((p2ω[i] - p_exact(cap2.C_ω[i][1], tend))^2 for i in idx2) / length(idx2))

        e = sqrt((e1^2 + e2^2) / 2)
        return e, step(grid[1])
    end

    tend = 0.2

    errs_be = Float64[]
    dts_be = [0.1, 0.05, 0.025]
    for dt in dts_be
        e, _ = run_error(81, dt, :BE, tend)
        push!(errs_be, e)
    end
    order_be = log(errs_be[end - 1] / errs_be[end]) / log(dts_be[end - 1] / dts_be[end])
    @test order_be > 0.85

    errs_cn = Float64[]
    dts_cn = [0.2, 0.1, 0.05]
    for dt in dts_cn
        e, _ = run_error(81, dt, :CN, tend)
        push!(errs_cn, e)
    end
    order_cn = log(errs_cn[end - 1] / errs_cn[end]) / log(dts_cn[end - 1] / dts_cn[end])
    @test order_cn > 0.45

    errs_h = Float64[]
    hs = Float64[]
    for n in (41, 81)
        e, h = run_error(n, 2e-4, :CN, 0.02)
        push!(errs_h, e)
        push!(hs, h)
    end
    order_h = log(errs_h[1] / errs_h[2]) / log(hs[1] / hs[2])
    @test order_h > 1.8
end
