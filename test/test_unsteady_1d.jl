@testset "unsteady Darcy 1D MMS (BE/CN + spatial)" begin
    λ = 1.3
    S = λ * pi^2

    p_exact(x, t) = exp(-t) * sin(pi * x)
    source(x, t) = 0.0

    function run_error(n, dt, scheme, tend)
        grid = (range(0.0, 1.0; length=n),)
        cap = assembled_capacity(full_moments(grid); bc=0.0)
        bc = BorderConditions(; left=Dirichlet(0.0), right=Dirichlet(0.0))
        ops = DiffusionOps(cap; periodic=periodic_flags(bc, 1))
        model = DarcyModelMono(cap, ops, λ; source=source, storage=S, bc_border=bc)

        u0 = [p_exact(cap.C_ω[i][1], 0.0) for i in 1:cap.ntotal]
        sol = solve_unsteady!(model, u0, (0.0, tend); dt=dt, scheme=scheme, save_history=false)

        lay = model.layout.offsets
        pω = sol.states[end][lay.ω]
        idx = active_indices(cap)
        err = sqrt(sum((pω[i] - p_exact(cap.C_ω[i][1], tend))^2 for i in idx) / length(idx))
        return err, step(grid[1])
    end

    errs_be = Float64[]
    dts_be = [0.1, 0.05, 0.025]
    for dt in dts_be
        e, _ = run_error(101, dt, :BE, 0.2)
        push!(errs_be, e)
    end
    order_be = log(errs_be[end - 1] / errs_be[end]) / log(dts_be[end - 1] / dts_be[end])
    @test order_be > 0.9

    errs_cn = Float64[]
    dts_cn = [0.2, 0.1, 0.05]
    for dt in dts_cn
        e, _ = run_error(101, dt, :CN, 0.2)
        push!(errs_cn, e)
    end
    order_cn = log(errs_cn[end - 1] / errs_cn[end]) / log(dts_cn[end - 1] / dts_cn[end])
    @test order_cn > 2.0

    errs_h = Float64[]
    hs = Float64[]
    for n in (41, 81)
        e, h = run_error(n, 2e-4, :CN, 0.02)
        push!(errs_h, e)
        push!(hs, h)
    end
    order_h = log(errs_h[1] / errs_h[2]) / log(hs[1] / hs[2])
    @test order_h > 1.4
end
