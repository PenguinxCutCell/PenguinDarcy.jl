@testset "moving mono free-boundary planar translation" begin
    U = 0.2
    h0 = 0.43
    tend = 0.05

    function run_case(dt; nx=41, ny=33)
        grid = CartesianGrid((0.0, 0.0), (1.0, 1.0), (nx, ny))
        y = collect(grid1d(grid, 2))
        xf0 = fill(h0, length(y))
        xf0[end] = xf0[1]

        tracker = HeightFunctionTracker(
            grid,
            xf0;
            axis=:x,
            periodic_transverse=true,
            damping=0.7,
            max_iter=25,
            tol_interface=1e-10,
            tol_update=1e-10,
        )

        bc = BorderConditions(
            ;
            left=Neumann(-U),
            right=Neumann(0.0),
            bottom=Periodic(),
            top=Periodic(),
        )

        model = MovingDarcyModelMono(
            tracker,
            1.0;
            source=(x, y, t) -> 0.0,
            bc_border=bc,
            p_ext=0.0,
            surface_tension=0.0,
        )

        sol = solve_unsteady_moving!(model, (0.0, tend); dt=dt, save_history=false)

        hvals = sol.interface_states[end]
        hnum = sum(hvals[1:(end - 1)]) / (length(hvals) - 1)
        hexact = h0 + U * tend
        herr = abs(hnum - hexact)

        asm = assemble_unsteady_moving!(model; t=tend)
        dmodel = asm.darcy_model
        lay = dmodel.layout.offsets
        pω = sol.system.x[lay.ω]

        idx = findall(i -> isfinite(dmodel.cap.buf.V[i]) && dmodel.cap.buf.V[i] > 0.0, 1:dmodel.cap.ntotal)
        pexact(x) = U * (hnum - x)
        perr = sqrt(sum((pω[i] - pexact(dmodel.cap.C_ω[i][1]))^2 for i in idx) / length(idx))

        mres = maximum(abs, interface_mass_residual(sol))
        return herr, perr, mres
    end

    e1, p1, m1 = run_case(0.02)
    e2, p2, m2 = run_case(0.01)
    e3, p3, m3 = run_case(0.005)

    @test e3 < 1e-6
    @test p3 < 1e-8
    @test maximum((m1, m2, m3)) < 2e-5

    order = log(e2 / e3) / log(0.01 / 0.005)
    @test order > 0.85

    @test p1 >= p2 >= p3
end
