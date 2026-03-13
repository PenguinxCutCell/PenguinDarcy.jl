@testset "moving free-boundary residual history" begin
    grid = CartesianGrid((0.0, 0.0), (1.0, 1.0), (41, 33))
    y = collect(grid1d(grid, 2))
    xf0 = fill(0.43, length(y))
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
        left=Neumann(-0.2),
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

    sol = solve_unsteady_moving!(model, (0.0, 0.01); dt=0.01, save_history=false)

    hist = sol.step_diagnostics[end].iteration
    @test length(hist.residual_inf) >= 2
    @test hist.residual_inf[end] < hist.residual_inf[1]
    @test hist.update_inf[end] < hist.update_inf[1]
    @test hist.converged

    # Must stop because tolerances were met, not because max_iter was exhausted.
    @test length(hist.residual_inf) < tracker.max_iter
    @test hist.residual_inf[end] < 10 * tracker.tol_interface
end
