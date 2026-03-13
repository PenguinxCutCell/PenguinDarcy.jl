@testset "moving mono hydrostatic free-surface equilibrium" begin
    ρ = 2.0
    g = 1.5

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
        max_iter=15,
        tol_interface=1e-10,
        tol_update=1e-10,
    )

    p_h(x, y, t=0.0) = ρ * g * x
    bc = BorderConditions(
        ;
        left=Dirichlet(p_h),
        right=Dirichlet(p_h),
        bottom=Periodic(),
        top=Periodic(),
    )

    model = MovingDarcyModelMono(
        tracker,
        1.2;
        source=(x, y, t) -> 0.0,
        ρ=ρ,
        gravity=(g, 0.0),
        bc_border=bc,
        p_ext=p_h,
        surface_tension=0.0,
    )

    sol = solve_unsteady_moving!(model, (0.0, 0.05); dt=0.01, save_history=false)

    drift = maximum(abs.(sol.interface_states[end] .- xf0))
    vmax = maximum(d.max_normal_velocity for d in sol.step_diagnostics)

    @test drift < 1e-10
    @test vmax < 5e-9
end
