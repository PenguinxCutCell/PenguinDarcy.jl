@testset "moving diphasic hydrostatic equilibrium" begin
    ρ = 1.7
    g = 2.0

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
        max_iter=20,
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

    model = MovingDarcyModelDiph(
        tracker,
        1.0,
        0.5;
        source=(x, y, t) -> (0.0, 0.0),
        bc_border=bc,
        bc_interface=DarcyContinuity(0.0, 0.0),
        ρ=(ρ, ρ),
        gravity=(g, 0.0),
        surface_tension=0.0,
    )

    sol = solve_unsteady_moving!(model, (0.0, 0.05); dt=0.01, save_history=false)

    drift = maximum(abs.(sol.interface_states[end] .- xf0))
    vmax = maximum(d.max_normal_velocity for d in sol.step_diagnostics)
    vjump = maximum(d.max_normal_velocity_jump for d in sol.step_diagnostics)

    @test drift < 1e-10
    @test vmax < 1e-9
    @test vjump < 1e-9

    asm = assemble_unsteady_moving!(model; t=0.05)
    dmodel = asm.darcy_model
    lay = dmodel.layout.offsets
    p1γ = sol.system.x[lay.γ1]
    p2γ = sol.system.x[lay.γ2]
    idxγ = findall(i -> isfinite(dmodel.cap1.buf.Γ[i]) && dmodel.cap1.buf.Γ[i] > 0.0, 1:dmodel.cap1.ntotal)
    @test maximum(abs.(p1γ[idxγ] .- p2γ[idxγ])) < 1e-10
end
