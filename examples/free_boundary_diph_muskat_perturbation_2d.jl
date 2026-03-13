using PenguinDarcy
using CartesianGrids
using PenguinBCs

function run_case(label; ρ1, ρ2, g, λ1=1.0, λ2=1.0, tend=0.05, dt=0.01)
    grid = CartesianGrid((0.0, 0.0), (1.0, 1.0), (65, 49))
    y = collect(grid1d(grid, 2))

    h0 = 0.45
    a0 = 0.02
    xf0 = [h0 + a0 * sin(2 * pi * yy) for yy in y]
    xf0[end] = xf0[1]

    tracker = HeightFunctionTracker(
        grid,
        xf0;
        axis=:x,
        periodic_transverse=true,
        damping=0.6,
        max_iter=25,
        tol_interface=1e-9,
        tol_update=1e-9,
    )

    bc = BorderConditions(
        ;
        left=Dirichlet(0.0),
        right=Dirichlet(0.0),
        bottom=Periodic(),
        top=Periodic(),
    )

    model = MovingDarcyModelDiph(
        tracker,
        λ1,
        λ2;
        source=(x, y, t) -> (0.0, 0.0),
        bc_border=bc,
        bc_interface=DarcyContinuity(),
        ρ=(ρ1, ρ2),
        gravity=(g, 0.0),
    )

    sol = solve_unsteady_moving!(model, (0.0, tend); dt=dt, save_history=true)

    amps = Float64[]
    for xf in sol.interface_states
        vals = xf[1:(end - 1)]
        push!(amps, (maximum(vals) - minimum(vals)) / 2)
    end

    vjump = maximum(d.max_normal_velocity_jump for d in sol.step_diagnostics)
    vols1 = [d.phase_volume[1] for d in sol.step_diagnostics]
    vols2 = [d.phase_volume[2] for d in sol.step_diagnostics]
    iters = [length(d.iteration.residual_inf) for d in sol.step_diagnostics]

    println("\n=== ", label, " ===")
    println("density pair (ρ1, ρ2) = (", ρ1, ", ", ρ2, ")")
    println("amplitude history = ", amps)
    println("amplitude change  = ", amps[end] - amps[1])
    println("phase-1 volume range = ", (minimum(vols1), maximum(vols1)))
    println("phase-2 volume range = ", (minimum(vols2), maximum(vols2)))
    println("max normal-velocity jump = ", vjump)
    println("nonlinear iterations per step = ", iters)

    return sol
end

stable = run_case("stable", ρ1=1.0, ρ2=2.0, g=1.0)
unstable = run_case("unstable", ρ1=2.0, ρ2=1.0, g=1.0)

println("\nStable final amplitude   : ", (maximum(stable.interface_states[end][1:end-1]) - minimum(stable.interface_states[end][1:end-1])) / 2)
println("Unstable final amplitude : ", (maximum(unstable.interface_states[end][1:end-1]) - minimum(unstable.interface_states[end][1:end-1])) / 2)
