using PenguinDarcy
using CartesianGrids
using PenguinBCs

ρ = 2.0
g = 1.5

grid = CartesianGrid((0.0, 0.0), (1.0, 1.0), (81, 65))
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
)

sol = solve_unsteady_moving!(model, (0.0, 0.05); dt=0.01, save_history=true)

drift = maximum(abs.(sol.interface_states[end] .- xf0))
max_vn = maximum(d.max_normal_velocity for d in sol.step_diagnostics)
hydro_res = maximum(abs, interface_mass_residual(sol))

println("max |V_n| over steps: ", max_vn)
println("max interface drift: ", drift)
println("hydrostatic mass-balance residual (max): ", hydro_res)
println("nonlinear iterations per step: ", [length(d.iteration.residual_inf) for d in sol.step_diagnostics])
