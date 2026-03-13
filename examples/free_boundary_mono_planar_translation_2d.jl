using PenguinDarcy
using CartesianGrids
using PenguinBCs

U = 0.2
h0 = 0.43
tend = 0.05

grid = CartesianGrid((0.0, 0.0), (1.0, 1.0), (81, 65))
y = collect(grid1d(grid, 2))
xf0 = fill(h0, length(y))
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

sol = solve_unsteady_moving!(model, (0.0, tend); dt=0.01, save_history=true)

println("time  h_num        h_exact      abs_error")
for (k, t) in enumerate(sol.times)
    hvals = sol.interface_states[k]
    hnum = sum(hvals[1:(end - 1)]) / (length(hvals) - 1)
    hexact = h0 + U * t
    println(rpad(string(round(t, digits=4)), 6), " ",
            rpad(string(round(hnum, digits=10)), 13), " ",
            rpad(string(round(hexact, digits=10)), 13), " ",
            round(abs(hnum - hexact), digits=12))
end

asm = assemble_unsteady_moving!(model; t=tend)
dmodel = asm.darcy_model
lay = dmodel.layout.offsets
pω = sol.system.x[lay.ω]
idx = findall(i -> isfinite(dmodel.cap.buf.V[i]) && dmodel.cap.buf.V[i] > 0.0, 1:dmodel.cap.ntotal)
hend = sum(model.tracker.xf[1:(end - 1)]) / (length(model.tracker.xf) - 1)
pexact(x) = U * (hend - x)
perr = sqrt(sum((pω[i] - pexact(dmodel.cap.C_ω[i][1]))^2 for i in idx) / length(idx))

println("\nfinal pressure RMS error: ", perr)
println("final phase volume: ", sol.step_diagnostics[end].phase_volume)
println("max |global mass residual| over steps: ", maximum(abs, interface_mass_residual(sol)))
println("nonlinear iterations per step: ", [length(d.iteration.residual_inf) for d in sol.step_diagnostics])
