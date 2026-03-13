using PenguinDarcy
using CartesianGrids
using PenguinBCs
using LinearAlgebra

q = 0.15
λ1 = 1.0
λ2 = 0.4
J = 0.05
h0 = 0.43
tend = 0.03

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
    left=Dirichlet(1.0),
    right=Neumann(q),
    bottom=Periodic(),
    top=Periodic(),
)

model = MovingDarcyModelDiph(
    tracker,
    λ1,
    λ2;
    source=(x, y, t) -> (0.0, 0.0),
    bc_border=bc,
    bc_interface=DarcyContinuity(J, 0.0),
)

sol = solve_unsteady_moving!(model, (0.0, tend); dt=0.01, save_history=false)

hnum = sum(model.tracker.xf[1:(end - 1)]) / (length(model.tracker.xf) - 1)
hexact = h0 + q * tend
println("interface position: num = ", hnum, ", exact = ", hexact, ", abs err = ", abs(hnum - hexact))

asm = assemble_unsteady_moving!(model; t=tend)
dmodel = asm.darcy_model
lay = dmodel.layout.offsets

p1ω = sol.system.x[lay.ω1]
p2ω = sol.system.x[lay.ω2]
idx1 = findall(i -> isfinite(dmodel.cap1.buf.V[i]) && dmodel.cap1.buf.V[i] > 0.0, 1:dmodel.cap1.ntotal)
idx2 = findall(i -> isfinite(dmodel.cap2.buf.V[i]) && dmodel.cap2.buf.V[i] > 0.0, 1:dmodel.cap2.ntotal)

x1 = [dmodel.cap1.C_ω[i][1] for i in idx1]
y1 = [p1ω[i] for i in idx1]
x2 = [dmodel.cap2.C_ω[i][1] for i in idx2]
y2 = [p2ω[i] for i in idx2]

c1 = [ones(length(x1)) x1] \ y1
c2 = [ones(length(x2)) x2] \ y2
println("phase-1 slope dp/dx: ", c1[2], " (expected ", -q / λ1, ")")
println("phase-2 slope dp/dx: ", c2[2], " (expected ", -q / λ2, ")")

idxγ = findall(i -> isfinite(dmodel.cap1.buf.Γ[i]) && dmodel.cap1.buf.Γ[i] > 0.0, 1:dmodel.cap1.ntotal)
p1γ = sol.system.x[lay.γ1]
p2γ = sol.system.x[lay.γ2]
println("max pressure-jump error: ", maximum(abs.(p1γ[idxγ] .- p2γ[idxγ] .- J)))
println("max normal-velocity jump: ", sol.step_diagnostics[end].max_normal_velocity_jump)
println("phase volumes: ", sol.step_diagnostics[end].phase_volume)
println("nonlinear iterations per step: ", [length(d.iteration.residual_inf) for d in sol.step_diagnostics])
