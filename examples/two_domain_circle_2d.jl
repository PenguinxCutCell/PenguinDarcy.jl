using CartesianGeometry: geometric_moments, nan
using CartesianOperators
using PenguinBCs
using PenguinDarcy

nx, ny = 81, 81
grid = (range(0.0, 1.0; length=nx), range(0.0, 1.0; length=ny))

radius = 0.22
center = (0.5, 0.5)
ϕ1(x, y) = hypot(x - center[1], y - center[2]) - radius
ϕ2(x, y) = -ϕ1(x, y)

cap1 = assembled_capacity(geometric_moments(ϕ1, grid, Float64, nan; method=:vofijul); bc=0.0)
cap2 = assembled_capacity(geometric_moments(ϕ2, grid, Float64, nan; method=:vofijul); bc=0.0)

p_far(x, y, t=0.0) = 1.0 - x
bc = BorderConditions(
    ; left=Dirichlet(p_far), right=Dirichlet(p_far),
    bottom=Dirichlet(p_far), top=Dirichlet(p_far),
)

ops1 = DiffusionOps(cap1; periodic=periodic_flags(bc, 2))
ops2 = DiffusionOps(cap2; periodic=periodic_flags(bc, 2))

model = DarcyModelDiph(
    cap1,
    ops1,
    1.0,
    cap2,
    ops2,
    0.7;
    source=((x, y, t) -> (0.0, 0.0)),
    bc_border=bc,
    bc_interface=DarcyContinuity(),
)

sys = solve_steady!(model)
lay = model.layout.offsets
p1ω = sys.x[lay.ω1]
p2ω = sys.x[lay.ω2]

iq = interface_discharge(model, sys.x)
mb = mass_balance(model, sys.x)

println("pressure extrema phase1: min=", minimum(p1ω), ", max=", maximum(p1ω))
println("pressure extrema phase2: min=", minimum(p2ω), ", max=", maximum(p2ω))
println("integrated distributed source: ", integrated_source(model; t=0.0))
println("integrated well rate: ", integrated_well_rate(model; t=0.0))
println("left boundary discharge phase1: ", boundary_discharge(model, sys.x, :left; phase=1))
println("right boundary discharge phase1: ", boundary_discharge(model, sys.x, :right; phase=1))
println("left boundary discharge phase2: ", boundary_discharge(model, sys.x, :left; phase=2))
println("right boundary discharge phase2: ", boundary_discharge(model, sys.x, :right; phase=2))
println("bottom boundary discharge phase1: ", boundary_discharge(model, sys.x, :bottom; phase=1))
println("top boundary discharge phase1: ", boundary_discharge(model, sys.x, :top; phase=1))
println("interface discharge totals: phase1=", iq.total_phase1, ", phase2=", iq.total_phase2, ", balance=", iq.balance)
println("global mass-balance residual: ", mb.total.imbalance)
