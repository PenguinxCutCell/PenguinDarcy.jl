using CartesianGeometry: geometric_moments, nan
using CartesianOperators
using PenguinBCs
using PenguinDarcy

L = 1.0
ξ = 0.37
U0 = 1.0
UL = 0.0
λ1 = 1.5
λ2 = 0.7

grid = (range(0.0, L; length=161),)

moms1 = geometric_moments(x -> x - ξ, grid, Float64, nan; method=:vofijul)
moms2 = geometric_moments(x -> -(x - ξ), grid, Float64, nan; method=:vofijul)
cap1 = assembled_capacity(moms1; bc=0.0)
cap2 = assembled_capacity(moms2; bc=0.0)

bc = BorderConditions(; left=Dirichlet(U0), right=Dirichlet(UL))
ops1 = DiffusionOps(cap1; periodic=periodic_flags(bc, 1))
ops2 = DiffusionOps(cap2; periodic=periodic_flags(bc, 1))

model = DarcyModelDiph(
    cap1,
    ops1,
    λ1,
    cap2,
    ops2,
    λ2;
    source=((x, t) -> (0.0, 0.0)),
    bc_border=bc,
    bc_interface=DarcyContinuity(),
)

sys = solve_steady!(model)
lay = model.layout.offsets
p1ω = sys.x[lay.ω1]
p2ω = sys.x[lay.ω2]

q_exact = (U0 - UL) / (ξ / λ1 + (L - ξ) / λ2)
p1_exact(x) = U0 - (q_exact / λ1) * x
p2_exact(x) = UL + (q_exact / λ2) * (L - x)

idx1 = [i for i in 1:cap1.ntotal if isfinite(cap1.buf.V[i]) && cap1.buf.V[i] > 0.0]
idx2 = [i for i in 1:cap2.ntotal if isfinite(cap2.buf.V[i]) && cap2.buf.V[i] > 0.0]

err1 = sqrt(sum((p1ω[i] - p1_exact(cap1.C_ω[i][1]))^2 for i in idx1) / length(idx1))
err2 = sqrt(sum((p2ω[i] - p2_exact(cap2.C_ω[i][1]))^2 for i in idx2) / length(idx2))

iq = interface_discharge(model, sys.x)
mb = mass_balance(model, sys.x)

println("pressure extrema phase1: min=", minimum(p1ω), ", max=", maximum(p1ω))
println("pressure extrema phase2: min=", minimum(p2ω), ", max=", maximum(p2ω))
println("L2 errors: phase1=", err1, ", phase2=", err2)
println("integrated distributed source: ", integrated_source(model; t=0.0))
println("integrated well rate: ", integrated_well_rate(model; t=0.0))
println("left boundary discharge phase1: ", boundary_discharge(model, sys.x, :left; phase=1))
println("right boundary discharge phase2: ", boundary_discharge(model, sys.x, :right; phase=2))
println("interface discharge totals: phase1=", iq.total_phase1, ", phase2=", iq.total_phase2, ", balance=", iq.balance)
println("global mass-balance residual: ", mb.total.imbalance)
