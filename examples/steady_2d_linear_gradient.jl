using CartesianGeometry: geometric_moments, nan
using CartesianOperators
using PenguinBCs
using PenguinDarcy

grid = (range(0.0, 1.0; length=41), range(0.0, 1.0; length=41))
moms = geometric_moments((x, y, t=0.0) -> -1.0, grid, Float64, nan; method=:vofijul)
cap = assembled_capacity(moms; bc=0.0)

a, b = 1.2, -0.6
λ = 1.7
p_exact(x, y, t=0.0) = a * x + b * y

bc = BorderConditions(
    ; left=Dirichlet(p_exact), right=Dirichlet(p_exact),
    bottom=Dirichlet(p_exact), top=Dirichlet(p_exact),
)
ops = DiffusionOps(cap; periodic=periodic_flags(bc, 2))
model = DarcyModelMono(cap, ops, λ; source=(x, y, t) -> 0.0, bc_border=bc)

sys = solve_steady!(model)
lay = model.layout.offsets
pω = sys.x[lay.ω]
vel = recover_velocity(model, sys.x)
mb = compute_mass_balance(model, sys.x)

println("pressure extrema: min=", minimum(pω), ", max=", maximum(pω))
println("velocity mean x: ", sum(vel.x) / length(vel.x), " (exact ", -λ * a, ")")
println("velocity mean y: ", sum(vel.y) / length(vel.y), " (exact ", -λ * b, ")")
println("left discharge: ", boundary_discharge(model, sys.x, :left))
println("right discharge: ", boundary_discharge(model, sys.x, :right))
println("bottom discharge: ", boundary_discharge(model, sys.x, :bottom))
println("top discharge: ", boundary_discharge(model, sys.x, :top))
println("mass imbalance: ", mb.imbalance)
