using CartesianGeometry: geometric_moments, nan
using CartesianOperators
using PenguinBCs
using PenguinDarcy

grid = (range(0.0, 1.0; length=65),)
moms = geometric_moments((x, t=0.0) -> -1.0, grid, Float64, nan; method=:vofijul)
cap = assembled_capacity(moms; bc=0.0)

bc = BorderConditions(; left=Dirichlet(1.0), right=Dirichlet(0.0))
ops = DiffusionOps(cap; periodic=periodic_flags(bc, 1))

model = DarcyModelMono(cap, ops, 2.0; source=(x, t) -> 0.0, bc_border=bc)
sys = solve_steady!(model)

lay = model.layout.offsets
pω = sys.x[lay.ω]
vel = recover_velocity(model, sys.x)
mb = compute_mass_balance(model, sys.x)

println("pressure extrema: min=", minimum(pω), ", max=", maximum(pω))
println("velocity mean x: ", sum(vel.x) / length(vel.x))
println("left discharge: ", boundary_discharge(model, sys.x, :left))
println("right discharge: ", boundary_discharge(model, sys.x, :right))
println("mass source integral: ", mb.source_integral)
println("mass imbalance: ", mb.imbalance)
