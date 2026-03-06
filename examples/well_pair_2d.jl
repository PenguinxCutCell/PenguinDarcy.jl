using CartesianGeometry: geometric_moments, nan
using CartesianOperators
using PenguinBCs
using PenguinDarcy

full_moments(grid) = geometric_moments((args...) -> -1.0, grid, Float64, nan; method=:vofijul)

grid = (range(0.0, 1.0; length=41), range(0.0, 1.0; length=41))
cap = assembled_capacity(full_moments(grid); bc=0.0)

bc = BorderConditions(
    ; left=Dirichlet(0.0), right=Dirichlet(0.0),
    bottom=Dirichlet(0.0), top=Dirichlet(0.0),
)
ops = DiffusionOps(cap; periodic=periodic_flags(bc, 2))

Q = 1.25
wells = [
    PointWell((0.30, 0.50), +Q; radius=0.08, phase=1),
    PointWell((0.70, 0.50), -Q; radius=0.08, phase=1),
]

model = DarcyModelMono(cap, ops, 1.0; source=(x, y, t) -> 0.0, wells=wells, bc_border=bc)
sys = solve_steady!(model)
lay = model.layout.offsets
pω = sys.x[lay.ω]

mb = mass_balance(model, sys.x)

println("pressure extrema: min=", minimum(pω), ", max=", maximum(pω))
println("integrated distributed source: ", integrated_source(model; t=0.0))
println("integrated well rate: ", integrated_well_rate(model; t=0.0))
println("left boundary discharge: ", boundary_discharge(model, sys.x, :left))
println("right boundary discharge: ", boundary_discharge(model, sys.x, :right))
println("bottom boundary discharge: ", boundary_discharge(model, sys.x, :bottom))
println("top boundary discharge: ", boundary_discharge(model, sys.x, :top))
println("global mass-balance residual: ", mb.imbalance)
