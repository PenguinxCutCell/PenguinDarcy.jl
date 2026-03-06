using SparseArrays
using CartesianGeometry: geometric_moments, nan
using CartesianOperators
using PenguinBCs
using PenguinDarcy

nx, ny = 41, 41
grid = (range(0.0, 1.0; length=nx), range(0.0, 1.0; length=ny))

radius = 0.18
center = (0.5, 0.5)
body(x, y) = -(hypot(x - center[1], y - center[2]) - radius)

moms = geometric_moments(body, grid, Float64, nan; method=:vofijul)
cap = assembled_capacity(moms; bc=0.0)

p_far(x, y, t=0.0) = 1.0 - x
bc = BorderConditions(
    ; left=Dirichlet(p_far), right=Dirichlet(p_far),
    bottom=Dirichlet(p_far), top=Dirichlet(p_far),
)
ops = DiffusionOps(cap; periodic=periodic_flags(bc, 2))
model = DarcyModelMono(
    cap,
    ops,
    1.0;
    source=(x, y, t) -> 0.0,
    bc_border=bc,
    bc_interface=Robin(0.0, 1.0, 0.0),
)

sys = solve_steady!(model)
lay = model.layout.offsets
pω = sys.x[lay.ω]
pγ = sys.x[lay.γ]

λf = face_mobility_values(model; t=0.0)
Wλ = spdiagm(0 => model.ops.Winv.nzval .* λf)
qγ = model.ops.H' * (Wλ * (model.ops.G * pω + model.ops.H * pγ))
idxγ = findall(i -> isfinite(cap.buf.Γ[i]) && cap.buf.Γ[i] > 0.0, 1:cap.ntotal)

mb = compute_mass_balance(model, sys.x)

println("pressure extrema: min=", minimum(pω), ", max=", maximum(pω))
println("max impermeable interface residual |q|: ", maximum(abs.(qγ[idxγ])))
println("left discharge: ", boundary_discharge(model, sys.x, :left))
println("right discharge: ", boundary_discharge(model, sys.x, :right))
println("bottom discharge: ", boundary_discharge(model, sys.x, :bottom))
println("top discharge: ", boundary_discharge(model, sys.x, :top))
println("mass source integral: ", mb.source_integral)
println("mass imbalance: ", mb.imbalance)
