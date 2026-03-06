using CartesianGeometry: geometric_moments, nan
using CartesianOperators
using PenguinBCs
using PenguinDarcy

full_moments(grid) = geometric_moments((args...) -> -1.0, grid, Float64, nan; method=:vofijul)

ρ = 2.3
g = 1.7
λ = 1.9

grid = (range(0.0, 1.0; length=81),)
cap = assembled_capacity(full_moments(grid); bc=0.0)

p_eq(x, t=0.0) = ρ * g * x
bc = BorderConditions(; left=Dirichlet(p_eq), right=Dirichlet(p_eq))
ops = DiffusionOps(cap; periodic=periodic_flags(bc, 1))

model = DarcyModelMono(
    cap,
    ops,
    λ;
    source=(x, t) -> 0.0,
    ρ=ρ,
    gravity=(g,),
    bc_border=bc,
)

sys = solve_steady!(model)
lay = model.layout.offsets
pω = sys.x[lay.ω]
vel = recover_velocity(model, sys.x)
LI = LinearIndices(cap.nnodes)
idx_int = Int[]
for I in CartesianIndices(cap.nnodes)
    i = LI[I]
    if I[1] > 1 && I[1] < cap.nnodes[1] && cap.buf.V[i] > 0.0
        push!(idx_int, i)
    end
end

mb = mass_balance(model, sys.x)

println("pressure extrema: min=", minimum(pω), ", max=", maximum(pω))
println("max |u_x| interior: ", maximum(abs.(vel.x[idx_int])))
println("integrated distributed source: ", integrated_source(model; t=0.0))
println("integrated well rate: ", integrated_well_rate(model; t=0.0))
println("left boundary discharge: ", boundary_discharge(model, sys.x, :left))
println("right boundary discharge: ", boundary_discharge(model, sys.x, :right))
println("global mass-balance residual: ", mb.imbalance)
