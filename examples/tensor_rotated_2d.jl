using LinearAlgebra
using CartesianGeometry: geometric_moments, nan
using CartesianOperators
using PenguinBCs
using PenguinDarcy

full_moments(grid) = geometric_moments((args...) -> -1.0, grid, Float64, nan; method=:vofijul)

θ = pi / 5
λ1 = 2.5
λ2 = 0.6
R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
Λ = R * Diagonal([λ1, λ2]) * transpose(R)

p_exact(x, y, t=0.0) = sin(pi * x) * sin(pi * y)
source(x, y, t=0.0) = begin
    pxx = -pi^2 * sin(pi * x) * sin(pi * y)
    pyy = pxx
    pxy = pi^2 * cos(pi * x) * cos(pi * y)
    -(Λ[1, 1] * pxx + (Λ[1, 2] + Λ[2, 1]) * pxy + Λ[2, 2] * pyy)
end

grid = (range(0.0, 1.0; length=97), range(0.0, 1.0; length=97))
cap = assembled_capacity(full_moments(grid); bc=0.0)

bc = BorderConditions(
    ; left=Dirichlet(p_exact), right=Dirichlet(p_exact),
    bottom=Dirichlet(p_exact), top=Dirichlet(p_exact),
)
ops = DiffusionOps(cap; periodic=periodic_flags(bc, 2))

model = DarcyModelMono(cap, ops, Λ; source=source, bc_border=bc)
sys = solve_steady!(model)
lay = model.layout.offsets
pω = sys.x[lay.ω]

idx = [i for i in 1:cap.ntotal if isfinite(cap.buf.V[i]) && cap.buf.V[i] > 0.0]
err = sqrt(sum((pω[i] - p_exact(cap.C_ω[i]...))^2 for i in idx) / length(idx))

Tf = PenguinDarcy.tensor_flux_operator(model, model.ops, model.cap; t=0.0, phase=1)
nt = cap.ntotal
cross12 = norm(Tf[1:nt, nt + 1:2 * nt])
cross21 = norm(Tf[nt + 1:2 * nt, 1:nt])

mb = mass_balance(model, sys.x)

println("pressure extrema: min=", minimum(pω), ", max=", maximum(pω))
println("L2 pressure error: ", err)
println("cross-coupling norms: T12=", cross12, ", T21=", cross21)
println("integrated distributed source: ", integrated_source(model; t=0.0))
println("integrated well rate: ", integrated_well_rate(model; t=0.0))
println("left boundary discharge: ", boundary_discharge(model, sys.x, :left))
println("right boundary discharge: ", boundary_discharge(model, sys.x, :right))
println("bottom boundary discharge: ", boundary_discharge(model, sys.x, :bottom))
println("top boundary discharge: ", boundary_discharge(model, sys.x, :top))
println("global mass-balance residual: ", mb.imbalance)
