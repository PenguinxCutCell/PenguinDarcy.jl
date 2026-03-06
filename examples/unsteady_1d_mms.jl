using CartesianGeometry: geometric_moments, nan
using CartesianOperators
using PenguinBCs
using PenguinDarcy

λ = 1.3
S = λ * pi^2
p_exact(x, t) = exp(-t) * sin(pi * x)
source(x, t) = 0.0

function run_case(dt, scheme)
    grid = (range(0.0, 1.0; length=101),)
    moms = geometric_moments((x, t=0.0) -> -1.0, grid, Float64, nan; method=:vofijul)
    cap = assembled_capacity(moms; bc=0.0)
    bc = BorderConditions(; left=Dirichlet(0.0), right=Dirichlet(0.0))
    ops = DiffusionOps(cap; periodic=periodic_flags(bc, 1))
    model = DarcyModelMono(cap, ops, λ; source=source, storage=S, bc_border=bc)

    u0 = [p_exact(cap.C_ω[i][1], 0.0) for i in 1:cap.ntotal]
    sol = solve_unsteady!(model, u0, (0.0, 0.2); dt=dt, scheme=scheme, save_history=false)

    lay = model.layout.offsets
    pω = sol.states[end][lay.ω]
    idx = [i for i in 1:cap.ntotal if cap.buf.V[i] > 0.0]
    err = sqrt(sum((pω[i] - p_exact(cap.C_ω[i][1], 0.2))^2 for i in idx) / length(idx))
    return err
end

for (scheme, dts) in ((:BE, (0.1, 0.05, 0.025)), (:CN, (0.2, 0.1, 0.05)))
    errs = Float64[]
    for dt in dts
        e = run_case(dt, scheme)
        push!(errs, e)
        println("scheme=", scheme, ", dt=", dt, ", L2 error=", e)
    end
    order = log(errs[end - 1] / errs[end]) / log(2.0)
    println("scheme=", scheme, ", last-step observed temporal order=", order)
end
