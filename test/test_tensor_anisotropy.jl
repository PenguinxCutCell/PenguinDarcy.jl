@testset "mono Darcy full-tensor anisotropy MMS" begin
    Λ = [2.0 0.4; 0.4 1.0]

    p_exact(x, y, t=0.0) = sin(pi * x) * sin(pi * y)
    source(x, y, t=0.0) = pi^2 * (3.0 * sin(pi * x) * sin(pi * y) - 0.8 * cos(pi * x) * cos(pi * y))

    px_exact(x, y) = pi * cos(pi * x) * sin(pi * y)
    py_exact(x, y) = pi * sin(pi * x) * cos(pi * y)
    ux_exact(x, y) = -(Λ[1, 1] * px_exact(x, y) + Λ[1, 2] * py_exact(x, y))
    uy_exact(x, y) = -(Λ[2, 1] * px_exact(x, y) + Λ[2, 2] * py_exact(x, y))

    function face_coord(cap, I, d)
        LI = LinearIndices(cap.nnodes)
        lin = LI[I]
        cω = cap.C_ω
        if I[d] == 1
            return (d == 1 ? cap.xyz[1][1] : cω[lin][1], d == 2 ? cap.xyz[2][1] : cω[lin][2])
        end
        Iminus = CartesianIndex(ntuple(k -> (k == d ? I[k] - 1 : I[k]), 2))
        linm = LI[Iminus]
        return (
            (cω[lin][1] + cω[linm][1]) / 2,
            (cω[lin][2] + cω[linm][2]) / 2,
        )
    end

    function run_case(n)
        grid = (range(0.0, 1.0; length=n), range(0.0, 1.0; length=n))
        cap = assembled_capacity(full_moments(grid); bc=0.0)

        bc = BorderConditions(
            ; left=Dirichlet(p_exact), right=Dirichlet(p_exact),
            bottom=Dirichlet(p_exact), top=Dirichlet(p_exact),
        )
        ops = DiffusionOps(cap; periodic=periodic_flags(bc, 2))
        model = DarcyModelMono(cap, ops, Λ; source=source, bc_border=bc)

        # Ensure off-diagonal tensor couplings are present in the discrete flux map.
        Tf = PenguinDarcy.tensor_flux_operator(model, model.ops, model.cap; t=0.0, phase=1)
        nt = cap.ntotal
        @test norm(Tf[1:nt, nt + 1:2 * nt]) > 1e-10
        @test norm(Tf[nt + 1:2 * nt, 1:nt]) > 1e-10

        sys = solve_steady!(model)
        lay = model.layout.offsets
        pω = sys.x[lay.ω]
        idx = active_indices(cap)

        errp = sqrt(sum((pω[i] - p_exact(cap.C_ω[i]...))^2 for i in idx) / length(idx))

        vel = recover_velocity(model, sys.x)
        LI = LinearIndices(cap.nnodes)
        errux = 0.0
        erruy = 0.0
        nx = 0
        ny = 0

        for I in CartesianIndices(cap.nnodes)
            i = LI[I]
            halo = any(d -> I[d] == cap.nnodes[d], 1:2)
            halo && continue
            cap.buf.V[i] > 0.0 || continue

            if I[1] > 1
                x, y = face_coord(cap, I, 1)
                errux += (vel.x[i] - ux_exact(x, y))^2
                nx += 1
            end
            if I[2] > 1
                x, y = face_coord(cap, I, 2)
                erruy += (vel.y[i] - uy_exact(x, y))^2
                ny += 1
            end
        end

        erru = sqrt(errux / max(nx, 1) + erruy / max(ny, 1))
        return errp, erru, step(grid[1])
    end

    errs_p = Float64[]
    errs_u = Float64[]
    hs = Float64[]
    for n in (33, 65)
        ep, eu, h = run_case(n)
        push!(errs_p, ep)
        push!(errs_u, eu)
        push!(hs, h)
    end

    order_p = log(errs_p[1] / errs_p[2]) / log(hs[1] / hs[2])
    order_u = log(errs_u[1] / errs_u[2]) / log(hs[1] / hs[2])

    @test order_p > 0.95
    @test order_u > 0.45
end
