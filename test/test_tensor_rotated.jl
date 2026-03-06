@testset "mono Darcy rotated anisotropy tensor" begin
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

    function run_case(n)
        grid = (range(0.0, 1.0; length=n), range(0.0, 1.0; length=n))
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
        idx = active_indices(cap)
        err = sqrt(sum((pω[i] - p_exact(cap.C_ω[i]...))^2 for i in idx) / length(idx))

        Tf = PenguinDarcy.tensor_flux_operator(model, model.ops, model.cap; t=0.0, phase=1)
        nt = cap.ntotal
        cross12 = norm(Tf[1:nt, nt + 1:2 * nt])
        cross21 = norm(Tf[nt + 1:2 * nt, 1:nt])

        return err, step(grid[1]), cross12, cross21
    end

    e1, h1, c12_1, c21_1 = run_case(49)
    e2, h2, c12_2, c21_2 = run_case(97)

    @test c12_1 > 1e-10
    @test c21_1 > 1e-10
    @test c12_2 > 1e-10
    @test c21_2 > 1e-10

    order = log(e1 / e2) / log(h1 / h2)
    @test order > 0.9
end
