@testset "mono Darcy conservative wells" begin
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

    @test abs(integrated_well_rate(model; t=0.0)) < 1e-12

    sys = solve_steady!(model)
    mb = mass_balance(model, sys.x)

    @test abs(mb.source_integral) < 1e-12
    @test abs(mb.well_rate) < 1e-12
    @test abs(mb.boundary_discharge) < 1e-8
    @test abs(mb.imbalance) < 1e-8
end
