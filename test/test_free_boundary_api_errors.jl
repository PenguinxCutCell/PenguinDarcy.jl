@testset "moving free-boundary API errors" begin
    grid = CartesianGrid((0.0, 0.0), (1.0, 1.0), (21, 17))
    y = collect(grid1d(grid, 2))
    xf = fill(0.4, length(y))
    xf[end] = xf[1]

    @test_throws ArgumentError HeightFunctionTracker(grid, xf; axis=:z)

    xf_bad = copy(xf)
    xf_bad[3] = 1.2
    @test_throws ArgumentError HeightFunctionTracker(grid, xf_bad; axis=:x, periodic_transverse=true)

    tracker = HeightFunctionTracker(grid, xf; axis=:x, periodic_transverse=true)

    bc_nonperiodic = BorderConditions(; left=Neumann(0.0), right=Neumann(0.0))
    @test_throws ArgumentError MovingDarcyModelMono(tracker, 1.0; bc_border=bc_nonperiodic)

    bc_periodic = BorderConditions(
        ;
        left=Neumann(0.0),
        right=Neumann(0.0),
        bottom=Periodic(),
        top=Periodic(),
    )

    @test_throws ArgumentError MovingDarcyModelMono(tracker, 1.0; bc_border=bc_periodic, storage=1.0)

    @test_throws ArgumentError MovingDarcyModelDiph(
        tracker,
        1.0,
        1.0;
        bc_border=bc_periodic,
        bc_interface=DarcyContinuity(0.0, 1e-4),
    )
end
