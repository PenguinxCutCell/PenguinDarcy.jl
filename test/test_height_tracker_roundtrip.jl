@testset "height-function tracker graph/sdf roundtrip" begin
    grid = CartesianGrid((0.0, 0.0), (1.0, 1.0), (65, 49))

    y = collect(grid1d(grid, 2))
    xf = [0.45 + 0.06 * sin(2 * pi * yy) for yy in y]
    xf[end] = xf[1]

    tx = HeightFunctionTracker(grid, xf; axis=:x, periodic_transverse=true)
    ϕx = PenguinDarcy.rebuild_signed_distance_or_geometry(tx)
    xf_rt = GlobalHeightFunctions.xf_from_sdf(ϕx, grid; axis=:x)
    GlobalHeightFunctions.ensure_periodic!(xf_rt)

    @test maximum(abs.(xf_rt .- tx.xf)) < 5e-3
    @test abs(xf_rt[end] - xf_rt[1]) < 1e-12

    x = collect(grid1d(grid, 1))
    yf = [0.55 + 0.04 * cos(2 * pi * xx) for xx in x]
    yf[end] = yf[1]

    ty = HeightFunctionTracker(grid, yf; axis=:y, periodic_transverse=true)
    ϕy = PenguinDarcy.rebuild_signed_distance_or_geometry(ty)
    yf_rt = GlobalHeightFunctions.xf_from_sdf(ϕy, grid; axis=:y)
    GlobalHeightFunctions.ensure_periodic!(yf_rt)

    @test maximum(abs.(yf_rt .- ty.xf)) < 5e-3
    @test abs(yf_rt[end] - yf_rt[1]) < 1e-12
end
