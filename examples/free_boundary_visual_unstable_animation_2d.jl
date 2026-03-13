using PenguinDarcy
using CartesianGrids
using PenguinBCs
using CairoMakie

function solve_unstable_case(; tend=0.10, dt=0.01)
    grid = CartesianGrid((0.0, 0.0), (1.0, 1.0), (81, 65))
    y = collect(grid1d(grid, 2))

    h0 = 0.45
    a0 = 0.025
    xf0 = [h0 + a0 * sin(2 * pi * yy) for yy in y]
    xf0[end] = xf0[1]

    tracker = HeightFunctionTracker(
        grid,
        xf0;
        axis=:x,
        periodic_transverse=true,
        damping=0.6,
        max_iter=25,
        tol_interface=1e-9,
        tol_update=1e-9,
    )

    bc = BorderConditions(
        ;
        left=Dirichlet(0.0),
        right=Dirichlet(0.0),
        bottom=Periodic(),
        top=Periodic(),
    )

    # Unstable Muskat-like configuration for axis=:x and gravity along +x.
    model = MovingDarcyModelDiph(
        tracker,
        1.0,
        1.0;
        source=(x, y, t) -> (0.0, 0.0),
        bc_border=bc,
        bc_interface=DarcyContinuity(),
        ρ=(2.0, 1.0),
        gravity=(1.0, 0.0),
    )

    return solve_unsteady_moving!(model, (0.0, tend); dt=dt, save_history=true), y
end

left_polygon(xf::AbstractVector, ys::AbstractVector) = begin
    pts = Point2f[Point2f(0.0f0, Float32(ys[1]))]
    append!(pts, [Point2f(Float32(xf[j]), Float32(ys[j])) for j in eachindex(ys)])
    push!(pts, Point2f(0.0f0, Float32(ys[end])))
    return pts
end

right_polygon(xf::AbstractVector, ys::AbstractVector) = begin
    pts = Point2f[Point2f(1.0f0, Float32(ys[1])), Point2f(1.0f0, Float32(ys[end]))]
    append!(pts, [Point2f(Float32(xf[j]), Float32(ys[j])) for j in reverse(eachindex(ys))])
    return pts
end

function make_visuals(sol, y; outdir=joinpath(@__DIR__, "output", "free_boundary_visual"))
    mkpath(outdir)

    ys = y[1:(end - 1)]
    interfaces = [xf[1:(end - 1)] for xf in sol.interface_states]
    times = sol.times
    amps = [(maximum(xf) - minimum(xf)) / 2 for xf in interfaces]
    nframes = length(interfaces)

    frame_idx = Observable(1)

    left_obs = lift(i -> left_polygon(interfaces[i], ys), frame_idx)
    right_obs = lift(i -> right_polygon(interfaces[i], ys), frame_idx)
    x_obs = lift(i -> interfaces[i], frame_idx)
    t_obs = lift(i -> times[i], frame_idx)
    a_obs = lift(i -> amps[i], frame_idx)

    fig = Figure(size=(1400, 900), backgroundcolor=RGBf(0.04, 0.06, 0.10))

    ax = Axis(
        fig[1, 1],
        title=lift((tt, aa) -> "Unstable Free-Boundary Darcy (Muskat-like)   t=$(round(tt, digits=3))   amp=$(round(aa, digits=4))", t_obs, a_obs),
        xlabel="x",
        ylabel="y",
        limits=(0.0, 1.0, 0.0, 1.0),
        backgroundcolor=RGBf(0.08, 0.11, 0.16),
    )

    poly!(ax, left_obs, color=RGBAf(0.12, 0.42, 0.90, 0.92), strokecolor=:transparent)
    poly!(ax, right_obs, color=RGBAf(0.95, 0.48, 0.16, 0.88), strokecolor=:transparent)
    lines!(ax, x_obs, ys, color=RGBAf(0.97, 0.98, 1.0, 0.98), linewidth=3)

    ax2 = Axis(
        fig[2, 1],
        xlabel="time",
        ylabel="amplitude",
        title="Interface amplitude history",
        backgroundcolor=RGBf(0.08, 0.11, 0.16),
    )
    lines!(ax2, times, amps, color=RGBAf(0.9, 0.9, 0.25, 1.0), linewidth=3)
    scatter!(ax2, lift(i -> [times[i]], frame_idx), lift(i -> [amps[i]], frame_idx), color=RGBAf(1.0, 0.2, 0.2, 1.0), markersize=12)

    still_ids = unique(round.(Int, LinRange(1, nframes, min(nframes, 8))))
    still_paths = String[]
    for i in still_ids
        frame_idx[] = i
        out = joinpath(outdir, "frame_$(lpad(i, 3, '0')).png")
        save(out, fig)
        push!(still_paths, out)
    end

    gif_path = joinpath(outdir, "free_boundary_unstable.gif")
    record(fig, gif_path, 1:nframes; framerate=12) do i
        frame_idx[] = i
    end

    return still_paths, gif_path
end

sol, y = solve_unstable_case()
stills, gif_path = make_visuals(sol, y)

println("Saved still images:")
for p in stills
    println(" - ", p)
end
println("Saved animation:")
println(" - ", gif_path)
