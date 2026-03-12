using PenguinSolverCore

@testset "DarcyCoupledModelMono SolverCore adapter" begin
    grid = (range(0.0, 1.0; length=65),)
    cap = assembled_capacity(full_moments(grid); bc=0.0)
    bc = BorderConditions(; left=Dirichlet(1.0), right=Dirichlet(0.0))
    ops = DiffusionOps(cap; periodic=periodic_flags(bc, 1))

    λ0 = 2.0
    α = 0.5
    callback_calls = Ref(0)

    mobility_cb = function (base_model, concentration, t)
        callback_calls[] += 1
        lay = base_model.layout.offsets
        cω = concentration[lay.ω]
        cavg = sum(cω) / length(cω)
        return λ0 * (1 + α * cavg)
    end

    base = DarcyModelMono(cap, ops, λ0; source=(x, t) -> 0.0, bc_border=bc)
    coupled_model = DarcyCoupledModelMono(base; mobility_update! = mobility_cb)

    block = CoupledBlock(:darcy, coupled_model; init=nothing, cache=Dict{Symbol,Any}())

    @test block.state.pressure isa Vector{Float64}
    @test block.state.velocity isa NamedTuple

    advance_steady!(block)
    vel0 = get_coupling_field(block, Val(:velocity))

    @test hasproperty(vel0, :x)
    @test length(vel0.x) == cap.ntotal

    cin = fill(0.2, length(block.state.pressure))
    set_coupling_field!(block, Val(:concentration), cin)
    @test get_coupling_field(block, Val(:concentration)) == cin

    advance_steady!(block)
    @test callback_calls[] >= 1

    vel1 = get_coupling_field(block, Val(:velocity))
    λ_expected = λ0 * (1 + α * 0.2)

    LI = LinearIndices(cap.nnodes)
    idx_v = Int[]
    for I in CartesianIndices(cap.nnodes)
        if I[1] > 1 && I[1] < cap.nnodes[1]
            i = LI[I]
            cap.buf.V[i] > 0.0 && push!(idx_v, i)
        end
    end

    @test maximum(abs.(vel1.x[idx_v] .- λ_expected)) < 3e-6

    advance_unsteady!(block, 0.0, 0.1)
    @test block.state.last_time ≈ 0.1
    @test haskey(block.cache, :darcy_system)

    s = PenguinSolverCore.block_summary(block)
    @test occursin("DarcyCoupledModelMono", s)
end
