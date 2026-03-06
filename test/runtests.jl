using Test
using LinearAlgebra
using SparseArrays

using CartesianGeometry: geometric_moments, nan
using CartesianOperators
using PenguinBCs
using PenguinDarcy

full_moments(grid) = geometric_moments((args...) -> -1.0, grid, Float64, nan; method=:vofijul)

function physical_indices(nnodes::NTuple{N,Int}) where {N}
    LI = LinearIndices(nnodes)
    idx = Int[]
    for I in CartesianIndices(nnodes)
        if all(d -> I[d] < nnodes[d], 1:N)
            push!(idx, LI[I])
        end
    end
    return idx
end

function active_indices(cap)
    LI = LinearIndices(cap.nnodes)
    idx = Int[]
    for I in CartesianIndices(cap.nnodes)
        i = LI[I]
        halo = any(d -> I[d] == cap.nnodes[d], 1:length(cap.nnodes))
        halo && continue
        v = cap.buf.V[i]
        if isfinite(v) && v > 0.0
            push!(idx, i)
        end
    end
    return idx
end

interface_indices(cap) = findall(i -> isfinite(cap.buf.Γ[i]) && cap.buf.Γ[i] > 0.0, 1:cap.ntotal)

include("test_steady_1d.jl")
include("test_steady_2d.jl")
include("test_unsteady_1d.jl")
include("test_embedded_boundary.jl")
