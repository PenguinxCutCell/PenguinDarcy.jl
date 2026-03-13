"""
Return `(N, T)` for the tracker geometry space dimension and numeric type.
"""
function tracker_dimension_type(::AbstractDarcyTracker)
    throw(MethodError(tracker_dimension_type, (AbstractDarcyTracker,)))
end

"""
Return a copy of the tracker state (graph position array/scalar).
"""
function tracker_state(::AbstractDarcyTracker)
    throw(MethodError(tracker_state, (AbstractDarcyTracker,)))
end

"""
Return interface positions in graph form (`x = h(...)` or analogous axis choice).
"""
interface_positions(tracker::AbstractDarcyTracker) = tracker_state(tracker)

"""
Rebuild a signed-distance-like level-set field for the current tracker state.
"""
function rebuild_signed_distance_or_geometry(::AbstractDarcyTracker)
    throw(MethodError(rebuild_signed_distance_or_geometry, (AbstractDarcyTracker,)))
end

"""
Rebuild one-phase capacities from the current tracker state.
"""
function rebuild_capacities(::AbstractDarcyTracker; kwargs...)
    throw(MethodError(rebuild_capacities, (AbstractDarcyTracker,)))
end

"""
Return interface normals sampled on the tracker graph nodes.
"""
function interface_normals(::AbstractDarcyTracker)
    throw(MethodError(interface_normals, (AbstractDarcyTracker,)))
end

"""
Return curvature sampled on the tracker graph nodes.
"""
function interface_curvature(::AbstractDarcyTracker)
    throw(MethodError(interface_curvature, (AbstractDarcyTracker,)))
end

"""
Update tracker state to a new graph and store optional iteration metadata.
"""
function update_interface!(::AbstractDarcyTracker, state_new; kwargs...)
    throw(MethodError(update_interface!, (AbstractDarcyTracker, typeof(state_new))))
end

"""
Predict next-step interface state from previous tracker history.
"""
function predictor_interface_state(::AbstractDarcyTracker, dt)
    throw(MethodError(predictor_interface_state, (AbstractDarcyTracker, typeof(dt))))
end

"""
Check graph validity constraints for a tracker state.
"""
function check_graph_validity(::AbstractDarcyTracker)
    throw(MethodError(check_graph_validity, (AbstractDarcyTracker,)))
end

"""
Return the running tracker-side iteration history buffer.
"""
function interface_iteration_history(::AbstractDarcyTracker)
    throw(MethodError(interface_iteration_history, (AbstractDarcyTracker,)))
end

tracker_periodic_transverse(::AbstractDarcyTracker) = false

function tracker_transverse_dims(tracker::AbstractDarcyTracker)
    N, _ = tracker_dimension_type(tracker)
    return ntuple(i -> i + 1, N - 1)
end

@inline function _transverse_dims(axis::Int, N::Int)
    return ntuple(i -> i < axis ? i : i + 1, N - 1)
end

@inline function _axis_symbol(axis::Int)
    if axis == 1
        return :x
    elseif axis == 2
        return :y
    elseif axis == 3
        return :z
    end
    return Symbol("axis", axis)
end

@inline function _clamp_or_wrap(v::T, a::T, b::T, periodic::Bool) where {T}
    if periodic
        L = b - a
        L == zero(T) && return a
        return a + mod(v - a, L)
    end
    return clamp(v, a, b)
end

function _nearest_index(nodes::AbstractVector{T}, x::T) where {T}
    n = length(nodes)
    n >= 1 || throw(ArgumentError("nodes must be non-empty"))
    i = searchsortedlast(nodes, x)
    i = clamp(i, 1, n)
    if i == n
        return n
    elseif i == 1
        return abs(nodes[2] - x) < abs(nodes[1] - x) ? 2 : 1
    end
    return abs(nodes[i + 1] - x) < abs(nodes[i] - x) ? (i + 1) : i
end

function _interp_pwlinear(nodes::AbstractVector{T}, vals::AbstractVector{T}, x::T; periodic::Bool=false) where {T}
    length(nodes) == length(vals) || throw(DimensionMismatch("nodes/vals length mismatch"))
    n = length(nodes)
    n == 1 && return vals[1]

    xq = _clamp_or_wrap(x, nodes[1], nodes[end], periodic)
    if xq <= nodes[1]
        return vals[1]
    elseif xq >= nodes[end]
        return vals[end]
    end

    i = searchsortedlast(nodes, xq)
    i = clamp(i, 1, n - 1)
    x0 = nodes[i]
    x1 = nodes[i + 1]
    θ = x1 == x0 ? zero(T) : (xq - x0) / (x1 - x0)
    return (one(T) - θ) * vals[i] + θ * vals[i + 1]
end

function _interp_bilinear(
    xnodes::AbstractVector{T},
    ynodes::AbstractVector{T},
    vals::AbstractMatrix{T},
    x::T,
    y::T;
    periodic_x::Bool=false,
    periodic_y::Bool=false,
) where {T}
    nx = length(xnodes)
    ny = length(ynodes)
    size(vals) == (nx, ny) || throw(DimensionMismatch("vals shape must be ($(nx), $(ny))"))

    xq = _clamp_or_wrap(x, xnodes[1], xnodes[end], periodic_x)
    yq = _clamp_or_wrap(y, ynodes[1], ynodes[end], periodic_y)

    ix = clamp(searchsortedlast(xnodes, xq), 1, nx - 1)
    iy = clamp(searchsortedlast(ynodes, yq), 1, ny - 1)

    x0, x1 = xnodes[ix], xnodes[ix + 1]
    y0, y1 = ynodes[iy], ynodes[iy + 1]
    θx = x1 == x0 ? zero(T) : (xq - x0) / (x1 - x0)
    θy = y1 == y0 ? zero(T) : (yq - y0) / (y1 - y0)

    v00 = vals[ix, iy]
    v10 = vals[ix + 1, iy]
    v01 = vals[ix, iy + 1]
    v11 = vals[ix + 1, iy + 1]

    vx0 = (one(T) - θx) * v00 + θx * v10
    vx1 = (one(T) - θx) * v01 + θx * v11
    return (one(T) - θy) * vx0 + θy * vx1
end

function _fill_missing_by_nearest!(a::AbstractVector{T}, mask::AbstractVector{Bool}) where {T}
    n = length(a)
    n == length(mask) || throw(DimensionMismatch("mask length mismatch"))
    all(mask) && return a

    idx = findall(mask)
    isempty(idx) && throw(ArgumentError("cannot fill missing values: no valid samples"))

    for i in 1:n
        mask[i] && continue
        best = idx[1]
        bestd = abs(i - best)
        for j in idx
            d = abs(i - j)
            if d < bestd
                bestd = d
                best = j
            end
        end
        a[i] = a[best]
    end
    return a
end

function _fill_missing_by_nearest!(a::AbstractMatrix{T}, mask::AbstractMatrix{Bool}) where {T}
    size(a) == size(mask) || throw(DimensionMismatch("mask shape mismatch"))
    all(mask) && return a

    valid = findall(mask)
    isempty(valid) && throw(ArgumentError("cannot fill missing values: no valid samples"))

    for I in CartesianIndices(a)
        mask[I] && continue
        best = valid[1]
        bestd = sum(abs.(Tuple(I) .- Tuple(best)))
        for J in valid
            d = sum(abs.(Tuple(I) .- Tuple(J)))
            if d < bestd
                bestd = d
                best = J
            end
        end
        a[I] = a[best]
    end
    return a
end
