mutable struct HeightFunctionTracker{N,T,AT} <: AbstractDarcyTracker
    grid::CartesianGrid{N,T}
    axis::Int
    xf::AT
    xf_prev::AT
    periodic_transverse::Bool
    interp::Symbol
    curvature_smoothing::Int
    damping::T
    max_iter::Int
    tol_interface::T
    tol_update::T
    dt_prev::T
    history::Vector{InterfaceIterationHistory{T}}
end

@inline function _tracker_tshape(grid::CartesianGrid{N}, axis::Int) where {N}
    return N == 1 ? () : ntuple(k -> grid.n[k < axis ? k : k + 1], N - 1)
end

function _init_tracker_xf(
    xf0,
    grid::CartesianGrid{N,T},
    axis::Int,
) where {N,T}
    tshape = _tracker_tshape(grid, axis)
    if N == 1
        if xf0 isa Number
            out = Array{T,0}(undef)
            out[] = convert(T, xf0)
            return out
        elseif xf0 isa AbstractArray
            out = Array{T,0}(undef)
            out[] = convert(T, xf0[])
            return out
        end
        throw(ArgumentError("xf0 must be scalar for 1D trackers"))
    end

    if xf0 isa Number
        return fill(convert(T, xf0), tshape...)
    elseif xf0 isa AbstractArray
        size(xf0) == tshape || throw(DimensionMismatch("xf0 shape $(size(xf0)) must match transverse shape $tshape"))
        return convert(Array{T,N-1}, xf0)
    elseif xf0 isa Function
        out = Array{T,N-1}(undef, tshape...)
        tdims = _transverse_dims(axis, N)
        xnodes = ntuple(k -> collect(T, grid1d(grid, tdims[k])), N - 1)
        if N == 2
            for j in eachindex(out)
                s = xnodes[1][j]
                if applicable(xf0, s)
                    out[j] = convert(T, xf0(s))
                elseif applicable(xf0, s, zero(T))
                    out[j] = convert(T, xf0(s, zero(T)))
                else
                    throw(ArgumentError("xf0 function must accept transverse coordinates"))
                end
            end
        elseif N == 3
            for I in CartesianIndices(out)
                y = xnodes[1][I[1]]
                z = xnodes[2][I[2]]
                if applicable(xf0, y, z)
                    out[I] = convert(T, xf0(y, z))
                elseif applicable(xf0, y, z, zero(T))
                    out[I] = convert(T, xf0(y, z, zero(T)))
                else
                    throw(ArgumentError("xf0 function must accept transverse coordinates"))
                end
            end
        else
            throw(ArgumentError("HeightFunctionTracker supports N <= 3"))
        end
        return out
    end

    throw(ArgumentError("unsupported xf0 initializer of type $(typeof(xf0))"))
end

function HeightFunctionTracker(
    grid::CartesianGrid{N,T},
    xf0;
    axis::Union{Int,Symbol}=1,
    periodic_transverse::Bool=false,
    interp::Symbol=:linear,
    curvature_smoothing::Int=0,
    damping::Real=0.5,
    max_iter::Int=20,
    tol_interface::Real=1e-8,
    tol_update::Real=1e-8,
) where {N,T}
    axis_idx = GlobalHeightFunctions.axis_to_index(axis, Val(N))
    xf = _init_tracker_xf(xf0, grid, axis_idx)
    periodic_transverse && GlobalHeightFunctions.ensure_periodic!(xf)
    xf_prev = copy(xf)

    tr = HeightFunctionTracker{N,T,typeof(xf)}(
        grid,
        axis_idx,
        xf,
        xf_prev,
        periodic_transverse,
        interp,
        curvature_smoothing,
        convert(T, damping),
        max_iter,
        convert(T, tol_interface),
        convert(T, tol_update),
        zero(T),
        InterfaceIterationHistory{T}[],
    )
    check_graph_validity(tr)
    return tr
end

tracker_dimension_type(::HeightFunctionTracker{N,T}) where {N,T} = (N, T)
tracker_state(tr::HeightFunctionTracker) = copy(tr.xf)
interface_iteration_history(tr::HeightFunctionTracker) = copy(tr.history)
tracker_periodic_transverse(tr::HeightFunctionTracker) = tr.periodic_transverse
tracker_transverse_dims(tr::HeightFunctionTracker{N}) where {N} = _transverse_dims(tr.axis, N)

function _tracker_eval_graph(
    tr::HeightFunctionTracker{N,T},
    xf,
    x::NTuple{N,T},
) where {N,T}
    ia = tr.axis
    if N == 1
        return x[1] - xf[]
    elseif N == 2
        td = _transverse_dims(ia, N)[1]
        snodes = collect(T, grid1d(tr.grid, td))
        h = _interp_pwlinear(snodes, vec(xf), x[td]; periodic=tr.periodic_transverse)
        return x[ia] - h
    elseif N == 3
        td = _transverse_dims(ia, N)
        ynodes = collect(T, grid1d(tr.grid, td[1]))
        znodes = collect(T, grid1d(tr.grid, td[2]))
        h = _interp_bilinear(
            ynodes,
            znodes,
            xf,
            x[td[1]],
            x[td[2]];
            periodic_x=tr.periodic_transverse,
            periodic_y=tr.periodic_transverse,
        )
        return x[ia] - h
    end
    throw(ArgumentError("HeightFunctionTracker supports N <= 3"))
end

function rebuild_signed_distance_or_geometry(tr::HeightFunctionTracker{N,T}) where {N,T}
    return GlobalHeightFunctions.phi_from_xf(
        tr.xf,
        tr.grid;
        axis=tr.axis,
        interp=tr.interp,
        periodic=tr.periodic_transverse,
    )
end

function rebuild_capacities(
    tr::HeightFunctionTracker{N,T};
    geom_method::Symbol=:vofijul,
    bc::T=zero(T),
    two_phase::Bool=false,
) where {N,T}
    check_graph_validity(tr)
    xyz = ntuple(d -> grid1d(tr.grid, d), N)

    f = (args...) -> begin
        x = ntuple(d -> convert(T, args[d]), N)
        _tracker_eval_graph(tr, tr.xf, x)
    end

    moms1 = geometric_moments(f, xyz, T, nan; method=geom_method)
    cap1 = assembled_capacity(moms1; bc=bc)
    if !two_phase
        return cap1
    end

    f2 = (args...) -> begin
        x = ntuple(d -> convert(T, args[d]), N)
        -_tracker_eval_graph(tr, tr.xf, x)
    end
    moms2 = geometric_moments(f2, xyz, T, nan; method=geom_method)
    cap2 = assembled_capacity(moms2; bc=bc)
    return cap1, cap2
end

function predictor_interface_state(tr::HeightFunctionTracker{N,T}, dt::Real) where {N,T}
    dtt = convert(T, dt)
    dtt > zero(T) || throw(ArgumentError("dt must be positive"))

    xpred = copy(tr.xf)
    if tr.dt_prev > zero(T)
        ratio = dtt / tr.dt_prev
        xpred .= tr.xf .+ ratio .* (tr.xf .- tr.xf_prev)
    end
    tr.periodic_transverse && GlobalHeightFunctions.ensure_periodic!(xpred)
    return xpred
end

function check_graph_validity(tr::HeightFunctionTracker{N,T}) where {N,T}
    xlo = tr.grid.lc[tr.axis]
    xhi = tr.grid.hc[tr.axis]
    for v in tr.xf
        isfinite(v) || throw(ArgumentError("height-function tracker contains non-finite values"))
        (xlo <= v <= xhi) || throw(ArgumentError("height-function tracker leaves box bounds along axis $(_axis_symbol(tr.axis))"))
    end

    if tr.periodic_transverse
        # Keep periodic copies synchronized in the same convention as GHF.
        GlobalHeightFunctions.ensure_periodic!(tr.xf)
    end
    return true
end

function update_interface!(
    tr::HeightFunctionTracker{N,T},
    state_new;
    dt::Union{Nothing,Real}=nothing,
    iteration::Union{Nothing,InterfaceIterationHistory{T}}=nothing,
) where {N,T}
    size(state_new) == size(tr.xf) || throw(DimensionMismatch("new tracker state shape mismatch"))
    tr.xf_prev .= tr.xf
    tr.xf .= state_new
    tr.periodic_transverse && GlobalHeightFunctions.ensure_periodic!(tr.xf)
    check_graph_validity(tr)
    if dt !== nothing
        tr.dt_prev = convert(T, dt)
    end
    if iteration !== nothing
        push!(tr.history, iteration)
    end
    return tr
end

function _first_derivative_1d(vals::AbstractVector{T}, nodes::AbstractVector{T}; periodic::Bool=false) where {T}
    n = length(vals)
    n == length(nodes) || throw(DimensionMismatch("nodes/vals length mismatch"))
    out = zeros(T, n)
    n <= 1 && return out
    Δ = nodes[2] - nodes[1]

    if periodic
        nphys = n > 1 ? n - 1 : n
        for i in 1:nphys
            im = i == 1 ? nphys : i - 1
            ip = i == nphys ? 1 : i + 1
            out[i] = (vals[ip] - vals[im]) / (2 * Δ)
        end
        out[end] = out[1]
        return out
    end

    out[1] = (vals[2] - vals[1]) / Δ
    @inbounds for i in 2:(n - 1)
        out[i] = (vals[i + 1] - vals[i - 1]) / (2 * Δ)
    end
    out[n] = (vals[n] - vals[n - 1]) / Δ
    return out
end

function _second_derivative_1d(vals::AbstractVector{T}, nodes::AbstractVector{T}; periodic::Bool=false) where {T}
    n = length(vals)
    n == length(nodes) || throw(DimensionMismatch("nodes/vals length mismatch"))
    out = zeros(T, n)
    n <= 2 && return out
    Δ = nodes[2] - nodes[1]
    Δ2 = Δ * Δ

    if periodic
        nphys = n > 1 ? n - 1 : n
        for i in 1:nphys
            im = i == 1 ? nphys : i - 1
            ip = i == nphys ? 1 : i + 1
            out[i] = (vals[ip] - 2 * vals[i] + vals[im]) / Δ2
        end
        out[end] = out[1]
        return out
    end

    out[1] = (vals[3] - 2 * vals[2] + vals[1]) / Δ2
    @inbounds for i in 2:(n - 1)
        out[i] = (vals[i + 1] - 2 * vals[i] + vals[i - 1]) / Δ2
    end
    out[n] = (vals[n] - 2 * vals[n - 1] + vals[n - 2]) / Δ2
    return out
end

function _derivative_dim1(vals::AbstractMatrix{T}, nodes::AbstractVector{T}; periodic::Bool=false) where {T}
    ny, nz = size(vals)
    out = similar(vals)
    for j in 1:nz
        out[:, j] .= _first_derivative_1d(view(vals, :, j), nodes; periodic=periodic)
    end
    return out
end

function _derivative_dim2(vals::AbstractMatrix{T}, nodes::AbstractVector{T}; periodic::Bool=false) where {T}
    ny, nz = size(vals)
    out = similar(vals)
    for i in 1:ny
        out[i, :] .= _first_derivative_1d(view(vals, i, :), nodes; periodic=periodic)
    end
    return out
end

function _second_derivative_dim1(vals::AbstractMatrix{T}, nodes::AbstractVector{T}; periodic::Bool=false) where {T}
    ny, nz = size(vals)
    out = similar(vals)
    for j in 1:nz
        out[:, j] .= _second_derivative_1d(view(vals, :, j), nodes; periodic=periodic)
    end
    return out
end

function _second_derivative_dim2(vals::AbstractMatrix{T}, nodes::AbstractVector{T}; periodic::Bool=false) where {T}
    ny, nz = size(vals)
    out = similar(vals)
    for i in 1:ny
        out[i, :] .= _second_derivative_1d(view(vals, i, :), nodes; periodic=periodic)
    end
    return out
end

function _smooth_curvature!(κ::AbstractVector{T}, niter::Int, periodic::Bool) where {T}
    niter <= 0 && return κ
    buf = similar(κ)
    n = length(κ)
    for _ in 1:niter
        if periodic
            nphys = n > 1 ? n - 1 : n
            for i in 1:nphys
                im = i == 1 ? nphys : i - 1
                ip = i == nphys ? 1 : i + 1
                buf[i] = (κ[im] + 2 * κ[i] + κ[ip]) / 4
            end
            buf[end] = buf[1]
        else
            buf[1] = κ[1]
            for i in 2:(n - 1)
                buf[i] = (κ[i - 1] + 2 * κ[i] + κ[i + 1]) / 4
            end
            buf[end] = κ[end]
        end
        κ .= buf
    end
    return κ
end

function _smooth_curvature!(κ::AbstractMatrix{T}, niter::Int, periodic::Bool) where {T}
    niter <= 0 && return κ
    ny, nz = size(κ)
    buf = similar(κ)
    for _ in 1:niter
        for I in CartesianIndices(κ)
            i, j = Tuple(I)
            im = periodic ? (i == 1 ? ny - 1 : i - 1) : max(i - 1, 1)
            ip = periodic ? (i == ny ? 2 : i + 1) : min(i + 1, ny)
            jm = periodic ? (j == 1 ? nz - 1 : j - 1) : max(j - 1, 1)
            jp = periodic ? (j == nz ? 2 : j + 1) : min(j + 1, nz)
            if periodic && (i == ny || j == nz)
                continue
            end
            buf[I] = (κ[i, j] * 4 + κ[im, j] + κ[ip, j] + κ[i, jm] + κ[i, jp]) / 8
        end
        if periodic
            ny > 1 && (buf[end, :] .= buf[1, :])
            nz > 1 && (buf[:, end] .= buf[:, 1])
        end
        κ .= buf
    end
    return κ
end

function _graph_first_derivatives(tr::HeightFunctionTracker{2,T}) where {T}
    td = _transverse_dims(tr.axis, 2)[1]
    snodes = collect(T, grid1d(tr.grid, td))
    h1 = _first_derivative_1d(vec(tr.xf), snodes; periodic=tr.periodic_transverse)
    return (h1,)
end

function _graph_first_derivatives(tr::HeightFunctionTracker{3,T}) where {T}
    td = _transverse_dims(tr.axis, 3)
    ynodes = collect(T, grid1d(tr.grid, td[1]))
    znodes = collect(T, grid1d(tr.grid, td[2]))
    hy = _derivative_dim1(tr.xf, ynodes; periodic=tr.periodic_transverse)
    hz = _derivative_dim2(tr.xf, znodes; periodic=tr.periodic_transverse)
    return hy, hz
end

function interface_normals(tr::HeightFunctionTracker{N,T}) where {N,T}
    if N == 1
        n = Array{T,0}(undef)
        n[] = one(T)
        return (axis_component=n, components=(n,))
    elseif N == 2
        (h1,) = _graph_first_derivatives(tr)
        den = sqrt.(one(T) .+ h1 .^ 2)
        na = one(T) ./ den

        comps = ntuple(d -> begin
            if d == tr.axis
                return na
            else
                return -h1 ./ den
            end
        end, 2)
        return (axis_component=na, components=comps)
    elseif N == 3
        hy, hz = _graph_first_derivatives(tr)
        den = sqrt.(one(T) .+ hy .^ 2 .+ hz .^ 2)
        na = one(T) ./ den
        td = _transverse_dims(tr.axis, 3)
        comps = ntuple(d -> begin
            if d == tr.axis
                return na
            elseif d == td[1]
                return -hy ./ den
            else
                return -hz ./ den
            end
        end, 3)
        return (axis_component=na, components=comps)
    end

    throw(ArgumentError("HeightFunctionTracker supports N <= 3"))
end

function interface_curvature(tr::HeightFunctionTracker{N,T}) where {N,T}
    if N == 1
        κ = Array{T,0}(undef)
        κ[] = zero(T)
        return κ
    elseif N == 2
        td = _transverse_dims(tr.axis, 2)[1]
        snodes = collect(T, grid1d(tr.grid, td))
        h = vec(tr.xf)
        h1 = _first_derivative_1d(h, snodes; periodic=tr.periodic_transverse)
        h2 = _second_derivative_1d(h, snodes; periodic=tr.periodic_transverse)
        κ = -h2 ./ (one(T) .+ h1 .^ 2) .^ (T(3) / T(2))
        _smooth_curvature!(κ, tr.curvature_smoothing, tr.periodic_transverse)
        return κ
    elseif N == 3
        td = _transverse_dims(tr.axis, 3)
        ynodes = collect(T, grid1d(tr.grid, td[1]))
        znodes = collect(T, grid1d(tr.grid, td[2]))

        h = tr.xf
        hy = _derivative_dim1(h, ynodes; periodic=tr.periodic_transverse)
        hz = _derivative_dim2(h, znodes; periodic=tr.periodic_transverse)
        hyy = _second_derivative_dim1(h, ynodes; periodic=tr.periodic_transverse)
        hzz = _second_derivative_dim2(h, znodes; periodic=tr.periodic_transverse)
        hyz = _derivative_dim1(hz, ynodes; periodic=tr.periodic_transverse)

        den = (one(T) .+ hy .^ 2 .+ hz .^ 2) .^ (T(3) / T(2))
        num = (one(T) .+ hz .^ 2) .* hyy .- 2 .* hy .* hz .* hyz .+ (one(T) .+ hy .^ 2) .* hzz
        κ = -num ./ den
        _smooth_curvature!(κ, tr.curvature_smoothing, tr.periodic_transverse)
        return κ
    end

    throw(ArgumentError("HeightFunctionTracker supports N <= 3"))
end
