@inline _inf_norm(v::Number) = abs(v)
@inline _inf_norm(v) = maximum(abs, v)

function interface_normal_velocity(
    model::DarcyModelMono{N,T},
    state;
    t::Real=zero(T),
) where {N,T}
    tt = convert(T, t)
    flux = recover_flux(model, state; t=tt)
    qn = model.ops.H' * flux.faces
    vn = zeros(T, length(qn))
    @inbounds for i in eachindex(qn)
        Γi = model.cap.buf.Γ[i]
        if isfinite(Γi) && Γi > zero(T)
            vn[i] = qn[i] / Γi
        end
    end
    return vn
end

function interface_normal_velocity(
    model::DarcyModelDiph{N,T},
    state;
    t::Real=zero(T),
) where {N,T}
    tt = convert(T, t)
    flux = recover_flux(model, state; t=tt)
    q1 = model.ops1.H' * flux.phase1.faces
    q2 = model.ops2.H' * flux.phase2.faces
    vn1 = zeros(T, length(q1))
    vn2 = zeros(T, length(q2))
    mismatch = zeros(T, length(q1))

    @inbounds for i in eachindex(q1)
        Γi = model.cap1.buf.Γ[i]
        if isfinite(Γi) && Γi > zero(T)
            vn1[i] = q1[i] / Γi
            vn2[i] = q2[i] / Γi
            # Same normal orientation (phase-1 normal): u1·n - u2·n = u1·n + u2·n2.
            mismatch[i] = vn1[i] + vn2[i]
        end
    end

    return (phase1=vn1, phase2=vn2, mismatch=mismatch)
end

function interface_normal_velocity(
    tracker::HeightFunctionTracker{N,T},
    cap::AssembledCapacity{N,T},
    qn_cell::AbstractVector{T},
) where {N,T}
    length(qn_cell) == cap.ntotal || throw(DimensionMismatch("qn_cell length must match cap.ntotal"))
    ia = tracker.axis

    if N == 1
        out = Array{T,0}(undef)
        num = zero(T)
        den = zero(T)
        @inbounds for i in 1:cap.ntotal
            Γi = cap.buf.Γ[i]
            if isfinite(Γi) && Γi > zero(T)
                num += qn_cell[i]
                den += Γi
            end
        end
        den > zero(T) || throw(ArgumentError("no active interface cells were found for the current tracker state"))
        out[] = num / den
        return out
    elseif N == 2
        td = _transverse_dims(ia, N)[1]
        snodes = collect(T, grid1d(tracker.grid, td))
        n = length(snodes)
        num = zeros(T, n)
        den = zeros(T, n)
        mask = falses(n)

        @inbounds for i in 1:cap.ntotal
            Γi = cap.buf.Γ[i]
            if !(isfinite(Γi) && Γi > zero(T))
                continue
            end
            s = cap.C_γ[i][td]
            j = _nearest_index(snodes, s)
            num[j] += qn_cell[i]
            den[j] += Γi
            mask[j] = true
        end

        any(mask) || throw(ArgumentError("no active interface cells were found for the current tracker state"))

        vn = similar(num)
        @inbounds for j in 1:n
            if mask[j]
                vn[j] = num[j] / den[j]
            else
                vn[j] = zero(T)
            end
        end
        _fill_missing_by_nearest!(vn, mask)
        tracker.periodic_transverse && GlobalHeightFunctions.ensure_periodic!(vn)
        return vn
    elseif N == 3
        td = _transverse_dims(ia, N)
        ynodes = collect(T, grid1d(tracker.grid, td[1]))
        znodes = collect(T, grid1d(tracker.grid, td[2]))

        ny = length(ynodes)
        nz = length(znodes)
        num = zeros(T, ny, nz)
        den = zeros(T, ny, nz)
        mask = falses(ny, nz)

        @inbounds for i in 1:cap.ntotal
            Γi = cap.buf.Γ[i]
            if !(isfinite(Γi) && Γi > zero(T))
                continue
            end
            y = cap.C_γ[i][td[1]]
            z = cap.C_γ[i][td[2]]
            iy = _nearest_index(ynodes, y)
            iz = _nearest_index(znodes, z)
            num[iy, iz] += qn_cell[i]
            den[iy, iz] += Γi
            mask[iy, iz] = true
        end

        any(mask) || throw(ArgumentError("no active interface cells were found for the current tracker state"))

        vn = similar(num)
        @inbounds for I in CartesianIndices(vn)
            if mask[I]
                vn[I] = num[I] / den[I]
            else
                vn[I] = zero(T)
            end
        end
        _fill_missing_by_nearest!(vn, mask)
        tracker.periodic_transverse && GlobalHeightFunctions.ensure_periodic!(vn)
        return vn
    end

    throw(ArgumentError("HeightFunctionTracker supports N <= 3"))
end

function _interface_porosity_values(
    interface_porosity,
    tracker::HeightFunctionTracker{N,T},
    t::T,
) where {N,T}
    shp = size(tracker.xf)
    if interface_porosity isa Number
        return fill(convert(T, interface_porosity), shp...)
    elseif interface_porosity isa Function
        ia = tracker.axis
        td = _transverse_dims(ia, N)
        vals = similar(tracker.xf)
        if N == 1
            x = tracker.xf[]
            if applicable(interface_porosity, x, t)
                vals[] = convert(T, interface_porosity(x, t))
            elseif applicable(interface_porosity, x)
                vals[] = convert(T, interface_porosity(x))
            else
                throw(ArgumentError("interface_porosity callback must accept interface coordinates"))
            end
            return vals
        elseif N == 2
            snodes = collect(T, grid1d(tracker.grid, td[1]))
            for j in eachindex(vals)
                xa = tracker.xf[j]
                coords = ia == 1 ? (xa, snodes[j]) : (snodes[j], xa)
                if applicable(interface_porosity, coords..., t)
                    vals[j] = convert(T, interface_porosity(coords..., t))
                elseif applicable(interface_porosity, coords...)
                    vals[j] = convert(T, interface_porosity(coords...))
                else
                    throw(ArgumentError("interface_porosity callback must accept interface coordinates"))
                end
            end
            return vals
        elseif N == 3
            ynodes = collect(T, grid1d(tracker.grid, td[1]))
            znodes = collect(T, grid1d(tracker.grid, td[2]))
            for I in CartesianIndices(vals)
                xa = tracker.xf[I]
                yy = ynodes[I[1]]
                zz = znodes[I[2]]
                coords = if ia == 1
                    (xa, yy, zz)
                elseif ia == 2
                    (yy, xa, zz)
                else
                    (yy, zz, xa)
                end
                if applicable(interface_porosity, coords..., t)
                    vals[I] = convert(T, interface_porosity(coords..., t))
                elseif applicable(interface_porosity, coords...)
                    vals[I] = convert(T, interface_porosity(coords...))
                else
                    throw(ArgumentError("interface_porosity callback must accept interface coordinates"))
                end
            end
            return vals
        end
    end

    throw(ArgumentError("interface_porosity must be a positive scalar or callback"))
end

function _graph_rhs_from_normal_velocity(
    tracker::HeightFunctionTracker{N,T},
    vn,
    interface_porosity,
    t::T,
) where {N,T}
    normals = interface_normals(tracker)
    na = normals.axis_component
    ϕ = _interface_porosity_values(interface_porosity, tracker, t)

    out = similar(tracker.xf)
    @inbounds for i in eachindex(out)
        ϕi = ϕ[i]
        ϕi > zero(T) || throw(ArgumentError("interface_porosity must remain positive on all interface points"))
        nai = na[i]
        abs(nai) > sqrt(eps(T)) || throw(ArgumentError("graph kinematics failure: n·e_axis is too small"))
        out[i] = vn[i] / (ϕi * nai)
    end
    return out
end

function _interface_residual(x_guess, x_old, dt, rhs)
    return x_guess .- x_old .- dt .* rhs
end

function _clamp_tracker_state!(xnew, tracker::HeightFunctionTracker{N,T}) where {N,T}
    xlo = tracker.grid.lc[tracker.axis]
    xhi = tracker.grid.hc[tracker.axis]
    @inbounds for i in eachindex(xnew)
        xnew[i] = clamp(xnew[i], xlo, xhi)
    end
    tracker.periodic_transverse && GlobalHeightFunctions.ensure_periodic!(xnew)
    return xnew
end
