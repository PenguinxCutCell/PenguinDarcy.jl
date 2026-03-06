abstract type AbstractWell end

struct PointWell{N,T,R} <: AbstractWell
    x0::NTuple{N,T}
    rate::R
    radius::T
    phase::Int
end

function PointWell(x0::NTuple{N,T}, rate; radius::Real=zero(T), phase::Int=1) where {N,T<:Real}
    return PointWell{N,T,typeof(rate)}(x0, rate, convert(T, radius), phase)
end

struct CellWell{M,R} <: AbstractWell
    mask_or_indices::M
    rate::R
    phase::Int
end

CellWell(mask_or_indices, rate; phase::Int=1) = CellWell{typeof(mask_or_indices),typeof(rate)}(mask_or_indices, rate, phase)

function _eval_rate(rate, x::SVector{N,T}, t::T) where {N,T}
    if rate isa Number
        return convert(T, rate)
    elseif rate isa Function
        if applicable(rate, x..., t)
            return convert(T, rate(x..., t))
        elseif applicable(rate, x...,)
            return convert(T, rate(x...))
        elseif applicable(rate, t)
            return convert(T, rate(t))
        elseif applicable(rate)
            return convert(T, rate())
        end
    end
    throw(ArgumentError("well rate must be numeric or callable with (x..., t), (x...), (t), or ()"))
end

function _active_cell_indices(cap::AssembledCapacity{N,T}) where {N,T}
    LI = LinearIndices(cap.nnodes)
    idx = Int[]
    for I in CartesianIndices(cap.nnodes)
        i = LI[I]
        halo = any(d -> I[d] == cap.nnodes[d], 1:N)
        halo && continue
        Vi = cap.buf.V[i]
        if isfinite(Vi) && Vi > zero(T)
            push!(idx, i)
        end
    end
    return idx
end

function assemble_sources!(
    rhs::AbstractVector{T},
    model::DarcyModelMono{N,T},
    cap::AssembledCapacity{N,T};
    t::T,
    phase::Int=1,
) where {N,T}
    length(rhs) == cap.ntotal || throw(DimensionMismatch("rhs length must be cap.ntotal"))
    sval = _source_values_mono(cap, model.source, t)
    @inbounds for i in eachindex(rhs)
        Vi = cap.buf.V[i]
        if isfinite(Vi) && Vi > zero(T)
            rhs[i] += Vi * sval[i]
        end
    end
    return rhs
end

function assemble_sources!(
    rhs::AbstractVector{T},
    model::DarcyModelDiph{N,T},
    cap::AssembledCapacity{N,T};
    t::T,
    phase::Int,
) where {N,T}
    length(rhs) == cap.ntotal || throw(DimensionMismatch("rhs length must be cap.ntotal"))
    source = phase == 1 ? model.source1 : model.source2
    sval = _source_values_mono(cap, source, t)
    @inbounds for i in eachindex(rhs)
        Vi = cap.buf.V[i]
        if isfinite(Vi) && Vi > zero(T)
            rhs[i] += Vi * sval[i]
        end
    end
    return rhs
end

function _cellwell_indices(w::CellWell, cap::AssembledCapacity{N,T}, t::T) where {N,T}
    active = _active_cell_indices(cap)
    if w.mask_or_indices isa AbstractVector{Bool}
        m = w.mask_or_indices
        length(m) == cap.ntotal || throw(DimensionMismatch("CellWell boolean mask length must be cap.ntotal"))
        return [i for i in active if m[i]]
    elseif w.mask_or_indices isa AbstractVector{<:Integer}
        return [i for i in w.mask_or_indices if (1 <= i <= cap.ntotal && i in active)]
    elseif w.mask_or_indices isa Function
        return [i for i in active if begin
            x = cap.C_ω[i]
            if applicable(w.mask_or_indices, x..., t)
                Bool(w.mask_or_indices(x..., t))
            elseif applicable(w.mask_or_indices, x...)
                Bool(w.mask_or_indices(x...))
            else
                throw(ArgumentError("CellWell selector function must accept (x...) or (x..., t)"))
            end
        end]
    end
    throw(ArgumentError("unsupported CellWell mask_or_indices type $(typeof(w.mask_or_indices))"))
end

function assemble_wells!(rhs::AbstractVector{T}, wells::AbstractVector, cap::AssembledCapacity{N,T}; t::T, phase::Int=1) where {N,T}
    length(rhs) == cap.ntotal || throw(DimensionMismatch("rhs length must be cap.ntotal"))
    active = _active_cell_indices(cap)

    for w in wells
        w isa AbstractWell || throw(ArgumentError("well collection must contain AbstractWell elements, got $(typeof(w))"))
        getfield(w, :phase) == phase || continue
        if w isa PointWell
            x0 = SVector{N,T}(w.x0)
            Q = _eval_rate(w.rate, x0, t)
            candidates = Int[]
            weights = T[]

            if w.radius > zero(T)
                for i in active
                    r = norm(cap.C_ω[i] - x0)
                    if r <= w.radius
                        push!(candidates, i)
                        push!(weights, max(w.radius - r, sqrt(eps(T))))
                    end
                end
            end

            if isempty(candidates)
                # Fallback: nearest active cell, still conservative.
                dist_min = T(Inf)
                i_min = first(active)
                for i in active
                    r = norm(cap.C_ω[i] - x0)
                    if r < dist_min
                        dist_min = r
                        i_min = i
                    end
                end
                candidates = [i_min]
                weights = [one(T)]
            end

            sw = sum(weights)
            sw > zero(T) || throw(ArgumentError("invalid PointWell weights (sum <= 0)"))
            for (i, wgt) in zip(candidates, weights)
                rhs[i] += Q * (wgt / sw)
            end
        elseif w isa CellWell
            idx = _cellwell_indices(w, cap, t)
            isempty(idx) && continue
            xref = cap.C_ω[idx[1]]
            Q = _eval_rate(w.rate, xref, t)
            qshare = Q / length(idx)
            for i in idx
                rhs[i] += qshare
            end
        else
            throw(ArgumentError("unsupported well type $(typeof(w))"))
        end
    end

    return rhs
end

function integrated_well_rate(wells::AbstractVector, cap::AssembledCapacity{N,T}; t::T, phase::Int=1) where {N,T}
    rhs = zeros(T, cap.ntotal)
    assemble_wells!(rhs, wells, cap; t=t, phase=phase)
    return sum(rhs)
end

integrated_well_rate(model::DarcyModelMono{N,T}; t::Real=zero(T)) where {N,T} =
    integrated_well_rate(model.wells, model.cap; t=convert(T, t), phase=1)

function integrated_well_rate(model::DarcyModelDiph{N,T}; t::Real=zero(T)) where {N,T}
    tt = convert(T, t)
    q1 = integrated_well_rate(model.wells1, model.cap1; t=tt, phase=1)
    q2 = integrated_well_rate(model.wells2, model.cap2; t=tt, phase=2)
    return (phase1=q1, phase2=q2, total=q1 + q2)
end

function integrated_source(model::DarcyModelMono{N,T}; t::Real=zero(T)) where {N,T}
    tt = convert(T, t)
    rhs = zeros(T, model.cap.ntotal)
    assemble_sources!(rhs, model, model.cap; t=tt, phase=1)
    return sum(rhs)
end

function integrated_source(model::DarcyModelDiph{N,T}; t::Real=zero(T)) where {N,T}
    tt = convert(T, t)
    rhs1 = zeros(T, model.cap1.ntotal)
    rhs2 = zeros(T, model.cap2.ntotal)
    assemble_sources!(rhs1, model, model.cap1; t=tt, phase=1)
    assemble_sources!(rhs2, model, model.cap2; t=tt, phase=2)
    s1 = sum(rhs1)
    s2 = sum(rhs2)
    return (phase1=s1, phase2=s2, total=s1 + s2)
end
