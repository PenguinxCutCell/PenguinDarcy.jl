function _eval_density(ρ, x::SVector{N,T}, t::T) where {N,T}
    if ρ isa Number
        return convert(T, ρ)
    elseif ρ isa Function
        if applicable(ρ, x..., t)
            return convert(T, ρ(x..., t))
        elseif applicable(ρ, x...)
            return convert(T, ρ(x...))
        elseif applicable(ρ, t)
            return convert(T, ρ(t))
        elseif applicable(ρ)
            return convert(T, ρ())
        end
    end
    throw(ArgumentError("density ρ must be numeric or callable with (x..., t), (x...), (t), or ()"))
end

function _eval_gravity(gravity, x::SVector{N,T}, t::T) where {N,T}
    if gravity isa NTuple{N}
        return SVector{N,T}(ntuple(d -> _eval_fun_or_const(gravity[d], x, t), N))
    elseif gravity isa AbstractVector
        length(gravity) == N || throw(DimensionMismatch("gravity vector length must be $N"))
        return SVector{N,T}(ntuple(d -> convert(T, gravity[d]), N))
    elseif gravity isa Function
        g = if applicable(gravity, x..., t)
            gravity(x..., t)
        elseif applicable(gravity, x...)
            gravity(x...)
        elseif applicable(gravity, t)
            gravity(t)
        elseif applicable(gravity)
            gravity()
        else
            throw(ArgumentError("gravity callback must accept (x..., t), (x...), (t), or ()"))
        end
        length(g) == N || throw(DimensionMismatch("gravity callback must return length-$N vector/tuple"))
        return SVector{N,T}(ntuple(d -> convert(T, g[d]), N))
    end
    throw(ArgumentError("gravity must be NTuple, vector, or callable"))
end

function gravity_vector(model::DarcyModelMono{N,T}, t::Real=zero(T)) where {N,T}
    tt = convert(T, t)
    xref = model.cap.C_ω[1]
    if model.variable === :head
        return SVector{N,T}(ntuple(d -> _eval_fun_or_const(model.gravity[d], xref, tt), N))
    end
    ρv = _eval_density(model.ρ, xref, tt)
    gv = _eval_gravity(model.gravity, xref, tt)
    return ρv .* gv
end

function gravity_vector(model::DarcyModelDiph{N,T}, t::Real=zero(T); phase::Int=1) where {N,T}
    tt = convert(T, t)
    cap = phase == 1 ? model.cap1 : model.cap2
    ρ = phase == 1 ? model.ρ1 : model.ρ2
    xref = cap.C_ω[1]
    if model.variable === :head
        return SVector{N,T}(ntuple(d -> _eval_fun_or_const(model.gravity[d], xref, tt), N))
    end
    ρv = _eval_density(ρ, xref, tt)
    gv = _eval_gravity(model.gravity, xref, tt)
    return ρv .* gv
end

function face_bodyforce_values(
    model::DarcyModelMono{N,T},
    ops::DiffusionOps{N,T},
    cap::AssembledCapacity{N,T};
    t::T,
) where {N,T}
    nt = cap.ntotal
    out = zeros(T, N * nt)
    model.variable === :head && return out

    LI = LinearIndices(cap.nnodes)
    cω = cap.C_ω
    @inbounds for d in 1:N
        offset = (d - 1) * nt
        xlow = convert(T, cap.xyz[d][1])
        for I in CartesianIndices(cap.nnodes)
            lin = LI[I]
            if any(k -> I[k] == cap.nnodes[k], 1:N)
                out[offset + lin] = zero(T)
                continue
            end
            xface = if I[d] == 1
                SVector{N,T}(ntuple(k -> (k == d ? xlow : cω[lin][k]), N))
            else
                Iminus = CartesianIndex(ntuple(k -> (k == d ? I[k] - 1 : I[k]), N))
                linm = LI[Iminus]
                SVector{N,T}(ntuple(k -> (cω[lin][k] + cω[linm][k]) / T(2), N))
            end
            ρv = _eval_density(model.ρ, xface, t)
            gv = _eval_gravity(model.gravity, xface, t)
            out[offset + lin] = ρv * gv[d]
        end
    end
    return out
end

function face_bodyforce_values(
    model::DarcyModelDiph{N,T},
    ops::DiffusionOps{N,T},
    cap::AssembledCapacity{N,T};
    t::T,
    phase::Int,
) where {N,T}
    nt = cap.ntotal
    out = zeros(T, N * nt)
    model.variable === :head && return out

    ρ = phase == 1 ? model.ρ1 : model.ρ2

    LI = LinearIndices(cap.nnodes)
    cω = cap.C_ω
    @inbounds for d in 1:N
        offset = (d - 1) * nt
        xlow = convert(T, cap.xyz[d][1])
        for I in CartesianIndices(cap.nnodes)
            lin = LI[I]
            if any(k -> I[k] == cap.nnodes[k], 1:N)
                out[offset + lin] = zero(T)
                continue
            end
            xface = if I[d] == 1
                SVector{N,T}(ntuple(k -> (k == d ? xlow : cω[lin][k]), N))
            else
                Iminus = CartesianIndex(ntuple(k -> (k == d ? I[k] - 1 : I[k]), N))
                linm = LI[Iminus]
                SVector{N,T}(ntuple(k -> (cω[lin][k] + cω[linm][k]) / T(2), N))
            end
            ρv = _eval_density(ρ, xface, t)
            gv = _eval_gravity(model.gravity, xface, t)
            out[offset + lin] = ρv * gv[d]
        end
    end
    return out
end
