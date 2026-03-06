function _normalize_coeff_mode(mode::Symbol)::Symbol
    if mode === :harmonic || mode === :arithmetic || mode === :face || mode === :cell
        return mode
    end
    throw(ArgumentError("unknown coeff_mode `$mode`; expected :harmonic, :arithmetic, :face, or :cell"))
end

function _eval_fun_or_const(v, x::SVector{N,T}, t::T) where {N,T}
    if v isa Number
        return convert(T, v)
    elseif v isa Function
        if applicable(v, x..., t)
            return convert(T, v(x..., t))
        elseif applicable(v, x...)
            return convert(T, v(x...))
        end
    end
    throw(ArgumentError("callback/value must be numeric, (x...), or (x..., t)"))
end

function _source_values_mono(cap::AssembledCapacity{N,T}, source, t::T) where {N,T}
    out = Vector{T}(undef, cap.ntotal)
    @inbounds for i in eachindex(out)
        out[i] = _eval_fun_or_const(source, cap.C_ω[i], t)
    end
    return out
end

function _storage_values_mono(cap::AssembledCapacity{N,T}, storage, t::T) where {N,T}
    out = Vector{T}(undef, cap.ntotal)
    if storage isa Number
        fill!(out, convert(T, storage))
        return out
    end

    @inbounds for i in eachindex(out)
        out[i] = _eval_fun_or_const(storage, cap.C_ω[i], t)
    end
    return out
end

function _harmonic_mean(a::T, b::T) where {T}
    den = a + b
    return iszero(den) ? zero(T) : (T(2) * a * b) / den
end

function _sample_mobility(cap::AssembledCapacity{N,T}, λ, t::T) where {N,T}
    out = Vector{T}(undef, cap.ntotal)
    @inbounds for i in eachindex(out)
        out[i] = convert(T, eval_coeff(λ, cap.C_ω[i], t, i))
    end
    return out
end

function face_mobility_values(
    cap::AssembledCapacity{N,T},
    λ,
    t::T,
    mode::Symbol,
) where {N,T}
    mode_eff = _normalize_coeff_mode(mode)
    nt = cap.ntotal
    vals = Vector{T}(undef, N * nt)
    LI = LinearIndices(cap.nnodes)
    cω = cap.C_ω

    if mode_eff === :face
        @inbounds for d in 1:N
            offset = (d - 1) * nt
            xlow = convert(T, cap.xyz[d][1])
            for I in CartesianIndices(cap.nnodes)
                lin = LI[I]
                if any(k -> I[k] == cap.nnodes[k], 1:N)
                    vals[offset + lin] = zero(T)
                    continue
                end
                xface = if I[d] == 1
                    SVector{N,T}(ntuple(k -> (k == d ? xlow : cω[lin][k]), N))
                else
                    Iminus = CartesianIndex(ntuple(k -> (k == d ? I[k] - 1 : I[k]), N))
                    linm = LI[Iminus]
                    SVector{N,T}(ntuple(k -> (cω[lin][k] + cω[linm][k]) / T(2), N))
                end
                vals[offset + lin] = convert(T, eval_coeff(λ, xface, t, lin))
            end
        end
        return vals
    end

    kcell = _sample_mobility(cap, λ, t)
    @inbounds for d in 1:N
        offset = (d - 1) * nt
        for I in CartesianIndices(cap.nnodes)
            lin = LI[I]
            if any(k -> I[k] == cap.nnodes[k], 1:N)
                vals[offset + lin] = zero(T)
                continue
            end
            ki = kcell[lin]
            if mode_eff === :cell
                vals[offset + lin] = ki
            elseif I[d] == 1
                vals[offset + lin] = ki
            else
                Iminus = CartesianIndex(ntuple(k -> (k == d ? I[k] - 1 : I[k]), N))
                kim = kcell[LI[Iminus]]
                vals[offset + lin] = mode_eff === :harmonic ? _harmonic_mean(ki, kim) : (ki + kim) / T(2)
            end
        end
    end
    return vals
end

function _diag_tensor_face_values(model, ops::DiffusionOps{N,T}, cap::AssembledCapacity{N,T}; t::T, phase::Int=1) where {N,T}
    parts = ntuple(d -> face_tensor_values(model, ops, cap, d, d; t=t, phase=phase), N)
    return join_face_stack(parts)
end

function face_mobility_values(model::DarcyModelMono{N,T}; t::T=zero(T)) where {N,T}
    mob = _wrap_mobility(model.λ, N)
    if mob isa ScalarMobility
        return face_mobility_values(model.cap, model.λ, t, model.coeff_mode)
    end
    return _diag_tensor_face_values(model, model.ops, model.cap; t=t, phase=1)
end

function face_mobility_values(model::DarcyModelDiph{N,T}; t::T=zero(T), phase::Int=1) where {N,T}
    if phase == 1
        mob = _wrap_mobility(model.λ1, N)
        if mob isa ScalarMobility
            return face_mobility_values(model.cap1, model.λ1, t, model.coeff_mode)
        end
        return _diag_tensor_face_values(model, model.ops1, model.cap1; t=t, phase=1)
    elseif phase == 2
        mob = _wrap_mobility(model.λ2, N)
        if mob isa ScalarMobility
            return face_mobility_values(model.cap2, model.λ2, t, model.coeff_mode)
        end
        return _diag_tensor_face_values(model, model.ops2, model.cap2; t=t, phase=2)
    end
    throw(ArgumentError("phase must be 1 or 2"))
end
