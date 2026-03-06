function _extract_mono_pressure_state(model::DarcyModelMono{N,T}, state) where {N,T}
    lay = model.layout.offsets
    nt = model.cap.ntotal
    nsys = maximum((last(lay.ω), last(lay.γ)))

    if length(state) == nsys
        pω = Vector{T}(state[lay.ω])
        pγ = Vector{T}(state[lay.γ])
        return pω, pγ
    elseif length(state) == nt
        pω = Vector{T}(state)
        pγ = zeros(T, nt)
        return pω, pγ
    end

    throw(DimensionMismatch("state length must be $nt (ω block) or $nsys (full ω+γ system)"))
end

function _extract_diph_pressure_state(model::DarcyModelDiph{N,T}, state) where {N,T}
    lay = model.layout.offsets
    nt = model.cap1.ntotal
    nsys = maximum((last(lay.ω1), last(lay.γ1), last(lay.ω2), last(lay.γ2)))

    if length(state) == nsys
        p1ω = Vector{T}(state[lay.ω1])
        p1γ = Vector{T}(state[lay.γ1])
        p2ω = Vector{T}(state[lay.ω2])
        p2γ = Vector{T}(state[lay.γ2])
        return p1ω, p1γ, p2ω, p2γ
    elseif length(state) == 2 * nt
        s = Vector{T}(state)
        p1ω = s[1:nt]
        p2ω = s[(nt + 1):(2 * nt)]
        return p1ω, zeros(T, nt), p2ω, zeros(T, nt)
    end

    throw(DimensionMismatch("state length must be $(2 * nt) (ω1+ω2) or $nsys (full system)"))
end

@inline function _face_component(uf::AbstractVector{T}, nt::Int, d::Int) where {T}
    return Vector{T}(view(uf, (d - 1) * nt + 1:d * nt))
end

function _faces_named_tuple(uf::Vector{T}, nt::Int, ::Val{1}) where {T}
    return (faces=uf, x=_face_component(uf, nt, 1))
end

function _faces_named_tuple(uf::Vector{T}, nt::Int, ::Val{2}) where {T}
    return (faces=uf, x=_face_component(uf, nt, 1), y=_face_component(uf, nt, 2))
end

function _faces_named_tuple(uf::Vector{T}, nt::Int, ::Val{3}) where {T}
    return (faces=uf, x=_face_component(uf, nt, 1), y=_face_component(uf, nt, 2), z=_face_component(uf, nt, 3))
end

function recover_flux(model::DarcyModelMono{N,T}, state; t::Real=zero(T)) where {N,T}
    tt = convert(T, t)
    pω, pγ = _extract_mono_pressure_state(model, state)
    gfaces = gradient(model.ops, pω, pγ)
    Tf = tensor_flux_operator(model, model.ops, model.cap; t=tt, phase=1)
    bfaces = face_bodyforce_values(model, model.ops, model.cap; t=tt)
    uf = -Tf * gfaces + Tf * bfaces
    return _faces_named_tuple(Vector{T}(uf), model.cap.ntotal, Val(N))
end

function recover_flux(model::DarcyModelDiph{N,T}, state; t::Real=zero(T)) where {N,T}
    tt = convert(T, t)
    p1ω, p1γ, p2ω, p2γ = _extract_diph_pressure_state(model, state)

    g1 = gradient(model.ops1, p1ω, p1γ)
    T1 = tensor_flux_operator(model, model.ops1, model.cap1; t=tt, phase=1)
    b1 = face_bodyforce_values(model, model.ops1, model.cap1; t=tt, phase=1)
    u1 = -T1 * g1 + T1 * b1

    g2 = gradient(model.ops2, p2ω, p2γ)
    T2 = tensor_flux_operator(model, model.ops2, model.cap2; t=tt, phase=2)
    b2 = face_bodyforce_values(model, model.ops2, model.cap2; t=tt, phase=2)
    u2 = -T2 * g2 + T2 * b2

    return (
        phase1=_faces_named_tuple(Vector{T}(u1), model.cap1.ntotal, Val(N)),
        phase2=_faces_named_tuple(Vector{T}(u2), model.cap2.ntotal, Val(N)),
    )
end

recover_velocity(model::Union{DarcyModelMono,DarcyModelDiph}, state; kwargs...) = recover_flux(model, state; kwargs...)
