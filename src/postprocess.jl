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
    λfaces = face_mobility_values(model; t=tt)
    uf = -λfaces .* gfaces
    return _faces_named_tuple(uf, model.cap.ntotal, Val(N))
end

recover_velocity(model::DarcyModelMono, state; kwargs...) = recover_flux(model, state; kwargs...)
