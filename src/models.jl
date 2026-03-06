struct DarcyModelMono{N,T,LT,ST,ΣT,IT}
    ops::DiffusionOps{N,T}
    cap::AssembledCapacity{N,T}
    λ::LT
    source::ST
    storage::ΣT
    bc_border::BorderConditions
    bc_interface::IT
    layout::UnknownLayout
    coeff_mode::Symbol
end

function DarcyModelMono(
    cap::AssembledCapacity{N,T},
    ops::DiffusionOps{N,T},
    λ;
    source=((args...) -> zero(T)),
    storage=one(T),
    bc_border::BorderConditions=BorderConditions(),
    bc_interface::Union{Nothing,PenguinBCs.Robin}=nothing,
    layout::UnknownLayout=layout_mono(cap.ntotal),
    coeff_mode::Symbol=:harmonic,
) where {N,T}
    coeff_mode_eff = _normalize_coeff_mode(coeff_mode)
    return DarcyModelMono{N,T,typeof(λ),typeof(source),typeof(storage),typeof(bc_interface)}(
        ops, cap, λ, source, storage, bc_border, bc_interface, layout, coeff_mode_eff
    )
end
