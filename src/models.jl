@inline function _normalize_variable(variable::Symbol)
    variable === :pressure && return variable
    variable === :head && return variable
    throw(ArgumentError("unknown variable `$variable`; expected :pressure or :head"))
end

struct DarcyModelMono{N,T,ΛT,ST,ΣT,RT,GT,WT,IT}
    ops::DiffusionOps{N,T}
    cap::AssembledCapacity{N,T}
    λ::ΛT
    source::ST
    storage::ΣT
    ρ::RT
    gravity::GT
    wells::WT
    bc_border::BorderConditions
    bc_interface::IT
    layout::UnknownLayout
    coeff_mode::Symbol
    variable::Symbol
end

function DarcyModelMono(
    cap::AssembledCapacity{N,T},
    ops::DiffusionOps{N,T},
    λ;
    source=((args...) -> zero(T)),
    storage=one(T),
    ρ=one(T),
    gravity=ntuple(_ -> zero(T), N),
    wells::AbstractVector=Any[],
    bc_border::BorderConditions=BorderConditions(),
    bc_interface::Union{Nothing,PenguinBCs.Robin}=nothing,
    layout::UnknownLayout=layout_mono(cap.ntotal),
    coeff_mode::Symbol=:harmonic,
    variable::Symbol=:pressure,
) where {N,T}
    coeff_mode_eff = _normalize_coeff_mode(coeff_mode)
    variable_eff = _normalize_variable(variable)
    return DarcyModelMono{N,T,typeof(λ),typeof(source),typeof(storage),typeof(ρ),typeof(gravity),typeof(wells),typeof(bc_interface)}(
        ops,
        cap,
        λ,
        source,
        storage,
        ρ,
        gravity,
        wells,
        bc_border,
        bc_interface,
        layout,
        coeff_mode_eff,
        variable_eff,
    )
end
