struct DarcyModelDiph{N,T,Λ1T,S1T,Σ1T,R1T,W1T,Λ2T,S2T,Σ2T,R2T,W2T,GT,IT}
    ops1::DiffusionOps{N,T}
    cap1::AssembledCapacity{N,T}
    λ1::Λ1T
    source1::S1T
    storage1::Σ1T
    ρ1::R1T
    wells1::W1T

    ops2::DiffusionOps{N,T}
    cap2::AssembledCapacity{N,T}
    λ2::Λ2T
    source2::S2T
    storage2::Σ2T
    ρ2::R2T
    wells2::W2T

    gravity::GT
    bc_border::BorderConditions
    bc_interface::IT
    layout::UnknownLayout
    coeff_mode::Symbol
    variable::Symbol
end

function _split_diph_source(source, ::Type{T}) where {T}
    if source isa Tuple && length(source) == 2
        return source[1], source[2]
    elseif source isa Function
        s1 = (args...) -> begin
            val = applicable(source, args...) ? source(args...) : source(args[1:(end - 1)]...)
            val[1]
        end
        s2 = (args...) -> begin
            val = applicable(source, args...) ? source(args...) : source(args[1:(end - 1)]...)
            val[2]
        end
        return s1, s2
    end
    throw(ArgumentError("diph source must be a 2-tuple or a function returning a 2-tuple"))
end

function _split_scalar_or_tuple(v)
    if v isa Tuple
        length(v) == 2 || throw(ArgumentError("expected 2-tuple"))
        return v[1], v[2]
    end
    return v, v
end

function DarcyModelDiph(
    cap1::AssembledCapacity{N,T},
    ops1::DiffusionOps{N,T},
    λ1,
    cap2::AssembledCapacity{N,T},
    ops2::DiffusionOps{N,T},
    λ2;
    source=((args...) -> (zero(T), zero(T))),
    storage=one(T),
    ρ=one(T),
    gravity=ntuple(_ -> zero(T), N),
    wells1::AbstractVector=Any[],
    wells2::AbstractVector=Any[],
    bc_border::BorderConditions=BorderConditions(),
    bc_interface::AbstractDarcyInterfaceBC=DarcyContinuity(),
    layout::UnknownLayout=layout_diph(cap1.ntotal),
    coeff_mode::Symbol=:harmonic,
    variable::Symbol=:pressure,
) where {N,T}
    cap1.ntotal == cap2.ntotal || throw(ArgumentError("cap1 and cap2 must have identical ntotal"))
    cap1.nnodes == cap2.nnodes || throw(ArgumentError("cap1 and cap2 must have identical nnodes"))

    s1, s2 = _split_diph_source(source, T)
    Σ1, Σ2 = _split_scalar_or_tuple(storage)
    ρ1, ρ2 = _split_scalar_or_tuple(ρ)

    coeff_mode_eff = _normalize_coeff_mode(coeff_mode)
    variable_eff = _normalize_variable(variable)

    return DarcyModelDiph{N,T,typeof(λ1),typeof(s1),typeof(Σ1),typeof(ρ1),typeof(wells1),typeof(λ2),typeof(s2),typeof(Σ2),typeof(ρ2),typeof(wells2),typeof(gravity),typeof(bc_interface)}(
        ops1,
        cap1,
        λ1,
        s1,
        Σ1,
        ρ1,
        wells1,
        ops2,
        cap2,
        λ2,
        s2,
        Σ2,
        ρ2,
        wells2,
        gravity,
        bc_border,
        bc_interface,
        layout,
        coeff_mode_eff,
        variable_eff,
    )
end
