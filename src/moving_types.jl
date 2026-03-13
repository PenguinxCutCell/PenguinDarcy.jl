abstract type AbstractDarcyTracker end

"""
Nonlinear interface-coupling history for one implicit moving-interface step.
"""
struct InterfaceIterationHistory{T}
    residual_inf::Vector{T}
    update_inf::Vector{T}
    velocity_jump_inf::Vector{T}
    damping::Vector{T}
    converged::Bool
end

"""
Per-time-step diagnostics for moving Darcy free-boundary solves.
"""
struct MovingStepDiagnostics{T,VT}
    t::T
    dt::T
    phase_volume::VT
    volume_change::VT
    interface_measure::T
    global_mass_residual::T
    interface_kinematic_residual::T
    max_normal_velocity::T
    max_normal_velocity_jump::T
    capillary_jump_min::T
    capillary_jump_max::T
    interface_cfl::T
    iteration::InterfaceIterationHistory{T}
end

"""
Return object for `solve_unsteady_moving!`.
"""
struct MovingDarcySolution{T,MT,DT}
    model::MT
    times::Vector{T}
    pressure_states::Vector{Vector{T}}
    interface_states::Vector{Any}
    step_diagnostics::Vector{DT}
    system::LinearSystem{T}
end

struct MovingDarcyModelMono{N,T,TR,LT,ST,ΣT,RT,GT,WT,PET,STT,PT}
    tracker::TR
    λ::LT
    source::ST
    storage::ΣT
    ρ::RT
    gravity::GT
    wells::WT
    bc_border::BorderConditions
    p_ext::PET
    surface_tension::STT
    interface_porosity::PT
    coeff_mode::Symbol
    variable::Symbol
    geom_method::Symbol
end

struct MovingDarcyModelDiph{N,T,TR,Λ1T,S1T,Σ1T,R1T,W1T,Λ2T,S2T,Σ2T,R2T,W2T,GT,IT,STT,PT}
    tracker::TR

    λ1::Λ1T
    source1::S1T
    storage1::Σ1T
    ρ1::R1T
    wells1::W1T

    λ2::Λ2T
    source2::S2T
    storage2::Σ2T
    ρ2::R2T
    wells2::W2T

    gravity::GT
    bc_border::BorderConditions
    bc_interface::IT
    surface_tension::STT
    interface_porosity::PT
    coeff_mode::Symbol
    variable::Symbol
    geom_method::Symbol
end

@inline _is_zero_storage_value(v::Number) = iszero(v)
@inline _is_zero_storage_value(v::Tuple) = all(_is_zero_storage_value, v)
_is_zero_storage_value(v) = false

function _check_moving_storage(storage)
    _is_zero_storage_value(storage) && return storage
    throw(ArgumentError("moving Darcy models currently support quasi-steady solves only; storage must be zero for moving geometry"))
end

function _check_interface_porosity(interface_porosity)
    if interface_porosity isa Number
        interface_porosity > 0 || throw(ArgumentError("interface_porosity must be positive"))
    end
    return interface_porosity
end

function _check_tracker_periodicity_compat(tracker::AbstractDarcyTracker, bc_border::BorderConditions)
    N, _ = tracker_dimension_type(tracker)
    pflags = periodic_flags(bc_border, N)
    if tracker_periodic_transverse(tracker)
        for d in tracker_transverse_dims(tracker)
            pflags[d] || throw(ArgumentError("tracker periodic_transverse=true requires periodic boundary conditions on both sides of transverse dimension $d"))
        end
    end
    return nothing
end

function MovingDarcyModelMono(
    tracker::AbstractDarcyTracker,
    λ;
    source=((args...) -> 0.0),
    storage=0.0,
    ρ=1.0,
    gravity=nothing,
    wells::AbstractVector=Any[],
    bc_border::BorderConditions=BorderConditions(),
    p_ext=0.0,
    surface_tension=0.0,
    interface_porosity=1.0,
    coeff_mode::Symbol=:harmonic,
    variable::Symbol=:pressure,
    geom_method::Symbol=:vofijul,
)
    N, T = tracker_dimension_type(tracker)
    gravity_eff = gravity === nothing ? ntuple(_ -> zero(T), N) : gravity
    coeff_mode_eff = _normalize_coeff_mode(coeff_mode)
    variable_eff = _normalize_variable(variable)
    storage_eff = _check_moving_storage(storage)
    porosity_eff = _check_interface_porosity(interface_porosity)
    _check_tracker_periodicity_compat(tracker, bc_border)

    return MovingDarcyModelMono{N,T,typeof(tracker),typeof(λ),typeof(source),typeof(storage_eff),typeof(ρ),typeof(gravity_eff),typeof(wells),typeof(p_ext),typeof(surface_tension),typeof(porosity_eff)}(
        tracker,
        λ,
        source,
        storage_eff,
        ρ,
        gravity_eff,
        wells,
        bc_border,
        p_ext,
        surface_tension,
        porosity_eff,
        coeff_mode_eff,
        variable_eff,
        geom_method,
    )
end

function MovingDarcyModelDiph(
    tracker::AbstractDarcyTracker,
    λ1,
    λ2;
    source=((args...) -> (0.0, 0.0)),
    storage=0.0,
    ρ=1.0,
    gravity=nothing,
    wells1::AbstractVector=Any[],
    wells2::AbstractVector=Any[],
    bc_border::BorderConditions=BorderConditions(),
    bc_interface::AbstractDarcyInterfaceBC=DarcyContinuity(),
    surface_tension=0.0,
    interface_porosity=1.0,
    coeff_mode::Symbol=:harmonic,
    variable::Symbol=:pressure,
    geom_method::Symbol=:vofijul,
)
    N, T = tracker_dimension_type(tracker)
    s1, s2 = _split_diph_source(source, T)
    Σ1, Σ2 = _split_scalar_or_tuple(storage)
    ρ1, ρ2 = _split_scalar_or_tuple(ρ)

    gravity_eff = gravity === nothing ? ntuple(_ -> zero(T), N) : gravity
    coeff_mode_eff = _normalize_coeff_mode(coeff_mode)
    variable_eff = _normalize_variable(variable)
    _check_moving_storage((Σ1, Σ2))
    porosity_eff = _check_interface_porosity(interface_porosity)
    _check_tracker_periodicity_compat(tracker, bc_border)

    if bc_interface isa DarcyContinuity
        iszero(bc_interface.flux_jump) || throw(ArgumentError("moving two-phase Darcy currently requires flux_jump = 0 (interfacial mass transfer is not implemented)"))
    elseif bc_interface isa DarcyMembrane
        throw(ArgumentError("DarcyMembrane interface law is not implemented yet for moving two-phase Darcy free boundaries"))
    end

    return MovingDarcyModelDiph{N,T,typeof(tracker),typeof(λ1),typeof(s1),typeof(Σ1),typeof(ρ1),typeof(wells1),typeof(λ2),typeof(s2),typeof(Σ2),typeof(ρ2),typeof(wells2),typeof(gravity_eff),typeof(bc_interface),typeof(surface_tension),typeof(porosity_eff)}(
        tracker,
        λ1,
        s1,
        Σ1,
        ρ1,
        wells1,
        λ2,
        s2,
        Σ2,
        ρ2,
        wells2,
        gravity_eff,
        bc_border,
        bc_interface,
        surface_tension,
        porosity_eff,
        coeff_mode_eff,
        variable_eff,
        geom_method,
    )
end

interface_iteration_history(model::Union{MovingDarcyModelMono,MovingDarcyModelDiph}) = interface_iteration_history(model.tracker)

function update_interface!(model::Union{MovingDarcyModelMono,MovingDarcyModelDiph}, state_new; kwargs...)
    return update_interface!(model.tracker, state_new; kwargs...)
end
