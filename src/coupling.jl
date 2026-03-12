mutable struct DarcyCoupledModelMono
    base_model::Any
    mobility_update!::Any
    metadata::Dict{Symbol,Any}
    quasi_steady::Bool
end

function DarcyCoupledModelMono(base_model::DarcyModelMono;
                               mobility_update! = nothing,
                               metadata::AbstractDict{Symbol,<:Any}=Dict{Symbol,Any}(),
                               quasi_steady::Bool=true)
    return DarcyCoupledModelMono(base_model, mobility_update!, Dict{Symbol,Any}(metadata), quasi_steady)
end

mutable struct DarcyCoupledStateMono
    pressure
    velocity
    incoming_concentration
    last_time
    last_mobility
end

_darcy_scalar_type(model::DarcyModelMono{N,T}) where {N,T} = T

function _empty_velocity(model::DarcyModelMono{1,T}) where {T}
    nt = model.cap.ntotal
    return (faces=zeros(T, nt), x=zeros(T, nt))
end

function _empty_velocity(model::DarcyModelMono{2,T}) where {T}
    nt = model.cap.ntotal
    return (faces=zeros(T, 2 * nt), x=zeros(T, nt), y=zeros(T, nt))
end

function _empty_velocity(model::DarcyModelMono{3,T}) where {T}
    nt = model.cap.ntotal
    return (faces=zeros(T, 3 * nt), x=zeros(T, nt), y=zeros(T, nt), z=zeros(T, nt))
end

function _darcy_system_size(model::DarcyModelMono)
    lay = model.layout.offsets
    return maximum((last(lay.ω), last(lay.γ)))
end

function _init_get(init, key::Symbol, default)
    if init === nothing
        return default
    elseif init isa NamedTuple
        return hasproperty(init, key) ? getproperty(init, key) : default
    elseif init isa AbstractDict
        return get(init, key, default)
    end
    return default
end

function _without_keys(kwargs::Base.Pairs, drop::Tuple{Vararg{Symbol}})
    kept = Symbol[]
    vals = Any[]
    for key in keys(kwargs)
        key in drop && continue
        push!(kept, key)
        push!(vals, kwargs[key])
    end
    return NamedTuple{Tuple(kept)}(Tuple(vals))
end

function _sanitize_scalar_field(data)
    if data isa AbstractArray
        clean = copy(data)
        @inbounds for i in eachindex(clean)
            if !isfinite(clean[i])
                clean[i] = zero(eltype(clean))
            end
        end
        return clean
    end
    return deepcopy(data)
end

function _rebuild_with_mobility(model::DarcyModelMono, λ)
    return DarcyModelMono(
        model.cap,
        model.ops,
        λ;
        source=model.source,
        storage=model.storage,
        ρ=model.ρ,
        gravity=model.gravity,
        wells=model.wells,
        bc_border=model.bc_border,
        bc_interface=model.bc_interface,
        layout=model.layout,
        coeff_mode=model.coeff_mode,
        variable=model.variable,
    )
end

function _evaluate_mobility_callback(callback, base_model::DarcyModelMono, concentration, t)
    if applicable(callback, base_model, concentration, t)
        return callback(base_model, concentration, t)
    elseif applicable(callback, concentration, t)
        return callback(concentration, t)
    elseif applicable(callback, concentration)
        return callback(concentration)
    end
    throw(ArgumentError(
        "mobility_update! must accept (base_model, concentration, t), (concentration, t), or (concentration)."))
end

function _maybe_update_mobility!(coupled::DarcyCoupledModelMono, state::DarcyCoupledStateMono, t)
    callback = coupled.mobility_update!
    (callback === nothing || state.incoming_concentration === nothing) && return nothing

    result = _evaluate_mobility_callback(callback, coupled.base_model, state.incoming_concentration, t)
    if result !== nothing
        coupled.base_model = _rebuild_with_mobility(coupled.base_model, result)
    end

    state.last_mobility = coupled.base_model.λ
    coupled.metadata[:last_mobility] = state.last_mobility
    return state.last_mobility
end

function PenguinSolverCore.initialize_state(model::DarcyCoupledModelMono, init)
    base = model.base_model
    T = _darcy_scalar_type(base)
    nsys = _darcy_system_size(base)

    pressure = _init_get(init, :pressure, zeros(T, nsys))
    velocity = _init_get(init, :velocity, _empty_velocity(base))
    concentration = _init_get(init, :concentration, zeros(T, nsys))
    last_time = _init_get(init, :time, nothing)
    last_mobility = _init_get(init, :mobility, base.λ)

    return DarcyCoupledStateMono(pressure, velocity, concentration, last_time, last_mobility)
end

function PenguinSolverCore.copy_state(state::DarcyCoupledStateMono)
    return DarcyCoupledStateMono(
        deepcopy(state.pressure),
        deepcopy(state.velocity),
        deepcopy(state.incoming_concentration),
        deepcopy(state.last_time),
        deepcopy(state.last_mobility),
    )
end

function PenguinSolverCore.advance_steady!(block::CoupledBlock{M,S,C}; kwargs...) where {M<:DarcyCoupledModelMono,S<:DarcyCoupledStateMono,C}
    coupled = block.model
    state = block.state
    T = _darcy_scalar_type(coupled.base_model)

    tval = haskey(kwargs, :t) ? convert(T, kwargs[:t]) : zero(T)
    _maybe_update_mobility!(coupled, state, tval)

    solve_kwargs = _without_keys(kwargs, (:t,))
    sys = solve_steady!(coupled.base_model; t=tval, solve_kwargs...)

    state.pressure = copy(sys.x)
    state.velocity = recover_velocity(coupled.base_model, state.pressure; t=tval)
    state.last_time = tval

    if block.cache isa AbstractDict
        block.cache[:darcy_system] = sys
    end

    return block
end

function PenguinSolverCore.advance_unsteady!(block::CoupledBlock{M,S,C}, t, dt; kwargs...) where {M<:DarcyCoupledModelMono,S<:DarcyCoupledStateMono,C}
    coupled = block.model
    coupled.quasi_steady || throw(ArgumentError(
        "DarcyCoupledModelMono currently supports quasi-steady stepping only."))

    T = _darcy_scalar_type(coupled.base_model)
    t_eval = convert(T, t) + convert(T, dt)
    _maybe_update_mobility!(coupled, block.state, t_eval)

    solve_kwargs = _without_keys(kwargs, (:t, :dt))
    sys = solve_steady!(coupled.base_model; t=t_eval, solve_kwargs...)

    block.state.pressure = copy(sys.x)
    block.state.velocity = recover_velocity(coupled.base_model, block.state.pressure; t=t_eval)
    block.state.last_time = t_eval

    if block.cache isa AbstractDict
        block.cache[:darcy_system] = sys
        block.cache[:quasi_steady] = true
    end

    return block
end

function PenguinSolverCore.get_coupling_field(block::CoupledBlock{M,S,C}, ::Val{:velocity}) where {M<:DarcyCoupledModelMono,S<:DarcyCoupledStateMono,C}
    return block.state.velocity
end

function PenguinSolverCore.get_coupling_field(block::CoupledBlock{M,S,C}, ::Val{:concentration}) where {M<:DarcyCoupledModelMono,S<:DarcyCoupledStateMono,C}
    return _sanitize_scalar_field(block.state.incoming_concentration)
end

function PenguinSolverCore.set_coupling_field!(block::CoupledBlock{M,S,C}, ::Val{:concentration}, data) where {M<:DarcyCoupledModelMono,S<:DarcyCoupledStateMono,C}
    block.state.incoming_concentration = _sanitize_scalar_field(data)
    return block
end

function PenguinSolverCore.set_coupling_field!(block::CoupledBlock{M,S,C}, ::Val{:velocity}, data) where {M<:DarcyCoupledModelMono,S<:DarcyCoupledStateMono,C}
    block.state.velocity = deepcopy(data)
    return block
end

function PenguinSolverCore.block_summary(block::CoupledBlock{M,S,C}) where {M<:DarcyCoupledModelMono,S<:DarcyCoupledStateMono,C}
    qs = block.model.quasi_steady ? "quasi-steady" : "transient"
    return "DarcyCoupledModelMono(name=$(block.name), mode=$qs)"
end
