function solve_steady!(model::DarcyModelMono{N,T}; t::Real=zero(T), method::Symbol=:direct, kwargs...) where {N,T}
    tt = convert(T, t)
    n = maximum((last(model.layout.offsets.ω), last(model.layout.offsets.γ)))
    sys = LinearSystem(spzeros(T, n, n), zeros(T, n))
    assemble_steady_mono!(sys, model, tt)
    solve!(sys; method=method, kwargs...)
    return sys
end

function solve_steady!(model::DarcyModelDiph{N,T}; t::Real=zero(T), method::Symbol=:direct, kwargs...) where {N,T}
    tt = convert(T, t)
    lay = model.layout.offsets
    n = maximum((last(lay.ω1), last(lay.γ1), last(lay.ω2), last(lay.γ2)))
    sys = LinearSystem(spzeros(T, n, n), zeros(T, n))
    assemble_steady_diph!(sys, model, tt)
    solve!(sys; method=method, kwargs...)
    return sys
end

function solve_unsteady!(
    model::DarcyModelMono{N,T},
    u0,
    tspan::Tuple{<:Real,<:Real};
    dt::Real,
    scheme=:BE,
    method::Symbol=:direct,
    save_history::Bool=true,
    kwargs...,
) where {N,T}
    t0 = convert(T, tspan[1])
    tend = convert(T, tspan[2])
    tend >= t0 || throw(ArgumentError("tspan must satisfy tend >= t0"))

    dt0 = convert(T, dt)
    dt0 > zero(T) || throw(ArgumentError("dt must be positive"))
    θ = _theta_from_scheme(T, scheme)

    u = _init_unsteady_state_mono(model, u0)
    lay = model.layout.offsets
    nsys = maximum((last(lay.ω), last(lay.γ)))

    times = T[t0]
    states = Vector{Vector{T}}()
    save_history && push!(states, copy(u))

    sys = LinearSystem(spzeros(T, nsys, nsys), zeros(T, nsys); x=copy(u))
    tol = sqrt(eps(T)) * max(one(T), abs(t0), abs(tend))
    t = t0

    while t < tend - tol
        dt_step = min(dt0, tend - t)
        assemble_unsteady_mono!(sys, model, u, t, dt_step, θ)
        solve!(sys; method=method, reuse_factorization=false, kwargs...)
        u .= sys.x
        t += dt_step
        push!(times, t)
        save_history && push!(states, copy(u))
    end

    if !save_history
        states = [copy(u)]
        times = T[t]
    end
    return (times=times, states=states, system=sys, reused_constant_operator=false)
end

function solve_unsteady!(
    model::DarcyModelDiph{N,T},
    u0,
    tspan::Tuple{<:Real,<:Real};
    dt::Real,
    scheme=:BE,
    method::Symbol=:direct,
    save_history::Bool=true,
    kwargs...,
) where {N,T}
    t0 = convert(T, tspan[1])
    tend = convert(T, tspan[2])
    tend >= t0 || throw(ArgumentError("tspan must satisfy tend >= t0"))

    dt0 = convert(T, dt)
    dt0 > zero(T) || throw(ArgumentError("dt must be positive"))
    θ = _theta_from_scheme(T, scheme)

    u = _init_unsteady_state_diph(model, u0)
    lay = model.layout.offsets
    nsys = maximum((last(lay.ω1), last(lay.γ1), last(lay.ω2), last(lay.γ2)))

    times = T[t0]
    states = Vector{Vector{T}}()
    save_history && push!(states, copy(u))

    sys = LinearSystem(spzeros(T, nsys, nsys), zeros(T, nsys); x=copy(u))
    tol = sqrt(eps(T)) * max(one(T), abs(t0), abs(tend))
    t = t0

    while t < tend - tol
        dt_step = min(dt0, tend - t)
        assemble_unsteady_diph!(sys, model, u, t, dt_step, θ)
        solve!(sys; method=method, reuse_factorization=false, kwargs...)
        u .= sys.x
        t += dt_step
        push!(times, t)
        save_history && push!(states, copy(u))
    end

    if !save_history
        states = [copy(u)]
        times = T[t]
    end
    return (times=times, states=states, system=sys, reused_constant_operator=false)
end
