function _stalled_damping(prev_res::T, cur_res::T, ω::T) where {T}
    if cur_res > prev_res * (one(T) - sqrt(eps(T)))
        return max(convert(T, 0.1), ω / 2)
    end
    return ω
end

function _solve_moving_step_mono!(
    model::MovingDarcyModelMono{N,T},
    t::T,
    dt::T;
    method::Symbol=:direct,
    kwargs...,
) where {N,T}
    tracker = model.tracker
    x_old = copy(tracker.xf)
    x_guess = predictor_interface_state(tracker, dt)
    _clamp_tracker_state!(x_guess, tracker)

    residual_hist = T[]
    update_hist = T[]
    jump_hist = T[]
    damping_hist = T[]

    converged = false
    ω = tracker.damping

    sys_last = nothing
    dmodel_last = nothing
    asm_last = nothing
    vn_last = nothing

    for k in 1:tracker.max_iter
        tracker.xf .= x_guess
        tracker.periodic_transverse && GlobalHeightFunctions.ensure_periodic!(tracker.xf)
        check_graph_validity(tracker)

        asm = assemble_unsteady_moving!(model; t=t + dt)
        dmodel = asm.darcy_model
        sys = solve_steady!(dmodel; t=t + dt, method=method, kwargs...)

        flux = recover_flux(dmodel, sys.x; t=t + dt)
        qn = dmodel.ops.H' * flux.faces
        vn_graph = interface_normal_velocity(tracker, dmodel.cap, qn)

        rhs = _graph_rhs_from_normal_velocity(tracker, vn_graph, model.interface_porosity, t + dt)
        R = _interface_residual(x_guess, x_old, dt, rhs)
        δx = -ω .* R

        x_new = x_guess .+ δx
        _clamp_tracker_state!(x_new, tracker)

        res = _inf_norm(R)
        upd = _inf_norm(δx)

        push!(residual_hist, res)
        push!(update_hist, upd)
        push!(jump_hist, zero(T))
        push!(damping_hist, ω)

        sys_last = sys
        dmodel_last = dmodel
        asm_last = asm
        vn_last = vn_graph

        if k > 1
            ω = _stalled_damping(residual_hist[end - 1], res, ω)
        end

        x_guess .= x_new

        if res <= tracker.tol_interface && upd <= tracker.tol_update
            converged = true
            break
        end
    end

    tracker.xf .= x_old
    tracker.periodic_transverse && GlobalHeightFunctions.ensure_periodic!(tracker.xf)

    iter = InterfaceIterationHistory(residual_hist, update_hist, jump_hist, damping_hist, converged)
    update_interface!(tracker, x_guess; dt=dt, iteration=iter)

    final_res = isempty(residual_hist) ? typemax(T) : residual_hist[end]
    final_upd = isempty(update_hist) ? typemax(T) : update_hist[end]
    vmax = vn_last === nothing ? zero(T) : _inf_norm(vn_last)

    return (
        system=sys_last,
        darcy_model=dmodel_last,
        assembled=asm_last,
        vn_graph=vn_last,
        iteration=iter,
        residual=final_res,
        update=final_upd,
        max_normal_velocity=vmax,
        max_normal_velocity_jump=zero(T),
    )
end

function _solve_moving_step_diph!(
    model::MovingDarcyModelDiph{N,T},
    t::T,
    dt::T;
    method::Symbol=:direct,
    kwargs...,
) where {N,T}
    tracker = model.tracker
    x_old = copy(tracker.xf)
    x_guess = predictor_interface_state(tracker, dt)
    _clamp_tracker_state!(x_guess, tracker)

    residual_hist = T[]
    update_hist = T[]
    jump_hist = T[]
    damping_hist = T[]

    converged = false
    ω = tracker.damping

    sys_last = nothing
    dmodel_last = nothing
    asm_last = nothing
    vn1_last = nothing
    vn2_last = nothing

    for k in 1:tracker.max_iter
        tracker.xf .= x_guess
        tracker.periodic_transverse && GlobalHeightFunctions.ensure_periodic!(tracker.xf)
        check_graph_validity(tracker)

        asm = assemble_unsteady_moving!(model; t=t + dt)
        dmodel = asm.darcy_model
        sys = solve_steady!(dmodel; t=t + dt, method=method, kwargs...)

        flux = recover_flux(dmodel, sys.x; t=t + dt)
        q1 = dmodel.ops1.H' * flux.phase1.faces
        q2 = dmodel.ops2.H' * flux.phase2.faces

        vn1_graph = interface_normal_velocity(tracker, dmodel.cap1, q1)
        vn2_graph = interface_normal_velocity(tracker, dmodel.cap2, q2)
        vjump = _inf_norm(vn1_graph .+ vn2_graph)

        rhs = _graph_rhs_from_normal_velocity(tracker, vn1_graph, model.interface_porosity, t + dt)
        R = _interface_residual(x_guess, x_old, dt, rhs)
        δx = -ω .* R

        x_new = x_guess .+ δx
        _clamp_tracker_state!(x_new, tracker)

        res = _inf_norm(R)
        upd = _inf_norm(δx)

        push!(residual_hist, res)
        push!(update_hist, upd)
        push!(jump_hist, vjump)
        push!(damping_hist, ω)

        sys_last = sys
        dmodel_last = dmodel
        asm_last = asm
        vn1_last = vn1_graph
        vn2_last = vn2_graph

        if k > 1
            ω = _stalled_damping(residual_hist[end - 1], res, ω)
        end

        x_guess .= x_new

        if res <= tracker.tol_interface && upd <= tracker.tol_update
            converged = true
            break
        end
    end

    tracker.xf .= x_old
    tracker.periodic_transverse && GlobalHeightFunctions.ensure_periodic!(tracker.xf)

    iter = InterfaceIterationHistory(residual_hist, update_hist, jump_hist, damping_hist, converged)
    update_interface!(tracker, x_guess; dt=dt, iteration=iter)

    final_res = isempty(residual_hist) ? typemax(T) : residual_hist[end]
    final_upd = isempty(update_hist) ? typemax(T) : update_hist[end]
    vmax = vn1_last === nothing ? zero(T) : _inf_norm(vn1_last)
    vmax_jump = isempty(jump_hist) ? zero(T) : jump_hist[end]

    return (
        system=sys_last,
        darcy_model=dmodel_last,
        assembled=asm_last,
        vn1_graph=vn1_last,
        vn2_graph=vn2_last,
        iteration=iter,
        residual=final_res,
        update=final_upd,
        max_normal_velocity=vmax,
        max_normal_velocity_jump=vmax_jump,
    )
end

function _initial_moving_state(
    model::MovingDarcyModelMono{N,T};
    t::T,
    method::Symbol=:direct,
    kwargs...,
) where {N,T}
    asm0 = assemble_unsteady_moving!(model; t=t)
    sys0 = solve_steady!(asm0.darcy_model; t=t, method=method, kwargs...)
    return asm0, sys0
end

function _initial_moving_state(
    model::MovingDarcyModelDiph{N,T};
    t::T,
    method::Symbol=:direct,
    kwargs...,
) where {N,T}
    asm0 = assemble_unsteady_moving!(model; t=t)
    sys0 = solve_steady!(asm0.darcy_model; t=t, method=method, kwargs...)
    return asm0, sys0
end

function solve_unsteady_moving!(
    model::MovingDarcyModelMono{N,T},
    tspan::Tuple{<:Real,<:Real};
    dt::Real,
    method::Symbol=:direct,
    save_history::Bool=true,
    kwargs...,
) where {N,T}
    t0 = convert(T, tspan[1])
    tend = convert(T, tspan[2])
    tend >= t0 || throw(ArgumentError("tspan must satisfy tend >= t0"))

    dt0 = convert(T, dt)
    dt0 > zero(T) || throw(ArgumentError("dt must be positive"))

    tracker = model.tracker
    check_graph_validity(tracker)

    asm0, sys0 = _initial_moving_state(model; t=t0, method=method, kwargs...)
    cap0 = asm0.darcy_model.cap
    vol_old = _phase_volume(cap0)

    times = T[t0]
    states = Vector{Vector{T}}()
    interfaces = Vector{Any}()
    diags = Vector{MovingStepDiagnostics{T,T}}()

    if save_history
        push!(states, copy(sys0.x))
        push!(interfaces, copy(tracker.xf))
    end

    sys_last = sys0
    t = t0
    tol = sqrt(eps(T)) * max(one(T), abs(t0), abs(tend))

    while t < tend - tol
        dt_step = min(dt0, tend - t)
        step = _solve_moving_step_mono!(model, t, dt_step; method=method, kwargs...)
        t += dt_step

        dmodel = step.darcy_model
        cap = dmodel.cap
        vol_new = _phase_volume(cap)
        dV = vol_new - vol_old
        imeas = _interface_measure(cap)
        mass_res = _moving_mass_residual_mono(dmodel, step.system.x, vol_old, vol_new, dt_step, t)
        cfl = _interface_cfl_recommendation(tracker, dt_step, step.max_normal_velocity)

        diag = MovingStepDiagnostics(
            t,
            dt_step,
            vol_new,
            dV,
            imeas,
            mass_res,
            step.residual,
            step.max_normal_velocity,
            step.max_normal_velocity_jump,
            step.assembled.capillary_jump_min,
            step.assembled.capillary_jump_max,
            cfl,
            step.iteration,
        )
        push!(diags, diag)

        push!(times, t)
        if save_history
            push!(states, copy(step.system.x))
            push!(interfaces, copy(tracker.xf))
        end

        vol_old = vol_new
        sys_last = step.system
    end

    if !save_history
        states = [copy(sys_last.x)]
        interfaces = [copy(tracker.xf)]
        times = T[t]
    end

    return MovingDarcySolution{T,typeof(model),MovingStepDiagnostics{T,T}}(
        model,
        times,
        states,
        interfaces,
        diags,
        sys_last,
    )
end

function solve_unsteady_moving!(
    model::MovingDarcyModelDiph{N,T},
    tspan::Tuple{<:Real,<:Real};
    dt::Real,
    method::Symbol=:direct,
    save_history::Bool=true,
    kwargs...,
) where {N,T}
    t0 = convert(T, tspan[1])
    tend = convert(T, tspan[2])
    tend >= t0 || throw(ArgumentError("tspan must satisfy tend >= t0"))

    dt0 = convert(T, dt)
    dt0 > zero(T) || throw(ArgumentError("dt must be positive"))

    tracker = model.tracker
    check_graph_validity(tracker)

    asm0, sys0 = _initial_moving_state(model; t=t0, method=method, kwargs...)
    cap10 = asm0.darcy_model.cap1
    cap20 = asm0.darcy_model.cap2
    vol1_old = _phase_volume(cap10)
    vol2_old = _phase_volume(cap20)

    times = T[t0]
    states = Vector{Vector{T}}()
    interfaces = Vector{Any}()
    diags = Vector{MovingStepDiagnostics{T,NTuple{2,T}}}()

    if save_history
        push!(states, copy(sys0.x))
        push!(interfaces, copy(tracker.xf))
    end

    sys_last = sys0
    t = t0
    tol = sqrt(eps(T)) * max(one(T), abs(t0), abs(tend))

    while t < tend - tol
        dt_step = min(dt0, tend - t)
        step = _solve_moving_step_diph!(model, t, dt_step; method=method, kwargs...)
        t += dt_step

        dmodel = step.darcy_model
        cap1 = dmodel.cap1
        cap2 = dmodel.cap2
        vol1_new = _phase_volume(cap1)
        vol2_new = _phase_volume(cap2)
        imeas = _interface_measure(cap1)

        mass_res = _moving_mass_residual_diph(
            dmodel,
            step.system.x,
            vol1_old,
            vol1_new,
            vol2_old,
            vol2_new,
            dt_step,
            t,
        )
        cfl = _interface_cfl_recommendation(tracker, dt_step, step.max_normal_velocity)

        diag = MovingStepDiagnostics(
            t,
            dt_step,
            (vol1_new, vol2_new),
            (vol1_new - vol1_old, vol2_new - vol2_old),
            imeas,
            mass_res.total,
            step.residual,
            step.max_normal_velocity,
            step.max_normal_velocity_jump,
            step.assembled.capillary_jump_min,
            step.assembled.capillary_jump_max,
            cfl,
            step.iteration,
        )
        push!(diags, diag)

        push!(times, t)
        if save_history
            push!(states, copy(step.system.x))
            push!(interfaces, copy(tracker.xf))
        end

        vol1_old = vol1_new
        vol2_old = vol2_new
        sys_last = step.system
    end

    if !save_history
        states = [copy(sys_last.x)]
        interfaces = [copy(tracker.xf)]
        times = T[t]
    end

    return MovingDarcySolution{T,typeof(model),MovingStepDiagnostics{T,NTuple{2,T}}}(
        model,
        times,
        states,
        interfaces,
        diags,
        sys_last,
    )
end
