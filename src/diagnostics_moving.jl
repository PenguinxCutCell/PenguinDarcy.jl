function _phase_volume(cap::AssembledCapacity{N,T}) where {N,T}
    vol = zero(T)
    @inbounds for i in 1:cap.ntotal
        Vi = cap.buf.V[i]
        if isfinite(Vi) && Vi > zero(T)
            vol += Vi
        end
    end
    return vol
end

function _interface_measure(cap::AssembledCapacity{N,T}) where {N,T}
    meas = zero(T)
    @inbounds for i in 1:cap.ntotal
        Γi = cap.buf.Γ[i]
        if isfinite(Γi) && Γi > zero(T)
            meas += Γi
        end
    end
    return meas
end

function _sum_outer_boundary_discharge(model::DarcyModelMono{N,T}, state, t::T) where {N,T}
    q = zero(T)
    for side in _all_sides(N)
        if get(model.bc_border.borders, side, nothing) isa Periodic
            continue
        end
        q += boundary_discharge(model, state, side; t=t)
    end
    return q
end

function _sum_outer_boundary_discharge(model::DarcyModelDiph{N,T}, state, t::T) where {N,T}
    q1 = zero(T)
    q2 = zero(T)
    for side in _all_sides(N)
        if get(model.bc_border.borders, side, nothing) isa Periodic
            continue
        end
        q1 += boundary_discharge(model, state, side; t=t, phase=1)
        q2 += boundary_discharge(model, state, side; t=t, phase=2)
    end
    return q1, q2
end

function _moving_mass_residual_mono(
    model::DarcyModelMono{N,T},
    state,
    vol_old::T,
    vol_new::T,
    dt::T,
    t::T,
) where {N,T}
    src = integrated_source(model; t=t)
    qw = integrated_well_rate(model; t=t)
    qbd = _sum_outer_boundary_discharge(model, state, t)
    return (vol_new - vol_old) / dt - (src + qw - qbd)
end

function _moving_mass_residual_diph(
    model::DarcyModelDiph{N,T},
    state,
    vol1_old::T,
    vol1_new::T,
    vol2_old::T,
    vol2_new::T,
    dt::T,
    t::T,
) where {N,T}
    src = integrated_source(model; t=t)
    qw = integrated_well_rate(model; t=t)
    qbd1, qbd2 = _sum_outer_boundary_discharge(model, state, t)

    res1 = (vol1_new - vol1_old) / dt - (src.phase1 + qw.phase1 - qbd1)
    res2 = (vol2_new - vol2_old) / dt - (src.phase2 + qw.phase2 - qbd2)
    return (phase1=res1, phase2=res2, total=res1 + res2)
end

function _interface_cfl_recommendation(
    tracker::HeightFunctionTracker{N,T},
    dt::T,
    vmax::T,
) where {N,T}
    Δmin = minimum(meshsize(tracker.grid))
    Δmin > zero(T) || return zero(T)
    return dt * vmax / Δmin
end

interface_mass_residual(diag::MovingStepDiagnostics) = diag.global_mass_residual
interface_iteration_history(diag::MovingStepDiagnostics) = diag.iteration

function interface_mass_residual(sol::MovingDarcySolution)
    return [interface_mass_residual(d) for d in sol.step_diagnostics]
end

function interface_iteration_history(sol::MovingDarcySolution)
    return [interface_iteration_history(d) for d in sol.step_diagnostics]
end
