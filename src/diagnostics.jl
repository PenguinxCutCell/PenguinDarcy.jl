function _all_sides(N::Int)
    if N == 1
        return (:left, :right)
    elseif N == 2
        return (:left, :right, :bottom, :top)
    elseif N == 3
        return (:left, :right, :bottom, :top, :backward, :forward)
    end
    throw(ArgumentError("unsupported dimension N=$N; expected 1, 2, or 3"))
end

function _phase_boundary_discharge(
    cap::AssembledCapacity{N,T},
    pω::AbstractVector{T},
    λ,
    ρ,
    gravity,
    variable::Symbol,
    bc_border::BorderConditions,
    side::Symbol,
    t::T,
) where {N,T}
    side_bc = get(bc_border.borders, side, Neumann(0.0))
    side_bc isa Periodic && return zero(T)

    d, is_high, normal_sign = side_info(side, N)
    xyz_d = cap.xyz[d]
    length(xyz_d) >= 2 || throw(ArgumentError("need at least 2 grid nodes in each dimension"))
    Δd = abs(xyz_d[2] - xyz_d[1])
    δ = Δd / T(2)
    x_d = is_high ? xyz_d[end] : xyz_d[1]

    LI = LinearIndices(cap.nnodes)
    out = zero(T)
    for I in each_boundary_cell(cap.nnodes, side)
        lin = LI[I]
        Aface = cap.buf.A[d][lin]
        if !(isfinite(Aface) && !iszero(Aface))
            continue
        end

        Cω = cap.C_ω[lin]
        xface = SVector{N,T}(ntuple(k -> (k == d ? x_d : Cω[k]), N))
        λnn = _normal_mobility_at_face(λ, N, d, xface, t, lin)
        pP = pω[lin]
        bn = if variable === :head
            zero(T)
        else
            ρv = _eval_density(ρ, xface, t)
            gv = _eval_gravity(gravity, xface, t)
            convert(T, normal_sign) * (ρv * gv[d])
        end

        qn = if side_bc isa Dirichlet
            pB = convert(T, eval_bc(side_bc.value, xface, t))
            -λnn * ((pB - pP) / δ - bn)
        elseif side_bc isa Neumann
            convert(T, eval_bc(side_bc.value, xface, t))
        elseif side_bc isa Robin
            α = convert(T, eval_bc(side_bc.α, xface, t))
            β = convert(T, eval_bc(side_bc.β, xface, t))
            g = convert(T, eval_bc(side_bc.value, xface, t))
            den = α * δ - β * λnn
            iszero(den) && throw(ArgumentError("degenerate Darcy Robin boundary at side `$side`: α*δ - β*λ = 0"))
            λnn * (α * (pP + δ * bn) - g) / den
        else
            throw(ArgumentError("unsupported boundary type $(typeof(side_bc))"))
        end

        out += qn * Aface
    end

    return out
end

function boundary_discharge(
    model::DarcyModelMono{N,T},
    state,
    side::Symbol;
    t::Real=zero(T),
) where {N,T}
    tt = convert(T, t)
    pω, _ = _extract_mono_pressure_state(model, state)
    return _phase_boundary_discharge(
        model.cap,
        pω,
        model.λ,
        model.ρ,
        model.gravity,
        model.variable,
        model.bc_border,
        side,
        tt,
    )
end

function boundary_discharge(
    model::DarcyModelDiph{N,T},
    state,
    side::Symbol;
    t::Real=zero(T),
    phase::Int=1,
) where {N,T}
    tt = convert(T, t)
    p1ω, _, p2ω, _ = _extract_diph_pressure_state(model, state)
    if phase == 1
        return _phase_boundary_discharge(model.cap1, p1ω, model.λ1, model.ρ1, model.gravity, model.variable, model.bc_border, side, tt)
    elseif phase == 2
        return _phase_boundary_discharge(model.cap2, p2ω, model.λ2, model.ρ2, model.gravity, model.variable, model.bc_border, side, tt)
    end
    throw(ArgumentError("phase must be 1 or 2"))
end

function interface_discharge(model::DarcyModelMono{N,T}, state; t::Real=zero(T)) where {N,T}
    flux = recover_flux(model, state; t=t)
    q = model.ops.H' * flux.faces
    idxγ = _interface_mask(model.cap)
    total = zero(T)
    for i in eachindex(q)
        idxγ[i] && (total += q[i])
    end
    return (per_cell=q, total=total)
end

function interface_discharge(model::DarcyModelDiph{N,T}, state; t::Real=zero(T)) where {N,T}
    flux = recover_flux(model, state; t=t)
    q1 = model.ops1.H' * flux.phase1.faces
    q2 = model.ops2.H' * flux.phase2.faces
    idxγ = _interface_mask(model.cap1)

    t1 = zero(T)
    t2 = zero(T)
    bal = zero(T)
    for i in eachindex(q1)
        idxγ[i] || continue
        t1 += q1[i]
        t2 += q2[i]
        bal += q1[i] + q2[i]
    end

    return (phase1=q1, phase2=q2, total_phase1=t1, total_phase2=t2, balance=bal)
end

function mass_balance(
    model::DarcyModelMono{N,T},
    state;
    t::Real=zero(T),
    p_prev=nothing,
    dt=nothing,
) where {N,T}
    tt = convert(T, t)
    src = integrated_source(model; t=tt)
    qw = integrated_well_rate(model; t=tt)

    qbd = zero(T)
    per_side = Dict{Symbol,T}()
    for side in _all_sides(N)
        if get(model.bc_border.borders, side, nothing) isa Periodic
            continue
        end
        q = boundary_discharge(model, state, side; t=tt)
        per_side[side] = q
        qbd += q
    end

    storage_rate = zero(T)
    if !(p_prev === nothing || dt === nothing)
        dtt = convert(T, dt)
        dtt > zero(T) || throw(ArgumentError("dt must be positive when computing storage_rate"))
        pω, _ = _extract_mono_pressure_state(model, state)
        pω0, _ = _extract_mono_pressure_state(model, p_prev)
        Sω = _storage_values_mono(model.cap, model.storage, tt)
        @inbounds for i in eachindex(pω)
            Vi = model.cap.buf.V[i]
            if isfinite(Vi) && Vi > zero(T)
                storage_rate += Sω[i] * Vi * (pω[i] - pω0[i]) / dtt
            end
        end
    end

    residual = storage_rate + qbd - src - qw
    return (
        storage_rate=storage_rate,
        source_integral=src,
        well_rate=qw,
        boundary_discharge=qbd,
        residual=residual,
        imbalance=residual,
        per_side=per_side,
    )
end

function mass_balance(
    model::DarcyModelDiph{N,T},
    state;
    t::Real=zero(T),
    p_prev=nothing,
    dt=nothing,
) where {N,T}
    tt = convert(T, t)
    src = integrated_source(model; t=tt)
    qw = integrated_well_rate(model; t=tt)

    qbd1 = zero(T)
    qbd2 = zero(T)
    per_side_1 = Dict{Symbol,T}()
    per_side_2 = Dict{Symbol,T}()
    for side in _all_sides(N)
        if get(model.bc_border.borders, side, nothing) isa Periodic
            continue
        end
        q1 = boundary_discharge(model, state, side; t=tt, phase=1)
        q2 = boundary_discharge(model, state, side; t=tt, phase=2)
        per_side_1[side] = q1
        per_side_2[side] = q2
        qbd1 += q1
        qbd2 += q2
    end

    srate1 = zero(T)
    srate2 = zero(T)
    if !(p_prev === nothing || dt === nothing)
        dtt = convert(T, dt)
        dtt > zero(T) || throw(ArgumentError("dt must be positive when computing storage_rate"))
        p1ω, _, p2ω, _ = _extract_diph_pressure_state(model, state)
        p10, _, p20, _ = _extract_diph_pressure_state(model, p_prev)
        S1 = _storage_values_mono(model.cap1, model.storage1, tt)
        S2 = _storage_values_mono(model.cap2, model.storage2, tt)
        @inbounds for i in eachindex(p1ω)
            V1 = model.cap1.buf.V[i]
            if isfinite(V1) && V1 > zero(T)
                srate1 += S1[i] * V1 * (p1ω[i] - p10[i]) / dtt
            end
            V2 = model.cap2.buf.V[i]
            if isfinite(V2) && V2 > zero(T)
                srate2 += S2[i] * V2 * (p2ω[i] - p20[i]) / dtt
            end
        end
    end

    iq = interface_discharge(model, state; t=tt)
    res1 = srate1 + qbd1 + iq.total_phase1 - src.phase1 - qw.phase1
    res2 = srate2 + qbd2 + iq.total_phase2 - src.phase2 - qw.phase2

    total_storage = srate1 + srate2
    total_boundary = qbd1 + qbd2
    total_source = src.total
    total_well = qw.total
    residual = total_storage + total_boundary - total_source - total_well

    return (
        phase1=(storage_rate=srate1, source_integral=src.phase1, well_rate=qw.phase1, boundary_discharge=qbd1, interface_discharge=iq.total_phase1, residual=res1, imbalance=res1, per_side=per_side_1),
        phase2=(storage_rate=srate2, source_integral=src.phase2, well_rate=qw.phase2, boundary_discharge=qbd2, interface_discharge=iq.total_phase2, residual=res2, imbalance=res2, per_side=per_side_2),
        total=(storage_rate=total_storage, source_integral=total_source, well_rate=total_well, boundary_discharge=total_boundary, residual=residual, imbalance=residual),
        interface_balance=iq.balance,
    )
end

compute_mass_balance(args...; kwargs...) = mass_balance(args...; kwargs...)
