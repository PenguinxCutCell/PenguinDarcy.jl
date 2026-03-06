function _interface_coupling_darcy(
    cap1::AssembledCapacity{N,T},
    cap2::AssembledCapacity{N,T},
    ic::AbstractDarcyInterfaceBC,
    t::T,
) where {N,T}
    nt = cap1.ntotal
    gp = zeros(T, nt)
    gq = zeros(T, nt)
    σ = zeros(T, nt)

    tol = sqrt(eps(T))
    @inbounds for i in 1:nt
        γ1 = cap1.buf.Γ[i]
        γ2 = cap2.buf.Γ[i]
        active1 = isfinite(γ1) && γ1 > zero(T)
        active2 = isfinite(γ2) && γ2 > zero(T)
        if active1 != active2
            throw(ArgumentError("cap1/cap2 interface masks differ at index $i"))
        end
        active1 || continue
        if abs(γ1 - γ2) > tol * max(one(T), abs(γ1), abs(γ2))
            throw(ArgumentError("cap1/cap2 interface measures differ at index $i"))
        end

        x = cap1.C_γ[i]
        if ic isa DarcyContinuity
            gp[i] = _eval_fun_or_const(ic.pressure_jump, x, t)
            gq[i] = _eval_fun_or_const(ic.flux_jump, x, t)
        elseif ic isa DarcyMembrane
            σ[i] = _eval_fun_or_const(ic.σ, x, t)
            gq[i] = _eval_fun_or_const(ic.flux_jump, x, t)
        else
            throw(ArgumentError("unsupported Darcy interface BC type $(typeof(ic))"))
        end
    end

    return gp, gq, σ
end

function _tensor_core_phase_diph(model::DarcyModelDiph{N,T}, phase::Int, t::T) where {N,T}
    cap = phase == 1 ? model.cap1 : model.cap2
    ops = phase == 1 ? model.ops1 : model.ops2

    Tf = tensor_flux_operator(model, ops, cap; t=t, phase=phase)
    Q = Tf * ops.Winv

    G = ops.G
    H = ops.H
    K = G' * Q * G
    C = G' * Q * H
    J = H' * Q * G
    L = H' * Q * H

    bfaces = face_bodyforce_values(model, ops, cap; t=t, phase=phase)
    ub = Tf * bfaces
    bωg = G' * ub
    bγg = H' * ub

    return K, C, J, L, bωg, bγg
end

function assemble_steady_diph!(sys::LinearSystem{T}, model::DarcyModelDiph{N,T}, t::T) where {N,T}
    nt = model.cap1.ntotal
    lay = model.layout.offsets
    nsys = maximum((last(lay.ω1), last(lay.γ1), last(lay.ω2), last(lay.γ2)))

    K1, C1, J1, L1, bωg1, bγg1 = _tensor_core_phase_diph(model, 1, t)
    K2, C2, J2, L2, bωg2, bγg2 = _tensor_core_phase_diph(model, 2, t)

    rhs1 = zeros(T, nt)
    rhs2 = zeros(T, nt)
    w1 = zeros(T, nt)
    w2 = zeros(T, nt)
    assemble_sources!(rhs1, model, model.cap1; t=t, phase=1)
    assemble_sources!(rhs2, model, model.cap2; t=t, phase=2)
    assemble_wells!(w1, model.wells1, model.cap1; t=t, phase=1)
    assemble_wells!(w2, model.wells2, model.cap2; t=t, phase=2)

    gp, gq, σ = _interface_coupling_darcy(model.cap1, model.cap2, model.bc_interface, t)

    maskγ = _interface_mask(model.cap1)
    mγ = T[maskγ[i] ? one(T) : zero(T) for i in 1:nt]
    Iγ = model.cap1.Γ
    Im = spdiagm(0 => mγ)
    Iσ = spdiagm(0 => σ)

    A11 = K1
    A12 = C1
    A13 = spzeros(T, nt, nt)
    A14 = spzeros(T, nt, nt)

    # γ1 row: pressure continuity (or membrane law)
    A21 = spzeros(T, nt, nt)
    A22 = Iγ
    A23 = spzeros(T, nt, nt)
    A24 = -Iγ

    A31 = spzeros(T, nt, nt)
    A32 = spzeros(T, nt, nt)
    A33 = K2
    A34 = C2

    # γ2 row: flux continuity
    A41 = -Im * J1
    A42 = -Im * L1
    A43 = -Im * J2
    A44 = -Im * L2

    b1 = rhs1 .+ w1 .+ bωg1
    b2 = Iγ * gp
    b3 = rhs2 .+ w2 .+ bωg2
    b4 = Im * gq .- Im * bγg1 .- Im * bγg2

    if model.bc_interface isa DarcyMembrane
        # p1 - p2 - σ*q1 = σ*flux_jump
        A21 = Im * (Iσ * J1)
        A22 = Iγ + Im * (Iσ * L1)
        b2 = Im * (Iσ * gq) .+ Im * (Iσ * bγg1)
    end

    A, b = if _is_canonical_diph_layout(lay, nt)
        ([A11 A12 A13 A14; A21 A22 A23 A24; A31 A32 A33 A34; A41 A42 A43 A44], vcat(b1, b2, b3, b4))
    else
        Awork = spzeros(T, nsys, nsys)
        bwork = zeros(T, nsys)

        _insert_block!(Awork, lay.ω1, lay.ω1, A11)
        _insert_block!(Awork, lay.ω1, lay.γ1, A12)
        _insert_block!(Awork, lay.ω1, lay.ω2, A13)
        _insert_block!(Awork, lay.ω1, lay.γ2, A14)

        _insert_block!(Awork, lay.γ1, lay.ω1, A21)
        _insert_block!(Awork, lay.γ1, lay.γ1, A22)
        _insert_block!(Awork, lay.γ1, lay.ω2, A23)
        _insert_block!(Awork, lay.γ1, lay.γ2, A24)

        _insert_block!(Awork, lay.ω2, lay.ω1, A31)
        _insert_block!(Awork, lay.ω2, lay.γ1, A32)
        _insert_block!(Awork, lay.ω2, lay.ω2, A33)
        _insert_block!(Awork, lay.ω2, lay.γ2, A34)

        _insert_block!(Awork, lay.γ2, lay.ω1, A41)
        _insert_block!(Awork, lay.γ2, lay.γ1, A42)
        _insert_block!(Awork, lay.γ2, lay.ω2, A43)
        _insert_block!(Awork, lay.γ2, lay.γ2, A44)

        _insert_vec!(bwork, lay.ω1, b1)
        _insert_vec!(bwork, lay.γ1, b2)
        _insert_vec!(bwork, lay.ω2, b3)
        _insert_vec!(bwork, lay.γ2, b4)
        (Awork, bwork)
    end

    sys.A = A
    sys.b = b
    length(sys.x) == nsys || (sys.x = zeros(T, nsys))
    sys.cache = nothing

    lay1 = UnknownLayout(nt, (ω=lay.ω1, γ=lay.γ1))
    lay2 = UnknownLayout(nt, (ω=lay.ω2, γ=lay.γ2))
    _apply_box_bc_darcy!(
        sys.A,
        sys.b,
        model.cap1,
        model.ops1,
        model.λ1,
        model.bc_border;
        ρ=model.ρ1,
        gravity=model.gravity,
        variable=model.variable,
        t=t,
        layout=lay1,
    )
    _apply_box_bc_darcy!(
        sys.A,
        sys.b,
        model.cap2,
        model.ops2,
        model.λ2,
        model.bc_border;
        ρ=model.ρ2,
        gravity=model.gravity,
        variable=model.variable,
        t=t,
        layout=lay2,
    )

    active_rows = _diph_row_activity(model.cap1, model.cap2, lay)
    sys.A, sys.b = _apply_row_identity_constraints!(sys.A, sys.b, active_rows)
    return sys
end

function _init_unsteady_state_diph(model::DarcyModelDiph{N,T}, u0) where {N,T}
    lay = model.layout.offsets
    nt = model.cap1.ntotal
    nsys = maximum((last(lay.ω1), last(lay.γ1), last(lay.ω2), last(lay.γ2)))
    u = zeros(T, nsys)
    if length(u0) == nsys
        u .= Vector{T}(u0)
    elseif length(u0) == 2 * nt
        u0v = Vector{T}(u0)
        u[lay.ω1] .= u0v[1:nt]
        u[lay.ω2] .= u0v[(nt + 1):(2 * nt)]
        # Use phase-local traces as a consistent default for CN/θ-schemes.
        u[lay.γ1] .= u[lay.ω1]
        u[lay.γ2] .= u[lay.ω2]
    else
        throw(DimensionMismatch("u0 length must be $(2 * nt) (ω1+ω2) or $nsys (full system)"))
    end
    return u
end

function assemble_unsteady_diph!(
    sys::LinearSystem{T},
    model::DarcyModelDiph{N,T},
    uⁿ,
    t::T,
    dt::T,
    scheme,
) where {N,T}
    θ = _theta_from_scheme(T, scheme)
    assemble_steady_diph!(sys, model, t + θ * dt)

    lay = model.layout.offsets
    nt = model.cap1.ntotal
    nsys = maximum((last(lay.ω1), last(lay.γ1), last(lay.ω2), last(lay.γ2)))

    ufull = if length(uⁿ) == nsys
        Vector{T}(uⁿ)
    elseif length(uⁿ) == 2 * nt
        v = zeros(T, nsys)
        u0v = Vector{T}(uⁿ)
        v[lay.ω1] .= u0v[1:nt]
        v[lay.ω2] .= u0v[(nt + 1):(2 * nt)]
        v
    else
        v = zeros(T, nsys)
        v[lay.ω1] .= Vector{T}(uⁿ[lay.ω1])
        v[lay.ω2] .= Vector{T}(uⁿ[lay.ω2])
        v
    end

    if θ != one(T)
        Aω1_prev = sys.A[lay.ω1, :]
        Aω2_prev = sys.A[lay.ω2, :]
        corr1 = Aω1_prev * ufull
        corr2 = Aω2_prev * ufull
        _scale_rows!(sys.A, lay.ω1, θ)
        _scale_rows!(sys.A, lay.ω2, θ)
        _insert_vec!(sys.b, lay.ω1, (-(one(T) - θ)) .* corr1)
        _insert_vec!(sys.b, lay.ω2, (-(one(T) - θ)) .* corr2)
    end

    S1 = _storage_values_mono(model.cap1, model.storage1, t + θ * dt)
    S2 = _storage_values_mono(model.cap2, model.storage2, t + θ * dt)
    M1 = (S1 .* model.cap1.buf.V) ./ dt
    M2 = (S2 .* model.cap2.buf.V) ./ dt
    rows = vcat(collect(lay.ω1), collect(lay.ω2))
    vals = vcat(M1, M2)
    sys.A = sys.A + sparse(rows, rows, vals, nsys, nsys)

    _insert_vec!(sys.b, lay.ω1, M1 .* Vector{T}(ufull[lay.ω1]))
    _insert_vec!(sys.b, lay.ω2, M2 .* Vector{T}(ufull[lay.ω2]))

    active_rows = _diph_row_activity(model.cap1, model.cap2, lay)
    sys.A, sys.b = _apply_row_identity_constraints!(sys.A, sys.b, active_rows)
    sys.cache = nothing
    return sys
end

function PenguinSolverCore.assemble!(sys::LinearSystem{T}, model::DarcyModelDiph{N,T}, t, dt) where {N,T}
    assemble_steady_diph!(sys, model, convert(T, t))
end
