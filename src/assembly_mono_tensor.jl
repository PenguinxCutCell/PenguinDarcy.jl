function _interface_mask(cap::AssembledCapacity{N,T}) where {N,T}
    m = Vector{Bool}(undef, cap.ntotal)
    @inbounds for i in eachindex(m)
        γ = cap.buf.Γ[i]
        m[i] = isfinite(γ) && γ > zero(T)
    end
    return m
end

function _interface_diagonals_mono(cap::AssembledCapacity{N,T}, ic::Union{Nothing,PenguinBCs.Robin}, t::T) where {N,T}
    α = zeros(T, cap.ntotal)
    β = zeros(T, cap.ntotal)
    g = zeros(T, cap.ntotal)
    ic === nothing && return α, β, g

    mask = _interface_mask(cap)
    @inbounds for i in eachindex(mask)
        mask[i] || continue
        x = cap.C_γ[i]
        α[i] = convert(T, eval_bc(ic.α, x, t))
        β[i] = convert(T, eval_bc(ic.β, x, t))
        g[i] = convert(T, eval_bc(ic.value, x, t))
    end
    return α, β, g
end

function _insert_block!(A::SparseMatrixCSC{T,Int}, rows::UnitRange{Int}, cols::UnitRange{Int}, B::SparseMatrixCSC{T,Int}) where {T}
    size(B, 1) == length(rows) || throw(DimensionMismatch("block rows do not match target range"))
    size(B, 2) == length(cols) || throw(DimensionMismatch("block cols do not match target range"))
    @inbounds for j in 1:size(B, 2)
        for p in nzrange(B, j)
            i = B.rowval[p]
            A[rows[i], cols[j]] = A[rows[i], cols[j]] + B.nzval[p]
        end
    end
    return A
end

function _insert_vec!(b::Vector{T}, rows::UnitRange{Int}, v::Vector{T}) where {T}
    length(v) == length(rows) || throw(DimensionMismatch("vector block length mismatch"))
    @inbounds for i in eachindex(v)
        b[rows[i]] += v[i]
    end
    return b
end

function _cell_activity_masks(cap::AssembledCapacity{N,T}) where {N,T}
    nt = cap.ntotal
    activeω = BitVector(undef, nt)
    activeγ = BitVector(undef, nt)
    LI = LinearIndices(cap.nnodes)
    for I in CartesianIndices(cap.nnodes)
        lin = LI[I]
        halo = any(d -> I[d] == cap.nnodes[d], 1:N)
        if halo
            activeω[lin] = false
            activeγ[lin] = false
            continue
        end
        v = cap.buf.V[lin]
        γ = cap.buf.Γ[lin]
        activeω[lin] = isfinite(v) && v > zero(T)
        activeγ[lin] = isfinite(γ) && γ > zero(T)
    end
    return activeω, activeγ
end

function _mono_row_activity(cap::AssembledCapacity{N,T}, lay) where {N,T}
    activeω, activeγ = _cell_activity_masks(cap)
    nsys = maximum((last(lay.ω), last(lay.γ)))
    active = falses(nsys)
    @inbounds for i in 1:cap.ntotal
        active[lay.ω[i]] = activeω[i]
        active[lay.γ[i]] = activeγ[i]
    end
    return active
end

function _diph_row_activity(cap1::AssembledCapacity{N,T}, cap2::AssembledCapacity{N,T}, lay) where {N,T}
    activeω1, activeγ1 = _cell_activity_masks(cap1)
    activeω2, activeγ2 = _cell_activity_masks(cap2)
    nsys = maximum((last(lay.ω1), last(lay.γ1), last(lay.ω2), last(lay.γ2)))
    active = falses(nsys)
    @inbounds for i in 1:cap1.ntotal
        active[lay.ω1[i]] = activeω1[i]
        active[lay.γ1[i]] = activeγ1[i]
        active[lay.ω2[i]] = activeω2[i]
        active[lay.γ2[i]] = activeγ2[i]
    end
    return active
end

function _apply_row_identity_constraints!(
    A::SparseMatrixCSC{T,Int},
    b::Vector{T},
    active_rows::BitVector,
) where {T}
    n = size(A, 1)
    size(A, 2) == n || throw(ArgumentError("row-identity constraints require square matrix"))
    length(b) == n || throw(ArgumentError("rhs length mismatch"))
    length(active_rows) == n || throw(ArgumentError("active row mask length mismatch"))

    p = Vector{T}(undef, n)
    @inbounds for i in 1:n
        ai = active_rows[i]
        p[i] = ai ? zero(T) : one(T)
        ai || (b[i] = zero(T))
    end

    @inbounds for j in 1:size(A, 2)
        aj = active_rows[j]
        for k in nzrange(A, j)
            if !(aj && active_rows[A.rowval[k]])
                A.nzval[k] = zero(T)
            end
        end
    end
    dropzeros!(A)

    Aout = A + spdiagm(0 => p)
    return Aout, b
end

function _scale_rows!(A::SparseMatrixCSC{T,Int}, rows::UnitRange{Int}, α::T) where {T}
    α == one(T) && return A
    r1 = first(rows)
    r2 = last(rows)
    @inbounds for j in 1:size(A, 2)
        for p in nzrange(A, j)
            i = A.rowval[p]
            if r1 <= i <= r2
                A.nzval[p] *= α
            end
        end
    end
    return A
end

function _is_canonical_mono_layout(lay, nt::Int)
    return lay.ω == (1:nt) && lay.γ == ((nt + 1):(2 * nt))
end

function _is_canonical_diph_layout(lay, nt::Int)
    return lay.ω1 == (1:nt) &&
           lay.γ1 == ((nt + 1):(2 * nt)) &&
           lay.ω2 == ((2 * nt + 1):(3 * nt)) &&
           lay.γ2 == ((3 * nt + 1):(4 * nt))
end

function _side_pairs(N::Int)
    if N == 1
        return ((:left, :right),)
    elseif N == 2
        return ((:left, :right), (:bottom, :top))
    elseif N == 3
        return ((:left, :right), (:bottom, :top), (:backward, :forward))
    end
    throw(ArgumentError("unsupported dimension N=$N; expected 1, 2, or 3"))
end

function _diag_add!(A::SparseMatrixCSC{T,Int}, row::Int, v::T) where {T}
    A[row, row] = A[row, row] + v
    return A
end

function _normal_mobility_at_face(λ, N::Int, d::Int, xface::SVector{M,T}, t::T, row_lin::Int) where {M,T}
    mob = _wrap_mobility(λ, N)
    return _mobility_entry(mob, d, d, xface, t, row_lin)
end

"""
Apply outer box BCs for Darcy pressure equation with flux semantics:
- Dirichlet: p = p_D
- Neumann:   u⋅n = q_n
- Robin:     α*p + β*(u⋅n) = g
where u = -Λ(∇p - ρg). In v0.2 we use the normal mobility component for the
boundary elimination formula.
"""
function _apply_box_bc_darcy!(
    A::SparseMatrixCSC{T,Int},
    b::Vector{T},
    cap::AssembledCapacity{N,T},
    ops::DiffusionOps{N,T},
    λ,
    bc_border::BorderConditions;
    ρ=one(T),
    gravity=ntuple(_ -> zero(T), N),
    variable::Symbol=:pressure,
    t::T=zero(T),
    layout::UnknownLayout=layout_mono(cap.ntotal),
) where {N,T}
    length(b) >= last(layout.offsets.ω) || throw(ArgumentError("rhs vector does not contain ω block"))
    ops.nnodes == cap.nnodes || throw(ArgumentError("ops and cap have inconsistent nnodes"))
    validate_borderconditions!(bc_border, N)

    pairs = _side_pairs(N)
    LI = LinearIndices(cap.nnodes)

    for pair in pairs
        for side in pair
            side_bc = get(bc_border.borders, side, Neumann(0.0))
            side_bc isa Periodic && continue

            d, is_high, normal_sign = side_info(side, N)
            xyz_d = cap.xyz[d]
            length(xyz_d) >= 2 || throw(ArgumentError("need at least 2 grid nodes in each dimension"))
            Δd = abs(xyz_d[2] - xyz_d[1])
            δ = Δd / T(2)
            x_d = is_high ? xyz_d[end] : xyz_d[1]

            for I in each_boundary_cell(cap.nnodes, side)
                row_lin = LI[I]
                row = layout.offsets.ω[row_lin]

                Aface = cap.buf.A[d][row_lin]
                if !isfinite(Aface) || iszero(Aface)
                    continue
                end

                Cω = cap.C_ω[row_lin]
                xface = SVector{N,T}(ntuple(k -> (k == d ? x_d : Cω[k]), N))
                λnn = _normal_mobility_at_face(λ, N, d, xface, t, row_lin)
                a = λnn * Aface / δ
                bn = if variable === :head
                    zero(T)
                else
                    ρv = _eval_density(ρ, xface, t)
                    gv = _eval_gravity(gravity, xface, t)
                    convert(T, normal_sign) * (ρv * gv[d])
                end

                if side_bc isa Dirichlet
                    p_b = convert(T, eval_bc(side_bc.value, xface, t))
                    _diag_add!(A, row, a)
                    p_eff = is_high ? (p_b - δ * bn) : p_b
                    b[row] += a * p_eff
                elseif side_bc isa Neumann
                    # Neumann value is Darcy normal flux q_n = u⋅n.
                    qn = convert(T, eval_bc(side_bc.value, xface, t))
                    iszero(qn) && continue
                    b[row] += -qn * Aface
                elseif side_bc isa Robin
                    α = convert(T, eval_bc(side_bc.α, xface, t))
                    β = convert(T, eval_bc(side_bc.β, xface, t))
                    g = convert(T, eval_bc(side_bc.value, xface, t))
                    den = α * δ - β * λnn

                    iszero(den) && throw(ArgumentError("degenerate Darcy Robin boundary at side `$side`: α*δ - β*λnn = 0"))

                    _diag_add!(A, row, a * (α * δ) / den)
                    g_eff = is_high ? (g - α * δ * bn) : g
                    b[row] += a * (δ * g_eff) / den
                else
                    throw(ArgumentError("unsupported boundary type $(typeof(side_bc))"))
                end
            end
        end
    end

    return A, b
end

function _theta_from_scheme(::Type{T}, scheme) where {T}
    if scheme isa Symbol
        if scheme === :BE
            return one(T)
        elseif scheme === :CN
            return convert(T, 0.5)
        end
        throw(ArgumentError("unknown scheme `$scheme`; expected :BE or :CN"))
    elseif scheme isa Real
        return convert(T, scheme)
    end
    throw(ArgumentError("scheme must be a Symbol (:BE/:CN) or a numeric theta"))
end

function _init_unsteady_state_mono(model::DarcyModelMono{N,T}, u0) where {N,T}
    lay = model.layout.offsets
    nt = model.cap.ntotal
    nsys = maximum((last(lay.ω), last(lay.γ)))
    u = zeros(T, nsys)
    if length(u0) == nsys
        u .= Vector{T}(u0)
    elseif length(u0) == nt
        uω = Vector{T}(u0)
        u[lay.ω] .= uω
        # Use a consistent first guess for interface traces when only ω is provided.
        u[lay.γ] .= uω
    else
        throw(DimensionMismatch("u0 length must be $nt (ω block) or $nsys (full system)"))
    end
    return u
end

function _tensor_core_mono(
    model::DarcyModelMono{N,T},
    t::T,
) where {N,T}
    cap = model.cap
    ops = model.ops
    Tf = tensor_flux_operator(model, ops, cap; t=t, phase=1)
    Q = Tf * ops.Winv

    G = ops.G
    H = ops.H
    K = G' * Q * G
    C = G' * Q * H
    J = H' * Q * G
    L = H' * Q * H

    bfaces = face_bodyforce_values(model, ops, cap; t=t)
    ub = Tf * bfaces
    bωg = G' * ub
    bγg = H' * ub

    return K, C, J, L, bωg, bγg
end

function assemble_steady_mono!(sys::LinearSystem{T}, model::DarcyModelMono{N,T}, t::T) where {N,T}
    nt = model.cap.ntotal
    lay = model.layout.offsets
    nsys = maximum((last(lay.ω), last(lay.γ)))

    K, C, J, L, bωg, bγg = _tensor_core_mono(model, t)
    α, β, gγ = _interface_diagonals_mono(model.cap, model.bc_interface, t)

    rhs_source = zeros(T, nt)
    rhs_well = zeros(T, nt)
    assemble_sources!(rhs_source, model, model.cap; t=t, phase=1)
    assemble_wells!(rhs_well, model.wells, model.cap; t=t, phase=1)

    Iβ = spdiagm(0 => β)
    Iα = spdiagm(0 => α)
    Iγ = model.cap.Γ

    A11 = K
    A12 = C
    # Darcy interface Robin uses α*p + β*(u⋅n)=g and u⋅n = -q + q_body.
    A21 = -Iβ * J
    A22 = -Iβ * L + Iα * Iγ
    if model.bc_interface === nothing
        A12 = spzeros(T, nt, nt)
        A21 = spzeros(T, nt, nt)
        A22 = spdiagm(0 => ones(T, nt))
    end
    b1 = rhs_source .+ rhs_well .+ bωg
    b2 = Iγ * gγ + Iβ * bγg

    A, b = if _is_canonical_mono_layout(lay, nt)
        ([A11 A12; A21 A22], vcat(b1, b2))
    else
        Awork = spzeros(T, nsys, nsys)
        bwork = zeros(T, nsys)
        _insert_block!(Awork, lay.ω, lay.ω, A11)
        _insert_block!(Awork, lay.ω, lay.γ, A12)
        _insert_block!(Awork, lay.γ, lay.ω, A21)
        _insert_block!(Awork, lay.γ, lay.γ, A22)
        _insert_vec!(bwork, lay.ω, b1)
        _insert_vec!(bwork, lay.γ, b2)
        (Awork, bwork)
    end

    sys.A = A
    sys.b = b
    length(sys.x) == nsys || (sys.x = zeros(T, nsys))
    sys.cache = nothing

    _apply_box_bc_darcy!(
        sys.A,
        sys.b,
        model.cap,
        model.ops,
        model.λ,
        model.bc_border;
        ρ=model.ρ,
        gravity=model.gravity,
        variable=model.variable,
        t=t,
        layout=model.layout,
    )
    active_rows = _mono_row_activity(model.cap, lay)
    sys.A, sys.b = _apply_row_identity_constraints!(sys.A, sys.b, active_rows)
    return sys
end

function assemble_unsteady_mono!(
    sys::LinearSystem{T},
    model::DarcyModelMono{N,T},
    uⁿ,
    t::T,
    dt::T,
    scheme,
) where {N,T}
    θ = _theta_from_scheme(T, scheme)
    assemble_steady_mono!(sys, model, t + θ * dt)

    lay = model.layout.offsets
    nt = model.cap.ntotal
    nsys = maximum((last(lay.ω), last(lay.γ)))

    ufull = if length(uⁿ) == nsys
        Vector{T}(uⁿ)
    elseif length(uⁿ) == nt
        v = zeros(T, nsys)
        v[lay.ω] .= Vector{T}(uⁿ)
        v
    else
        v = zeros(T, nsys)
        v[lay.ω] .= Vector{T}(uⁿ[lay.ω])
        v
    end

    if θ != one(T)
        Aω_prev = sys.A[lay.ω, :]
        corr = Aω_prev * ufull
        _scale_rows!(sys.A, lay.ω, θ)
        _insert_vec!(sys.b, lay.ω, (-(one(T) - θ)) .* corr)
    end

    Sω = _storage_values_mono(model.cap, model.storage, t + θ * dt)
    M = (Sω .* model.cap.buf.V) ./ dt
    sys.A = sys.A + sparse(lay.ω, lay.ω, M, nsys, nsys)

    uω = Vector{T}(ufull[lay.ω])
    _insert_vec!(sys.b, lay.ω, M .* uω)
    active_rows = _mono_row_activity(model.cap, lay)
    sys.A, sys.b = _apply_row_identity_constraints!(sys.A, sys.b, active_rows)
    sys.cache = nothing
    return sys
end

function PenguinSolverCore.assemble!(sys::LinearSystem{T}, model::DarcyModelMono{N,T}, t, dt) where {N,T}
    assemble_steady_mono!(sys, model, convert(T, t))
end
