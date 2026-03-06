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

function _side_component_name(d::Int)
    if d == 1
        return :x
    elseif d == 2
        return :y
    elseif d == 3
        return :z
    end
    throw(ArgumentError("unsupported face direction d=$d"))
end

function boundary_discharge(
    model::DarcyModelMono{N,T},
    state,
    side::Symbol;
    t::Real=zero(T),
) where {N,T}
    tt = convert(T, t)
    pω, _ = _extract_mono_pressure_state(model, state)
    d, is_high, _ = side_info(side, N)
    side_bc = get(model.bc_border.borders, side, Neumann(0.0))
    side_bc isa Periodic && return zero(T)

    xyz_d = model.cap.xyz[d]
    length(xyz_d) >= 2 || throw(ArgumentError("need at least 2 grid nodes in each dimension"))
    Δd = abs(xyz_d[2] - xyz_d[1])
    δ = Δd / T(2)
    x_d = is_high ? xyz_d[end] : xyz_d[1]

    cap = model.cap
    LI = LinearIndices(cap.nnodes)
    out = zero(T)
    for I in each_boundary_cell(cap.nnodes, side)
        lin = LI[I]
        Aface = cap.buf.A[d][lin]
        if isfinite(Aface) && !iszero(Aface)
            Cω = cap.C_ω[lin]
            xface = SVector{N,T}(ntuple(k -> (k == d ? x_d : Cω[k]), N))
            λ_face = convert(T, eval_coeff(model.λ, xface, tt, lin))
            pP = pω[lin]

            qn = if side_bc isa Dirichlet
                pB = convert(T, eval_bc(side_bc.value, xface, tt))
                -λ_face * (pB - pP) / δ
            elseif side_bc isa Neumann
                convert(T, eval_bc(side_bc.value, xface, tt))
            elseif side_bc isa Robin
                α = convert(T, eval_bc(side_bc.α, xface, tt))
                β = convert(T, eval_bc(side_bc.β, xface, tt))
                g = convert(T, eval_bc(side_bc.value, xface, tt))
                den = α * δ - β * λ_face
                iszero(den) && throw(ArgumentError("degenerate Darcy Robin boundary at side `$side`: α*δ - β*λ = 0"))
                λ_face * (α * pP - g) / den
            else
                throw(ArgumentError("unsupported boundary type $(typeof(side_bc))"))
            end

            out += qn * Aface
        end
    end
    return out
end

function compute_mass_balance(
    model::DarcyModelMono{N,T},
    state;
    t::Real=zero(T),
) where {N,T}
    tt = convert(T, t)
    sω = _source_values_mono(model.cap, model.source, tt)

    source_integral = zero(T)
    @inbounds for i in eachindex(sω)
        Vi = model.cap.buf.V[i]
        if isfinite(Vi) && Vi > zero(T)
            source_integral += Vi * sω[i]
        end
    end

    per_side = Dict{Symbol,T}()
    total_boundary_discharge = zero(T)
    for side in _all_sides(N)
        if get(model.bc_border.borders, side, nothing) isa Periodic
            continue
        end
        q = boundary_discharge(model, state, side; t=tt)
        per_side[side] = q
        total_boundary_discharge += q
    end

    return (
        source_integral=source_integral,
        boundary_discharge=total_boundary_discharge,
        imbalance=total_boundary_discharge - source_integral,
        per_side=per_side,
    )
end
