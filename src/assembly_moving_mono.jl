function _interface_nodes_scalar_field(
    coeff,
    tracker::HeightFunctionTracker{N,T},
    t::T,
) where {N,T}
    out = similar(tracker.xf)
    ia = tracker.axis

    if coeff isa Number
        fill!(out, convert(T, coeff))
        return out
    elseif !(coeff isa Function)
        throw(ArgumentError("coefficient must be scalar or callback"))
    end

    td = _transverse_dims(ia, N)
    if N == 1
        x = tracker.xf[]
        pt = SVector{1,T}(x)
        out[] = _eval_fun_or_const(coeff, pt, t)
        return out
    elseif N == 2
        snodes = collect(T, grid1d(tracker.grid, td[1]))
        for j in eachindex(out)
            xa = tracker.xf[j]
            coords = ia == 1 ? (xa, snodes[j]) : (snodes[j], xa)
            pt = SVector{2,T}(coords)
            out[j] = _eval_fun_or_const(coeff, pt, t)
        end
        return out
    elseif N == 3
        ynodes = collect(T, grid1d(tracker.grid, td[1]))
        znodes = collect(T, grid1d(tracker.grid, td[2]))
        for I in CartesianIndices(out)
            xa = tracker.xf[I]
            yy = ynodes[I[1]]
            zz = znodes[I[2]]
            coords = if ia == 1
                (xa, yy, zz)
            elseif ia == 2
                (yy, xa, zz)
            else
                (yy, zz, xa)
            end
            pt = SVector{3,T}(coords)
            out[I] = _eval_fun_or_const(coeff, pt, t)
        end
        return out
    end

    throw(ArgumentError("HeightFunctionTracker supports N <= 3"))
end

function _interp_tracker_field(
    tracker::HeightFunctionTracker{N,T},
    field,
    x::SVector{N,T},
) where {N,T}
    ia = tracker.axis
    if N == 1
        return field[]
    elseif N == 2
        td = _transverse_dims(ia, N)[1]
        snodes = collect(T, grid1d(tracker.grid, td))
        return _interp_pwlinear(snodes, vec(field), x[td]; periodic=tracker.periodic_transverse)
    elseif N == 3
        td = _transverse_dims(ia, N)
        ynodes = collect(T, grid1d(tracker.grid, td[1]))
        znodes = collect(T, grid1d(tracker.grid, td[2]))
        return _interp_bilinear(
            ynodes,
            znodes,
            field,
            x[td[1]],
            x[td[2]];
            periodic_x=tracker.periodic_transverse,
            periodic_y=tracker.periodic_transverse,
        )
    end
    throw(ArgumentError("HeightFunctionTracker supports N <= 3"))
end

function _mono_interface_pressure_callback(
    model::MovingDarcyModelMono{N,T},
    κfield,
) where {N,T}
    tracker = model.tracker
    return function (args...)
        if length(args) == N + 1
            x = SVector{N,T}(ntuple(d -> convert(T, args[d]), N))
            tt = convert(T, args[end])
        elseif length(args) == N
            x = SVector{N,T}(ntuple(d -> convert(T, args[d]), N))
            tt = zero(T)
        else
            throw(ArgumentError("interface pressure callback expects (x..., t) or (x...)"))
        end
        pext = _eval_fun_or_const(model.p_ext, x, tt)
        σv = _eval_fun_or_const(model.surface_tension, x, tt)
        κv = _interp_tracker_field(tracker, κfield, x)
        return pext + σv * κv
    end
end

function assemble_unsteady_moving!(
    model::MovingDarcyModelMono{N,T};
    t::Real=zero(T),
) where {N,T}
    tt = convert(T, t)
    cap = rebuild_capacities(model.tracker; geom_method=model.geom_method, bc=zero(T), two_phase=false)
    ops = DiffusionOps(cap; periodic=periodic_flags(model.bc_border, N))

    κ = interface_curvature(model.tracker)
    σ = _interface_nodes_scalar_field(model.surface_tension, model.tracker, tt)
    jump = σ .* κ
    jump_min = minimum(jump)
    jump_max = maximum(jump)

    bc_if = Robin(1.0, 0.0, _mono_interface_pressure_callback(model, κ))

    dmodel = DarcyModelMono(
        cap,
        ops,
        model.λ;
        source=model.source,
        storage=zero(T),
        ρ=model.ρ,
        gravity=model.gravity,
        wells=model.wells,
        bc_border=model.bc_border,
        bc_interface=bc_if,
        coeff_mode=model.coeff_mode,
        variable=model.variable,
    )

    return (darcy_model=dmodel, curvature=κ, capillary_jump=jump, capillary_jump_min=jump_min, capillary_jump_max=jump_max)
end
