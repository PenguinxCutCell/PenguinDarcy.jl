function _diph_pressure_jump_callback(
    model::MovingDarcyModelDiph{N,T},
    κfield,
) where {N,T}
    tracker = model.tracker
    base = model.bc_interface
    base isa DarcyContinuity || throw(ArgumentError("moving two-phase assembly currently supports DarcyContinuity only"))

    return function (args...)
        if length(args) == N + 1
            x = SVector{N,T}(ntuple(d -> convert(T, args[d]), N))
            tt = convert(T, args[end])
        elseif length(args) == N
            x = SVector{N,T}(ntuple(d -> convert(T, args[d]), N))
            tt = zero(T)
        else
            throw(ArgumentError("pressure-jump callback expects (x..., t) or (x...)"))
        end

        jbase = _eval_fun_or_const(base.pressure_jump, x, tt)
        σv = _eval_fun_or_const(model.surface_tension, x, tt)
        κv = _interp_tracker_field(tracker, κfield, x)
        return jbase + σv * κv
    end
end

function assemble_unsteady_moving!(
    model::MovingDarcyModelDiph{N,T};
    t::Real=zero(T),
) where {N,T}
    tt = convert(T, t)
    model.bc_interface isa DarcyContinuity || throw(ArgumentError("moving two-phase Darcy currently supports DarcyContinuity only"))

    cap1, cap2 = rebuild_capacities(model.tracker; geom_method=model.geom_method, bc=zero(T), two_phase=true)
    ops1 = DiffusionOps(cap1; periodic=periodic_flags(model.bc_border, N))
    ops2 = DiffusionOps(cap2; periodic=periodic_flags(model.bc_border, N))

    κ = interface_curvature(model.tracker)
    σ = _interface_nodes_scalar_field(model.surface_tension, model.tracker, tt)
    jump_cap = σ .* κ

    bc_if = DarcyContinuity(
        _diph_pressure_jump_callback(model, κ),
        0.0,
    )

    dmodel = DarcyModelDiph(
        cap1,
        ops1,
        model.λ1,
        cap2,
        ops2,
        model.λ2;
        source=(model.source1, model.source2),
        storage=(zero(T), zero(T)),
        ρ=(model.ρ1, model.ρ2),
        gravity=model.gravity,
        wells1=model.wells1,
        wells2=model.wells2,
        bc_border=model.bc_border,
        bc_interface=bc_if,
        coeff_mode=model.coeff_mode,
        variable=model.variable,
    )

    return (
        darcy_model=dmodel,
        curvature=κ,
        capillary_jump=jump_cap,
        capillary_jump_min=minimum(jump_cap),
        capillary_jump_max=maximum(jump_cap),
    )
end
