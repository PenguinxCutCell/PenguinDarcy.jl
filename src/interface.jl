abstract type AbstractDarcyInterfaceBC end

struct DarcyContinuity <: AbstractDarcyInterfaceBC
    pressure_jump
    flux_jump
end

DarcyContinuity() = DarcyContinuity(0.0, 0.0)

struct DarcyMembrane <: AbstractDarcyInterfaceBC
    σ
    flux_jump
end

DarcyMembrane(σ) = DarcyMembrane(σ, 0.0)
