module PenguinDarcy

using LinearAlgebra
using SparseArrays
using StaticArrays

using CartesianGeometry: GeometricMoments, geometric_moments, nan
using CartesianGrids: CartesianGrid, grid1d, meshsize
using CartesianOperators
using GlobalHeightFunctions
using PenguinBCs
using PenguinSolverCore

export DarcyModelMono
export DarcyModelDiph
export DarcyCoupledModelMono
export DarcyContinuity, DarcyMembrane
export PointWell, CellWell
export assemble_steady_mono!, assemble_unsteady_mono!
export assemble_steady_diph!, assemble_unsteady_diph!
export solve_steady!, solve_unsteady!
export AbstractDarcyTracker, HeightFunctionTracker
export MovingDarcyModelMono, MovingDarcyModelDiph
export assemble_unsteady_moving!, solve_unsteady_moving!
export interface_normal_velocity, update_interface!
export interface_mass_residual, interface_iteration_history
export recover_velocity, recover_flux
export mass_balance, compute_mass_balance, boundary_discharge, interface_discharge
export face_mobility_values
export integrated_source, integrated_well_rate

include("models.jl")
include("interface.jl")
include("models_diph.jl")
include("coupling.jl")
include("wells.jl")
include("gravity.jl")
include("tensors.jl")
include("coeffs.jl")
include("assembly_mono_tensor.jl")
include("assembly_diph.jl")
include("solve.jl")
include("postprocess_tensor.jl")
include("diagnostics.jl")
include("moving_types.jl")
include("tracker.jl")
include("tracker_height_function.jl")
include("kinematics.jl")
include("assembly_moving_mono.jl")
include("assembly_moving_diph.jl")
include("diagnostics_moving.jl")
include("solve_moving.jl")

end
