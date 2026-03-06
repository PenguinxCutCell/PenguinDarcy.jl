module PenguinDarcy

using LinearAlgebra
using SparseArrays
using StaticArrays

using CartesianGeometry: GeometricMoments, geometric_moments, nan
using CartesianGrids: CartesianGrid, grid1d
using CartesianOperators
using PenguinBCs
using PenguinSolverCore

export DarcyModelMono
export assemble_steady_mono!, assemble_unsteady_mono!
export solve_steady!, solve_unsteady!
export recover_velocity, recover_flux
export compute_mass_balance, boundary_discharge
export face_mobility_values

include("models.jl")
include("coeffs.jl")
include("assembly.jl")
include("solve.jl")
include("postprocess.jl")
include("diagnostics.jl")

end
