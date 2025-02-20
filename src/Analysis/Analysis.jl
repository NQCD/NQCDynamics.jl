"""
Analysis functions common enough to be included in the main package.
"""
module Analysis
using Reexport: @reexport

# Diatomic analysis functions
include("diatomic.jl")
export Diatomic

# Allow re-processing of simulation data.
include("postprocess.jl")
export Postprocess

# Rigid rotator energies
include("rigid_rotator.jl")
export RigidRotator

end
