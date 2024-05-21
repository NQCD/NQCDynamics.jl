"""
Analysis functions common enough to be included in the main package.
"""
module Analysis
using Reexport: @reexport

# Diatomic analysis functions @alexsp32
include("diatomic.jl")
export Diatomic

# Allow re-processing of simulation data.
include("postprocess.jl")
export Postprocess

# Rebinding of quantise_diatomic under Analysis.
using NQCDynamics: InitialConditions
@reexport InitialConditions.QuantisedDiatomic

include("rigid_rotator.jl")
export RigidRotator

end
