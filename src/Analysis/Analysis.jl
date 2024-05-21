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
function quantise_diatomic(sim, v, r; args...)
    InitialConditions.QuantisedDiatomic.quantise_diatomic(sim, v, r; args...)
end
export quantise_diatomic

end
