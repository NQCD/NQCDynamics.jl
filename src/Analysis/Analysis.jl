"""
Analysis functions common enough to be included in the main package.
"""
module Analysis
using Reexport: @reexport

# Diatomic analysis functions @alexsp32
include("diatomic.jl")
export Diatomic

# Rebinding of quantise_diatomic under Analysis. 
using NQCDynamics: InitialConditions
@reexport using .InitialConditions: quantise_diatomic

end