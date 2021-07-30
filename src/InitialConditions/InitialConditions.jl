"""
    InitialConditions
    
Functions and types for generating initial conditions for simulations.
"""
module InitialConditions

using Reexport
using ..NonadiabaticMolecularDynamics

include("nuclear_distributions.jl")

include("MetropolisHastings.jl")
export MetropolisHastings
include("QuantisedDiatomic.jl")
export QuantisedDiatomic

include("advancedmh_sampling.jl")

end # module