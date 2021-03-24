"""
    InitialConditions
    
Functions and types for generating initial conditions for simulations.
"""
module InitialConditions

using Reexport
using ..NonadiabaticMolecularDynamics

include("nuclear_distributions.jl")

include("MetropolisHastings.jl")
@reexport using .MetropolisHastings
include("QuantisedDiatomic.jl")
@reexport using .QuantisedDiatomic

end # module