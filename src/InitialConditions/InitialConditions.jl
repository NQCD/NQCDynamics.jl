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
include("Diatomic/Diatomic.jl")
@reexport using .Diatomic

end # module