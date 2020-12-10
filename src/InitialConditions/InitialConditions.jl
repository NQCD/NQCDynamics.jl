"""
    InitialConditions
    
Functions and types for generating initial conditions for simulations.
"""
module InitialConditions

using Reexport

include("MetropolisHastings.jl")
@reexport using .MetropolisHastings

end # module