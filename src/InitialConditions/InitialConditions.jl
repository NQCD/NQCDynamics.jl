"""
    InitialConditions
    
Functions and types for generating initial conditions for simulations.
"""
module InitialConditions

using Reexport
using ..NonadiabaticMolecularDynamics

include("distributions/dynamical_distribution.jl")
include("distributions/boltzmann_velocity.jl")

include("MetropolisHastings.jl")
export MetropolisHastings
include("QuantisedDiatomic.jl")
export QuantisedDiatomic

include("advancedmh_sampling.jl")

end # module