"""
    InitialConditions
    
Functions and types for generating initial conditions for simulations.
"""
module InitialConditions

include("QuantisedAtomic.jl")
export QuantisedDiatomic

include("QuantisedDiatomic.jl")
export QuantisedDiatomic

include("ThermalMonteCarlo.jl")
export ThermalMonteCarlo

include("MetropolisHastings.jl")
export MetropolisHastings

end # module
