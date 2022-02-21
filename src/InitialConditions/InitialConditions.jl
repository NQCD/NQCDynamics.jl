"""
    InitialConditions
    
Functions and types for generating initial conditions for simulations.
"""
module InitialConditions


include("QuantisedAtomic.jl")
export QuantisedAtomic

include("QuantisedDiatomic/QuantisedDiatomic.jl")
export QuantisedDiatomic

include("ThermalMonteCarlo.jl")
export ThermalMonteCarlo

include("MetropolisHastings.jl")
export MetropolisHastings

end # module
