"""
    InitialConditions
    
Functions and types for generating initial conditions for simulations.
"""
module InitialConditions


include("QuantisedDiatomic/QuantisedDiatomic.jl")
export QuantisedDiatomic

include("ConfigureAtomic.jl")
export ConfigureAtomic

include("ThermalMonteCarlo.jl")
export ThermalMonteCarlo

include("MetropolisHastings.jl")
export MetropolisHastings

end # module
