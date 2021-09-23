"""
    InitialConditions
    
Functions and types for generating initial conditions for simulations.
"""
module InitialConditions

include("distributions/dynamical_distribution.jl")
export DynamicalDistribution

include("distributions/boltzmann_velocity.jl")
export BoltzmannVelocityDistribution

include("distributions/harmonic_wigner.jl")
export MomentumHarmonicWigner
export PositionHarmonicWigner
export VelocityHarmonicWigner

include("QuantisedDiatomic.jl")
export QuantisedDiatomic

include("ThermalMonteCarlo.jl")
export ThermalMonteCarlo

include("MetropolisHastings.jl")
export MetropolisHastings

end # module
