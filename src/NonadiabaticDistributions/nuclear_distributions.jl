
abstract type NuclearDistribution <: NonadiabaticDistribution end

abstract type PositionDistribution <: NuclearDistribution end
abstract type VelocityDistribution <: NuclearDistribution end

include("dynamical_distribution.jl")
export DynamicalDistribution

include("boltzmann_velocity.jl")
export BoltzmannVelocityDistribution

include("harmonic_wigner.jl")
export MomentumHarmonicWigner
export PositionHarmonicWigner
export VelocityHarmonicWigner

include("harmonic_ring_polymer.jl")
export HarmonicRingPolymer
