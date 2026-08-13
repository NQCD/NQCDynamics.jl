module NQCDistributions

using Distributions: Sampleable, Univariate
using Random: Random, AbstractRNG, SamplerTrivial, Xoshiro
using RingPolymerArrays: RingPolymerArray, eachbead
using UnitfulAtomic: austrip
using Distributions: Normal
using RingPolymerArrays: NormalModeTransformation, transform_from_normal_modes!
using LinearAlgebra: mul!, diagind
using ComponentArrays: ComponentVector

include("sampleable_components.jl")

include("nuclear/boltzmann.jl")
export VelocityBoltzmann

include("nuclear/harmonic_wigner.jl")
export VelocityHarmonicWigner
export MomentumHarmonicWigner
export PositionHarmonicWigner

include("nuclear/harmonic_ring_polymer.jl")
export PositionHarmonicRingPolymer

include("nuclear/free_ring_polymer.jl")
export PositionFreeRingPolymer

include("electronic/electronic.jl")
export Adiabatic, Diabatic
export PureState
export MixedState
export FermiDiracState

include("dynamical_distribution.jl")
export DynamicalDistribution

include("product_distribution.jl")

end
