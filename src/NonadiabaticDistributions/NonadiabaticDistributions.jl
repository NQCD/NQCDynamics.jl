
"""
    NonadiabaticDistributions

Module containing types for representing distributions of combined nuclear and electronic configurations.
"""
module NonadiabaticDistributions

"""
    NonadiabaticDistribution

Top level type for all nonadiabatic distributions.

These distributions contain a both electronic and nuclear degrees of freedom.
"""
abstract type NonadiabaticDistribution end

include("nuclear_distributions.jl")

include("electronic_distributions.jl")
export Diabatic
export Adiabatic
export SingleState
export ElectronicPopulation

struct CombinedDistribution{N<:NuclearDistribution, E<:ElectronicDistribution} <: NonadiabaticDistribution
    nuclear::N
    electronic::E
end

Base.:*(N::NuclearDistribution, E::ElectronicDistribution) = CombinedDistribution(N, E)
Base.:*(E::ElectronicDistribution, N::NuclearDistribution) = CombinedDistribution(N, E)

Base.lastindex(dist::CombinedDistribution) = lastindex(dist.nuclear)

end # module
