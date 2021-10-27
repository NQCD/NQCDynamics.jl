
module NonadiabaticDistributions

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

maxindex(dist::CombinedDistribution) = maxindex(dist.nuclear)

end # module
