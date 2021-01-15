using Test
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.InitialConditions
using Random
using Distributions

a = [rand(Float64, 3, 2) for i=1:10]
d = PhasespaceDistribution(a, a, (3,2))
@test eltype(a[1]) == eltype(d)
@test (3,2) == size(d)
@test rand(d) isa Vector{<:Matrix}

d = PhasespaceDistribution(a, Normal(), (3,2))
@test rand(d) isa Vector{<:Matrix}

d = PhasespaceDistribution(1, Normal(), (3,2))
@test rand(d) isa Vector{<:Matrix}