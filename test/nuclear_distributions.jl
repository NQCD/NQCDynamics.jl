using Test
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.InitialConditions
using Random

a = [rand(Float64, 3, 2) for i=1:10]
d = PositionDistribution(a)
@test rand(d) isa Matrix
blank = zero(a[1])
@test blank == zero(a[1])
rand!(d, blank)
@test blank != zero(a[1])

d = PhasespaceDistribution(a, a)
@test rand(d) isa Phasespace
