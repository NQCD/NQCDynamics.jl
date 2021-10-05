
using Test, NonadiabaticMolecularDynamics.InitialConditions

ω = 3
β = 5
m = 10
mom = MomentumHarmonicWigner(ω, β)
pos = PositionHarmonicWigner(ω, β)
vel = VelocityHarmonicWigner(ω, β, m)

norm = tanh(β*ω/2) / π
@test norm ≈ 1 / (mom.σ * pos.σ * 2π)
@test norm ≈ 1 / (vel.σ * pos.σ * 2π) / m

vec = MomentumHarmonicWigner.([1, 2, 3, 4, 5], β)

d = DynamicalDistribution(vec, vec, (1, 1, 1))
@test_throws ErrorException rand(d)
d = DynamicalDistribution(vec, vec, (2, 5, 2))
rand(d)
