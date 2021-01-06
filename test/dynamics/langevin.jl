using Test
using NonadiabaticMolecularDynamics
using Unitful

@test Dynamics.Langevin{Float64}(1.0) isa Dynamics.Langevin
