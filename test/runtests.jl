using Test
using NonadiabaticMolecularDynamics

const tests = [
    "atoms"
    "cells"
    "calculators"
    "simulations"
    "ring_polymers"
    "monte_carlo"
    "phasespace"
    "dynamics/dynamics"
    "io/io"
    "model"
    "nuclear_distributions"
]

for t in tests
    @testset "Test $t" begin
        include("$t.jl")
    end
end
