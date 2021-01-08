using Test
using NonadiabaticMolecularDynamics

const tests = [
    "atoms"
    "cells"
    "calculators"
    "ring_polymers"
    "monte_carlo"
    "phasespace"
    "dynamics/langevin"
    "dynamics/mdef"
    "dynamics/fssh"
    "dynamics/fermionic"
    "dynamics/ensembles"
    "io/io"
    "model"
    "nuclear_distributions"
]

for t in tests
    @testset "Test $t" begin
        include("$t.jl")
    end
end
