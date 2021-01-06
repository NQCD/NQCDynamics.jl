using Test
using NonadiabaticMolecularDynamics
using Aqua

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

Aqua.test_ambiguities(NonadiabaticMolecularDynamics)
# Aqua.test_unbound_args(NonadiabaticMolecularDynamics)
Aqua.test_undefined_exports(NonadiabaticMolecularDynamics)
Aqua.test_stale_deps(NonadiabaticMolecularDynamics, ignore=[:Documenter])
Aqua.test_project_toml_formatting(NonadiabaticMolecularDynamics)
