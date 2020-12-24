using Test
using NonadiabaticMolecularDynamics
using Aqua

Aqua.test_ambiguities(NonadiabaticMolecularDynamics)
Aqua.test_unbound_args(NonadiabaticMolecularDynamics)
Aqua.test_undefined_exports(NonadiabaticMolecularDynamics)
Aqua.test_stale_deps(NonadiabaticMolecularDynamics, ignore=[:DifferentialEquations])
Aqua.test_project_toml_formatting(NonadiabaticMolecularDynamics)

const tests = [
    "atoms"
    "cells"
    "calculators"
    "ring_polymers"
    "monte_carlo"
    "phasespace"
    "dynamics/langevin"
    "dynamics/mdef"
    "io/io"
    "model"
]

for t in tests
    @testset "Test $t" begin
        include("$t.jl")
    end
end