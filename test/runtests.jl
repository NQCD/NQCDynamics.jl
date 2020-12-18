using Test
using NonadiabaticMolecularDynamics

const tests = [
    "model"
    "systems"
    "io/io"
    "dynamics"
    "ring_polymers"
    "dynamics/langevin"
    "dynamics/mdef"
    "monte_carlo"
]

for t in tests
    @testset "Test $t" begin
        include("$t.jl")
    end
end