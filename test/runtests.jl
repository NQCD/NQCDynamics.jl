using Test
using SafeTestsets
using Documenter
using NonadiabaticMolecularDynamics

doctestsetup = quote
    using NonadiabaticMolecularDynamics
    using Plots
    using Symbolics
end 
DocMeta.setdocmeta!(NonadiabaticMolecularDynamics, :DocTestSetup, doctestsetup; recursive=true)
doctest(NonadiabaticMolecularDynamics)

@time @safetestset "Calculator Tests" begin include("calculators.jl") end
@time @safetestset "Simulation Tests" begin include("simulations.jl") end
@time @safetestset "Ring Polymer Tests" begin include("ring_polymers.jl") end
@time @safetestset "RingPolymerArrays Tests" begin include("ring_polymer_array.jl") end
@time @safetestset "Monte Carlo Tests" begin include("monte_carlo.jl") end
@time @safetestset "AdvancedMH Sampling Tests" begin include("advancedmh_sampling.jl") end
@time @safetestset "Dynamical Variables Tests" begin include("dynamical_variables.jl") end
include("dynamics/dynamics.jl")
@time @safetestset "IO Tests" begin include("io/io.jl") end
@time @safetestset "Distribution Tests" begin include("nuclear_distributions.jl") end
@time @safetestset "Ensemble Tests" begin include("ensembles.jl") end
