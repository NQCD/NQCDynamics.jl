using Test
using SafeTestsets

@safetestset "Atom Tests" begin include("atoms.jl") end
@safetestset "Cell Tests" begin include("cells.jl") end
@safetestset "Calculator Tests" begin include("calculators.jl") end
@safetestset "Simulation Tests" begin include("simulations.jl") end
@safetestset "Ring Polymer Tests" begin include("ring_polymers.jl") end
@safetestset "Monte Carlo Tests" begin include("monte_carlo.jl") end
@safetestset "Quantised Diatomic Tests" begin include("quantised_diatomic.jl") end
@safetestset "Dynamical Variables Tests" begin include("dynamical_variables.jl") end
@safetestset "Dynamics Tests" begin include("dynamics/dynamics.jl") end
@safetestset "IO Tests" begin include("io/io.jl") end
@safetestset "Model Tests" begin include("model.jl") end
@safetestset "Distribution Tests" begin include("nuclear_distributions.jl") end
