using Test
using SafeTestsets

@time @safetestset "Atom Tests" begin include("atoms.jl") end
@time @safetestset "Cell Tests" begin include("cells.jl") end
@time @safetestset "Calculator Tests" begin include("calculators.jl") end
@time @safetestset "Simulation Tests" begin include("simulations.jl") end
@time @safetestset "Ring Polymer Tests" begin include("ring_polymers.jl") end
@time @safetestset "RingPolymerArrays Tests" begin include("ring_polymer_array.jl") end
@time @safetestset "Monte Carlo Tests" begin include("monte_carlo.jl") end
@time @safetestset "Dynamical Variables Tests" begin include("dynamical_variables.jl") end
include("dynamics/dynamics.jl")
@time @safetestset "IO Tests" begin include("io/io.jl") end
@time @safetestset "Model Tests" begin include("model.jl") end
@time @safetestset "Distribution Tests" begin include("nuclear_distributions.jl") end
