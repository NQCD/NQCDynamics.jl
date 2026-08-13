using SafeTestsets, Test

@safetestset "SampleableComponents" begin include("sampleable_components.jl") end
@safetestset "DynamicalDistribution" begin include("dynamical_distribution.jl") end
@safetestset "ElectronicDistributions" begin include("electronic.jl") end

@testset "Nuclear distributions" begin
    @safetestset "HarmonicWigner" begin include("harmonic_wigner.jl") end
    @safetestset "PositionHarmonicRingPolymer" begin include("harmonic_ring_polymer.jl") end
    @safetestset "PositionFreeRingPolymer" begin include("free_ring_polymer.jl") end
    @safetestset "Boltzmann" begin include("boltzmann.jl") end
end
