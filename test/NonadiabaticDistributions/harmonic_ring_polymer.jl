using Test
using NQCDynamics

ω = 3.4
β = 11.2
m = 1
n_beads = 400
s = HarmonicRingPolymer(ω, β, m, n_beads)
configs = (reshape(rand(s), (1,1,n_beads)) for _ in 1:5000)

sim = RingPolymerSimulation(Atoms(1), Harmonic(m=m, ω=ω), n_beads, temperature=1/β)

@testset "Energy expectation" begin
    quantum = ω/2 * coth(ω*β/2)
    estimate = Estimators.@estimate total_energy(sim, configs)
    @test quantum ≈ estimate rtol=1e-2
end

@testset "Radius of gyration expectation" begin
    T = 1/β
    quantum = sqrt(T / ω^2 * (ω / (2T) * coth(ω/2T) - 1))
    estimate = Estimators.@estimate radius_of_gyration(sim, configs)
    @test quantum ≈ estimate[1] rtol=1e-2
end

@testset "select_item" begin
    atoms = 10
    beads = 4
    s = HarmonicRingPolymer.(1:atoms, 1, 1, beads)

    r = NonadiabaticDistributions.select_item(s, 0, (3, atoms, beads))
    @test r isa Array
    @test_throws ErrorException("Sample size does not match distribution.") NonadiabaticDistributions.select_item(s, 0, (3, 1, beads))
    @test_throws ErrorException("Sample size does not match distribution.") NonadiabaticDistributions.select_item(s, 0, (3, 1, 1))
    @test_throws ErrorException("`size` (3, 1) must have 3 dimensions.") NonadiabaticDistributions.select_item(s, 0, (3, 1))
end

@testset "modified centre" begin
    centre = 57
    s = HarmonicRingPolymer(ω, β, m, n_beads; centre)
    @test mean(rand(s)) ≈ 57 atol=0.1
end
