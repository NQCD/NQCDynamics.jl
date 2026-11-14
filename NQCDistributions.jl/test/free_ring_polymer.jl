using Test
using NQCDynamics
using NQCDistributions
using Statistics: mean

ω = 3.4
β = 11.2
m = 1.5
n_beads = 400
s = PositionFreeRingPolymer{Float64}(β, m, (1,1,n_beads))
configs = [reshape(rand(s), (1,1,n_beads)) for _ in 1:10000]

sim = RingPolymerSimulation(Atoms(m), Free(), n_beads, temperature=1/β)

@testset "Energy expectation" begin
    estimate = Estimators.@estimate total_energy(sim, configs)
    @test 1 / 2β ≈ estimate rtol=1e-2
end

@testset "Radius of gyration expectation" begin
    T = 1/β
    analytic = sqrt(β / 12m) # https://www.cmt.york.ac.uk/cmd/mjg_thesis.pdf
    estimate = Estimators.@estimate radius_of_gyration(sim, configs)
    @test analytic ≈ estimate[1] rtol=1e-1
end

@testset "modified centre" begin
    centre = 57
    s = PositionFreeRingPolymer{Float64}(β, m, (1,1, n_beads); centre)
    @test mean(rand(s)) ≈ 57 atol=0.5
end

@testset "DynamicalDistribution" begin
    d = DynamicalDistribution(s, s, (1, 1, n_beads))
    rand(d)
    d[1]
end