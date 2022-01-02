
using Test
using NonadiabaticMolecularDynamics
using ComponentArrays: ComponentVector
using StatsBase: mean
using Distributions: Normal

model = Harmonic()
sim = Simulation(Atoms(1), model; temperature=5)
rp_sim = RingPolymerSimulation(Atoms(1), model, 10; temperature=5)

vector = [ComponentVector(r=fill(1, size(sim)), v=randn(size(sim))) for _=1:100]
rp_vector = [ComponentVector(r=fill(1, size(rp_sim)), v=randn(size(rp_sim))) for _=1:100]

dist = DynamicalDistribution(Normal(), [randn(size(sim)) for _=1:100], size(sim))

@testset "@estimate" begin
    avg = mean(Estimators.potential_energy.(sim, vector))
    est = Estimators.@estimate potential_energy(sim, vector)
    @test avg ≈ est

    Estimators.@estimate potential_energy(sim, dist)
end


@testset "potential_energy" begin
    function f(u)
        r = DynamicsUtils.get_positions(u)
        sum(model.m .* model.ω.^2 .* r.^2) / 2
    end
    u = rand(vector)
    @test Estimators.potential_energy(sim, u) ≈ f(u)
    u = rand(rp_vector)
    @test Estimators.potential_energy(rp_sim, u) ≈ f(u) / nbeads(rp_sim)
end

@testset "kinetic_energy" begin
    function f(u)
        v = DynamicsUtils.get_velocities(u)
        sum(v.^2) / 2
    end
    u = rand(vector)
    @test Estimators.kinetic_energy(sim, u) ≈ f(u)
end

