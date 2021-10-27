using NonadiabaticMolecularDynamics
using Test
using RecursiveArrayTools
using OrdinaryDiffEq
using ComponentArrays

atoms = Atoms([:H])
model = NonadiabaticModels.Harmonic()
sim = Simulation(atoms, model)

positions = [randn(1, length(atoms)) for i=1:10]
velocities = [randn(1, length(atoms)) for i=1:10]
distribution = DynamicalDistribution(positions, velocities, (1, 1))

u = rand(distribution)
prob = DynamicsMethods.create_problem(u, (0.0, 1.0), sim)
@testset "OrderedSelection" begin 
    selector = Ensembles.OrderedSelection(distribution, 1:10)
    new_prob = selector(prob, 3, false)
    u = NonadiabaticDistributions.pick(distribution, 3)
    @test new_prob.u0 == ArrayPartition(u.v, u.r)
    selector = Ensembles.OrderedSelection(distribution, 6:10)
    new_prob = selector(prob, 3, false)
    u = NonadiabaticDistributions.pick(distribution, 8)
    @test new_prob.u0 == ArrayPartition(u.v, u.r)
end

@testset "RandomSelection" begin
    selector = Ensembles.RandomSelection(distribution)
    new_prob = selector(prob, 3, false)
end

@testset "SelectWithCallbacks" begin
    selector = Ensembles.RandomSelection(distribution)
    Ensembles.SelectWithCallbacks(selector, CallbackSet(), (:position,), 10)
end

@testset "SumReduction" begin
    prob = ODEProblem((u,p,t) -> 1.01u, 0.5, (0.0,1.0))
    prob_func(prob, i, repeat) = remake(prob,u0=rand())
    output_func(sol,i) = last(sol), false
    reduction = Ensembles.SumReduction()
    ensemble = EnsembleProblem(prob, prob_func=prob_func, output_func=output_func, reduction=reduction, u_init=reduction.u_init)
    sol = solve(ensemble, Tsit5(), trajectories=1000, batch_size=20)
    @test sol.u / 1000 ≈ 0.5 * exp(1.01) rtol=1e-1
end

@testset "run_trajectories" begin
    tspan = (0.0, 10.0)
    out = Ensembles.run_trajectories(sim, tspan, distribution;
        output=(:position), dt=1, trajectories=10)
end

