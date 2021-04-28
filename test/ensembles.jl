using NonadiabaticMolecularDynamics
using Test
using RecursiveArrayTools
using OrdinaryDiffEq

atoms = Atoms([:H])
model = Models.Harmonic()
sim = Simulation(atoms, model, Dynamics.Classical(); DoFs=1)

positions = [randn(1, length(atoms)) for i=1:10]
velocities = [randn(1, length(atoms)) for i=1:10]
distribution = InitialConditions.DynamicalDistribution(positions, velocities, (1, 1))

# solution = Dynamics.run_ensemble(distribution, (0.0, 1e3), sim; trajectories=10, dt=1)

# @testset "run_dissociation_ensemble" begin
#     atoms = Atoms([:H, :H])
#     model = Models.Free()
#     sim = Simulation(atoms, model; DoFs=1)
#     positions = [zeros(1, length(atoms)) for i=1:100]
#     velocities = [randn(1, length(atoms)) for i=1:100]
#     distribution = InitialConditions.DynamicalDistribution(positions, velocities, (1, 2))

#     solution = Dynamics.run_dissociation_ensemble(distribution, (0.0, 1e2), sim; trajectories=100, dt=1)
#     @test solution.u < 1.0
# end

prob = Dynamics.create_problem(ArrayPartition(rand(distribution)...), (0.0, 1.0), sim)
@testset "OrderedSelection" begin 
    selector = Ensembles.OrderedSelection(distribution)
    new_prob = selector(prob, 3, false)
    @test new_prob.u0 == ArrayPartition(InitialConditions.pick(distribution, 3)...)
end

@testset "RandomSelection" begin
    selector = Ensembles.RandomSelection(distribution)
    new_prob = selector(prob, 3, false)
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

@testset "OutputFinal" begin
    output = Ensembles.OutputFinal()
    @test output(1:10, 1) == (10, false)
end

@testset "OutputDissociation" begin
    output = Ensembles.OutputDissociation(1.0, [1, 2])
    r = [1 0; 0 0; 0 0]
    v = zeros(3, 2)
    u = ArrayPartition(v, r)
    sol = [u]
    @test output(sol, 1) == (0, false)
    r = [2 0; 0 0; 0 0]
    u = ArrayPartition(v, r)
    sol = [u]
    @test output(sol, 1) == (1, false)
end

@testset "OutputDiabaticPopulation" begin
    sim = Simulation{FSSH}(atoms, Models.TullyModelTwo(); DoFs=1)
    output = Ensembles.OutputDiabaticPopulation(sim)
end

