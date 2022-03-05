
using NQCDynamics
using Test
using RecursiveArrayTools
using DiffEqBase: CallbackSet

atoms = Atoms([:H])
model = NQCModels.Harmonic()
sim = Simulation(atoms, model)

positions = [randn(1, length(atoms)) for i=1:10]
velocities = [randn(1, length(atoms)) for i=1:10]
distribution = DynamicalDistribution(positions, velocities, (1, 1))

u = rand(distribution)
prob = DynamicsMethods.create_problem(u, (0.0, 1.0), sim)
@testset "OrderedSelection" begin 
    selector = Ensembles.OrderedSelection(distribution, 1:10)
    new_prob = selector(prob, 3, false)
    u = getindex(distribution, 3)
    @test new_prob.u0 == ArrayPartition(u.v, u.r)
    selector = Ensembles.OrderedSelection(distribution, 6:10)
    new_prob = selector(prob, 3, false)
    u = getindex(distribution, 8)
    @test new_prob.u0 == ArrayPartition(u.v, u.r)
end

@testset "RandomSelection" begin
    selector = Ensembles.RandomSelection(distribution)
    new_prob = selector(prob, 3, false)
end
