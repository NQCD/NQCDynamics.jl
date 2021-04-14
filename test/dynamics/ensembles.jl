using NonadiabaticMolecularDynamics
using Test

atoms = Atoms{Float64}([:H])
model = Models.Harmonic()
sim = Simulation(atoms, model, Dynamics.Classical(); DoFs=1)

positions = [randn(1, length(atoms)) for i=1:10]
velocities = [randn(1, length(atoms)) for i=1:10]
distribution = InitialConditions.DynamicalDistribution(positions, velocities, (1, 1))

solution = Dynamics.run_ensemble(distribution, (0.0, 1e3), sim; trajectories=10, dt=1)

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