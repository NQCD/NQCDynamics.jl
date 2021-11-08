using NonadiabaticMolecularDynamics
using Test
using RecursiveArrayTools: ArrayPartition
using ComponentArrays: ComponentVector

atoms = Atoms([:H])
model = NonadiabaticModels.Harmonic()
sim = Simulation(atoms, model)

positions = [randn(1, length(atoms)) for i=1:10]
velocities = [randn(1, length(atoms)) for i=1:10]
distribution = DynamicalDistribution(positions, velocities, (1, 1))
tspan = (0.0, 10.0)

@testset "run_trajectories" begin
    out = Ensembles.run_trajectories(sim, tspan, distribution;
        output=(:position), dt=1, trajectories=10)
    @test length(out) == 10
    @test out[1].position isa Vector{<:Matrix}
    @test out[1].t isa Vector{<:Real}
end

@testset "run_ensemble" begin
    out = Ensembles.run_ensemble(sim, tspan, distribution;
        output=Ensembles.OutputFinal(), dt=1, trajectories=10, reduction=:append)
    @test out isa Vector{<:ArrayPartition}
end

@testset "run_ensemble, reduction=$reduction" for reduction ∈ (:sum, :mean)
    out = Ensembles.run_ensemble(sim, tspan, distribution;
        output=Ensembles.OutputFinal(), dt=1, trajectories=10, reduction=reduction)
    @test out isa ComponentVector
end
