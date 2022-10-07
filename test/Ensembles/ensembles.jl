using NQCDynamics
using Test
using RecursiveArrayTools: ArrayPartition

atoms = Atoms([:H])
model = NQCModels.Harmonic()
sim = Simulation(atoms, model)

positions = [randn(1, length(atoms)) for i=1:10]
velocities = [randn(1, length(atoms)) for i=1:10]
distribution = DynamicalDistribution(positions, velocities, (1, 1))
u0 = rand(distribution)
tspan = (0.0, 10.0)

@testset "run_ensemble" for reduction in (
    AppendReduction(), MeanReduction(), SumReduction(), FileReduction("test.h5")
)
    out = run_ensemble(sim, tspan, distribution; reduction,
        output=(PositionOutput, VelocityOutput, HamiltonianOutput), dt=1, trajectories=10)
    if (reduction isa MeanReduction) || (reduction isa SumReduction)
        reduction isa SumReduction && @test out[:Time] == 0.0:10.0:100.0
        reduction isa MeanReduction && @test out[:Time] == 0.0:10.0
        @test out[:PositionOutput] isa Vector{<:AbstractMatrix}
        @test out[:VelocityOutput] isa Vector{<:AbstractMatrix}
        @test out[:HamiltonianOutput] isa Vector{<:Number}
    elseif (reduction isa AppendReduction)
        @test all(x -> x == 0.0:10.0, (traj[:Time] for traj in out))
        @test all(x -> x isa Vector{<:AbstractMatrix}, (traj[:PositionOutput] for traj in out))
        @test all(x -> x isa Vector{<:AbstractMatrix}, (traj[:VelocityOutput] for traj in out))
        @test all(x -> x isa Vector{<:Number}, (traj[:HamiltonianOutput] for traj in out))
    end
end

@testset "run_ensemble" begin
    out = run_ensemble(sim, tspan, distribution;
        output=Ensembles.OutputFinal(), dt=1, trajectories=10, reduction=AppendReduction())
    final = [o[:OutputFinal] for o in out]
    @test final isa Vector{<:ArrayPartition}
end

@testset "run_ensemble, reduction=$reduction" for reduction ∈ (SumReduction(), MeanReduction())
    out = run_ensemble(sim, tspan, distribution;
        output=Ensembles.OutputFinal(), dt=1, trajectories=10, reduction=reduction)
    @test out[:OutputFinal] isa ArrayPartition
end
