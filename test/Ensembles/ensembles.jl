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

@testset "run_dynamics" for reduction in (
    AppendReduction(), MeanReduction(), SumReduction(), FileReduction("test.h5")
)
    out = run_dynamics(sim, tspan, distribution; reduction,
        output=(OutputPosition, OutputVelocity, OutputTotalEnergy), dt=1, trajectories=10)
    if (reduction isa MeanReduction) || (reduction isa SumReduction)
        reduction isa SumReduction && @test out[:Time] == 0.0:10.0:100.0
        reduction isa MeanReduction && @test out[:Time] == 0.0:10.0
        @test out[:OutputPosition] isa Vector{<:AbstractMatrix}
        @test out[:OutputVelocity] isa Vector{<:AbstractMatrix}
        @test out[:OutputTotalEnergy] isa Vector{<:Number}
    elseif (reduction isa AppendReduction)
        @test all(x -> x == 0.0:10.0, (traj[:Time] for traj in out))
        @test all(x -> x isa Vector{<:AbstractMatrix}, (traj[:OutputPosition] for traj in out))
        @test all(x -> x isa Vector{<:AbstractMatrix}, (traj[:OutputVelocity] for traj in out))
        @test all(x -> x isa Vector{<:Number}, (traj[:OutputTotalEnergy] for traj in out))
    end
end

@testset "run_dynamics" begin
    out = run_dynamics(sim, tspan, distribution;
        output=OutputFinal, dt=1, trajectories=10, reduction=AppendReduction())
    final = [o[:OutputFinal] for o in out]
    @test final isa Vector{<:ArrayPartition}
end

@testset "run_dynamics, reduction=$reduction" for reduction ∈ (SumReduction(), MeanReduction())
    out = run_dynamics(sim, tspan, distribution;
        output=OutputFinal, dt=1, trajectories=10, reduction=reduction)
    @test out[:OutputFinal] isa ArrayPartition
end
