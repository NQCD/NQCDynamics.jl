using NQCDynamics
using Test
using RecursiveArrayTools: ArrayPartition
using ComponentArrays: ComponentVector

atoms = Atoms([:H])
model = NQCModels.Harmonic()
sim = Simulation(atoms, model)

positions = [randn(1, length(atoms)) for i=1:10]
velocities = [randn(1, length(atoms)) for i=1:10]
distribution = DynamicalDistribution(positions, velocities, (1, 1))
u0 = rand(distribution)
tspan = (0.0, 10.0)

@testset "run_ensemble" for reduction in (:append, :sum, :mean)
    out = run_ensemble(sim, tspan, distribution; reduction,
        output=(:position, :velocity, :hamiltonian), dt=1, trajectories=10)
    if (reduction === :sum) || (reduction === :mean)
        @test out.t == 0.0:10.0
        @test out.position isa Vector{<:AbstractMatrix}
        @test out.velocity isa Vector{<:AbstractMatrix}
        @test out.hamiltonian isa Vector{<:Number}
    else
        @test all(x -> x == 0.0:10.0, (traj.t for traj in out))
        @test all(x -> x isa Vector{<:AbstractMatrix}, (traj.position for traj in out))
        @test all(x -> x isa Vector{<:AbstractMatrix}, (traj.velocity for traj in out))
        @test all(x -> x isa Vector{<:Number}, (traj.hamiltonian for traj in out))
    end
end

@testset "run_ensemble" begin
    out = run_ensemble(sim, tspan, distribution;
        output=Ensembles.OutputFinal(), dt=1, trajectories=10, reduction=:append)
    @test out isa Vector{<:ArrayPartition}
end

@testset "run_ensemble, reduction=$reduction" for reduction ∈ (:sum, :mean)
    out = run_ensemble(sim, tspan, distribution;
        output=Ensembles.OutputFinal(), dt=1, trajectories=10, reduction=reduction, u_init=u0)
    @test out isa ComponentVector
end
