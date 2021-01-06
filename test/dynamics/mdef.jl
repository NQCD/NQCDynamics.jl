using Test
using NonadiabaticMolecularDynamics
using Unitful
using DifferentialEquations

@test Dynamics.MDEF{Float64}(10, 3) isa Dynamics.MDEF
atoms = Atoms{Float64}([:H, :C])
mdef = Dynamics.MDEF{Float64}(length(atoms), 1)
sim = Simulation(atoms, Models.FrictionHarmonic(), mdef; temperature=10u"K", DoFs=1)

R = zeros(sim.DoFs, length(sim.atoms)) 
P = zeros(sim.DoFs, length(sim.atoms)) 
u = Phasespace(R, P)
du = zero(u)

n = sim.DoFs*length(sim.atoms)*2
Dynamics.set_force!(du, u, sim)

@testset "random_force!" begin
    # Test that only the bottom left of the matrix is filled
    blank = zeros(n, n)
    Dynamics.random_force!(blank, u, sim, 1.0)
    @test all(blank[1:end, 1:n÷2] .== 0)
    @test all(blank[1:n÷2, 1:end] .== 0)
    @test all(blank[n÷2+1:end, n÷2+1:end] .!= 0)
end

prob = SDEProblem(Dynamics.motion!, Dynamics.random_force!, u, (0.0, 1.0), sim; noise_rate_prototype=zeros(n, n))
sol = Dynamics.run_trajectory(u, (0.0, 1.0), sim)