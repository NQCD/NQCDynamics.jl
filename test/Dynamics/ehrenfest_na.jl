using Test
using NQCDynamics
using Statistics: var

kT = 9.5e-4
M = 30 # number of bath states
Γ = 6.4e-3
W = 6Γ / 2 # bandwidth  parameter

basemodel = MiaoSubotnik(;Γ)
bath = TrapezoidalRule(M, -W, W)
model = AndersonHolstein(basemodel, bath)
atoms = Atoms(2000)
r = randn(1,1)
v = randn(1,1)
n_electrons = M ÷ 2

sim = Simulation{EhrenfestNA}(atoms, model)
u = DynamicsVariables(sim, v, r)

@testset "create_problem" begin
    sim = Simulation{EhrenfestNA}(atoms, model)
    DynamicsMethods.create_problem(u, (0.0, 1.0), sim)
end

@testset "Energy conservation" begin
    sim = Simulation{EhrenfestNA}(atoms, model)
    u = DynamicsVariables(sim, zeros(1,1), hcat(300.0))
    tspan = (0.0, 40000.0)
    dt = 200.0
    output = (OutputTotalEnergy, OutputKineticEnergy, OutputPotentialEnergy, OutputPosition, OutputVelocity, OutputQuantumSubsystem)
    traj1 = run_dynamics(sim, tspan, u; dt, output)
    @test isapprox(var(traj1[:OutputTotalEnergy]), 0; atol=1e-6)
end
