using Test
using NQCDynamics
using LinearAlgebra
using Random
using Distributions
using NQCDynamics: DynamicsMethods, DynamicsUtils, Calculators
using NQCDynamics.DynamicsMethods: SurfaceHoppingMethods
using OrdinaryDiffEq
using ComponentArrays

kT = 9.5e-4
M = 30 # number of bath states
Γ = 6.4e-3
W = 6Γ / 2 # bandwidth  parameter
model = MiaoSubotnik(M=M, W=W, Γ=Γ)
atoms = Atoms(2000)
r = randn(1,1)
v = randn(1,1)
n_electrons = M ÷ 2

sim = Simulation{AdiabaticIESH}(atoms, model; n_electrons=n_electrons)
u = DynamicsVariables(sim, v, r)
sim.method.state .= u.state

@testset "compute_overlap!" begin
    S = zeros(n_electrons, n_electrons)
    state = collect(1:n_electrons)
    ψ = DynamicsUtils.get_quantum_subsystem(u)
    SurfaceHoppingMethods.compute_overlap!(sim, S, ψ, state)
    @test S == I # Check overlap is identity

    # All electrons in state 1, so they all have unity overlap
    # with first electron wavefunction, which is in state 1.
    state = ones(Int, n_electrons)
    SurfaceHoppingMethods.compute_overlap!(sim, S, ψ, state)
    @test all(S[:,1] .== 1) # Check only first column ones
    @test all(S[:,2:end] .== 0)

    @testset "adiabatic vs diabatic" begin
        sim = Simulation{DiabaticIESH}(atoms, model; n_electrons=n_electrons)
        u = DynamicsVariables(sim, v, r)
        Sdiabatic = zeros(n_electrons, n_electrons)
        ψ = DynamicsUtils.get_quantum_subsystem(u)
        SurfaceHoppingMethods.compute_overlap!(sim, Sdiabatic, ψ, u.state)

        sim = Simulation{AdiabaticIESH}(atoms, model; n_electrons=n_electrons)
        u = DynamicsVariables(sim, v, r)
        Sadiabatic = zeros(n_electrons, n_electrons)
        ψ = DynamicsUtils.get_quantum_subsystem(u)
        SurfaceHoppingMethods.compute_overlap!(sim, Sadiabatic, ψ, u.state)

        @test Sdiabatic ≈ Sadiabatic
    end
end

@testset "calculate_Akj" begin
    S = zeros(n_electrons, n_electrons)
    state = collect(1:n_electrons)
    ψ = DynamicsUtils.get_quantum_subsystem(u)
    Akj = SurfaceHoppingMethods.calculate_Akj(sim, S, ψ, 1.0, state)
    @test Akj ≈ 1.0
end

@testset "evaluate_hopping_probability!" begin
    Calculators.update_electronics!(sim.calculator, hcat(20.0))
    z = deepcopy(u)
    rand!(DynamicsUtils.get_quantum_subsystem(z))
    SurfaceHoppingMethods.evaluate_hopping_probability!(sim, z, 1.0)
end

@testset "execute_hop!" begin
    problem = ODEProblem(DynamicsMethods.motion!, u, (0.0, 1.0), sim)
    integrator = init(problem, Tsit5(), callback=SurfaceHoppingMethods.HoppingCallback)

    initial_state = collect(1:n_electrons)
    final_state = collect(1:n_electrons)
    final_state[end] = final_state[end] + 1

    Calculators.update_electronics!(integrator.p.calculator, get_positions(integrator.u))
    DynamicsUtils.get_velocities(integrator.u) .= 2 # Set high momentum to ensure successful hop
    integrator.u.state .= initial_state
    integrator.p.method.new_state .= final_state
    eigs = SurfaceHoppingMethods.get_hopping_eigenvalues(integrator.p)
    new_state, old_state = SurfaceHoppingMethods.unpack_states(sim)
    ΔE = SurfaceHoppingMethods.calculate_potential_energy_change(eigs, new_state, old_state)

    KE_initial = DynamicsUtils.classical_kinetic_energy(integrator.p, get_velocities(integrator.u))
    H_initial = DynamicsUtils.classical_hamiltonian(integrator.p, integrator.u)
    @test integrator.u.state == initial_state

    SurfaceHoppingMethods.execute_hop!(integrator)

    KE_final = DynamicsUtils.classical_kinetic_energy(integrator.p, get_velocities(integrator.u))
    H_final = DynamicsUtils.classical_hamiltonian(integrator.p, integrator.u)
    ΔKE = KE_final - KE_initial

    @test integrator.u.state == final_state # Check state has changed
    @test H_final ≈ H_initial
    @test ΔKE ≈ -ΔE rtol=1e-3 # Test for energy conservation
end

@testset "algorithm comparison" begin
    sim = Simulation{AdiabaticIESH}(atoms, model; n_electrons=n_electrons)
    u = DynamicsVariables(sim, zeros(1,1), randn(1,1))
    tspan = (0.0, 100.0)
    dt = 1.0
    traj1 = run_trajectory(u, tspan, sim; dt, algorithm=DynamicsMethods.IntegrationAlgorithms.VerletwithElectronics(), output=:u)
    traj2 = run_trajectory(u, tspan, sim; algorithm=Tsit5(), output=:u, saveat=tspan[1]:dt:tspan[2])
    @test traj1.u ≈ traj2.u atol=1e-2
end
