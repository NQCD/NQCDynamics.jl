using Test
using NQCDynamics
using LinearAlgebra
using Distributions
using NQCDynamics: DynamicsMethods, DynamicsUtils, Calculators
using NQCDynamics.DynamicsMethods: SurfaceHoppingMethods
using OrdinaryDiffEq
using ComponentArrays
using Unitful, UnitfulAtomic

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

sim = Simulation{AdiabaticIESH}(atoms, model)
u = DynamicsVariables(sim, v, r)
sim.method.state .= u.state
SurfaceHoppingMethods.set_unoccupied_states!(sim)

@testset "DynamicsVariables" begin
    model = AndersonHolstein(basemodel, bath; fermi_level=0.0u"eV")
    sim = Simulation{AdiabaticIESH}(atoms, model)
    distribution = NQCDistributions.FermiDiracState(0.0u"eV", 0u"K")
    u = DynamicsVariables(sim, v, r, distribution)
    @test u.state == 1:n_electrons

    model = AndersonHolstein(basemodel, bath; fermi_level=0.1u"eV")
    sim = Simulation{AdiabaticIESH}(atoms, model)
    distribution = NQCDistributions.FermiDiracState(0.0u"eV", 300u"K")
    @test_throws ErrorException DynamicsVariables(sim, v, r, distribution)

    model = AndersonHolstein(basemodel, bath; fermi_level=0.1u"eV")
    sim = Simulation{AdiabaticIESH}(atoms, model)
    distribution = NQCDistributions.FermiDiracState(0.1u"eV", 300u"K")
    avg = zeros(nstates(model))
    samples = 5000
    for i=1:samples
        u = DynamicsVariables(sim, v, r, distribution)
        avg[u.state] .+= 1
    end
    avg ./= samples
    eigs = Calculators.get_eigen(sim.calculator, r)
    occupations = NQCDistributions.fermi.(eigs.values, distribution.fermi_level, distribution.β)       

    @test avg ≈ occupations atol=0.2
end

@testset "set_unoccupied_states!" begin
    sim = Simulation{AdiabaticIESH}(atoms, model)
    for (occupied, unoccupied) in zip([1:15, 6:20], [16:31, vcat(1:5, 21:31)])
        sim.method.state .= occupied
        SurfaceHoppingMethods.set_unoccupied_states!(sim)
        @test sim.method.unoccupied == unoccupied
    end
end

@testset "create_problem" begin
    sim = Simulation{AdiabaticIESH}(atoms, model)
    DynamicsMethods.create_problem(u, (0.0, 1.0), sim)
    @test sim.method.state == 1:15
    @test sim.method.unoccupied == 16:31
end

@testset "evaluate_v_dot_d!" begin
    d = Calculators.get_nonadiabatic_coupling(sim.calculator, r)
    SurfaceHoppingMethods.evaluate_v_dot_d!(sim, v, d)
end

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
    ψ = DynamicsUtils.get_quantum_subsystem(z)
    rand!(ψ)
    @views for i in axes(ψ,2)
        normalize!(ψ[:,i])
    end
    SurfaceHoppingMethods.evaluate_hopping_probability!(sim, z, 1.0, rand())
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
    eigs = DynamicsUtils.get_hopping_eigenvalues(integrator.p, DynamicsUtils.get_positions(integrator.u))
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
    sim = Simulation{AdiabaticIESH}(atoms, model)
    u = DynamicsVariables(sim, randn(1,1), randn(1,1))
    tspan = (0.0, 100.0)
    dt = 1.0
    output = (OutputPosition, OutputVelocity, OutputQuantumSubsystem, OutputTotalEnergy)
    @time traj1 = run_dynamics(sim, tspan, u; dt, algorithm=DynamicsMethods.IntegrationAlgorithms.VerletwithElectronics(), output)
    @time traj2 = run_dynamics(sim, tspan, u; algorithm=Tsit5(), saveat=tspan[1]:dt:tspan[2], output, reltol=1e-6, abstol=1e-6)
    @time traj3 = run_dynamics(sim, tspan, u; dt, algorithm=DynamicsMethods.IntegrationAlgorithms.VerletwithElectronics2(MagnusMidpoint(krylov=false), dt=dt/5), output)

    # We cannot compare these when hopping is happening since the trajectories will be different. 
end

@testset "DecoherenceCorrectionEDC" begin
    sim = Simulation{AdiabaticIESH}(atoms, model; decoherence=SurfaceHoppingMethods.DecoherenceCorrectionEDC())
    u = DynamicsVariables(sim, zeros(1,1), 100randn(1,1))
    tspan = (0.0, 100000.0)
    dt = 100.0
    traj = run_dynamics(sim, tspan, u; dt, output=(OutputQuantumSubsystem, OutputSurfaceHops))
    @test traj[:OutputSurfaceHops] > 0 # Ensure some hops occur
    norms = zeros(length(traj[:Time]))
    for i in eachindex(norms)
        ψ = traj[:OutputQuantumSubsystem][i]
        @views for j in axes(ψ,2)
            norms[i] += norm(ψ[:,j])
        end
    end
    @test all(i->isapprox(i, n_electrons; rtol=1e-3), norms) # Test that decoherence conserves norm
end
