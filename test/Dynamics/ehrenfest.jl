using Test
using NQCDynamics
using NQCDynamics: DynamicsMethods, Calculators
using NQCDynamics.DynamicsMethods.EhrenfestMethods
using OrdinaryDiffEq
using Statistics: var

@test Ehrenfest{Float64}(2) isa Ehrenfest
atoms = Atoms(:H)

@testset "Ehrenfest" begin
    sim = Simulation{Ehrenfest}(atoms, NQCModels.DoubleWell())

    r = zeros(size(sim)) 
    v = randn(size(sim)) 
    u = DynamicsVariables(sim, v, r, PureState(1, Adiabatic()))
    du = zero(u)

    @test DynamicsUtils.get_quantum_subsystem(u) ≈ Complex.([1 0; 0 0])

    DynamicsMethods.motion!(du, u, sim, 0.0)

    Calculators.evaluate_nonadiabatic_coupling!(sim.calculator, r)
    @test sim.calculator.nonadiabatic_coupling ≈ -sim.calculator.nonadiabatic_coupling'

    problem = ODEProblem(DynamicsMethods.motion!, u, (0.0, 1.0), sim)
    integrator = init(problem, Tsit5())

    σ = DynamicsUtils.get_quantum_subsystem(u)
    dσ = zero(σ)
    DynamicsUtils.set_quantum_derivative!(dσ, v, σ, sim)

    @testset "get_diabatic_population" begin
        population = Estimators.diabatic_population(sim, u)
        @test population ≈ [0.5, 0.5]
    end

    @testset "Algorithm comparison" begin
        atoms = Atoms(2000)
        sim = Simulation{Ehrenfest}(atoms, NQCModels.TullyModelTwo())
        v = hcat(100 / 2000)
        r = hcat(-10.0)
        output = (OutputTotalEnergy, OutputKineticEnergy, OutputPotentialEnergy, OutputPosition, OutputVelocity, OutputQuantumSubsystem, OutputAdiabaticPopulation)

        dt = 2e-3
        u = DynamicsVariables(sim, v, r, PureState(1, Adiabatic()))
        traj1 = run_dynamics(sim, (0.0, 500.0), u; output, dt)
        @test isapprox(var(traj1[:OutputTotalEnergy]), 0; atol=1e-6)

        u = DynamicsVariables(sim, v, r, PureState(1, Adiabatic()))
        traj2 = run_dynamics(sim, (0.0, 500.0), u; algorithm=Vern9(), output, reltol=1e-12, abstol=1e-12, saveat=dt)
        @test isapprox(var(traj2[:OutputTotalEnergy]), 0; atol=1e-6)

        @test traj1[:OutputKineticEnergy] ≈ traj2[:OutputKineticEnergy] rtol=1e-3
        @test traj1[:OutputPotentialEnergy] ≈ traj2[:OutputPotentialEnergy] rtol=1e-1
        @test traj1[:OutputVelocity] ≈ traj2[:OutputVelocity] rtol=1e-3
        @test traj1[:OutputPosition] ≈ traj2[:OutputPosition] rtol=1e-3
        @test traj1[:OutputQuantumSubsystem] ≈ traj2[:OutputQuantumSubsystem] rtol=1e-2
    end

    @testset "FermiDirac Ehrenfest algorithm comparison" begin
        kT = 9.5e-4
        M = 30 # number of bath states
        Γ = 6.4e-3
        W = 6Γ / 2 # bandwidth  parameter

        basemodel = MiaoSubotnik(;Γ)
        bath = TrapezoidalRule(M, -W, W)
        model = AndersonHolstein(basemodel, bath)
        atoms = Atoms(2000)

        sim = Simulation{Ehrenfest}(atoms, model)
        v = zeros(1,1)
        r = hcat(21.0)
        u = DynamicsVariables(sim, v, r, FermiDiracState(0.0, 0.0))
        tspan = (0.0, 2000.0)
        dt = 10.0
        output = (OutputTotalEnergy, OutputKineticEnergy, OutputPotentialEnergy, OutputPosition, OutputVelocity, OutputQuantumSubsystem)
        traj1 = run_dynamics(sim, tspan, u; dt, output, algorithm=Vern9(), abstol=1e-15, reltol=1e-15, saveat=dt)
        @test isapprox(var(traj1[:OutputTotalEnergy]), 0; atol=1e-6)
    
        u = DynamicsVariables(sim, v, r, FermiDiracState(0.0, 0.0))
        traj2 = run_dynamics(sim, tspan, u; dt, output) # default algorithm is fixed timestep
        @test isapprox(var(traj2[:OutputTotalEnergy]), 0; atol=1e-6)

        @test traj1[:OutputKineticEnergy] ≈ traj2[:OutputKineticEnergy] rtol=1e-3
        @test traj1[:OutputPotentialEnergy] ≈ traj2[:OutputPotentialEnergy] rtol=1e-3
        @test traj1[:OutputVelocity] ≈ traj2[:OutputVelocity] rtol=1e-3
        @test traj1[:OutputPosition] ≈ traj2[:OutputPosition] rtol=1e-3
        @test traj1[:OutputQuantumSubsystem] ≈ traj2[:OutputQuantumSubsystem] rtol=1e-2
    end
end

@testset "Ehrenfest RPMD" begin
    atoms = Atoms(2000)
    sim = RingPolymerSimulation{Ehrenfest}(atoms, NQCModels.TullyModelTwo(), 10)
    v = fill(100 / 2000, 1, 1, 10)
    r = fill(-10.0, 1, 1, 10)
    u = DynamicsVariables(sim, v, r, PureState(1, Adiabatic()))
    solution = run_dynamics(sim, (0.0, 500.0), u, output=OutputTotalEnergy, dt=0.01)
    @test isapprox(var(solution[:OutputTotalEnergy]), 0; atol=1e-6)
end
