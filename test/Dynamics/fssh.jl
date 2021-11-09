using Test
using NonadiabaticMolecularDynamics
using OrdinaryDiffEq
using StaticArrays
using Random: seed!
using NonadiabaticMolecularDynamics: DynamicsMethods, DynamicsUtils, Calculators
using NonadiabaticMolecularDynamics.SurfaceHoppingMethods: SurfaceHoppingMethods
using NonadiabaticMolecularDynamics.DynamicsUtils: get_positions, get_velocities

@test SurfaceHoppingMethods.FSSH{Float64}(2, :standard) isa SurfaceHoppingMethods.FSSH
atoms = Atoms(1)

@testset "FSSH" begin
    sim = Simulation{FSSH}(atoms, NonadiabaticModels.DoubleWell())
    sim.method.state = 1

    r = zeros(size(sim)) 
    v = randn(size(sim)) 
    u = DynamicsVariables(sim, v, r, SingleState(1, Adiabatic()))
    du = zero(u)

    @test DynamicsUtils.get_quantum_subsystem(u) ≈ Complex.([1 0; 0 0])

    DynamicsMethods.motion!(du, u, sim, 0.0)

    Calculators.evaluate_nonadiabatic_coupling!(sim.calculator)
    @test sim.calculator.nonadiabatic_coupling ≈ -sim.calculator.nonadiabatic_coupling'

    problem = ODEProblem(DynamicsMethods.motion!, u, (0.0, 1.0), sim)
    integrator = init(problem, Tsit5(), callback=SurfaceHoppingMethods.HoppingCallback)
    DynamicsMethods.SurfaceHoppingMethods.evaluate_hopping_probability!(sim, u, get_proposed_dt(integrator))

    σ = DynamicsUtils.get_quantum_subsystem(u)
    dσ = zero(σ)
    DynamicsUtils.set_quantum_derivative!(dσ, v, σ, sim)

    @testset "get_diabatic_population" begin
        population = Estimators.diabatic_population(sim, u)
        @test population ≈ [0.5, 0.5]
    end

    @testset "select_new_state" begin
        sim.method.hopping_probability .= [0, 1.0]
        u.state = 1
        @test 2 == SurfaceHoppingMethods.select_new_state(sim, u)
        sim.method.hopping_probability .= [0, 1.0]
        u.state = 2
        @test 2 == SurfaceHoppingMethods.select_new_state(sim, u)
        sim.method.hopping_probability .= [1.0, 0.0]
        u.state = 2
        @test 1 == SurfaceHoppingMethods.select_new_state(sim, u)
    end

    @testset "rescale_velocity!" begin
        get_velocities(integrator.u) .= 0.0 # Set momentum to zero to force frustrated hop
        integrator.p.method.new_state = 2
        cont = SurfaceHoppingMethods.rescale_velocity!(integrator.p, integrator.u)
        @test cont == false # Check hop is rejected 

        get_velocities(integrator.u) .= 1e5 # Set momentum to big to force hop
        cont = SurfaceHoppingMethods.rescale_velocity!(integrator.p, integrator.u)
        @test cont == true # Check hop is accepted
    end

    @testset "calculate_potential_energy_change" begin
        integrator.p.calculator.eigenvalues = SVector{2}([0.9, -0.3])
        eigs = SurfaceHoppingMethods.get_hopping_eigenvalues(integrator.p)
        @test 1.2 ≈ SurfaceHoppingMethods.calculate_potential_energy_change(eigs, 1, 2)
        @test -1.2 ≈ SurfaceHoppingMethods.calculate_potential_energy_change(eigs, 2, 1)
    end

    @testset "execute_hop!" begin
        Calculators.update_electronics!(integrator.p.calculator, get_positions(integrator.u))
        DynamicsUtils.get_velocities(integrator.u) .= 2 # Set high momentum to ensure successful hop
        integrator.u.state = 1
        integrator.p.method.new_state = 2
        eigs = SurfaceHoppingMethods.get_hopping_eigenvalues(integrator.p)
        ΔE = SurfaceHoppingMethods.calculate_potential_energy_change(eigs, 2, 1)

        KE_initial = DynamicsUtils.classical_kinetic_energy(integrator.p, get_velocities(integrator.u))
        H_initial = DynamicsUtils.classical_hamiltonian(integrator.p, integrator.u)
        @test integrator.u.state == 1

        SurfaceHoppingMethods.execute_hop!(integrator)

        KE_final = DynamicsUtils.classical_kinetic_energy(integrator.p, get_velocities(integrator.u))
        H_final = DynamicsUtils.classical_hamiltonian(integrator.p, integrator.u)
        ΔKE = KE_final - KE_initial

        @test integrator.u.state == 2 # Check state has changed
        @test H_final ≈ H_initial
        @test ΔKE ≈ -ΔE rtol=1e-3 # Test for energy conservation
    end

    @testset "run_trajectory" begin
        atoms = Atoms(2000)
        sim = Simulation{FSSH}(atoms, NonadiabaticModels.TullyModelTwo())
        v = hcat(100 / 2000)
        r = hcat(-10.0)
        u = DynamicsVariables(sim, v, r, SingleState(1, Adiabatic()))
        solution = run_trajectory(u, (0.0, 500.0), sim, output=(:hamiltonian, :state), reltol=1e-6)
        @test solution.hamiltonian[1] ≈ solution.hamiltonian[end] rtol=1e-2
    end
end

@testset "RPSH" begin
    sim = RingPolymerSimulation{FSSH}(atoms, NonadiabaticModels.DoubleWell(), 5)

    r = RingPolymerArray(zeros(size(sim)))
    v = RingPolymerArray(randn(size(sim)))
    u = DynamicsVariables(sim, v, r, SingleState(1, Adiabatic()))
    sim.method.state = u.state

    problem = ODEProblem(DynamicsMethods.motion!, u, (0.0, 1.0), sim)
    integrator = init(problem, Tsit5(), callback=SurfaceHoppingMethods.HoppingCallback)
    SurfaceHoppingMethods.evaluate_hopping_probability!(sim, u, get_proposed_dt(integrator))

    @testset "rescale_velocity!" begin
        get_velocities(integrator.u) .= 0.0 # Set momentum to zero to force frustrated hop
        integrator.p.method.new_state = 2
        cont = SurfaceHoppingMethods.rescale_velocity!(integrator.p, integrator.u)
        @test cont == false # Check hop is rejected 

        get_velocities(integrator.u) .= 1e5 # Set momentum to big to force hop
        cont = SurfaceHoppingMethods.rescale_velocity!(integrator.p, integrator.u)
        @test cont == true # Check hop is accepted
    end

    @testset "calculate_potential_energy_change" begin
        integrator.p.calculator.centroid_eigenvalues = SVector{2}([0.9, -0.3])
        eigs = SurfaceHoppingMethods.get_hopping_eigenvalues(integrator.p)
        @test 1.2 ≈ SurfaceHoppingMethods.calculate_potential_energy_change(eigs, 1, 2)
        @test -1.2 ≈ SurfaceHoppingMethods.calculate_potential_energy_change(eigs, 2, 1)
    end

    @testset "execute_hop!" begin
        Calculators.update_electronics!(integrator.p.calculator, get_positions(integrator.u))
        get_velocities(integrator.u) .= 5 # Set high momentum to ensure successful hop
        integrator.u.state = 1
        integrator.p.method.new_state = 2
        eigs = SurfaceHoppingMethods.get_hopping_eigenvalues(integrator.p)
        ΔE = SurfaceHoppingMethods.calculate_potential_energy_change(eigs, 2, 1)

        KE_initial = DynamicsUtils.classical_kinetic_energy(integrator.p, get_velocities(integrator.u))
        H_initial = DynamicsUtils.classical_hamiltonian(integrator.p, integrator.u)
        potential = H_initial - KE_initial
        @test integrator.u.state == 1

        SurfaceHoppingMethods.execute_hop!(integrator)

        KE_final = DynamicsUtils.classical_kinetic_energy(integrator.p, get_velocities(integrator.u))
        H_final = DynamicsUtils.classical_hamiltonian(integrator.p, integrator.u)
        potential = H_final - KE_final

        @test integrator.u.state == 2 # Check state has changed
        @test H_initial ≈ H_final # Check ring polymer hamiltonian conserved
        ΔKE = (KE_final - KE_initial) / length(sim.beads)
        @test ΔKE ≈ -sum(ΔE) rtol=1e-3 # Test for energy conservation
    end

    @testset "run_trajectory" begin
        atoms = Atoms(2000)
        sim = RingPolymerSimulation{FSSH}(atoms, TullyModelTwo(), 5; temperature=0.01)
        v = fill(100 / 2000, size(sim))
        r = fill(-10.0, size(sim)) .+ randn(1,1,5)
        u = DynamicsVariables(sim, v, r, SingleState(1, Adiabatic()))
        seed!(1)
        solution = run_trajectory(u, (0.0, 1000.0), sim, output=(:hamiltonian), dt=0.1)
        seed!(1)
        solution2 = run_trajectory(u, (0.0, 1000.0), sim, output=(:hamiltonian), algorithm=Tsit5(), abstol=1e-10, reltol=1e-10, saveat=0:0.1:1000.0)
        # Ring polymer Hamiltonian is not strictly conserved during hoppping
        @test solution.hamiltonian[1] ≈ solution.hamiltonian[end] rtol=1e-2
    end

end
