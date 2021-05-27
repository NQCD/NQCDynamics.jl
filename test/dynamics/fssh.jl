using Test
using NonadiabaticMolecularDynamics
using OrdinaryDiffEq

@test Dynamics.FSSH{Float64}(2) isa Dynamics.FSSH
atoms = Atoms(1)

@testset "FSSH" begin
    sim = Simulation(atoms, Models.DoubleWell(), Dynamics.FSSH{Float64}(2); DoFs=1)
    sim.method.state = 1

    r = zeros(sim.DoFs, length(sim.atoms)) 
    v = rand(sim.DoFs, length(sim.atoms)) 
    u = SurfaceHoppingVariables(v, r, 2, 1)
    du = zero(u)

    @test Dynamics.get_density_matrix(u) ≈ Complex.([1 0; 0 0])

    Dynamics.motion!(du, u, sim, 0.0)

    Calculators.evaluate_nonadiabatic_coupling!(sim.calculator)
    @test sim.calculator.nonadiabatic_coupling ≈ -sim.calculator.nonadiabatic_coupling'

    problem = ODEProblem(Dynamics.motion!, u, (0.0, 1.0), sim)
    integrator = init(problem, Tsit5(), callback=Dynamics.HoppingCallback)
    Dynamics.evaluate_hopping_probability!(sim, u, get_proposed_dt(integrator))

    σ = get_density_matrix(u)
    dσ = zero(σ)
    Dynamics.set_density_matrix_derivative!(dσ, v, σ, sim)

    @testset "get_diabatic_population" begin
        population = Dynamics.get_diabatic_population(sim, u)
        @test population ≈ [0.5, 0.5]
    end

    @testset "select_new_state" begin
        sim.method.hopping_probability .= [0, 1.0]
        u.state = 1
        @test 2 == Dynamics.select_new_state(sim, u)
        sim.method.hopping_probability .= [0, 1.0]
        u.state = 2
        @test 2 == Dynamics.select_new_state(sim, u)
        sim.method.hopping_probability .= [1.0, 0.0]
        u.state = 2
        @test 1 == Dynamics.select_new_state(sim, u)
    end

    @testset "rescale_velocity!" begin
        get_velocities(integrator.u) .= 0.0 # Set momentum to zero to force frustrated hop
        integrator.p.method.new_state = 2
        cont = Dynamics.rescale_velocity!(integrator.p, integrator.u)
        @test cont == false # Check hop is rejected 

        get_velocities(integrator.u) .= 1e5 # Set momentum to big to force hop
        cont = Dynamics.rescale_velocity!(integrator.p, integrator.u)
        @test cont == true # Check hop is accepted
    end

    @testset "calculate_potential_energy_change" begin
        integrator.p.calculator.eigenvalues .= [0.9, -0.3]
        @test 1.2 ≈ Dynamics.calculate_potential_energy_change(integrator.p.calculator, 1, 2)
        @test -1.2 ≈ Dynamics.calculate_potential_energy_change(integrator.p.calculator, 2, 1)
    end

    @testset "execute_hop!" begin
        Calculators.update_electronics!(integrator.p.calculator, get_positions(integrator.u))
        get_velocities(integrator.u) .= 2 # Set high momentum to ensure successful hop
        integrator.u.state = 1
        integrator.p.method.new_state = 2
        ΔE = Dynamics.calculate_potential_energy_change(integrator.p.calculator, 2, 1)

        KE_initial = evaluate_kinetic_energy(atoms.masses, get_velocities(integrator.u))
        H_initial = evaluate_hamiltonian(integrator.p, integrator.u)
        @test integrator.u.state == 1

        Dynamics.execute_hop!(integrator)

        KE_final = evaluate_kinetic_energy(atoms.masses, get_velocities(integrator.u))
        H_final = evaluate_hamiltonian(integrator.p, integrator.u)
        ΔKE = KE_final - KE_initial

        @test integrator.u.state == 2 # Check state has changed
        @test H_final ≈ H_initial
        @test ΔKE ≈ -ΔE rtol=1e-3 # Test for energy conservation
    end

    @testset "run_trajectory" begin
        atoms = Atoms(2000)
        sim = Simulation{FSSH}(atoms, Models.TullyModelTwo(); DoFs=1)
        v = hcat(100 / 2000)
        r = hcat(-10.0)
        u = SurfaceHoppingVariables(v, r, 2, 1)
        solution = Dynamics.run_trajectory(u, (0.0, 500.0), sim, output=(:hamiltonian, :state), reltol=1e-6)
        @test solution.hamiltonian[1] ≈ solution.hamiltonian[end]
    end
end

@testset "RPSH" begin
    sim = RingPolymerSimulation{FSSH}(atoms, Models.DoubleWell(), 5; DoFs=1)

    r = RingPolymerArray(zeros(sim.DoFs, length(sim.atoms), 5))
    v = RingPolymerArray(rand(sim.DoFs, length(sim.atoms), 5))
    u = SurfaceHoppingVariables(v, r, 2, 1)
    sim.method.state = u.state

    problem = ODEProblem(Dynamics.motion!, u, (0.0, 1.0), sim)
    integrator = init(problem, Tsit5(), callback=Dynamics.HoppingCallback)
    Dynamics.evaluate_hopping_probability!(sim, u, get_proposed_dt(integrator))

    @testset "rescale_velocity!" begin
        get_velocities(integrator.u) .= 0.0 # Set momentum to zero to force frustrated hop
        integrator.p.method.new_state = 2
        cont = Dynamics.rescale_velocity!(integrator.p, integrator.u)
        @test cont == false # Check hop is rejected 

        get_velocities(integrator.u) .= 1e5 # Set momentum to big to force hop
        cont = Dynamics.rescale_velocity!(integrator.p, integrator.u)
        @test cont == true # Check hop is accepted
    end

    @testset "calculate_potential_energy_change" begin
        integrator.p.calculator.eigenvalues .= [[0.9, -0.3] for i=1:5]
        @test fill(1.2, 5) ≈ Dynamics.calculate_potential_energy_change(integrator.p.calculator, 1, 2)
        @test fill(-1.2, 5) ≈ Dynamics.calculate_potential_energy_change(integrator.p.calculator, 2, 1)
    end

    @testset "execute_hop!" begin
        Calculators.update_electronics!(integrator.p.calculator, get_positions(integrator.u))
        get_velocities(integrator.u) .= 5 # Set high momentum to ensure successful hop
        integrator.u.state = 1
        integrator.p.method.new_state = 2
        ΔE = Dynamics.calculate_potential_energy_change(integrator.p.calculator, 2, 1)

        KE_initial = evaluate_kinetic_energy(atoms.masses, get_velocities(integrator.u))
        H_initial = evaluate_hamiltonian(integrator.p, integrator.u)
        potential = H_initial - KE_initial
        @test integrator.u.state == 1

        Dynamics.execute_hop!(integrator)

        KE_final = evaluate_kinetic_energy(atoms.masses, get_velocities(integrator.u))
        H_final = evaluate_hamiltonian(integrator.p, integrator.u)
        potential = H_final - KE_final

        @test integrator.u.state == 2 # Check state has changed
        @test H_initial ≈ H_final # Check ring polymer hamiltonian conserved
        ΔKE = KE_final - KE_initial
        @test ΔKE ≈ -sum(ΔE) rtol=1e-3 # Test for energy conservation
    end

    @testset "run_trajectory" begin
        atoms = Atoms(2000)
        sim = RingPolymerSimulation{FSSH}(atoms, Models.TullyModelTwo(), 5; DoFs=1, temperature=0.01)
        v = RingPolymerArray(fill(100 / 2000, 1, 1, 5))
        r = RingPolymerArray(fill(-10.0, 1, 1, 5)) .+ randn(1,1,5)
        u = SurfaceHoppingVariables(v, r, 2, 1)
        solution = Dynamics.run_trajectory(u, (0.0, 1000.0), sim, output=(:hamiltonian), reltol=1e-10)
        @test solution.hamiltonian[1] ≈ solution.hamiltonian[end] rtol=1e-3
    end

end
