using Test
using NonadiabaticMolecularDynamics
using OrdinaryDiffEq

@test Dynamics.FSSH{Float64}(2) isa Dynamics.FSSH
atoms = Atoms([:H])

@testset "FSSH" begin
    sim = Simulation(atoms, Models.DoubleWell(), Dynamics.FSSH{Float64}(2); DoFs=1)

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
        get_velocities(integrator.u) .= 1e4 # Set high momentum to ensure successful hop
        integrator.p.method.new_state = 2
        ΔE = Dynamics.calculate_potential_energy_change(integrator.p.calculator, 2, 1)

        KE_initial = sum(get_velocities(integrator.u).^2 .* integrator.p.atoms.masses')/2
        @test integrator.u.state == 1
        Dynamics.execute_hop!(integrator)
        KE_final = sum(get_velocities(integrator.u).^2 .* integrator.p.atoms.masses')/2
        @test integrator.u.state == 2 # Check state has changed
        @test KE_final - KE_initial ≈ ΔE rtol=1e-3 # Test for energy conservation
    end

    @testset "run_trajectory" begin
        solution = Dynamics.run_trajectory(u, (0.0, 1.0), sim)
    end
end

@testset "RPSH" begin
    sim = RingPolymerSimulation{FSSH}(atoms, Models.DoubleWell(), 5; DoFs=1)

    r = RingPolymerArray(zeros(sim.DoFs, length(sim.atoms), 5))
    v = RingPolymerArray(rand(sim.DoFs, length(sim.atoms), 5))
    u = SurfaceHoppingVariables(v, r, 2, 1)

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
        @test 1.2 ≈ Dynamics.calculate_potential_energy_change(integrator.p.calculator, 1, 2)
        @test -1.2 ≈ Dynamics.calculate_potential_energy_change(integrator.p.calculator, 2, 1)
    end

    @testset "execute_hop!" begin
        get_velocities(integrator.u) .= 1e4 # Set high momentum to ensure successful hop
        integrator.p.method.new_state = 2
        ΔE = Dynamics.calculate_potential_energy_change(integrator.p.calculator, 2, 1)

        KE_initial = sum(get_velocities(integrator.u).^2 .* integrator.p.atoms.masses')/2 ./ 5
        @test integrator.u.state == 1
        Dynamics.execute_hop!(integrator)
        KE_final = sum(get_velocities(integrator.u).^2 .* integrator.p.atoms.masses')/2 ./ 5
        @test integrator.u.state == 2 # Check state has changed
        @test KE_final - KE_initial ≈ ΔE rtol=1e-3 # Test for energy conservation
    end

    @testset "run_trajectory" begin

        sim = RingPolymerSimulation{FSSH}(atoms, Models.DoubleWell(), 5; DoFs=1)

        r = RingPolymerArray(zeros(sim.DoFs, length(sim.atoms), 5))
        v = RingPolymerArray(rand(sim.DoFs, length(sim.atoms), 5))
        u = SurfaceHoppingVariables(v, r, 2, 1)

        solution = Dynamics.run_trajectory(u, (0.0, 1.0), sim)
    end

end
