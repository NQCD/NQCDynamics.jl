using Test
using NonadiabaticMolecularDynamics
using OrdinaryDiffEq

@test Dynamics.ehrenfest_basic{Float64}(2) isa Dynamics.ehrenfest_basic
atoms = Atoms(1)

@testset "Ehrenfest" begin
    sim = Simulation(atoms, Models.DoubleWell(), Dynamics.ehrenfest_basic{Float64}(2); DoFs=1)
    #sim.method.state = 1

    r = zeros(sim.DoFs, length(sim.atoms)) 
    v = rand(sim.DoFs, length(sim.atoms)) 
    u = EhrenfestVariables(v, r, 2, 1)
    du = zero(u)

    @test Dynamics.get_density_matrix(u) ≈ Complex.([1 0; 0 0])

    Dynamics.motion!(du, u, sim, 0.0)

    Calculators.evaluate_nonadiabatic_coupling!(sim.calculator)
    @test sim.calculator.nonadiabatic_coupling ≈ -sim.calculator.nonadiabatic_coupling'

    problem = ODEProblem(Dynamics.motion!, u, (0.0, 1.0), sim)
    integrator = init(problem, Tsit5(), callback=Dynamics.HoppingCallback)

    σ = get_density_matrix(u)
    dσ = zero(σ)
    Dynamics.set_density_matrix_derivative!(dσ, v, σ, sim)

    @testset "get_population" begin
        population = Dynamics.get_population(sim, u)
        @test population ≈ [0.5, 0.5]
    end

    @testset "calculate_potential_energy_change" begin
        integrator.p.calculator.eigenvalues .= [0.9, -0.3]
        @test 1.2 ≈ Dynamics.calculate_potential_energy_change(integrator.p.calculator, 1, 2)
        @test -1.2 ≈ Dynamics.calculate_potential_energy_change(integrator.p.calculator, 2, 1)
    end

    @testset "run_trajectory" begin
        atoms = Atoms(2000)
        sim = Simulation{ehrenfest_basic}(atoms, Models.TullyModelTwo(); DoFs=1)
        v = hcat(100 / 2000)
        r = hcat(-10.0)
        u = EhrenfestVariables(v, r, 2, 1)
        solution = Dynamics.run_trajectory(u, (0.0, 500.0), sim, output=(:hamiltonian, :state), reltol=1e-6)
        @test solution.hamiltonian[1] ≈ solution.hamiltonian[end]
    end
end
