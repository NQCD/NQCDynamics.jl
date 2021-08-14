using Test
using NonadiabaticMolecularDynamics
using OrdinaryDiffEq

@test Dynamics.Ehrenfest{Float64}(2) isa Dynamics.Ehrenfest
atoms = Atoms(:H)

@testset "Ehrenfest" begin
    sim = Simulation(atoms, NonadiabaticModels.DoubleWell(), Dynamics.Ehrenfest{Float64}(2); DoFs=1)

    r = zeros(sim.DoFs, length(sim.atoms)) 
    v = rand(sim.DoFs, length(sim.atoms)) 
    u = DynamicsVariables(sim, v, r, 1; type=:adiabatic)
    du = zero(u)

    @test Dynamics.get_quantum_subsystem(u) ≈ Complex.([1 0; 0 0])

    Dynamics.motion!(du, u, sim, 0.0)

    Calculators.evaluate_nonadiabatic_coupling!(sim.calculator)
    @test sim.calculator.nonadiabatic_coupling ≈ -sim.calculator.nonadiabatic_coupling'

    problem = ODEProblem(Dynamics.motion!, u, (0.0, 1.0), sim)
    integrator = init(problem, Tsit5(), callback=Dynamics.HoppingCallback)

    σ = get_quantum_subsystem(u)
    dσ = zero(σ)
    Dynamics.set_quantum_derivative!(dσ, v, σ, sim)

    @testset "get_diabatic_population" begin
        population = Dynamics.get_diabatic_population(sim, u)
        @test population ≈ [0.5, 0.5]
    end


    sim1 = Simulation{Ehrenfest}(atoms, NonadiabaticModels.TullyModelOne(); DoFs=1)
    r = fill(-5.0, sim1.DoFs, length(sim1.atoms))
    v1 = fill(8.9, sim1.DoFs, length(sim1.atoms)) ./ sim1.atoms.masses[1]
    z1 = DynamicsVariables(sim, v1, r, 1; type=:adiabatic)
    solution1 = Dynamics.run_trajectory(z1, (0.0, 2500.0), sim1; output=(:population), reltol=1e-6)

    @testset "run_trajectory" begin
        atoms = Atoms(2000)
        sim = Simulation{Ehrenfest}(atoms, NonadiabaticModels.TullyModelTwo(); DoFs=1)
        v = hcat(100 / 2000)
        r = hcat(-10.0)
        u = DynamicsVariables(sim, v, r, 1; type=:adiabatic)
        solution = Dynamics.run_trajectory(u, (0.0, 500.0), sim, output=(:hamiltonian), reltol=1e-6)
        @test solution.hamiltonian[1] ≈ solution.hamiltonian[end] rtol=1e-2
    end
end

@testset "Ehrenfest RPMD" begin
    atoms = Atoms(2000)
    sim = RingPolymerSimulation{Ehrenfest}(atoms, NonadiabaticModels.TullyModelTwo(), 10; DoFs=1)
    v = fill(100 / 2000, 1, 1, 10)
    r = fill(-10.0, 1, 1, 10)
    u = DynamicsVariables(sim, v, r, 1; type=:adiabatic)
    solution = Dynamics.run_trajectory(u, (0.0, 500.0), sim, output=(:hamiltonian), reltol=1e-6)
    @test solution.hamiltonian[1] ≈ solution.hamiltonian[end] rtol=1e-2
end
