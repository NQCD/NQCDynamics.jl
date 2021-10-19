using Test
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics: DynamicsMethods, Calculators
using NonadiabaticMolecularDynamics.DynamicsMethods.EhrenfestMethods
using OrdinaryDiffEq

@test Ehrenfest{Float64}(2) isa Ehrenfest
atoms = Atoms(:H)

@testset "Ehrenfest" begin
    sim = Simulation{Ehrenfest}(atoms, NonadiabaticModels.DoubleWell())

    r = zeros(size(sim)) 
    v = randn(size(sim)) 
    u = DynamicsVariables(sim, v, r, 1; type=:adiabatic)
    du = zero(u)

    @test DynamicsUtils.get_quantum_subsystem(u) ≈ Complex.([1 0; 0 0])

    DynamicsMethods.motion!(du, u, sim, 0.0)

    Calculators.evaluate_nonadiabatic_coupling!(sim.calculator)
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


    sim1 = Simulation{Ehrenfest}(atoms, NonadiabaticModels.TullyModelOne())
    r = fill(-5.0, size(sim))
    v1 = fill(8.9, size(sim)) ./ sim1.atoms.masses[1]
    z1 = DynamicsVariables(sim, v1, r, 1; type=:adiabatic)
    solution1 = run_trajectory(z1, (0.0, 2500.0), sim1; output=(:population), reltol=1e-6)

    @testset "run_trajectory" begin
        atoms = Atoms(2000)
        sim = Simulation{Ehrenfest}(atoms, NonadiabaticModels.TullyModelTwo())
        v = hcat(100 / 2000)
        r = hcat(-10.0)
        u = DynamicsVariables(sim, v, r, 1; type=:adiabatic)
        solution = run_trajectory(u, (0.0, 500.0), sim, output=(:hamiltonian), reltol=1e-6)
        @test solution.hamiltonian[1] ≈ solution.hamiltonian[end] rtol=1e-2
    end
end

@testset "Ehrenfest RPMD" begin
    atoms = Atoms(2000)
    sim = RingPolymerSimulation{Ehrenfest}(atoms, NonadiabaticModels.TullyModelTwo(), 10)
    v = fill(100 / 2000, 1, 1, 10)
    r = fill(-10.0, 1, 1, 10)
    u = DynamicsVariables(sim, v, r, 1; type=:adiabatic)
    solution = run_trajectory(u, (0.0, 500.0), sim, output=(:hamiltonian), dt=1)
    @test solution.hamiltonian[1] ≈ solution.hamiltonian[end] rtol=1e-2
end
