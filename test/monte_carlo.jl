using Test
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.InitialConditions
using Unitful
using StatsBase

atoms = Atoms.AtomicParameters(Atoms.PeriodicCell(hcat(1)), fill(:H, 3))
model = Models.Analytic.Harmonic(1.0, 1.0, 0.1)
sys = System(atoms, model, 1; temperature=100u"K")

@testset "propose_move!" begin
    Rᵢ = zeros(n_DoF(sys), n_atoms(sys))
    Rₚ = zero(Rᵢ)
    Δ = Dict([(:H, 1.0)])
    @test Rᵢ == Rₚ
    MetropolisHastings.propose_move!(sys, Rᵢ, Rₚ, Δ, n_DoF(sys), n_atoms(sys), AnalyticWeights([0]))
    @test Rᵢ != Rₚ
end

@testset "write_output!" begin
    Rₚ = fill(0.1, n_DoF(sys), n_atoms(sys))
    output = InitialConditions.MetropolisHastings.MonteCarloOutput{Float64}(size(Rₚ), 20)
    MetropolisHastings.write_output!(output, Rₚ, 1.0, 10)
    @test all(output.R[10] .== 0.1)
end

@testset "acceptance_probability" begin
    e1 = 0.0
    e2 = -0.001
    @test MetropolisHastings.acceptance_probability(e2, e1, 1.0) == 1
    e1 = 0.0
    e2 = 0.001
    @test 1 > MetropolisHastings.acceptance_probability(e2, e1, 1.0) > 0
    e1 = 0.0
    e2 = Inf
    @test MetropolisHastings.acceptance_probability(e2, e1, 1.0) == 0
end

@testset "apply_cell_boundaries!" begin
    cell = Atoms.PeriodicCell([1 0 0; 0 1 0; 0 0 1])
    n_atoms = 4
    R = rand(3, n_atoms)
    A = copy(R)
    MetropolisHastings.apply_cell_boundaries!(cell, A, n_atoms, 3)
    @test R == A # Check unchanged when inside cell
    A = copy(R)
    A += rand(3, n_atoms) # Move atoms out of cell
    MetropolisHastings.apply_cell_boundaries!(cell, A, n_atoms, 3)
    @test all(A .<= 1) # Check they're all back in
end

@testset "run_monte_carlo_sampling" begin
    Δ = Dict([(:H, 1.0)])
    run_monte_carlo_sampling(sys, zeros(1, 3), Δ; passes=5000)
end
