using Test
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.InitialConditions
using Unitful

atoms = Atoms.AtomicParameters(Atoms.PeriodicCell(hcat(1)), [:H, :H, :C])
model = Models.Analytic.Harmonic(1.0, 1.0, 0.0)

Δ = Dict([(:H, 0.1), (:C, 0.1)])
sys = System{MonteCarlo}(atoms, model, 100u"K", Δ, 1; passes=1000)

@testset "propose_move!" begin
    Rᵢ = zeros(n_DoF(sys), n_atoms(sys))
    Rₚ = zero(Rᵢ)
    @test Rᵢ == Rₚ
    MetropolisHastings.propose_move!(sys, Rᵢ, Rₚ)
    @test Rᵢ != Rₚ
end

@testset "assess_proposal!" begin
    Rᵢ = randn(n_DoF(sys), n_atoms(sys))
    Rₚ = zero(Rᵢ)
    # Accept the move
    sys.dynamics.Eᵢ = 1
    output = MetropolisHastings.MonteCarloOutput{Float64}(Rᵢ, 1)
    MetropolisHastings.assess_proposal!(sys, Rᵢ, Rₚ, output, 1)
    @test output.R[1] == Rₚ
    # Reject the move
    sys.dynamics.Eᵢ = -1
    MetropolisHastings.assess_proposal!(sys, Rᵢ, Rₚ, output, 1)
    @test output.R[1] == Rᵢ
end

@testset "write_output!" begin
    Rₚ = fill(0.1, n_DoF(sys), n_atoms(sys))
    output = InitialConditions.MetropolisHastings.MonteCarloOutput{Float64}(Rₚ, 20)
    MetropolisHastings.write_output!(output, Rₚ, 1.0, 10)
    Rₚ .+= 1
    MetropolisHastings.write_output!(output, Rₚ, 1.1, 11)
    @test all(output.R[10] .== 0.1)
    @test output.energy[10] == 1.0
    @test all(output.R[11] .== 1.1)
    @test output.energy[11] == 1.1
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
    MetropolisHastings.apply_cell_boundaries!(cell, A, n_atoms)
    @test R == A # Check unchanged when inside cell
    A += 2rand(3, n_atoms) # Move atoms out of cell
    MetropolisHastings.apply_cell_boundaries!(cell, A, n_atoms)
    @test all(0 .<= A .<= 1) # Check they're all back in
    A -= 2rand(3, n_atoms) # Move atoms out of cell
    MetropolisHastings.apply_cell_boundaries!(cell, A, n_atoms)
    @test all(0 .<= A .<= 1) # Check they're all back in
end

@testset "run_monte_carlo_sampling" begin
    R0 = rand(n_DoF(sys), n_atoms(sys))
    out = InitialConditions.run_monte_carlo_sampling(sys, R0)
    @test !(out.R[1] ≈ out.R[10])
    @test !(out.energy[1] ≈ out.energy[20])
end
