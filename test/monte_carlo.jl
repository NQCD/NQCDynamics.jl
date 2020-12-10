using Test
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.InitialConditions
using Unitful

atoms = Atoms.AtomicParameters(Atoms.InfiniteCell(), fill(:H, 3))
model = Models.Analytic.Harmonic(1.0, 1.0, 0.1)
sys = System(atoms, model, 1; temperature=100u"K")

@testset "propose_move!" begin
    Rᵢ = zeros(n_DoF(sys), n_atoms(sys))
    Rₚ = zero(Rᵢ)
    @test Rᵢ == Rₚ
    MetropolisHastings.propose_move!(Rᵢ, Rₚ, 1.0, n_DoF(sys), n_atoms(sys))
    @test Rᵢ != Rₚ
end

@testset "write_output!" begin
    Rₚ = fill(0.1, n_DoF(sys), n_atoms(sys))
    output = zeros(n_DoF(sys), n_atoms(sys))
    MetropolisHastings.write_output!(output, Rₚ)
    @test all(output .== 0.1)
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

@testset "run_monte_carlo_sampling" begin
    run_monte_carlo_sampling(sys, zeros(1, 3); passes=5000, Δ=0.1)
end
