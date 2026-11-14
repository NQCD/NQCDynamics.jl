using Test
using NQCModels
using NQCModels.BathDiscretisations: lorentzian_integral
using LinearAlgebra: Hermitian, diagind

const M = 20
const N = 200

schemes = (
    TrapezoidalRule(M, -1.0, 1.0),
    ShenviGaussLegendre(M, -1.0, 1.0),
    WindowedTrapezoidalRule(M, -1.0, 1.0, -0.2, 0.2),
)

@testset "lorentzian_integral" begin
    @test lorentzian_integral((-Inf, Inf), 0.1) ≈ 1.0
    @test lorentzian_integral((0.0, Inf), 0.1) ≈ 0.5
end

@testset "$(nameof(typeof(disc)))" for disc in schemes
    @testset "widebandbath (default)" begin
        bath = discrete_bath(disc)
        @test bath isa DiscreteBath
        @test NQCModels.nstates(bath) == M
        @test bath.bathstates == disc.bathstates
        @test bath.bathcoupling == disc.bathcoupling
        @test all(==(1.0), bath.bathdegeneracy)
        @test bath.bathtype == :widebandbath
        @test bath.discretisationtype == disc.discretisationtype
        @test bath.N == M
        @test bath.M == M
    end

    @testset "widebandbath with N electronic states" begin
        bath = discrete_bath(disc, widebandbath(disc, N))
        @test bath.N == N
        @test bath.M == M
        @test sum(bath.bathdegeneracy) ≈ N
    end

    @testset "lorentzianbath" begin
        bath = discrete_bath(disc, lorentzianbath(disc, N; W = 4.5))
        @test bath.bathtype == :lorentzianbath
        @test bath.N == N
        @test length(bath.bathfunction) == M
        @test all(>(0), bath.bathdegeneracy)
    end

    @testset "Hamiltonian filling" begin
        # Wide band DiscreteBath must reproduce the discretisation scheme on its own
        out_scheme = Hermitian(zeros(M + 1, M + 1))
        out_bath = Hermitian(zeros(M + 1, M + 1))
        bath = discrete_bath(disc)
        fillbathstates!(out_scheme, disc)
        fillbathstates!(out_bath, bath)
        @test out_bath[diagind(out_bath)[2:end]] == disc.bathstates
        fillbathcoupling!(out_scheme, 0.5, disc)
        fillbathcoupling!(out_bath, 0.5, bath)
        @test out_bath ≈ out_scheme

        # Degeneracy scales the couplings
        bath_N = discrete_bath(disc, widebandbath(disc, N))
        out = Hermitian(zeros(M + 1, M + 1))
        fillbathcoupling!(out, 0.5, bath_N, 2.0)
        @test out[1, 2:end] ≈ disc.bathcoupling .* 0.5 .* bath_N.bathdegeneracy .* 2.0
    end

    @testset "AndersonHolstein with DiscreteBath" begin
        impurity = ErpenbeckThoss(Γ = 2.0)
        r = hcat(0.5)
        @test potential(AndersonHolstein(impurity, discrete_bath(disc)), r) ≈
              potential(AndersonHolstein(impurity, disc), r)
    end
end

@testset "Lorentzian degeneracy bounded by N" begin
    disc = TrapezoidalRule(M, -1.0, 1.0)
    bath = discrete_bath(disc, lorentzianbath(disc, N; W = 4.5))
    @test sum(bath.bathdegeneracy) < N
end

@testset "Schemes not yet migrated to BathDiscretisationScheme" begin
    # Known issues carried over from NQCD/NQCModels.jl#73
    @test_broken discrete_bath(ReferenceGaussLegendre(M, -1.0, 1.0)) isa DiscreteBath
    @test_broken discrete_bath(FullGaussLegendre(M, -1.0, 1.0)) isa DiscreteBath
    @test_broken discrete_bath(GapGaussLegendre(M, -1.0, 1.0, 0.1)) isa DiscreteBath
    @test_broken discrete_bath(GapTrapezoidalRule(M, -1.0, 1.0, 0.1)) isa DiscreteBath
end
