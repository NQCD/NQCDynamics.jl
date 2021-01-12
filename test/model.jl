using Test
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.Models
using LinearAlgebra
using FiniteDiff

function finite_difference_gradient(model, R)
    function f(x)
        out = [0.0]
        potential!(model, out, x)
        out[1]
    end
    FiniteDiff.finite_difference_gradient(f, R)
end

R = rand(3, 10)

@testset "Harmonic" begin
    model = Harmonic()
    potential!(model, [0.0], R)
    D = zero(R)
    derivative!(model, D, R)
    @test finite_difference_gradient(model, R) ≈ D
end

@testset "DiatomicHarmonic" begin
    model = DiatomicHarmonic()
    R = rand(3, 2)
    V = [0.0]
    potential!(model, V, R)

    D = zero(R)
    derivative!(model, D, R)
    @test finite_difference_gradient(model, R) ≈ D
end

@testset "Free" begin
    model = Free()
    potential!(model, [0.0], R)
    D = zero(R)
    derivative!(model, D, R)
    @test finite_difference_gradient(model, R) ≈ D
end

@testset "DoubleWell" begin
    R = rand(1,1)
    model = DoubleWell()
    potential!(model, Hermitian(zeros(model.n_states, model.n_states)), R)
    derivative!(model, [Hermitian(zeros(model.n_states, model.n_states))]', R)
end

@testset "TullyModelOne" begin
    R = rand(1,1)
    model = TullyModelOne()
    potential!(model, Hermitian(zeros(model.n_states, model.n_states)), R)
    derivative!(model, [Hermitian(zeros(model.n_states, model.n_states))]', R)
end

@testset "TullyModelTwo" begin
    R = rand(1,1)
    model = TullyModelTwo()
    potential!(model, Hermitian(zeros(model.n_states, model.n_states)), R)
    derivative!(model, [Hermitian(zeros(model.n_states, model.n_states))]', R)
end

@testset "FrictionHarmonic" begin
    R = rand(1,3)
    model = FrictionHarmonic()
    potential!(model, [0.0], R)

    D = zero(R)
    derivative!(model, D, R)
    @test finite_difference_gradient(model, R) ≈ D

    friction!(model, zero(R), R)
end

@testset "ScatteringAndersonHolstein" begin
    R = rand(1,1)
    model = ScatteringAndersonHolstein()
    potential!(model, Hermitian(zeros(model.n_states, model.n_states)), R)
    derivative!(model, [Hermitian(zeros(model.n_states, model.n_states))]', R)
end

@testset "PdH" begin
    model = PdH([:Pd, :Pd, :H], PeriodicCell([10 0 0; 0 10 0; 0 0 10]), 10.0)
    R = rand(3, 3) * 10
    V = [0.0]
    D = zero(R)
    potential!(model, V, R)
    derivative!(model, D, R)
    @test finite_difference_gradient(model, R) ≈ D
end