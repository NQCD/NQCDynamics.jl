using Test
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.Models
using LinearAlgebra

R = rand(3, 10)

@testset "Harmonic" begin
    model = Harmonic()
    potential!(model, [0.0], R)
    derivative!(model, zero(R), R)
end

@testset "Free" begin
    model = Free()
    potential!(model, [0.0], R)
    derivative!(model, zero(R), R)
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
    derivative!(model, zero(R), R)
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
    h = 1e-4
    V1 = [0.0]
    for i=1:length(R)
        R[i] += h
        potential!(model, V1, R)
        @test D[i] ≈ (V1 - V)[1] / h rtol=1e-1
        R[i] -= h
    end
end